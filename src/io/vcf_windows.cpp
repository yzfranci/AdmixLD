#include "vcf_windows.hpp"

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <cmath>
#include <cstdlib>
#include <limits>
#include <iostream>

WindowMatrix load_windows_from_vcf(
	const std::string& vcf_path,
	const VcfLoadOptions& opt
) {
	WindowMatrix out;

	htsFile* fp = bcf_open(vcf_path.c_str(), "r");
	if (!fp) {
		throw std::runtime_error("Failed to open VCF: " + vcf_path);
	}

	bcf_hdr_t* hdr = bcf_hdr_read(fp);
	if (!hdr) {
		bcf_close(fp);
		throw std::runtime_error("Failed to read VCF header: " + vcf_path);
	}

	int nsamples = bcf_hdr_nsamples(hdr);
	out.sample_names.reserve(nsamples);
	for (int i = 0; i < nsamples; ++i)
		out.sample_names.emplace_back(hdr->samples[i]);

	int ds_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "DS");
	int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
	bool has_ds = (ds_id >= 0);
	bool has_gt = (gt_id >= 0);
	if (!has_ds && !has_gt) {
		bcf_hdr_destroy(hdr);
		bcf_close(fp);
		throw std::runtime_error("VCF has neither DS nor GT FORMAT");
	}

	int max_windows_load = opt.max_windows;
	if (max_windows_load <= 0)
		max_windows_load = 50000;

	out.X = Eigen::MatrixXf(nsamples, max_windows_load);
	const float NA = std::numeric_limits<float>::quiet_NaN();

	out.meta.chrom.reserve(max_windows_load);
	out.meta.start.reserve(max_windows_load);
	out.meta.end.reserve(max_windows_load);

	bcf1_t* rec = bcf_init();
	int nwin = 0;

	while (bcf_read(fp, hdr, rec) == 0 && nwin < max_windows_load) {
		bcf_unpack(rec, BCF_UN_STR);

		const char* chrom = bcf_hdr_id2name(hdr, rec->rid);
		int start = rec->pos + 1;

		int end = start;
		int32_t* end_ptr = nullptr;
		int nend = 0;
		if (bcf_get_info_int32(hdr, rec, "END", &end_ptr, &nend) > 0 && nend > 0)
			end = end_ptr[0];
		free(end_ptr);

		out.meta.chrom.emplace_back(chrom);
		out.meta.start.push_back(start);
		out.meta.end.push_back(end);

		if (has_ds) {
			float* ds = nullptr;
			int nds = 0;
			int ret = bcf_get_format_float(hdr, rec, "DS", &ds, &nds);
			if (ret <= 0) {
				out.X.col(nwin).setConstant(NA);
			} else {
				for (int i = 0; i < nsamples; ++i) {
					float v = ds[i];
					if (bcf_float_is_missing(v) || bcf_float_is_vector_end(v))
						out.X(i, nwin) = NA;
					else
						out.X(i, nwin) = v;
				}
			}
			free(ds);
		} else {
			int32_t* gt = nullptr;
			int ngt = 0;
			int ret = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
			if (ret <= 0 || ngt < 2 * nsamples) {
				out.X.col(nwin).setConstant(NA);
			} else {
				for (int i = 0; i < nsamples; ++i) {
					int g0 = gt[2 * i];
					int g1 = gt[2 * i + 1];
					if (g0 == bcf_gt_missing || g1 == bcf_gt_missing) {
						out.X(i, nwin) = NA;
						continue;
					}
					int dosage = bcf_gt_allele(g0) + bcf_gt_allele(g1);
					out.X(i, nwin) = (float)dosage;
				}
			}
			free(gt);
		}

		++nwin;
	}

	bcf_destroy(rec);
	bcf_hdr_destroy(hdr);
	bcf_close(fp);

	out.X.conservativeResize(Eigen::NoChange, nwin);
	return out;
}
