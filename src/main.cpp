#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <Eigen/Dense>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>

static void usage() {
	std::cerr
		<< "adfinder (load windows into matrix)\n"
		<< "Usage:\n"
		<< "  adfinder --vcf input.vcf[.gz] --out output_prefix [--min-abs-r 0.2]\n"
		<< "  --intra        Scan intrachromosomal pairs (default: interchrom only)\n"
		<< "  --max-windows N  Load at most N windows (default: all)\n"
		<< "  --hi FILE       User provided hybrid index file (TSV: sample<TAB>hi)\n"
		<< "  --block-size INT       Block size for processing (default: 1024)\n";

}

static bool load_hi_file(
	const std::string& path,
	const bcf_hdr_t* hdr,
	Eigen::VectorXf& h_out
) {
	std::ifstream in(path);
	if (!in) {
		std::cerr << "Error: cannot open HI file: " << path << "\n";
		return false;
	}

	int nsamples = bcf_hdr_nsamples(hdr);

	std::unordered_map<std::string, float> m;
	m.reserve((size_t)nsamples * 2);

	std::string line;
	while (std::getline(in, line)) {
		if (line.size() == 0) continue;
		if (line[0] == '#') continue;

		std::istringstream ss(line);
		std::string sample;
		std::string hi_str;

		if (!(ss >> sample >> hi_str))
			continue;

		try {
			float hi = std::stof(hi_str);
			m[sample] = hi;
		} catch (...) {
			continue;
		}
	}

	h_out.resize(nsamples);

	int missing = 0;
	for (int i = 0; i < nsamples; ++i) {
		const char* s = hdr->samples[i];
		auto it = m.find(s);
		if (it == m.end()) {
			h_out(i) = std::numeric_limits<float>::quiet_NaN();
			++missing;
		} else {
			h_out(i) = it->second;
		}
	}

	if (missing > 0) {
		std::cerr << "Error: HI file missing " << missing << " / " << nsamples
				  << " samples (must contain all samples).\n";
		return false;
	}

	return true;
}

static bool write_hi_file(
	const std::string& path,
	const bcf_hdr_t* hdr,
	const Eigen::VectorXf& h
) {
	std::ofstream out(path);
	if (!out) {
		std::cerr << "Error: cannot write HI file: " << path << "\n";
		return false;
	}

	int nsamples = bcf_hdr_nsamples(hdr);
	out << "sample\thi\n";
	for (int i = 0; i < nsamples; ++i) {
		out << hdr->samples[i] << "\t" << h(i) << "\n";
	}
	return true;
}


int main(int argc, char** argv) {
	std::string vcf_path;
	std::string out;
	std::string hi_path;
	double min_abs_r = 0.2;
	bool intra = false;
	int max_windows = -1;	// -1 = all
	int block_size = 1024;

	for (int i = 1; i < argc; ++i) {
		std::string a = argv[i];
		if (a == "--help" || a == "-h") {
			usage();
			return 0;
		} else if (a == "--vcf" && i + 1 < argc) {
			vcf_path = argv[++i];
		} else if (a == "--out" && i + 1 < argc) {
			out = argv[++i];
		} else if (a == "--min-abs-r" && i + 1 < argc) {
			min_abs_r = std::stod(argv[++i]);
		} else if (a == "--intra") {
			intra = true;
		} else if (a == "--max-windows" && i + 1 < argc) {
			max_windows = std::stoi(argv[++i]);
		} else if (a == "--hi" && i + 1 < argc) {
			hi_path = argv[++i];
		} else if (a == "--block-size" && i + 1 < argc) {
			block_size = std::stoi(argv[++i]);
		} else {
			std::cerr << "Unknown/invalid arg: " << a << "\n";
			usage();
			return 2;
		}
	}


	if (vcf_path.empty() || out.empty()) {
		std::cerr << "Error: --vcf and --out are required.\n";
		usage();
		return 2;
	}

	htsFile* fp = bcf_open(vcf_path.c_str(), "r");
	if (!fp) {
		std::cerr << "Error: failed to open VCF/BCF: " << vcf_path << "\n";
		return 1;
	}

	bcf_hdr_t* hdr = bcf_hdr_read(fp);
	if (!hdr) {
		std::cerr << "Error: failed to read VCF header: " << vcf_path << "\n";
		bcf_close(fp);
		return 1;
	}

	const int nsamples = bcf_hdr_nsamples(hdr);
	std::cout << "adfinder VCF OK\n";
	std::cout << "  vcf        = " << vcf_path << "\n";
	std::cout << "  out        = " << out << "\n";
	std::cout << "  min_abs_r  = " << min_abs_r << "\n";
	std::cout << "  nsamples   = " << nsamples << "\n";

	int to_print = (nsamples < 5) ? nsamples : 5;
	for (int i = 0; i < to_print; ++i) {
		std::cout << "  sample[" << i << "] = " << hdr->samples[i] << "\n";
	}

	// ---- Decide DS vs GT (DS preferred) ----
	const int ds_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "DS");
	const int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
	const bool has_ds = (ds_id >= 0);
	const bool has_gt = (gt_id >= 0);

	if (!has_ds && !has_gt) {
		std::cerr << "Error: VCF has neither DS nor GT FORMAT field.\n";
		bcf_hdr_destroy(hdr);
		bcf_close(fp);
		return 1;
	}

	std::cout << "  using FORMAT/" << (has_ds ? "DS" : "GT") << "\n";

	// ---- Read first N windows into matrix X (nsamples × nwin) ----
	int max_windows_load = max_windows;
	if (max_windows_load <= 0)
		max_windows_load = 50000;	// effectively "all"


	std::vector<std::string> chroms;
	std::vector<int> starts;
	std::vector<int> ends;
	chroms.reserve(max_windows_load);
	starts.reserve(max_windows_load);
	ends.reserve(max_windows_load);

	// Ancestry matrix
	Eigen::MatrixXf X(nsamples, max_windows_load);

	const float NA = std::numeric_limits<float>::quiet_NaN();

	bcf1_t* rec = bcf_init();
	int nwin = 0;

	while (bcf_read(fp, hdr, rec) == 0 && nwin < max_windows_load) {

		bcf_unpack(rec, BCF_UN_STR);

		const char* chrom = bcf_hdr_id2name(hdr, rec->rid);
		const int start = rec->pos + 1;

		int end = start;
		int32_t* end_ptr = nullptr;
		int nend = 0;
		if (bcf_get_info_int32(hdr, rec, "END", &end_ptr, &nend) > 0 && nend > 0) {
			end = end_ptr[0];
		}
		free(end_ptr);

		chroms.emplace_back(chrom);
		starts.push_back(start);
		ends.push_back(end);

		if (has_ds) {
			float* ds = nullptr;
			int nds = 0;
			int ret = bcf_get_format_float(hdr, rec, "DS", &ds, &nds);
			if (ret <= 0) {
				X.col(nwin).setConstant(NA);
			} else {
				for (int i = 0; i < nsamples; ++i) {
					float v = ds[i];
					if (bcf_float_is_missing(v) || bcf_float_is_vector_end(v)) {
						X(i, nwin) = NA;
					} else {
						X(i, nwin) = v;
					}
				}
			}
			free(ds);
		} else {
			int32_t* gt = nullptr;
			int ngt = 0;
			int ret = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt);
			if (ret <= 0 || ngt < 2 * nsamples) {
				X.col(nwin).setConstant(NA);
			} else {
				for (int i = 0; i < nsamples; ++i) {
					int g0 = gt[2 * i];
					int g1 = gt[2 * i + 1];

					if (g0 == bcf_gt_missing || g1 == bcf_gt_missing) {
						X(i, nwin) = NA;
						continue;
					}
					int dosage = bcf_gt_allele(g0) + bcf_gt_allele(g1);
					X(i, nwin) = static_cast<float>(dosage);
				}
			}
			free(gt);
		}

		++nwin;
	}

	bcf_destroy(rec);

	// ---- Group window indices by chromosome ----
	std::unordered_map<std::string, std::vector<int>> windows_by_chr;
	windows_by_chr.reserve(64);

	std::vector<std::string> chr_order;
	chr_order.reserve(64);

	for (int w = 0; w < nwin; ++w) {
		auto it = windows_by_chr.find(chroms[w]);
		if (it == windows_by_chr.end()) {
			windows_by_chr[chroms[w]] = std::vector<int>();
			windows_by_chr[chroms[w]].reserve(1024);
			chr_order.push_back(chroms[w]);
		}
		windows_by_chr[chroms[w]].push_back(w);
	}

	std::cout << "Chromosomes in loaded windows: " << chr_order.size() << "\n";


	X.conservativeResize(Eigen::NoChange, nwin);

	std::cout << "Loaded " << nwin << " windows into matrix X ("
			  << nsamples << " × " << nwin << ")\n";

	int pr_s = std::min(3, nsamples);
	int pr_w = std::min(5, nwin);

	std::cout << "Preview (first " << pr_s << " samples × first "
			  << pr_w << " windows):\n";

	for (int i = 0; i < pr_s; ++i) {
		std::cout << "  " << hdr->samples[i] << " : ";
		for (int w = 0; w < pr_w; ++w) {
			float v = X(i, w);
			if (std::isnan(v)) std::cout << "NA";
			else std::cout << v;
			if (w + 1 < pr_w) std::cout << "\t";
		}
		std::cout << "\n";
	}

	std::cout << "First " << pr_w << " windows:\n";
	for (int w = 0; w < pr_w; ++w) {
		std::cout << "  w" << w << " "
				  << chroms[w] << ":" << starts[w] << "-" << ends[w] << "\n";
	}

	// ---- Compute hybrid index h (mean ancestry per individual) ----
	// For GT dosage in {0,1,2}, convert to proportion by dividing by 2.
	// Ignore NA entries.
	Eigen::VectorXf h;

	if (!hi_path.empty()) {
		std::cout << "Using HI from file: " << hi_path << "\n";
		if (!load_hi_file(hi_path, hdr, h)) {
			bcf_hdr_destroy(hdr);
			bcf_close(fp);
			return 1;
		}
	} else {
		std::cout << "Computing HI from X\n";
		h.resize(nsamples);

		for (int i = 0; i < nsamples; ++i) {
			double sum = 0.0;
			int count = 0;
			for (int w = 0; w < nwin; ++w) {
				float v = X(i, w);
				if (!std::isnan(v)) {
					sum += (double)v;
					++count;
				}
			}
			if (count == 0) {
				h(i) = std::numeric_limits<float>::quiet_NaN();
			} else {
				h(i) = (float)((sum / count) / 2.0);
			}
		}
	}

	// Always write the HI we used (computed or loaded)
	std::string hi_out_path = out + ".hi.tsv";
	if (!write_hi_file(hi_out_path, hdr, h)) {
		bcf_hdr_destroy(hdr);
		bcf_close(fp);
		return 1;
	}
	std::cout << "Wrote HI to: " << hi_out_path << "\n";


	// Summary stats (ignore NaNs)
	double hsum = 0.0;
	int hcount = 0;
	double hmin = 1e300;
	double hmax = -1e300;

	for (int i = 0; i < nsamples; ++i) {
		float v = h(i);
		if (!std::isnan(v)) {
			hsum += v;
			++hcount;
			if (v < hmin) hmin = v;
			if (v > hmax) hmax = v;
		}
	}

	std::cout << "Hybrid index (from X):\n";
	if (hcount == 0) {
		std::cout << "  (all missing)\n";
	} else {
		std::cout << "  n   = " << hcount << "\n";
		std::cout << "  min = " << hmin << "\n";
		std::cout << "  mean= " << (hsum / hcount) << "\n";
		std::cout << "  max = " << hmax << "\n";
	}

	// Preview first 5 h values
	int ph = std::min(15, nsamples);
	std::cout << "First " << ph << " HI values:\n";
	for (int i = 0; i < ph; ++i) {
		std::cout << "  " << hdr->samples[i] << "\t" << h(i) << "\n";
	}

	// ---- Residualize each window on HI and z-score => Z (nsamples × nwin) ----
	// Model per window: x_iw = a_w + b_w * h_i + e_iw
	// Z_iw = (e_iw - mean(e_.w)) / sd(e_.w)
	//
	// Note: For simplicity, if a window has ANY missing values, we mark it invalid (all NaN).

	Eigen::MatrixXf Z(nsamples, nwin);
	Z.setConstant(std::numeric_limits<float>::quiet_NaN());

	double hmean = 0.0;
	int hcount2 = 0;
	for (int i = 0; i < nsamples; ++i) {
		float hv = h(i);
		if (!std::isnan(hv)) {
			hmean += hv;
			++hcount2;
		}
	}
	if (hcount2 == 0) {
		std::cerr << "Error: HI is all missing; cannot residualize.\n";
		return 1;
	}
	hmean /= hcount2;

	// Centered HI vector
	Eigen::VectorXf hc(nsamples);
	for (int i = 0; i < nsamples; ++i) {
		float hv = h(i);
		if (std::isnan(hv)) hc(i) = std::numeric_limits<float>::quiet_NaN();
		else hc(i) = hv - static_cast<float>(hmean);
	}

	// Denominator: sum(hc^2) over non-missing HI
	double den = 0.0;
	for (int i = 0; i < nsamples; ++i) {
		float v = hc(i);
		if (!std::isnan(v)) den += static_cast<double>(v) * static_cast<double>(v);
	}
	if (den <= 0.0) {
		std::cerr << "Error: HI variance is zero; cannot residualize.\n";
		return 1;
	}

	int n_valid_windows = 0;

	for (int w = 0; w < nwin; ++w) {
		// Check missingness in this window (and HI)
		bool ok = true;

		// Compute mean(x) over individuals with non-missing x and HI
		double xmean = 0.0;
		int xcount = 0;

		for (int i = 0; i < nsamples; ++i) {
			float x = X(i, w);
			float hci = hc(i);
			if (std::isnan(x) || std::isnan(hci)) {
				ok = false;
				break;
			}
			xmean += x;
			++xcount;
		}

		if (!ok || xcount < 3) {
			// Leave Z(:,w) as NaN (invalid)
			continue;
		}

		xmean /= xcount;

		// Compute slope b = cov(xc, hc) / var(hc)
		double num = 0.0;
		for (int i = 0; i < nsamples; ++i) {
			double xc = static_cast<double>(X(i, w)) - xmean;
			double hcv = static_cast<double>(hc(i));
			num += xc * hcv;
		}
		double b = num / den;

		// Residuals: r = (x - xmean) - b * hc
		// Then z-score residuals
		double rmean = 0.0;
		for (int i = 0; i < nsamples; ++i) {
			double r = (static_cast<double>(X(i, w)) - xmean) - b * static_cast<double>(hc(i));
			rmean += r;
		}
		rmean /= nsamples;

		double rss = 0.0;
		for (int i = 0; i < nsamples; ++i) {
			double r = (static_cast<double>(X(i, w)) - xmean) - b * static_cast<double>(hc(i));
			double d = r - rmean;
			rss += d * d;
		}

		double sd = std::sqrt(rss / (nsamples - 1));
		if (!(sd > 0.0) || !std::isfinite(sd)) {
			// Constant residuals => invalid window
			continue;
		}

		for (int i = 0; i < nsamples; ++i) {
			double r = (static_cast<double>(X(i, w)) - xmean) - b * static_cast<double>(hc(i));
			Z(i, w) = static_cast<float>((r - rmean) / sd);
		}

		++n_valid_windows;
	}

	std::cout << "Residualization complete:\n";
	std::cout << "  valid_windows = " << n_valid_windows << " / " << nwin << "\n";

	// Quick sanity: print correlation between window 0 and 1 (if both valid)
	auto is_window_valid = [&](int w) -> bool {
		for (int i = 0; i < nsamples; ++i) {
			if (std::isnan(Z(i, w))) return false;
		}
		return true;
	};

	if (nwin >= 2 && is_window_valid(0) && is_window_valid(1)) {
		double dot = 0.0;
		for (int i = 0; i < nsamples; ++i) dot += static_cast<double>(Z(i, 0)) * static_cast<double>(Z(i, 1));
		double r = dot / (nsamples - 1);
		std::cout << "  sanity r(Z_w0, Z_w1) = " << r << "\n";
	} else {
		std::cout << "  sanity r(Z_w0, Z_w1) skipped (invalid window or <2 windows)\n";
	}

	// ---- Scan pairs using BLOCK matrix multiply (with window IDs; no END columns) ----
	std::string out_path = out + ".hits.tsv";
	std::ofstream of(out_path);
	if (!of) {
		std::cerr << "Error: cannot write to " << out_path << "\n";
		return 1;
	}

	of << "wA\tchrA\tstartA\twB\tchrB\tstartB\tr\tn\n";

	long long tested = 0;
	long long kept = 0;

	std::cout << "Scan mode: " << (intra ? "intrachromosomal" : "interchromosomal") << "\n";
	std::cout << "Block size: " << block_size << "\n";

	const float denom = 1.0f / (float)(nsamples - 1);

	// Reusable buffers (allocated to max block size once)
	Eigen::MatrixXf A(nsamples, block_size);
	Eigen::MatrixXf B(nsamples, block_size);
	Eigen::MatrixXf R(block_size, block_size);

	auto emit_hit = [&](int a, int b, float r) {
		of << a << "\t" << chroms[a] << "\t" << starts[a] << "\t"
		<< b << "\t" << chroms[b] << "\t" << starts[b] << "\t"
		<< r << "\t" << nsamples << "\n";
		++kept;
	};

	if (intra) {
		// Intrachromosomal: within each chromosome, scan all pairs (blockwise within same chr)
		for (const auto& chr : chr_order) {
			const auto& idx = windows_by_chr[chr];
			int m = (int)idx.size();

			for (int i0 = 0; i0 < m; i0 += block_size) {
				int b1 = std::min(block_size, m - i0);

				// Load block A
				for (int k = 0; k < b1; ++k)
					A.col(k) = Z.col(idx[i0 + k]);

				for (int j0 = i0; j0 < m; j0 += block_size) {
					int b2 = std::min(block_size, m - j0);

					// Load block B
					for (int k = 0; k < b2; ++k)
						B.col(k) = Z.col(idx[j0 + k]);

					// Compute correlations for this block
					R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
					R.topLeftCorner(b1, b2) *= denom;

					// Emit only upper triangle when j0 == i0 (avoid duplicates)
					for (int ia = 0; ia < b1; ++ia) {
						int a = idx[i0 + ia];

						int jb_start = 0;
						if (j0 == i0)
							jb_start = ia + 1;

						for (int ib = jb_start; ib < b2; ++ib) {
							int b = idx[j0 + ib];
							float r = R(ia, ib);

							++tested;

							if (std::fabs(r) >= min_abs_r)
								emit_hit(a, b, r);
						}
					}
				}
			}
		}
	} else {
		// Interchromosomal: across chromosome pairs, full rectangle (blockwise)
		int C = (int)chr_order.size();

		for (int c1 = 0; c1 < C; ++c1) {
			const auto& chr1 = chr_order[c1];
			const auto& idx1 = windows_by_chr[chr1];
			int m1 = (int)idx1.size();

			for (int c2 = c1 + 1; c2 < C; ++c2) {
				const auto& chr2 = chr_order[c2];
				const auto& idx2 = windows_by_chr[chr2];
				int m2 = (int)idx2.size();

				for (int i0 = 0; i0 < m1; i0 += block_size) {
					int b1 = std::min(block_size, m1 - i0);

					// Load block A
					for (int k = 0; k < b1; ++k)
						A.col(k) = Z.col(idx1[i0 + k]);

					for (int j0 = 0; j0 < m2; j0 += block_size) {
						int b2 = std::min(block_size, m2 - j0);

						// Load block B
						for (int k = 0; k < b2; ++k)
							B.col(k) = Z.col(idx2[j0 + k]);

						// Compute correlations for this block
						R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
						R.topLeftCorner(b1, b2) *= denom;

						// Threshold + write
						for (int ia = 0; ia < b1; ++ia) {
							int a = idx1[i0 + ia];
							for (int ib = 0; ib < b2; ++ib) {
								int b = idx2[j0 + ib];
								float r = R(ia, ib);

								++tested;

								if (std::fabs(r) >= min_abs_r)
									emit_hit(a, b, r);
							}
						}
					}
				}
			}
		}
	}

	of.close();

	std::cout << "Scan complete:\n";
	std::cout << "  tested_pairs = " << tested << "\n";
	std::cout << "  kept_pairs   = " << kept << " (|r| >= " << min_abs_r << ")\n";
	std::cout << "  wrote        = " << out_path << "\n";

	bcf_hdr_destroy(hdr);
	bcf_close(fp);
	return 0;
}
