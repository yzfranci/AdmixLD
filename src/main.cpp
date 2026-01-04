#include <Eigen/Dense>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <cstdint>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>

#include "io/vcf_windows.hpp"
#include "io/bed.hpp"
#include "core/hybrid_index.hpp"
#include "core/residualize.hpp"
#include "core/scan_blocks.hpp"
#include "config.hpp"

static void usage() {
	std::cerr
		<< "adfinder (load windows into matrix)\n"
		<< "Usage:\n"
		<< "  adfinder --vcf input.vcf[.gz] --out output_prefix\n"
		<< "  --min-abs-r FLOAT      Keep pairs with |r| >= value (default)\n"
		<< "  --min-neg-r FLOAT      Keep pairs with r <= -value (asymmetric)\n"
		<< "  --min-pos-r FLOAT      Keep pairs with r >= value (asymmetric)\n"
		<< "  --intra            Scan intrachromosomal pairs (default: interchrom only)\n"
		<< "  --max-dist INT     Intra only: max end-distance (bp) between window pairs\n"
		<< "  --max-windows N    Load at most N windows (default cap if <=0)\n"
		<< "  --hi FILE          User provided hybrid index file (TSV: sample<TAB>hi)\n"
		<< "  --full-hi           Compute hybrid index using the full VCF (ignore --chr/--bed)\n"
		<< "  --block-size INT   Block size for processing (default: 1024)\n"
		<< "  --distrib              Write empirical scan r distribution summary\n"
		<< "  --distrib-sample INT   Reservoir sample size for distrib summary (default: 200000)\n"
		<< "  --permute N            Run N interchrom chr-block permutations (summary stats)\n"
		<< "  --permute-sample INT   Reservoir sample size for percentile estimates (default: 200000)\n"
		<< "  --seed INT             RNG seed for permutations (default: 1)\n"
		<< "  --chr STR              Keep only this chromosome (repeatable)\n"
		<< "  --bed FILE             Keep windows whose pos is within BED intervals (chr start end; no header)\n"
		<< "  --target-chr STR        Scan one target window/pos vs all others (target chromosome)\n"
		<< "  --target-pos INT        Scan one target window/pos vs all others (target position; matches single pos column)\n";
}

int main(int argc, char** argv) {
	std::string vcf_path;
	std::string out;
	std::string hi_path;
	std::vector<std::string> keep_chrs;
	std::string bed_path;
	std::string target_chr;

	bool full_hi = false;
	double min_abs_r = 0;
	bool intra = false;
	int max_dist = -1;
	int max_windows = -1;	// <=0 means loader default cap
	int block_size = ADFINDER_DEFAULT_BLOCK_SIZE;
	int n_perm = 0;
	int perm_sample = 200000;
	uint64_t seed = 1;
	bool distrib = false;
	int distrib_sample = 200000;
	bool has_target = false;
	int target_pos = -1;
	bool has_min_neg_r = false;
	bool has_min_pos_r = false;
	float min_neg_r = 0.0f;
	float min_pos_r = 0.0f;


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

		} else if (a == "--max-dist" && i + 1 < argc) {
			max_dist = std::stoi(argv[++i]);	

		} else if (a == "--max-windows" && i + 1 < argc) {
			max_windows = std::stoi(argv[++i]);

		} else if (a == "--hi" && i + 1 < argc) {
			hi_path = argv[++i];

		} else if (a == "--full-hi") {
			full_hi = true;

		} else if (a == "--block-size" && i + 1 < argc) {
			block_size = std::stoi(argv[++i]);

		} else if (a == "--permute" && i + 1 < argc) {
			n_perm = std::stoi(argv[++i]);

		} else if (a == "--permute-sample" && i + 1 < argc) {
			perm_sample = std::stoi(argv[++i]);

		} else if (a == "--seed" && i + 1 < argc) {
			seed = (uint64_t)std::stoull(argv[++i]);

		} else if (a == "--distrib") {
			distrib = true;

		} else if (a == "--distrib-sample" && i + 1 < argc) {
			distrib_sample = std::stoi(argv[++i]);
		
		} else if (a == "--chr" && i + 1 < argc) {
			keep_chrs.push_back(argv[++i]);

		} else if (a == "--bed" && i + 1 < argc) {
			bed_path = argv[++i];

		} else if (a == "--target-chr" && i + 1 < argc) {
			has_target = true;
			target_chr = argv[++i];

		} else if (a == "--target-pos" && i + 1 < argc) {
			has_target = true;
			target_pos = std::stoi(argv[++i]);

		} else if (a == "--min-neg-r" && i + 1 < argc) {
			min_neg_r = std::stof(argv[++i]);
			has_min_neg_r = true;

		} else if (a == "--min-pos-r" && i + 1 < argc) {
			min_pos_r = std::stof(argv[++i]);
			has_min_pos_r = true;

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

	// ---- Set asymetric threshold if min_neg_r and/or min_pos_r are set ----
	bool use_asym = false;

	if (has_min_neg_r || has_min_pos_r) {
		use_asym = true;

		if (!has_min_neg_r)
			min_neg_r = 0.0f;
		if (!has_min_pos_r)
			min_pos_r = 0.0f;

		if (min_neg_r < 0.0f || min_pos_r < 0.0f) {
			std::cerr << "Error: --min-neg-r and --min-pos-r must be >= 0\n";
			return 2;
		}
	}


	// ---- load windows via module ----
	VcfLoadOptions vopt;
	vopt.max_windows = max_windows;

	WindowMatrix wm = load_windows_from_vcf(vcf_path, vopt);

	const int nsamples = (int)wm.sample_names.size();

	Eigen::MatrixXf X = wm.X;
	Eigen::MatrixXf X_full = X;

	std::vector<std::string> chroms = wm.meta.chrom;
	std::vector<int> starts = wm.meta.start;
	std::vector<int> ends = wm.meta.end;
	int nwin = (int)chroms.size();

	// ---- Optional filters (chr / bed) applied BEFORE residualization ----
	std::unordered_map<std::string, bool> chr_keep_map;
	if (!keep_chrs.empty()) {
		chr_keep_map.reserve(keep_chrs.size() * 2 + 1);
		for (const auto& c : keep_chrs)
			chr_keep_map[c] = true;
	}

	std::unordered_map<std::string, std::vector<BedInterval>> bed_by_chr;
	bool use_bed = false;
	if (!bed_path.empty()) {
		use_bed = true;
		if (!read_bed(bed_path, bed_by_chr))
			return 1;
		std::cout << "Using BED filter: " << bed_path << "\n";
	}

	std::vector<int> keep_idx;
	keep_idx.reserve((size_t)nwin);

	for (int w = 0; w < nwin; ++w) {
		const std::string& chr = chroms[w];
		int pos = ends[w];

		if (!keep_chrs.empty()) {
			if (chr_keep_map.find(chr) == chr_keep_map.end())
				continue;
		}

		if (use_bed) {
			if (!bed_contains(bed_by_chr, chr, pos))
				continue;
		}

		keep_idx.push_back(w);
	}

	if (keep_idx.empty()) {
		std::cerr << "Error: no windows remain after filtering.\n";
		return 1;
	}

	if ((int)keep_idx.size() != nwin) {
		std::cout << "Filtering windows:\n";
		std::cout << "  before = " << nwin << "\n";
		std::cout << "  after  = " << (int)keep_idx.size() << "\n";

		int nwin2 = (int)keep_idx.size();

		Eigen::MatrixXf Xf(nsamples, nwin2);
		std::vector<std::string> chroms_f((size_t)nwin2);
		std::vector<int> starts_f((size_t)nwin2);
		std::vector<int> ends_f((size_t)nwin2);

		for (int j = 0; j < nwin2; ++j) {
			int w = keep_idx[j];
			Xf.col(j) = X.col(w);
			chroms_f[j] = chroms[w];
			starts_f[j] = starts[w];
			ends_f[j] = ends[w];
		}

		X.swap(Xf);
		chroms.swap(chroms_f);
		starts.swap(starts_f);
		ends.swap(ends_f);
		nwin = nwin2;
	}


	std::cout << "adfinder VCF OK\n";
	std::cout << "  vcf        = " << vcf_path << "\n";
	std::cout << "  out        = " << out << "\n";
	std::cout << "  nsamples   = " << nsamples << "\n";
	std::cout << "  windows    = " << nwin << "\n";

	int to_print = std::min(5, nsamples);
	for (int i = 0; i < to_print; ++i)
		std::cout << "  sample[" << i << "] = " << wm.sample_names[i] << "\n";
	
	std::cout << std::flush;


	// ---- Group window indices by chromosome ----
	std::vector<std::string> chr_order;
	auto windows_by_chr = group_by_chr(chroms, chr_order);
	std::cout << "Chromosomes in loaded windows: " << chr_order.size() << "\n";

	// Ensure within-chromosome ordering by position (required for --max-dist upper_bound)
	for (auto& kv : windows_by_chr) {
		auto& idx = kv.second;
		std::sort(idx.begin(), idx.end(),
			[&](int a, int b) { return ends[a] < ends[b]; }
		);
	}

	// ---- Compute / load hybrid index ----
	Eigen::VectorXf h;

	if (full_hi && !hi_path.empty()) {
		std::cerr << "Warning: --full-hi ignored because --hi FILE was provided\n";
	}

	if (!hi_path.empty()) {
		std::cout << "Using HI from file: " << hi_path << "\n";
		if (!load_hi_tsv(hi_path, wm.sample_names, h))
			return 1;
	} else if (full_hi) {
		std::cout << "Computing HI from FULL VCF (pre-filter)\n";
		h = compute_hi_from_X(X_full);
	} else {
		std::cout << "Computing HI from FILTERED windows (if filters such as --chr or --bed are used)\n";
		h = compute_hi_from_X(X);
	}


	// Always write the HI we used (computed or loaded)
	std::string hi_out_path = out + ".hi.tsv";
	if (!write_hi_tsv(hi_out_path, wm.sample_names, h))
		return 1;
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

	int ph = std::min(15, nsamples);
	std::cout << "First " << ph << " HI values:\n";
	for (int i = 0; i < ph; ++i)
		std::cout << "  " << wm.sample_names[i] << "\t" << h(i) << "\n";

	std::cout << std::flush;

	// ---- Resolve target window index ----
	int target_w = -1;
	if (has_target) {
		if (target_chr.empty() || target_pos < 0) {
			std::cerr << "Error: --target-chr and --target-pos must both be provided.\n";
			return 1;
		}

		for (int w = 0; w < nwin; ++w) {
			if (chroms[w] == target_chr && ends[w] == target_pos) {
				target_w = w;
				break;
			}
		}

		if (target_w < 0) {
			std::cerr << "Error: target window not found after filtering: "
				<< target_chr << ":" << target_pos << "\n";
			return 1;
		}

		std::cout << "Target window: w=" << target_w
			<< " chr=" << target_chr
			<< " pos=" << target_pos << "\n";
	}

	// ---- Residualize each window on HI and z-score => Z (nsamples × nwin) ----
	int n_valid_windows = 0;
	Eigen::MatrixXf Z;

	try {
		Z = residualize_and_zscore(X, h, n_valid_windows);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return 1;
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
		for (int i = 0; i < nsamples; ++i)
			dot += (double)Z(i, 0) * (double)Z(i, 1);
		double r = dot / (nsamples - 1);
		std::cout << "  sanity r(Z_w0, Z_w1) = " << r << "\n";
	} else {
		std::cout << "  sanity r(Z_w0, Z_w1) skipped (invalid window or <2 windows)\n";
	}

	std::cout << std::flush;

	// ---- Permutation test ----
	if (n_perm > 0) {
		if (intra) {
			std::cerr << "Error: --permute is currently implemented for interchrom scans only (omit --intra).\n";
			return 1;
		}

		ScanOptions popt;
		popt.intra = intra;
		popt.max_dist = max_dist;
		popt.block_size = block_size;
		popt.nsamples = nsamples;
		popt.min_abs_r = (float)min_abs_r;
		popt.use_asym = false;
		popt.min_neg_r = 0.0f;
		popt.min_pos_r = 0.0f;

		std::vector<PermSummary> summ;

		std::cout << "Permutation test (interchrom, chr-block):\n";
		std::cout << "  n_perm         = " << n_perm << "\n";
		std::cout << "  seed           = " << seed << "\n";
		std::cout << "  perm_sample    = " << perm_sample << "\n";

		std::cout << std::flush;

		if (!permute_interchrom_summary_chrblock(Z, windows_by_chr, chr_order, popt, seed, n_perm, perm_sample, summ))
			return 1;

		std::string perm_path = out + ".perm.summary.tsv";
		std::ofstream pf(perm_path);
		if (!pf) {
			std::cerr << "Error: cannot write to " << perm_path << "\n";
			return 1;
		}

		pf << "rep\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\n";

		for (int r = 0; r < (int)summ.size(); ++r) {
			pf << r
				<< "\t" << summ[r].max_r
				<< "\t" << summ[r].p99
				<< "\t" << summ[r].p95
				<< "\t" << summ[r].p75
				<< "\t" << summ[r].median
				<< "\t" << summ[r].p25
				<< "\t" << summ[r].p05
				<< "\t" << summ[r].p01
				<< "\t" << summ[r].min_r
				<< "\n";
		}
		pf.close();

		std::cout << "  wrote permutation summaries: " << perm_path << "\n";

		// If permutations were requested, stop here (do not run hits scan)
		std::cout << "Permutation-only mode (--permute was set). Skipping hits scan.\n";
		return 0;
	}


	// ---- Scan pairs using BLOCK matrix multiply ----
	std::string out_path = out + ".hits.tsv";

	ScanOptions opt;
	opt.intra = intra;
	opt.max_dist = max_dist;
	opt.block_size = block_size;
	opt.nsamples = nsamples;
	opt.min_abs_r = (float)min_abs_r;
	opt.use_asym = use_asym;
	opt.min_neg_r = min_neg_r;
	opt.min_pos_r = min_pos_r;


	long long tested = 0;
	long long kept = 0;

	std::cout << "Scan mode: " << (intra ? "intrachromosomal" : "interchromosomal") << "\n";
	std::cout << "Block size: " << block_size << "\n";

	std::string distrib_path;
	if (distrib)
		distrib_path = out + ".scan.summary.tsv";

	if (target_w >= 0) {
		if (!scan_target_write_hits(
			Z, chroms, ends, windows_by_chr, chr_order,
			opt, target_w, out_path, tested, kept,
			distrib_path, distrib_sample, seed
		))
			return 1;
	} else {
		if (!scan_blocks_write_hits(
			Z, chroms, ends, windows_by_chr, chr_order,
			opt, out_path, tested, kept,
			distrib_path, distrib_sample, seed
		))
			return 1;
	}

	std::cout << "Scan complete:\n";
	std::cout << "  tested_pairs = " << tested << "\n";
	std::cout << "  kept_pairs   = " << kept << " (|r| >= " << min_abs_r << ")\n";
	std::cout << "  wrote        = " << out_path << "\n";
	if (distrib)
		std::cout << "  wrote        = " << distrib_path << "\n";
	if (use_asym) {
		std::cout << "  r filter    = ";
		if (has_min_neg_r)
			std::cout << "r <= -" << min_neg_r << " ";
		if (has_min_pos_r)
			std::cout << "r >= " << min_pos_r;
		std::cout << "\n";
	} else {
		std::cout << "  |r| filter  = " << min_abs_r << "\n";
	}

	return 0;
}
