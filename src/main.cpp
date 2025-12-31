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

#include "io/vcf_windows.hpp"
#include "core/hybrid_index.hpp"
#include "core/residualize.hpp"
#include "core/scan_blocks.hpp"
#include "config.hpp"

static void usage() {
	std::cerr
		<< "adfinder (load windows into matrix)\n"
		<< "Usage:\n"
		<< "  adfinder --vcf input.vcf[.gz] --out output_prefix [--min-abs-r 0.2]\n"
		<< "  --intra            Scan intrachromosomal pairs (default: interchrom only)\n"
		<< "  --max-windows N    Load at most N windows (default cap if <=0)\n"
		<< "  --hi FILE          User provided hybrid index file (TSV: sample<TAB>hi)\n"
		<< "  --block-size INT   Block size for processing (default: 1024)\n"
		<< "  --permute N            Run N interchrom chr-block permutations (summary stats)\n"
		<< "  --permute-sample INT   Reservoir sample size for percentile estimates (default: 200000)\n"
		<< "  --seed INT             RNG seed for permutations (default: 1)\n";
}

int main(int argc, char** argv) {
	std::string vcf_path;
	std::string out;
	std::string hi_path;

	double min_abs_r = 0;
	bool intra = false;
	int max_windows = -1;	// <=0 means loader default cap
	int block_size = ADFINDER_DEFAULT_BLOCK_SIZE;
	int n_perm = 0;
	int perm_sample = 200000;
	uint64_t seed = 1;


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

		} else if (a == "--permute" && i + 1 < argc) {
			n_perm = std::stoi(argv[++i]);

		} else if (a == "--permute-sample" && i + 1 < argc) {
			perm_sample = std::stoi(argv[++i]);

		} else if (a == "--seed" && i + 1 < argc) {
			seed = (uint64_t)std::stoull(argv[++i]);

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

	// ---- NEW: load windows via module ----
	VcfLoadOptions vopt;
	vopt.max_windows = max_windows;

	WindowMatrix wm = load_windows_from_vcf(vcf_path, vopt);

	const int nsamples = (int)wm.sample_names.size();
	const int nwin = (int)wm.X.cols();

	Eigen::MatrixXf X = wm.X;

	const std::vector<std::string>& chroms = wm.meta.chrom;
	const std::vector<int>& starts = wm.meta.start;
	const std::vector<int>& ends = wm.meta.end;

	std::cout << "adfinder VCF OK\n";
	std::cout << "  vcf        = " << vcf_path << "\n";
	std::cout << "  out        = " << out << "\n";
	std::cout << "  min_abs_r  = " << min_abs_r << "\n";
	std::cout << "  nsamples   = " << nsamples << "\n";
	std::cout << "  windows    = " << nwin << "\n";

	int to_print = std::min(5, nsamples);
	for (int i = 0; i < to_print; ++i)
		std::cout << "  sample[" << i << "] = " << wm.sample_names[i] << "\n";

	// ---- Group window indices by chromosome ----
		std::vector<std::string> chr_order;
		auto windows_by_chr = group_by_chr(chroms, chr_order);
		std::cout << "Chromosomes in loaded windows: " << chr_order.size() << "\n";

	// ---- Compute / load hybrid index ----
	Eigen::VectorXf h;

	if (!hi_path.empty()) {
		std::cout << "Using HI from file: " << hi_path << "\n";
		if (!load_hi_tsv(hi_path, wm.sample_names, h))
			return 1;
	} else {
		std::cout << "Computing HI from X\n";
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

	// ---- Permutation test ----
	if (n_perm > 0) {
		if (intra) {
			std::cerr << "Error: --permute is currently implemented for interchrom scans only (omit --intra).\n";
			return 1;
		}

		ScanOptions popt;
		popt.intra = false;
		popt.block_size = block_size;
		popt.min_abs_r = (float)min_abs_r;
		popt.nsamples = nsamples;

		std::vector<PermSummary> summ;

		std::cout << "Permutation test (interchrom, chr-block):\n";
		std::cout << "  n_perm         = " << n_perm << "\n";
		std::cout << "  seed           = " << seed << "\n";
		std::cout << "  perm_sample    = " << perm_sample << "\n";

		if (!permute_interchrom_summary_chrblock(Z, windows_by_chr, chr_order, popt, seed, n_perm, perm_sample, summ))
			return 1;

		std::string perm_path = out + ".perm.summary.tsv";
		std::ofstream pf(perm_path);
		if (!pf) {
			std::cerr << "Error: cannot write to " << perm_path << "\n";
			return 1;
		}

		pf << "rep\tmax_r\tp99\tp95\tmedian\tp05\tp01\tmin_r\n";
		for (int r = 0; r < (int)summ.size(); ++r) {
			pf << r
				<< "\t" << summ[r].max_r
				<< "\t" << summ[r].p99
				<< "\t" << summ[r].p95
				<< "\t" << summ[r].median
				<< "\t" << summ[r].p05
				<< "\t" << summ[r].p01
				<< "\t" << summ[r].min_r
				<< "\n";
		}
		pf.close();

		std::cout << "  wrote permutation summaries: " << perm_path << "\n";
	}


	// ---- Scan pairs using BLOCK matrix multiply ----
	std::string out_path = out + ".hits.tsv";

	ScanOptions opt;
	opt.intra = intra;
	opt.block_size = block_size;
	opt.min_abs_r = (float)min_abs_r;
	opt.nsamples = nsamples;

	long long tested = 0;
	long long kept = 0;

	std::cout << "Scan mode: " << (intra ? "intrachromosomal" : "interchromosomal") << "\n";
	std::cout << "Block size: " << block_size << "\n";

	if (!scan_blocks_write_hits(Z, chroms, starts, windows_by_chr, chr_order, opt, out_path, tested, kept))
		return 1;

	std::cout << "Scan complete:\n";
	std::cout << "  tested_pairs = " << tested << "\n";
	std::cout << "  kept_pairs   = " << kept << " (|r| >= " << min_abs_r << ")\n";
	std::cout << "  wrote        = " << out_path << "\n";

	return 0;
}
