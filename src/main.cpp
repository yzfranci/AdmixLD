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
#include <unordered_set>

#include "io/vcf_windows.hpp"
#include "io/msp_windows.hpp"
#include "io/bed.hpp"
#include "io/sample_vector.hpp"
#include "core/hybrid_index.hpp"
#include "core/residualize.hpp"
#include "core/scan_blocks.hpp"
#include "core/scan_dynamic_hi.hpp"
#include "config.hpp"

/*
	AdmixLD main entry point

	Pipeline:
	1) Parse/validate CLI
	2) Load ancestry blocks from VCF
	3) Apply optional block filters (--chr/--no-chr/--bed)
	4) Build chromosome index (block indices by chromosome)
	5) Compute/load hybrid index (full-genome by default; filtered if --compute-hi)
	6) Optional modes: HI-only, permutations, sample-geno scans
	7) Residualize blocks on HI and scan for LD hits
*/

static void usage() {
	std::cerr
		<< "admixld (load blocks into matrix)\n"
		<< "Usage:\n"
		<< "  admixld --vcf input.vcf[.gz] --out output_prefix\n"
		<< "  admixld --msp input.msp.tsv   --out output_prefix\n"
		<< "  --min-abs-r FLOAT      Keep pairs with |r| >= value (default)\n"
		<< "  --min-neg-r FLOAT      Keep pairs with r <= -value (asymmetric)\n"
		<< "  --min-pos-r FLOAT      Keep pairs with r >= value (asymmetric)\n"
		<< "  --intra                Scan intrachromosomal pairs (default: interchrom only)\n"
		<< "  --phased               Use phased haplotypes (2n rows); requires --intra; MSP/VCF GT only\n"
		<< "  --max-dist INT         Intra only: max end-distance (bp) between block pairs\n"
		<< "  --max-windows N        Load at most N blocks from the VCF (used solely for test runs and debugging)\n"
		<< "  --min-callrate FLOAT   Minimum block call rate; values <1.0 enable within-block mean imputation (beta; might introduce bias).\n"
		<< "  --cov FILE             Covariate file (TSV: for HI: sample<TAB>hi; For more covariates (PCA etc...) sample<TAB>cov1 [cov2 ...]; header optional); overrides computed HI\n"
		<< "  --hi-mode STR          HI correction: global (default) | excl-focus (LOCO/LOCO2; incompatible with --cov)\n"
		<< "  --compute-hi           Compute HI using FILTERED blocks only, write out.hi.tsv, and exit (ignored with --cov)\n"
		<< "  --unweighted-hi        Use unweighted HI (mean(dosage)/2; legacy behavior; ignored with --cov)\n"
		<< "  --pos-is-start         For weighted HI: interpret VCF block pos as START (ignored with --cov)\n"
		<< "  --tile-size INT		 tile size for processing (default: 1024)\n"
		<< "  --threads INT          Number of OpenMP threads for scan steps (default: 1)\n"
		<< "  --distrib              Write empirical scan distribution summary\n"
		<< "  --distrib-sample INT   Reservoir sample size for distribution summary (default: 200000)\n"
		<< "  --permute N            Run N interchrom full-shuffle permutations (summary stats)\n"
		<< "  --permute-sample INT   Reservoir sample size for percentile estimates (default: 200000)\n"
		<< "  --seed INT             RNG seed for permutations (default: 1)\n"
		<< "  --chr STR              Keep only this chromosome (repeatable)\n"
		<< "  --no-chr STR           Exclude this chromosome (repeatable; opposite of --chr)\n"
		<< "  --bed FILE             Keep blocks whose position is within BED intervals (chr start end; no header)\n"
		<< "  --target-chr STR       Scan one target block/pos vs all others (target chromosome)\n"
		<< "  --target-pos INT       Scan one target block/pos vs all others (target position; matches single position)\n"
		<< "  --sample-geno FILE     Per-sample numeric mito haplotype/genotype/trait TSV: sample<TAB>value\n"
	<< "  --keep-indv FILE       Keep only samples listed in FILE (one ID per line)\n";
}

struct CliOptions {
	std::string vcf_path;
	std::string msp_path;
	std::string out;
	std::string cov_path;

	std::vector<std::string> keep_chrs;
	std::vector<std::string> drop_chrs;

	std::string bed_path;

	bool has_target = false;
	std::string target_chr;
	int target_pos = -1;

	std::string sample_geno_path;
	std::string keep_indv_path;

	std::string hi_mode = "global";	// global | excl-focus

	bool compute_hi_only = false;
	bool unweighted_hi = false;
	bool pos_is_start = false;

	bool intra = false;
	bool phased = false;
	int max_dist = -1;
	int max_windows = -1;	// <=0 means loader default cap
	double min_callrate = 1.0;	// fraction in [0,1]

	int tile_size = ADMIXLD_DEFAULT_TILE_SIZE;
	int threads = 1;

	bool distrib = false;
	int distrib_sample = 200000;

	int n_perm = 0;
	int perm_sample = 200000;
	uint64_t seed = 1;

	double min_abs_r = 0.0;
	bool has_min_neg_r = false;
	bool has_min_pos_r = false;
	float min_neg_r = 0.0f;
	float min_pos_r = 0.0f;
};

static bool parse_args(int argc, char** argv, CliOptions& opt) {
	for (int i = 1; i < argc; ++i) {
		std::string a = argv[i];

		if (a == "--help" || a == "-h") {
			usage();
			std::exit(0);

		} else if (a == "--vcf" && i + 1 < argc) {
			opt.vcf_path = argv[++i];

		} else if (a == "--msp" && i + 1 < argc) {
			opt.msp_path = argv[++i];

		} else if (a == "--out" && i + 1 < argc) {
			opt.out = argv[++i];

		} else if (a == "--min-abs-r" && i + 1 < argc) {
			opt.min_abs_r = std::stod(argv[++i]);

		} else if (a == "--intra") {
			opt.intra = true;

		} else if (a == "--phased") {
			opt.phased = true;

		} else if (a == "--max-dist" && i + 1 < argc) {
			opt.max_dist = std::stoi(argv[++i]);

		} else if (a == "--max-windows" && i + 1 < argc) {
			opt.max_windows = std::stoi(argv[++i]);

		} else if (a == "--min-callrate" && i + 1 < argc) {
			opt.min_callrate = std::stod(argv[++i]);

		} else if (a == "--cov" && i + 1 < argc) {
			opt.cov_path = argv[++i];

		} else if (a == "--hi-mode" && i + 1 < argc) {
			opt.hi_mode = argv[++i];

		} else if (a == "--compute-hi") {
			opt.compute_hi_only = true;

		} else if (a == "--unweighted-hi") {
			opt.unweighted_hi = true;

		} else if (a == "--pos-is-start") {
			opt.pos_is_start = true;

		} else if (a == "--tile-size" && i + 1 < argc) {
			opt.tile_size = std::stoi(argv[++i]);

		} else if (a == "--threads" && i + 1 < argc) {
			opt.threads = std::stoi(argv[++i]);

		} else if (a == "--permute" && i + 1 < argc) {
			opt.n_perm = std::stoi(argv[++i]);

		} else if (a == "--permute-sample" && i + 1 < argc) {
			opt.perm_sample = std::stoi(argv[++i]);

		} else if (a == "--seed" && i + 1 < argc) {
			opt.seed = (uint64_t)std::stoull(argv[++i]);

		} else if (a == "--distrib") {
			opt.distrib = true;

		} else if (a == "--distrib-sample" && i + 1 < argc) {
			opt.distrib_sample = std::stoi(argv[++i]);

		} else if (a == "--chr" && i + 1 < argc) {
			opt.keep_chrs.push_back(argv[++i]);

		} else if (a == "--no-chr" && i + 1 < argc) {
			opt.drop_chrs.push_back(argv[++i]);

		} else if (a == "--bed" && i + 1 < argc) {
			opt.bed_path = argv[++i];

		} else if (a == "--target-chr" && i + 1 < argc) {
			opt.has_target = true;
			opt.target_chr = argv[++i];

		} else if (a == "--target-pos" && i + 1 < argc) {
			opt.has_target = true;
			opt.target_pos = std::stoi(argv[++i]);

		} else if (a == "--sample-geno" && i + 1 < argc) {
			opt.sample_geno_path = argv[++i];

		} else if (a == "--keep-indv" && i + 1 < argc) {
			opt.keep_indv_path = argv[++i];

		} else if (a == "--min-neg-r" && i + 1 < argc) {
			opt.min_neg_r = std::stof(argv[++i]);
			opt.has_min_neg_r = true;

		} else if (a == "--min-pos-r" && i + 1 < argc) {
			opt.min_pos_r = std::stof(argv[++i]);
			opt.has_min_pos_r = true;

		} else {
			std::cerr << "Unknown/invalid arg: " << a << "\n";
			usage();
			return false;
		}
	}
	return true;
}

static int validate_options(const CliOptions& opt) {
	if (opt.vcf_path.empty() == opt.msp_path.empty()) {
		std::cerr << "Error: exactly one of --vcf or --msp is required.\n";
		usage();
		return 2;
	}

	if (opt.out.empty()) {
		std::cerr << "Error: --out is required.\n";
		usage();
		return 2;
	}

	if (opt.hi_mode != "global" && opt.hi_mode != "excl-focus") {
		std::cerr << "Error: --hi-mode must be global or excl-focus\n";
		return 2;
	}

	if (opt.threads < 1) {
		std::cerr << "Error: --threads must be >= 1\n";
		return 2;
	}

	if (!opt.cov_path.empty() && opt.hi_mode == "excl-focus") {
		std::cerr << "Error: --hi-mode excl-focus is incompatible with --cov (LOCO requires internally computed per-chromosome HI)\n";
		return 2;
	}

	if (opt.phased && !opt.intra) {
		std::cerr << "Error: --phased requires --intra (phased LD is only meaningful for intrachromosomal pairs)\n";
		return 2;
	}

	if (opt.phased && opt.hi_mode == "excl-focus") {
		std::cerr << "Error: --phased is incompatible with --hi-mode excl-focus\n";
		return 2;
	}

	if (opt.phased && !opt.sample_geno_path.empty()) {
		std::cerr << "Error: --phased is incompatible with --sample-geno\n";
		return 2;
	}

	if (!opt.cov_path.empty()) {
		if (opt.compute_hi_only)
			std::cerr << "Warning: --compute-hi is ignored when --cov is provided\n";
		if (opt.unweighted_hi)
			std::cerr << "Warning: --unweighted-hi is ignored when --cov is provided\n";
		if (opt.pos_is_start)
			std::cerr << "Warning: --pos-is-start is ignored when --cov is provided\n";
	}

	if (opt.compute_hi_only && opt.hi_mode != "global") {
		std::cerr << "Error: --compute-hi requires --hi-mode global (excl-focus is only relevant during scanning)\n";
		return 2;
	}

	if ((opt.has_min_neg_r && opt.min_neg_r < 0.0f) || (opt.has_min_pos_r && opt.min_pos_r < 0.0f)) {
		std::cerr << "Error: --min-neg-r and --min-pos-r must be >= 0\n";
		return 2;
	}

	if (opt.has_target) {
		if (opt.target_chr.empty() || opt.target_pos < 0) {
			std::cerr << "Error: --target-chr and --target-pos must both be provided.\n";
			return 2;
		}
	}

	if (opt.min_callrate < 0.0 || opt.min_callrate > 1.0) {
		std::cerr << "Error: --min-callrate must be in [0,1]\n";
		return 2;
	}

	return 0;
}

static bool read_bed_if_needed(
	const CliOptions& opt,
	std::unordered_map<std::string, std::vector<BedInterval>>& bed_by_chr,
	bool& use_bed
) {
	use_bed = false;
	if (opt.bed_path.empty())
		return true;

	use_bed = true;
	if (!read_bed(opt.bed_path, bed_by_chr))
		return false;

	std::cout << "Using BED filter: " << opt.bed_path << "\n";
	return true;
}

static void build_chr_maps(
	const CliOptions& opt,
	std::unordered_map<std::string, bool>& chr_keep_map,
	std::unordered_map<std::string, bool>& chr_drop_map
) {
	chr_keep_map.clear();
	chr_drop_map.clear();

	if (!opt.keep_chrs.empty()) {
		chr_keep_map.reserve(opt.keep_chrs.size() * 2 + 1);
		for (const auto& c : opt.keep_chrs)
			chr_keep_map[c] = true;
	}

	if (!opt.drop_chrs.empty()) {
		chr_drop_map.reserve(opt.drop_chrs.size() * 2 + 1);
		for (const auto& c : opt.drop_chrs)
			chr_drop_map[c] = true;
	}
}

static bool apply_block_filters(
	Eigen::MatrixXf& X,
	std::vector<std::string>& chroms,
	std::vector<int>& pos,
	std::vector<int>& pos_start,
	const std::vector<std::string>& sample_names,
	const CliOptions& opt
) {
	const int nrows = (int)X.rows();	// may be 2*nsamples in phased mode
	int nwin = (int)chroms.size();

	// Filters are applied before residualization so downstream matrices match the scanned block set.
	std::unordered_map<std::string, bool> chr_keep_map;
	std::unordered_map<std::string, bool> chr_drop_map;
	build_chr_maps(opt, chr_keep_map, chr_drop_map);

	std::unordered_map<std::string, std::vector<BedInterval>> bed_by_chr;
	bool use_bed = false;
	if (!read_bed_if_needed(opt, bed_by_chr, use_bed))
		return false;

	std::vector<int> keep_idx;
	keep_idx.reserve((size_t)nwin);

	for (int w = 0; w < nwin; ++w) {
		const std::string& chr = chroms[w];
		int p = pos[w];

		if (!opt.keep_chrs.empty()) {
			if (chr_keep_map.find(chr) == chr_keep_map.end())
				continue;
		}

		if (!opt.drop_chrs.empty()) {
			if (chr_drop_map.find(chr) != chr_drop_map.end())
				continue;
		}

		if (use_bed) {
			if (!bed_contains(bed_by_chr, chr, p))
				continue;
		}

		keep_idx.push_back(w);
	}

	if (!opt.drop_chrs.empty()) {
		std::cout << "Excluded chromosomes:\n";
		for (const auto& c : opt.drop_chrs)
			std::cout << "  - " << c << "\n";
	}

	if (keep_idx.empty()) {
		std::cerr << "Error: no blocks remain after filtering.\n";
		return false;
	}

	if ((int)keep_idx.size() == nwin)
		return true;

	std::cout << "Filtering blocks:\n";
	std::cout << "  before = " << nwin << "\n";
	std::cout << "  after  = " << (int)keep_idx.size() << "\n";

	int nwin2 = (int)keep_idx.size();
	Eigen::MatrixXf Xf(nrows, nwin2);
	std::vector<std::string> chroms_f((size_t)nwin2);
	std::vector<int> pos_f((size_t)nwin2);

	for (int j = 0; j < nwin2; ++j) {
		int w = keep_idx[j];
		Xf.col(j) = X.col(w);
		chroms_f[j] = chroms[w];
		pos_f[j] = pos[w];
	}

	X.swap(Xf);
	chroms.swap(chroms_f);
	pos.swap(pos_f);

	if (!pos_start.empty()) {
		std::vector<int> pos_start_f((size_t)nwin2);
		for (int j = 0; j < nwin2; ++j)
			pos_start_f[j] = pos_start[(size_t)keep_idx[j]];
		pos_start.swap(pos_start_f);
	}

	return true;
}

static bool apply_callrate_filter(
	Eigen::MatrixXf& X,
	std::vector<std::string>& chroms,
	std::vector<int>& pos,
	std::vector<int>& pos_start,
	const std::vector<std::string>& sample_names,
	double min_callrate
) {
	const int nrows = (int)X.rows();	// may be 2*nsamples in phased mode
	int nblocks = (int)chroms.size();

	if (nblocks == 0)
		return true;

	// Default is 1.0, which matches current strict behavior.
	// If set < 1.0, we keep blocks called in at least that fraction of rows.
	if (min_callrate >= 1.0)
		return true;

	const int min_called = (int)std::ceil(min_callrate * (double)nrows);

	std::vector<int> keep_idx;
	keep_idx.reserve((size_t)nblocks);

	for (int w = 0; w < nblocks; ++w) {
		int called = 0;
		for (int i = 0; i < nrows; ++i) {
			float v = X(i, w);
			if (!std::isnan(v))
				++called;
		}

		if (called >= min_called)
			keep_idx.push_back(w);
	}

	if (keep_idx.empty()) {
		std::cerr << "Error: no blocks remain after --min-callrate filtering.\n";
		return false;
	}

	if ((int)keep_idx.size() == nblocks) {
		std::cout << "Call-rate filter: kept all blocks (min_callrate=" << min_callrate << ")\n";
		return true;
	}

	std::cout << "Call-rate filter (--min-callrate):\n";
	std::cout << "  min_callrate = " << min_callrate << "\n";
	std::cout << "  min_called   = " << min_called << " / " << nrows << "\n";
	std::cout << "  before       = " << nblocks << "\n";
	std::cout << "  after        = " << (int)keep_idx.size() << "\n";

	const int nblocks2 = (int)keep_idx.size();
	Eigen::MatrixXf Xf(nrows, nblocks2);
	std::vector<std::string> chroms_f((size_t)nblocks2);
	std::vector<int> pos_f((size_t)nblocks2);

	for (int j = 0; j < nblocks2; ++j) {
		int w = keep_idx[j];
		Xf.col(j) = X.col(w);
		chroms_f[j] = chroms[w];
		pos_f[j] = pos[w];
	}

	X.swap(Xf);
	chroms.swap(chroms_f);
	pos.swap(pos_f);

	if (!pos_start.empty()) {
		std::vector<int> pos_start_f((size_t)nblocks2);
		for (int j = 0; j < nblocks2; ++j)
			pos_start_f[j] = pos_start[(size_t)keep_idx[j]];
		pos_start.swap(pos_start_f);
	}

	return true;
}

static void mean_impute_missing_per_block(Eigen::MatrixXf& X) {
	const int nsamples = (int)X.rows();
	const int nblocks = (int)X.cols();

	for (int w = 0; w < nblocks; ++w) {
		double sum = 0.0;
		int cnt = 0;

		for (int i = 0; i < nsamples; ++i) {
			float v = X(i, w);
			if (!std::isnan(v)) {
				sum += (double)v;
				++cnt;
			}
		}

		if (cnt == 0)
			continue;

		float mu = (float)(sum / (double)cnt);

		for (int i = 0; i < nsamples; ++i) {
			if (std::isnan(X(i, w)))
				X(i, w) = mu;
		}
	}
}

static void build_blocks_by_chr_sorted(
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	std::unordered_map<std::string, std::vector<int>>& blocks_by_chr,
	std::vector<std::string>& chr_order
) {
	chr_order.clear();
	blocks_by_chr = group_by_chr(chroms, chr_order);

	for (auto& kv : blocks_by_chr) {
		auto& idx = kv.second;
		std::sort(idx.begin(), idx.end(),
			[&](int a, int b) { return pos[a] < pos[b]; }
		);
	}
}

static bool compute_hi(
	Eigen::VectorXf& h,
	const Eigen::MatrixXf& X,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::vector<int>& pos_start,
	const Eigen::MatrixXf& X_full,
	const std::vector<std::string>& chroms_full,
	const std::vector<int>& pos_full,
	const std::vector<int>& pos_start_full,
	const CliOptions& opt
) {
	// Default: compute HI on full (unfiltered) set for stable ancestry covariate.
	// --compute-hi: compute HI on the filtered block set, then exit (utility mode).
	const Eigen::MatrixXf& X_hi = opt.compute_hi_only ? X : X_full;
	const std::vector<std::string>& chroms_hi = opt.compute_hi_only ? chroms : chroms_full;
	const std::vector<int>& pos_hi = opt.compute_hi_only ? pos : pos_full;
	const std::vector<int>& pos_start_hi = opt.compute_hi_only ? pos_start : pos_start_full;

	if (opt.compute_hi_only)
		std::cout << "Computing HI from FILTERED blocks (--chr/--no-chr/--bed applied)\n";
	else
		std::cout << "Computing HI from full block set (independent of --chr/--no-chr/--bed)\n";

	const int n_diploid = opt.phased ? (int)X_hi.rows() / 2 : (int)X_hi.rows();
	if (opt.unweighted_hi) {
		h = compute_hi_from_X(X_hi, opt.phased, n_diploid);
	} else {
		h = compute_hi_from_X_weighted(
			X_hi,
			chroms_hi,
			pos_hi,
			opt.pos_is_start,
			pos_start_hi,
			opt.phased,
			n_diploid
		);
	}

	return true;
}

static bool write_hi_and_print_summary(
	const std::string& hi_out_path,
	const WindowMatrix& wm,
	const Eigen::VectorXf& h,
	const CliOptions& opt
) {
	std::cout << "HI mode: "
		<< (opt.unweighted_hi ? "unweighted" : "weighted")
		<< (opt.unweighted_hi ? "" : (opt.pos_is_start ? " (pos=start)" : " (pos=end)"))
		<< "\n";

	if (!write_hi_tsv(hi_out_path, wm.sample_names, h))
		return false;

	std::cout << "Wrote HI to: " << hi_out_path << "\n";

	const int nsamples = (int)wm.sample_names.size();
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

	std::cout << "Hybrid index (from blocks):\n";
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
	return true;
}

static bool build_hi_components_if_needed(
	HiComponentsWeighted& hc_full,
	bool& has_hc_full,
	const Eigen::MatrixXf& X_full,
	const std::vector<std::string>& chroms_full,
	const std::vector<int>& pos_full,
	const std::vector<int>& pos_start_full,
	const CliOptions& opt
) {
	has_hc_full = false;

	if (opt.hi_mode != "excl-focus")
		return true;

	// excl-focus requires weighted components so we can subtract focus-chrom contributions on the fly.
	if (opt.unweighted_hi) {
		std::cerr << "Error: --hi-mode excl-focus currently requires weighted HI (omit --unweighted-hi)\n";
		return false;
	}

	hc_full = build_hi_components_weighted(
		X_full,
		chroms_full,
		pos_full,
		opt.pos_is_start,
		pos_start_full
	);
	has_hc_full = true;
	return true;
}

static int resolve_target_block(
	const CliOptions& opt,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos
) {
	if (!opt.has_target)
		return -1;

	for (int w = 0; w < (int)chroms.size(); ++w) {
		if (chroms[w] == opt.target_chr && pos[w] == opt.target_pos)
			return w;
	}

	return -2;
}

static bool residualize_blocks(
	Eigen::MatrixXf& Z,
	const Eigen::MatrixXf& X,
	const Eigen::MatrixXf& H,
	int& n_valid_blocks
) {
	try {
		Z = residualize_and_zscore(X, H, n_valid_blocks);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return false;
	}
	return true;
}

static void print_residualization_sanity(const Eigen::MatrixXf& Z) {
	const int nrows = (int)Z.rows();
	auto is_block_valid = [&](int w) -> bool {
		for (int i = 0; i < nrows; ++i) {
			if (std::isnan(Z(i, w))) return false;
		}
		return true;
	};

	const int nblocks = (int)Z.cols();
	if (nblocks >= 2 && is_block_valid(0) && is_block_valid(1)) {
		double dot = 0.0;
		for (int i = 0; i < nrows; ++i)
			dot += (double)Z(i, 0) * (double)Z(i, 1);
		double r = dot / (nrows - 1);
		std::cout << "  sanity r(Z_b0, Z_b1) = " << r << "\n";
	} else {
		std::cout << "  sanity r(Z_b0, Z_b1) skipped (invalid block or <2 blocks)\n";
	}
}

static void write_perm_summary_tsv(
	const std::string& path,
	const std::vector<PermSummary>& summ
) {
	std::ofstream pf(path);
	pf << "rep\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\tmean_r2\tsd_r2\n";

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
			<< "\t" << summ[r].mean
			<< "\t" << summ[r].sd
			<< "\t" << summ[r].mean_r2
			<< "\t" << summ[r].sd_r2
			<< "\n";
	}
}

static bool apply_sample_filter(WindowMatrix& wm, const std::string& keep_indv_path) {
	std::ifstream f(keep_indv_path);
	if (!f) {
		std::cerr << "Error: cannot open --keep-indv file: " << keep_indv_path << "\n";
		return false;
	}

	std::unordered_set<std::string> keep_set;
	std::string line;
	while (std::getline(f, line)) {
		while (!line.empty() && (line.back() == '\r' || line.back() == ' ' || line.back() == '\t'))
			line.pop_back();
		if (!line.empty())
			keep_set.insert(line);
	}

	if (keep_set.empty()) {
		std::cerr << "Error: --keep-indv file has no valid IDs.\n";
		return false;
	}

	const int n_dip = wm.nsamples_diploid;
	std::vector<int> keep_idx;
	keep_idx.reserve((size_t)n_dip);

	for (int i = 0; i < n_dip; ++i) {
		if (keep_set.count(wm.sample_names[i]))
			keep_idx.push_back(i);
	}

	if (keep_idx.empty()) {
		std::cerr << "Error: no samples remain after --keep-indv filtering.\n";
		return false;
	}

	const int n_keep = (int)keep_idx.size();

	std::cout << "Sample filter (--keep-indv):\n";
	std::cout << "  requested = " << keep_set.size() << "\n";
	std::cout << "  matched   = " << n_keep << " / " << n_dip << "\n";

	if (n_keep == n_dip)
		return true;

	const int nwin = (int)wm.X.cols();

	if (!wm.phased) {
		Eigen::MatrixXf Xf(n_keep, nwin);
		for (int i = 0; i < n_keep; ++i)
			Xf.row(i) = wm.X.row(keep_idx[i]);
		wm.X = std::move(Xf);
	} else {
		Eigen::MatrixXf Xf(2 * n_keep, nwin);
		for (int i = 0; i < n_keep; ++i) {
			Xf.row(2 * i)     = wm.X.row(2 * keep_idx[i]);
			Xf.row(2 * i + 1) = wm.X.row(2 * keep_idx[i] + 1);
		}
		wm.X = std::move(Xf);
	}

	std::vector<std::string> names_f;
	names_f.reserve((size_t)n_keep);
	for (int i : keep_idx)
		names_f.push_back(wm.sample_names[i]);
	wm.sample_names = std::move(names_f);
	wm.nsamples_diploid = n_keep;

	return true;
}

static bool make_asym_flags(
	const CliOptions& opt,
	bool& use_asym,
	float& min_neg_r,
	float& min_pos_r
) {
	use_asym = false;
	min_neg_r = opt.min_neg_r;
	min_pos_r = opt.min_pos_r;

	if (opt.has_min_neg_r || opt.has_min_pos_r) {
		use_asym = true;

		if (!opt.has_min_neg_r)
			min_neg_r = 0.0f;
		if (!opt.has_min_pos_r)
			min_pos_r = 0.0f;

		if (min_neg_r < 0.0f || min_pos_r < 0.0f) {
			std::cerr << "Error: --min-neg-r and --min-pos-r must be >= 0\n";
			return false;
		}
	}

	return true;
}

int main(int argc, char** argv) {
	// Stage 1: Parse + validate CLI
	CliOptions cli;
	if (!parse_args(argc, argv, cli))
		return 2;

	if (int rc = validate_options(cli); rc != 0)
		return rc;

	const bool sample_geno_mode = !cli.sample_geno_path.empty();
	if (sample_geno_mode && cli.intra)
		std::cout << "Note: --intra ignored in --sample-geno mode\n";
	if (sample_geno_mode && cli.has_target)
		std::cout << "Note: --target-* ignored in --sample-geno mode\n";

	// Stage 2: Load blocks from VCF or MSP
	WindowMatrix wm;
	{
		const int max_win = cli.max_windows;
		if (!cli.vcf_path.empty()) {
			VcfLoadOptions vopt;
			vopt.max_windows = max_win;
			vopt.phased = cli.phased;
			wm = load_windows_from_vcf(cli.vcf_path, vopt);
		} else {
			if (cli.pos_is_start)
				std::cerr << "Warning: --pos-is-start has no effect with --msp (exact block lengths are used).\n";
			MspLoadOptions mopt;
			mopt.max_windows = max_win;
			mopt.phased = cli.phased;
			wm = load_windows_from_msp(cli.msp_path, mopt);
		}
	}

	// Stage 2b: Apply optional sample filter (before X/X_full split)
	if (!cli.keep_indv_path.empty()) {
		if (!apply_sample_filter(wm, cli.keep_indv_path))
			return 1;
	}

	const int nsamples_diploid = wm.nsamples_diploid;
	const int nsamples = (int)wm.sample_names.size();	// diploid sample count for HI output / reporting

	Eigen::MatrixXf X = wm.X;
	Eigen::MatrixXf X_full = X;

	std::vector<std::string> chroms = wm.meta.chrom;
	std::vector<int> pos = wm.meta.pos;
	std::vector<int> pos_start = wm.meta.pos_start;

	// Copies for full-genome HI + LOCO components (must be defined BEFORE any filtering).
	std::vector<std::string> chroms_full = chroms;
	std::vector<int> pos_full = pos;
	std::vector<int> pos_start_full = pos_start;

	// Stage 3a: Apply optional block filters
	if (!apply_block_filters(X, chroms, pos, pos_start, wm.sample_names, cli))
		return 1;

	// Stage 3b: Apply optional callrate filters
	if (!apply_callrate_filter(X, chroms, pos, pos_start, wm.sample_names, cli.min_callrate))
		return 1;

	if (cli.min_callrate < 1.0) {
		std::cout << "Missing-data handling: mean-imputing within blocks (enabled by --min-callrate < 1.0)\n";
		mean_impute_missing_per_block(X);
	}

	const int nblocks = (int)chroms.size();

	std::cout << "admixld input OK\n";
	if (!cli.vcf_path.empty())
		std::cout << "  vcf        = " << cli.vcf_path << "\n";
	else
		std::cout << "  msp        = " << cli.msp_path << "\n";
	std::cout << "  out        = " << cli.out << "\n";
	std::cout << "  nsamples   = " << nsamples_diploid
		<< (cli.phased ? " (phased: " + std::to_string(2 * nsamples_diploid) + " haplotype rows)" : "") << "\n";
	std::cout << "  blocks     = " << nblocks << "\n";
	std::cout << "HI correction mode: " << cli.hi_mode << "\n";

	int to_print = std::min(5, nsamples);
	for (int i = 0; i < to_print; ++i)
		std::cout << "  sample[" << i << "] = " << wm.sample_names[i] << "\n";
	std::cout << std::flush;

	// Stage 4: Build chromosome index for blocks
	std::unordered_map<std::string, std::vector<int>> blocks_by_chr;
	std::vector<std::string> chr_order;
	build_blocks_by_chr_sorted(chroms, pos, blocks_by_chr, chr_order);
	std::cout << "Chromosomes in loaded blocks: " << chr_order.size() << "\n";

	// Stage 5: Compute HI from blocks or load covariate matrix from --cov
	Eigen::MatrixXf H;

	if (!cli.cov_path.empty()) {
		if (!load_cov_tsv(cli.cov_path, wm.sample_names, H))
			return 1;
	} else {
		Eigen::VectorXf h;
		if (!compute_hi(h, X, chroms, pos, pos_start, X_full, chroms_full, pos_full, pos_start_full, cli))
			return 1;

		std::string hi_out_path = cli.out + ".hi.tsv";
		if (!write_hi_and_print_summary(hi_out_path, wm, h, cli))
			return 1;

		if (cli.compute_hi_only) {
			std::cout << "--compute-hi set: exiting after HI computation (no scans)\n";
			return 0;
		}

		H.resize(h.size(), 1);
		H.col(0) = h;
	}

	// In phased mode, duplicate H rows so H.rows() == X.rows() == 2*nsamples_diploid
	if (cli.phased) {
		const int k = (int)H.cols();
		Eigen::MatrixXf H2(2 * nsamples_diploid, k);
		for (int i = 0; i < nsamples_diploid; ++i) {
			H2.row(2 * i)     = H.row(i);
			H2.row(2 * i + 1) = H.row(i);
		}
		H = std::move(H2);
	}

	// Stage 6: Build HI components for excl-focus mode (full genome)
	HiComponentsWeighted hc_full;
	bool has_hc_full = false;
	if (!build_hi_components_if_needed(hc_full, has_hc_full, X_full, chroms_full, pos_full, pos_start_full, cli))
		return 2;

	// Stage 7: Resolve target block (after filtering)
	int target_w = -1;
	if (cli.has_target) {
		int tw = resolve_target_block(cli, chroms, pos);
		if (tw == -2) {
			std::cerr << "Error: target block not found after filtering: "
				<< cli.target_chr << ":" << cli.target_pos << "\n";
			return 1;
		}
		target_w = tw;
		std::cout << "Target block: w=" << target_w
			<< " chr=" << cli.target_chr
			<< " pos=" << cli.target_pos << "\n";
	}

	// Stage 8: Sample-geno mode (trait/mito vs all blocks)
	Eigen::VectorXf gZ;
	int g_valid = 0;

	Eigen::VectorXf g;
	bool g_loaded = false;

	if (sample_geno_mode) {
		if (!load_sample_vector_tsv(cli.sample_geno_path, wm.sample_names, g))
			return 1;
		g_loaded = true;

		if (cli.hi_mode == "global") {
			try {
				gZ = residualize_and_zscore_vector(g, H, g_valid);
			} catch (const std::exception& e) {
				std::cerr << "Error: sample-geno residualization failed: "
					<< e.what() << "\n";
				return 1;
			}

			std::cout << "Sample-geno residualization complete:\n";
			std::cout << "  valid_samples = " << g_valid << " / " << nsamples << "\n";
		} else {
			std::cout << "Sample-geno loaded (excl-focus mode: residualize per chromosome during scan)\n";
		}
	}

	// Stage 9: Residualize blocks on HI and z-score => Z (nsamples × nblocks)
	Eigen::MatrixXf Z;
	int n_valid_blocks = 0;

	if (!residualize_blocks(Z, X, H, n_valid_blocks))
		return 1;

	std::cout << "Residualization complete:\n";
	std::cout << "  valid_blocks = " << n_valid_blocks << " / " << nblocks << "\n";
	print_residualization_sanity(Z);
	std::cout << std::flush;

	// Stage 10: Scan options for permutations (symmetric only)
	ScanOptions popt;
	popt.intra = cli.intra;
	popt.max_dist = cli.max_dist;
	popt.tile_size = cli.tile_size;
	popt.nsamples = (int)Z.rows();
	popt.threads = cli.threads;
	popt.min_abs_r = (float)cli.min_abs_r;
	popt.use_asym = false;
	popt.min_neg_r = 0.0f;
	popt.min_pos_r = 0.0f;

	// Stage 10a: Sample-geno permutation test
	if (sample_geno_mode && cli.n_perm > 0) {
		std::vector<PermSummary> summ;

		std::cout << "Permutation test (sample-geno vs blocks):\n";
		std::cout << "  n_perm      = " << cli.n_perm << "\n";
		std::cout << "  seed        = " << cli.seed << "\n";
		std::cout << "  sample_size = " << cli.perm_sample << "\n";

		if (cli.hi_mode == "global") {
			if (!permute_sample_vector_summary(
				Z, gZ,
				blocks_by_chr, chr_order,
				popt,
				cli.seed, cli.n_perm, cli.perm_sample,
				summ
			))
				return 1;
		} else {
			if (!g_loaded) {
				std::cerr << "Error: sample-geno vector was not loaded\n";
				return 1;
			}
			if (!has_hc_full) {
				std::cerr << "Error: HI components missing for excl-focus mode\n";
				return 1;
			}
			if (!permute_sample_vector_summary_excl_focus(
				X,
				g,
				chroms,
				pos,
				blocks_by_chr,
				chr_order,
				hc_full,
				popt,
				cli.seed,
				cli.n_perm,
				cli.perm_sample,
				summ
			))
				return 1;
		}

		std::string perm_path = cli.out + ".samplegeno.perm.summary.tsv";
		write_perm_summary_tsv(perm_path, summ);
		std::cout << "  wrote permutation summaries: " << perm_path << "\n";
		return 0;
	}

	// Stage 10b: Block/block permutation test (interchrom chr-block)
	if (cli.n_perm > 0) {
		if (cli.intra) {
			std::cerr << "Error: --permute is currently implemented for interchrom scans only (omit --intra).\n";
			return 1;
		}

		std::vector<PermSummary> summ;

		std::cout << "Permutation test (interchrom, chr-block):\n";
		std::cout << "  n_perm      = " << cli.n_perm << "\n";
		std::cout << "  seed        = " << cli.seed << "\n";
		std::cout << "  perm_sample = " << cli.perm_sample << "\n";
		std::cout << std::flush;

		if (cli.hi_mode == "global") {
			if (!permute_interchrom_summary_chrblock(
				Z,
				blocks_by_chr, chr_order,
				popt,
				cli.seed, cli.n_perm, cli.perm_sample,
				summ
			))
				return 1;
		} else {
			if (!has_hc_full) {
				std::cerr << "Error: HI components missing for excl-focus mode\n";
				return 1;
			}
			if (!permute_interchrom_summary_chrblock_excl_focus(
				X,
				chroms,
				pos,
				blocks_by_chr,
				chr_order,
				hc_full,
				popt,
				cli.seed,
				cli.n_perm,
				cli.perm_sample,
				summ
			))
				return 1;
		}

		std::string perm_path = cli.out + ".perm.summary.tsv";
		write_perm_summary_tsv(perm_path, summ);
		std::cout << "  wrote permutation summaries: " << perm_path << "\n";
		std::cout << "Permutation mode (--permute was set). Skipping hits scan.\n";
		return 0;
	}

	// Stage 11: Final scan options (including asymmetric filtering)
	bool use_asym = false;
	float min_neg_r = 0.0f;
	float min_pos_r = 0.0f;

	if (!make_asym_flags(cli, use_asym, min_neg_r, min_pos_r))
		return 2;

	ScanOptions opt;
	opt.intra = cli.intra;
	opt.max_dist = cli.max_dist;
	opt.tile_size = cli.tile_size;
	opt.nsamples = (int)Z.rows();
	opt.threads = cli.threads;
	opt.min_abs_r = (float)cli.min_abs_r;
	opt.use_asym = use_asym;
	opt.min_neg_r = min_neg_r;
	opt.min_pos_r = min_pos_r;

	// Stage 11a: Sample-geno scan overrides block/block scan
	if (sample_geno_mode) {
		std::string out_path = cli.out + ".samplegeno.hits.tsv";
		std::string distrib_path;
		if (cli.distrib)
			distrib_path = cli.out + ".samplegeno.scan.summary.tsv";

		long long tested = 0;
		long long kept = 0;

		std::cout << "Scan mode: sample-geno vs all blocks\n";

		if (cli.hi_mode == "global") {
			if (!scan_vector_vs_windows_write_hits(
				Z, gZ,
				chroms, pos,
				blocks_by_chr, chr_order,
				opt,
				out_path,
				tested,
				kept,
				distrib_path,
				cli.distrib_sample,
				cli.seed
			))
				return 1;

		} else {
			if (!g_loaded) {
				std::cerr << "Error: sample-geno vector was not loaded\n";
				return 1;
			}
			if (!has_hc_full) {
				std::cerr << "Error: HI components missing for excl-focus mode\n";
				return 1;
			}

			if (!scan_vector_vs_windows_write_hits_excl_focus(
				X,
				g,
				chroms,
				pos,
				blocks_by_chr,
				chr_order,
				hc_full,
				opt,
				out_path,
				tested,
				kept,
				distrib_path,
				cli.distrib_sample,
				cli.seed
			))
				return 1;
		}

		std::cout << "Sample-geno scan complete:\n";
		std::cout << "  tested_pairs = " << tested << "\n";
		std::cout << "  kept_pairs   = " << kept << "\n";
		std::cout << "  wrote        = " << out_path << "\n";
		if (cli.distrib)
			std::cout << "  wrote        = " << distrib_path << "\n";

		return 0;
	}

	// Stage 11b: Block/block scan
	std::string out_path = cli.out + ".hits.tsv";
	std::string distrib_path;
	if (cli.distrib)
		distrib_path = cli.out + ".scan.summary.tsv";

	long long tested = 0;
	long long kept = 0;

	std::cout << "Scan mode: " << (cli.intra ? "intrachromosomal" : "interchromosomal") << "\n";
	std::cout << "Block size: " << cli.tile_size << "\n";

	if (cli.hi_mode == "global") {
		if (target_w >= 0) {
			if (!scan_target_write_hits(
				Z, chroms, pos, blocks_by_chr, chr_order,
				opt, target_w, out_path, tested, kept,
				distrib_path, cli.distrib_sample, cli.seed
			))
				return 1;
		} else {
			if (!scan_blocks_write_hits(
				Z, chroms, pos, blocks_by_chr, chr_order,
				opt, out_path, tested, kept,
				distrib_path, cli.distrib_sample, cli.seed
			))
				return 1;
		}
	} else {
		if (target_w >= 0) {
			if (!scan_target_write_hits_excl_focus(
				X,
				chroms,
				pos,
				blocks_by_chr,
				chr_order,
				hc_full,
				opt,
				target_w,
				out_path,
				tested,
				kept,
				distrib_path,
				cli.distrib_sample,
				cli.seed
			))
				return 1;
		} else {
			if (!scan_blocks_write_hits_excl_focus(
				X,
				chroms,
				pos,
				blocks_by_chr,
				chr_order,
				hc_full,
				opt,
				out_path,
				tested,
				kept,
				distrib_path,
				cli.distrib_sample,
				cli.seed
			))
				return 1;
		}
	}

	std::cout << "Scan complete:\n";
	std::cout << "  tested_pairs = " << tested << "\n";
	std::cout << "  kept_pairs   = " << kept << " (|r| >= " << cli.min_abs_r << ")\n";
	std::cout << "  wrote        = " << out_path << "\n";
	if (cli.distrib)
		std::cout << "  wrote        = " << distrib_path << "\n";

	if (use_asym) {
		std::cout << "  r filter    = ";
		if (cli.has_min_neg_r)
			std::cout << "r <= -" << min_neg_r << " ";
		if (cli.has_min_pos_r)
			std::cout << "r >= " << min_pos_r;
		std::cout << "\n";
	} else {
		std::cout << "  |r| filter  = " << cli.min_abs_r << "\n";
	}

	return 0;
}
