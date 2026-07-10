#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "../io/ref_freq.hpp"

// mean(dosage)/2 per sample, ignoring NaNs.
// Phased mode: X is 2*nsamples x nwin (rows 2i/2i+1 = sample i's haplotypes).
Eigen::VectorXf compute_hi_from_X(
	const Eigen::MatrixXf& X,
	bool phased = false,
	int nsamples_diploid = 0
);

// Tract-length-weighted HI. pos_start non-empty (MSP) => exact lengths from
// pos_start/pos; otherwise lengths are inferred from consecutive pos (see
// pos_is_start).
Eigen::VectorXf compute_hi_from_X_weighted(
	const Eigen::MatrixXf& X,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	bool pos_is_start,
	const std::vector<int>& pos_start = {},
	bool phased = false,
	int nsamples_diploid = 0
);

struct HiComponentsWeighted {
	Eigen::VectorXf wsum_total;	// per sample: sum(v * L) over non-missing
	Eigen::VectorXf wtot_total;	// per sample: sum(L) over non-missing

	std::vector<std::string> chr_order;		// stable insertion order
	std::vector<Eigen::VectorXf> wsum_chr;	// [c] per sample
	std::vector<Eigen::VectorXf> wtot_chr;	// [c] per sample
};

// Per-chromosome weighted components (for LOCO/LOCO2).
// unweighted=true forces weight 1 per marker instead of tract length.
HiComponentsWeighted build_hi_components_weighted(
	const Eigen::MatrixXf& X,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	bool pos_is_start,
	const std::vector<int>& pos_start = {},
	bool unweighted = false
);

// Weighted HI excluding none / one / two chromosomes (empty string = none).
Eigen::VectorXf hi_from_components_weighted_excluding(
	const HiComponentsWeighted& hc,
	const std::string& exclude_chrA = "",
	const std::string& exclude_chrB = ""
);

// Per-chromosome frequency-based HI components (LOCO with --ref-freq).
struct HiComponentsFreq {
	Eigen::VectorXf num_total;	// per sample
	Eigen::VectorXf den_total;	// per sample

	std::vector<std::string> chr_order;		// stable insertion order
	std::vector<Eigen::VectorXf> num_chr;	// [c] per sample
	std::vector<Eigen::VectorXf> den_chr;	// [c] per sample
};

// Frequency-based HI. freqs must have length == X.cols(); markers with NaN
// freqs are skipped.
Eigen::VectorXf compute_hi_freq(
	const Eigen::MatrixXf& X,
	const std::vector<MarkerFreq>& freqs,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	bool pos_is_start,
	const std::vector<int>& pos_start = {},
	bool phased = false,
	int nsamples_diploid = 0
);

// Build per-chromosome frequency-based HI components for LOCO.
HiComponentsFreq build_hi_components_freq(
	const Eigen::MatrixXf& X,
	const std::vector<MarkerFreq>& freqs,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	bool pos_is_start,
	const std::vector<int>& pos_start = {}
);

// Compute frequency-based HI from components excluding none / one / two chromosomes.
Eigen::VectorXf hi_from_components_freq_excluding(
	const HiComponentsFreq& hc,
	const std::string& exclude_chrA = "",
	const std::string& exclude_chrB = ""
);

// Load a single-column HI file from TSV: sample<TAB>hi (header allowed).
// Exactly one value column is required; returns false if more columns are found.
// Requires all samples to be present; no missing values allowed.
bool load_cov_tsv(
	const std::string& path,
	const std::vector<std::string>& sample_names,
	Eigen::VectorXf& h_out
);

// Write hybrid index to TSV: sample<TAB>hi
bool write_hi_tsv(
	const std::string& path,
	const std::vector<std::string>& sample_names,
	const Eigen::VectorXf& h
);
