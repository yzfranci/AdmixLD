#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "../io/ref_freq.hpp"

// Compute hybrid index from X as mean(dosage)/2 per sample, ignoring NaNs.
// Diploid mode (phased=false): X is nsamples x nwin, values in [0,2].
// Phased mode (phased=true):   X is 2*nsamples x nwin, values in {0,1};
//   rows 2i and 2i+1 are the two haplotypes of sample i.
//   nsamples_diploid must equal nsamples when phased=true.
Eigen::VectorXf compute_hi_from_X(
	const Eigen::MatrixXf& X,
	bool phased = false,
	int nsamples_diploid = 0
);

// Weighted HI using tract lengths per chromosome.
// If pos_start is non-empty (MSP input): exact lengths computed as pos[w] - pos_start[w].
// Otherwise (VCF input): lengths inferred from pos alone:
//   pos_is_start=false (default): pos is END, len[k] = pos[k] - pos[k-1]
//   pos_is_start=true:            pos is START, len[k] = pos[k+1] - pos[k]
// Uses weights only for non-missing genotypes in X.
// phased/nsamples_diploid: same semantics as compute_hi_from_X.
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

// Build per-chromosome weighted components (for LOCO/LOCO2).
// Uses the same tract-length logic as compute_hi_from_X_weighted.
HiComponentsWeighted build_hi_components_weighted(
	const Eigen::MatrixXf& X,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	bool pos_is_start,
	const std::vector<int>& pos_start = {}
);

// Compute weighted HI from components excluding none / one / two chromosomes.
// If exclude_chrA is empty => exclude none.
// If exclude_chrB is empty => exclude only A (if A provided).
Eigen::VectorXf hi_from_components_weighted_excluding(
	const HiComponentsWeighted& hc,
	const std::string& exclude_chrA = "",
	const std::string& exclude_chrB = ""
);

// Struct for frequency-based HI components (LOCO with --ref-freq).
// Accumulates delta^2-weighted numerator and denominator per chromosome.
// num = Σ((d_iw/2 - p2_w)*delta_w*L_w), den = Σ(delta_w^2*L_w) for non-NaN d_iw.
struct HiComponentsFreq {
	Eigen::VectorXf num_total;	// per sample
	Eigen::VectorXf den_total;	// per sample

	std::vector<std::string> chr_order;		// stable insertion order
	std::vector<Eigen::VectorXf> num_chr;	// [c] per sample
	std::vector<Eigen::VectorXf> den_chr;	// [c] per sample
};

// Compute frequency-based HI: HI_i = Σ((d/2-p2)*delta*L) / Σ(delta^2*L).
// freqs must have length == X.cols(); markers with NaN freqs are skipped.
// Uses the same tract-length logic as compute_hi_from_X_weighted.
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
