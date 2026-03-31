#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

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

// Weighted HI using block lengths per chromosome.
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
// Uses the same block-length logic as compute_hi_from_X_weighted.
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

// Load covariate matrix from TSV: sample<TAB>cov1 [<TAB>cov2 ...] (header allowed).
// Number of covariate columns is inferred from the first data row.
// Requires all samples to be present; no missing values allowed.
bool load_cov_tsv(
	const std::string& path,
	const std::vector<std::string>& sample_names,
	Eigen::MatrixXf& H_out
);

// Write hybrid index to TSV: sample<TAB>hi
bool write_hi_tsv(
	const std::string& path,
	const std::vector<std::string>& sample_names,
	const Eigen::VectorXf& h
);
