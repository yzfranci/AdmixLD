#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

// Compute hybrid index from X (dosage 0/1/2) as mean(dosage)/2 per sample, ignoring NaNs
Eigen::VectorXf compute_hi_from_X(
	const Eigen::MatrixXf& X
);

// Weighted HI using inferred block lengths per chromosome.
// If pos_is_start=false (default): pos is END of block
//		len[0] = pos0
//		len[k] = pos[k] - pos[k-1]
// If pos_is_start=true: pos is START of block
//		len[k] = pos[k+1] - pos[k] (except last, uses previous len or 1)
// Uses weights only for non-missing genotypes in X.
Eigen::VectorXf compute_hi_from_X_weighted(
	const Eigen::MatrixXf& X,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	bool pos_is_start
);

struct HiComponentsWeighted {
	Eigen::VectorXf wsum_total;	// per sample: sum(v * L) over non-missing
	Eigen::VectorXf wtot_total;	// per sample: sum(L) over non-missing

	std::vector<std::string> chr_order;		// stable insertion order
	std::vector<Eigen::VectorXf> wsum_chr;	// [c] per sample
	std::vector<Eigen::VectorXf> wtot_chr;	// [c] per sample
};

// Build per-chromosome weighted components (for LOCO/LOCO2).
// Uses the same block-length inference as compute_hi_from_X_weighted.
HiComponentsWeighted build_hi_components_weighted(
	const Eigen::MatrixXf& X,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	bool pos_is_start
);

// Compute weighted HI from components excluding none / one / two chromosomes.
// If exclude_chrA is empty => exclude none.
// If exclude_chrB is empty => exclude only A (if A provided).
Eigen::VectorXf hi_from_components_weighted_excluding(
	const HiComponentsWeighted& hc,
	const std::string& exclude_chrA = "",
	const std::string& exclude_chrB = ""
);

// Load hybrid index from TSV: sample<TAB>hi (header allowed). Requires all samples.
bool load_hi_tsv(
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
