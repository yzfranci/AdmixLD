#pragma once

#include <Eigen/Dense>
#include <vector>
#include "../io/ref_freq.hpp"

// Plug-in residualization of marker matrix X on hybrid index h, then z-score.
// Without freqs: r(r,w) = X(r,w) - dos_scale * h(sample_idx)
// With freqs:    r(r,w) = X(r,w) - dos_scale * (p2_w + delta_w * h(sample_idx))
// phased: dos_scale=1, sample_idx = row/2; diploid: dos_scale=2, sample_idx = row.
// h has nsamples_diploid entries regardless of phased mode.
Eigen::MatrixXf residualize_and_zscore(
	const Eigen::MatrixXf& X,
	const Eigen::VectorXf& h,
	bool phased,
	int& n_valid_windows,
	const std::vector<MarkerFreq>* freqs = nullptr
);

// Plug-in residualization of a haploid 0/1 vector on HI: r_i = y_i - h_i.
// NaN entries in y are propagated. Values must be 0 or 1 (validated at load time).
Eigen::VectorXf residualize_and_zscore_vector(
	const Eigen::VectorXf& y,
	const Eigen::VectorXf& h,
	int& n_valid
);

// Per-chromosome plug-in residualization for excl-focus mode (phased not needed: blocked).
// win_idx: column indices into X. freqs: optional per-marker frequencies indexed by column.
Eigen::MatrixXf residualize_and_zscore_subset(
	const Eigen::MatrixXf& X,
	const Eigen::VectorXf& h,
	const std::vector<int>& win_idx,
	int& n_valid_windows,
	const std::vector<MarkerFreq>* freqs = nullptr
);
