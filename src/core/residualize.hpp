#pragma once

#include <Eigen/Dense>

// Residualize each window on covariate matrix H (nsamples x k) and z-score residuals.
// Returns Z (nsamples × nwin) with NaNs for invalid windows (any NaN in X for that window).
// H must have no missing values. n_valid_windows is set to the number of valid windows.
Eigen::MatrixXf residualize_and_zscore(
	const Eigen::MatrixXf& X,
	const Eigen::MatrixXf& H,
	int& n_valid_windows
);

Eigen::VectorXf residualize_and_zscore_vector(
	const Eigen::VectorXf& y,
	const Eigen::MatrixXf& H,
	int& n_valid
);

Eigen::MatrixXf residualize_and_zscore_subset(
	const Eigen::MatrixXf& X,
	const Eigen::VectorXf& h,
	const std::vector<int>& win_idx,
	int& n_valid_windows
);
