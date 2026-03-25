#include "residualize.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <limits>

Eigen::MatrixXf residualize_and_zscore(
	const Eigen::MatrixXf& X,
	const Eigen::MatrixXf& H,
	int& n_valid_windows
) {
	const int nsamples = (int)X.rows();
	const int nwin     = (int)X.cols();
	const int k        = (int)H.cols();

	Eigen::MatrixXf Z(nsamples, nwin);
	Z.setConstant(std::numeric_limits<float>::quiet_NaN());

	// Center each covariate column (H has no missing values by contract)
	Eigen::MatrixXd Hc(nsamples, k);
	for (int j = 0; j < k; ++j) {
		double col_mean = 0.0;
		for (int i = 0; i < nsamples; ++i)
			col_mean += (double)H(i, j);
		col_mean /= nsamples;
		for (int i = 0; i < nsamples; ++i)
			Hc(i, j) = (double)H(i, j) - col_mean;
	}

	// Precompute (H^T H)^{-1} H^T  (k × nsamples) — reused for every window
	Eigen::MatrixXd HtH = Hc.transpose() * Hc;
	Eigen::LDLT<Eigen::MatrixXd> ldlt(HtH);
	if (ldlt.info() != Eigen::Success || !ldlt.isPositive())
		throw std::runtime_error("Covariate matrix is singular or not positive-definite; cannot residualize.");
	Eigen::MatrixXd HtH_inv_Ht = ldlt.solve(Hc.transpose());  // k × nsamples

	n_valid_windows = 0;

	for (int w = 0; w < nwin; ++w) {
		// Skip window if any genotype is missing
		bool ok = true;
		double xmean = 0.0;
		for (int i = 0; i < nsamples; ++i) {
			float x = X(i, w);
			if (std::isnan(x)) { ok = false; break; }
			xmean += x;
		}
		if (!ok || nsamples < k + 2)
			continue;
		xmean /= nsamples;

		// Center X for this window
		Eigen::VectorXd xc(nsamples);
		for (int i = 0; i < nsamples; ++i)
			xc(i) = (double)X(i, w) - xmean;

		// Regression coefficients and residuals
		Eigen::VectorXd b = HtH_inv_Ht * xc;
		Eigen::VectorXd r = xc - Hc * b;

		// Z-score residuals
		double rmean = r.mean();
		double rss   = (r.array() - rmean).square().sum();
		double sd    = std::sqrt(rss / (nsamples - 1));
		if (!(sd > 0.0) || !std::isfinite(sd))
			continue;

		for (int i = 0; i < nsamples; ++i)
			Z(i, w) = (float)((r(i) - rmean) / sd);

		++n_valid_windows;
	}

	return Z;
}


Eigen::VectorXf residualize_and_zscore_vector(
	const Eigen::VectorXf& y,
	const Eigen::MatrixXf& H,
	int& n_valid
) {
	const int n = (int)y.size();
	const int k = (int)H.cols();

	// Collect valid indices: non-NaN in y (H has no missing values by contract)
	std::vector<int> idx;
	idx.reserve(n);
	for (int i = 0; i < n; ++i)
		if (!std::isnan(y(i))) idx.push_back(i);
	n_valid = (int)idx.size();

	if (n_valid < k + 2)
		throw std::runtime_error("residualize_and_zscore_vector: too few valid samples");

	// Build subsets
	Eigen::VectorXd yv(n_valid);
	Eigen::MatrixXd Hv(n_valid, k);
	for (int ii = 0; ii < n_valid; ++ii) {
		yv(ii) = (double)y(idx[ii]);
		for (int j = 0; j < k; ++j)
			Hv(ii, j) = (double)H(idx[ii], j);
	}

	// Center y and each covariate column using means from valid samples
	double ymean = yv.mean();
	Eigen::VectorXd yc = yv.array() - ymean;

	Eigen::VectorXd col_means(k);
	Eigen::MatrixXd Hvc(n_valid, k);
	for (int j = 0; j < k; ++j) {
		col_means(j) = Hv.col(j).mean();
		Hvc.col(j) = Hv.col(j).array() - col_means(j);
	}

	// Fit b = (Hvc^T Hvc)^{-1} Hvc^T yc
	Eigen::LDLT<Eigen::MatrixXd> ldlt(Hvc.transpose() * Hvc);
	if (ldlt.info() != Eigen::Success || !ldlt.isPositive())
		throw std::runtime_error("residualize_and_zscore_vector: singular covariate matrix");
	Eigen::VectorXd b = ldlt.solve(Hvc.transpose() * yc);

	// Compute residuals for all n samples (NaN where y is NaN)
	Eigen::VectorXf r(n);
	for (int i = 0; i < n; ++i) {
		if (std::isnan(y(i))) {
			r(i) = std::numeric_limits<float>::quiet_NaN();
			continue;
		}
		Eigen::VectorXd hc_i(k);
		for (int j = 0; j < k; ++j)
			hc_i(j) = (double)H(i, j) - col_means(j);
		r(i) = (float)(((double)y(i) - ymean) - hc_i.dot(b));
	}

	// Z-score valid residuals
	double rmean = 0.0;
	int cnt = 0;
	for (int i = 0; i < n; ++i)
		if (!std::isnan(r(i))) { rmean += r(i); ++cnt; }
	rmean /= cnt;

	double var = 0.0;
	for (int i = 0; i < n; ++i)
		if (!std::isnan(r(i))) { double d = r(i) - rmean; var += d * d; }
	var /= (cnt - 1);
	double sd = std::sqrt(var);
	if (sd == 0.0)
		throw std::runtime_error("residualize_and_zscore_vector: zero variance after residualization");

	for (int i = 0; i < n; ++i)
		if (!std::isnan(r(i)))
			r(i) = (float)(((double)r(i) - rmean) / sd);

	return r;
}

Eigen::MatrixXf residualize_and_zscore_subset(
	const Eigen::MatrixXf& X,
	const Eigen::VectorXf& h,
	const std::vector<int>& win_idx,
	int& n_valid_windows
) {
	const int nsamples = (int)X.rows();
	const int m = (int)win_idx.size();

	Eigen::MatrixXf Z(nsamples, m);
	Z.setConstant(std::numeric_limits<float>::quiet_NaN());

	// Compute mean HI (ignore NaNs)
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
		throw std::runtime_error("HI is all missing; cannot residualize.");
	}
	hmean /= hcount2;

	// Centered HI vector
	Eigen::VectorXf hc(nsamples);
	for (int i = 0; i < nsamples; ++i) {
		float hv = h(i);
		if (std::isnan(hv))
			hc(i) = std::numeric_limits<float>::quiet_NaN();
		else
			hc(i) = hv - (float)hmean;
	}

	// Denominator: sum(hc^2) over non-missing HI
	double den = 0.0;
	for (int i = 0; i < nsamples; ++i) {
		float v = hc(i);
		if (!std::isnan(v))
			den += (double)v * (double)v;
	}
	if (den <= 0.0) {
		throw std::runtime_error("HI variance is zero; cannot residualize.");
	}

	n_valid_windows = 0;

	for (int j = 0; j < m; ++j) {
		int w = win_idx[(size_t)j];

		bool ok = true;

		// Mean of X in this window (requires no missing X and no missing HI)
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

		if (!ok || xcount < 3)
			continue;

		xmean /= xcount;

		// Slope b = cov(xc, hc) / var(hc)
		double num = 0.0;
		for (int i = 0; i < nsamples; ++i) {
			double xc = (double)X(i, w) - xmean;
			double hcv = (double)hc(i);
			num += xc * hcv;
		}
		double b = num / den;

		// Residual mean
		double rmean = 0.0;
		for (int i = 0; i < nsamples; ++i) {
			double r = ((double)X(i, w) - xmean) - b * (double)hc(i);
			rmean += r;
		}
		rmean /= nsamples;

		// Residual sd
		double rss = 0.0;
		for (int i = 0; i < nsamples; ++i) {
			double r = ((double)X(i, w) - xmean) - b * (double)hc(i);
			double d = r - rmean;
			rss += d * d;
		}

		double sd = std::sqrt(rss / (nsamples - 1));
		if (!(sd > 0.0) || !std::isfinite(sd))
			continue;

		for (int i = 0; i < nsamples; ++i) {
			double r = ((double)X(i, w) - xmean) - b * (double)hc(i);
			Z(i, j) = (float)((r - rmean) / sd);
		}

		++n_valid_windows;
	}

	return Z;
}
