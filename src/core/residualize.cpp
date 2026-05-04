#include "residualize.hpp"

#include <cmath>
#include <stdexcept>
#include <limits>

Eigen::MatrixXf residualize_and_zscore(
	const Eigen::MatrixXf& X,
	const Eigen::VectorXf& h,
	bool phased,
	int& n_valid_windows,
	const std::vector<MarkerFreq>* freqs
) {
	const int nrows = (int)X.rows();
	const int nwin  = (int)X.cols();
	const float dos_scale = phased ? 1.0f : 2.0f;

	Eigen::MatrixXf Z(nrows, nwin);
	Z.setConstant(std::numeric_limits<float>::quiet_NaN());
	n_valid_windows = 0;

	for (int w = 0; w < nwin; ++w) {
		float offset = 0.0f;
		float scale  = dos_scale;
		if (freqs) {
			const auto& f = (*freqs)[w];
			const float delta = f.p1 - f.p2;
			offset = dos_scale * f.p2;
			scale  = dos_scale * delta;
		}

		bool ok = true;
		for (int r = 0; r < nrows; ++r)
			if (std::isnan(X(r, w))) { ok = false; break; }
		if (!ok) continue;

		// Plug-in residual: r(row,w) = X(row,w) - offset - scale * h(sample_idx)
		double rmean = 0.0;
		for (int r = 0; r < nrows; ++r) {
			const int si = phased ? r / 2 : r;
			rmean += (double)X(r, w) - (double)offset - (double)scale * (double)h(si);
		}
		rmean /= nrows;

		double rss = 0.0;
		for (int r = 0; r < nrows; ++r) {
			const int si = phased ? r / 2 : r;
			const double res = (double)X(r, w) - (double)offset - (double)scale * (double)h(si) - rmean;
			rss += res * res;
		}

		const double sd = std::sqrt(rss / (nrows - 1));
		if (!(sd > 0.0) || !std::isfinite(sd))
			continue;

		for (int r = 0; r < nrows; ++r) {
			const int si = phased ? r / 2 : r;
			const double res = (double)X(r, w) - (double)offset - (double)scale * (double)h(si) - rmean;
			Z(r, w) = (float)(res / sd);
		}
		++n_valid_windows;
	}

	return Z;
}

Eigen::VectorXf residualize_and_zscore_vector(
	const Eigen::VectorXf& y,
	const Eigen::VectorXf& h,
	int& n_valid
) {
	const int n = (int)y.size();

	// Plug-in residual: r_i = y_i - h_i (scale=1, haploid 0/1)
	Eigen::VectorXf r(n);
	n_valid = 0;
	double rmean = 0.0;

	for (int i = 0; i < n; ++i) {
		if (std::isnan(y(i))) {
			r(i) = std::numeric_limits<float>::quiet_NaN();
		} else {
			r(i) = y(i) - h(i);
			rmean += (double)r(i);
			++n_valid;
		}
	}

	if (n_valid < 3)
		throw std::runtime_error("residualize_and_zscore_vector: too few valid samples");

	rmean /= n_valid;

	double rss = 0.0;
	for (int i = 0; i < n; ++i)
		if (!std::isnan(r(i))) {
			const double d = (double)r(i) - rmean;
			rss += d * d;
		}

	const double sd = std::sqrt(rss / (n_valid - 1));
	if (!(sd > 0.0) || !std::isfinite(sd))
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
	int& n_valid_windows,
	const std::vector<MarkerFreq>* freqs
) {
	const int nsamples = (int)X.rows();
	const int m = (int)win_idx.size();

	Eigen::MatrixXf Z(nsamples, m);
	Z.setConstant(std::numeric_limits<float>::quiet_NaN());
	n_valid_windows = 0;

	for (int j = 0; j < m; ++j) {
		const int w = win_idx[j];

		float offset = 0.0f, scale = 2.0f;
		if (freqs) {
			const auto& f = (*freqs)[w];
			const float delta = f.p1 - f.p2;
			offset = 2.0f * f.p2;
			scale  = 2.0f * delta;
		}

		bool ok = true;
		for (int i = 0; i < nsamples; ++i)
			if (std::isnan(X(i, w))) { ok = false; break; }
		if (!ok || nsamples < 3) continue;

		double rmean = 0.0;
		for (int i = 0; i < nsamples; ++i)
			rmean += (double)X(i, w) - (double)offset - (double)scale * (double)h(i);
		rmean /= nsamples;

		double rss = 0.0;
		for (int i = 0; i < nsamples; ++i) {
			const double res = (double)X(i, w) - (double)offset - (double)scale * (double)h(i) - rmean;
			rss += res * res;
		}

		const double sd = std::sqrt(rss / (nsamples - 1));
		if (!(sd > 0.0) || !std::isfinite(sd)) continue;

		for (int i = 0; i < nsamples; ++i) {
			const double res = (double)X(i, w) - (double)offset - (double)scale * (double)h(i) - rmean;
			Z(i, j) = (float)(res / sd);
		}
		++n_valid_windows;
	}

	return Z;
}
