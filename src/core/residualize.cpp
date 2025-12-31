#include "residualize.hpp"

#include <cmath>
#include <iostream>
#include <limits>

Eigen::MatrixXf residualize_and_zscore(
	const Eigen::MatrixXf& X,
	const Eigen::VectorXf& h,
	int& n_valid_windows
) {
	int nsamples = (int)X.rows();
	int nwin = (int)X.cols();

	Eigen::MatrixXf Z(nsamples, nwin);
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

	for (int w = 0; w < nwin; ++w) {
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
			Z(i, w) = (float)((r - rmean) / sd);
		}

		++n_valid_windows;
	}

	return Z;
}
