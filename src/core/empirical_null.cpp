#include "empirical_null.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

double norm_cdf(double x) {
	return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

static double norm_pdf(double x, double mu, double sigma) {
	if (sigma <= 0.0)
		return 0.0;
	const double z = (x - mu) / sigma;
	static const double inv_sqrt_2pi = 0.3989422804014327;
	return (inv_sqrt_2pi / sigma) * std::exp(-0.5 * z * z);
}

// Peter Acklam's rational approximation to the inverse standard normal CDF.
double norm_ppf(double p) {
	if (p <= 0.0)
		return -std::numeric_limits<double>::infinity();
	if (p >= 1.0)
		return std::numeric_limits<double>::infinity();

	static const double a[6] = {
		-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
		1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00
	};
	static const double b[5] = {
		-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
		6.680131188771972e+01, -1.328068155288572e+01
	};
	static const double c[6] = {
		-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
		-2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00
	};
	static const double d[4] = {
		7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
		3.754408661907416e+00
	};

	const double p_low = 0.02425;
	const double p_high = 1.0 - p_low;

	double q, r, x;

	if (p < p_low) {
		q = std::sqrt(-2.0 * std::log(p));
		x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
			((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
	} else if (p <= p_high) {
		q = p - 0.5;
		r = q * q;
		x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
			(((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
	} else {
		q = std::sqrt(-2.0 * std::log(1.0 - p));
		x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
			((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
	}

	return x;
}

static double median_inplace(std::vector<double>& v) {
	if (v.empty())
		return std::numeric_limits<double>::quiet_NaN();

	const size_t n = v.size();
	const size_t mid = n / 2;
	std::nth_element(v.begin(), v.begin() + (long)mid, v.end());
	double hi = v[mid];

	if (n % 2 == 1)
		return hi;

	std::nth_element(v.begin(), v.begin() + (long)(mid - 1), v.begin() + (long)mid);
	double lo = v[mid - 1];
	return 0.5 * (lo + hi);
}

// Monotonize a Storey q-value sequence sorted by ascending p: running min
// from the largest p down, so the result is non-decreasing in k.
static void monotonize_q(std::vector<double>& q) {
	for (size_t k = q.size(); k >= 2; --k)
		q[k - 2] = std::min(q[k - 2], q[k - 1]);
	for (auto& x : q)
		x = std::min(std::max(x, 0.0), 1.0);
}

static void fit_tail(
	const std::vector<double>& p_raw,
	double lambda_cut,
	double target_fdr,
	double& pi0_out,
	std::vector<double>& p_sorted,
	std::vector<double>& q_sorted,
	double& p_thresh
) {
	const int R = (int)p_raw.size();

	int n_above = 0;
	for (double p : p_raw)
		if (p > lambda_cut)
			++n_above;

	double pi0 = 1.0;
	if (R > 0 && lambda_cut < 1.0)
		pi0 = (double)n_above / ((double)R * (1.0 - lambda_cut));
	pi0 = std::min(std::max(pi0, 0.0), 1.0);
	pi0_out = pi0;

	p_sorted = p_raw;
	std::sort(p_sorted.begin(), p_sorted.end());

	q_sorted.assign((size_t)R, 1.0);
	for (int k = 1; k <= R; ++k)
		q_sorted[(size_t)(k - 1)] = pi0 * (double)R * p_sorted[(size_t)(k - 1)] / (double)k;
	monotonize_q(q_sorted);

	// Largest k with q_sorted[k-1] <= target_fdr (q_sorted is non-decreasing).
	auto it = std::upper_bound(q_sorted.begin(), q_sorted.end(), target_fdr);
	size_t n_hits = (size_t)std::distance(q_sorted.begin(), it);
	p_thresh = (n_hits > 0) ? p_sorted[n_hits - 1] : -1.0;
}

NullFit fit_empirical_null(
	const std::vector<double>& z_reservoir,
	long long n_pairs,
	int n_samples,
	double target_fdr,
	double lambda_cut,
	int nbins
) {
	NullFit fit;
	fit.n_pairs = n_pairs;
	fit.reservoir_size = (int)z_reservoir.size();

	if (z_reservoir.empty())
		return fit;

	std::vector<double> zs = z_reservoir;
	fit.mu0 = median_inplace(zs);

	std::vector<double> absdev(z_reservoir.size());
	for (size_t i = 0; i < z_reservoir.size(); ++i)
		absdev[i] = std::fabs(z_reservoir[i] - fit.mu0);
	fit.sigma0 = 1.4826 * median_inplace(absdev);

	if (!(fit.sigma0 > 0.0)) {
		fit.ok = false;
		return fit;
	}

	fit.lambda = (n_samples > 3)
		? fit.sigma0 * std::sqrt((double)(n_samples - 3))
		: std::numeric_limits<double>::quiet_NaN();

	const int R = (int)z_reservoir.size();
	std::vector<double> zstar(z_reservoir.size());
	std::vector<double> p_pos(z_reservoir.size()), p_neg(z_reservoir.size());
	for (int i = 0; i < R; ++i) {
		zstar[(size_t)i] = (z_reservoir[(size_t)i] - fit.mu0) / fit.sigma0;
		p_pos[(size_t)i] = 1.0 - norm_cdf(zstar[(size_t)i]);
		p_neg[(size_t)i] = norm_cdf(zstar[(size_t)i]);
	}

	double p_thresh_pos = -1.0, p_thresh_neg = -1.0;
	fit_tail(p_pos, lambda_cut, target_fdr, fit.pi0_pos, fit.p_sorted_pos, fit.q_sorted_pos, p_thresh_pos);
	fit_tail(p_neg, lambda_cut, target_fdr, fit.pi0_neg, fit.p_sorted_neg, fit.q_sorted_neg, p_thresh_neg);

	fit.zstar_thresh_pos = (p_thresh_pos >= 0.0)
		? norm_ppf(1.0 - p_thresh_pos)
		: std::numeric_limits<double>::infinity();
	fit.zstar_thresh_neg = (p_thresh_neg >= 0.0)
		? norm_ppf(p_thresh_neg)
		: -std::numeric_limits<double>::infinity();

	// Smoothed marginal density of z (for local_fdr), and two-sided pi0
	// (Efron density-ratio at the null mode).
	double zmin = *std::min_element(z_reservoir.begin(), z_reservoir.end());
	double zmax = *std::max_element(z_reservoir.begin(), z_reservoir.end());
	if (zmax <= zmin) {
		zmax = zmin + 1.0;
		zmin -= 1.0;
	}
	nbins = std::max(nbins, 10);

	std::vector<double> counts((size_t)nbins, 0.0);
	const double bw = (zmax - zmin) / (double)nbins;
	for (double z : z_reservoir) {
		int b = (int)((z - zmin) / bw);
		b = std::min(std::max(b, 0), nbins - 1);
		counts[(size_t)b] += 1.0;
	}

	// Gaussian smoothing of the histogram (kernel sigma = 1.5 bins).
	const double smooth_sigma = 1.5;
	const int half_win = std::max(1, (int)std::ceil(4.0 * smooth_sigma));
	std::vector<double> kernel((size_t)(2 * half_win + 1));
	double ksum = 0.0;
	for (int k = -half_win; k <= half_win; ++k) {
		double w = std::exp(-0.5 * (double)(k * k) / (smooth_sigma * smooth_sigma));
		kernel[(size_t)(k + half_win)] = w;
		ksum += w;
	}
	for (auto& w : kernel)
		w /= ksum;

	std::vector<double> smoothed((size_t)nbins, 0.0);
	for (int i = 0; i < nbins; ++i) {
		double acc = 0.0;
		for (int k = -half_win; k <= half_win; ++k) {
			int j = i + k;
			j = std::min(std::max(j, 0), nbins - 1);	// clamp edges
			acc += counts[(size_t)j] * kernel[(size_t)(k + half_win)];
		}
		smoothed[(size_t)i] = acc;
	}

	double total = 0.0;
	for (double c : counts)
		total += c;
	fit.hist_lo = zmin;
	fit.hist_hi = zmax;
	fit.hist_bw = bw;
	fit.hist_density.assign((size_t)nbins, 0.0);
	if (total > 0.0 && bw > 0.0)
		for (int i = 0; i < nbins; ++i)
			fit.hist_density[(size_t)i] = smoothed[(size_t)i] / (total * bw);

	int b0 = (int)((fit.mu0 - zmin) / bw);
	b0 = std::min(std::max(b0, 0), nbins - 1);
	double f_mu0 = fit.hist_density.empty() ? 0.0 : fit.hist_density[(size_t)b0];
	double f0_mu0 = norm_pdf(fit.mu0, fit.mu0, fit.sigma0);
	fit.pi0_2sided = (f0_mu0 > 0.0) ? std::min(std::max(f_mu0 / f0_mu0, 0.0), 1.0) : 1.0;

	fit.ok = true;
	return fit;
}

double qvalue_for_p(const NullFit& fit, bool tail_pos, double p) {
	const std::vector<double>& ps = tail_pos ? fit.p_sorted_pos : fit.p_sorted_neg;
	const std::vector<double>& qs = tail_pos ? fit.q_sorted_pos : fit.q_sorted_neg;

	if (ps.empty())
		return 1.0;

	auto it = std::lower_bound(ps.begin(), ps.end(), p);
	if (it == ps.begin())
		return qs.front();
	if (it == ps.end())
		return 1.0;

	size_t hi = (size_t)std::distance(ps.begin(), it);
	size_t lo = hi - 1;
	double p_lo = ps[lo], p_hi = ps[hi];
	if (p_hi <= p_lo)
		return qs[lo];

	double t = (p - p_lo) / (p_hi - p_lo);
	return qs[lo] + t * (qs[hi] - qs[lo]);
}

double local_fdr_for_z(const NullFit& fit, double z) {
	if (fit.hist_density.empty() || fit.hist_bw <= 0.0)
		return 1.0;

	int b = (int)((z - fit.hist_lo) / fit.hist_bw);
	b = std::min(std::max(b, 0), (int)fit.hist_density.size() - 1);
	double f = fit.hist_density[(size_t)b];
	double f0 = norm_pdf(z, fit.mu0, fit.sigma0);

	if (f <= 0.0)
		return 1.0;

	return std::min(std::max(fit.pi0_2sided * f0 / f, 0.0), 1.0);
}
