#include "empirical_null.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
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

// e_raw holds the reservoir's extremeness values for one tail (zstar or
// -zstar); p = Phi(-e) keeps full precision far out in either tail.
static void fit_tail(
	const std::vector<double>& e_raw,
	double lambda_cut,
	double target_fdr,
	double& pi0_out,
	NullTail& tail
) {
	const int R = (int)e_raw.size();

	tail.e_sorted = e_raw;
	std::sort(tail.e_sorted.begin(), tail.e_sorted.end(), std::greater<double>());

	std::vector<double> p((size_t)R);
	int n_above = 0;
	for (int k = 0; k < R; ++k) {
		p[(size_t)k] = norm_cdf(-tail.e_sorted[(size_t)k]);
		if (p[(size_t)k] > lambda_cut)
			++n_above;
	}

	double pi0 = 1.0;
	if (R > 0 && lambda_cut < 1.0)
		pi0 = (double)n_above / ((double)R * (1.0 - lambda_cut));
	pi0 = std::min(std::max(pi0, 0.0), 1.0);
	pi0_out = pi0;

	tail.q_sorted.assign((size_t)R, 1.0);
	for (int k = 1; k <= R; ++k)
		tail.q_sorted[(size_t)(k - 1)] = pi0 * (double)R * p[(size_t)(k - 1)] / (double)k;
	monotonize_q(tail.q_sorted);

	// Largest k with q_sorted[k-1] <= target_fdr (q_sorted is non-decreasing).
	auto it = std::upper_bound(tail.q_sorted.begin(), tail.q_sorted.end(), target_fdr);
	size_t n_hits = (size_t)std::distance(tail.q_sorted.begin(), it);
	tail.e_thresh = (n_hits > 0)
		? tail.e_sorted[n_hits - 1]
		: std::numeric_limits<double>::infinity();
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
	fit.target_fdr = target_fdr;

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
	std::vector<double> e_pos((size_t)R), e_neg((size_t)R);
	for (int i = 0; i < R; ++i) {
		e_pos[(size_t)i] = (z_reservoir[(size_t)i] - fit.mu0) / fit.sigma0;
		e_neg[(size_t)i] = -e_pos[(size_t)i];
	}

	fit_tail(e_pos, lambda_cut, target_fdr, fit.pi0_pos, fit.tail_pos);
	fit_tail(e_neg, lambda_cut, target_fdr, fit.pi0_neg, fit.tail_neg);

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

double reservoir_qvalue(const NullTail& tail, double e) {
	auto it = std::lower_bound(tail.e_sorted.begin(), tail.e_sorted.end(), e,
		[](double x, double value) { return x > value; });
	if (it == tail.e_sorted.end())
		return 1.0;
	return tail.q_sorted[(size_t)std::distance(tail.e_sorted.begin(), it)];
}

std::mt19937_64 fdr_block_rng(uint64_t seed, uint64_t block_id) {
	// splitmix64 finalizer, so neighbouring blocks get unrelated streams.
	uint64_t x = seed + 0x9e3779b97f4a7c15ULL * (block_id + 1);
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
	return std::mt19937_64(x ^ (x >> 31));
}

double local_fdr_for_z(const NullFit& fit, double z) {
	if (fit.hist_density.empty() || fit.hist_bw <= 0.0)
		return 1.0;

	// Clamp before the int cast: z is +-inf for pairs with |r| = 1.
	const double last = (double)fit.hist_density.size() - 1.0;
	const double bpos = (z - fit.hist_lo) / fit.hist_bw;
	const int b = std::isnan(bpos) ? 0 : (int)std::min(std::max(bpos, 0.0), last);
	double f = fit.hist_density[(size_t)b];
	double f0 = norm_pdf(z, fit.mu0, fit.sigma0);

	if (f <= 0.0)
		return 1.0;

	return std::min(std::max(fit.pi0_2sided * f0 / f, 0.0), 1.0);
}
