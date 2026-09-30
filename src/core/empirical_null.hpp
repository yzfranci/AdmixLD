#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

// Empirical-null recalibration + Storey q-value / Efron local-fdr helpers,
// operating on a reservoir sample of Fisher z = atanh(r) values drawn from a
// single chromosome-pair block.

// One tail of the fitted null. Pairs are ordered by their extremeness e,
// which is zstar in the positive tail and -zstar in the negative tail, so
// both tails share the one-sided p-value p = Phi(-e).
struct NullTail {
	// Reservoir e values (descending) and their monotonized Storey q-values
	// (non-decreasing), used for pairs outside the exactly ranked tail.
	std::vector<double> e_sorted;
	std::vector<double> q_sorted;

	// Least extreme reservoir e whose q-value is <= the target FDR
	// (+inf when none is).
	double e_thresh = std::numeric_limits<double>::infinity();
};

struct NullFit {
	bool ok = false;

	long long n_pairs = 0;	// exact pair count in the block (m)
	int reservoir_size = 0;	// R
	double target_fdr = 0.0;

	double mu0 = 0.0;
	double sigma0 = 0.0;
	double lambda = 0.0;	// inflation factor vs theoretical null

	double pi0_pos = 1.0;
	double pi0_neg = 1.0;
	double pi0_2sided = 1.0;

	NullTail tail_pos, tail_neg;

	// Smoothed marginal density of z (histogram + Gaussian smoothing),
	// for local_fdr evaluation.
	double hist_lo = 0.0, hist_hi = 0.0, hist_bw = 1.0;
	std::vector<double> hist_density;	// smoothed density per bin
};

// Fit an empirical null from a reservoir of z = atanh(r) values sampled
// (uniformly, i.i.d.) from a chromosome-pair block of n_pairs total tests.
// n_samples is the number of individuals used in the correlation (for lambda).
NullFit fit_empirical_null(
	const std::vector<double>& z_reservoir,
	long long n_pairs,
	int n_samples,
	double target_fdr,
	double lambda_cut,
	int nbins = 120
);

// Reservoir q-value for a pair of extremeness e: the q-value of the most
// extreme reservoir point that is no more extreme than the pair (1 if none).
double reservoir_qvalue(const NullTail& tail, double e);

// Efron-style local fdr at a given z, using the two-sided pi0 and the
// smoothed marginal density fit.
double local_fdr_for_z(const NullFit& fit, double z);

double norm_cdf(double x);

// RNG for one block's calibration reservoir, so the sample depends only on
// --seed and the block, not on which thread happens to run the block.
std::mt19937_64 fdr_block_rng(uint64_t seed, uint64_t block_id);

// Pass-2 hit caller for one block. The max_exact most extreme pairs of each
// tail are kept, and once the block is done they get Storey q-values from
// their exact rank among all n_pairs pairs. A pair pushed out of (or never
// entering) that set falls back to the reservoir's q-value curve and is
// emitted right away if it passes.
//
// emit(a, b, r, z, zstar, pvalue, qvalue, local_fdr) is called once per hit.
class FdrBlockCaller {
public:
	FdrBlockCaller(const NullFit& fit, int max_exact)
		: fit_(fit), max_exact_(std::max(max_exact, 0)) {}

	long long hits_pos = 0;
	long long hits_neg = 0;

	template <typename EmitFn>
	void add(int a, int b, float r, EmitFn&& emit) {
		const double zstar = (fisher_z(r) - fit_.mu0) / fit_.sigma0;
		if (std::isnan(zstar))
			return;

		const bool pos = zstar >= 0.0;
		Candidate c{pos ? zstar : -zstar, a, b, r};
		std::vector<Candidate>& top = pos ? top_pos_ : top_neg_;

		if ((int)top.size() < max_exact_) {
			top.push_back(c);
			std::push_heap(top.begin(), top.end(), more_extreme);
			return;
		}
		if (!top.empty() && c.e > top.front().e) {
			std::pop_heap(top.begin(), top.end(), more_extreme);
			std::swap(c, top.back());
			std::push_heap(top.begin(), top.end(), more_extreme);
		}

		const NullTail& tail = pos ? fit_.tail_pos : fit_.tail_neg;
		if (c.e >= tail.e_thresh)
			write(c, pos, reservoir_qvalue(tail, c.e), emit);
	}

	template <typename EmitFn>
	void finish(EmitFn&& emit) {
		finish_tail(top_pos_, true, emit);
		finish_tail(top_neg_, false, emit);
	}

private:
	struct Candidate {
		double e;
		int a;
		int b;
		float r;
	};

	// Same clamp as the calibration pass, so zstar matches the reservoir's.
	static double fisher_z(float r) {
		const float rc = std::min(std::max(r, -0.999999999f), 0.999999999f);
		return std::atanh((double)rc);
	}

	// As a heap comparator this keeps the least extreme candidate at the front.
	static bool more_extreme(const Candidate& x, const Candidate& y) {
		return x.e > y.e;
	}

	template <typename EmitFn>
	void write(const Candidate& c, bool pos, double qvalue, EmitFn& emit) {
		const double z = fisher_z(c.r);
		emit(c.a, c.b, c.r, z, pos ? c.e : -c.e, norm_cdf(-c.e), qvalue, local_fdr_for_z(fit_, z));
		if (pos) ++hits_pos; else ++hits_neg;
	}

	template <typename EmitFn>
	void finish_tail(std::vector<Candidate>& top, bool pos, EmitFn& emit) {
		if (top.empty())
			return;

		std::sort(top.begin(), top.end(), more_extreme);

		const NullTail& tail = pos ? fit_.tail_pos : fit_.tail_neg;
		const double pi0 = pos ? fit_.pi0_pos : fit_.pi0_neg;
		const size_t K = top.size();

		// Monotonize from the least extreme kept pair up; the reservoir curve
		// stands in for the pairs beyond the kept set.
		std::vector<double> q(K);
		double running = reservoir_qvalue(tail, top[K - 1].e);
		for (size_t k = K; k-- > 0;) {
			const double raw = pi0 * (double)fit_.n_pairs * norm_cdf(-top[k].e) / (double)(k + 1);
			running = std::min(running, raw);
			q[k] = std::min(std::max(running, 0.0), 1.0);
		}

		std::vector<size_t> hits;
		for (size_t k = 0; k < K && q[k] <= fit_.target_fdr; ++k)
			hits.push_back(k);
		std::sort(hits.begin(), hits.end(), [&](size_t x, size_t y) {
			return top[x].a != top[y].a ? top[x].a < top[y].a : top[x].b < top[y].b;
		});
		for (size_t k : hits)
			write(top[k], pos, q[k], emit);

		top.clear();
	}

	const NullFit& fit_;
	int max_exact_;
	std::vector<Candidate> top_pos_, top_neg_;
};
