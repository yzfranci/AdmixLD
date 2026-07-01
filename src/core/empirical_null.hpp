#pragma once

#include <cstdint>
#include <vector>

// Empirical-null recalibration + Storey q-value / Efron local-fdr helpers,
// operating on a reservoir sample of Fisher z = atanh(r) values drawn from a
// single chromosome-pair block.

struct NullFit {
	bool ok = false;

	long long n_pairs = 0;	// exact pair count in the block (m)
	int reservoir_size = 0;	// R

	double mu0 = 0.0;
	double sigma0 = 0.0;
	double lambda = 0.0;	// inflation factor vs theoretical null

	double pi0_pos = 1.0;
	double pi0_neg = 1.0;
	double pi0_2sided = 1.0;

	// Critical zstar thresholds for the requested target FDR.
	// A pair is a positive-tail hit if zstar >= zstar_thresh_pos,
	// a negative-tail hit if zstar <= zstar_thresh_neg.
	// Set to +inf / -inf when no threshold achieves the target (no hits).
	double zstar_thresh_pos = 0.0;
	double zstar_thresh_neg = 0.0;

	// Sorted (p, q) arrays per tail, for q-value interpolation in pass 2.
	std::vector<double> p_sorted_pos, q_sorted_pos;
	std::vector<double> p_sorted_neg, q_sorted_neg;

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

// q-value for an observed one-sided p-value, interpolated from the fitted
// tail's sorted (p, q) curve. tail_pos selects which tail's curve to use.
double qvalue_for_p(const NullFit& fit, bool tail_pos, double p);

// Efron-style local fdr at a given z, using the two-sided pi0 and the
// smoothed marginal density fit.
double local_fdr_for_z(const NullFit& fit, double z);

double norm_cdf(double x);
double norm_ppf(double p);
