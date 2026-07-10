#include "scan_dynamic_hi.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include <cstdio>

#ifdef ADMIXLD_HAS_OPENMP
#include <omp.h>
#endif

#include "residualize.hpp"
#include "empirical_null.hpp"

static float quantile_from_sorted(
	const std::vector<float>& v,
	double p
) {
	if (v.empty())
		return std::numeric_limits<float>::quiet_NaN();

	double x = p * (v.size() - 1);
	size_t i0 = (size_t)std::floor(x);
	size_t i1 = (size_t)std::ceil(x);
	if (i0 == i1)
		return v[i0];

	double t = x - (double)i0;
	return (float)((1.0 - t) * v[i0] + t * v[i1]);
}

static bool keep_hit_r(
	float r,
	const ScanOptions& opt
) {
	if (!std::isfinite(r))
		return false;

	if (opt.use_asym) {
		bool keep = false;
		if (opt.min_neg_r > 0.0f && r <= -opt.min_neg_r)
			keep = true;
		if (opt.min_pos_r > 0.0f && r >= opt.min_pos_r)
			keep = true;
		return keep;
	}

	return (std::fabs(r) >= opt.min_abs_r);
}

struct ThreadDistrib {
	long long tested_pairs = 0;
	float min_r = std::numeric_limits<float>::infinity();
	float max_r = -std::numeric_limits<float>::infinity();
	std::mt19937_64 rng;
	std::vector<float> sample;
};

static void consider_distrib(
	ThreadDistrib& d,
	float r,
	int sample_cap
) {
	if (!std::isfinite(r))
		return;

	if (r < d.min_r) d.min_r = r;
	if (r > d.max_r) d.max_r = r;

	++d.tested_pairs;

	if ((int)d.sample.size() < sample_cap) {
		d.sample.push_back(r);
	} else {
		std::uniform_int_distribution<long long> U(0, d.tested_pairs - 1);
		long long j = U(d.rng);
		if (j < sample_cap)
			d.sample[(size_t)j] = r;
	}
}

static void merge_part_files(
	const std::string& out_path,
	const std::vector<std::string>& part_paths,
	const std::string& header
) {
	std::ofstream of(out_path);
	if (!of)
		throw std::runtime_error("cannot write to " + out_path);

	of << header;

	for (size_t t = 0; t < part_paths.size(); ++t) {
		std::ifstream pf(part_paths[t]);
		if (!pf)
			continue;

		std::string line;
		while (std::getline(pf, line))
			of << line << "\n";
	}
}

// Concatenates two TSV files into out_path, keeping only the first header.
static bool concat_tsv_drop_second_header(
	const std::string& out_path,
	const std::string& first_path,
	const std::string& second_path
) {
	std::ofstream of(out_path);
	if (!of) {
		std::cerr << "Error: cannot write to " << out_path << "\n";
		return false;
	}
	{
		std::ifstream f1(first_path);
		std::string line;
		while (std::getline(f1, line))
			of << line << "\n";
	}
	{
		std::ifstream f2(second_path);
		std::string line;
		bool first = true;
		while (std::getline(f2, line)) {
			if (first) { first = false; continue; }
			of << line << "\n";
		}
	}
	return true;
}

static void compute_mean_sd_from_sample(
	const std::vector<float>& v,
	float& mean_out,
	float& sd_out
) {
	mean_out = std::numeric_limits<float>::quiet_NaN();
	sd_out = std::numeric_limits<float>::quiet_NaN();

	long long n = 0;
	double mu = 0.0;
	double m2 = 0.0;

	for (float x : v) {
		if (!std::isfinite(x))
			continue;

		++n;
		double dx = (double)x - mu;
		mu += dx / (double)n;
		double dx2 = (double)x - mu;
		m2 += dx * dx2;
	}

	if (n <= 0)
		return;

	mean_out = (float)mu;
	if (n > 1) {
		double var = m2 / (double)(n - 1);
		if (var < 0.0) var = 0.0;
		sd_out = (float)std::sqrt(var);
	}
}

template<typename HC, typename HiExcludeFn>
static bool scan_markers_write_hits_excl_focus_T(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HC& hc_full,
	HiExcludeFn hi_excluding,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path,
	const std::vector<MarkerFreq>* freqs_scan = nullptr
) {
	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;

	tested_pairs = 0;
	kept_pairs = 0;

	const bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
	nthreads = opt.threads;
	if (nthreads < 1) nthreads = 1;
	omp_set_num_threads(nthreads);
	#endif

	if (nthreads > (int)chr_order_scan.size())
		nthreads = (int)chr_order_scan.size();
	if (nthreads < 1)
		nthreads = 1;

	std::vector<std::string> part_paths((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t) {
		part_paths[(size_t)t] = out_path + ".part." + std::to_string(t);
		std::remove(part_paths[(size_t)t].c_str());
	}

	std::vector<long long> tested_t((size_t)nthreads, 0);
	std::vector<long long> kept_t((size_t)nthreads, 0);

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	const float denom = 1.0f / (float)(nsamples - 1);

	auto do_intra_chr = [&](int tid, const std::string& chr) {
		const auto& idx = windows_by_chr_scan.at(chr);
		const int m = (int)idx.size();
		if (m < 2)
			return;

		Eigen::VectorXf h = hi_excluding(hc_full, chr, std::string(""));

		int n_valid = 0;
		Eigen::MatrixXf Zc = residualize_and_zscore_subset(X_scan, h, idx, n_valid, freqs_scan);

		std::ofstream ofp(part_paths[(size_t)tid], std::ios::out | std::ios::app);
		if (!ofp) {
			std::cerr << "Error: cannot write to " << part_paths[(size_t)tid] << "\n";
			return;
		}

		Eigen::MatrixXf A(nsamples, tile_size);
		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::MatrixXf R(tile_size, tile_size);

		for (int i0 = 0; i0 < m; i0 += tile_size) {
			const int b1 = std::min(tile_size, m - i0);

			for (int k = 0; k < b1; ++k)
				A.col(k) = Zc.col(i0 + k);

			for (int j0 = i0; j0 < m; j0 += tile_size) {
				const int b2 = std::min(tile_size, m - j0);

				for (int k = 0; k < b2; ++k)
					B.col(k) = Zc.col(j0 + k);

				R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
				R.topLeftCorner(b1, b2) *= denom;

				for (int ia = 0; ia < b1; ++ia) {
					const int a = idx[i0 + ia];
					const int posA = pos_scan[a];

					int jb_start = 0;
					if (j0 == i0)
						jb_start = ia + 1;

					int jb_end = b2;

					if (opt.max_dist >= 0) {
						const int limit = posA + opt.max_dist;

						auto begin_it = idx.begin() + j0;
						auto end_it = idx.begin() + j0 + b2;

						auto ub = std::upper_bound(begin_it, end_it, limit,
							[&](int value, int widx) {
								return value < pos_scan[widx];
							}
						);

						jb_end = (int)std::distance(begin_it, ub);
					}

					if (jb_end <= jb_start)
						continue;

					for (int ib = jb_start; ib < jb_end; ++ib) {
						const int b = idx[j0 + ib];

						if (opt.min_dist >= 0 && (pos_scan[b] - posA) < opt.min_dist)
							continue;

						const float r = R(ia, ib);

						if (!std::isfinite(r))
							continue;

						++tested_t[(size_t)tid];
						if (do_distrib)
							consider_distrib(td[(size_t)tid], r, distrib_sample);

						if (keep_hit_r(r, opt)) {
							ofp << a << "\t" << chroms_scan[a] << "\t" << pos_scan[a] << "\t"
								<< b << "\t" << chroms_scan[b] << "\t" << pos_scan[b] << "\t"
								<< r << "\t" << nsamples << "\n";
							++kept_t[(size_t)tid];
						}
					}
				}
			}
		}
	};

	auto do_inter_pair = [&](int tid, const std::string& chr1, const std::string& chr2) {
		const auto& idx1 = windows_by_chr_scan.at(chr1);
		const auto& idx2 = windows_by_chr_scan.at(chr2);
		const int m1 = (int)idx1.size();
		const int m2 = (int)idx2.size();
		if (m1 == 0 || m2 == 0)
			return;

		Eigen::VectorXf h = hi_excluding(hc_full, chr1, chr2);

		int n_valid1 = 0;
		int n_valid2 = 0;
		Eigen::MatrixXf Z1 = residualize_and_zscore_subset(X_scan, h, idx1, n_valid1, freqs_scan);
		Eigen::MatrixXf Z2 = residualize_and_zscore_subset(X_scan, h, idx2, n_valid2, freqs_scan);

		std::ofstream ofp(part_paths[(size_t)tid], std::ios::out | std::ios::app);
		if (!ofp) {
			std::cerr << "Error: cannot write to " << part_paths[(size_t)tid] << "\n";
			return;
		}

		Eigen::MatrixXf A(nsamples, tile_size);
		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::MatrixXf R(tile_size, tile_size);

		for (int i0 = 0; i0 < m1; i0 += tile_size) {
			const int b1 = std::min(tile_size, m1 - i0);

			for (int k = 0; k < b1; ++k)
				A.col(k) = Z1.col(i0 + k);

			for (int j0 = 0; j0 < m2; j0 += tile_size) {
				const int b2 = std::min(tile_size, m2 - j0);

				for (int k = 0; k < b2; ++k)
					B.col(k) = Z2.col(j0 + k);

				R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
				R.topLeftCorner(b1, b2) *= denom;

				for (int ia = 0; ia < b1; ++ia) {
					const int a = idx1[i0 + ia];
					for (int ib = 0; ib < b2; ++ib) {
						const int b = idx2[j0 + ib];
						const float r = R(ia, ib);

						if (!std::isfinite(r))
							continue;

						++tested_t[(size_t)tid];
						if (do_distrib)
							consider_distrib(td[(size_t)tid], r, distrib_sample);

						if (keep_hit_r(r, opt)) {
							ofp << a << "\t" << chroms_scan[a] << "\t" << pos_scan[a] << "\t"
								<< b << "\t" << chroms_scan[b] << "\t" << pos_scan[b] << "\t"
								<< r << "\t" << nsamples << "\n";
							++kept_t[(size_t)tid];
						}
					}
				}
			}
		}
	};

	if (opt.intra) {
		const int C = (int)chr_order_scan.size();

		#ifdef ADMIXLD_HAS_OPENMP
		#pragma omp parallel for schedule(dynamic)
		#endif
		for (int c = 0; c < C; ++c) {
			int tid = 0;
			#ifdef ADMIXLD_HAS_OPENMP
			tid = omp_get_thread_num();
			if (tid >= nthreads) tid = tid % nthreads;
			#endif

			do_intra_chr(tid, chr_order_scan[c]);
		}
	} else {
		const int C = (int)chr_order_scan.size();

		std::vector<std::pair<int,int>> jobs;
		jobs.reserve((size_t)C * (size_t)(C - 1) / 2);
		for (int c1 = 0; c1 < C; ++c1)
			for (int c2 = c1 + 1; c2 < C; ++c2)
				jobs.push_back({c1, c2});

		#ifdef ADMIXLD_HAS_OPENMP
		#pragma omp parallel for schedule(dynamic)
		#endif
		for (int j = 0; j < (int)jobs.size(); ++j) {
			int tid = 0;
			#ifdef ADMIXLD_HAS_OPENMP
			tid = omp_get_thread_num();
			if (tid >= nthreads) tid = tid % nthreads;
			#endif

			int c1 = jobs[(size_t)j].first;
			int c2 = jobs[(size_t)j].second;
			do_inter_pair(tid, chr_order_scan[c1], chr_order_scan[c2]);
		}
	}

	for (int t = 0; t < nthreads; ++t) {
		tested_pairs += tested_t[(size_t)t];
		kept_pairs += kept_t[(size_t)t];
	}

	try {
		merge_part_files(
			out_path,
			part_paths,
			"wA\tchrA\tposA\twB\tchrB\tposB\tr\tn\n"
		);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return false;
	}

	for (int t = 0; t < nthreads; ++t)
		std::remove(part_paths[(size_t)t].c_str());

	if (do_distrib) {
		long long total_tested = 0;
		float global_min = std::numeric_limits<float>::infinity();
		float global_max = -std::numeric_limits<float>::infinity();

		std::vector<float> dsample;
		dsample.reserve((size_t)distrib_sample * (size_t)nthreads);

		for (int t = 0; t < nthreads; ++t) {
			const auto& d = td[(size_t)t];
			if (d.tested_pairs == 0)
				continue;

			total_tested += d.tested_pairs;
			if (d.min_r < global_min) global_min = d.min_r;
			if (d.max_r > global_max) global_max = d.max_r;

			dsample.insert(dsample.end(), d.sample.begin(), d.sample.end());
		}

		if (!distrib_path.empty()) {
			std::ofstream df(distrib_path);
			if (!df) {
				std::cerr << "Error: cannot write to " << distrib_path << "\n";
				return false;
			}

			df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\tmean_r2\tsd_r2\n";

			if (total_tested == 0 || dsample.empty()) {
				df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
			} else {
				if ((int)dsample.size() > distrib_sample) {
					std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
					std::shuffle(dsample.begin(), dsample.end(), rng);
					dsample.resize((size_t)distrib_sample);
				}

				float mean = std::numeric_limits<float>::quiet_NaN();
				float sd = std::numeric_limits<float>::quiet_NaN();
				compute_mean_sd_from_sample(dsample, mean, sd);

				float mean_r2 = std::numeric_limits<float>::quiet_NaN();
				float sd_r2 = std::numeric_limits<float>::quiet_NaN();
				{
					double mu2 = 0.0, m2 = 0.0;
					long long n2 = 0;
					for (float x : dsample) {
						if (!std::isfinite(x)) continue;
						++n2;
						double r2 = (double)x * (double)x;
						double dr2 = r2 - mu2;
						mu2 += dr2 / (double)n2;
						m2 += dr2 * (r2 - mu2);
					}
					if (n2 > 0) mean_r2 = (float)mu2;
					if (n2 > 1) sd_r2 = (float)std::sqrt(m2 / (double)(n2 - 1));
				}

				std::sort(dsample.begin(), dsample.end());

				float p01 = quantile_from_sorted(dsample, 0.01);
				float p05 = quantile_from_sorted(dsample, 0.05);
				float p25 = quantile_from_sorted(dsample, 0.25);
				float med = quantile_from_sorted(dsample, 0.50);
				float p75 = quantile_from_sorted(dsample, 0.75);
				float p95 = quantile_from_sorted(dsample, 0.95);
				float p99 = quantile_from_sorted(dsample, 0.99);

				df << total_tested
					<< "\t" << global_max
					<< "\t" << p99
					<< "\t" << p95
					<< "\t" << p75
					<< "\t" << med
					<< "\t" << p25
					<< "\t" << p05
					<< "\t" << p01
					<< "\t" << global_min
					<< "\t" << mean
					<< "\t" << sd
					<< "\t" << mean_r2
					<< "\t" << sd_r2
					<< "\n";
			}
		} else if (!dsample.empty() && (int)dsample.size() > distrib_sample) {
			std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
			std::shuffle(dsample.begin(), dsample.end(), rng);
			dsample.resize((size_t)distrib_sample);
		}

		if (!reservoir_path.empty() && !dsample.empty()) {
			std::ofstream rf(reservoir_path);
			if (!rf) {
				std::cerr << "Error: cannot write reservoir to " << reservoir_path << "\n";
				return false;
			}
			rf << "r\n";
			for (float x : dsample)
				rf << x << "\n";
		}
	}

	return true;
}

bool scan_markers_write_hits_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	auto fn = [](const HiComponentsWeighted& hc, const std::string& a, const std::string& b) {
		return hi_from_components_weighted_excluding(hc, a, b);
	};
	return scan_markers_write_hits_excl_focus_T(
		X_scan, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, out_path, tested_pairs, kept_pairs,
		distrib_path, distrib_sample, distrib_seed, reservoir_path
	);
}

bool scan_markers_write_hits_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::vector<MarkerFreq>& freqs_scan,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	auto fn = [](const HiComponentsFreq& hc, const std::string& a, const std::string& b) {
		return hi_from_components_freq_excluding(hc, a, b);
	};
	return scan_markers_write_hits_excl_focus_T(
		X_scan, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, out_path, tested_pairs, kept_pairs,
		distrib_path, distrib_sample, distrib_seed, reservoir_path, &freqs_scan
	);
}

// Same (chromosome, distance-bin) block structure as
// scan_markers_write_hits_fdr_intra, but each block is residualized
// excluding its own chromosome from the HI; recomputed per (chr,bin) job.
template<typename HC, typename HiExcludeFn>
static bool scan_markers_write_hits_excl_focus_fdr_intra_T(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HC& hc_full,
	HiExcludeFn hi_excluding,
	const ScanOptions& opt,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	uint64_t seed,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path,
	const std::vector<MarkerFreq>* freqs_scan
) {
	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;
	const float denom = 1.0f / (float)(nsamples - 1);
	const int rsample = opt.fdr_sample > 0 ? opt.fdr_sample : 200000;
	const int n_bins = opt.fdr_intra_bins > 0 ? opt.fdr_intra_bins : 8;
	const long long real_floor = opt.min_dist >= 0 ? (long long)opt.min_dist : 0;

	tested_pairs = 0;
	kept_pairs = 0;

	const bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	long long effective_min = std::max(real_floor, 1LL);
	long long effective_max = opt.max_dist >= 0 ? (long long)opt.max_dist : 0;
	if (effective_max <= 0) {
		for (const auto& chr : chr_order_scan) {
			const auto& idx = windows_by_chr_scan.at(chr);
			if (idx.size() < 2)
				continue;
			long long span = (long long)pos_scan[idx.back()] - (long long)pos_scan[idx.front()];
			if (span > effective_max)
				effective_max = span;
		}
	}
	if (effective_max <= effective_min)
		effective_max = effective_min + 1;

	std::vector<long long> edges((size_t)n_bins + 1);
	{
		const double log_lo = std::log((double)effective_min);
		const double log_hi = std::log((double)effective_max);
		for (int i = 0; i <= n_bins; ++i) {
			const double t = (double)i / (double)n_bins;
			edges[(size_t)i] = (long long)std::llround(std::exp(log_lo + t * (log_hi - log_lo)));
		}
		for (int i = 1; i <= n_bins; ++i)
			if (edges[(size_t)i] <= edges[(size_t)(i - 1)])
				edges[(size_t)i] = edges[(size_t)(i - 1)] + 1;
	}

	std::cout << "[scan-fdr-loco-intra] target_fdr=" << opt.fdr_target
		  << " fdr_sample=" << rsample
		  << " fdr_min_pairs=" << opt.fdr_min_pairs
		  << " fdr_lambda=" << opt.fdr_lambda_cut
		  << " bins=" << n_bins
		  << " min_dist=" << real_floor
		  << " max_dist=" << effective_max
		  << " tile_size=" << tile_size << "\n";

	const int C = (int)chr_order_scan.size();
	std::vector<std::pair<int,int>> jobs;	// (chr_index, bin_index)
	jobs.reserve((size_t)C * (size_t)n_bins);
	for (int c = 0; c < C; ++c)
		for (int bidx = 0; bidx < n_bins; ++bidx)
			jobs.push_back({c, bidx});

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
	nthreads = opt.threads;
	if (nthreads < 1) nthreads = 1;
	omp_set_num_threads(nthreads);
	#endif
	if (nthreads > (int)jobs.size())
		nthreads = std::max(1, (int)jobs.size());

	std::vector<std::string> hit_part_paths((size_t)nthreads);
	std::vector<std::string> summary_part_paths((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t) {
		hit_part_paths[(size_t)t] = out_path + ".part." + std::to_string(t);
		summary_part_paths[(size_t)t] = summary_path + ".part." + std::to_string(t);
		std::remove(hit_part_paths[(size_t)t].c_str());
		std::remove(summary_part_paths[(size_t)t].c_str());
	}

	std::vector<long long> tested_t((size_t)nthreads, 0);
	std::vector<long long> kept_t((size_t)nthreads, 0);
	std::vector<std::mt19937_64> rngs((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t)
		rngs[(size_t)t] = std::mt19937_64(seed + (uint64_t)t);

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	#ifdef ADMIXLD_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int j = 0; j < (int)jobs.size(); ++j) {
		int tid = 0;
		#ifdef ADMIXLD_HAS_OPENMP
		tid = omp_get_thread_num();
		if (tid >= nthreads) tid = tid % nthreads;
		#endif

		const int c = jobs[(size_t)j].first;
		const int bidx = jobs[(size_t)j].second;
		const auto& chr = chr_order_scan[c];
		const auto& idx = windows_by_chr_scan.at(chr);
		const int m = (int)idx.size();
		if (m < 2)
			continue;

		const long long bin_lo = (bidx == 0) ? real_floor : edges[(size_t)bidx];
		const long long bin_hi = edges[(size_t)bidx + 1];
		if (bin_lo >= bin_hi)
			continue;

		std::ofstream hfp(hit_part_paths[(size_t)tid], std::ios::out | std::ios::app);
		std::ofstream sfp(summary_part_paths[(size_t)tid], std::ios::out | std::ios::app);
		if (!hfp || !sfp) {
			#ifdef ADMIXLD_HAS_OPENMP
			#pragma omp critical
			#endif
			{
				std::cerr << "Error: cannot write part files for thread " << tid << "\n";
			}
			continue;
		}

		Eigen::VectorXf h = hi_excluding(hc_full, chr, std::string(""));

		int n_valid = 0;
		Eigen::MatrixXf Zc = residualize_and_zscore_subset(X_scan, h, idx, n_valid, freqs_scan);

		Eigen::MatrixXf A(nsamples, tile_size);
		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::MatrixXf R(tile_size, tile_size);

		auto for_each_pair = [&](auto&& cb) {
			for (int i0 = 0; i0 < m; i0 += tile_size) {
				const int b1 = std::min(tile_size, m - i0);
				for (int k = 0; k < b1; ++k)
					A.col(k) = Zc.col(i0 + k);

				for (int j0 = i0; j0 < m; j0 += tile_size) {
					const int b2 = std::min(tile_size, m - j0);
					for (int k = 0; k < b2; ++k)
						B.col(k) = Zc.col(j0 + k);

					R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
					R.topLeftCorner(b1, b2) *= denom;

					for (int ia = 0; ia < b1; ++ia) {
						const int a = idx[i0 + ia];
						const int posA = pos_scan[a];

						int jb_dedup_start = 0;
						if (j0 == i0)
							jb_dedup_start = ia + 1;

						auto begin_it = idx.begin() + j0;
						auto end_it = idx.begin() + j0 + b2;

						auto lo_it = std::lower_bound(begin_it, end_it, (long long)posA + bin_lo,
							[&](int widx, long long value) { return (long long)pos_scan[widx] < value; });
						auto hi_it = std::upper_bound(begin_it, end_it, (long long)posA + bin_hi - 1,
							[&](long long value, int widx) { return value < (long long)pos_scan[widx]; });

						const int jb_start = std::max(jb_dedup_start, (int)std::distance(begin_it, lo_it));
						const int jb_end = (int)std::distance(begin_it, hi_it);

						if (jb_end <= jb_start)
							continue;

						for (int ib = jb_start; ib < jb_end; ++ib) {
							const int b = idx[j0 + ib];
							const float r = R(ia, ib);
							if (!std::isfinite(r))
								continue;
							cb(a, b, r);
						}
					}
				}
			}
		};

		long long m_local = 0;
		std::vector<double> reservoir;
		reservoir.reserve((size_t)rsample);
		auto& rng = rngs[(size_t)tid];

		for_each_pair([&](int, int, float r) {
			++m_local;
			float rc = std::min(std::max(r, -0.999999999f), 0.999999999f);
			double z = std::atanh((double)rc);
			if ((int)reservoir.size() < rsample) {
				reservoir.push_back(z);
			} else {
				std::uniform_int_distribution<long long> U(0, m_local - 1);
				long long jr = U(rng);
				if (jr < rsample)
					reservoir[(size_t)jr] = z;
			}
			if (do_distrib)
				consider_distrib(td[(size_t)tid], r, distrib_sample);
		});

		NullFit fit = fit_empirical_null(reservoir, m_local, nsamples, opt.fdr_target, opt.fdr_lambda_cut);

		long long kept_local = 0;
		long long hits_pos = 0, hits_neg = 0;

		if (fit.ok && m_local >= opt.fdr_min_pairs) {
			for_each_pair([&](int a, int b, float r) {
				float rc = std::min(std::max(r, -0.999999999f), 0.999999999f);
				double z = std::atanh((double)rc);
				double zstar = (z - fit.mu0) / fit.sigma0;

				bool is_pos_hit = std::isfinite(fit.zstar_thresh_pos) && zstar >= fit.zstar_thresh_pos;
				bool is_neg_hit = !is_pos_hit && std::isfinite(fit.zstar_thresh_neg) && zstar <= fit.zstar_thresh_neg;

				if (!is_pos_hit && !is_neg_hit)
					return;

				double pvalue = is_pos_hit ? (1.0 - norm_cdf(zstar)) : norm_cdf(zstar);
				double qvalue = qvalue_for_p(fit, is_pos_hit, pvalue);
				double lfdr = local_fdr_for_z(fit, z);

				hfp << a << "\t" << chroms_scan[a] << "\t" << pos_scan[a] << "\t"
					<< b << "\t" << chroms_scan[b] << "\t" << pos_scan[b] << "\t"
					<< r << "\t" << nsamples << "\t"
					<< z << "\t" << zstar << "\t" << pvalue << "\t" << qvalue << "\t" << lfdr << "\n";

				++kept_local;
				if (is_pos_hit) ++hits_pos; else ++hits_neg;
			});
		}

		sfp << chr << "\t" << bin_lo << "\t" << bin_hi << "\t" << m_local << "\t"
			<< (fit.ok && m_local >= opt.fdr_min_pairs ? "ok" : "skipped") << "\t"
			<< fit.mu0 << "\t" << fit.sigma0 << "\t" << fit.lambda << "\t"
			<< fit.pi0_pos << "\t" << fit.pi0_neg << "\t"
			<< hits_pos << "\t" << hits_neg << "\n";

		tested_t[(size_t)tid] += m_local;
		kept_t[(size_t)tid] += kept_local;
	}

	for (int t = 0; t < nthreads; ++t) {
		tested_pairs += tested_t[(size_t)t];
		kept_pairs += kept_t[(size_t)t];
	}

	try {
		merge_part_files(
			out_path, hit_part_paths,
			"wA\tchrA\tposA\twB\tchrB\tposB\tr\tn\tz\tzstar\tpvalue\tqvalue\tlocal_fdr\n"
		);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return false;
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(hit_part_paths[(size_t)t].c_str());

	try {
		merge_part_files(
			summary_path, summary_part_paths,
			"chr\tdist_lo\tdist_hi\tn_pairs\tstatus\tmu0\tsigma0\tlambda\tpi0_pos\tpi0_neg\tn_hits_pos\tn_hits_neg\n"
		);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return false;
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(summary_part_paths[(size_t)t].c_str());

	if (do_distrib) {
		long long total_tested = 0;
		float global_min = std::numeric_limits<float>::infinity();
		float global_max = -std::numeric_limits<float>::infinity();

		std::vector<float> dsample;
		dsample.reserve((size_t)distrib_sample * (size_t)nthreads);

		for (int t = 0; t < nthreads; ++t) {
			const auto& d = td[(size_t)t];
			if (d.tested_pairs == 0)
				continue;

			total_tested += d.tested_pairs;
			if (d.min_r < global_min) global_min = d.min_r;
			if (d.max_r > global_max) global_max = d.max_r;

			dsample.insert(dsample.end(), d.sample.begin(), d.sample.end());
		}

		if ((int)dsample.size() > distrib_sample) {
			std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
			std::shuffle(dsample.begin(), dsample.end(), rng);
			dsample.resize((size_t)distrib_sample);
		}

		if (!distrib_path.empty()) {
			std::ofstream df(distrib_path);
			if (!df) {
				std::cerr << "Error: cannot write to " << distrib_path << "\n";
				return false;
			}

			df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\tmean_r2\tsd_r2\n";

			if (total_tested == 0 || dsample.empty()) {
				df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
			} else {
				float mean = std::numeric_limits<float>::quiet_NaN();
				float sd = std::numeric_limits<float>::quiet_NaN();
				compute_mean_sd_from_sample(dsample, mean, sd);

				float mean_r2 = std::numeric_limits<float>::quiet_NaN();
				float sd_r2 = std::numeric_limits<float>::quiet_NaN();
				{
					double mu2 = 0.0, m2v = 0.0;
					long long n2 = 0;
					for (float x : dsample) {
						if (!std::isfinite(x)) continue;
						++n2;
						double r2 = (double)x * (double)x;
						double dr2 = r2 - mu2;
						mu2 += dr2 / (double)n2;
						m2v += dr2 * (r2 - mu2);
					}
					if (n2 > 0) mean_r2 = (float)mu2;
					if (n2 > 1) sd_r2 = (float)std::sqrt(m2v / (double)(n2 - 1));
				}

				std::sort(dsample.begin(), dsample.end());

				float p01 = quantile_from_sorted(dsample, 0.01);
				float p05 = quantile_from_sorted(dsample, 0.05);
				float p25 = quantile_from_sorted(dsample, 0.25);
				float med = quantile_from_sorted(dsample, 0.50);
				float p75 = quantile_from_sorted(dsample, 0.75);
				float p95 = quantile_from_sorted(dsample, 0.95);
				float p99 = quantile_from_sorted(dsample, 0.99);

				df << total_tested
					<< "\t" << global_max
					<< "\t" << p99
					<< "\t" << p95
					<< "\t" << p75
					<< "\t" << med
					<< "\t" << p25
					<< "\t" << p05
					<< "\t" << p01
					<< "\t" << global_min
					<< "\t" << mean
					<< "\t" << sd
					<< "\t" << mean_r2
					<< "\t" << sd_r2
					<< "\n";
			}
		}

		if (!reservoir_path.empty() && !dsample.empty()) {
			std::ofstream rf(reservoir_path);
			if (!rf) {
				std::cerr << "Error: cannot write reservoir to " << reservoir_path << "\n";
				return false;
			}
			rf << "r\n";
			for (float x : dsample)
				rf << x << "\n";
		}
	}

	return true;
}

// Per chromosome-pair block, runs the calibrate+hit-call two-pass procedure
// against that block's LOCO-residualized Z1/Z2. Delegates to the intra
// variant above when opt.intra is set.
template<typename HC, typename HiExcludeFn>
static bool scan_markers_write_hits_excl_focus_fdr_T(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HC& hc_full,
	HiExcludeFn hi_excluding,
	const ScanOptions& opt,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	uint64_t seed,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path,
	const std::vector<MarkerFreq>* freqs_scan = nullptr
) {
	if (opt.intra) {
		return scan_markers_write_hits_excl_focus_fdr_intra_T(
			X_scan, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
			hc_full, hi_excluding, opt, out_path, summary_path, tested_pairs, kept_pairs, seed,
			distrib_path, distrib_sample, distrib_seed, reservoir_path, freqs_scan
		);
	}

	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;
	const float denom = 1.0f / (float)(nsamples - 1);
	const int rsample = opt.fdr_sample > 0 ? opt.fdr_sample : 200000;

	std::cout << "[scan-fdr-loco] target_fdr=" << opt.fdr_target
		  << " fdr_sample=" << rsample
		  << " fdr_min_pairs=" << opt.fdr_min_pairs
		  << " fdr_lambda=" << opt.fdr_lambda_cut
		  << " tile_size=" << tile_size << "\n";

	tested_pairs = 0;
	kept_pairs = 0;

	const bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	const int C = (int)chr_order_scan.size();
	std::vector<std::pair<int,int>> jobs;
	jobs.reserve((size_t)C * (size_t)(C - 1) / 2);
	for (int c1 = 0; c1 < C; ++c1)
		for (int c2 = c1 + 1; c2 < C; ++c2)
			jobs.push_back({c1, c2});

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
	nthreads = opt.threads;
	if (nthreads < 1) nthreads = 1;
	omp_set_num_threads(nthreads);
	#endif
	if (nthreads > (int)jobs.size())
		nthreads = std::max(1, (int)jobs.size());

	std::vector<std::string> hit_part_paths((size_t)nthreads);
	std::vector<std::string> summary_part_paths((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t) {
		hit_part_paths[(size_t)t] = out_path + ".part." + std::to_string(t);
		summary_part_paths[(size_t)t] = summary_path + ".part." + std::to_string(t);
		std::remove(hit_part_paths[(size_t)t].c_str());
		std::remove(summary_part_paths[(size_t)t].c_str());
	}

	std::vector<long long> tested_t((size_t)nthreads, 0);
	std::vector<long long> kept_t((size_t)nthreads, 0);
	std::vector<std::mt19937_64> rngs((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t)
		rngs[(size_t)t] = std::mt19937_64(seed + (uint64_t)t);

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	#ifdef ADMIXLD_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int j = 0; j < (int)jobs.size(); ++j) {
		int tid = 0;
		#ifdef ADMIXLD_HAS_OPENMP
		tid = omp_get_thread_num();
		if (tid >= nthreads) tid = tid % nthreads;
		#endif

		const int c1 = jobs[(size_t)j].first;
		const int c2 = jobs[(size_t)j].second;
		const std::string& chr1 = chr_order_scan[c1];
		const std::string& chr2 = chr_order_scan[c2];
		const auto& idx1 = windows_by_chr_scan.at(chr1);
		const auto& idx2 = windows_by_chr_scan.at(chr2);
		const int m1 = (int)idx1.size();
		const int m2 = (int)idx2.size();
		if (m1 == 0 || m2 == 0)
			continue;

		std::ofstream hfp(hit_part_paths[(size_t)tid], std::ios::out | std::ios::app);
		std::ofstream sfp(summary_part_paths[(size_t)tid], std::ios::out | std::ios::app);
		if (!hfp || !sfp) {
			#ifdef ADMIXLD_HAS_OPENMP
			#pragma omp critical
			#endif
			{
				std::cerr << "Error: cannot write part files for thread " << tid << "\n";
			}
			continue;
		}

		// Block-specific LOCO HI (excludes both chr1 and chr2); Z1/Z2 don't
		// change between the calibration and hit-calling passes below.
		Eigen::VectorXf h = hi_excluding(hc_full, chr1, chr2);

		int n_valid1 = 0, n_valid2 = 0;
		Eigen::MatrixXf Z1 = residualize_and_zscore_subset(X_scan, h, idx1, n_valid1, freqs_scan);
		Eigen::MatrixXf Z2 = residualize_and_zscore_subset(X_scan, h, idx2, n_valid2, freqs_scan);

		Eigen::MatrixXf A(nsamples, tile_size);
		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::MatrixXf R(tile_size, tile_size);

		auto for_each_pair = [&](auto&& cb) {
			for (int i0 = 0; i0 < m1; i0 += tile_size) {
				const int b1 = std::min(tile_size, m1 - i0);
				for (int k = 0; k < b1; ++k)
					A.col(k) = Z1.col(i0 + k);

				for (int j0 = 0; j0 < m2; j0 += tile_size) {
					const int b2 = std::min(tile_size, m2 - j0);
					for (int k = 0; k < b2; ++k)
						B.col(k) = Z2.col(j0 + k);

					R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
					R.topLeftCorner(b1, b2) *= denom;

					for (int ia = 0; ia < b1; ++ia) {
						const int a = idx1[i0 + ia];
						for (int ib = 0; ib < b2; ++ib) {
							const int b = idx2[j0 + ib];
							const float r = R(ia, ib);
							if (!std::isfinite(r))
								continue;
							cb(a, b, r);
						}
					}
				}
			}
		};

		// Pass 1: exact pair count + reservoir sample of z = atanh(r), plus
		// the global raw-r distrib reservoir (always runs, even for blocks
		// later skipped in pass 2).
		long long m_local = 0;
		std::vector<double> reservoir;
		reservoir.reserve((size_t)rsample);
		auto& rng = rngs[(size_t)tid];

		for_each_pair([&](int, int, float r) {
			++m_local;
			float rc = std::min(std::max(r, -0.999999999f), 0.999999999f);
			double z = std::atanh((double)rc);
			if ((int)reservoir.size() < rsample) {
				reservoir.push_back(z);
			} else {
				std::uniform_int_distribution<long long> U(0, m_local - 1);
				long long jr = U(rng);
				if (jr < rsample)
					reservoir[(size_t)jr] = z;
			}

			if (do_distrib)
				consider_distrib(td[(size_t)tid], r, distrib_sample);
		});

		NullFit fit = fit_empirical_null(reservoir, m_local, nsamples, opt.fdr_target, opt.fdr_lambda_cut);

		long long kept_local = 0;
		long long hits_pos = 0, hits_neg = 0;

		// Pass 2: apply the fitted per-tail zstar thresholds and write hits.
		if (fit.ok && m_local >= opt.fdr_min_pairs) {
			for_each_pair([&](int a, int b, float r) {
				float rc = std::min(std::max(r, -0.999999999f), 0.999999999f);
				double z = std::atanh((double)rc);
				double zstar = (z - fit.mu0) / fit.sigma0;

				bool is_pos_hit = std::isfinite(fit.zstar_thresh_pos) && zstar >= fit.zstar_thresh_pos;
				bool is_neg_hit = !is_pos_hit && std::isfinite(fit.zstar_thresh_neg) && zstar <= fit.zstar_thresh_neg;

				if (!is_pos_hit && !is_neg_hit)
					return;

				double pvalue = is_pos_hit ? (1.0 - norm_cdf(zstar)) : norm_cdf(zstar);
				double qvalue = qvalue_for_p(fit, is_pos_hit, pvalue);
				double lfdr = local_fdr_for_z(fit, z);

				hfp << a << "\t" << chroms_scan[a] << "\t" << pos_scan[a] << "\t"
					<< b << "\t" << chroms_scan[b] << "\t" << pos_scan[b] << "\t"
					<< r << "\t" << nsamples << "\t"
					<< z << "\t" << zstar << "\t" << pvalue << "\t" << qvalue << "\t" << lfdr << "\n";

				++kept_local;
				if (is_pos_hit) ++hits_pos; else ++hits_neg;
			});
		}

		sfp << chr1 << "\t" << chr2 << "\t" << m_local << "\t"
			<< (fit.ok && m_local >= opt.fdr_min_pairs ? "ok" : "skipped") << "\t"
			<< fit.mu0 << "\t" << fit.sigma0 << "\t" << fit.lambda << "\t"
			<< fit.pi0_pos << "\t" << fit.pi0_neg << "\t"
			<< hits_pos << "\t" << hits_neg << "\n";

		tested_t[(size_t)tid] += m_local;
		kept_t[(size_t)tid] += kept_local;
	}

	for (int t = 0; t < nthreads; ++t) {
		tested_pairs += tested_t[(size_t)t];
		kept_pairs += kept_t[(size_t)t];
	}

	try {
		merge_part_files(
			out_path, hit_part_paths,
			"wA\tchrA\tposA\twB\tchrB\tposB\tr\tn\tz\tzstar\tpvalue\tqvalue\tlocal_fdr\n"
		);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return false;
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(hit_part_paths[(size_t)t].c_str());

	try {
		merge_part_files(
			summary_path, summary_part_paths,
			"chrA\tchrB\tn_pairs\tstatus\tmu0\tsigma0\tlambda\tpi0_pos\tpi0_neg\tn_hits_pos\tn_hits_neg\n"
		);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return false;
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(summary_part_paths[(size_t)t].c_str());

	if (do_distrib) {
		long long total_tested = 0;
		float global_min = std::numeric_limits<float>::infinity();
		float global_max = -std::numeric_limits<float>::infinity();

		std::vector<float> dsample;
		dsample.reserve((size_t)distrib_sample * (size_t)nthreads);

		for (int t = 0; t < nthreads; ++t) {
			const auto& d = td[(size_t)t];
			if (d.tested_pairs == 0)
				continue;

			total_tested += d.tested_pairs;
			if (d.min_r < global_min) global_min = d.min_r;
			if (d.max_r > global_max) global_max = d.max_r;

			dsample.insert(dsample.end(), d.sample.begin(), d.sample.end());
		}

		if ((int)dsample.size() > distrib_sample) {
			std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
			std::shuffle(dsample.begin(), dsample.end(), rng);
			dsample.resize((size_t)distrib_sample);
		}

		if (!distrib_path.empty()) {
			std::ofstream df(distrib_path);
			if (!df) {
				std::cerr << "Error: cannot write to " << distrib_path << "\n";
				return false;
			}

			df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\tmean_r2\tsd_r2\n";

			if (total_tested == 0 || dsample.empty()) {
				df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
			} else {
				float mean = std::numeric_limits<float>::quiet_NaN();
				float sd = std::numeric_limits<float>::quiet_NaN();
				compute_mean_sd_from_sample(dsample, mean, sd);

				float mean_r2 = std::numeric_limits<float>::quiet_NaN();
				float sd_r2 = std::numeric_limits<float>::quiet_NaN();
				{
					double mu2 = 0.0, m2v = 0.0;
					long long n2 = 0;
					for (float x : dsample) {
						if (!std::isfinite(x)) continue;
						++n2;
						double r2 = (double)x * (double)x;
						double dr2 = r2 - mu2;
						mu2 += dr2 / (double)n2;
						m2v += dr2 * (r2 - mu2);
					}
					if (n2 > 0) mean_r2 = (float)mu2;
					if (n2 > 1) sd_r2 = (float)std::sqrt(m2v / (double)(n2 - 1));
				}

				std::sort(dsample.begin(), dsample.end());

				float p01 = quantile_from_sorted(dsample, 0.01);
				float p05 = quantile_from_sorted(dsample, 0.05);
				float p25 = quantile_from_sorted(dsample, 0.25);
				float med = quantile_from_sorted(dsample, 0.50);
				float p75 = quantile_from_sorted(dsample, 0.75);
				float p95 = quantile_from_sorted(dsample, 0.95);
				float p99 = quantile_from_sorted(dsample, 0.99);

				df << total_tested
					<< "\t" << global_max
					<< "\t" << p99
					<< "\t" << p95
					<< "\t" << p75
					<< "\t" << med
					<< "\t" << p25
					<< "\t" << p05
					<< "\t" << p01
					<< "\t" << global_min
					<< "\t" << mean
					<< "\t" << sd
					<< "\t" << mean_r2
					<< "\t" << sd_r2
					<< "\n";
			}
		}

		if (!reservoir_path.empty() && !dsample.empty()) {
			std::ofstream rf(reservoir_path);
			if (!rf) {
				std::cerr << "Error: cannot write reservoir to " << reservoir_path << "\n";
				return false;
			}
			rf << "r\n";
			for (float x : dsample)
				rf << x << "\n";
		}
	}

	return true;
}

bool scan_markers_write_hits_excl_focus_fdr(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	uint64_t seed,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	auto fn = [](const HiComponentsWeighted& hc, const std::string& a, const std::string& b) {
		return hi_from_components_weighted_excluding(hc, a, b);
	};
	return scan_markers_write_hits_excl_focus_fdr_T(
		X_scan, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, out_path, summary_path, tested_pairs, kept_pairs, seed,
		distrib_path, distrib_sample, distrib_seed, reservoir_path
	);
}

bool scan_markers_write_hits_excl_focus_fdr(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::vector<MarkerFreq>& freqs_scan,
	uint64_t seed,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	auto fn = [](const HiComponentsFreq& hc, const std::string& a, const std::string& b) {
		return hi_from_components_freq_excluding(hc, a, b);
	};
	return scan_markers_write_hits_excl_focus_fdr_T(
		X_scan, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, out_path, summary_path, tested_pairs, kept_pairs, seed,
		distrib_path, distrib_sample, distrib_seed, reservoir_path, &freqs_scan
	);
}

template<typename HC, typename HiExcludeFn>
static bool scan_target_write_hits_excl_focus_T(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HC& hc_full,
	HiExcludeFn hi_excluding,
	const ScanOptions& opt,
	int target_w,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path,
	const std::vector<MarkerFreq>* freqs_scan = nullptr
) {
	const bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;
	const float denom = 1.0f / (float)(nsamples - 1);

	tested_pairs = 0;
	kept_pairs = 0;

	const std::string& tchr = chroms_scan[target_w];
	const int tpos = pos_scan[target_w];

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
	nthreads = opt.threads;
	if (nthreads < 1) nthreads = 1;
	omp_set_num_threads(nthreads);
	#endif

	if (nthreads > (int)chr_order_scan.size())
		nthreads = (int)chr_order_scan.size();
	if (nthreads < 1)
		nthreads = 1;

	std::vector<std::string> part_paths((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t) {
		part_paths[(size_t)t] = out_path + ".part." + std::to_string(t);
		std::remove(part_paths[(size_t)t].c_str());
	}

	std::vector<long long> tested_t((size_t)nthreads, 0);
	std::vector<long long> kept_t((size_t)nthreads, 0);

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	const int C = (int)chr_order_scan.size();

	#ifdef ADMIXLD_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int c = 0; c < C; ++c) {
		int tid = 0;
		#ifdef ADMIXLD_HAS_OPENMP
		tid = omp_get_thread_num();
		if (tid >= nthreads) tid = tid % nthreads;
		#endif

		const auto& chr = chr_order_scan[c];

		if (!opt.intra && chr == tchr)
			continue;

		const auto& idx = windows_by_chr_scan.at(chr);
		const int m = (int)idx.size();
		if (m == 0)
			continue;

		std::ofstream ofp(part_paths[(size_t)tid], std::ios::out | std::ios::app);
		if (!ofp) {
			#ifdef ADMIXLD_HAS_OPENMP
			#pragma omp critical
			#endif
			{
				std::cerr << "Error: cannot write to " << part_paths[(size_t)tid] << "\n";
			}
			continue;
		}

		const bool same_chr = (chr == tchr);

		Eigen::VectorXf h;
		if (same_chr)
			h = hi_excluding(hc_full, tchr, std::string(""));
		else
			h = hi_excluding(hc_full, tchr, chr);

		auto local_consider = [&](float r) {
			if (!do_distrib)
				return;
			consider_distrib(td[(size_t)tid], r, distrib_sample);
		};

		if (same_chr) {
			int n_valid = 0;
			Eigen::MatrixXf Zc = residualize_and_zscore_subset(X_scan, h, idx, n_valid, freqs_scan);

			int j_target = -1;
			for (int j = 0; j < m; ++j) {
				if (idx[j] == target_w) {
					j_target = j;
					break;
				}
			}
			if (j_target < 0)
				continue;

			Eigen::VectorXf v = Zc.col(j_target);

			Eigen::MatrixXf B(nsamples, tile_size);
			Eigen::RowVectorXf R(tile_size);

			long long tested_local = 0;
			long long kept_local = 0;

			for (int j0 = 0; j0 < m; j0 += tile_size) {
				const int b2 = std::min(tile_size, m - j0);

				for (int k = 0; k < b2; ++k)
					B.col(k) = Zc.col(j0 + k);

				R.head(b2).noalias() = v.transpose() * B.leftCols(b2);
				R.head(b2) *= denom;

				for (int ib = 0; ib < b2; ++ib) {
					const int b = idx[j0 + ib];

					if (b == target_w)
						continue;

					if (opt.max_dist >= 0 || opt.min_dist >= 0) {
						int dpos = pos_scan[b] - tpos;
						if (dpos < 0) dpos = -dpos;
						if (opt.max_dist >= 0 && dpos > opt.max_dist)
							continue;
						if (opt.min_dist >= 0 && dpos < opt.min_dist)
							continue;
					}

					const float r = R(ib);
					if (!std::isfinite(r))
						continue;

					++tested_local;
					local_consider(r);

					if (keep_hit_r(r, opt)) {
						ofp << target_w << "\t" << chroms_scan[target_w] << "\t" << pos_scan[target_w] << "\t"
							<< b << "\t" << chroms_scan[b] << "\t" << pos_scan[b] << "\t"
							<< r << "\t" << nsamples << "\n";
						++kept_local;
					}
				}
			}

			tested_t[(size_t)tid] += tested_local;
			kept_t[(size_t)tid] += kept_local;
		} else {
			std::vector<int> tvec(1);
			tvec[0] = target_w;

			int n_valid_t = 0;
			Eigen::MatrixXf Zt = residualize_and_zscore_subset(X_scan, h, tvec, n_valid_t, freqs_scan);
			Eigen::VectorXf v = Zt.col(0);

			int n_valid_b = 0;
			Eigen::MatrixXf Zb = residualize_and_zscore_subset(X_scan, h, idx, n_valid_b, freqs_scan);

			Eigen::MatrixXf B(nsamples, tile_size);
			Eigen::RowVectorXf R(tile_size);

			long long tested_local = 0;
			long long kept_local = 0;

			for (int j0 = 0; j0 < m; j0 += tile_size) {
				const int b2 = std::min(tile_size, m - j0);

				for (int k = 0; k < b2; ++k)
					B.col(k) = Zb.col(j0 + k);

				R.head(b2).noalias() = v.transpose() * B.leftCols(b2);
				R.head(b2) *= denom;

				for (int ib = 0; ib < b2; ++ib) {
					const int b = idx[j0 + ib];

					const float r = R(ib);
					if (!std::isfinite(r))
						continue;

					++tested_local;
					local_consider(r);

					if (keep_hit_r(r, opt)) {
						ofp << target_w << "\t" << chroms_scan[target_w] << "\t" << pos_scan[target_w] << "\t"
							<< b << "\t" << chroms_scan[b] << "\t" << pos_scan[b] << "\t"
							<< r << "\t" << nsamples << "\n";
						++kept_local;
					}
				}
			}

			tested_t[(size_t)tid] += tested_local;
			kept_t[(size_t)tid] += kept_local;
		}
	}

	for (int t = 0; t < nthreads; ++t) {
		tested_pairs += tested_t[(size_t)t];
		kept_pairs += kept_t[(size_t)t];
	}

	{
		std::ofstream of(out_path);
		if (!of) {
			std::cerr << "Error: cannot write to " << out_path << "\n";
			return false;
		}

		of << "wA\tchrA\tposA\twB\tchrB\tposB\tr\tn\n";

		for (int t = 0; t < nthreads; ++t) {
			std::ifstream pf(part_paths[(size_t)t]);
			if (!pf)
				continue;

			std::string line;
			while (std::getline(pf, line))
				of << line << "\n";
		}
	}

	for (int t = 0; t < nthreads; ++t)
		std::remove(part_paths[(size_t)t].c_str());

	if (do_distrib) {
		long long total_tested = 0;
		float global_min = std::numeric_limits<float>::infinity();
		float global_max = -std::numeric_limits<float>::infinity();

		std::vector<float> dsample;
		dsample.reserve((size_t)distrib_sample * (size_t)nthreads);

		for (int t = 0; t < nthreads; ++t) {
			const auto& d = td[(size_t)t];
			if (d.tested_pairs == 0)
				continue;

			total_tested += d.tested_pairs;
			if (d.min_r < global_min) global_min = d.min_r;
			if (d.max_r > global_max) global_max = d.max_r;

			dsample.insert(dsample.end(), d.sample.begin(), d.sample.end());
		}

		if (!distrib_path.empty()) {
			std::ofstream df(distrib_path);
			if (!df) {
				std::cerr << "Error: cannot write to " << distrib_path << "\n";
				return false;
			}

			df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\tmean_r2\tsd_r2\n";

			if (total_tested == 0 || dsample.empty()) {
				df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
			} else {
				if ((int)dsample.size() > distrib_sample) {
					std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
					std::shuffle(dsample.begin(), dsample.end(), rng);
					dsample.resize((size_t)distrib_sample);
				}

				float mean = std::numeric_limits<float>::quiet_NaN();
				float sd = std::numeric_limits<float>::quiet_NaN();
				compute_mean_sd_from_sample(dsample, mean, sd);

				float mean_r2 = std::numeric_limits<float>::quiet_NaN();
				float sd_r2 = std::numeric_limits<float>::quiet_NaN();
				{
					double mu2 = 0.0, m2 = 0.0;
					long long n2 = 0;
					for (float x : dsample) {
						if (!std::isfinite(x)) continue;
						++n2;
						double r2 = (double)x * (double)x;
						double dr2 = r2 - mu2;
						mu2 += dr2 / (double)n2;
						m2 += dr2 * (r2 - mu2);
					}
					if (n2 > 0) mean_r2 = (float)mu2;
					if (n2 > 1) sd_r2 = (float)std::sqrt(m2 / (double)(n2 - 1));
				}

				std::sort(dsample.begin(), dsample.end());

				float p01 = quantile_from_sorted(dsample, 0.01);
				float p05 = quantile_from_sorted(dsample, 0.05);
				float p25 = quantile_from_sorted(dsample, 0.25);
				float med = quantile_from_sorted(dsample, 0.50);
				float p75 = quantile_from_sorted(dsample, 0.75);
				float p95 = quantile_from_sorted(dsample, 0.95);
				float p99 = quantile_from_sorted(dsample, 0.99);

				df << total_tested
					<< "\t" << global_max
					<< "\t" << p99
					<< "\t" << p95
					<< "\t" << p75
					<< "\t" << med
					<< "\t" << p25
					<< "\t" << p05
					<< "\t" << p01
					<< "\t" << global_min
					<< "\t" << mean
					<< "\t" << sd
					<< "\t" << mean_r2
					<< "\t" << sd_r2
					<< "\n";
			}
		} else if (!dsample.empty() && (int)dsample.size() > distrib_sample) {
			std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
			std::shuffle(dsample.begin(), dsample.end(), rng);
			dsample.resize((size_t)distrib_sample);
		}

		if (!reservoir_path.empty() && !dsample.empty()) {
			std::ofstream rf(reservoir_path);
			if (!rf) {
				std::cerr << "Error: cannot write reservoir to " << reservoir_path << "\n";
				return false;
			}
			rf << "r\n";
			for (float x : dsample)
				rf << x << "\n";
		}
	}

	return true;
}

bool scan_target_write_hits_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	int target_w,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	auto fn = [](const HiComponentsWeighted& hc, const std::string& a, const std::string& b) {
		return hi_from_components_weighted_excluding(hc, a, b);
	};
	return scan_target_write_hits_excl_focus_T(
		X_scan, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, target_w, out_path, tested_pairs, kept_pairs,
		distrib_path, distrib_sample, distrib_seed, reservoir_path
	);
}

bool scan_target_write_hits_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const ScanOptions& opt,
	int target_w,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::vector<MarkerFreq>& freqs_scan,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	auto fn = [](const HiComponentsFreq& hc, const std::string& a, const std::string& b) {
		return hi_from_components_freq_excluding(hc, a, b);
	};
	return scan_target_write_hits_excl_focus_T(
		X_scan, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, target_w, out_path, tested_pairs, kept_pairs,
		distrib_path, distrib_sample, distrib_seed, reservoir_path, &freqs_scan
	);
}

template<typename HC, typename HiExcludeFn>
static bool scan_vector_vs_windows_write_hits_excl_focus_T(
	const Eigen::MatrixXf& X_scan,
	const Eigen::VectorXf& v_raw,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HC& hc_full,
	HiExcludeFn hi_excluding,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path,
	const std::vector<MarkerFreq>* freqs_scan = nullptr
) {
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	const int nsamples = opt.nsamples;
	if (v_raw.size() != nsamples) {
		std::cerr << "Error: scan_vector_vs_windows_write_hits_excl_focus: v_raw size != nsamples\n";
		return false;
	}

	tested_pairs = 0;
	kept_pairs = 0;

	const bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	const float denom = 1.0f / (float)(nsamples - 1);

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
	nthreads = opt.threads;
	if (nthreads < 1) nthreads = 1;
	omp_set_num_threads(nthreads);
	#endif

	if (nthreads > (int)chr_order_scan.size())
		nthreads = (int)chr_order_scan.size();
	if (nthreads < 1)
		nthreads = 1;

	std::vector<std::string> part_paths((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t) {
		part_paths[(size_t)t] = out_path + ".part." + std::to_string(t);
		std::remove(part_paths[(size_t)t].c_str());
	}

	std::vector<long long> tested_t((size_t)nthreads, 0);
	std::vector<long long> kept_t((size_t)nthreads, 0);

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	const int C = (int)chr_order_scan.size();

	#ifdef ADMIXLD_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int c = 0; c < C; ++c) {
		int tid = 0;
		#ifdef ADMIXLD_HAS_OPENMP
		tid = omp_get_thread_num();
		if (tid >= nthreads) tid = tid % nthreads;
		#endif

		const auto& chr = chr_order_scan[c];
		const auto& idx = windows_by_chr_scan.at(chr);
		const int m = (int)idx.size();
		if (m == 0)
			continue;

		Eigen::VectorXf h = hi_excluding(hc_full, chr, std::string(""));

		int v_valid = 0;
		Eigen::VectorXf vZ;
		try {
			vZ = residualize_and_zscore_vector(v_raw, h, v_valid);
		} catch (...) {
			continue;
		}

		int n_valid = 0;
		Eigen::MatrixXf Zc;
		try {
			Zc = residualize_and_zscore_subset(X_scan, h, idx, n_valid, freqs_scan);
		} catch (...) {
			continue;
		}

		std::ofstream ofp(part_paths[(size_t)tid], std::ios::out | std::ios::app);
		if (!ofp) {
			#ifdef ADMIXLD_HAS_OPENMP
			#pragma omp critical
			#endif
			{
				std::cerr << "Error: cannot write to " << part_paths[(size_t)tid] << "\n";
			}
			continue;
		}

		long long tested_local = 0;
		long long kept_local = 0;

		for (int j = 0; j < m; ++j) {
			const int w = idx[j];

			double dot = 0.0;
			for (int i = 0; i < nsamples; ++i) {
				float zi = Zc(i, j);
				float vi = vZ(i);
				if (!std::isnan(zi) && !std::isnan(vi))
					dot += (double)zi * (double)vi;
			}

			float r = (float)(dot * denom);
			if (!std::isfinite(r))
				continue;

			++tested_local;
			if (do_distrib)
				consider_distrib(td[(size_t)tid], r, distrib_sample);

			if (keep_hit_r(r, opt)) {
				ofp << "sample_haplo"
					<< "\t" << w
					<< "\t" << chroms_scan[w]
					<< "\t" << pos_scan[w]
					<< "\t" << r
					<< "\t" << nsamples
					<< "\n";
				++kept_local;
			}
		}

		tested_t[(size_t)tid] += tested_local;
		kept_t[(size_t)tid] += kept_local;
	}

	for (int t = 0; t < nthreads; ++t) {
		tested_pairs += tested_t[(size_t)t];
		kept_pairs += kept_t[(size_t)t];
	}

	{
		std::ofstream of(out_path);
		if (!of) {
			std::cerr << "Error: cannot write to " << out_path << "\n";
			return false;
		}

		of << "tag\twB\tchrB\tposB\tr\tn\n";

		for (int t = 0; t < nthreads; ++t) {
			std::ifstream pf(part_paths[(size_t)t]);
			if (!pf)
				continue;

			std::string line;
			while (std::getline(pf, line))
				of << line << "\n";
		}
	}

	for (int t = 0; t < nthreads; ++t)
		std::remove(part_paths[(size_t)t].c_str());

	if (do_distrib) {
		long long total_tested = 0;
		float global_min = std::numeric_limits<float>::infinity();
		float global_max = -std::numeric_limits<float>::infinity();

		std::vector<float> dsample;
		dsample.reserve((size_t)distrib_sample * (size_t)nthreads);

		for (int t = 0; t < nthreads; ++t) {
			const auto& d = td[(size_t)t];
			if (d.tested_pairs == 0)
				continue;

			total_tested += d.tested_pairs;
			if (d.min_r < global_min) global_min = d.min_r;
			if (d.max_r > global_max) global_max = d.max_r;

			dsample.insert(dsample.end(), d.sample.begin(), d.sample.end());
		}

		if (!distrib_path.empty()) {
			std::ofstream df(distrib_path);
			if (!df) {
				std::cerr << "Error: cannot write distrib summary to " << distrib_path << "\n";
				return false;
			}

			df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\tmean_r2\tsd_r2\n";

			if (total_tested == 0 || dsample.empty()) {
				df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
			} else {
				if ((int)dsample.size() > distrib_sample) {
					std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
					std::shuffle(dsample.begin(), dsample.end(), rng);
					dsample.resize((size_t)distrib_sample);
				}

				float mean = std::numeric_limits<float>::quiet_NaN();
				float sd = std::numeric_limits<float>::quiet_NaN();
				compute_mean_sd_from_sample(dsample, mean, sd);

				float mean_r2 = std::numeric_limits<float>::quiet_NaN();
				float sd_r2 = std::numeric_limits<float>::quiet_NaN();
				{
					double mu2 = 0.0, m2 = 0.0;
					long long n2 = 0;
					for (float x : dsample) {
						if (!std::isfinite(x)) continue;
						++n2;
						double r2 = (double)x * (double)x;
						double dr2 = r2 - mu2;
						mu2 += dr2 / (double)n2;
						m2 += dr2 * (r2 - mu2);
					}
					if (n2 > 0) mean_r2 = (float)mu2;
					if (n2 > 1) sd_r2 = (float)std::sqrt(m2 / (double)(n2 - 1));
				}

				std::sort(dsample.begin(), dsample.end());

				float p01 = quantile_from_sorted(dsample, 0.01);
				float p05 = quantile_from_sorted(dsample, 0.05);
				float p25 = quantile_from_sorted(dsample, 0.25);
				float med = quantile_from_sorted(dsample, 0.50);
				float p75 = quantile_from_sorted(dsample, 0.75);
				float p95 = quantile_from_sorted(dsample, 0.95);
				float p99 = quantile_from_sorted(dsample, 0.99);

				df << total_tested
					<< "\t" << global_max
					<< "\t" << p99
					<< "\t" << p95
					<< "\t" << p75
					<< "\t" << med
					<< "\t" << p25
					<< "\t" << p05
					<< "\t" << p01
					<< "\t" << global_min
					<< "\t" << mean
					<< "\t" << sd
					<< "\t" << mean_r2
					<< "\t" << sd_r2
					<< "\n";
			}
		} else if (!dsample.empty() && (int)dsample.size() > distrib_sample) {
			std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
			std::shuffle(dsample.begin(), dsample.end(), rng);
			dsample.resize((size_t)distrib_sample);
		}

		if (!reservoir_path.empty() && !dsample.empty()) {
			std::ofstream rf(reservoir_path);
			if (!rf) {
				std::cerr << "Error: cannot write reservoir to " << reservoir_path << "\n";
				return false;
			}
			rf << "r\n";
			for (float x : dsample)
				rf << x << "\n";
		}
	}

	return true;
}

bool scan_vector_vs_windows_write_hits_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const Eigen::VectorXf& v_raw,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	auto fn = [](const HiComponentsWeighted& hc, const std::string& a, const std::string& b) {
		return hi_from_components_weighted_excluding(hc, a, b);
	};
	return scan_vector_vs_windows_write_hits_excl_focus_T(
		X_scan, v_raw, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, out_path, tested_pairs, kept_pairs,
		distrib_path, distrib_sample, distrib_seed, reservoir_path
	);
}

bool scan_vector_vs_windows_write_hits_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const Eigen::VectorXf& v_raw,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::vector<MarkerFreq>& freqs_scan,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	auto fn = [](const HiComponentsFreq& hc, const std::string& a, const std::string& b) {
		return hi_from_components_freq_excluding(hc, a, b);
	};
	return scan_vector_vs_windows_write_hits_excl_focus_T(
		X_scan, v_raw, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, out_path, tested_pairs, kept_pairs,
		distrib_path, distrib_sample, distrib_seed, reservoir_path, &freqs_scan
	);
}

// Shared calibrate+hit-call pass for a single-vector-vs-one-chromosome FDR
// block, given that block's precomputed (b, r) pairs (bs/rs, index-aligned).
// Used by both the LOCO target and LOCO sample-haplo FDR scans below.
template <typename RowWriterFn>
static void run_vector_block_fdr(
	const std::vector<int>& bs,
	const std::vector<float>& rs,
	int nsamples,
	const ScanOptions& opt,
	int rsample,
	std::mt19937_64& rng,
	bool do_distrib,
	ThreadDistrib* td_ptr,
	int distrib_sample,
	std::ofstream& hfp,
	std::ofstream& sfp,
	const std::string& chr_label,
	long long& tested_out,
	long long& kept_out,
	RowWriterFn write_row
) {
	const long long m_local = (long long)bs.size();

	std::vector<double> reservoir;
	reservoir.reserve((size_t)std::min<long long>(rsample, m_local));

	for (long long k = 0; k < m_local; ++k) {
		const float r = rs[(size_t)k];
		const float rc = std::min(std::max(r, -0.999999999f), 0.999999999f);
		const double z = std::atanh((double)rc);

		if ((int)reservoir.size() < rsample) {
			reservoir.push_back(z);
		} else {
			std::uniform_int_distribution<long long> U(0, k);
			long long jr = U(rng);
			if (jr < rsample)
				reservoir[(size_t)jr] = z;
		}

		if (do_distrib)
			consider_distrib(*td_ptr, r, distrib_sample);
	}

	NullFit fit = fit_empirical_null(reservoir, m_local, nsamples, opt.fdr_target, opt.fdr_lambda_cut);

	long long kept_local = 0;
	long long hits_pos = 0, hits_neg = 0;

	if (fit.ok && m_local >= opt.fdr_min_pairs) {
		for (long long k = 0; k < m_local; ++k) {
			const int b = bs[(size_t)k];
			const float r = rs[(size_t)k];
			const float rc = std::min(std::max(r, -0.999999999f), 0.999999999f);
			const double z = std::atanh((double)rc);
			const double zstar = (z - fit.mu0) / fit.sigma0;

			bool is_pos_hit = std::isfinite(fit.zstar_thresh_pos) && zstar >= fit.zstar_thresh_pos;
			bool is_neg_hit = !is_pos_hit && std::isfinite(fit.zstar_thresh_neg) && zstar <= fit.zstar_thresh_neg;

			if (!is_pos_hit && !is_neg_hit)
				continue;

			double pvalue = is_pos_hit ? (1.0 - norm_cdf(zstar)) : norm_cdf(zstar);
			double qvalue = qvalue_for_p(fit, is_pos_hit, pvalue);
			double lfdr = local_fdr_for_z(fit, z);

			write_row(hfp, b, r, z, zstar, pvalue, qvalue, lfdr);

			++kept_local;
			if (is_pos_hit) ++hits_pos; else ++hits_neg;
		}
	}

	sfp << chr_label << "\t" << m_local << "\t"
		<< (fit.ok && m_local >= opt.fdr_min_pairs ? "ok" : "skipped") << "\t"
		<< fit.mu0 << "\t" << fit.sigma0 << "\t" << fit.lambda << "\t"
		<< fit.pi0_pos << "\t" << fit.pi0_neg << "\t"
		<< hits_pos << "\t" << hits_neg << "\n";

	tested_out += m_local;
	kept_out += kept_local;
}

// Distance-stratified intrachromosomal variant: HI excludes the target's
// own chromosome entirely, computed once and reused across every bin.
template<typename HC, typename HiExcludeFn>
static bool scan_target_write_hits_excl_focus_fdr_intra_T(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const HC& hc_full,
	HiExcludeFn hi_excluding,
	const ScanOptions& opt,
	int target_w,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	uint64_t seed,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path,
	const std::vector<MarkerFreq>* freqs_scan
) {
	const int nsamples = opt.nsamples;
	const int rsample = opt.fdr_sample > 0 ? opt.fdr_sample : 200000;
	const int n_bins = opt.fdr_intra_bins > 0 ? opt.fdr_intra_bins : 8;
	const long long real_floor = opt.min_dist >= 0 ? (long long)opt.min_dist : 0;

	tested_pairs = 0;
	kept_pairs = 0;

	const std::string& tchr = chroms_scan[target_w];
	const int tpos = pos_scan[target_w];
	const auto& idx = windows_by_chr_scan.at(tchr);
	const int m = (int)idx.size();

	const bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	long long effective_min = std::max(real_floor, 1LL);
	long long effective_max = opt.max_dist >= 0 ? (long long)opt.max_dist : 0;
	if (effective_max <= 0 && m >= 2)
		effective_max = (long long)pos_scan[idx.back()] - (long long)pos_scan[idx.front()];
	if (effective_max <= effective_min)
		effective_max = effective_min + 1;

	std::vector<long long> edges((size_t)n_bins + 1);
	{
		const double log_lo = std::log((double)effective_min);
		const double log_hi = std::log((double)effective_max);
		for (int i = 0; i <= n_bins; ++i) {
			const double t = (double)i / (double)n_bins;
			edges[(size_t)i] = (long long)std::llround(std::exp(log_lo + t * (log_hi - log_lo)));
		}
		for (int i = 1; i <= n_bins; ++i)
			if (edges[(size_t)i] <= edges[(size_t)(i - 1)])
				edges[(size_t)i] = edges[(size_t)(i - 1)] + 1;
	}

	std::cout << "[scan-fdr-loco-target-intra] target_fdr=" << opt.fdr_target
		  << " fdr_sample=" << rsample
		  << " fdr_min_pairs=" << opt.fdr_min_pairs
		  << " fdr_lambda=" << opt.fdr_lambda_cut
		  << " bins=" << n_bins
		  << " min_dist=" << real_floor
		  << " max_dist=" << effective_max << "\n";

	Eigen::VectorXf h = hi_excluding(hc_full, tchr, std::string(""));

	std::vector<int> tvec(1, target_w);
	int n_valid_t = 0;
	Eigen::MatrixXf Zt = residualize_and_zscore_subset(X_scan, h, tvec, n_valid_t, freqs_scan);
	Eigen::VectorXf v = Zt.col(0);

	int n_valid_b = 0;
	Eigen::MatrixXf Zc = residualize_and_zscore_subset(X_scan, h, idx, n_valid_b, freqs_scan);

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
	nthreads = opt.threads;
	if (nthreads < 1) nthreads = 1;
	omp_set_num_threads(nthreads);
	#endif
	if (nthreads > n_bins)
		nthreads = std::max(1, n_bins);

	std::vector<std::string> hit_part_paths((size_t)nthreads);
	std::vector<std::string> summary_part_paths((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t) {
		hit_part_paths[(size_t)t] = out_path + ".part." + std::to_string(t);
		summary_part_paths[(size_t)t] = summary_path + ".part." + std::to_string(t);
		std::remove(hit_part_paths[(size_t)t].c_str());
		std::remove(summary_part_paths[(size_t)t].c_str());
	}

	std::vector<long long> tested_t((size_t)nthreads, 0);
	std::vector<long long> kept_t((size_t)nthreads, 0);
	std::vector<std::mt19937_64> rngs((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t)
		rngs[(size_t)t] = std::mt19937_64(seed + (uint64_t)t);

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	#ifdef ADMIXLD_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int bidx = 0; bidx < n_bins; ++bidx) {
		int tid = 0;
		#ifdef ADMIXLD_HAS_OPENMP
		tid = omp_get_thread_num();
		if (tid >= nthreads) tid = tid % nthreads;
		#endif

		const long long bin_lo = (bidx == 0) ? real_floor : edges[(size_t)bidx];
		const long long bin_hi = edges[(size_t)bidx + 1];
		if (bin_lo >= bin_hi || m < 1)
			continue;

		std::ofstream hfp(hit_part_paths[(size_t)tid], std::ios::out | std::ios::app);
		std::ofstream sfp(summary_part_paths[(size_t)tid], std::ios::out | std::ios::app);
		if (!hfp || !sfp) {
			#ifdef ADMIXLD_HAS_OPENMP
			#pragma omp critical
			#endif
			{
				std::cerr << "Error: cannot write part files for thread " << tid << "\n";
			}
			continue;
		}

		std::vector<int> bs;
		std::vector<float> rs;
		bs.reserve((size_t)m);
		rs.reserve((size_t)m);

		for (int k = 0; k < m; ++k) {
			const int b = idx[k];
			if (b == target_w)
				continue;

			long long d = (long long)pos_scan[b] - (long long)tpos;
			if (d < 0) d = -d;
			if (d < bin_lo || d >= bin_hi)
				continue;

			double dot = 0.0;
			for (int i = 0; i < nsamples; ++i)
				dot += (double)v(i) * (double)Zc(i, k);

			const float r = (float)(dot / (double)(nsamples - 1));
			if (!std::isfinite(r))
				continue;

			bs.push_back(b);
			rs.push_back(r);
		}

		auto write_row = [&](std::ofstream& out, int b, float r, double z, double zstar, double pvalue, double qvalue, double lfdr) {
			out << target_w << "\t" << tchr << "\t" << tpos << "\t"
				<< b << "\t" << chroms_scan[b] << "\t" << pos_scan[b] << "\t"
				<< r << "\t" << nsamples << "\t"
				<< z << "\t" << zstar << "\t" << pvalue << "\t" << qvalue << "\t" << lfdr << "\n";
		};

		long long tested_local = 0;
		long long kept_local = 0;

		const std::string bin_label = tchr + "\t" + std::to_string(bin_lo) + "\t" + std::to_string(bin_hi);

		run_vector_block_fdr(
			bs, rs, nsamples, opt, rsample, rngs[(size_t)tid],
			do_distrib, do_distrib ? &td[(size_t)tid] : nullptr, distrib_sample,
			hfp, sfp, bin_label, tested_local, kept_local, write_row
		);

		tested_t[(size_t)tid] += tested_local;
		kept_t[(size_t)tid] += kept_local;
	}

	for (int t = 0; t < nthreads; ++t) {
		tested_pairs += tested_t[(size_t)t];
		kept_pairs += kept_t[(size_t)t];
	}

	try {
		merge_part_files(
			out_path, hit_part_paths,
			"wA\tchrA\tposA\twB\tchrB\tposB\tr\tn\tz\tzstar\tpvalue\tqvalue\tlocal_fdr\n"
		);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return false;
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(hit_part_paths[(size_t)t].c_str());

	try {
		merge_part_files(
			summary_path, summary_part_paths,
			"chr\tdist_lo\tdist_hi\tn_pairs\tstatus\tmu0\tsigma0\tlambda\tpi0_pos\tpi0_neg\tn_hits_pos\tn_hits_neg\n"
		);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return false;
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(summary_part_paths[(size_t)t].c_str());

	if (do_distrib) {
		long long total_tested = 0;
		float global_min = std::numeric_limits<float>::infinity();
		float global_max = -std::numeric_limits<float>::infinity();

		std::vector<float> dsample;
		dsample.reserve((size_t)distrib_sample * (size_t)nthreads);

		for (int t = 0; t < nthreads; ++t) {
			const auto& d = td[(size_t)t];
			if (d.tested_pairs == 0)
				continue;

			total_tested += d.tested_pairs;
			if (d.min_r < global_min) global_min = d.min_r;
			if (d.max_r > global_max) global_max = d.max_r;

			dsample.insert(dsample.end(), d.sample.begin(), d.sample.end());
		}

		if ((int)dsample.size() > distrib_sample) {
			std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
			std::shuffle(dsample.begin(), dsample.end(), rng);
			dsample.resize((size_t)distrib_sample);
		}

		if (!distrib_path.empty()) {
			std::ofstream df(distrib_path);
			if (!df) {
				std::cerr << "Error: cannot write to " << distrib_path << "\n";
				return false;
			}

			df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\tmean_r2\tsd_r2\n";

			if (total_tested == 0 || dsample.empty()) {
				df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
			} else {
				float mean = std::numeric_limits<float>::quiet_NaN();
				float sd = std::numeric_limits<float>::quiet_NaN();
				compute_mean_sd_from_sample(dsample, mean, sd);

				float mean_r2 = std::numeric_limits<float>::quiet_NaN();
				float sd_r2 = std::numeric_limits<float>::quiet_NaN();
				{
					double mu2 = 0.0, m2v = 0.0;
					long long n2 = 0;
					for (float x : dsample) {
						if (!std::isfinite(x)) continue;
						++n2;
						double r2 = (double)x * (double)x;
						double dr2 = r2 - mu2;
						mu2 += dr2 / (double)n2;
						m2v += dr2 * (r2 - mu2);
					}
					if (n2 > 0) mean_r2 = (float)mu2;
					if (n2 > 1) sd_r2 = (float)std::sqrt(m2v / (double)(n2 - 1));
				}

				std::sort(dsample.begin(), dsample.end());

				float p01 = quantile_from_sorted(dsample, 0.01);
				float p05 = quantile_from_sorted(dsample, 0.05);
				float p25 = quantile_from_sorted(dsample, 0.25);
				float med = quantile_from_sorted(dsample, 0.50);
				float p75 = quantile_from_sorted(dsample, 0.75);
				float p95 = quantile_from_sorted(dsample, 0.95);
				float p99 = quantile_from_sorted(dsample, 0.99);

				df << total_tested
					<< "\t" << global_max
					<< "\t" << p99
					<< "\t" << p95
					<< "\t" << p75
					<< "\t" << med
					<< "\t" << p25
					<< "\t" << p05
					<< "\t" << p01
					<< "\t" << global_min
					<< "\t" << mean
					<< "\t" << sd
					<< "\t" << mean_r2
					<< "\t" << sd_r2
					<< "\n";
			}
		}

		if (!reservoir_path.empty() && !dsample.empty()) {
			std::ofstream rf(reservoir_path);
			if (!rf) {
				std::cerr << "Error: cannot write reservoir to " << reservoir_path << "\n";
				return false;
			}
			rf << "r\n";
			for (float x : dsample)
				rf << x << "\n";
		}
	}

	return true;
}

// Per chromosome block, runs the calibrate+hit-call two-pass procedure
// against that block's LOCO-residualized target vector / marker matrix.
// Delegates to the intra variant above when opt.intra is set.
template<typename HC, typename HiExcludeFn>
static bool scan_target_write_hits_excl_focus_fdr_T(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HC& hc_full,
	HiExcludeFn hi_excluding,
	const ScanOptions& opt,
	int target_w,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	uint64_t seed,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path,
	const std::vector<MarkerFreq>* freqs_scan = nullptr
) {
	// --intra is additive: the interchromosomal portion below always runs;
	// when opt.intra is set its output is merged with the intrachromosomal
	// distance-binned portion at the end of this function.
	const std::string inter_hits = opt.intra ? out_path + ".fdrpart.inter" : out_path;
	const std::string inter_summary = opt.intra ? summary_path + ".fdrpart.inter" : summary_path;

	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;
	const float denom = 1.0f / (float)(nsamples - 1);
	const int rsample = opt.fdr_sample > 0 ? opt.fdr_sample : 200000;

	std::cout << "[scan-fdr-loco-target] target_fdr=" << opt.fdr_target
		  << " fdr_sample=" << rsample
		  << " fdr_min_pairs=" << opt.fdr_min_pairs
		  << " fdr_lambda=" << opt.fdr_lambda_cut << "\n";

	tested_pairs = 0;
	kept_pairs = 0;

	const std::string& tchr = chroms_scan[target_w];
	const int tpos = pos_scan[target_w];

	const bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
	nthreads = opt.threads;
	if (nthreads < 1) nthreads = 1;
	omp_set_num_threads(nthreads);
	#endif
	if (nthreads > (int)chr_order_scan.size())
		nthreads = std::max(1, (int)chr_order_scan.size());

	std::vector<std::string> hit_part_paths((size_t)nthreads);
	std::vector<std::string> summary_part_paths((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t) {
		hit_part_paths[(size_t)t] = inter_hits + ".part." + std::to_string(t);
		summary_part_paths[(size_t)t] = inter_summary + ".part." + std::to_string(t);
		std::remove(hit_part_paths[(size_t)t].c_str());
		std::remove(summary_part_paths[(size_t)t].c_str());
	}

	std::vector<long long> tested_t((size_t)nthreads, 0);
	std::vector<long long> kept_t((size_t)nthreads, 0);
	std::vector<std::mt19937_64> rngs((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t)
		rngs[(size_t)t] = std::mt19937_64(seed + (uint64_t)t);

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	const int C = (int)chr_order_scan.size();

	#ifdef ADMIXLD_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int c = 0; c < C; ++c) {
		int tid = 0;
		#ifdef ADMIXLD_HAS_OPENMP
		tid = omp_get_thread_num();
		if (tid >= nthreads) tid = tid % nthreads;
		#endif

		const auto& chr = chr_order_scan[c];
		if (chr == tchr)
			continue;

		const auto& idx = windows_by_chr_scan.at(chr);
		const int m = (int)idx.size();
		if (m == 0)
			continue;

		std::ofstream hfp(hit_part_paths[(size_t)tid], std::ios::out | std::ios::app);
		std::ofstream sfp(summary_part_paths[(size_t)tid], std::ios::out | std::ios::app);
		if (!hfp || !sfp) {
			#ifdef ADMIXLD_HAS_OPENMP
			#pragma omp critical
			#endif
			{
				std::cerr << "Error: cannot write part files for thread " << tid << "\n";
			}
			continue;
		}

		// Block-specific LOCO HI (excludes both the target's own chromosome
		// and this block's chromosome); Zt/Zb don't change between passes.
		Eigen::VectorXf h = hi_excluding(hc_full, tchr, chr);

		std::vector<int> tvec(1, target_w);
		int n_valid_t = 0;
		Eigen::MatrixXf Zt = residualize_and_zscore_subset(X_scan, h, tvec, n_valid_t, freqs_scan);
		Eigen::VectorXf v = Zt.col(0);

		int n_valid_b = 0;
		Eigen::MatrixXf Zb = residualize_and_zscore_subset(X_scan, h, idx, n_valid_b, freqs_scan);

		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::RowVectorXf R(tile_size);

		std::vector<int> bs;
		std::vector<float> rs;
		bs.reserve((size_t)m);
		rs.reserve((size_t)m);

		for (int j0 = 0; j0 < m; j0 += tile_size) {
			const int b2 = std::min(tile_size, m - j0);
			for (int k = 0; k < b2; ++k)
				B.col(k) = Zb.col(j0 + k);

			R.head(b2).noalias() = v.transpose() * B.leftCols(b2);
			R.head(b2) *= denom;

			for (int ib = 0; ib < b2; ++ib) {
				const int b = idx[j0 + ib];
				const float r = R(ib);
				if (!std::isfinite(r))
					continue;
				bs.push_back(b);
				rs.push_back(r);
			}
		}

		auto write_row = [&](std::ofstream& out, int b, float r, double z, double zstar, double pvalue, double qvalue, double lfdr) {
			out << target_w << "\t" << tchr << "\t" << tpos << "\t"
				<< b << "\t" << chroms_scan[b] << "\t" << pos_scan[b] << "\t"
				<< r << "\t" << nsamples << "\t"
				<< z << "\t" << zstar << "\t" << pvalue << "\t" << qvalue << "\t" << lfdr << "\n";
		};

		long long tested_local = 0;
		long long kept_local = 0;

		run_vector_block_fdr(
			bs, rs, nsamples, opt, rsample, rngs[(size_t)tid],
			do_distrib, do_distrib ? &td[(size_t)tid] : nullptr, distrib_sample,
			hfp, sfp, chr + "\tNA\tNA", tested_local, kept_local, write_row
		);

		tested_t[(size_t)tid] += tested_local;
		kept_t[(size_t)tid] += kept_local;
	}

	for (int t = 0; t < nthreads; ++t) {
		tested_pairs += tested_t[(size_t)t];
		kept_pairs += kept_t[(size_t)t];
	}

	try {
		merge_part_files(
			inter_hits, hit_part_paths,
			"wA\tchrA\tposA\twB\tchrB\tposB\tr\tn\tz\tzstar\tpvalue\tqvalue\tlocal_fdr\n"
		);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return false;
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(hit_part_paths[(size_t)t].c_str());

	try {
		merge_part_files(
			inter_summary, summary_part_paths,
			"chr\tdist_lo\tdist_hi\tn_pairs\tstatus\tmu0\tsigma0\tlambda\tpi0_pos\tpi0_neg\tn_hits_pos\tn_hits_neg\n"
		);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return false;
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(summary_part_paths[(size_t)t].c_str());

	if (do_distrib) {
		long long total_tested = 0;
		float global_min = std::numeric_limits<float>::infinity();
		float global_max = -std::numeric_limits<float>::infinity();

		std::vector<float> dsample;
		dsample.reserve((size_t)distrib_sample * (size_t)nthreads);

		for (int t = 0; t < nthreads; ++t) {
			const auto& d = td[(size_t)t];
			if (d.tested_pairs == 0)
				continue;

			total_tested += d.tested_pairs;
			if (d.min_r < global_min) global_min = d.min_r;
			if (d.max_r > global_max) global_max = d.max_r;

			dsample.insert(dsample.end(), d.sample.begin(), d.sample.end());
		}

		if ((int)dsample.size() > distrib_sample) {
			std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
			std::shuffle(dsample.begin(), dsample.end(), rng);
			dsample.resize((size_t)distrib_sample);
		}

		if (!distrib_path.empty()) {
			std::ofstream df(distrib_path);
			if (!df) {
				std::cerr << "Error: cannot write to " << distrib_path << "\n";
				return false;
			}

			df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\tmean_r2\tsd_r2\n";

			if (total_tested == 0 || dsample.empty()) {
				df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
			} else {
				float mean = std::numeric_limits<float>::quiet_NaN();
				float sd = std::numeric_limits<float>::quiet_NaN();
				compute_mean_sd_from_sample(dsample, mean, sd);

				float mean_r2 = std::numeric_limits<float>::quiet_NaN();
				float sd_r2 = std::numeric_limits<float>::quiet_NaN();
				{
					double mu2 = 0.0, m2v = 0.0;
					long long n2 = 0;
					for (float x : dsample) {
						if (!std::isfinite(x)) continue;
						++n2;
						double r2 = (double)x * (double)x;
						double dr2 = r2 - mu2;
						mu2 += dr2 / (double)n2;
						m2v += dr2 * (r2 - mu2);
					}
					if (n2 > 0) mean_r2 = (float)mu2;
					if (n2 > 1) sd_r2 = (float)std::sqrt(m2v / (double)(n2 - 1));
				}

				std::sort(dsample.begin(), dsample.end());

				float p01 = quantile_from_sorted(dsample, 0.01);
				float p05 = quantile_from_sorted(dsample, 0.05);
				float p25 = quantile_from_sorted(dsample, 0.25);
				float med = quantile_from_sorted(dsample, 0.50);
				float p75 = quantile_from_sorted(dsample, 0.75);
				float p95 = quantile_from_sorted(dsample, 0.95);
				float p99 = quantile_from_sorted(dsample, 0.99);

				df << total_tested
					<< "\t" << global_max
					<< "\t" << p99
					<< "\t" << p95
					<< "\t" << p75
					<< "\t" << med
					<< "\t" << p25
					<< "\t" << p05
					<< "\t" << p01
					<< "\t" << global_min
					<< "\t" << mean
					<< "\t" << sd
					<< "\t" << mean_r2
					<< "\t" << sd_r2
					<< "\n";
			}
		}

		if (!reservoir_path.empty() && !dsample.empty()) {
			std::ofstream rf(reservoir_path);
			if (!rf) {
				std::cerr << "Error: cannot write reservoir to " << reservoir_path << "\n";
				return false;
			}
			rf << "r\n";
			for (float x : dsample)
				rf << x << "\n";
		}
	}

	if (!opt.intra)
		return true;

	// Intrachromosomal portion (--intra is additive, see comment at the top
	// of this function): distance-binned same-chromosome blocks, merged
	// into the interchromosomal output above.
	const std::string intra_hits = out_path + ".fdrpart.intra";
	const std::string intra_summary = summary_path + ".fdrpart.intra";

	long long tested_intra = 0, kept_intra = 0;
	if (!scan_target_write_hits_excl_focus_fdr_intra_T(
		X_scan, chroms_scan, pos_scan, windows_by_chr_scan,
		hc_full, hi_excluding, opt, target_w, intra_hits, intra_summary, tested_intra, kept_intra, seed,
		"", 0, distrib_seed, "", freqs_scan
	))
		return false;

	bool ok = concat_tsv_drop_second_header(out_path, inter_hits, intra_hits)
		&& concat_tsv_drop_second_header(summary_path, inter_summary, intra_summary);

	std::remove(inter_hits.c_str());
	std::remove(inter_summary.c_str());
	std::remove(intra_hits.c_str());
	std::remove(intra_summary.c_str());

	if (!ok)
		return false;

	tested_pairs += tested_intra;
	kept_pairs += kept_intra;
	return true;
}

bool scan_target_write_hits_excl_focus_fdr(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	int target_w,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	uint64_t seed,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	auto fn = [](const HiComponentsWeighted& hc, const std::string& a, const std::string& b) {
		return hi_from_components_weighted_excluding(hc, a, b);
	};
	return scan_target_write_hits_excl_focus_fdr_T(
		X_scan, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, target_w, out_path, summary_path, tested_pairs, kept_pairs, seed,
		distrib_path, distrib_sample, distrib_seed, reservoir_path
	);
}

bool scan_target_write_hits_excl_focus_fdr(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const ScanOptions& opt,
	int target_w,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::vector<MarkerFreq>& freqs_scan,
	uint64_t seed,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	auto fn = [](const HiComponentsFreq& hc, const std::string& a, const std::string& b) {
		return hi_from_components_freq_excluding(hc, a, b);
	};
	return scan_target_write_hits_excl_focus_fdr_T(
		X_scan, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, target_w, out_path, summary_path, tested_pairs, kept_pairs, seed,
		distrib_path, distrib_sample, distrib_seed, reservoir_path, &freqs_scan
	);
}

// Per chromosome block, HI excludes only that block's chromosome; the
// vector and block's markers are residualized against that HI, then the
// calibrate+hit-call two-pass procedure runs against the resulting (b, r) pairs.
template<typename HC, typename HiExcludeFn>
static bool scan_vector_vs_windows_write_hits_excl_focus_fdr_T(
	const Eigen::MatrixXf& X_scan,
	const Eigen::VectorXf& v_raw,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HC& hc_full,
	HiExcludeFn hi_excluding,
	const ScanOptions& opt,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	uint64_t seed,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path,
	const std::vector<MarkerFreq>* freqs_scan = nullptr
) {
	// --intra has no meaning for a sample-level vector (no chromosome of its
	// own to be "intra" relative to); silently ignored here, same as the
	// non-FDR sample-haplo scan.

	const int nsamples = opt.nsamples;
	if (v_raw.size() != nsamples) {
		std::cerr << "Error: scan_vector_vs_windows_write_hits_excl_focus_fdr: v_raw size != nsamples\n";
		return false;
	}

	const int rsample = opt.fdr_sample > 0 ? opt.fdr_sample : 200000;

	std::cout << "[scan-fdr-loco-samplehaplo] target_fdr=" << opt.fdr_target
		  << " fdr_sample=" << rsample
		  << " fdr_min_pairs=" << opt.fdr_min_pairs
		  << " fdr_lambda=" << opt.fdr_lambda_cut << "\n";

	tested_pairs = 0;
	kept_pairs = 0;

	const bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
	nthreads = opt.threads;
	if (nthreads < 1) nthreads = 1;
	omp_set_num_threads(nthreads);
	#endif
	if (nthreads > (int)chr_order_scan.size())
		nthreads = std::max(1, (int)chr_order_scan.size());

	std::vector<std::string> hit_part_paths((size_t)nthreads);
	std::vector<std::string> summary_part_paths((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t) {
		hit_part_paths[(size_t)t] = out_path + ".part." + std::to_string(t);
		summary_part_paths[(size_t)t] = summary_path + ".part." + std::to_string(t);
		std::remove(hit_part_paths[(size_t)t].c_str());
		std::remove(summary_part_paths[(size_t)t].c_str());
	}

	std::vector<long long> tested_t((size_t)nthreads, 0);
	std::vector<long long> kept_t((size_t)nthreads, 0);
	std::vector<std::mt19937_64> rngs((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t)
		rngs[(size_t)t] = std::mt19937_64(seed + (uint64_t)t);

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	const int C = (int)chr_order_scan.size();

	#ifdef ADMIXLD_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int c = 0; c < C; ++c) {
		int tid = 0;
		#ifdef ADMIXLD_HAS_OPENMP
		tid = omp_get_thread_num();
		if (tid >= nthreads) tid = tid % nthreads;
		#endif

		const auto& chr = chr_order_scan[c];
		const auto& idx = windows_by_chr_scan.at(chr);
		const int m = (int)idx.size();
		if (m == 0)
			continue;

		Eigen::VectorXf h = hi_excluding(hc_full, chr, std::string(""));

		int v_valid = 0;
		Eigen::VectorXf vZ;
		try {
			vZ = residualize_and_zscore_vector(v_raw, h, v_valid);
		} catch (...) {
			continue;
		}

		int n_valid = 0;
		Eigen::MatrixXf Zc;
		try {
			Zc = residualize_and_zscore_subset(X_scan, h, idx, n_valid, freqs_scan);
		} catch (...) {
			continue;
		}

		std::ofstream hfp(hit_part_paths[(size_t)tid], std::ios::out | std::ios::app);
		std::ofstream sfp(summary_part_paths[(size_t)tid], std::ios::out | std::ios::app);
		if (!hfp || !sfp) {
			#ifdef ADMIXLD_HAS_OPENMP
			#pragma omp critical
			#endif
			{
				std::cerr << "Error: cannot write part files for thread " << tid << "\n";
			}
			continue;
		}

		// Per-sample NaN-aware dot product: unlike the target scan's marker
		// columns (all-or-nothing NaN), vZ can have genuine per-sample
		// missingness carried over from the raw --sample-haplo vector.
		std::vector<int> bs;
		std::vector<float> rs;
		bs.reserve((size_t)m);
		rs.reserve((size_t)m);

		const float denom = 1.0f / (float)(nsamples - 1);
		for (int j = 0; j < m; ++j) {
			const int w = idx[j];

			double dot = 0.0;
			for (int i = 0; i < nsamples; ++i) {
				const float zi = Zc(i, j);
				const float vi = vZ(i);
				if (!std::isnan(zi) && !std::isnan(vi))
					dot += (double)zi * (double)vi;
			}

			const float r = (float)(dot * denom);
			if (!std::isfinite(r))
				continue;

			bs.push_back(w);
			rs.push_back(r);
		}

		auto write_row = [&](std::ofstream& out, int b, float r, double z, double zstar, double pvalue, double qvalue, double lfdr) {
			out << "sample_haplo" << "\t" << b << "\t" << chroms_scan[b] << "\t" << pos_scan[b] << "\t"
				<< r << "\t" << nsamples << "\t"
				<< z << "\t" << zstar << "\t" << pvalue << "\t" << qvalue << "\t" << lfdr << "\n";
		};

		long long tested_local = 0;
		long long kept_local = 0;

		run_vector_block_fdr(
			bs, rs, nsamples, opt, rsample, rngs[(size_t)tid],
			do_distrib, do_distrib ? &td[(size_t)tid] : nullptr, distrib_sample,
			hfp, sfp, chr + "\tNA\tNA", tested_local, kept_local, write_row
		);

		tested_t[(size_t)tid] += tested_local;
		kept_t[(size_t)tid] += kept_local;
	}

	for (int t = 0; t < nthreads; ++t) {
		tested_pairs += tested_t[(size_t)t];
		kept_pairs += kept_t[(size_t)t];
	}

	try {
		merge_part_files(
			out_path, hit_part_paths,
			"tag\twB\tchrB\tposB\tr\tn\tz\tzstar\tpvalue\tqvalue\tlocal_fdr\n"
		);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return false;
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(hit_part_paths[(size_t)t].c_str());

	try {
		merge_part_files(
			summary_path, summary_part_paths,
			"chr\tdist_lo\tdist_hi\tn_pairs\tstatus\tmu0\tsigma0\tlambda\tpi0_pos\tpi0_neg\tn_hits_pos\tn_hits_neg\n"
		);
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return false;
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(summary_part_paths[(size_t)t].c_str());

	if (do_distrib) {
		long long total_tested = 0;
		float global_min = std::numeric_limits<float>::infinity();
		float global_max = -std::numeric_limits<float>::infinity();

		std::vector<float> dsample;
		dsample.reserve((size_t)distrib_sample * (size_t)nthreads);

		for (int t = 0; t < nthreads; ++t) {
			const auto& d = td[(size_t)t];
			if (d.tested_pairs == 0)
				continue;

			total_tested += d.tested_pairs;
			if (d.min_r < global_min) global_min = d.min_r;
			if (d.max_r > global_max) global_max = d.max_r;

			dsample.insert(dsample.end(), d.sample.begin(), d.sample.end());
		}

		if ((int)dsample.size() > distrib_sample) {
			std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
			std::shuffle(dsample.begin(), dsample.end(), rng);
			dsample.resize((size_t)distrib_sample);
		}

		if (!distrib_path.empty()) {
			std::ofstream df(distrib_path);
			if (!df) {
				std::cerr << "Error: cannot write to " << distrib_path << "\n";
				return false;
			}

			df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\tmean_r2\tsd_r2\n";

			if (total_tested == 0 || dsample.empty()) {
				df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
			} else {
				float mean = std::numeric_limits<float>::quiet_NaN();
				float sd = std::numeric_limits<float>::quiet_NaN();
				compute_mean_sd_from_sample(dsample, mean, sd);

				float mean_r2 = std::numeric_limits<float>::quiet_NaN();
				float sd_r2 = std::numeric_limits<float>::quiet_NaN();
				{
					double mu2 = 0.0, m2v = 0.0;
					long long n2 = 0;
					for (float x : dsample) {
						if (!std::isfinite(x)) continue;
						++n2;
						double r2 = (double)x * (double)x;
						double dr2 = r2 - mu2;
						mu2 += dr2 / (double)n2;
						m2v += dr2 * (r2 - mu2);
					}
					if (n2 > 0) mean_r2 = (float)mu2;
					if (n2 > 1) sd_r2 = (float)std::sqrt(m2v / (double)(n2 - 1));
				}

				std::sort(dsample.begin(), dsample.end());

				float p01 = quantile_from_sorted(dsample, 0.01);
				float p05 = quantile_from_sorted(dsample, 0.05);
				float p25 = quantile_from_sorted(dsample, 0.25);
				float med = quantile_from_sorted(dsample, 0.50);
				float p75 = quantile_from_sorted(dsample, 0.75);
				float p95 = quantile_from_sorted(dsample, 0.95);
				float p99 = quantile_from_sorted(dsample, 0.99);

				df << total_tested
					<< "\t" << global_max
					<< "\t" << p99
					<< "\t" << p95
					<< "\t" << p75
					<< "\t" << med
					<< "\t" << p25
					<< "\t" << p05
					<< "\t" << p01
					<< "\t" << global_min
					<< "\t" << mean
					<< "\t" << sd
					<< "\t" << mean_r2
					<< "\t" << sd_r2
					<< "\n";
			}
		}

		if (!reservoir_path.empty() && !dsample.empty()) {
			std::ofstream rf(reservoir_path);
			if (!rf) {
				std::cerr << "Error: cannot write reservoir to " << reservoir_path << "\n";
				return false;
			}
			rf << "r\n";
			for (float x : dsample)
				rf << x << "\n";
		}
	}

	return true;
}

bool scan_vector_vs_windows_write_hits_excl_focus_fdr(
	const Eigen::MatrixXf& X_scan,
	const Eigen::VectorXf& v_raw,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	uint64_t seed,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	auto fn = [](const HiComponentsWeighted& hc, const std::string& a, const std::string& b) {
		return hi_from_components_weighted_excluding(hc, a, b);
	};
	return scan_vector_vs_windows_write_hits_excl_focus_fdr_T(
		X_scan, v_raw, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, out_path, summary_path, tested_pairs, kept_pairs, seed,
		distrib_path, distrib_sample, distrib_seed, reservoir_path
	);
}

bool scan_vector_vs_windows_write_hits_excl_focus_fdr(
	const Eigen::MatrixXf& X_scan,
	const Eigen::VectorXf& v_raw,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::vector<MarkerFreq>& freqs_scan,
	uint64_t seed,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	auto fn = [](const HiComponentsFreq& hc, const std::string& a, const std::string& b) {
		return hi_from_components_freq_excluding(hc, a, b);
	};
	return scan_vector_vs_windows_write_hits_excl_focus_fdr_T(
		X_scan, v_raw, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, out_path, summary_path, tested_pairs, kept_pairs, seed,
		distrib_path, distrib_sample, distrib_seed, reservoir_path, &freqs_scan
	);
}
