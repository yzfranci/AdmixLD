#include "scan_markers.hpp"
#include "empirical_null.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <iomanip>
#include <cstdio>

#ifdef ADMIXLD_HAS_OPENMP
#include <omp.h>
#endif


std::unordered_map<std::string, std::vector<int>> group_by_chr(
	const std::vector<std::string>& chroms,
	std::vector<std::string>& chr_order
) {
	std::unordered_map<std::string, std::vector<int>> windows_by_chr;
	windows_by_chr.reserve(64);

	chr_order.clear();
	chr_order.reserve(64);

	int nwin = (int)chroms.size();
	for (int w = 0; w < nwin; ++w) {
		auto it = windows_by_chr.find(chroms[w]);
		if (it == windows_by_chr.end()) {
			windows_by_chr[chroms[w]] = std::vector<int>();
			windows_by_chr[chroms[w]].reserve(1024);
			chr_order.push_back(chroms[w]);
		}
		windows_by_chr[chroms[w]].push_back(w);
	}

	return windows_by_chr;
}

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

bool scan_markers_write_hits(
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;

	std::cout << "[scan] intra=" << (opt.intra ? 1 : 0)
		  << " max_dist=" << opt.max_dist
		  << " tile_size=" << opt.tile_size
		  << "\n";

	tested_pairs = 0;
	kept_pairs = 0;

	const bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
		nthreads = opt.threads;
		omp_set_num_threads(nthreads);
	#endif

	// One temp file per thread
	std::vector<std::string> part_paths((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t) {
		part_paths[(size_t)t] = out_path + ".part." + std::to_string(t);
		std::remove(part_paths[(size_t)t].c_str());
	}

	// Per-thread counters
	std::vector<long long> tested_t((size_t)nthreads, 0);
	std::vector<long long> kept_t((size_t)nthreads, 0);

	// Per-thread distrib tracking
	struct ThreadDistrib {
		long long tested_pairs = 0;
		float min_r = std::numeric_limits<float>::infinity();
		float max_r = -std::numeric_limits<float>::infinity();
		std::mt19937_64 rng;
		std::vector<float> sample;
	};

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	auto keep_hit = [&](float r) -> bool {
		if (opt.use_asym) {
			bool keep = false;
			if (opt.min_neg_r > 0.0f && r <= -opt.min_neg_r)
				keep = true;
			if (opt.min_pos_r > 0.0f && r >= opt.min_pos_r)
				keep = true;
			return keep;
		}
		return (std::fabs(r) >= opt.min_abs_r);
	};

	const float denom = 1.0f / (float)(nsamples - 1);

	if (opt.intra) {
		const int C = (int)chr_order.size();

		#ifdef ADMIXLD_HAS_OPENMP
		#pragma omp parallel for schedule(dynamic)
		#endif
		for (int c = 0; c < C; ++c) {
			int tid = 0;
			#ifdef ADMIXLD_HAS_OPENMP
			tid = omp_get_thread_num();
			#endif

			const auto& chr = chr_order[c];
			const auto& idx = windows_by_chr.at(chr);
			const int m = (int)idx.size();

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

			Eigen::MatrixXf A(nsamples, tile_size);
			Eigen::MatrixXf B(nsamples, tile_size);
			Eigen::MatrixXf R(tile_size, tile_size);

			long long tested_local = 0;
			long long kept_local = 0;

			auto consider_distrib = [&](float r) {
				if (!do_distrib)
					return;

				auto& d = td[(size_t)tid];

				if (r < d.min_r) d.min_r = r;
				if (r > d.max_r) d.max_r = r;

				++d.tested_pairs;

				if ((int)d.sample.size() < distrib_sample) {
					d.sample.push_back(r);
				} else {
					std::uniform_int_distribution<long long> U(0, d.tested_pairs - 1);
					long long j = U(d.rng);
					if (j < distrib_sample)
						d.sample[(size_t)j] = r;
				}
			};

			for (int i0 = 0; i0 < m; i0 += tile_size) {
				const int b1 = std::min(tile_size, m - i0);

				for (int k = 0; k < b1; ++k)
					A.col(k) = Z.col(idx[i0 + k]);

				for (int j0 = i0; j0 < m; j0 += tile_size) {
					const int b2 = std::min(tile_size, m - j0);

					for (int k = 0; k < b2; ++k)
						B.col(k) = Z.col(idx[j0 + k]);

					R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
					R.topLeftCorner(b1, b2) *= denom;

					for (int ia = 0; ia < b1; ++ia) {
						const int a = idx[i0 + ia];
						const int posA = pos[a];

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
									return value < pos[widx];
								}
							);

							jb_end = (int)std::distance(begin_it, ub);
						}

						if (jb_end <= jb_start)
							continue;

						for (int ib = jb_start; ib < jb_end; ++ib) {
							const int b = idx[j0 + ib];

							if (opt.min_dist >= 0 && (pos[b] - posA) < opt.min_dist)
								continue;

							const float r = R(ia, ib);

							++tested_local;
							consider_distrib(r);

							if (keep_hit(r)) {
								ofp << a << "\t" << chroms[a] << "\t" << pos[a] << "\t"
									<< b << "\t" << chroms[b] << "\t" << pos[b] << "\t"
									<< r << "\t" << nsamples << "\n";
								++kept_local;
							}
						}
					}
				}
			}

			tested_t[(size_t)tid] += tested_local;
			kept_t[(size_t)tid] += kept_local;
		}
	} else {
		const int C = (int)chr_order.size();

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
			#endif

			const int c1 = jobs[(size_t)j].first;
			const int c2 = jobs[(size_t)j].second;

			const auto& chr1 = chr_order[c1];
			const auto& chr2 = chr_order[c2];
			const auto& idx1 = windows_by_chr.at(chr1);
			const auto& idx2 = windows_by_chr.at(chr2);
			const int m1 = (int)idx1.size();
			const int m2 = (int)idx2.size();

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

			Eigen::MatrixXf A(nsamples, tile_size);
			Eigen::MatrixXf B(nsamples, tile_size);
			Eigen::MatrixXf R(tile_size, tile_size);

			long long tested_local = 0;
			long long kept_local = 0;

			auto consider_distrib = [&](float r) {
				if (!do_distrib)
					return;

				auto& d = td[(size_t)tid];

				if (r < d.min_r) d.min_r = r;
				if (r > d.max_r) d.max_r = r;

				++d.tested_pairs;

				if ((int)d.sample.size() < distrib_sample) {
					d.sample.push_back(r);
				} else {
					std::uniform_int_distribution<long long> U(0, d.tested_pairs - 1);
					long long j2 = U(d.rng);
					if (j2 < distrib_sample)
						d.sample[(size_t)j2] = r;
				}
			};

			for (int i0 = 0; i0 < m1; i0 += tile_size) {
				const int b1 = std::min(tile_size, m1 - i0);

				for (int k = 0; k < b1; ++k)
					A.col(k) = Z.col(idx1[i0 + k]);

				for (int j0 = 0; j0 < m2; j0 += tile_size) {
					const int b2 = std::min(tile_size, m2 - j0);

					for (int k = 0; k < b2; ++k)
						B.col(k) = Z.col(idx2[j0 + k]);

					R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
					R.topLeftCorner(b1, b2) *= denom;

					for (int ia = 0; ia < b1; ++ia) {
						const int a = idx1[i0 + ia];
						for (int ib = 0; ib < b2; ++ib) {
							const int b = idx2[j0 + ib];
							const float r = R(ia, ib);

							++tested_local;
							consider_distrib(r);

							if (keep_hit(r)) {
								ofp << a << "\t" << chroms[a] << "\t" << pos[a] << "\t"
									<< b << "\t" << chroms[b] << "\t" << pos[b] << "\t"
									<< r << "\t" << nsamples << "\n";
								++kept_local;
							}
						}
					}
				}
			}

			tested_t[(size_t)tid] += tested_local;
			kept_t[(size_t)tid] += kept_local;
		}
	}

	// Aggregate counts
	for (int t = 0; t < nthreads; ++t) {
		tested_pairs += tested_t[(size_t)t];
		kept_pairs += kept_t[(size_t)t];
	}

	// Merge part files into final output
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

	// Cleanup part files
	for (int t = 0; t < nthreads; ++t)
		std::remove(part_paths[(size_t)t].c_str());

	// Write distrib summary and/or reservoir
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

		double mu = 0.0, m2_acc = 0.0, mu_r2 = 0.0, m2_r2 = 0.0;
		long long n = 0;
		for (float x : dsample) {
			if (!std::isfinite(x))
				continue;
			++n;
			double dx = (double)x - mu;
			mu += dx / (double)n;
			m2_acc += dx * ((double)x - mu);
			double r2 = (double)x * (double)x;
			double dr2 = r2 - mu_r2;
			mu_r2 += dr2 / (double)n;
			m2_r2 += dr2 * (r2 - mu_r2);
		}

		float mean = (n > 0) ? (float)mu : std::numeric_limits<float>::quiet_NaN();
		float sd = (n > 1) ? (float)std::sqrt(m2_acc / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();
		float mean_r2 = (n > 0) ? (float)mu_r2 : std::numeric_limits<float>::quiet_NaN();
		float sd_r2 = (n > 1) ? (float)std::sqrt(m2_r2 / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();

		std::sort(dsample.begin(), dsample.end());

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

// Blocks are (chromosome, distance-bin) instead of chromosome-pairs; bin
// edges are log-spaced between --min-dist and --max-dist, shared across all
// chromosomes. Otherwise the same two-pass calibrate/hit-call structure as
// the interchromosomal scan below.
static bool scan_markers_write_hits_fdr_intra(
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
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

	// Global log-spaced bin edges, shared across all chromosomes.
	long long effective_min = std::max(real_floor, 1LL);
	long long effective_max = opt.max_dist >= 0 ? (long long)opt.max_dist : 0;
	if (effective_max <= 0) {
		for (const auto& chr : chr_order) {
			const auto& idx = windows_by_chr.at(chr);
			if (idx.size() < 2)
				continue;
			long long span = (long long)pos[idx.back()] - (long long)pos[idx.front()];
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

	std::cout << "[scan-fdr-intra] target_fdr=" << opt.fdr_target
		  << " fdr_sample=" << rsample
		  << " fdr_min_pairs=" << opt.fdr_min_pairs
		  << " fdr_lambda=" << opt.fdr_lambda_cut
		  << " bins=" << n_bins
		  << " min_dist=" << real_floor
		  << " max_dist=" << effective_max
		  << " tile_size=" << tile_size << "\n";

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
		nthreads = opt.threads;
		omp_set_num_threads(nthreads);
	#endif

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

	struct ThreadDistrib {
		long long tested_pairs = 0;
		float min_r = std::numeric_limits<float>::infinity();
		float max_r = -std::numeric_limits<float>::infinity();
		std::mt19937_64 rng;
		std::vector<float> sample;
	};

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	const int C = (int)chr_order.size();
	std::vector<std::pair<int,int>> jobs;	// (chr_index, bin_index)
	jobs.reserve((size_t)C * (size_t)n_bins);
	for (int c = 0; c < C; ++c)
		for (int bidx = 0; bidx < n_bins; ++bidx)
			jobs.push_back({c, bidx});

	#ifdef ADMIXLD_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int j = 0; j < (int)jobs.size(); ++j) {
		int tid = 0;
		#ifdef ADMIXLD_HAS_OPENMP
		tid = omp_get_thread_num();
		#endif

		const int c = jobs[(size_t)j].first;
		const int bidx = jobs[(size_t)j].second;
		const auto& chr = chr_order[c];
		const auto& idx = windows_by_chr.at(chr);
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

		Eigen::MatrixXf A(nsamples, tile_size);
		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::MatrixXf R(tile_size, tile_size);

		// Only visits pairs with distance in [bin_lo, bin_hi), via
		// lower/upper_bound on the sorted position array.
		auto for_each_pair = [&](auto&& cb) {
			for (int i0 = 0; i0 < m; i0 += tile_size) {
				const int b1 = std::min(tile_size, m - i0);
				for (int k = 0; k < b1; ++k)
					A.col(k) = Z.col(idx[i0 + k]);

				for (int j0 = i0; j0 < m; j0 += tile_size) {
					const int b2 = std::min(tile_size, m - j0);
					for (int k = 0; k < b2; ++k)
						B.col(k) = Z.col(idx[j0 + k]);

					R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
					R.topLeftCorner(b1, b2) *= denom;

					for (int ia = 0; ia < b1; ++ia) {
						const int a = idx[i0 + ia];
						const int posA = pos[a];

						int jb_dedup_start = 0;
						if (j0 == i0)
							jb_dedup_start = ia + 1;

						auto begin_it = idx.begin() + j0;
						auto end_it = idx.begin() + j0 + b2;

						auto lo_it = std::lower_bound(begin_it, end_it, (long long)posA + bin_lo,
							[&](int widx, long long value) { return (long long)pos[widx] < value; });
						auto hi_it = std::upper_bound(begin_it, end_it, (long long)posA + bin_hi - 1,
							[&](long long value, int widx) { return value < (long long)pos[widx]; });

						const int jb_start = std::max(jb_dedup_start, (int)std::distance(begin_it, lo_it));
						const int jb_end = (int)std::distance(begin_it, hi_it);

						if (jb_end <= jb_start)
							continue;

						for (int ib = jb_start; ib < jb_end; ++ib) {
							const int b = idx[j0 + ib];
							cb(a, b, R(ia, ib));
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

			if (do_distrib) {
				auto& d = td[(size_t)tid];

				if (r < d.min_r) d.min_r = r;
				if (r > d.max_r) d.max_r = r;

				++d.tested_pairs;

				if ((int)d.sample.size() < distrib_sample) {
					d.sample.push_back(r);
				} else {
					std::uniform_int_distribution<long long> Ud(0, d.tested_pairs - 1);
					long long jd = Ud(d.rng);
					if (jd < distrib_sample)
						d.sample[(size_t)jd] = r;
				}
			}
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

				hfp << a << "\t" << chroms[a] << "\t" << pos[a] << "\t"
					<< b << "\t" << chroms[b] << "\t" << pos[b] << "\t"
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

	{
		std::ofstream of(out_path);
		if (!of) {
			std::cerr << "Error: cannot write to " << out_path << "\n";
			return false;
		}
		of << "wA\tchrA\tposA\twB\tchrB\tposB\tr\tn\tz\tzstar\tpvalue\tqvalue\tlocal_fdr\n";
		for (int t = 0; t < nthreads; ++t) {
			std::ifstream pf(hit_part_paths[(size_t)t]);
			if (!pf)
				continue;
			std::string line;
			while (std::getline(pf, line))
				of << line << "\n";
		}
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(hit_part_paths[(size_t)t].c_str());

	{
		std::ofstream of(summary_path);
		if (!of) {
			std::cerr << "Error: cannot write to " << summary_path << "\n";
			return false;
		}
		of << "chr\tdist_lo\tdist_hi\tn_pairs\tstatus\tmu0\tsigma0\tlambda\tpi0_pos\tpi0_neg\tn_hits_pos\tn_hits_neg\n";
		for (int t = 0; t < nthreads; ++t) {
			std::ifstream pf(summary_part_paths[(size_t)t]);
			if (!pf)
				continue;
			std::string line;
			while (std::getline(pf, line))
				of << line << "\n";
		}
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

		double mu = 0.0, m2_acc = 0.0, mu_r2 = 0.0, m2_r2 = 0.0;
		long long n = 0;
		for (float x : dsample) {
			if (!std::isfinite(x))
				continue;
			++n;
			double dx = (double)x - mu;
			mu += dx / (double)n;
			m2_acc += dx * ((double)x - mu);
			double r2 = (double)x * (double)x;
			double dr2 = r2 - mu_r2;
			mu_r2 += dr2 / (double)n;
			m2_r2 += dr2 * (r2 - mu_r2);
		}

		float mean = (n > 0) ? (float)mu : std::numeric_limits<float>::quiet_NaN();
		float sd = (n > 1) ? (float)std::sqrt(m2_acc / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();
		float mean_r2 = (n > 0) ? (float)mu_r2 : std::numeric_limits<float>::quiet_NaN();
		float sd_r2 = (n > 1) ? (float)std::sqrt(m2_r2 / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();

		std::sort(dsample.begin(), dsample.end());

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

bool scan_markers_write_hits_fdr(
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
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
	if (opt.intra) {
		return scan_markers_write_hits_fdr_intra(
			Z, chroms, pos, windows_by_chr, chr_order, opt,
			out_path, summary_path, tested_pairs, kept_pairs, seed,
			distrib_path, distrib_sample, distrib_seed, reservoir_path
		);
	}

	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;
	const float denom = 1.0f / (float)(nsamples - 1);
	const int rsample = opt.fdr_sample > 0 ? opt.fdr_sample : 200000;

	std::cout << "[scan-fdr] target_fdr=" << opt.fdr_target
		  << " fdr_sample=" << rsample
		  << " fdr_min_pairs=" << opt.fdr_min_pairs
		  << " fdr_lambda=" << opt.fdr_lambda_cut
		  << " tile_size=" << tile_size << "\n";

	tested_pairs = 0;
	kept_pairs = 0;

	const bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
		nthreads = opt.threads;
		omp_set_num_threads(nthreads);
	#endif

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

	struct ThreadDistrib {
		long long tested_pairs = 0;
		float min_r = std::numeric_limits<float>::infinity();
		float max_r = -std::numeric_limits<float>::infinity();
		std::mt19937_64 rng;
		std::vector<float> sample;
	};

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	const int C = (int)chr_order.size();
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
		#endif

		const int c1 = jobs[(size_t)j].first;
		const int c2 = jobs[(size_t)j].second;
		const auto& chr1 = chr_order[c1];
		const auto& chr2 = chr_order[c2];
		const auto& idx1 = windows_by_chr.at(chr1);
		const auto& idx2 = windows_by_chr.at(chr2);
		const int m1 = (int)idx1.size();
		const int m2 = (int)idx2.size();

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

		Eigen::MatrixXf A(nsamples, tile_size);
		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::MatrixXf R(tile_size, tile_size);

		// Tile loop over this chromosome-pair block; invokes cb(a, b, r) for
		// every pair. Called twice per block (calibration pass, hit pass).
		auto for_each_pair = [&](auto&& cb) {
			for (int i0 = 0; i0 < m1; i0 += tile_size) {
				const int b1 = std::min(tile_size, m1 - i0);
				for (int k = 0; k < b1; ++k)
					A.col(k) = Z.col(idx1[i0 + k]);

				for (int j0 = 0; j0 < m2; j0 += tile_size) {
					const int b2 = std::min(tile_size, m2 - j0);
					for (int k = 0; k < b2; ++k)
						B.col(k) = Z.col(idx2[j0 + k]);

					R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
					R.topLeftCorner(b1, b2) *= denom;

					for (int ia = 0; ia < b1; ++ia) {
						const int a = idx1[i0 + ia];
						for (int ib = 0; ib < b2; ++ib) {
							const int b = idx2[j0 + ib];
							cb(a, b, R(ia, ib));
						}
					}
				}
			}
		};

		// Pass 1: exact pair count + reservoir sample of z = atanh(r).
		// Also collects the global raw-r distrib reservoir (if requested),
		// since this pass always runs (even for blocks later skipped in pass 2).
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

			if (do_distrib) {
				auto& d = td[(size_t)tid];

				if (r < d.min_r) d.min_r = r;
				if (r > d.max_r) d.max_r = r;

				++d.tested_pairs;

				if ((int)d.sample.size() < distrib_sample) {
					d.sample.push_back(r);
				} else {
					std::uniform_int_distribution<long long> Ud(0, d.tested_pairs - 1);
					long long jd = Ud(d.rng);
					if (jd < distrib_sample)
						d.sample[(size_t)jd] = r;
				}
			}
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

				hfp << a << "\t" << chroms[a] << "\t" << pos[a] << "\t"
					<< b << "\t" << chroms[b] << "\t" << pos[b] << "\t"
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

	{
		std::ofstream of(out_path);
		if (!of) {
			std::cerr << "Error: cannot write to " << out_path << "\n";
			return false;
		}
		of << "wA\tchrA\tposA\twB\tchrB\tposB\tr\tn\tz\tzstar\tpvalue\tqvalue\tlocal_fdr\n";
		for (int t = 0; t < nthreads; ++t) {
			std::ifstream pf(hit_part_paths[(size_t)t]);
			if (!pf)
				continue;
			std::string line;
			while (std::getline(pf, line))
				of << line << "\n";
		}
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(hit_part_paths[(size_t)t].c_str());

	{
		std::ofstream of(summary_path);
		if (!of) {
			std::cerr << "Error: cannot write to " << summary_path << "\n";
			return false;
		}
		of << "chrA\tchrB\tn_pairs\tstatus\tmu0\tsigma0\tlambda\tpi0_pos\tpi0_neg\tn_hits_pos\tn_hits_neg\n";
		for (int t = 0; t < nthreads; ++t) {
			std::ifstream pf(summary_part_paths[(size_t)t]);
			if (!pf)
				continue;
			std::string line;
			while (std::getline(pf, line))
				of << line << "\n";
		}
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

		double mu = 0.0, m2_acc = 0.0, mu_r2 = 0.0, m2_r2 = 0.0;
		long long n = 0;
		for (float x : dsample) {
			if (!std::isfinite(x))
				continue;
			++n;
			double dx = (double)x - mu;
			mu += dx / (double)n;
			m2_acc += dx * ((double)x - mu);
			double r2 = (double)x * (double)x;
			double dr2 = r2 - mu_r2;
			mu_r2 += dr2 / (double)n;
			m2_r2 += dr2 * (r2 - mu_r2);
		}

		float mean = (n > 0) ? (float)mu : std::numeric_limits<float>::quiet_NaN();
		float sd = (n > 1) ? (float)std::sqrt(m2_acc / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();
		float mean_r2 = (n > 0) ? (float)mu_r2 : std::numeric_limits<float>::quiet_NaN();
		float sd_r2 = (n > 1) ? (float)std::sqrt(m2_r2 / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();

		std::sort(dsample.begin(), dsample.end());

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

// Shared core for the single-vector-vs-chromosome FDR scans. Each
// chromosome in chr_order (except skip_chrom) is one calibration+hit-calling
// block: v vs that chromosome's markers. write_row emits caller-specific columns.
template <typename RowWriterFn>
static bool scan_vector_fdr_core(
	const Eigen::VectorXf& v,
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const std::string& skip_chrom,
	const ScanOptions& opt,
	const std::string& out_path,
	const std::string& summary_path,
	const std::string& hits_header,
	long long& tested_pairs,
	long long& kept_pairs,
	uint64_t seed,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path,
	RowWriterFn write_row
) {
	const int nsamples = opt.nsamples;
	const int rsample = opt.fdr_sample > 0 ? opt.fdr_sample : 200000;

	std::cout << "[scan-fdr] target_fdr=" << opt.fdr_target
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
		omp_set_num_threads(nthreads);
	#endif

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

	struct ThreadDistrib {
		long long tested_pairs = 0;
		float min_r = std::numeric_limits<float>::infinity();
		float max_r = -std::numeric_limits<float>::infinity();
		std::mt19937_64 rng;
		std::vector<float> sample;
	};

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	const int C = (int)chr_order.size();

	#ifdef ADMIXLD_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int c = 0; c < C; ++c) {
		int tid = 0;
		#ifdef ADMIXLD_HAS_OPENMP
		tid = omp_get_thread_num();
		#endif

		const auto& chr = chr_order[c];
		if (!skip_chrom.empty() && chr == skip_chrom)
			continue;

		const auto& idx = windows_by_chr.at(chr);
		const int m2 = (int)idx.size();
		if (m2 == 0)
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

		// Compute r(b) once per marker in this block; reused for both the
		// calibration reservoir and the hit-calling pass (no need to
		// recompute, unlike the tiled marker-pair scan).
		std::vector<int> bs;
		std::vector<float> rs;
		bs.reserve((size_t)m2);
		rs.reserve((size_t)m2);

		for (int k = 0; k < m2; ++k) {
			const int b = idx[k];

			double dot = 0.0;
			int nvalid = 0;
			for (int i = 0; i < nsamples; ++i) {
				const float vi = v(i);
				const float zi = Z(i, b);
				if (!std::isnan(vi) && !std::isnan(zi)) {
					dot += (double)vi * (double)zi;
					++nvalid;
				}
			}
			if (nvalid < 3)
				continue;

			bs.push_back(b);
			rs.push_back((float)(dot / (double)(nvalid - 1)));
		}

		const long long m_local = (long long)bs.size();

		// Pass 1 (calibration reservoir) + distrib sampling.
		std::vector<double> reservoir;
		reservoir.reserve((size_t)std::min<long long>(rsample, m_local));
		auto& rng = rngs[(size_t)tid];

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

			if (do_distrib) {
				auto& d = td[(size_t)tid];

				if (r < d.min_r) d.min_r = r;
				if (r > d.max_r) d.max_r = r;

				++d.tested_pairs;

				if ((int)d.sample.size() < distrib_sample) {
					d.sample.push_back(r);
				} else {
					std::uniform_int_distribution<long long> Ud(0, d.tested_pairs - 1);
					long long jd = Ud(d.rng);
					if (jd < distrib_sample)
						d.sample[(size_t)jd] = r;
				}
			}
		}

		NullFit fit = fit_empirical_null(reservoir, m_local, nsamples, opt.fdr_target, opt.fdr_lambda_cut);

		long long kept_local = 0;
		long long hits_pos = 0, hits_neg = 0;

		// Pass 2: apply the fitted per-tail zstar thresholds and write hits.
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

		sfp << chr << "\tNA\tNA\t" << m_local << "\t"
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

	{
		std::ofstream of(out_path);
		if (!of) {
			std::cerr << "Error: cannot write to " << out_path << "\n";
			return false;
		}
		of << hits_header;
		for (int t = 0; t < nthreads; ++t) {
			std::ifstream pf(hit_part_paths[(size_t)t]);
			if (!pf)
				continue;
			std::string line;
			while (std::getline(pf, line))
				of << line << "\n";
		}
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(hit_part_paths[(size_t)t].c_str());

	{
		std::ofstream of(summary_path);
		if (!of) {
			std::cerr << "Error: cannot write to " << summary_path << "\n";
			return false;
		}
		of << "chr\tdist_lo\tdist_hi\tn_pairs\tstatus\tmu0\tsigma0\tlambda\tpi0_pos\tpi0_neg\tn_hits_pos\tn_hits_neg\n";
		for (int t = 0; t < nthreads; ++t) {
			std::ifstream pf(summary_part_paths[(size_t)t]);
			if (!pf)
				continue;
			std::string line;
			while (std::getline(pf, line))
				of << line << "\n";
		}
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

		double mu = 0.0, m2_acc = 0.0, mu_r2 = 0.0, m2_r2 = 0.0;
		long long n = 0;
		for (float x : dsample) {
			if (!std::isfinite(x))
				continue;
			++n;
			double dx = (double)x - mu;
			mu += dx / (double)n;
			m2_acc += dx * ((double)x - mu);
			double r2 = (double)x * (double)x;
			double dr2 = r2 - mu_r2;
			mu_r2 += dr2 / (double)n;
			m2_r2 += dr2 * (r2 - mu_r2);
		}

		float mean = (n > 0) ? (float)mu : std::numeric_limits<float>::quiet_NaN();
		float sd = (n > 1) ? (float)std::sqrt(m2_acc / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();
		float mean_r2 = (n > 0) ? (float)mu_r2 : std::numeric_limits<float>::quiet_NaN();
		float sd_r2 = (n > 1) ? (float)std::sqrt(m2_r2 / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();

		std::sort(dsample.begin(), dsample.end());

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

// Target marker vs its own chromosome's other markers, binned like
// scan_markers_write_hits_fdr_intra. r(b) is computed once per bin and
// reused for both the calibration and hit-calling passes.
static bool scan_target_write_hits_fdr_intra(
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
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
	const int nsamples = opt.nsamples;
	const int rsample = opt.fdr_sample > 0 ? opt.fdr_sample : 200000;
	const int n_bins = opt.fdr_intra_bins > 0 ? opt.fdr_intra_bins : 8;
	const long long real_floor = opt.min_dist >= 0 ? (long long)opt.min_dist : 0;

	tested_pairs = 0;
	kept_pairs = 0;

	const std::string& tchr = chroms[target_w];
	const int tpos = pos[target_w];
	const Eigen::VectorXf v = Z.col(target_w);

	const auto& idx = windows_by_chr.at(tchr);
	const int m = (int)idx.size();

	const bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	long long effective_min = std::max(real_floor, 1LL);
	long long effective_max = opt.max_dist >= 0 ? (long long)opt.max_dist : 0;
	if (effective_max <= 0 && m >= 2)
		effective_max = (long long)pos[idx.back()] - (long long)pos[idx.front()];
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

	std::cout << "[scan-fdr-target-intra] target_fdr=" << opt.fdr_target
		  << " fdr_sample=" << rsample
		  << " fdr_min_pairs=" << opt.fdr_min_pairs
		  << " fdr_lambda=" << opt.fdr_lambda_cut
		  << " bins=" << n_bins
		  << " min_dist=" << real_floor
		  << " max_dist=" << effective_max << "\n";

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
		nthreads = opt.threads;
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

	struct ThreadDistrib {
		long long tested_pairs = 0;
		float min_r = std::numeric_limits<float>::infinity();
		float max_r = -std::numeric_limits<float>::infinity();
		std::mt19937_64 rng;
		std::vector<float> sample;
	};

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

		// Distance is unsigned (target can have nearby markers on either side).
		std::vector<int> bs;
		std::vector<float> rs;
		bs.reserve((size_t)m);
		rs.reserve((size_t)m);

		for (int k = 0; k < m; ++k) {
			const int b = idx[k];
			if (b == target_w)
				continue;

			long long d = (long long)pos[b] - (long long)tpos;
			if (d < 0) d = -d;
			if (d < bin_lo || d >= bin_hi)
				continue;

			double dot = 0.0;
			for (int i = 0; i < nsamples; ++i)
				dot += (double)v(i) * (double)Z(i, b);

			const float r = (float)(dot / (double)(nsamples - 1));
			if (!std::isfinite(r))
				continue;

			bs.push_back(b);
			rs.push_back(r);
		}

		const long long m_local = (long long)bs.size();

		std::vector<double> reservoir;
		reservoir.reserve((size_t)std::min<long long>(rsample, m_local));
		auto& rng = rngs[(size_t)tid];

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

			if (do_distrib) {
				auto& d = td[(size_t)tid];
				if (r < d.min_r) d.min_r = r;
				if (r > d.max_r) d.max_r = r;
				++d.tested_pairs;
				if ((int)d.sample.size() < distrib_sample) {
					d.sample.push_back(r);
				} else {
					std::uniform_int_distribution<long long> Ud(0, d.tested_pairs - 1);
					long long jd = Ud(d.rng);
					if (jd < distrib_sample)
						d.sample[(size_t)jd] = r;
				}
			}
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

				hfp << target_w << "\t" << tchr << "\t" << tpos << "\t"
					<< b << "\t" << chroms[b] << "\t" << pos[b] << "\t"
					<< r << "\t" << nsamples << "\t"
					<< z << "\t" << zstar << "\t" << pvalue << "\t" << qvalue << "\t" << lfdr << "\n";

				++kept_local;
				if (is_pos_hit) ++hits_pos; else ++hits_neg;
			}
		}

		sfp << tchr << "\t" << bin_lo << "\t" << bin_hi << "\t" << m_local << "\t"
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

	{
		std::ofstream of(out_path);
		if (!of) {
			std::cerr << "Error: cannot write to " << out_path << "\n";
			return false;
		}
		of << "wA\tchrA\tposA\twB\tchrB\tposB\tr\tn\tz\tzstar\tpvalue\tqvalue\tlocal_fdr\n";
		for (int t = 0; t < nthreads; ++t) {
			std::ifstream pf(hit_part_paths[(size_t)t]);
			if (!pf)
				continue;
			std::string line;
			while (std::getline(pf, line))
				of << line << "\n";
		}
	}
	for (int t = 0; t < nthreads; ++t)
		std::remove(hit_part_paths[(size_t)t].c_str());

	{
		std::ofstream of(summary_path);
		if (!of) {
			std::cerr << "Error: cannot write to " << summary_path << "\n";
			return false;
		}
		of << "chr\tdist_lo\tdist_hi\tn_pairs\tstatus\tmu0\tsigma0\tlambda\tpi0_pos\tpi0_neg\tn_hits_pos\tn_hits_neg\n";
		for (int t = 0; t < nthreads; ++t) {
			std::ifstream pf(summary_part_paths[(size_t)t]);
			if (!pf)
				continue;
			std::string line;
			while (std::getline(pf, line))
				of << line << "\n";
		}
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

		double mu = 0.0, m2_acc = 0.0, mu_r2 = 0.0, m2_r2 = 0.0;
		long long n = 0;
		for (float x : dsample) {
			if (!std::isfinite(x))
				continue;
			++n;
			double dx = (double)x - mu;
			mu += dx / (double)n;
			m2_acc += dx * ((double)x - mu);
			double r2 = (double)x * (double)x;
			double dr2 = r2 - mu_r2;
			mu_r2 += dr2 / (double)n;
			m2_r2 += dr2 * (r2 - mu_r2);
		}

		float mean = (n > 0) ? (float)mu : std::numeric_limits<float>::quiet_NaN();
		float sd = (n > 1) ? (float)std::sqrt(m2_acc / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();
		float mean_r2 = (n > 0) ? (float)mu_r2 : std::numeric_limits<float>::quiet_NaN();
		float sd_r2 = (n > 1) ? (float)std::sqrt(m2_r2 / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();

		std::sort(dsample.begin(), dsample.end());

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

// --intra is additive: runs both the interchromosomal blocks
// (scan_vector_fdr_core) and the intrachromosomal distance-binned blocks
// (scan_target_write_hits_fdr_intra), then concatenates their output.
bool scan_target_write_hits_fdr(
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
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
	const Eigen::VectorXf v = Z.col(target_w);
	const std::string tchr = chroms[target_w];
	const int tpos = pos[target_w];
	const int nsamples = opt.nsamples;

	auto write_row = [&](std::ofstream& hfp, int b, float r, double z, double zstar, double pvalue, double qvalue, double lfdr) {
		hfp << target_w << "\t" << tchr << "\t" << tpos << "\t"
			<< b << "\t" << chroms[b] << "\t" << pos[b] << "\t"
			<< r << "\t" << nsamples << "\t"
			<< z << "\t" << zstar << "\t" << pvalue << "\t" << qvalue << "\t" << lfdr << "\n";
	};

	const std::string inter_hits = opt.intra ? out_path + ".fdrpart.inter" : out_path;
	const std::string inter_summary = opt.intra ? summary_path + ".fdrpart.inter" : summary_path;

	long long tested_inter = 0, kept_inter = 0;
	if (!scan_vector_fdr_core(
		v, Z, chroms, pos, windows_by_chr, chr_order, tchr,
		opt, inter_hits, inter_summary,
		"wA\tchrA\tposA\twB\tchrB\tposB\tr\tn\tz\tzstar\tpvalue\tqvalue\tlocal_fdr\n",
		tested_inter, kept_inter, seed,
		distrib_path, distrib_sample, distrib_seed, reservoir_path,
		write_row
	))
		return false;

	if (!opt.intra) {
		tested_pairs = tested_inter;
		kept_pairs = kept_inter;
		return true;
	}

	const std::string intra_hits = out_path + ".fdrpart.intra";
	const std::string intra_summary = summary_path + ".fdrpart.intra";

	long long tested_intra = 0, kept_intra = 0;
	if (!scan_target_write_hits_fdr_intra(
		Z, chroms, pos, windows_by_chr, opt, target_w,
		intra_hits, intra_summary, tested_intra, kept_intra, seed,
		"", 0, distrib_seed, ""
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

	tested_pairs = tested_inter + tested_intra;
	kept_pairs = kept_inter + kept_intra;
	return true;
}

bool scan_vector_vs_windows_write_hits_fdr(
	const Eigen::MatrixXf& Z,
	const Eigen::VectorXf& v,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
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
	// --intra has no meaning for a sample-level vector (no chromosome of its
	// own to be "intra" relative to); silently ignored here, same as the
	// non-FDR sample-haplo scan.
	const int nsamples = opt.nsamples;

	auto write_row = [&](std::ofstream& hfp, int b, float r, double z, double zstar, double pvalue, double qvalue, double lfdr) {
		hfp << "sample_haplo" << "\t" << b << "\t" << chroms[b] << "\t" << pos[b] << "\t"
			<< r << "\t" << nsamples << "\t"
			<< z << "\t" << zstar << "\t" << pvalue << "\t" << qvalue << "\t" << lfdr << "\n";
	};

	return scan_vector_fdr_core(
		v, Z, chroms, pos, windows_by_chr, chr_order, std::string(""),
		opt, out_path, summary_path,
		"tag\twB\tchrB\tposB\tr\tn\tz\tzstar\tpvalue\tqvalue\tlocal_fdr\n",
		tested_pairs, kept_pairs, seed,
		distrib_path, distrib_sample, distrib_seed, reservoir_path,
		write_row
	);
}

bool scan_target_write_hits(
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
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
	const bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;
	const float denom = 1.0f / (float)(nsamples - 1);

	tested_pairs = 0;
	kept_pairs = 0;

	const Eigen::VectorXf v = Z.col(target_w);

	const std::string& tchr = chroms[target_w];
	const int tpos = pos[target_w];

	int nthreads = 1;
	#ifdef ADMIXLD_HAS_OPENMP
		nthreads = opt.threads;
		omp_set_num_threads(nthreads);
	#endif

	std::vector<std::string> part_paths((size_t)nthreads);
	for (int t = 0; t < nthreads; ++t) {
		part_paths[(size_t)t] = out_path + ".part." + std::to_string(t);
		std::remove(part_paths[(size_t)t].c_str());
	}

	std::vector<long long> tested_t((size_t)nthreads, 0);
	std::vector<long long> kept_t((size_t)nthreads, 0);

	struct ThreadDistrib {
		long long tested_pairs = 0;
		float min_r = std::numeric_limits<float>::infinity();
		float max_r = -std::numeric_limits<float>::infinity();
		std::mt19937_64 rng;
		std::vector<float> sample;
	};

	std::vector<ThreadDistrib> td;
	if (do_distrib) {
		td.resize((size_t)nthreads);
		for (int t = 0; t < nthreads; ++t) {
			td[(size_t)t].rng = std::mt19937_64(distrib_seed + (uint64_t)t);
			td[(size_t)t].sample.reserve((size_t)distrib_sample);
		}
	}

	auto keep_hit = [&](float r) -> bool {
		if (opt.use_asym) {
			bool keep = false;
			if (opt.min_neg_r > 0.0f && r <= -opt.min_neg_r)
				keep = true;
			if (opt.min_pos_r > 0.0f && r >= opt.min_pos_r)
				keep = true;
			return keep;
		}
		return (std::fabs(r) >= opt.min_abs_r);
	};

	const int C = (int)chr_order.size();

	#ifdef ADMIXLD_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int c = 0; c < C; ++c) {
		int tid = 0;
		#ifdef ADMIXLD_HAS_OPENMP
		tid = omp_get_thread_num();
		#endif

		const auto& chr = chr_order[c];

		if (!opt.intra && chr == tchr)
			continue;

		const auto& idx = windows_by_chr.at(chr);
		const int m = (int)idx.size();

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

		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::RowVectorXf R(tile_size);

		long long tested_local = 0;
		long long kept_local = 0;

		auto consider_distrib = [&](float r) {
			if (!do_distrib)
				return;

			auto& d = td[(size_t)tid];

			if (r < d.min_r) d.min_r = r;
			if (r > d.max_r) d.max_r = r;

			++d.tested_pairs;

			if ((int)d.sample.size() < distrib_sample) {
				d.sample.push_back(r);
			} else {
				std::uniform_int_distribution<long long> U(0, d.tested_pairs - 1);
				long long j = U(d.rng);
				if (j < distrib_sample)
					d.sample[(size_t)j] = r;
			}
		};

		for (int j0 = 0; j0 < m; j0 += tile_size) {
			const int b2 = std::min(tile_size, m - j0);

			for (int k = 0; k < b2; ++k)
				B.col(k) = Z.col(idx[j0 + k]);

			R.head(b2).noalias() = v.transpose() * B.leftCols(b2);
			R.head(b2) *= denom;

			for (int ib = 0; ib < b2; ++ib) {
				const int b = idx[j0 + ib];

				if (b == target_w)
					continue;

				if (opt.intra && chr == tchr) {
					int d = pos[b] - tpos;
					if (d < 0) d = -d;
					if (opt.max_dist >= 0 && d > opt.max_dist)
						continue;
					if (opt.min_dist >= 0 && d < opt.min_dist)
						continue;
				}

				const float r = R(ib);

				++tested_local;
				consider_distrib(r);

				if (keep_hit(r)) {
					ofp << target_w << "\t" << chroms[target_w] << "\t" << pos[target_w] << "\t"
						<< b << "\t" << chroms[b] << "\t" << pos[b] << "\t"
						<< r << "\t" << nsamples << "\n";
					++kept_local;
				}
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

		if ((int)dsample.size() > distrib_sample) {
			std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
			std::shuffle(dsample.begin(), dsample.end(), rng);
			dsample.resize((size_t)distrib_sample);
		}

		double mu = 0.0, m2_acc = 0.0, mu_r2 = 0.0, m2_r2 = 0.0;
		long long n = 0;
		for (float x : dsample) {
			if (!std::isfinite(x))
				continue;
			++n;
			double dx = (double)x - mu;
			mu += dx / (double)n;
			m2_acc += dx * ((double)x - mu);
			double r2 = (double)x * (double)x;
			double dr2 = r2 - mu_r2;
			mu_r2 += dr2 / (double)n;
			m2_r2 += dr2 * (r2 - mu_r2);
		}

		float mean = (n > 0) ? (float)mu : std::numeric_limits<float>::quiet_NaN();
		float sd = (n > 1) ? (float)std::sqrt(m2_acc / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();
		float mean_r2 = (n > 0) ? (float)mu_r2 : std::numeric_limits<float>::quiet_NaN();
		float sd_r2 = (n > 1) ? (float)std::sqrt(m2_r2 / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();

		std::sort(dsample.begin(), dsample.end());

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

bool scan_vector_vs_windows_write_hits(
	const Eigen::MatrixXf& Z,
	const Eigen::VectorXf& v,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path,
	int distrib_sample,
	uint64_t distrib_seed,
	const std::string& reservoir_path
) {
	std::ofstream of(out_path);
	if (!of) {
		std::cerr << "Error: cannot write to " << out_path << "\n";
		return false;
	}

	const int nsamples = opt.nsamples;
	const int nwin = (int)Z.cols();
	const float denom = 1.0f / (float)(nsamples - 1);

	if (v.size() != nsamples) {
		std::cerr << "Error: scan_vector_vs_windows_write_hits: v size != nsamples\n";
		return false;
	}

	(void)nwin;

	of << "tag\twB\tchrB\tposB\tr\tn\n";

	bool do_distrib = (!distrib_path.empty() || !reservoir_path.empty());
	std::mt19937_64 drng(distrib_seed);

	struct DistribSummary {
		float min_r = std::numeric_limits<float>::infinity();
		float max_r = -std::numeric_limits<float>::infinity();
		long long tested_pairs = 0;
	};

	DistribSummary ds;
	std::vector<float> dsample;
	if (do_distrib) {
		if (distrib_sample <= 0)
			distrib_sample = 200000;
		dsample.reserve((size_t)distrib_sample);
	}

	auto consider_distrib = [&](float r) {
		if (!do_distrib)
			return;

		if (!std::isfinite(r))
			return;

		if (r < ds.min_r) ds.min_r = r;
		if (r > ds.max_r) ds.max_r = r;

		++ds.tested_pairs;

		if ((int)dsample.size() < distrib_sample) {
			dsample.push_back(r);
		} else {
			std::uniform_int_distribution<long long> U(0, ds.tested_pairs - 1);
			long long j = U(drng);
			if (j < distrib_sample)
				dsample[(size_t)j] = r;
		}
	};

	auto keep_hit = [&](float r) -> bool {
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
	};

	tested_pairs = 0;
	kept_pairs = 0;

	for (const auto& chr : chr_order) {
		const auto& idx = windows_by_chr.at(chr);

		for (int k = 0; k < (int)idx.size(); ++k) {
			int w = idx[k];

			double dot = 0.0;
			for (int i = 0; i < nsamples; ++i) {
				float zi = Z(i, w);
				float vi = v(i);
				if (!std::isnan(zi) && !std::isnan(vi))
					dot += (double)zi * (double)vi;
			}

			float r = (float)(dot * denom);

			++tested_pairs;
			consider_distrib(r);

			if (keep_hit(r)) {
				of << "sample_haplo"
					<< "\t" << w
					<< "\t" << chroms[w]
					<< "\t" << pos[w]
					<< "\t" << r
					<< "\t" << nsamples
					<< "\n";
				++kept_pairs;
			}
		}
	}

	of.close();

	if (do_distrib) {
		if ((int)dsample.size() > distrib_sample) {
			std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
			std::shuffle(dsample.begin(), dsample.end(), rng);
			dsample.resize((size_t)distrib_sample);
		}

		double mu = 0.0, m2_acc = 0.0, mu_r2 = 0.0, m2_r2 = 0.0;
		long long n = 0;
		for (float x : dsample) {
			if (!std::isfinite(x))
				continue;
			++n;
			double dx = (double)x - mu;
			mu += dx / (double)n;
			m2_acc += dx * ((double)x - mu);
			double r2 = (double)x * (double)x;
			double dr2 = r2 - mu_r2;
			mu_r2 += dr2 / (double)n;
			m2_r2 += dr2 * (r2 - mu_r2);
		}

		float mean = (n > 0) ? (float)mu : std::numeric_limits<float>::quiet_NaN();
		float sd = (n > 1) ? (float)std::sqrt(m2_acc / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();
		float mean_r2 = (n > 0) ? (float)mu_r2 : std::numeric_limits<float>::quiet_NaN();
		float sd_r2 = (n > 1) ? (float)std::sqrt(m2_r2 / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();

		std::sort(dsample.begin(), dsample.end());

		if (!distrib_path.empty()) {
			std::ofstream df(distrib_path);
			if (!df) {
				std::cerr << "Error: cannot write distrib summary to " << distrib_path << "\n";
				return false;
			}
			df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\tmean_r2\tsd_r2\n";
			if (ds.tested_pairs == 0 || dsample.empty()) {
				df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
			} else {
				float p01 = quantile_from_sorted(dsample, 0.01);
				float p05 = quantile_from_sorted(dsample, 0.05);
				float p25 = quantile_from_sorted(dsample, 0.25);
				float med = quantile_from_sorted(dsample, 0.50);
				float p75 = quantile_from_sorted(dsample, 0.75);
				float p95 = quantile_from_sorted(dsample, 0.95);
				float p99 = quantile_from_sorted(dsample, 0.99);
				df << ds.tested_pairs
					<< "\t" << ds.max_r
					<< "\t" << p99
					<< "\t" << p95
					<< "\t" << p75
					<< "\t" << med
					<< "\t" << p25
					<< "\t" << p05
					<< "\t" << p01
					<< "\t" << ds.min_r
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

