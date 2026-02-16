#include "scan_blocks.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <iomanip>
#include <cstdio>

#ifdef ADFINDER_HAS_OPENMP
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

bool scan_blocks_write_hits(
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
	uint64_t distrib_seed
) {
	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;

	std::cout << "[scan] intra=" << (opt.intra ? 1 : 0)
		  << " max_dist=" << opt.max_dist
		  << " tile_size=" << opt.tile_size
		  << "\n";

	tested_pairs = 0;
	kept_pairs = 0;

	const bool do_distrib = (!distrib_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	int nthreads = 1;
	#ifdef ADFINDER_HAS_OPENMP
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

		#ifdef ADFINDER_HAS_OPENMP
		#pragma omp parallel for schedule(dynamic)
		#endif
		for (int c = 0; c < C; ++c) {
			int tid = 0;
			#ifdef ADFINDER_HAS_OPENMP
			tid = omp_get_thread_num();
			#endif

			const auto& chr = chr_order[c];
			const auto& idx = windows_by_chr.at(chr);
			const int m = (int)idx.size();

			std::ofstream ofp(part_paths[(size_t)tid], std::ios::out | std::ios::app);
			if (!ofp) {
				#ifdef ADFINDER_HAS_OPENMP
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

		#ifdef ADFINDER_HAS_OPENMP
		#pragma omp parallel for schedule(dynamic)
		#endif
		for (int j = 0; j < (int)jobs.size(); ++j) {
			int tid = 0;
			#ifdef ADFINDER_HAS_OPENMP
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
				#ifdef ADFINDER_HAS_OPENMP
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

	// Write distrib summary (merged)
	if (do_distrib) {
		std::ofstream df(distrib_path);
		if (!df) {
			std::cerr << "Error: cannot write to " << distrib_path << "\n";
			return false;
		}

		df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\n";

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

		if (total_tested == 0 || dsample.empty()) {
			df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
		} else {
			if ((int)dsample.size() > distrib_sample) {
				std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
				std::shuffle(dsample.begin(), dsample.end(), rng);
				dsample.resize((size_t)distrib_sample);
			}

			// mean/sd on the same reservoir sample used for quantiles
			double sum = 0.0;
			double sumsq = 0.0;
			int n = 0;
			for (float x : dsample) {
				if (!std::isfinite(x))
					continue;
				sum += (double)x;
				sumsq += (double)x * (double)x;
				++n;
			}

			float mean = 0.0f;
			float sd = 0.0f;
			if (n > 0) {
				mean = (float)(sum / (double)n);
				if (n > 1) {
					double var = (sumsq - (sum * sum) / (double)n) / (double)(n - 1);
					if (var < 0.0) var = 0.0; // numeric guard
					sd = (float)std::sqrt(var);
				}
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
				<< "\n";
		}
	}

	return true;
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
	uint64_t distrib_seed
) {
	const bool do_distrib = (!distrib_path.empty());
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
	#ifdef ADFINDER_HAS_OPENMP
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

	#ifdef ADFINDER_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int c = 0; c < C; ++c) {
		int tid = 0;
		#ifdef ADFINDER_HAS_OPENMP
		tid = omp_get_thread_num();
		#endif

		const auto& chr = chr_order[c];

		if (!opt.intra && chr == tchr)
			continue;

		const auto& idx = windows_by_chr.at(chr);
		const int m = (int)idx.size();

		std::ofstream ofp(part_paths[(size_t)tid], std::ios::out | std::ios::app);
		if (!ofp) {
			#ifdef ADFINDER_HAS_OPENMP
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

				if (opt.intra && opt.max_dist >= 0 && chr == tchr) {
					int d = pos[b] - tpos;
					if (d < 0) d = -d;
					if (d > opt.max_dist)
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
		std::ofstream df(distrib_path);
		if (!df) {
			std::cerr << "Error: cannot write to " << distrib_path << "\n";
			return false;
		}

		df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\n";

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

		if (total_tested == 0 || dsample.empty()) {
			df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
		} else {
			if ((int)dsample.size() > distrib_sample) {
				std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
				std::shuffle(dsample.begin(), dsample.end(), rng);
				dsample.resize((size_t)distrib_sample);
			}

			double sum = 0.0;
			double sumsq = 0.0;
			for (float x : dsample) {
				sum += (double)x;
				sumsq += (double)x * (double)x;
			}
			double mean = sum / (double)dsample.size();
			double var = 0.0;
			if (dsample.size() > 1) {
				var = (sumsq - (double)dsample.size() * mean * mean) / (double)(dsample.size() - 1);
				if (var < 0.0) var = 0.0;
			}
			double sd = std::sqrt(var);

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
				<< "\t" << (float)mean
				<< "\t" << (float)sd
				<< "\n";
		}
	}

	return true;
}

static inline int hi_bin_index(float h, int nbins, float lo = 0.0f, float hi = 1.0f) {
	if (!(h == h))
		return -1;
	if (nbins <= 1)
		return 0;
	if (h <= lo)
		return 0;
	if (h >= hi)
		return nbins - 1;
	float t = (h - lo) / (hi - lo);
	int b = (int)std::floor(t * (float)nbins);
	if (b < 0) b = 0;
	if (b >= nbins) b = nbins - 1;
	return b;
}

bool permute_interchrom_summary_chrblock(
	const Eigen::MatrixXf& Z,
	const Eigen::VectorXf& h,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	uint64_t seed,
	int n_perm,
	int sample_size,
	int hi_bins,
	std::vector<PermSummary>& summaries_out
) {
	if (opt.intra) {
		std::cerr << "Error: HI-bin+cshift permutation summary is implemented for interchrom scans only.\n";
		return false;
	}
	if (n_perm <= 0) {
		summaries_out.clear();
		return true;
	}
	if (sample_size <= 0)
		sample_size = 200000;
	if (hi_bins <= 0)
		hi_bins = 20;

	int nsamples = opt.nsamples;
	int tile_size = opt.tile_size;

	if (h.size() != nsamples) {
		std::cerr << "Error: HI vector length != nsamples\n";
		return false;
	}

	summaries_out.clear();
	summaries_out.resize((size_t)n_perm);

	const float denom = 1.0f / (float)(nsamples - 1);

	int nthreads = 1;
	#ifdef ADFINDER_HAS_OPENMP
		nthreads = opt.threads;
		if (nthreads < 1)
			nthreads = 1;
		omp_set_dynamic(0);
		omp_set_num_threads(nthreads);
	#endif

	std::cout << "Permutation test (HI-bin + circular-shift): " << n_perm << " replicates";
	#ifdef ADFINDER_HAS_OPENMP
		std::cout << " using " << nthreads << " threads";
	#endif
	std::cout << "\n";
	std::cout << "  hi_bins     = " << hi_bins << "\n";
	std::cout << "  sample_size = " << sample_size << "\n";

	// Precompute sample membership per HI bin
	std::vector<std::vector<int>> bin_members((size_t)hi_bins);
	for (int i = 0; i < nsamples; ++i) {
		int b = hi_bin_index(h(i), hi_bins, 0.0f, 1.0f);
		if (b >= 0)
			bin_members[(size_t)b].push_back(i);
	}

	#ifdef ADFINDER_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
	#endif
	for (int rep = 0; rep < n_perm; ++rep) {
		std::mt19937_64 rng(seed + (uint64_t)rep);

		Eigen::MatrixXf A(nsamples, tile_size);
		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::MatrixXf R(tile_size, tile_size);

		// For each chromosome:
		// - p_chr maps output-row i -> source-row p_chr[i], permuted within HI bins
		// - shift_chr[source-row] is a random circular shift offset for that source individual
		std::unordered_map<std::string, std::vector<int>> p_by_chr;
		std::unordered_map<std::string, std::vector<int>> shift_by_chr;
		p_by_chr.reserve(chr_order.size());
		shift_by_chr.reserve(chr_order.size());

		for (const auto& chr : chr_order) {
			std::vector<int> p(nsamples);
			for (int i = 0; i < nsamples; ++i)
				p[i] = i;

			// Permute within HI bins
			for (int b = 0; b < hi_bins; ++b) {
				const auto& members = bin_members[(size_t)b];
				if ((int)members.size() <= 1)
					continue;

				std::vector<int> shuf = members;
				std::shuffle(shuf.begin(), shuf.end(), rng);

				for (int k = 0; k < (int)members.size(); ++k)
					p[members[(size_t)k]] = shuf[(size_t)k];
			}

			p_by_chr.emplace(chr, std::move(p));

			// Circular shifts are defined on the SOURCE individual index
			const auto& idx = windows_by_chr.at(chr);
			int m = (int)idx.size();

			std::vector<int> sh(nsamples, 0);
			if (m > 1) {
				std::uniform_int_distribution<int> U(0, m - 1);
				for (int s = 0; s < nsamples; ++s)
					sh[s] = U(rng);
			}
			shift_by_chr.emplace(chr, std::move(sh));
		}

		PermSummary s;
		s.min_r = std::numeric_limits<float>::infinity();
		s.max_r = -std::numeric_limits<float>::infinity();

		std::vector<float> sample;
		sample.reserve((size_t)sample_size);
		long long seen = 0;

		int C = (int)chr_order.size();
		for (int c1 = 0; c1 < C; ++c1) {
			const auto& chr1 = chr_order[c1];
			const auto& idx1 = windows_by_chr.at(chr1);
			const auto& p1 = p_by_chr.at(chr1);
			const auto& sh1 = shift_by_chr.at(chr1);
			int m1 = (int)idx1.size();

			for (int c2 = c1 + 1; c2 < C; ++c2) {
				const auto& chr2 = chr_order[c2];
				const auto& idx2 = windows_by_chr.at(chr2);
				const auto& p2 = p_by_chr.at(chr2);
				const auto& sh2 = shift_by_chr.at(chr2);
				int m2 = (int)idx2.size();

				for (int i0 = 0; i0 < m1; i0 += tile_size) {
					int b1 = std::min(tile_size, m1 - i0);

					for (int k = 0; k < b1; ++k) {
						int r1 = i0 + k;	// rank within chr1 (sorted)
						for (int i = 0; i < nsamples; ++i) {
							int src = p1[i];
							int rr = r1;
							if (m1 > 1)
								rr = (r1 + sh1[src]) % m1;
							int w = idx1[(size_t)rr];
							A(i, k) = Z(src, w);
						}
					}

					for (int j0 = 0; j0 < m2; j0 += tile_size) {
						int b2 = std::min(tile_size, m2 - j0);

						for (int k = 0; k < b2; ++k) {
							int r2 = j0 + k;	// rank within chr2 (sorted)
							for (int i = 0; i < nsamples; ++i) {
								int src = p2[i];
								int rr = r2;
								if (m2 > 1)
									rr = (r2 + sh2[src]) % m2;
								int w = idx2[(size_t)rr];
								B(i, k) = Z(src, w);
							}
						}

						R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
						R.topLeftCorner(b1, b2) *= denom;

						for (int ia = 0; ia < b1; ++ia) {
							for (int ib = 0; ib < b2; ++ib) {
								float r = R(ia, ib);

								if (r < s.min_r) s.min_r = r;
								if (r > s.max_r) s.max_r = r;

								++seen;
								if ((int)sample.size() < sample_size) {
									sample.push_back(r);
								} else {
									std::uniform_int_distribution<long long> U(0, seen - 1);
									long long j = U(rng);
									if (j < sample_size)
										sample[(size_t)j] = r;
								}
							}
						}
					}
				}
			}
		}
		// Mean + SD from the reservoir sample (consistent with percentile estimates)
		double mu = 0.0;
		double m2 = 0.0;
		long long n = 0;

		for (float x : sample) {
			if (std::isnan(x))
				continue;
			++n;
			double dx = (double)x - mu;
			mu += dx / (double)n;
			double dx2 = (double)x - mu;
			m2 += dx * dx2;
		}

		s.mean = (n > 0) ? (float)mu : std::numeric_limits<float>::quiet_NaN();
		s.sd = (n > 1) ? (float)std::sqrt(m2 / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();

		std::sort(sample.begin(), sample.end());
		s.p01 = quantile_from_sorted(sample, 0.01);
		s.p05 = quantile_from_sorted(sample, 0.05);
		s.p25 = quantile_from_sorted(sample, 0.25);
		s.median = quantile_from_sorted(sample, 0.50);
		s.p75 = quantile_from_sorted(sample, 0.75);
		s.p95 = quantile_from_sorted(sample, 0.95);
		s.p99 = quantile_from_sorted(sample, 0.99);

		summaries_out[(size_t)rep] = s;
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
	uint64_t distrib_seed
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

	bool do_distrib = !distrib_path.empty();
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
				of << "sample_geno"
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
		std::ofstream df(distrib_path);
		if (!df) {
			std::cerr << "Error: cannot write distrib summary to " << distrib_path << "\n";
			return false;
		}

		df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\tmean\tsd\n";

		if (ds.tested_pairs == 0 || dsample.empty()) {
			df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
			return true;
		}

		if ((int)dsample.size() > distrib_sample) {
			std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
			std::shuffle(dsample.begin(), dsample.end(), rng);
			dsample.resize((size_t)distrib_sample);
		}

		double sum = 0.0;
		double sumsq = 0.0;
		for (float x : dsample) {
			sum += (double)x;
			sumsq += (double)x * (double)x;
		}
		double mean = sum / (double)dsample.size();
		double var = 0.0;
		if (dsample.size() > 1) {
			var = (sumsq - (double)dsample.size() * mean * mean) / (double)(dsample.size() - 1);
			if (var < 0.0) var = 0.0;
		}
		double sd = std::sqrt(var);

		std::sort(dsample.begin(), dsample.end());

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
			<< "\t" << (float)mean
			<< "\t" << (float)sd
			<< "\n";
	}

	return true;
}

bool permute_sample_vector_summary(
	const Eigen::MatrixXf& Z,
	const Eigen::VectorXf& gZ,
	const Eigen::VectorXf& h,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	uint64_t seed,
	int n_perm,
	int sample_size,
	int hi_bins,
	std::vector<PermSummary>& summaries_out
) {
	if (n_perm <= 0) {
		summaries_out.clear();
		return true;
	}

	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;

	if (gZ.size() != nsamples) {
		std::cerr << "Error: gZ length does not match nsamples\n";
		return false;
	}
	if (h.size() != nsamples) {
		std::cerr << "Error: HI vector length != nsamples\n";
		return false;
	}

	if (sample_size <= 0)
		sample_size = 200000;
	if (hi_bins <= 0)
		hi_bins = 20;

	summaries_out.clear();
	summaries_out.resize((size_t)n_perm);

	const float denom = 1.0f / (float)(nsamples - 1);

	// Precompute membership per HI bin
	std::vector<std::vector<int>> bin_members((size_t)hi_bins);
	for (int i = 0; i < nsamples; ++i) {
		int b = hi_bin_index(h(i), hi_bins, 0.0f, 1.0f);
		if (b >= 0)
			bin_members[(size_t)b].push_back(i);
	}

	int nthreads = 1;
	#ifdef ADFINDER_HAS_OPENMP
		nthreads = opt.threads;
		if (nthreads < 1)
			nthreads = 1;
		omp_set_dynamic(0);
		omp_set_num_threads(nthreads);
	#endif

	std::cout << "Permutation test (sample-geno; HI-bin + circular-shift): " << n_perm << " replicates";
	#ifdef ADFINDER_HAS_OPENMP
		std::cout << " using " << nthreads << " threads";
	#endif
	std::cout << "\n";
	std::cout << "  hi_bins     = " << hi_bins << "\n";
	std::cout << "  sample_size = " << sample_size << "\n";

	#ifdef ADFINDER_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
	#endif
	for (int rep = 0; rep < n_perm; ++rep) {
		std::mt19937_64 rng(seed + (uint64_t)rep);

		// perm[i] = source individual row used for output-row i (permuted within HI bins)
		std::vector<int> perm(nsamples);
		for (int i = 0; i < nsamples; ++i)
			perm[i] = i;

		for (int b = 0; b < hi_bins; ++b) {
			const auto& members = bin_members[(size_t)b];
			if ((int)members.size() <= 1)
				continue;

			std::vector<int> shuf = members;
			std::shuffle(shuf.begin(), shuf.end(), rng);

			for (int k = 0; k < (int)members.size(); ++k)
				perm[members[(size_t)k]] = shuf[(size_t)k];
		}

		// For each chromosome, define per-source-individual circular shift in rank space [0..m-1]
		std::unordered_map<std::string, std::vector<int>> shift_by_chr;
		shift_by_chr.reserve(chr_order.size());

		for (const auto& chr : chr_order) {
			const auto& idx = windows_by_chr.at(chr);
			int m = (int)idx.size();

			std::vector<int> sh(nsamples, 0);
			if (m > 1) {
				std::uniform_int_distribution<int> U(0, m - 1);
				for (int s = 0; s < nsamples; ++s)
					sh[s] = U(rng);
			}
			shift_by_chr.emplace(chr, std::move(sh));
		}

		PermSummary s;
		s.min_r = std::numeric_limits<float>::infinity();
		s.max_r = -std::numeric_limits<float>::infinity();

		std::vector<float> sample;
		sample.reserve((size_t)sample_size);
		long long seen = 0;

		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::RowVectorXf R(tile_size);

		// Loop over chromosomes to apply circular shift within each chromosome
		for (const auto& chr : chr_order) {
			const auto& idx = windows_by_chr.at(chr);
			const auto& sh = shift_by_chr.at(chr);
			int m = (int)idx.size();

			for (int j0 = 0; j0 < m; j0 += tile_size) {
				int b2 = std::min(tile_size, m - j0);

				// Build B: each column is a window; rows are HI-binned permuted samples; columns are circular-shifted per source-individual
				for (int k = 0; k < b2; ++k) {
					int rrank = j0 + k;	// rank within this chromosome

					for (int i = 0; i < nsamples; ++i) {
						int src = perm[i];
						int rr = rrank;

						if (m > 1)
							rr = (rrank + sh[src]) % m;

						int w = idx[(size_t)rr];
						B(i, k) = Z(src, w);
					}
				}

				R.head(b2).noalias() = gZ.transpose() * B.leftCols(b2);
				R.head(b2) *= denom;

				for (int k = 0; k < b2; ++k) {
					float r = R(k);

					if (r < s.min_r) s.min_r = r;
					if (r > s.max_r) s.max_r = r;

					++seen;
					if ((int)sample.size() < sample_size) {
						sample.push_back(r);
					} else {
						std::uniform_int_distribution<long long> U(0, seen - 1);
						long long j = U(rng);
						if (j < sample_size)
							sample[(size_t)j] = r;
					}
				}
			}
		}

		// mean + sd from the reservoir sample (same philosophy as your quantiles)
		double mu = 0.0;
		double m2 = 0.0;
		long long n = 0;

		for (float x : sample) {
			if (std::isnan(x))
				continue;
			++n;
			double dx = (double)x - mu;
			mu += dx / (double)n;
			double dx2 = (double)x - mu;
			m2 += dx * dx2;
		}

		s.mean = (n > 0) ? (float)mu : std::numeric_limits<float>::quiet_NaN();
		s.sd = (n > 1) ? (float)std::sqrt(m2 / (double)(n - 1)) : std::numeric_limits<float>::quiet_NaN();

		std::sort(sample.begin(), sample.end());
		s.p01 = quantile_from_sorted(sample, 0.01);
		s.p05 = quantile_from_sorted(sample, 0.05);
		s.p25 = quantile_from_sorted(sample, 0.25);
		s.median = quantile_from_sorted(sample, 0.50);
		s.p75 = quantile_from_sorted(sample, 0.75);
		s.p95 = quantile_from_sorted(sample, 0.95);
		s.p99 = quantile_from_sorted(sample, 0.99);

		summaries_out[(size_t)rep] = s;
	}

	return true;
}
