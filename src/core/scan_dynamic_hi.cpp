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

#ifdef ADFINDER_HAS_OPENMP
#include <omp.h>
#endif

#include "residualize.hpp"

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

bool scan_blocks_write_hits_excl_focus(
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
	uint64_t distrib_seed
) {
	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;

	tested_pairs = 0;
	kept_pairs = 0;

	const bool do_distrib = (!distrib_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	int nthreads = 1;
	#ifdef ADFINDER_HAS_OPENMP
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

		Eigen::VectorXf h = hi_from_components_weighted_excluding(hc_full, chr);

		int n_valid = 0;
		Eigen::MatrixXf Zc = residualize_and_zscore_subset(X_scan, h, idx, n_valid);

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

		Eigen::VectorXf h = hi_from_components_weighted_excluding(hc_full, chr1, chr2);

		int n_valid1 = 0;
		int n_valid2 = 0;
		Eigen::MatrixXf Z1 = residualize_and_zscore_subset(X_scan, h, idx1, n_valid1);
		Eigen::MatrixXf Z2 = residualize_and_zscore_subset(X_scan, h, idx2, n_valid2);

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

		#ifdef ADFINDER_HAS_OPENMP
		#pragma omp parallel for schedule(dynamic)
		#endif
		for (int c = 0; c < C; ++c) {
			int tid = 0;
			#ifdef ADFINDER_HAS_OPENMP
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

		#ifdef ADFINDER_HAS_OPENMP
		#pragma omp parallel for schedule(dynamic)
		#endif
		for (int j = 0; j < (int)jobs.size(); ++j) {
			int tid = 0;
			#ifdef ADFINDER_HAS_OPENMP
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

		std::ofstream df(distrib_path);
		if (!df) {
			std::cerr << "Error: cannot write to " << distrib_path << "\n";
			return false;
		}

		df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\n";

		if (total_tested == 0 || dsample.empty()) {
			df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
		} else {
			if ((int)dsample.size() > distrib_sample) {
				std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
				std::shuffle(dsample.begin(), dsample.end(), rng);
				dsample.resize((size_t)distrib_sample);
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
				<< "\n";
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

	const std::string& tchr = chroms_scan[target_w];
	const int tpos = pos_scan[target_w];

	int nthreads = 1;
	#ifdef ADFINDER_HAS_OPENMP
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

	const int C = (int)chr_order_scan.size();

	#ifdef ADFINDER_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int c = 0; c < C; ++c) {
		int tid = 0;
		#ifdef ADFINDER_HAS_OPENMP
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
			#ifdef ADFINDER_HAS_OPENMP
			#pragma omp critical
			#endif
			{
				std::cerr << "Error: cannot write to " << part_paths[(size_t)tid] << "\n";
			}
			continue;
		}

		Eigen::VectorXf h;

		bool same_chr = (chr == tchr);
		if (same_chr) {
			h = hi_from_components_weighted_excluding(hc_full, tchr);
		} else {
			h = hi_from_components_weighted_excluding(hc_full, tchr, chr);
		}

		Eigen::VectorXf v;

		if (same_chr) {
			int n_valid = 0;
			Eigen::MatrixXf Zc = residualize_and_zscore_subset(X_scan, h, idx, n_valid);

			int j_target = -1;
			for (int j = 0; j < m; ++j) {
				if (idx[j] == target_w) {
					j_target = j;
					break;
				}
			}
			if (j_target < 0)
				continue;

			v = Zc.col(j_target);

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
					long long j2 = U(d.rng);
					if (j2 < distrib_sample)
						d.sample[(size_t)j2] = r;
				}
			};

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

					if (opt.max_dist >= 0) {
						int dpos = pos_scan[b] - tpos;
						if (dpos < 0) dpos = -dpos;
						if (dpos > opt.max_dist)
							continue;
					}

					const float r = R(ib);

					++tested_local;
					consider_distrib(r);

					if (keep_hit(r)) {
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
			Eigen::MatrixXf Zt = residualize_and_zscore_subset(X_scan, h, tvec, n_valid_t);
			v = Zt.col(0);

			int n_valid_b = 0;
			Eigen::MatrixXf Zb = residualize_and_zscore_subset(X_scan, h, idx, n_valid_b);

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
					long long j2 = U(d.rng);
					if (j2 < distrib_sample)
						d.sample[(size_t)j2] = r;
				}
			};

			for (int j0 = 0; j0 < m; j0 += tile_size) {
				const int b2 = std::min(tile_size, m - j0);

				for (int k = 0; k < b2; ++k)
					B.col(k) = Zb.col(j0 + k);

				R.head(b2).noalias() = v.transpose() * B.leftCols(b2);
				R.head(b2) *= denom;

				for (int ib = 0; ib < b2; ++ib) {
					const int b = idx[j0 + ib];

					const float r = R(ib);

					++tested_local;
					consider_distrib(r);

					if (keep_hit(r)) {
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
		std::ofstream df(distrib_path);
		if (!df) {
			std::cerr << "Error: cannot write to " << distrib_path << "\n";
			return false;
		}

		df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\n";

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
			df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
		} else {
			if ((int)dsample.size() > distrib_sample) {
				std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
				std::shuffle(dsample.begin(), dsample.end(), rng);
				dsample.resize((size_t)distrib_sample);
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
				<< "\n";
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
	uint64_t distrib_seed
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

	const bool do_distrib = (!distrib_path.empty());
	const float denom = 1.0f / (float)(nsamples - 1);

	int nthreads = 1;
	#ifdef ADFINDER_HAS_OPENMP
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

	auto consider_distrib = [&](int tid, float r) {
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

	const int C = (int)chr_order_scan.size();

	#ifdef ADFINDER_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int c = 0; c < C; ++c) {
		int tid = 0;
		#ifdef ADFINDER_HAS_OPENMP
		tid = omp_get_thread_num();
		if (tid >= nthreads) tid = tid % nthreads;
		#endif

		const auto& chr = chr_order_scan[c];
		const auto& idx = windows_by_chr_scan.at(chr);
		const int m = (int)idx.size();
		if (m == 0)
			continue;

		// LOCO HI excluding this chromosome
		Eigen::VectorXf h = hi_from_components_weighted_excluding(hc_full, chr);

		// Residualize v_raw against h (and z-score)
		int v_valid = 0;
		Eigen::VectorXf vZ;
		try {
			vZ = residualize_and_zscore_vector(v_raw, h, v_valid);
		} catch (const std::exception&) {
			continue;
		}

		// Residualize windows on this chromosome against h (and z-score)
		int n_valid = 0;
		Eigen::MatrixXf Zc;
		try {
			Zc = residualize_and_zscore_subset(X_scan, h, idx, n_valid);
		} catch (const std::exception&) {
			continue;
		}

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

		long long tested_local = 0;
		long long kept_local = 0;

		// Compute one r per window in this chromosome
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

			++tested_local;
			consider_distrib(tid, r);

			if (keep_hit(r)) {
				ofp << "sample_geno"
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

	// Merge part files
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

	// Distrib summary
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

		std::ofstream df(distrib_path);
		if (!df) {
			std::cerr << "Error: cannot write distrib summary to " << distrib_path << "\n";
			return false;
		}

		df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\n";

		if (total_tested == 0 || dsample.empty()) {
			df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
		} else {
			if ((int)dsample.size() > distrib_sample) {
				std::mt19937_64 rng(distrib_seed ^ 0x9e3779b97f4a7c15ULL);
				std::shuffle(dsample.begin(), dsample.end(), rng);
				dsample.resize((size_t)distrib_sample);
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

bool permute_interchrom_summary_chrblock_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	uint64_t seed,
	int n_perm,
	int sample_size,
	int hi_bins,
	std::vector<PermSummary>& summaries_out
) {
	(void)chroms_scan;
	(void)pos_scan;

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

	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;
	const float denom = 1.0f / (float)(nsamples - 1);

	summaries_out.clear();
	summaries_out.resize((size_t)n_perm);

	int nthreads = 1;
	#ifdef ADFINDER_HAS_OPENMP
	nthreads = opt.threads;
	if (nthreads < 1) nthreads = 1;
	omp_set_dynamic(0);
	omp_set_num_threads(nthreads);
	#endif

	std::cout << "Permutation test (excl-focus, interchrom; HI-bin + circular-shift): " << n_perm << " replicates";
	#ifdef ADFINDER_HAS_OPENMP
	std::cout << " using " << nthreads << " threads";
	#endif
	std::cout << "\n";
	std::cout << "  hi_bins     = " << hi_bins << "\n";
	std::cout << "  sample_size = " << sample_size << "\n";

	#ifdef ADFINDER_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int rep = 0; rep < n_perm; ++rep) {
		std::mt19937_64 rng(seed + (uint64_t)rep);

		PermSummary s;
		s.min_r = std::numeric_limits<float>::infinity();
		s.max_r = -std::numeric_limits<float>::infinity();

		// reservoir sample for quantiles + mean/sd
		std::vector<float> sample;
		sample.reserve((size_t)sample_size);
		long long seen = 0;

		// Welford mean/sd over the reservoir sample (after it is finalized)
		// (we compute after sampling so mean/sd match the sample used for quantiles)

		Eigen::MatrixXf A(nsamples, tile_size);
		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::MatrixXf R(tile_size, tile_size);

		const int C = (int)chr_order_scan.size();
		for (int c1 = 0; c1 < C; ++c1) {
			const auto& chr1 = chr_order_scan[c1];
			const auto& idx1 = windows_by_chr_scan.at(chr1);
			const int m1 = (int)idx1.size();
			if (m1 == 0)
				continue;

			for (int c2 = c1 + 1; c2 < C; ++c2) {
				const auto& chr2 = chr_order_scan[c2];
				const auto& idx2 = windows_by_chr_scan.at(chr2);
				const int m2 = (int)idx2.size();
				if (m2 == 0)
					continue;

				// LOCO2 HI for this chromosome pair (same as dynamic scan)
				Eigen::VectorXf h = hi_from_components_weighted_excluding(hc_full, chr1, chr2);

				// Bin membership for this pair-specific HI
				std::vector<std::vector<int>> bin_members((size_t)hi_bins);
				for (int i = 0; i < nsamples; ++i) {
					int b = hi_bin_index(h(i), hi_bins, 0.0f, 1.0f);
					if (b >= 0)
						bin_members[(size_t)b].push_back(i);
				}

				// Residualize windows using LOCO2 HI
				int n_valid1 = 0;
				int n_valid2 = 0;
				Eigen::MatrixXf Z1 = residualize_and_zscore_subset(X_scan, h, idx1, n_valid1);
				Eigen::MatrixXf Z2 = residualize_and_zscore_subset(X_scan, h, idx2, n_valid2);

				// Permutation within HI bins for each chromosome in the pair:
				// pX[i] = source individual row used for output-row i
				std::vector<int> p1(nsamples), p2(nsamples);
				for (int i = 0; i < nsamples; ++i) {
					p1[i] = i;
					p2[i] = i;
				}

				for (int b = 0; b < hi_bins; ++b) {
					const auto& members = bin_members[(size_t)b];
					if ((int)members.size() <= 1)
						continue;

					std::vector<int> shuf = members;
					std::shuffle(shuf.begin(), shuf.end(), rng);

					for (int k = 0; k < (int)members.size(); ++k) {
						p1[members[(size_t)k]] = shuf[(size_t)k];
					}

					// independent shuffle for chr2 (same bin partitions, different mapping)
					shuf = members;
					std::shuffle(shuf.begin(), shuf.end(), rng);
					for (int k = 0; k < (int)members.size(); ++k) {
						p2[members[(size_t)k]] = shuf[(size_t)k];
					}
				}

				// Circular shift per *source* individual, in rank space within each chromosome
				std::vector<int> sh1(nsamples, 0), sh2(nsamples, 0);
				if (m1 > 1) {
					std::uniform_int_distribution<int> U(0, m1 - 1);
					for (int src = 0; src < nsamples; ++src)
						sh1[src] = U(rng);
				}
				if (m2 > 1) {
					std::uniform_int_distribution<int> U(0, m2 - 1);
					for (int src = 0; src < nsamples; ++src)
						sh2[src] = U(rng);
				}

				for (int i0 = 0; i0 < m1; i0 += tile_size) {
					int b1 = std::min(tile_size, m1 - i0);

					for (int k = 0; k < b1; ++k) {
						int r1 = i0 + k;	// rank within chr1
						for (int i = 0; i < nsamples; ++i) {
							int src = p1[i];
							int rr = r1;
							if (m1 > 1)
								rr = (r1 + sh1[src]) % m1;
							A(i, k) = Z1(src, rr);
						}
					}

					for (int j0 = 0; j0 < m2; j0 += tile_size) {
						int b2 = std::min(tile_size, m2 - j0);

						for (int k = 0; k < b2; ++k) {
							int r2 = j0 + k;	// rank within chr2
							for (int i = 0; i < nsamples; ++i) {
								int src = p2[i];
								int rr = r2;
								if (m2 > 1)
									rr = (r2 + sh2[src]) % m2;
								B(i, k) = Z2(src, rr);
							}
						}

						R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
						R.topLeftCorner(b1, b2) *= denom;

						for (int ia = 0; ia < b1; ++ia) {
							for (int ib = 0; ib < b2; ++ib) {
								float r = R(ia, ib);
								if (!std::isfinite(r))
									continue;

								if (r < s.min_r) s.min_r = r;
								if (r > s.max_r) s.max_r = r;

								++seen;
								if ((int)sample.size() < sample_size) {
									sample.push_back(r);
								} else {
									std::uniform_int_distribution<long long> U(0, seen - 1);
									long long jj = U(rng);
									if (jj < sample_size)
										sample[(size_t)jj] = r;
								}
							}
						}
					}
				}
			}
		}

		// mean + sd from the reservoir sample (same sample used for quantiles)
		double mu = 0.0;
		double m2 = 0.0;
		long long n = 0;
		for (float x : sample) {
			if (!std::isfinite(x))
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

bool permute_sample_vector_summary_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const Eigen::VectorXf& g_raw,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	uint64_t seed,
	int n_perm,
	int sample_size,
	int hi_bins,
	std::vector<PermSummary>& summaries_out
) {
	(void)chroms_scan;
	(void)pos_scan;

	if (n_perm <= 0) {
		summaries_out.clear();
		return true;
	}

	const int nsamples = opt.nsamples;
	if (g_raw.size() != nsamples) {
		std::cerr << "Error: permute_sample_vector_summary_excl_focus: g_raw length does not match nsamples\n";
		return false;
	}

	if (sample_size <= 0)
		sample_size = 200000;
	if (hi_bins <= 0)
		hi_bins = 20;

	const int tile_size = opt.tile_size;
	const float denom = 1.0f / (float)(nsamples - 1);

	const int C = (int)chr_order_scan.size();

	struct ChrCache {
		const std::vector<int>* idx = nullptr;
		Eigen::VectorXf gZ;	// nsamples
		Eigen::MatrixXf Zc;	// nsamples x m (columns are rank within chr)
		Eigen::VectorXf h;	// nsamples (LOCO HI for this chr)
		int m = 0;
		bool ok = false;
	};

	std::vector<ChrCache> cache((size_t)C);

	for (int c = 0; c < C; ++c) {
		const auto& chr = chr_order_scan[c];
		auto it = windows_by_chr_scan.find(chr);
		if (it == windows_by_chr_scan.end())
			continue;

		cache[(size_t)c].idx = &it->second;
		cache[(size_t)c].m = (int)it->second.size();
		if (cache[(size_t)c].m <= 0)
			continue;

		// LOCO HI excluding this chromosome (store it for binning)
		cache[(size_t)c].h = hi_from_components_weighted_excluding(hc_full, chr);

		// Residualize phenotype vs LOCO HI
		int g_valid = 0;
		try {
			cache[(size_t)c].gZ = residualize_and_zscore_vector(g_raw, cache[(size_t)c].h, g_valid);
		} catch (...) {
			continue;
		}

		// Residualize windows on this chromosome vs LOCO HI
		int n_valid = 0;
		try {
			cache[(size_t)c].Zc = residualize_and_zscore_subset(X_scan, cache[(size_t)c].h, *cache[(size_t)c].idx, n_valid);
		} catch (...) {
			continue;
		}

		cache[(size_t)c].ok = true;
	}

	int nthreads = 1;
	#ifdef ADFINDER_HAS_OPENMP
	nthreads = opt.threads;
	if (nthreads < 1) nthreads = 1;
	omp_set_dynamic(0);
	omp_set_num_threads(nthreads);
	#endif

	summaries_out.clear();
	summaries_out.resize((size_t)n_perm);

	std::cout << "Permutation test (excl-focus, sample-geno; HI-bin + circular-shift): " << n_perm << " replicates";
	#ifdef ADFINDER_HAS_OPENMP
	std::cout << " using " << nthreads << " threads";
	#endif
	std::cout << "\n";
	std::cout << "  hi_bins     = " << hi_bins << "\n";
	std::cout << "  sample_size = " << sample_size << "\n";

	#ifdef ADFINDER_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int rep = 0; rep < n_perm; ++rep) {
		std::mt19937_64 rng(seed + (uint64_t)rep);

		PermSummary s;
		s.min_r = std::numeric_limits<float>::infinity();
		s.max_r = -std::numeric_limits<float>::infinity();

		std::vector<float> sample;
		sample.reserve((size_t)sample_size);
		long long seen = 0;

		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::RowVectorXf R(tile_size);

		for (int c = 0; c < C; ++c) {
			const auto& cc = cache[(size_t)c];
			if (!cc.ok || cc.m <= 0)
				continue;

			const Eigen::VectorXf& gZ = cc.gZ;
			const Eigen::MatrixXf& Zc = cc.Zc;
			const Eigen::VectorXf& h = cc.h;
			const int m = cc.m;

			// Build HI bins for this chromosome's LOCO HI
			std::vector<std::vector<int>> bin_members((size_t)hi_bins);
			for (int i = 0; i < nsamples; ++i) {
				int b = hi_bin_index(h(i), hi_bins, 0.0f, 1.0f);
				if (b >= 0)
					bin_members[(size_t)b].push_back(i);
			}

			// perm[i] = source sample row used for output-row i (within-bin permutation)
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

			// circular shift per *source* individual in rank space within this chromosome
			std::vector<int> sh(nsamples, 0);
			if (m > 1) {
				std::uniform_int_distribution<int> U(0, m - 1);
				for (int src = 0; src < nsamples; ++src)
					sh[src] = U(rng);
			}

			for (int j0 = 0; j0 < m; j0 += tile_size) {
				int b2 = std::min(tile_size, m - j0);

				for (int k = 0; k < b2; ++k) {
					int rrank = j0 + k;	// rank within chr
					for (int i = 0; i < nsamples; ++i) {
						int src = perm[i];
						int rr = rrank;
						if (m > 1)
							rr = (rrank + sh[src]) % m;
						B(i, k) = Zc(src, rr);
					}
				}

				R.head(b2).noalias() = gZ.transpose() * B.leftCols(b2);
				R.head(b2) *= denom;

				for (int k = 0; k < b2; ++k) {
					float r = R(k);
					if (!std::isfinite(r))
						continue;

					if (r < s.min_r) s.min_r = r;
					if (r > s.max_r) s.max_r = r;

					++seen;
					if ((int)sample.size() < sample_size) {
						sample.push_back(r);
					} else {
						std::uniform_int_distribution<long long> U(0, seen - 1);
						long long jj = U(rng);
						if (jj < sample_size)
							sample[(size_t)jj] = r;
					}
				}
			}
		}

		// mean + sd from the same reservoir sample used for quantiles
		double mu = 0.0;
		double m2 = 0.0;
		long long n = 0;
		for (float x : sample) {
			if (!std::isfinite(x))
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
