#include "scan_blocks.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <iomanip>

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
	std::ofstream of(out_path);
	if (!of) {
		std::cerr << "Error: cannot write to " << out_path << "\n";
		return false;
	}

	int nsamples = opt.nsamples;
	int block_size = opt.block_size;
	float min_abs_r = opt.min_abs_r;

	std::cout << "[scan] intra=" << (opt.intra ? 1 : 0)
		  << " max_dist=" << opt.max_dist
		  << " block_size=" << opt.block_size
		  << "\n";

	of << "wA\tchrA\tposA\twB\tchrB\tposB\tr\tn\n";

	tested_pairs = 0;
	kept_pairs = 0;

	bool do_distrib = (!distrib_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	ScanSummary ds;
	ds.tested_pairs = 0;
	ds.min_r = std::numeric_limits<float>::infinity();
	ds.max_r = -std::numeric_limits<float>::infinity();

	std::mt19937_64 drng(distrib_seed);
	std::vector<float> dsample;
	if (do_distrib) {
		dsample.reserve((size_t)distrib_sample);
	}

	auto consider_distrib = [&](float r) {
		if (!do_distrib)
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

	const float denom = 1.0f / (float)(nsamples - 1);

	Eigen::MatrixXf A(nsamples, block_size);
	Eigen::MatrixXf B(nsamples, block_size);
	Eigen::MatrixXf R(block_size, block_size);

	auto emit_hit = [&](int a, int b, float r) {
		of << a << "\t" << chroms[a] << "\t" << pos[a] << "\t"
		<< b << "\t" << chroms[b] << "\t" << pos[b] << "\t"
		<< r << "\t" << nsamples << "\n";
		++kept_pairs;
	};


	if (opt.intra) {
		for (const auto& chr : chr_order) {
			const auto& idx = windows_by_chr.at(chr);
			int m = (int)idx.size();

			for (int i0 = 0; i0 < m; i0 += block_size) {
				int b1 = std::min(block_size, m - i0);

				for (int k = 0; k < b1; ++k)
					A.col(k) = Z.col(idx[i0 + k]);

				for (int j0 = i0; j0 < m; j0 += block_size) {
					int b2 = std::min(block_size, m - j0);

					for (int k = 0; k < b2; ++k)
						B.col(k) = Z.col(idx[j0 + k]);

					R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
					R.topLeftCorner(b1, b2) *= denom;

					for (int ia = 0; ia < b1; ++ia) {
						int a = idx[i0 + ia];
						int posA = pos[a];

						int jb_start = 0;
						if (j0 == i0)
							jb_start = ia + 1;

						int jb_end = b2;

						// Distance filter (only for intra)
						if (opt.max_dist >= 0) {
							int limit = posA + opt.max_dist;

							// idx is assumed sorted by pos[] within chromosome.
							// Search within the j-block [j0, j0+b2) for first pos > limit.
							auto begin_it = idx.begin() + j0;
							auto end_it = idx.begin() + j0 + b2;

							auto ub = std::upper_bound(begin_it, end_it, limit,
								[&](int value, int widx) {
									return value < pos[widx];
								}
							);

							jb_end = (int)std::distance(begin_it, ub);

							// jb_end is relative to begin_it; convert to [0,b2]
							// (distance already in that coordinate system)
						}

						if (jb_end <= jb_start)
							continue;

						for (int ib = jb_start; ib < jb_end; ++ib) {
							int b = idx[j0 + ib];
							float r = R(ia, ib);

							++tested_pairs;

							consider_distrib(r);

							bool keep = false;

							if (opt.use_asym) {
								if (opt.min_neg_r > 0.0f && r <= -opt.min_neg_r)
									keep = true;
								if (opt.min_pos_r > 0.0f && r >= opt.min_pos_r)
									keep = true;
							} else {
								if (std::fabs(r) >= opt.min_abs_r)
									keep = true;
							}

							if (keep)
								emit_hit(a, b, r);

						}
					}

				}
			}
		}
	} else {
		int C = (int)chr_order.size();

		for (int c1 = 0; c1 < C; ++c1) {
			const auto& chr1 = chr_order[c1];
			const auto& idx1 = windows_by_chr.at(chr1);
			int m1 = (int)idx1.size();

			for (int c2 = c1 + 1; c2 < C; ++c2) {
				const auto& chr2 = chr_order[c2];
				const auto& idx2 = windows_by_chr.at(chr2);
				int m2 = (int)idx2.size();

				for (int i0 = 0; i0 < m1; i0 += block_size) {
					int b1 = std::min(block_size, m1 - i0);

					for (int k = 0; k < b1; ++k)
						A.col(k) = Z.col(idx1[i0 + k]);

					for (int j0 = 0; j0 < m2; j0 += block_size) {
						int b2 = std::min(block_size, m2 - j0);

						for (int k = 0; k < b2; ++k)
							B.col(k) = Z.col(idx2[j0 + k]);

						R.topLeftCorner(b1, b2).noalias() = A.leftCols(b1).transpose() * B.leftCols(b2);
						R.topLeftCorner(b1, b2) *= denom;

						for (int ia = 0; ia < b1; ++ia) {
							int a = idx1[i0 + ia];
							for (int ib = 0; ib < b2; ++ib) {
								int b = idx2[j0 + ib];
								float r = R(ia, ib);

								++tested_pairs;

								consider_distrib(r);

								bool keep = false;

								if (opt.use_asym) {
									if (opt.min_neg_r > 0.0f && r <= -opt.min_neg_r)
										keep = true;
									if (opt.min_pos_r > 0.0f && r >= opt.min_pos_r)
										keep = true;
								} else {
									if (std::fabs(r) >= opt.min_abs_r)
										keep = true;
								}

								if (keep)
									emit_hit(a, b, r);
							}
						}
					}
				}
			}
		}
	}

	if (do_distrib) {
		std::ofstream df(distrib_path);
		if (!df) {
			std::cerr << "Error: cannot write to " << distrib_path << "\n";
			return false;
		}

		df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\n";

		if (ds.tested_pairs == 0 || dsample.empty()) {
			df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
		} else {
			std::sort(dsample.begin(), dsample.end());
			ds.p01 = quantile_from_sorted(dsample, 0.01);
			ds.p05 = quantile_from_sorted(dsample, 0.05);
			ds.p25 = quantile_from_sorted(dsample, 0.25);
			ds.median = quantile_from_sorted(dsample, 0.50);
			ds.p75 = quantile_from_sorted(dsample, 0.75);
			ds.p95 = quantile_from_sorted(dsample, 0.95);
			ds.p99 = quantile_from_sorted(dsample, 0.99);

			df << ds.tested_pairs
				<< "\t" << ds.max_r
				<< "\t" << ds.p99
				<< "\t" << ds.p95
				<< "\t" << ds.p75
				<< "\t" << ds.median
				<< "\t" << ds.p25
				<< "\t" << ds.p05
				<< "\t" << ds.p01
				<< "\t" << ds.min_r
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
	std::ofstream of(out_path);
	if (!of) {
		std::cerr << "Error: cannot write to " << out_path << "\n";
		return false;
	}

	bool do_distrib = (!distrib_path.empty());
	if (distrib_sample <= 0)
		distrib_sample = 200000;

	int nsamples = opt.nsamples;
	int block_size = opt.block_size;
	float min_abs_r = opt.min_abs_r;
	const float denom = 1.0f / (float)(nsamples - 1);

	ScanSummary ds;
	ds.tested_pairs = 0;
	ds.min_r = std::numeric_limits<float>::infinity();
	ds.max_r = -std::numeric_limits<float>::infinity();

	std::mt19937_64 drng(distrib_seed);
	std::vector<float> dsample;
	if (do_distrib)
		dsample.reserve((size_t)distrib_sample);

	auto consider_distrib = [&](float r) {
		if (!do_distrib)
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

	of << "wA\tchrA\tposA\twB\tchrB\tposB\tr\tn\n";

	tested_pairs = 0;
	kept_pairs = 0;

	Eigen::VectorXf v = Z.col(target_w);				// nsamples
	Eigen::MatrixXf B(nsamples, block_size);			// nsamples x b2
	Eigen::RowVectorXf R(block_size);					// 1 x b2

	const std::string& tchr = chroms[target_w];
	int tpos = pos[target_w];

	for (const auto& chr : chr_order) {
		const auto& idx = windows_by_chr.at(chr);
		int m = (int)idx.size();

		// If not intra, skip target chromosome entirely
		if (!opt.intra && chr == tchr)
			continue;

		for (int j0 = 0; j0 < m; j0 += block_size) {
			int b2 = std::min(block_size, m - j0);

			for (int k = 0; k < b2; ++k)
				B.col(k) = Z.col(idx[j0 + k]);

			R.head(b2).noalias() = v.transpose() * B.leftCols(b2);
			R.head(b2) *= denom;

			for (int ib = 0; ib < b2; ++ib) {
				int b = idx[j0 + ib];

				if (b == target_w)
					continue;

				// If intra + max_dist, enforce distance on same chr only
				if (opt.intra && opt.max_dist >= 0 && chroms[b] == tchr) {
					int d = pos[b] - tpos;
					if (d < 0) d = -d;
					if (d > opt.max_dist)
						continue;
				}

				float r = R(ib);
				++tested_pairs;
				consider_distrib(r);

				bool keep = false;

				if (opt.use_asym) {
					if (opt.min_neg_r > 0.0f && r <= -opt.min_neg_r)
						keep = true;
					if (opt.min_pos_r > 0.0f && r >= opt.min_pos_r)
						keep = true;
				} else {
					if (std::fabs(r) >= opt.min_abs_r)
						keep = true;
				}

				if (keep) {
					of << target_w << "\t" << chroms[target_w] << "\t" << pos[target_w] << "\t"
						<< b << "\t" << chroms[b] << "\t" << pos[b] << "\t"
						<< r << "\t" << nsamples << "\n";
					++kept_pairs;
				}
			}
		}
	}

	// Write distrib summary if requested
	if (do_distrib) {
		std::ofstream df(distrib_path);
		if (!df) {
			std::cerr << "Error: cannot write to " << distrib_path << "\n";
			return false;
		}

		df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\n";

		if (ds.tested_pairs == 0 || dsample.empty()) {
			df << 0 << "\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
		} else {
			std::sort(dsample.begin(), dsample.end());
			ds.p01 = quantile_from_sorted(dsample, 0.01);
			ds.p05 = quantile_from_sorted(dsample, 0.05);
			ds.p25 = quantile_from_sorted(dsample, 0.25);
			ds.median = quantile_from_sorted(dsample, 0.50);
			ds.p75 = quantile_from_sorted(dsample, 0.75);
			ds.p95 = quantile_from_sorted(dsample, 0.95);
			ds.p99 = quantile_from_sorted(dsample, 0.99);

			df << ds.tested_pairs
				<< "\t" << ds.max_r
				<< "\t" << ds.p99
				<< "\t" << ds.p95
				<< "\t" << ds.p75
				<< "\t" << ds.median
				<< "\t" << ds.p25
				<< "\t" << ds.p05
				<< "\t" << ds.p01
				<< "\t" << ds.min_r
				<< "\n";
		}
	}

	return true;
}

bool permute_interchrom_summary_chrblock(
	const Eigen::MatrixXf& Z,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	uint64_t seed,
	int n_perm,
	int sample_size,
	std::vector<PermSummary>& summaries_out
) {
	if (opt.intra) {
		std::cerr << "Error: chr-block permutation summary is implemented for interchrom scans only.\n";
		return false;
	}
	if (n_perm <= 0) {
		summaries_out.clear();
		return true;
	}
	if (sample_size <= 0)
		sample_size = 200000;

	int nsamples = opt.nsamples;
	int block_size = opt.block_size;

	summaries_out.clear();
	summaries_out.resize((size_t)n_perm);

	Eigen::MatrixXf A(nsamples, block_size);
	Eigen::MatrixXf B(nsamples, block_size);
	Eigen::MatrixXf R(block_size, block_size);

	const float denom = 1.0f / (float)(nsamples - 1);

	std::vector<int> base(nsamples);
	for (int i = 0; i < nsamples; ++i)
		base[i] = i;

	std::unordered_map<std::string, std::vector<int>> perm_by_chr;
	perm_by_chr.reserve(chr_order.size());
	for (const auto& chr : chr_order)
		perm_by_chr[chr] = std::vector<int>(nsamples);

	summaries_out.clear();
	summaries_out.resize((size_t)n_perm);

	int nthreads = 1;
	#ifdef ADFINDER_HAS_OPENMP
	nthreads = omp_get_max_threads();
	#endif

	std::cout << "Permutation test: " << n_perm << " replicates";
	#ifdef ADFINDER_HAS_OPENMP
	std::cout << " using " << nthreads << " threads";
	#endif
	std::cout << "\n";

	#pragma omp parallel for schedule(dynamic)
	for (int rep = 0; rep < n_perm; ++rep) {
		// Thread-local RNG
		std::mt19937_64 rng(seed + (uint64_t)rep);

		// Thread-local buffers (avoid races)
		Eigen::MatrixXf A(nsamples, block_size);
		Eigen::MatrixXf B(nsamples, block_size);
		Eigen::MatrixXf R(block_size, block_size);

		// Base vector + perm maps (thread-local)
		std::vector<int> base(nsamples);
		for (int i = 0; i < nsamples; ++i)
			base[i] = i;

		std::unordered_map<std::string, std::vector<int>> perm_by_chr;
		perm_by_chr.reserve(chr_order.size());
		for (const auto& chr : chr_order)
			perm_by_chr[chr] = base;

		// Build one permutation per chromosome (for this rep)
		for (const auto& chr : chr_order) {
			auto& perm = perm_by_chr[chr];
			std::shuffle(perm.begin(), perm.end(), rng);
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
			const auto& p1 = perm_by_chr.at(chr1);
			int m1 = (int)idx1.size();

			for (int c2 = c1 + 1; c2 < C; ++c2) {
				const auto& chr2 = chr_order[c2];
				const auto& idx2 = windows_by_chr.at(chr2);
				const auto& p2 = perm_by_chr.at(chr2);
				int m2 = (int)idx2.size();

				for (int i0 = 0; i0 < m1; i0 += block_size) {
					int b1 = std::min(block_size, m1 - i0);

					for (int k = 0; k < b1; ++k) {
						int w = idx1[i0 + k];
						for (int i = 0; i < nsamples; ++i)
							A(i, k) = Z(p1[i], w);
					}

					for (int j0 = 0; j0 < m2; j0 += block_size) {
						int b2 = std::min(block_size, m2 - j0);

						for (int k = 0; k < b2; ++k) {
							int w = idx2[j0 + k];
							for (int i = 0; i < nsamples; ++i)
								B(i, k) = Z(p2[i], w);
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

	// output header
	of << "tag\twB\tchrB\tposB\tr\tn\n";

	// optional distribution summary (reservoir sample)
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

	// Scan in chr blocks for cache friendliness, but we compute one r per window.
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
		std::sort(dsample.begin(), dsample.end());

		float p01 = quantile_from_sorted(dsample, 0.01);
		float p05 = quantile_from_sorted(dsample, 0.05);
		float p25 = quantile_from_sorted(dsample, 0.25);
		float med = quantile_from_sorted(dsample, 0.50);
		float p75 = quantile_from_sorted(dsample, 0.75);
		float p95 = quantile_from_sorted(dsample, 0.95);
		float p99 = quantile_from_sorted(dsample, 0.99);

		std::ofstream df(distrib_path);
		if (!df) {
			std::cerr << "Error: cannot write distrib summary to " << distrib_path << "\n";
			return false;
		}

		df << "tested_pairs\tmax_r\tp99\tp95\tp75\tmedian\tp25\tp05\tp01\tmin_r\n";
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
			<< "\n";
	}


	return true;
}


bool permute_sample_vector_summary(
	const Eigen::MatrixXf& Z,
	const Eigen::VectorXf& gZ,
	const ScanOptions& opt,
	uint64_t seed,
	int n_perm,
	int sample_size,
	std::vector<PermSummary>& summaries_out
) {
	if (n_perm <= 0) {
		summaries_out.clear();
		return true;
	}

	const int nsamples = opt.nsamples;
	const int nwin = (int)Z.cols();
	const float denom = 1.0f / (float)(nsamples - 1);

	if (gZ.size() != nsamples) {
		std::cerr << "Error: gZ length does not match nsamples\n";
		return false;
	}

	if (sample_size <= 0)
		sample_size = 200000;

	summaries_out.clear();
	summaries_out.resize((size_t)n_perm);

	#ifdef ADFINDER_HAS_OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int rep = 0; rep < n_perm; ++rep) {
		std::mt19937_64 rng(seed + (uint64_t)rep);

		// permutation of sample indices
		std::vector<int> perm(nsamples);
		for (int i = 0; i < nsamples; ++i)
			perm[i] = i;
		std::shuffle(perm.begin(), perm.end(), rng);

		PermSummary s;
		s.min_r = std::numeric_limits<float>::infinity();
		s.max_r = -std::numeric_limits<float>::infinity();

		std::vector<float> sample;
		sample.reserve((size_t)sample_size);

		long long seen = 0;

		// loop over windows
		for (int w = 0; w < nwin; ++w) {
			double dot = 0.0;

			for (int i = 0; i < nsamples; ++i) {
				float zi = Z(i, w);
				float gi = gZ(perm[i]);

				if (!std::isnan(zi) && !std::isnan(gi))
					dot += (double)zi * (double)gi;
			}

			float r = (float)(dot * denom);

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
