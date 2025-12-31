#include "scan_blocks.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <iomanip>

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

bool scan_blocks_write_hits(
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& starts,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs
) {
	std::ofstream of(out_path);
	if (!of) {
		std::cerr << "Error: cannot write to " << out_path << "\n";
		return false;
	}

	int nsamples = opt.nsamples;
	int block_size = opt.block_size;
	float min_abs_r = opt.min_abs_r;

	of << "wA\tchrA\tstartA\twB\tchrB\tstartB\tr\tn\n";

	tested_pairs = 0;
	kept_pairs = 0;

	const float denom = 1.0f / (float)(nsamples - 1);

	Eigen::MatrixXf A(nsamples, block_size);
	Eigen::MatrixXf B(nsamples, block_size);
	Eigen::MatrixXf R(block_size, block_size);

	auto emit_hit = [&](int a, int b, float r) {
		of << a << "\t" << chroms[a] << "\t" << starts[a] << "\t"
		   << b << "\t" << chroms[b] << "\t" << starts[b] << "\t"
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
						int jb_start = (j0 == i0) ? (ia + 1) : 0;

						for (int ib = jb_start; ib < b2; ++ib) {
							int b = idx[j0 + ib];
							float r = R(ia, ib);

							++tested_pairs;

							if (std::fabs(r) >= min_abs_r)
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

								if (std::fabs(r) >= min_abs_r)
									emit_hit(a, b, r);
							}
						}
					}
				}
			}
		}
	}

	return true;
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
	summaries_out.reserve((size_t)n_perm);

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

	auto t0 = std::chrono::steady_clock::now();
	auto last = t0;

	int every = 1;
	if (n_perm >= 1000) every = 10;
	else if (n_perm >= 200) every = 5;
	else if (n_perm >= 50) every = 2;

	std::cout << "Permuting (interchrom chr-block): " << n_perm << " replicates\n";
	for (int rep = 0; rep < n_perm; ++rep) {
		std::mt19937_64 rng(seed + (uint64_t)rep);

		for (const auto& chr : chr_order) {
			auto& perm = perm_by_chr[chr];
			perm = base;
			std::shuffle(perm.begin(), perm.end(), rng);
		}

		PermSummary s;
		s.min_r = std::numeric_limits<float>::infinity();
		s.max_r = -std::numeric_limits<float>::infinity();

		// Reservoir sample of r values (signed)
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

								// Reservoir sampling
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
		s.median = quantile_from_sorted(sample, 0.50);
		s.p95 = quantile_from_sorted(sample, 0.95);
		s.p99 = quantile_from_sorted(sample, 0.99);

		summaries_out.push_back(s);

	if ((rep + 1) % every == 0 || (rep + 1) == n_perm) {
		auto now = std::chrono::steady_clock::now();
		double elapsed = std::chrono::duration<double>(now - t0).count();
		double since_last = std::chrono::duration<double>(now - last).count();

		// Also print if it's been a while since last update
		if (since_last >= 5.0 || (rep + 1) == n_perm) {
			double rate = (elapsed > 0.0) ? (double)(rep + 1) / elapsed : 0.0;
			double eta = (rate > 0.0) ? ((double)n_perm - (rep + 1)) / rate : 0.0;

			std::cout << "  perm " << (rep + 1) << "/" << n_perm
					<< "  elapsed=" << std::fixed << std::setprecision(1) << elapsed << "s"
					<< "  rate=" << std::setprecision(2) << rate << "/s"
					<< "  ETA=" << std::setprecision(1) << eta << "s"
					<< "  max_r=" << std::setprecision(4) << s.max_r
					<< "\n";

			last = now;
		}
	}

	}

	return true;
}
