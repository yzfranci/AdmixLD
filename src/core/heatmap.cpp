#include "heatmap.hpp"
#include "residualize.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>

#ifdef ADMIXLD_HAS_OPENMP
#include <omp.h>
#endif

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

static inline uint64_t pack_bins(int binA, int binB) {
	return ((uint64_t)(uint32_t)binA << 32) | (uint64_t)(uint32_t)binB;
}

// Per (chrA-bin, chrB-bin) cell accumulator: exact running mean/SD (Welford)
// plus a bounded reservoir (Algorithm R) for an approximate upper quantile.
// Bounding the reservoir per cell (rather than storing every r that lands in
// the cell) is what keeps memory bounded regardless of scan size.
struct HeatmapCell {
	long long n = 0;
	double mean = 0.0;
	double m2 = 0.0;
	std::vector<float> reservoir;
};

static void heatmap_cell_add(HeatmapCell& c, float r, int reservoir_size, std::mt19937_64& rng) {
	++c.n;
	double delta = (double)r - c.mean;
	c.mean += delta / (double)c.n;
	double delta2 = (double)r - c.mean;
	c.m2 += delta * delta2;

	if (reservoir_size <= 0)
		return;

	if ((int)c.reservoir.size() < reservoir_size) {
		c.reservoir.push_back(r);
	} else {
		std::uniform_int_distribution<long long> U(0, c.n - 1);
		long long j = U(rng);
		if (j < reservoir_size)
			c.reservoir[(size_t)j] = r;
	}
}

// Writes every non-empty cell of one block (one chromosome for intra, one
// chromosome pair for inter) directly to that thread's part file, then the
// caller drops the map — memory is bounded by the cells touched by the
// current block, not by the whole scan.
static void write_heatmap_block(
	std::ostream& out,
	const std::string& chrA,
	const std::string& chrB,
	long long bin_size,
	double quantile,
	std::unordered_map<uint64_t, HeatmapCell>& cells
) {
	for (auto& kv : cells) {
		const int binA = (int)(kv.first >> 32);
		const int binB = (int)(uint32_t)(kv.first & 0xFFFFFFFFULL);
		HeatmapCell& c = kv.second;

		const float mean_r = (c.n > 0) ? (float)c.mean : std::numeric_limits<float>::quiet_NaN();
		const float sd_r = (c.n > 1) ? (float)std::sqrt(c.m2 / (double)(c.n - 1)) : std::numeric_limits<float>::quiet_NaN();

		std::sort(c.reservoir.begin(), c.reservoir.end());
		const float q_r = quantile_from_sorted(c.reservoir, quantile);

		const long long binA_start = (long long)binA * bin_size;
		const long long binB_start = (long long)binB * bin_size;

		out << chrA << "\t" << binA_start << "\t" << (binA_start + bin_size) << "\t"
			<< chrB << "\t" << binB_start << "\t" << (binB_start + bin_size) << "\t"
			<< c.n << "\t" << mean_r << "\t" << sd_r << "\t" << q_r << "\n";
	}
}

static void write_heatmap_header(std::ostream& out) {
	out << "chrA\tbinA_start\tbinA_end\tchrB\tbinB_start\tbinB_end\tn_pairs\tmean_r\tsd_r\tquantile_r\n";
}

bool scan_markers_write_heatmap(
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const HeatmapOptions& opt,
	const std::string& out_path,
	long long& tested_pairs
) {
	(void)chroms;

	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;
	const long long bin_size = opt.bin_size;

	tested_pairs = 0;

	std::cout << "[heatmap] intra=" << (opt.intra ? 1 : 0)
		  << " max_dist=" << opt.max_dist
		  << " bin_size=" << bin_size
		  << " quantile=" << opt.quantile
		  << " reservoir=" << opt.reservoir_size
		  << "\n";

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

	const float denom = 1.0f / (float)(nsamples - 1);

	auto bin_of = [&](int p) -> int { return (int)(p / bin_size); };

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

			Eigen::MatrixXf A(nsamples, tile_size);
			Eigen::MatrixXf B(nsamples, tile_size);
			Eigen::MatrixXf R(tile_size, tile_size);

			long long tested_local = 0;
			std::unordered_map<uint64_t, HeatmapCell> cells;
			std::mt19937_64 rng(opt.seed ^ 0xA17E4A70000ULL ^ (uint64_t)c);

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

						const int binA = bin_of(posA);

						for (int ib = jb_start; ib < jb_end; ++ib) {
							const int b = idx[j0 + ib];
							const int posB = pos[b];

							if (opt.min_dist >= 0 && (posB - posA) < opt.min_dist)
								continue;

							const float r = R(ia, ib);
							if (!std::isfinite(r))
								continue;

							++tested_local;

							const int binB = bin_of(posB);
							heatmap_cell_add(cells[pack_bins(binA, binB)], r, opt.reservoir_size, rng);
						}
					}
				}
			}

			tested_t[(size_t)tid] += tested_local;

			std::ofstream ofp(part_paths[(size_t)tid], std::ios::out | std::ios::app);
			if (!ofp)
				std::cerr << "Error: cannot write to " << part_paths[(size_t)tid] << "\n";
			else
				write_heatmap_block(ofp, chr, chr, bin_size, opt.quantile, cells);
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

			Eigen::MatrixXf A(nsamples, tile_size);
			Eigen::MatrixXf B(nsamples, tile_size);
			Eigen::MatrixXf R(tile_size, tile_size);

			long long tested_local = 0;
			std::unordered_map<uint64_t, HeatmapCell> cells;
			std::mt19937_64 rng(opt.seed ^ 0xA17E4A70000ULL ^ (uint64_t)j);

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
						const int binA = bin_of(pos[a]);

						for (int ib = 0; ib < b2; ++ib) {
							const int b = idx2[j0 + ib];
							const float r = R(ia, ib);
							if (!std::isfinite(r))
								continue;

							++tested_local;

							const int binB = bin_of(pos[b]);
							heatmap_cell_add(cells[pack_bins(binA, binB)], r, opt.reservoir_size, rng);
						}
					}
				}
			}

			tested_t[(size_t)tid] += tested_local;

			std::ofstream ofp(part_paths[(size_t)tid], std::ios::out | std::ios::app);
			if (!ofp)
				std::cerr << "Error: cannot write to " << part_paths[(size_t)tid] << "\n";
			else
				write_heatmap_block(ofp, chr1, chr2, bin_size, opt.quantile, cells);
		}
	}

	for (int t = 0; t < nthreads; ++t)
		tested_pairs += tested_t[(size_t)t];

	{
		std::ofstream of(out_path);
		if (!of) {
			std::cerr << "Error: cannot write to " << out_path << "\n";
			return false;
		}
		write_heatmap_header(of);

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

	return true;
}

template<typename HC, typename HiExcludeFn>
static bool scan_markers_write_heatmap_excl_focus_T(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HC& hc_full,
	HiExcludeFn hi_excluding,
	const HeatmapOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	const std::vector<MarkerFreq>* freqs_scan = nullptr
) {
	(void)chroms_scan;

	const int nsamples = opt.nsamples;
	const int tile_size = opt.tile_size;
	const long long bin_size = opt.bin_size;

	tested_pairs = 0;

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

	const float denom = 1.0f / (float)(nsamples - 1);

	auto bin_of = [&](int p) -> int { return (int)(p / bin_size); };

	auto do_intra_chr = [&](int tid, int block_salt, const std::string& chr) {
		const auto& idx = windows_by_chr_scan.at(chr);
		const int m = (int)idx.size();
		if (m < 2)
			return;

		Eigen::VectorXf h = hi_excluding(hc_full, chr, std::string(""));

		int n_valid = 0;
		Eigen::MatrixXf Zc = residualize_and_zscore_subset(X_scan, h, idx, n_valid, freqs_scan);

		Eigen::MatrixXf A(nsamples, tile_size);
		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::MatrixXf R(tile_size, tile_size);

		long long tested_local = 0;
		std::unordered_map<uint64_t, HeatmapCell> cells;
		std::mt19937_64 rng(opt.seed ^ 0xA17E4A70000ULL ^ (uint64_t)block_salt);

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

					const int binA = bin_of(posA);

					for (int ib = jb_start; ib < jb_end; ++ib) {
						const int b = idx[j0 + ib];
						const int posB = pos_scan[b];

						if (opt.min_dist >= 0 && (posB - posA) < opt.min_dist)
							continue;

						const float r = R(ia, ib);
						if (!std::isfinite(r))
							continue;

						++tested_local;

						const int binB = bin_of(posB);
						heatmap_cell_add(cells[pack_bins(binA, binB)], r, opt.reservoir_size, rng);
					}
				}
			}
		}

		tested_t[(size_t)tid] += tested_local;

		std::ofstream ofp(part_paths[(size_t)tid], std::ios::out | std::ios::app);
		if (!ofp)
			std::cerr << "Error: cannot write to " << part_paths[(size_t)tid] << "\n";
		else
			write_heatmap_block(ofp, chr, chr, bin_size, opt.quantile, cells);
	};

	auto do_inter_pair = [&](int tid, int block_salt, const std::string& chr1, const std::string& chr2) {
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

		Eigen::MatrixXf A(nsamples, tile_size);
		Eigen::MatrixXf B(nsamples, tile_size);
		Eigen::MatrixXf R(tile_size, tile_size);

		long long tested_local = 0;
		std::unordered_map<uint64_t, HeatmapCell> cells;
		std::mt19937_64 rng(opt.seed ^ 0xA17E4A70000ULL ^ (uint64_t)block_salt);

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
					const int binA = bin_of(pos_scan[a]);

					for (int ib = 0; ib < b2; ++ib) {
						const int b = idx2[j0 + ib];
						const float r = R(ia, ib);
						if (!std::isfinite(r))
							continue;

						++tested_local;

						const int binB = bin_of(pos_scan[b]);
						heatmap_cell_add(cells[pack_bins(binA, binB)], r, opt.reservoir_size, rng);
					}
				}
			}
		}

		tested_t[(size_t)tid] += tested_local;

		std::ofstream ofp(part_paths[(size_t)tid], std::ios::out | std::ios::app);
		if (!ofp)
			std::cerr << "Error: cannot write to " << part_paths[(size_t)tid] << "\n";
		else
			write_heatmap_block(ofp, chr1, chr2, bin_size, opt.quantile, cells);
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

			do_intra_chr(tid, c, chr_order_scan[c]);
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

			const int c1 = jobs[(size_t)j].first;
			const int c2 = jobs[(size_t)j].second;
			do_inter_pair(tid, j, chr_order_scan[c1], chr_order_scan[c2]);
		}
	}

	for (int t = 0; t < nthreads; ++t)
		tested_pairs += tested_t[(size_t)t];

	{
		std::ofstream of(out_path);
		if (!of) {
			std::cerr << "Error: cannot write to " << out_path << "\n";
			return false;
		}
		write_heatmap_header(of);

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

	return true;
}

bool scan_markers_write_heatmap_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const HeatmapOptions& opt,
	const std::string& out_path,
	long long& tested_pairs
) {
	auto fn = [](const HiComponentsWeighted& hc, const std::string& a, const std::string& b) {
		return hi_from_components_weighted_excluding(hc, a, b);
	};
	return scan_markers_write_heatmap_excl_focus_T(
		X_scan, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, out_path, tested_pairs
	);
}

bool scan_markers_write_heatmap_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const HeatmapOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	const std::vector<MarkerFreq>& freqs_scan
) {
	auto fn = [](const HiComponentsFreq& hc, const std::string& a, const std::string& b) {
		return hi_from_components_freq_excluding(hc, a, b);
	};
	return scan_markers_write_heatmap_excl_focus_T(
		X_scan, chroms_scan, pos_scan, windows_by_chr_scan, chr_order_scan,
		hc_full, fn, opt, out_path, tested_pairs, &freqs_scan
	);
}
