#include "scan_blocks.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

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
