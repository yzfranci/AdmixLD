#include "hybrid_index.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <algorithm>

Eigen::VectorXf compute_hi_from_X(
	const Eigen::MatrixXf& X
) {
	int nsamples = (int)X.rows();
	int nwin = (int)X.cols();

	Eigen::VectorXf h(nsamples);

	for (int i = 0; i < nsamples; ++i) {
		double sum = 0.0;
		int count = 0;

		for (int w = 0; w < nwin; ++w) {
			float v = X(i, w);
			if (!std::isnan(v)) {
				sum += (double)v;
				++count;
			}
		}

		if (count == 0) {
			h(i) = std::numeric_limits<float>::quiet_NaN();
		} else {
			h(i) = (float)((sum / count) / 2.0);
		}
	}

	return h;
}

static std::vector<float> infer_block_lengths(
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	bool pos_is_start
) {
	int nwin = (int)chroms.size();
	std::vector<float> len((size_t)nwin, 1.0f);
	if ((int)pos.size() != nwin)
		return len;

	// group indices by chr
	std::unordered_map<std::string, std::vector<int>> idx_by_chr;
	idx_by_chr.reserve(64);

	for (int w = 0; w < nwin; ++w)
		idx_by_chr[chroms[w]].push_back(w);

	for (auto& kv : idx_by_chr) {
		auto& idx = kv.second;

		// sort by pos ascending
		std::sort(idx.begin(), idx.end(),
			[&](int a, int b) { return pos[a] < pos[b]; }
		);

		int m = (int)idx.size();
		if (m <= 0)
			continue;

		if (!pos_is_start) {
			// pos = END
			{
				int w0 = idx[0];
				float L = (float)pos[w0];
				if (L < 1.0f) L = 1.0f;
				len[(size_t)w0] = L;
			}

			for (int j = 1; j < m; ++j) {
				int w_prev = idx[j - 1];
				int w = idx[j];
				float L = (float)(pos[w] - pos[w_prev]);
				if (L < 1.0f) L = 1.0f;
				len[(size_t)w] = L;
			}

		} else {
			// pos = START
			float last_L = 1.0f;

			for (int j = 0; j < m - 1; ++j) {
				int w = idx[j];
				int w_next = idx[j + 1];
				float L = (float)(pos[w_next] - pos[w]);
				if (L < 1.0f) L = 1.0f;
				len[(size_t)w] = L;
				last_L = L;
			}

			// last window: reuse last inferred length (or 1)
			int w_last = idx[m - 1];
			if (last_L < 1.0f) last_L = 1.0f;
			len[(size_t)w_last] = last_L;
		}
	}

	return len;
}

Eigen::VectorXf compute_hi_from_X_weighted(
	const Eigen::MatrixXf& X,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	bool pos_is_start
) {
	int nsamples = (int)X.rows();
	int nwin = (int)X.cols();

	Eigen::VectorXf h(nsamples);

	if ((int)chroms.size() != nwin || (int)pos.size() != nwin) {
		std::cerr << "Warning: weighted HI requested but chrom/pos size mismatch; falling back to unweighted HI.\n";
		return compute_hi_from_X(X);
	}

	std::vector<float> wlen = infer_block_lengths(chroms, pos, pos_is_start);

	for (int i = 0; i < nsamples; ++i) {
		double wsum = 0.0;
		double wtot = 0.0;

		for (int w = 0; w < nwin; ++w) {
			float v = X(i, w);
			if (!std::isnan(v)) {
				double L = (double)wlen[(size_t)w];
				wsum += (double)v * L;
				wtot += L;
			}
		}

		if (wtot <= 0.0) {
			h(i) = std::numeric_limits<float>::quiet_NaN();
		} else {
			h(i) = (float)((wsum / wtot) / 2.0);
		}
	}

	return h;
}

bool load_hi_tsv(
	const std::string& path,
	const std::vector<std::string>& sample_names,
	Eigen::VectorXf& h_out
) {
	std::ifstream in(path);
	if (!in) {
		std::cerr << "Error: cannot open HI file: " << path << "\n";
		return false;
	}

	int nsamples = (int)sample_names.size();

	std::unordered_map<std::string, float> m;
	m.reserve((size_t)nsamples * 2);

	std::string line;
	while (std::getline(in, line)) {
		if (line.size() == 0) continue;
		if (line[0] == '#') continue;

		std::istringstream ss(line);
		std::string sample;
		std::string hi_str;

		if (!(ss >> sample >> hi_str))
			continue;

		try {
			float hi = std::stof(hi_str);
			m[sample] = hi;
		} catch (...) {
			// Allows header lines like "sample hi"
			continue;
		}
	}

	h_out.resize(nsamples);

	int missing = 0;
	for (int i = 0; i < nsamples; ++i) {
		auto it = m.find(sample_names[i]);
		if (it == m.end()) {
			h_out(i) = std::numeric_limits<float>::quiet_NaN();
			++missing;
		} else {
			h_out(i) = it->second;
		}
	}

	if (missing > 0) {
		std::cerr << "Error: HI file missing " << missing << " / " << nsamples
				  << " samples (must contain all samples).\n";
		return false;
	}

	return true;
}

bool write_hi_tsv(
	const std::string& path,
	const std::vector<std::string>& sample_names,
	const Eigen::VectorXf& h
) {
	std::ofstream out(path);
	if (!out) {
		std::cerr << "Error: cannot write HI file: " << path << "\n";
		return false;
	}

	int nsamples = (int)sample_names.size();
	out << "sample\thi\n";
	for (int i = 0; i < nsamples; ++i) {
		out << sample_names[i] << "\t" << h(i) << "\n";
	}

	return true;
}
