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

// Compute exact block lengths from MSP spos/epos columns.
static std::vector<float> compute_block_lengths_from_range_msp(
	const std::vector<int>& pos_start,
	const std::vector<int>& pos_end
) {
	int nwin = (int)pos_start.size();
	std::vector<float> len((size_t)nwin, 1.0f);
	for (int w = 0; w < nwin; ++w) {
		float L = (float)(pos_end[w] - pos_start[w]);
		if (L < 1.0f) L = 1.0f;
		len[(size_t)w] = L;
	}
	return len;
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
	bool pos_is_start,
	const std::vector<int>& pos_start
) {
	int nsamples = (int)X.rows();
	int nwin = (int)X.cols();

	Eigen::VectorXf h(nsamples);

	if ((int)chroms.size() != nwin || (int)pos.size() != nwin) {
		std::cerr << "Warning: weighted HI requested but chrom/pos size mismatch; falling back to unweighted HI.\n";
		return compute_hi_from_X(X);
	}

	std::vector<float> wlen = (!pos_start.empty() && (int)pos_start.size() == nwin)
		? compute_block_lengths_from_range_msp(pos_start, pos)
		: infer_block_lengths(chroms, pos, pos_is_start);

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

HiComponentsWeighted build_hi_components_weighted(
	const Eigen::MatrixXf& X,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	bool pos_is_start,
	const std::vector<int>& pos_start
) {
	const int nsamples = (int)X.rows();
	const int nwin = (int)X.cols();

	HiComponentsWeighted hc;
	hc.wsum_total = Eigen::VectorXf::Zero(nsamples);
	hc.wtot_total = Eigen::VectorXf::Zero(nsamples);

	if ((int)chroms.size() != nwin || (int)pos.size() != nwin) {
		std::cerr << "Warning: build_hi_components_weighted chrom/pos size mismatch; using length=1 for all windows.\n";
	}

	// Stable chromosome order: first occurrence in chroms[]
	std::unordered_map<std::string, int> chr_to_idx;
	chr_to_idx.reserve(64);

	for (int w = 0; w < nwin; ++w) {
		const std::string& c = (w < (int)chroms.size()) ? chroms[w] : std::string("");
		if (c.empty())
			continue;

		auto it = chr_to_idx.find(c);
		if (it == chr_to_idx.end()) {
			int id = (int)hc.chr_order.size();
			chr_to_idx[c] = id;
			hc.chr_order.push_back(c);
		}
	}

	const int C = (int)hc.chr_order.size();
	hc.wsum_chr.resize((size_t)C);
	hc.wtot_chr.resize((size_t)C);
	for (int c = 0; c < C; ++c) {
		hc.wsum_chr[(size_t)c] = Eigen::VectorXf::Zero(nsamples);
		hc.wtot_chr[(size_t)c] = Eigen::VectorXf::Zero(nsamples);
	}

	// Weights per window (block lengths)
	std::vector<float> wlen((size_t)nwin, 1.0f);
	if ((int)chroms.size() == nwin && (int)pos.size() == nwin) {
		wlen = (!pos_start.empty() && (int)pos_start.size() == nwin)
			? compute_block_lengths_from_range_msp(pos_start, pos)
			: infer_block_lengths(chroms, pos, pos_is_start);
	}

	// Accumulate components
	for (int w = 0; w < nwin; ++w) {
		if (w >= (int)chroms.size())
			continue;

		const std::string& chr = chroms[w];
		auto it = chr_to_idx.find(chr);
		if (it == chr_to_idx.end())
			continue;

		const int cidx = it->second;
		const float Lf = (w < (int)wlen.size()) ? wlen[(size_t)w] : 1.0f;
		const float L = (Lf >= 1.0f) ? Lf : 1.0f;

		for (int i = 0; i < nsamples; ++i) {
			const float v = X(i, w);
			if (std::isnan(v))
				continue;

			hc.wsum_total(i) += v * L;
			hc.wtot_total(i) += L;

			hc.wsum_chr[(size_t)cidx](i) += v * L;
			hc.wtot_chr[(size_t)cidx](i) += L;
		}
	}

	return hc;
}

Eigen::VectorXf hi_from_components_weighted_excluding(
	const HiComponentsWeighted& hc,
	const std::string& exclude_chrA,
	const std::string& exclude_chrB
) {
	const int nsamples = (int)hc.wsum_total.size();

	Eigen::VectorXf wsum = hc.wsum_total;
	Eigen::VectorXf wtot = hc.wtot_total;

	auto subtract_chr = [&](const std::string& chr) {
		if (chr.empty())
			return;

		for (int c = 0; c < (int)hc.chr_order.size(); ++c) {
			if (hc.chr_order[(size_t)c] == chr) {
				wsum -= hc.wsum_chr[(size_t)c];
				wtot -= hc.wtot_chr[(size_t)c];
				return;
			}
		}
	};

	subtract_chr(exclude_chrA);
	if (!exclude_chrB.empty() && exclude_chrB != exclude_chrA)
		subtract_chr(exclude_chrB);

	Eigen::VectorXf h(nsamples);
	for (int i = 0; i < nsamples; ++i) {
		const float d = wtot(i);
		if (!(d > 0.0f) || !std::isfinite(d)) {
			h(i) = std::numeric_limits<float>::quiet_NaN();
		} else {
			h(i) = (wsum(i) / d) / 2.0f;
		}
	}

	return h;
}

bool load_cov_tsv(
	const std::string& path,
	const std::vector<std::string>& sample_names,
	Eigen::MatrixXf& H_out
) {
	std::ifstream in(path);
	if (!in) {
		std::cerr << "Error: cannot open covariate file: " << path << "\n";
		return false;
	}

	int nsamples = (int)sample_names.size();
	int ncols = -1;

	std::unordered_map<std::string, std::vector<float>> m;
	m.reserve((size_t)nsamples * 2);

	std::string line;
	while (std::getline(in, line)) {
		if (line.empty() || line[0] == '#') continue;

		std::istringstream ss(line);
		std::string sample;
		if (!(ss >> sample)) continue;

		std::vector<float> vals;
		std::string tok;
		while (ss >> tok) {
			try {
				vals.push_back(std::stof(tok));
			} catch (...) {
				vals.clear();
				break;  // non-numeric token: header line
			}
		}
		if (vals.empty()) continue;

		if (ncols == -1) {
			ncols = (int)vals.size();
		} else if ((int)vals.size() != ncols) {
			std::cerr << "Error: covariate file has inconsistent column count at sample: "
			          << sample << " (expected " << ncols << ", got " << vals.size() << ")\n";
			return false;
		}

		m[sample] = std::move(vals);
	}

	if (ncols <= 0) {
		std::cerr << "Error: no data rows found in covariate file: " << path << "\n";
		return false;
	}

	H_out.resize(nsamples, ncols);

	int missing = 0;
	for (int i = 0; i < nsamples; ++i) {
		auto it = m.find(sample_names[i]);
		if (it == m.end()) {
			++missing;
		} else {
			for (int j = 0; j < ncols; ++j)
				H_out(i, j) = it->second[j];
		}
	}

	if (missing > 0) {
		std::cerr << "Error: covariate file missing " << missing << " / " << nsamples
		          << " samples (must contain all samples).\n";
		return false;
	}

	std::cout << "Loaded covariate matrix: " << nsamples << " samples x " << ncols
	          << " covariate" << (ncols > 1 ? "s" : "") << " from: " << path << "\n";
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
