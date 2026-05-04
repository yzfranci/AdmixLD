#include "ref_freq.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

bool load_ref_freq_map(
	const std::string& path,
	std::unordered_map<std::string, MarkerFreq>& out_map
) {
	std::ifstream in(path);
	if (!in) {
		std::cerr << "Error: cannot open --ref-freq file: " << path << "\n";
		return false;
	}

	out_map.clear();
	out_map.reserve(1 << 18);

	std::string line;
	long long n_loaded = 0;

	while (std::getline(in, line)) {
		if (line.empty() || line[0] == '#')
			continue;

		std::istringstream ss(line);
		std::string chrom_str, pos_str, p1_str, p2_str;
		if (!(ss >> chrom_str >> pos_str >> p1_str >> p2_str))
			continue;

		// Header detection: pos field must be a valid integer.
		int p;
		try {
			p = std::stoi(pos_str);
		} catch (...) {
			continue;
		}

		float p1f, p2f;
		try {
			p1f = std::stof(p1_str);
			p2f = std::stof(p2_str);
		} catch (...) {
			continue;
		}

		std::string key = chrom_str + ":" + std::to_string(p);
		out_map[key] = MarkerFreq{p1f, p2f};
		++n_loaded;
	}

	std::cout << "Loaded --ref-freq: " << path << " (" << n_loaded << " entries)\n";
	return true;
}

void match_ref_freq(
	const std::unordered_map<std::string, MarkerFreq>& map,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	std::vector<MarkerFreq>& out
) {
	const int nwin = (int)chroms.size();
	out.resize((size_t)nwin);

	for (int w = 0; w < nwin; ++w) {
		std::string key = chroms[(size_t)w] + ":" + std::to_string(pos[(size_t)w]);
		auto it = map.find(key);
		out[(size_t)w] = (it != map.end()) ? it->second : MarkerFreq{};
	}
}

int apply_ref_freq_filter(
	Eigen::MatrixXf& X,
	std::vector<std::string>& chroms,
	std::vector<int>& pos,
	std::vector<int>& pos_start,
	std::vector<MarkerFreq>& freqs,
	float min_delta_afd
) {
	const int nrows = (int)X.rows();
	const int nwin = (int)chroms.size();

	std::vector<int> keep;
	keep.reserve((size_t)nwin);
	int dropped_unmatched = 0;
	int dropped_afd = 0;

	for (int w = 0; w < nwin; ++w) {
		const MarkerFreq& mf = freqs[(size_t)w];
		if (std::isnan(mf.p1) || std::isnan(mf.p2)) {
			++dropped_unmatched;
			continue;
		}
		if (min_delta_afd > 0.0f && std::fabs(mf.p1 - mf.p2) < min_delta_afd) {
			++dropped_afd;
			continue;
		}
		keep.push_back(w);
	}

	const int nkeep = (int)keep.size();

	if (nkeep < nwin) {
		std::cout << "Ref-freq filter:\n";
		std::cout << "  before              = " << nwin << "\n";
		std::cout << "  dropped (unmatched) = " << dropped_unmatched << "\n";
		if (min_delta_afd > 0.0f)
			std::cout << "  dropped (|afd|<" << min_delta_afd << ")  = " << dropped_afd << "\n";
		std::cout << "  after               = " << nkeep << "\n";
	}

	if (nkeep == nwin)
		return nkeep;

	Eigen::MatrixXf Xf(nrows, nkeep);
	std::vector<std::string> chroms_f((size_t)nkeep);
	std::vector<int> pos_f((size_t)nkeep);
	std::vector<MarkerFreq> freqs_f((size_t)nkeep);

	for (int j = 0; j < nkeep; ++j) {
		int w = keep[(size_t)j];
		Xf.col(j) = X.col(w);
		chroms_f[(size_t)j] = chroms[(size_t)w];
		pos_f[(size_t)j] = pos[(size_t)w];
		freqs_f[(size_t)j] = freqs[(size_t)w];
	}

	X.swap(Xf);
	chroms.swap(chroms_f);
	pos.swap(pos_f);
	freqs.swap(freqs_f);

	if (!pos_start.empty()) {
		std::vector<int> ps_f((size_t)nkeep);
		for (int j = 0; j < nkeep; ++j)
			ps_f[(size_t)j] = pos_start[(size_t)keep[(size_t)j]];
		pos_start.swap(ps_f);
	}

	return nkeep;
}

int polarize_X(
	Eigen::MatrixXf& X,
	std::vector<MarkerFreq>& freqs,
	bool phased
) {
	const int nwin = (int)X.cols();
	const float max_dos = phased ? 1.0f : 2.0f;
	int n_flipped = 0;

	for (int w = 0; w < nwin; ++w) {
		MarkerFreq& mf = freqs[(size_t)w];
		if (std::isnan(mf.p1) || std::isnan(mf.p2))
			continue;
		if (mf.p2 > mf.p1) {
			// Flip: dosage=max_dos now encodes pop-1 ancestry
			X.col(w) = (max_dos - X.col(w).array()).matrix();
			std::swap(mf.p1, mf.p2);
			++n_flipped;
		}
	}

	return n_flipped;
}
