#include "hybrid_index.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_map>

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
