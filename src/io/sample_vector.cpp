#include "sample_vector.hpp"

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <iostream>
#include <limits>
#include <cmath>

static bool is_missing_token(const std::string& s) {
	if (s.empty())
		return true;
	if (s == "." || s == "NA" || s == "NaN" || s == "nan" || s == "NAN")
		return true;
	return false;
}

bool load_sample_vector_tsv(
	const std::string& path,
	const std::vector<std::string>& sample_names,
	Eigen::VectorXf& y_out
) {
	std::ifstream in(path);
	if (!in) {
		std::cerr << "Error: cannot open --sample-haplo file: " << path << "\n";
		return false;
	}

	std::unordered_map<std::string, float> map;
	map.reserve(sample_names.size() * 2 + 1);

	std::string line;
	while (std::getline(in, line)) {
		if (line.empty() || line[0] == '#')
			continue;

		std::istringstream ss(line);
		std::string sname;
		std::string sval;

		if (!(ss >> sname >> sval))
			continue;

		if (is_missing_token(sval)) {
			map[sname] = std::numeric_limits<float>::quiet_NaN();
			continue;
		}

		try {
			float v = std::stof(sval);
			map[sname] = v;
		} catch (...) {
			std::cerr << "Error: non-numeric value for sample '" << sname
				<< "': '" << sval << "'.\n";
			return false;
		}
	}

	int n = (int)sample_names.size();
	y_out.resize(n);

	int found = 0;
	int missing = 0;

	for (int i = 0; i < n; ++i) {
		auto it = map.find(sample_names[i]);
		if (it == map.end()) {
			y_out(i) = std::numeric_limits<float>::quiet_NaN();
			++missing;
		} else {
			y_out(i) = it->second;
			if (std::isnan(y_out(i)))
				++missing;
			else
				++found;
		}
	}

	std::cout << "Loaded sample-haplo vector:\n";
	std::cout << "  file    = " << path << "\n";
	std::cout << "  found   = " << found << "\n";
	std::cout << "  missing = " << missing << "\n";

	if (found < 3) {
		std::cerr << "Error: too few non-missing sample-haplo values (" << found << ").\n";
		return false;
	}

	// Strict 0/1 validation
	for (int i = 0; i < n; ++i) {
		float v = y_out(i);
		if (!std::isnan(v) && v != 0.0f && v != 1.0f) {
			std::cerr << "Error: --sample-haplo requires strictly 0 or 1 values (use missing tokens for absent data).\n";
			std::cerr << "  Sample '" << sample_names[i] << "' has value: " << v << "\n";
			return false;
		}
	}

	return true;
}
