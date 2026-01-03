#include "bed.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

bool read_bed(
	const std::string& path,
	std::unordered_map<std::string, std::vector<BedInterval>>& by_chr
) {
	std::ifstream in(path);
	if (!in) {
		std::cerr << "Error: cannot open BED file: " << path << "\n";
		return false;
	}

	by_chr.clear();

	std::string line;
	while (std::getline(in, line)) {
		if (line.empty() || line[0] == '#')
			continue;

		std::istringstream ss(line);
		std::string chr;
		int start = 0;
		int end = 0;

		if (!(ss >> chr >> start >> end))
			continue;

		if (end <= start)
			continue;

		by_chr[chr].push_back(BedInterval{start, end});
	}

	// Sort and merge intervals per chromosome
	for (auto& kv : by_chr) {
		auto& v = kv.second;
		if (v.empty())
			continue;

		std::sort(v.begin(), v.end(),
			[](const BedInterval& a, const BedInterval& b) {
				if (a.start != b.start) return a.start < b.start;
				return a.end < b.end;
			}
		);

		std::vector<BedInterval> merged;
		merged.reserve(v.size());

		for (const auto& it : v) {
			if (merged.empty() || it.start > merged.back().end) {
				merged.push_back(it);
			} else if (it.end > merged.back().end) {
				merged.back().end = it.end;
			}
		}

		v.swap(merged);
	}

	return true;
}

bool bed_contains(
	const std::unordered_map<std::string, std::vector<BedInterval>>& by_chr,
	const std::string& chr,
	int pos
) {
	auto it = by_chr.find(chr);
	if (it == by_chr.end())
		return false;

	const auto& v = it->second;
	int lo = 0;
	int hi = (int)v.size();

	// BED intervals are half-open: [start, end)
	while (lo < hi) {
		int mid = lo + (hi - lo) / 2;
		if (pos < v[mid].start) {
			hi = mid;
		} else if (pos >= v[mid].end) {
			lo = mid + 1;
		} else {
			return true;
		}
	}
	return false;
}
