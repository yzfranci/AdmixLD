#pragma once

#include <string>
#include <unordered_map>
#include <vector>

struct BedInterval {
	int start = 0;
	int end = 0;
};

struct BedRow {
	std::string chr;
	int start = 0;
	int end = 0;
};

bool read_bed(
	const std::string& path,
	std::unordered_map<std::string, std::vector<BedInterval>>& by_chr
);

bool bed_contains(
	const std::unordered_map<std::string, std::vector<BedInterval>>& by_chr,
	const std::string& chr,
	int pos
);

// Rows kept as given (no sorting/merging), for per-window scans.
bool read_bed_rows(
	const std::string& path,
	std::vector<BedRow>& rows
);
