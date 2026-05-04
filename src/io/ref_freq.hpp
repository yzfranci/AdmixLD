#pragma once

#include <Eigen/Dense>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

struct MarkerFreq {
	float p1 = std::numeric_limits<float>::quiet_NaN();
	float p2 = std::numeric_limits<float>::quiet_NaN();
};

// Load parental allele frequencies from TSV into a key→MarkerFreq map.
// Key format: "chrom:pos". File format: chrom<TAB>pos<TAB>p1<TAB>p2.
// Header line auto-detected by checking whether pos field parses as integer.
bool load_ref_freq_map(
	const std::string& path,
	std::unordered_map<std::string, MarkerFreq>& out_map
);

// Match per-marker freqs from a pre-loaded map for a given marker list.
// Unmatched markers get NaN p1/p2. out has the same length as chroms/pos.
void match_ref_freq(
	const std::unordered_map<std::string, MarkerFreq>& map,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	std::vector<MarkerFreq>& out
);

// Drop markers with NaN freqs (unmatched) or |p1-p2| < min_delta_afd.
// Rebuilds X, chroms, pos, pos_start, freqs in-place.
// Returns number of markers remaining; prints a summary if any are dropped.
int apply_ref_freq_filter(
	Eigen::MatrixXf& X,
	std::vector<std::string>& chroms,
	std::vector<int>& pos,
	std::vector<int>& pos_start,
	std::vector<MarkerFreq>& freqs,
	float min_delta_afd
);

// Polarize X so dosage=2 (or 1 phased) always means pop-1 ancestry.
// For markers where freqs[w].p2 > freqs[w].p1: flips X.col(w) and swaps p1/p2.
// Returns count of flipped markers.
int polarize_X(
	Eigen::MatrixXf& X,
	std::vector<MarkerFreq>& freqs,
	bool phased
);
