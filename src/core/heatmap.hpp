#pragma once

#include <Eigen/Dense>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdint>

#include "hybrid_index.hpp"
#include "../io/ref_freq.hpp"

// --heatmap: bins marker pairs into a 2D grid of (chrA, binA) x (chrB, binB)
// cells by genomic position, reporting a running mean/SD (exact) and an
// approximate upper quantile (bounded reservoir per cell) instead of writing
// every tested pair. A standalone scan mode: no hit-calling, no --fdr, no
// --target-*/--sample-haplo.
struct HeatmapOptions {
	bool intra = false;
	int max_dist = -1;
	int min_dist = -1;
	int tile_size = 1024;
	int nsamples = 0;
	int threads = 1;

	long long bin_size = 1000000;
	double quantile = 0.99;
	int reservoir_size = 2000;	// per-cell reservoir cap for the approximate quantile
	uint64_t seed = 1;
};

// Global-HI heatmap scan.
bool scan_markers_write_heatmap(
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const HeatmapOptions& opt,
	const std::string& out_path,
	long long& tested_pairs
);

// LOCO (--hi-mode excl-focus) heatmap scan: each block's HI/residualized
// data is recomputed excluding the relevant chromosome(s), same as
// scan_markers_write_hits_excl_focus.
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
);

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
);
