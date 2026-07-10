#pragma once

#include <Eigen/Dense>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdint>

#include "scan_markers.hpp"
#include "hybrid_index.hpp"
#include "../io/ref_freq.hpp"

bool scan_markers_write_hits_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1,
	const std::string& reservoir_path = ""
);

bool scan_markers_write_hits_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::vector<MarkerFreq>& freqs_scan,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1,
	const std::string& reservoir_path = ""
);

// LOCO version of scan_markers_write_hits_fdr: per chromosome-pair block,
// HI/Z1/Z2 are recomputed excluding both chromosomes (once, reused across
// the calibration and hit-calling passes). Interchromosomal only.
bool scan_markers_write_hits_excl_focus_fdr(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	uint64_t seed = 1,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1,
	const std::string& reservoir_path = ""
);

bool scan_markers_write_hits_excl_focus_fdr(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::vector<MarkerFreq>& freqs_scan,
	uint64_t seed = 1,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1,
	const std::string& reservoir_path = ""
);

bool scan_target_write_hits_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	int target_w,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1,
	const std::string& reservoir_path = ""
);

bool scan_target_write_hits_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const ScanOptions& opt,
	int target_w,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::vector<MarkerFreq>& freqs_scan,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1,
	const std::string& reservoir_path = ""
);

bool scan_vector_vs_windows_write_hits_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const Eigen::VectorXf& v_raw,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1,
	const std::string& reservoir_path = ""
);

bool scan_vector_vs_windows_write_hits_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const Eigen::VectorXf& v_raw,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::vector<MarkerFreq>& freqs_scan,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1,
	const std::string& reservoir_path = ""
);

// LOCO version of scan_target_write_hits_fdr: HI/residualized vectors are
// recomputed excluding both the target's chromosome and the block's
// chromosome. Interchromosomal only.
bool scan_target_write_hits_excl_focus_fdr(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	int target_w,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	uint64_t seed = 1,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1,
	const std::string& reservoir_path = ""
);

bool scan_target_write_hits_excl_focus_fdr(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const ScanOptions& opt,
	int target_w,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::vector<MarkerFreq>& freqs_scan,
	uint64_t seed = 1,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1,
	const std::string& reservoir_path = ""
);

// LOCO version of scan_vector_vs_windows_write_hits_fdr: HI/residualized
// vectors are recomputed excluding the block's chromosome. Interchromosomal only.
bool scan_vector_vs_windows_write_hits_excl_focus_fdr(
	const Eigen::MatrixXf& X_scan,
	const Eigen::VectorXf& v_raw,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	uint64_t seed = 1,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1,
	const std::string& reservoir_path = ""
);

bool scan_vector_vs_windows_write_hits_excl_focus_fdr(
	const Eigen::MatrixXf& X_scan,
	const Eigen::VectorXf& v_raw,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsFreq& hc_full,
	const ScanOptions& opt,
	const std::string& out_path,
	const std::string& summary_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::vector<MarkerFreq>& freqs_scan,
	uint64_t seed = 1,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1,
	const std::string& reservoir_path = ""
);

