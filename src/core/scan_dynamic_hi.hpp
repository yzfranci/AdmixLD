#pragma once

#include <Eigen/Dense>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdint>

#include "scan_blocks.hpp"
#include "hybrid_index.hpp"

bool scan_blocks_write_hits_excl_focus(
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
	uint64_t distrib_seed = 1
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
	uint64_t distrib_seed = 1
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
	uint64_t distrib_seed = 1
);

bool permute_interchrom_summary_chrblock_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	uint64_t seed,
	int n_perm,
	int sample_size,
	int hi_bins,
	std::vector<PermSummary>& summaries_out
);

bool permute_sample_vector_summary_excl_focus(
	const Eigen::MatrixXf& X_scan,
	const Eigen::VectorXf& g_raw,
	const std::vector<std::string>& chroms_scan,
	const std::vector<int>& pos_scan,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr_scan,
	const std::vector<std::string>& chr_order_scan,
	const HiComponentsWeighted& hc_full,
	const ScanOptions& opt,
	uint64_t seed,
	int n_perm,
	int sample_size,
	int hi_bins,
	std::vector<PermSummary>& summaries_out
);
