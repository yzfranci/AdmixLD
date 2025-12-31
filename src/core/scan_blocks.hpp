#pragma once

#include <Eigen/Dense>
#include <string>
#include <unordered_map>
#include <vector>

struct ScanOptions {
	bool intra = false;
	int block_size = 1024;
	float min_abs_r = 0.2f;
	int nsamples = 0;	// used for output column 'n'
};

std::unordered_map<std::string, std::vector<int>> group_by_chr(
	const std::vector<std::string>& chroms,
	std::vector<std::string>& chr_order
);

bool scan_blocks_write_hits(
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& starts,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs
);

struct PermSummary {
	float min_r = 0.0f;
	float p01 = 0.0f;
	float p05 = 0.0f;
	float median = 0.0f;
	float p95 = 0.0f;
	float p99 = 0.0f;
	float max_r = 0.0f;
};

bool permute_interchrom_summary_chrblock(
	const Eigen::MatrixXf& Z,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	uint64_t seed,
	int n_perm,
	int sample_size,
	std::vector<PermSummary>& summaries_out
);