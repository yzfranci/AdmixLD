#pragma once

#include <Eigen/Dense>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <random>

struct ScanOptions {
	bool intra = false;
	int max_dist = -1;
	int tile_size = 1024;

	// legacy symmetric filter
	float min_abs_r = 0.0f;

	// asymmetric filters
	bool use_asym = false;
	float min_neg_r = 0.0f;	// threshold magnitude, applied as r <= -min_neg_r
	float min_pos_r = 0.0f;	// threshold applied as r >= min_pos_r

	int nsamples = 0;

	int threads = 1;
};

std::unordered_map<std::string, std::vector<int>> group_by_chr(
	const std::vector<std::string>& chroms,
	std::vector<std::string>& chr_order
);

struct ScanSummary {
	long long tested_pairs = 0;
	float min_r = 0.0f;
	float p01 = 0.0f;
	float p05 = 0.0f;
	float p25 = 0.0f;
	float median = 0.0f;
	float p75 = 0.0f;
	float p95 = 0.0f;
	float p99 = 0.0f;
	float max_r = 0.0f;
	float mean = 0.0f;
	float sd = 0.0f;
};


bool scan_markers_write_hits(
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1
);

struct PermSummary {
	float min_r = 0.0f;
	float p01 = 0.0f;
	float p05 = 0.0f;
	float p25 = 0.0f;
	float median = 0.0f;
	float p75 = 0.0f;
	float p95 = 0.0f;
	float p99 = 0.0f;
	float max_r = 0.0f;
	float mean = 0.0f;
	float sd = 0.0f;
	float mean_r2 = 0.0f;
	float sd_r2 = 0.0f;
};


bool permute_interchrom_summary_chrmarker(
	const Eigen::MatrixXf& Z,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	uint64_t seed,
	int n_perm,
	int sample_size,
	std::vector<PermSummary>& summaries_out
);

bool scan_target_write_hits(
	const Eigen::MatrixXf& Z,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	int target_w,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1
);

bool scan_vector_vs_windows_write_hits(
	const Eigen::MatrixXf& Z,
	const Eigen::VectorXf& v,
	const std::vector<std::string>& chroms,
	const std::vector<int>& pos,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	const std::string& out_path,
	long long& tested_pairs,
	long long& kept_pairs,
	const std::string& distrib_path = "",
	int distrib_sample = 200000,
	uint64_t distrib_seed = 1
);

bool permute_sample_vector_summary(
	const Eigen::MatrixXf& Z,
	const Eigen::VectorXf& gZ,
	const std::unordered_map<std::string, std::vector<int>>& windows_by_chr,
	const std::vector<std::string>& chr_order,
	const ScanOptions& opt,
	uint64_t seed,
	int n_perm,
	int sample_size,
	std::vector<PermSummary>& summaries_out
);
