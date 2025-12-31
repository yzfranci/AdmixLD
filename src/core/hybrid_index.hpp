#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

// Compute hybrid index from X (dosage 0/1/2) as mean(dosage)/2 per sample, ignoring NaNs
Eigen::VectorXf compute_hi_from_X(
	const Eigen::MatrixXf& X
);

// Load hybrid index from TSV: sample<TAB>hi (header allowed). Requires all samples.
bool load_hi_tsv(
	const std::string& path,
	const std::vector<std::string>& sample_names,
	Eigen::VectorXf& h_out
);

// Write hybrid index to TSV: sample<TAB>hi
bool write_hi_tsv(
	const std::string& path,
	const std::vector<std::string>& sample_names,
	const Eigen::VectorXf& h
);
