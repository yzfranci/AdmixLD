#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

// Load per-sample haplotype vector from TSV (sample<TAB>value).
// Non-NaN values must be exactly 0 or 1; any other value is an error.
// Missing values: ".", "NA", "NaN", "nan", "NAN" are treated as NaN.
bool load_sample_vector_tsv(
	const std::string& path,
	const std::vector<std::string>& sample_names,
	Eigen::VectorXf& y_out
);
