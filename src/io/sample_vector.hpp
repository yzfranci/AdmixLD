#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

bool load_sample_vector_tsv(
	const std::string& path,
	const std::vector<std::string>& sample_names,
	Eigen::VectorXf& y_out
);
