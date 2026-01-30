#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

struct WindowMeta {
	std::vector<std::string> chrom;
	std::vector<int> pos;
};

struct WindowMatrix {
	Eigen::MatrixXf X;	// nsamples × nwin (raw DS/GT dosage)
	WindowMeta meta;
	std::vector<std::string> sample_names;
};
