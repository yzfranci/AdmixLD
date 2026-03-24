#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

struct WindowMeta {
	std::vector<std::string> chrom;
	std::vector<int> pos;        // end position of block (or sole position for VCF)
	std::vector<int> pos_start;  // start position of block (MSP input only; empty for VCF)
};

struct WindowMatrix {
	Eigen::MatrixXf X;	// nsamples × nwin (raw DS/GT dosage)
	WindowMeta meta;
	std::vector<std::string> sample_names;
};
