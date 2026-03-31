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
	Eigen::MatrixXf X;	// nrows × nwin; nrows = nsamples (diploid) or 2*nsamples (phased)
	WindowMeta meta;
	std::vector<std::string> sample_names;	// always diploid sample names (length = nsamples_diploid)
	bool phased = false;
	int nsamples_diploid = 0;	// number of diploid samples; X.rows() == (phased ? 2 : 1) * nsamples_diploid
};
