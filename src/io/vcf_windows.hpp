#pragma once

#include <string>
#include "../core/types.hpp"

struct VcfLoadOptions {
	int max_windows = 50000;	// 0/negative = “all” (we cap internally)
};

WindowMatrix load_windows_from_vcf(
	const std::string& vcf_path,
	const VcfLoadOptions& opt
);
