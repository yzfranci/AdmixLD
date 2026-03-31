#pragma once

#include <string>
#include "../core/types.hpp"

struct MspLoadOptions {
	int max_windows = 50000;	// 0/negative = "all" (capped internally)
	bool phased = false;		// store haplotypes as separate rows (2n rows) instead of summing to dosage
};

WindowMatrix load_windows_from_msp(
	const std::string& msp_path,
	const MspLoadOptions& opt
);
