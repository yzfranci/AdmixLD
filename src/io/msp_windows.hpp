#pragma once

#include <string>
#include "../core/types.hpp"

struct MspLoadOptions {
	int max_windows = 50000;	// 0/negative = "all" (capped internally)
};

WindowMatrix load_windows_from_msp(
	const std::string& msp_path,
	const MspLoadOptions& opt
);
