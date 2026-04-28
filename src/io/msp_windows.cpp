#include "msp_windows.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

WindowMatrix load_windows_from_msp(
	const std::string& msp_path,
	const MspLoadOptions& opt
) {
	WindowMatrix out;

	std::ifstream in(msp_path);
	if (!in)
		throw std::runtime_error("Failed to open MSP file: " + msp_path);

	std::string line;

	// Line 1: "#Subpopulation order/codes: ..."
	if (!std::getline(in, line) || line.empty() || line[0] != '#')
		throw std::runtime_error("MSP file missing subpopulation header line: " + msp_path);

	// Line 2: "#chm\tspos\tepos\tsgpos\tegpos\tn snps\t[sample.0\tsample.1\t...]"
	if (!std::getline(in, line) || line.empty() || line[0] != '#')
		throw std::runtime_error("MSP file missing column header line: " + msp_path);

	// Strip leading '#' and split by tab
	std::vector<std::string> header_fields;
	{
		std::istringstream ss(line.substr(1));
		std::string tok;
		while (std::getline(ss, tok, '\t'))
			header_fields.push_back(tok);
	}

	// Fixed columns: chm(0) spos(1) epos(2) sgpos(3) egpos(4) "n snps"(5)
	const int N_FIXED = 6;
	if ((int)header_fields.size() < N_FIXED + 2)
		throw std::runtime_error("MSP file header has too few columns (expected at least 2 haplotype columns): " + msp_path);

	int n_hap_cols = (int)header_fields.size() - N_FIXED;
	if (n_hap_cols % 2 != 0)
		throw std::runtime_error("MSP file has an odd number of haplotype columns (expected pairs .0/.1): " + msp_path);

	int nsamples = n_hap_cols / 2;

	// Extract sample names from haplotype column headers (strip ".0" suffix from the first of each pair)
	out.sample_names.reserve((size_t)nsamples);
	for (int i = 0; i < nsamples; ++i) {
		std::string hap0 = header_fields[(size_t)(N_FIXED + 2 * i)];
		if (hap0.size() >= 2 && hap0.substr(hap0.size() - 2) == ".0")
			hap0 = hap0.substr(0, hap0.size() - 2);
		out.sample_names.push_back(hap0);
	}

	out.nsamples_diploid = nsamples;
	out.phased = opt.phased;

	const int max_windows_load = opt.max_windows;  // 0 = no limit

	const int nrows = opt.phased ? 2 * nsamples : nsamples;
	int capacity = (max_windows_load > 0) ? max_windows_load : 65536;
	out.X = Eigen::MatrixXf(nrows, capacity);
	const float NA = std::numeric_limits<float>::quiet_NaN();

	if (max_windows_load > 0) {
		out.meta.chrom.reserve((size_t)max_windows_load);
		out.meta.pos.reserve((size_t)max_windows_load);
		out.meta.pos_start.reserve((size_t)max_windows_load);
	}

	int nwin = 0;

	while (std::getline(in, line) && (max_windows_load <= 0 || nwin < max_windows_load)) {
		if (line.empty() || line[0] == '#')
			continue;

		std::vector<std::string> fields;
		{
			std::istringstream ss(line);
			std::string tok;
			while (std::getline(ss, tok, '\t'))
				fields.push_back(tok);
		}

		if ((int)fields.size() < N_FIXED + n_hap_cols)
			continue;	// skip malformed lines

		if (nwin == capacity) {
			capacity *= 2;
			out.X.conservativeResize(Eigen::NoChange, capacity);
		}

		const std::string& chrom = fields[0];
		int spos = std::stoi(fields[1]);
		int epos = std::stoi(fields[2]);
		// fields[3]=sgpos, fields[4]=egpos, fields[5]=n_snps — not used

		out.meta.chrom.emplace_back(chrom);
		out.meta.pos.push_back(epos);        // pos = end of tract (consistent with VCF marker-position field)
		out.meta.pos_start.push_back(spos);  // pos_start = start of tract

		for (int i = 0; i < nsamples; ++i) {
			const std::string& v0 = fields[(size_t)(N_FIXED + 2 * i)];
			const std::string& v1 = fields[(size_t)(N_FIXED + 2 * i + 1)];

			if ((v0 != "0" && v0 != "1") || (v1 != "0" && v1 != "1")) {
				if (opt.phased) {
					out.X(2 * i,     nwin) = NA;
					out.X(2 * i + 1, nwin) = NA;
				} else {
					out.X(i, nwin) = NA;
				}
			} else if (opt.phased) {
				out.X(2 * i,     nwin) = (float)std::stoi(v0);
				out.X(2 * i + 1, nwin) = (float)std::stoi(v1);
			} else {
				out.X(i, nwin) = (float)(std::stoi(v0) + std::stoi(v1));
			}
		}

		++nwin;
	}

	if (max_windows_load > 0 && nwin == max_windows_load)
		std::cerr << "Warning: reached --max-windows limit (" << max_windows_load << "); remaining MSP records were not read.\n";

	out.X.conservativeResize(Eigen::NoChange, nwin);
	return out;
}
