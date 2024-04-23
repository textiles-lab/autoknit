#include "Stitch.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

bool load_stitches(std::string const &filename, std::vector< Stitch > *into_) {
	assert(into_);
	auto &into = *into_;
	into.clear();

	std::ifstream file(filename);
	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		Stitch temp;

		int32_t in[2];
		int32_t out[2];

		if (!(iss >> temp.yarn >> temp.type >> temp.direction >> in[0] >> in[1] >> out[0] >> out[1] >> temp.at.x >> temp.at.y >> temp.at.z)) {
			std::cerr << "ERROR: Failed to read stitch from '" << line << "'" << std::endl;
			return false;
		}
		temp.in[0] = in[0];
		temp.in[1] = in[1];
		temp.out[0] = out[0];
		temp.out[1] = out[1];

		if (!temp.check_type()) {
			std::cerr << "ERROR: Stitch does not have proper in/out for type." << std::endl;
			std::cerr << "  line: '" << line << "'" << std::endl;
			return false;
		}


		into.emplace_back(temp);
	}
	for (auto const &s : into) {
		uint32_t idx = &s - &into[0];
		auto check_in = [&](uint32_t in_idx) -> bool {
			if (in_idx == -1U) {
				return true;
			} else {
				if (in_idx >= idx) return false;
				if (into[in_idx].find_out(idx) == -1U) return false;
				return true;
			}
		};
		if (!check_in(s.in[0]) || !check_in(s.in[1])) {
			std::cerr << "Stitch does not have proper 'in' array." << std::endl;
			return false;
		}
		auto check_out = [&](uint32_t out_idx) -> bool {
			if (out_idx == -1U) {
				return true;
			} else {
				if (out_idx >= into.size()) return false;
				if (out_idx <= idx) return false;
				if (into[out_idx].find_in(idx) == -1U) return false;
				return true;
			}
		};
		if (!check_out(s.out[0]) || !check_out(s.out[1])) {
			std::cerr << "Stitch does not have proper 'out' array." << std::endl;
			return false;
		}
	}
	return true;
}

void save_stitches(std::string const &filename, std::vector< Stitch > const &from) {
	std::ofstream file(filename);
	for (auto const &s : from) {
		file << s.yarn
		<< ' ' << s.type
		<< ' ' << s.direction
		<< ' ' << (int32_t)s.in[0]
		<< ' ' << (int32_t)s.in[1]
		<< ' ' << (int32_t)s.out[0]
		<< ' ' << (int32_t)s.out[1]
		<< ' ' << s.at.x << ' ' << s.at.y << ' ' << s.at.z << '\n';
	}
}
