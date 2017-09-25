#pragma once

#include <limits>
#include <vector>
#include <string>

struct BedNeedle {
	enum Bed : char {
		Front = 'f',
		FrontSliders = 'F',
		BackSliders = 'B',
		Back = 'b'
	} bed = Front;
	int32_t needle = 0;
	bool operator==(BedNeedle const &o) const {
		return bed == o.bed && needle == o.needle;
	}
};

struct Constraints {
	int32_t min_free = std::numeric_limits< int32_t >::min();
	int32_t max_free = std::numeric_limits< int32_t >::max();
	uint32_t max_racking = 4;
};

struct Transfer {
	BedNeedle from;
	BedNeedle to;
};

bool plan_transfers(
	Constraints const &constraints,
	std::vector< BedNeedle > const &from_ccw,
	std::vector< BedNeedle > const &to_ccw,
	std::vector< uint8_t > const &slack,
	std::vector< Transfer> *transfers,
	std::string *error = nullptr
);
