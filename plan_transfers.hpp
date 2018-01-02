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
	} bed;
	int32_t needle;
	bool operator==(BedNeedle const &o) const {
		return bed == o.bed && needle == o.needle;
	}
	BedNeedle(Bed bed_ = Front, int32_t needle_ = 0) : bed(bed_), needle(needle_) { }
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

typedef int32_t Slack;
constexpr const Slack SlackForNoYarn = std::numeric_limits< Slack >::max();

bool plan_transfers(
	Constraints const &constraints,
	std::vector< BedNeedle > const &from_ccw,
	std::vector< BedNeedle > const &to_ccw,
	std::vector< Slack > const &slack, //slack between stitch and its ccw-ward neighbor
	std::vector< Transfer> *transfers,
	std::string *error = nullptr
);
