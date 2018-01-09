#pragma once

#include <limits>
#include <vector>
#include <string>
#include <cassert>

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
	std::string to_string() const {
		if (bed == Front)             return "f"  + std::to_string(needle);
		else if (bed == FrontSliders) return "fs" + std::to_string(needle);
		else if (bed == BackSliders)  return "bs" + std::to_string(needle);
		else if (bed == Back)         return "b"  + std::to_string(needle);
		else assert(0 && "Should have handled all cases");
	}
};

template< >
struct std::hash< BedNeedle > {
	size_t operator()(BedNeedle const &bn) const {
		return (size_t(bn.bed) << (8*(sizeof(size_t)-sizeof(bn.bed)))) ^ size_t(bn.needle);
	};
};

struct Constraints {
	int32_t min_free = std::numeric_limits< int32_t >::min();
	int32_t max_free = std::numeric_limits< int32_t >::max();
	uint32_t max_racking = 4;
};

struct Transfer {
	Transfer() = default;
	Transfer(BedNeedle const &_from, BedNeedle const &_to) : from(_from), to(_to) { }
	BedNeedle from;
	BedNeedle to;
	std::string to_string() const {
		return from.to_string() + " -> " + to.to_string();
	}
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
