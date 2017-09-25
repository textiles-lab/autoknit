#pragma once

#include "Shape.hpp"

#include <vector>

struct ScheduleCost {
	uint32_t shape = 0; //for awkwardly-oriented cycles
	uint32_t roll = 0; //for stitches that need to change beds
	uint32_t shift = 0; //for stitches that need to move left/right

	ScheduleCost() = default;
	ScheduleCost(uint32_t shape_, uint32_t roll_, uint32_t shift_) : shape(shape_), roll(roll_), shift(shift_) { }

	bool operator<(ScheduleCost const &o) const {
		if (shape != o.shape) return shape < o.shape;
		else if (roll != o.roll) return roll < o.roll;
		else return shift < o.shift;
	}

	bool operator==(ScheduleCost const &o) const {
		return shape == o.shape && roll == o.roll && shift == o.shift;
	}

	ScheduleCost &operator+=(ScheduleCost const &o) {
		//NOTE: unsafe around ScheduleCost::max()
		shape += o.shape;
		roll += o.roll;
		shift += o.shift;
		return *this;
	}

	static ScheduleCost shape_cost( Shape const &shape );

	static ScheduleCost transfer_cost(
		uint32_t from_count, Shape from_shape,
		uint32_t to_count, Shape to_shape,
		std::vector< uint32_t > const &map
		//TODO: bridges?
	);

	static ScheduleCost zero() {
		ScheduleCost z;
		return z;
	};
	static ScheduleCost max() {
		ScheduleCost m(-1U, -1U, -1U);
		return m;
	}
};
