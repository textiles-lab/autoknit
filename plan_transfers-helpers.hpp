#pragma once

#include "plan_transfers.hpp"

#include <cassert>

struct NeedleRollGoal {
	NeedleRollGoal() = default;
	NeedleRollGoal(int32_t needle_, int32_t roll_, int32_t goal_) : needle(needle_), roll(roll_), goal(goal_) {
	}
	int32_t needle = 0;
	int32_t roll = 0;
	int32_t goal = 0;

	//slack available on both sides of the stitch:
	Slack left_slack = 0;
	Slack right_slack = 0;

	bool has_same_goal_as(NeedleRollGoal const &o) const {
		return needle == o.needle && roll == o.roll && goal == o.goal;
	}

	uint32_t penalty(int32_t min_free, int32_t max_free) const {
		uint32_t dir = penalty_dir(min_free, max_free);
		uint32_t rec = penalty_rec(min_free, max_free);
		assert(dir == rec);
		return dir;
	}
	//These should be the same, but for now both are called and checked to shake out any bugs:
	uint32_t penalty_dir(int32_t min_free, int32_t max_free) const {
		if (roll == 0) {
			return std::abs(goal - needle);
		} else if (roll < 0) {
			return (needle - min_free)
				+ 2 * -roll
				+ (roll < -1 ? (-roll - 1) * (max_free - min_free) : 0)
				+ (roll % 2 == 0 ? max_free - needle : needle - min_free);
		} else { //(roll > 0)
			return (max_free - needle)
				+ 2 * roll
				+ (roll > 1 ? (roll - 1) * (max_free - min_free) : 0)
				+ (roll % 2 == 0 ? needle - min_free : max_free - needle);
		}
	}
	uint32_t penalty_rec(int32_t min_free, int32_t max_free) const {
		if (roll == 0) {
			return std::abs(goal - needle);
		} else if (roll < 0) {
			NeedleRollGoal after = *this;
			after.needle = min_free;
			after.roll = -(roll + 1);
			return (needle - min_free) + 2 + after.penalty_rec(min_free, max_free);
		} else { //(roll > 0)
			NeedleRollGoal after = *this;
			after.needle = max_free;
			after.roll = -(roll - 1);
			return (max_free - needle) + 2 + after.penalty_rec(min_free, max_free);
		}
	}
};

//collapse 'top' onto 'bottom':
void best_collapse(
	Constraints const &constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const &top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const &bottom,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal > *to_top,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal > *to_bottom,
	std::vector< Transfer > *plan
);

//shift top/bottom and move to other bed:
void best_shift(
	Constraints const &constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const &top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const &bottom,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal > *to_top,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal > *to_bottom,
	std::vector< Transfer > *plan
);

//expand bottom to other bed:
void best_expand(
	Constraints const &constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const &top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const &bottom,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal > *to_top,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal > *to_bottom,
	std::vector< Transfer > *plan
);
