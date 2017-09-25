#include "plan_transfers-helpers.hpp"

#include <iostream>

void best_shift(
	Constraints const &constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const &top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const &bottom,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal > *to_top_,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal > *to_bottom_,
	std::vector< Transfer > *plan_
) {
	//Shift will change both beds:
	assert(top_bed != to_top_bed);
	assert(bottom_bed != to_bottom_bed);

	//make sure all output arrays exist:
	assert(to_top_);
	auto &to_top = *to_top_;
	assert(to_bottom_);
	auto &to_bottom = *to_bottom_;
	assert(plan_);
	auto &plan = *plan_;

	//Clear output:
	to_top.clear();
	to_bottom.clear();
	//NOTE: don't clear plan, only append.

	//------------
	std::cerr << "WARNING: shift not implemented." << std::endl;
	to_top = top;
	to_bottom = bottom;
	(void)plan;
	//------------

}

