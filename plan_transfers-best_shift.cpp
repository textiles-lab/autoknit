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

	int32_t best_ofs = 0;
	uint32_t best_penalty = std::numeric_limits< uint32_t >::max();

	//Find best offset to shift to:

	auto do_ofs = [&](int32_t ofs) {
		//would shifting to this offset move outside limits?
		if (!top.empty()) {
			if (top[0].needle + ofs < constraints.min_free) return;
			if (top.back().needle + ofs > constraints.max_free) return;
		}
		if (!bottom.empty()) {
			if (bottom[0].needle + ofs < constraints.min_free) return;
			if (bottom.back().needle + ofs > constraints.max_free) return;
		}
		uint32_t penalty = 0;
		for (auto const &nrg : top) {
			penalty += nrg.after_offset_and_roll(ofs, 0).penalty(constraints.min_free, constraints.max_free);
		}
		for (auto const &nrg : bottom) {
			penalty += nrg.after_offset_and_roll(ofs, 0).penalty(constraints.min_free, constraints.max_free);
		}
		if (penalty < best_penalty) {
			best_penalty = penalty;
			best_ofs = ofs;
		}
		//std::cout << "  offset " << ofs << " has penalty " << penalty << std::endl; //DEBUG
	};

	for (int32_t ofs = 0; ofs < int32_t(constraints.max_racking); ++ofs) {
		do_ofs(ofs);
		if (ofs != 0) do_ofs(-ofs);
	}
	std::cout << "Best offset: " << best_ofs << std::endl; //DEBUG

	assert(best_penalty < std::numeric_limits< uint32_t >::max());

	//shift by that offset:

	std::vector< Transfer > ops;
	for (auto const &nrg : top) {
		ops.emplace_back(BedNeedle(top_bed, nrg.needle), BedNeedle(to_top_bed, nrg.needle + best_ofs));
	}
	for (auto const &nrg : bottom) {
		ops.emplace_back(BedNeedle(bottom_bed, nrg.needle), BedNeedle(to_bottom_bed, nrg.needle + best_ofs));
	}

	run_transfers(constraints,
		top_bed, top,
		bottom_bed, bottom,
		ops,
		to_top_bed, &to_top,
		to_bottom_bed, &to_bottom
	);

	std::cout << "Before Shift:\n"; //DEBUG
	draw_beds(top_bed, top, bottom_bed, bottom); //DEBUG

	std::cout << "After Shift:\n"; //DEBUG
	draw_beds(to_top_bed, to_top, to_bottom_bed, to_bottom); //DEBUG



	plan.insert(plan.end(), ops.begin(), ops.end());


}

