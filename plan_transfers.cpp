#include "plan_transfers.hpp"
#include "plan_transfers-helpers.hpp"

#include <cassert>
#include <algorithm>
#include <cmath>

bool plan_transfers(
	Constraints const &constraints,
	std::vector< BedNeedle > const &from_,
	std::vector< BedNeedle > const &to_,
	std::vector< Slack > const &slack_,
	std::vector< Transfer> *transfers_,
	std::string *error
) {
	assert(constraints.min_free < constraints.max_free);
	assert(constraints.max_racking >= 1);

	assert(from_.size() == to_.size());
	assert(from_.size() == slack_.size());
	assert(transfers_);

	auto &transfers = *transfers_;
	transfers.clear();
	if (error) *error = "";

	if (from_.empty()) return true; //empty cycle has empty plan

	std::vector< BedNeedle > from = from_;
	std::vector< BedNeedle > to = to_;
	std::vector< Slack > slack = slack_;
	//eliminate any stacked needles in 'from':
	for (uint32_t i = 0; i < from.size(); /* later */) {
		if (from.size() == 1) break;
		uint32_t n = (i + 1 < from.size() ? i + 1 : 0);
		if (from[i] == from[n]) {
			assert(to[i] == to[n]);
			from.erase(from.begin() + i);
			to.erase(to.begin() + i);
			slack.erase(slack.begin() + i);
		} else {
			++i;
		}
	}
	const uint32_t Count = from.size();
	assert(from.size() == Count);
	assert(to.size() == Count);
	assert(slack.size() == Count);


	//PARANOIA: check consistency of inputs:
	auto assert_valid_layout = [&](std::vector< BedNeedle > const &cycle, bool require_zero_ofs) {
		Slack left_slack = SlackForNoYarn;
		int32_t left_offset = 0;
		Slack right_slack = SlackForNoYarn;
		int32_t right_offset = 0;
		for (uint32_t i = 0; i < Count; ++i) {
			//in the constrained region:
			assert(cycle[i].needle >= constraints.min_free);
			assert(cycle[i].needle <= constraints.max_free);
			//slack:
			uint32_t n = (i + 1 < Count ? i + 1 : 0);
			if (cycle[i].bed == cycle[n].bed) {
				assert(std::abs(cycle[n].needle - cycle[i].needle) <= slack[i]);
			} else if (cycle[i].bed == BedNeedle::Front && cycle[n].bed == BedNeedle::Back) {
				assert(left_slack == SlackForNoYarn);
				left_offset = cycle[n].needle - cycle[i].needle;
				left_slack = slack[i];
			} else if (cycle[i].bed == BedNeedle::Back && cycle[n].bed == BedNeedle::Front) {
				assert(right_slack == SlackForNoYarn);
				right_offset = cycle[i].needle - cycle[n].needle;
				right_slack = slack[i];
			} else {
				assert(0);
			}
		}
		//std::cout << left_offset << "/" << int32_t(left_slack) << " , " << right_offset << "/" << int32_t(right_slack) << "  max_racking is " << constraints.max_racking << std::endl; //DEBUG
		//there should be at least one valid racking:
		bool has_valid = false;
		bool has_valid_at_zero = false;
		for (int32_t r = -int32_t(constraints.max_racking); r <= int32_t(constraints.max_racking); ++r) {
			if (std::abs(left_offset + r) > left_slack) continue;
			if (std::abs(right_offset + r) > right_slack) continue;
			has_valid = true;
			if (r == 0) has_valid_at_zero = true;
		}
		if (require_zero_ofs) {
			assert(has_valid_at_zero && "Must have valid layout at zero racking");
		}
		assert(has_valid && "Must have valid racking");
	};
	assert_valid_layout(from, false);
	assert_valid_layout(to, true);

	//can't split loops:
	for (uint32_t i = 0; i < Count; ++i) {
		uint32_t n = (i + 1 < Count ? i + 1 : 0);
		if (from[i] == from[n]) {
			assert(to[i] == to[n]);
		}
	}

	auto assert_ccw = [](std::vector< BedNeedle > const &cycle) {
		bool has_front_back = false;
		bool has_back_front = false;
		for (uint32_t i = 0; i < cycle.size(); ++i) {
			uint32_t n = (i + 1 < cycle.size() ? i + 1 : 0);
			if (cycle[i].bed == BedNeedle::Front) {
				if (cycle[n].bed == BedNeedle::Front) {
					if (cycle[i].needle <= cycle[n].needle) {
						//great!
					} else {
						//single-bed cycle, e.g.:
						// . . . . .
						// . 3 1 2 .
						assert(!has_front_back);
						assert(!has_back_front);
						has_front_back = has_back_front = true;
					}
				} else { assert(cycle[n].bed == BedNeedle::Back);
					assert(!has_front_back);
					has_front_back = true;
				}
			} else { assert(cycle[i].bed == BedNeedle::Back);
				if (cycle[n].bed == BedNeedle::Back) {
					if (cycle[i].needle >= cycle[n].needle) {
						//great!
					} else {
						//single-bed cycle, e.g.:
						// . 2 1 3 .
						// . . . . .
						assert(!has_back_front);
						assert(!has_front_back);
						has_back_front = has_front_back = true;
					}

				} else { assert(cycle[n].bed == BedNeedle::Front);
					assert(!has_back_front);
					has_back_front = true;
				}
			}
		}
		assert(has_front_back == has_back_front); //might be a degenerate "all on one bed" cycle.
	};
	assert_ccw(from);
	assert_ccw(to);

	//end PARANOIA


	//(a) compute roll/goal for each from index.
	//how many bed swaps between p and n:
	auto swaps = [](BedNeedle const &p, BedNeedle const &n) {
		if (p.bed == n.bed) {
			if (p.bed == BedNeedle::Front) {
				return (p.needle <= n.needle ? 0 : 2);
			} else { assert(p.bed == BedNeedle::Back);
				return (p.needle >= n.needle ? 0 : 2);
			}
		} else {
			return 1;
		}
	};
	std::vector< int32_t > winding;
	winding.reserve(Count);
	winding.emplace_back(from[0].bed == to[0].bed ? 0 : 1);
	for (uint32_t i = 1; i < Count; ++i) {
		uint32_t p = i - 1;
		winding.emplace_back(
			winding.back()
			- swaps(from[p], from[i])
			+ swaps(to[p], to[i])
		);
	}
	assert(winding.size() == Count);

	minimize_winding(&winding);

	std::vector< NeedleRollGoal > rg;
	rg.reserve(Count);
	for (uint32_t i = 0; i < Count; ++i) {
		uint32_t p = (i > 0 ? i - 1 : Count - 1);
		uint32_t n = (i + 1 < Count ? i + 1 : 0);
		rg.emplace_back(
			from[i].needle,
			(from[i].bed == BedNeedle::Front ? winding[i] : -winding[i]),
			to[i].needle,
			(from[i].bed == BedNeedle::Front ? slack[p] : slack[i]), //left slack
			(from[i].bed == BedNeedle::Front ? slack[i] : slack[p]) //right slack
		);
		rg.back().can_stack_left = (to[i] == (from[i].bed == BedNeedle::Front ? to[p] : to[n]));
		rg.back().can_stack_right = (to[i] == (from[i].bed == BedNeedle::Front ? to[n] : to[p]));
		//can't stack stitches with themselves:
		if (i == n) {
			assert(i == p);
			rg.back().can_stack_left = rg.back().can_stack_right = false;
		}

		//check to make sure roll/goal actually matches desired behavior:
		assert((from[i].bed == to[i].bed) == (rg.back().roll % 2 == 0));
	}
	assert(rg.size() == Count);

	//PARANOIA: make sure roll/goal matches on any overlapped needles:
	for (uint32_t i = 0; i < Count; ++i) {
		uint32_t n = (i + 1 < Count ? i + 1 : 0);
		assert((from[i] == from[n]) == (rg[i].has_same_goal_as(rg[n])));
	}

	//(b) Get everything into left-to-right sorted beds:
	std::vector< NeedleRollGoal > front, back;
	for (uint32_t i = 0; i < Count; ++i) {
		if (from[i].bed == BedNeedle::Front) {
			front.emplace_back(rg[i]);
		} else { assert(from[i].bed == BedNeedle::Back);
			back.emplace_back(rg[i]);
		}
	}
	auto compare_needles = [](NeedleRollGoal const &a, NeedleRollGoal const &b) {
		return a.needle < b.needle;
	};
	std::sort(front.begin(), front.end(), compare_needles);
	std::sort(back.begin(), back.end(), compare_needles);

	//(c) Potentially shrink free range to avoid problems:
	Constraints shrunk_constraints = constraints;
	int32_t min_used = std::numeric_limits< int32_t >::max();
	int32_t max_used = std::numeric_limits< int32_t >::min();
	for (uint32_t i = 0; i < Count; ++i) {
		min_used = std::min(min_used, std::min(from[i].needle, to[i].needle));
		max_used = std::max(max_used, std::max(from[i].needle, to[i].needle));
	}
	int32_t margin = std::max(int32_t(shrunk_constraints.max_racking), max_used - min_used + 1);
	shrunk_constraints.min_free = std::max(shrunk_constraints.min_free, min_used - margin);
	shrunk_constraints.max_free = std::min(shrunk_constraints.max_free, max_used + margin);

	//(c) Iterate Collapse/Shift/Expand until finished:
	while (1) {
		//--------------------------
		//eliminate duplicate needles:
		auto remove_duplicates = [](std::vector< NeedleRollGoal > &bed) {
			for (uint32_t i = 0; i + 1 < bed.size(); /* later */) {
				if (bed[i].needle == bed[i+1].needle) {
					assert(bed[i].can_stack_left || bed[i+1].can_stack_right);
					assert(bed[i].has_same_goal_as(bed[i+1]));
					bed[i+1].left_slack = bed[i].left_slack; //track slack
					bed[i+1].can_stack_left = bed[i].can_stack_left;
					bed.erase(bed.begin() + i);
				} else {
					assert(bed[i].needle < bed[i+1].needle);
					++i;
				}
			}
		};
		remove_duplicates(front);
		remove_duplicates(back);

		//--------------------------
		//stop if solved:
		bool solved = true;
		for (auto const &nrg : front) {
			if (!(nrg.roll == 0 && nrg.needle == nrg.goal)) {
				solved = false;
				break;
			}
		}
		if (solved) {
			for (auto const &nrg : back) {
				if (!(nrg.roll == 0 && nrg.needle == nrg.goal)) {
					solved = false;
					break;
				}
			}
		}
		if (solved) break;

		//--------------------------
		//TODO, maybe: check for canonical configuration in cache, return solution therefrom.

		//--------------------------
		//try both collapse/stretch/expand options:

		auto penalty = [](Constraints const &cons, std::vector< NeedleRollGoal > const &front, std::vector< NeedleRollGoal > const &back) {
			uint32_t ret = 0;
			for (auto const &nrg : front) {
				ret += nrg.penalty(cons.min_free, cons.max_free);
			}
			for (auto const &nrg : back) {
				ret += nrg.penalty(cons.min_free, cons.max_free);
			}
			return ret;
		};
		std::vector< NeedleRollGoal > best_front = front;
		std::vector< NeedleRollGoal > best_back = back;
		std::vector< Transfer > best_plan;
		uint32_t best_penalty = std::numeric_limits< uint32_t >::max();
		uint32_t starting_penalty = penalty(shrunk_constraints, best_front, best_back);
		//DEBUG:
		std::cout << " ------- [penalty: " << starting_penalty << "] -------\n";
		draw_beds(BedNeedle::Back, back, BedNeedle::Front, front);
		std::cout << " --------------\n";
		{ //collapse-to-back:
			std::vector< Transfer > plan;
			std::vector< NeedleRollGoal > collapsed_top, collapsed_bottom;
			std::vector< NeedleRollGoal > shifted_top, shifted_bottom;
			std::vector< NeedleRollGoal > after_front, after_back;
			best_collapse( //collapse 'top' onto 'bottom':
				shrunk_constraints,
				BedNeedle::Front, front, //top
				BedNeedle::Back, back, //bottom
				BedNeedle::BackSliders, &collapsed_top, //after top
				BedNeedle::Back, &collapsed_bottom, //after bottom
				&plan
			);
			best_shift( //shift collapsed top/bottom to other bed:
				shrunk_constraints,
				BedNeedle::BackSliders, collapsed_top, //top
				BedNeedle::Back, collapsed_bottom, //bottom
				BedNeedle::Front, &shifted_top, //after top
				BedNeedle::FrontSliders, &shifted_bottom, //after bottom
				&plan
			);
			best_expand( //expand bottom of top/bottom back to other bed:
				shrunk_constraints,
				BedNeedle::Front, shifted_top, //top
				BedNeedle::FrontSliders, shifted_bottom, //bottom
				BedNeedle::Front, &after_front, //after top
				BedNeedle::Back, &after_back, //after bottom
				&plan
			);
			uint32_t p = penalty(shrunk_constraints, after_front, after_back);

			std::cout << "-- -- -- (penalty: " << p << ") -- -- --\n";
			draw_beds(BedNeedle::Back, after_back, BedNeedle::Front, after_front);
			std::cout << "-- -- -- -- -- --\n";
			if (p < best_penalty) {
				best_front = after_front;
				best_back = after_back;
				best_plan = plan;
				best_penalty = p;
			}
		}
		{ //collapse-to-front:
			std::vector< Transfer > plan;
			std::vector< NeedleRollGoal > collapsed_top, collapsed_bottom;
			std::vector< NeedleRollGoal > shifted_top, shifted_bottom;
			std::vector< NeedleRollGoal > after_front, after_back;
			best_collapse( //collapse 'top' onto 'bottom'
				shrunk_constraints,
				BedNeedle::Back, back, //top
				BedNeedle::Front, front, //bottom
				BedNeedle::FrontSliders, &collapsed_top, //after top
				BedNeedle::Front, &collapsed_bottom, //after bottom
				&plan
			);
			best_shift( //shift collapsed top/bottom to other bed:
				shrunk_constraints,
				BedNeedle::FrontSliders, collapsed_top, //top
				BedNeedle::Front, collapsed_bottom, //bottom
				BedNeedle::Back, &shifted_top, //after top
				BedNeedle::BackSliders, &shifted_bottom, //after bottom
				&plan
			);
			best_expand( //expand bottom of top/bottom back to other bed:
				shrunk_constraints,
				BedNeedle::Back, shifted_top, //top
				BedNeedle::BackSliders, shifted_bottom, //bottom
				BedNeedle::Back, &after_back, //after top
				BedNeedle::Front, &after_front, //after bottom
				&plan
			);
			uint32_t p = penalty(shrunk_constraints, after_front, after_back);
			std::cout << "-- -- -- (penalty: " << p << ") -- -- --\n";
			draw_beds(BedNeedle::Back, after_back, BedNeedle::Front, after_front);
			std::cout << "-- -- -- -- -- --\n";
			if (p < best_penalty) {
				best_front = after_front;
				best_back = after_back;
				best_plan = plan;
				best_penalty = p;
			}
		}
		if (!(best_penalty < starting_penalty)) {
			std::cout << "ERROR: penalty DID NOT DECREASE; you may be in for an infinite planning loop [...I think this happens because the code doesn't force zero-racking configurations after expand...]" << std::endl;
			//assert(best_penalty < starting_penalty);
		}

		transfers.insert(transfers.end(), best_plan.begin(), best_plan.end());
		front = best_front;
		back = best_back;

	} //while (not solved)

	//PARANOIA:
	//TODO: potentially, check plan

	return true;
}
