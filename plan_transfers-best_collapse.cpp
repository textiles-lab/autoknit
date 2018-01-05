#include "plan_transfers-helpers.hpp"

#include <iostream>
#include <map>
#include <unordered_map>

void best_collapse(
	Constraints const &constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const &top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const &bottom,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal > *to_top_,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal > *to_bottom_,
	std::vector< Transfer > *plan_
) {
	//Collapse won't change the bottom bed's location, but will change the top's:
	assert(top_bed != to_top_bed);
	assert(bottom_bed == to_bottom_bed);

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
	//collapse moves stitches from top to bottom, from the edges in.

	//if no stitches on top, nothing to move, so done:
	if (top.empty()) {
		return;
	}

	//Otherwise, A*-style search for optimal sequence of moves:

	//Search state:
	#pragma pack(push,1)
	struct State {
		uint32_t l; //index of left stitch on top
		uint32_t r; //index of right stitch on top
		int32_t l_prev_needle; //needle of stitch to the left of 'l'
		int32_t r_next_needle; //needle of stitch to the right of 'r'
		bool l_prev_bottom; //adjacent stitch is on the bottom bed to the left?
		bool r_next_bottom; //adjacent stitch is on the bottom bed to the right?
		bool operator==(State const &o) const {
			return l == o.l
			    && r == o.r
			    && l_prev_needle == o.l_prev_needle
			    && r_next_needle == o.r_next_needle
			    && l_prev_bottom == o.l_prev_bottom
			    && r_next_bottom == o.r_next_bottom;
		};
	};
	#pragma pack(pop)
	static_assert(sizeof(State) == 4*4+2*1, "collapse's State is packed");

	struct HashState {
		size_t operator()(State const &state) const {
			static std::hash< std::string > hash;
			return hash(std::string(reinterpret_cast< char const * >(&state), reinterpret_cast< char const * >(&state) + sizeof(state)));
		}
	};

	struct Cost {
		uint32_t penalty = 0;
		bool operator<(Cost const &o) const {
			return penalty < o.penalty;
		}
		bool operator==(Cost const &o) const {
			return penalty == o.penalty;
		}
	};

	std::map< Cost, const State * > todo;
	std::unordered_map< State, std::pair< const State *, Cost >, HashState > best_source;

	auto queue_state = [&](State const &state, const State *from, Cost const &cost) {
		auto ret = best_source.insert(std::make_pair(state, std::make_pair(from, cost)));
		if (cost < ret.first->second.second) {
			ret.first->second = std::make_pair(from, cost);
			ret.second = true;
		}
		if (ret.second) {
			todo.insert(std::make_pair( cost, &ret.first->first ));
		}
	};

	//OKAY! This is getting into special case hell, which is best avoided.
	// So, things to think about:
	// the min/max offset range is limited by bridges as long as *any* stitch has been moved
	// bottom empty, no stitches:
	//   can move either left or right, roll or no roll
	// bottom empty, only left stitches moved:
	//   left moves are "usual", right can still roll or no roll, as long as it is to the right of the last left stitch.
	// bottom empty, only right stitches moved:
	//   right moves are "usual", left can still roll or no roll, as long as it is to the left of the last right stitch moved
	// These two bottom cases seem to be pretty much the same as "usual" moves, with modified roll bounds.

	// Can probably unify these cases by using move and roll ranges, with limited special-casing when setting the ranges up:
	//  l--o--o--r
	// [      ] <- range that is okay to "move" into (limits: previously moved stitch, if any; space required -- for early-out)
	// [   ]   <- range that it is okay to "roll" into (limits: previously rolled stitch, if any; bottom bed stitches, if any)
	//     p-- ... ---n
	   

	//Special-case for first expansion when bottom is empty;
	// this is done so that in general there are already two stitches on the bottom:
	auto special_expand_state = [&](State const &state, Cost const &cost) {
		assert(bottom.empty() && l == 0 && r + 1 == top.size());

		//option (0): only one stitch:
		if (state.l == state.r) {
			assert(0 && "TODO: the only-one-stitch initial expand case");
		}
		//option (1): move l and then move r:
		if (state.l != state.r)
		for (int32_t ofs_l = -int32_t(constraints.max_racking); ofs_l <= int32_t(constraints.max_racking); ++ofs_l) {
			for (int32_t ofs_r = -int32_t(constraints.max_racking); ofs_r <= int32_t(constraints.max_racking); ++ofs_r) {
			//move l:
				if (constraints.min_free <= top[l].needle + ofs && top[l].needle + ofs <= constraints.max_free) {
					State next_state;
					next_state.l = state.l + 1;
					next_state.r = state.r;
					next_state.r_next_needle = next_state.l_prev_needle = top[l].needle + ofs;

					{ //roll to bottom:
						next_state.r_next_bottom = next_state.l_prev_bottom = true;
						Cost next_cost = cost;
						NeedleRollGoal nrg = top[l];
						nrg.needle += ofs;
						nrg.roll = -(nrg.roll + 1);
						next_cost.penalty += nrg.penalty(constraints.min_free, constraints.max_free);
						queue_state(next_state, &state, next_cost);
					}

					{ //move to top:
						next_state.r_next_bottom = next_state.l_prev_bottom = false;
						Cost next_cost = cost;
						NeedleRollGoal nrg = top[l];
						nrg.needle += ofs;
						next_cost.penalty += nrg.penalty(constraints.min_free, constraints.max_free);
						queue_state(next_state, &state, next_cost);
					}
				}

				//move r:
				if (constraints.min_free <= top[r].needle + ofs && top[r].needle + ofs <= constraints.max_free) {
					State next_state;
					next_state.l = state.l;
					next_state.r = state.r - 1;
					next_state.r_next_needle = next_state.l_prev_needle = top[r].needle + ofs;

					{ //roll to bottom:
						next_state.r_next_bottom = next_state.l_prev_bottom = true;
						Cost next_cost = cost;
						NeedleRollGoal nrg = top[r];
						nrg.needle += ofs;
						nrg.roll = -(nrg.roll - 1);
						next_cost.penalty += nrg.penalty(constraints.min_free, constraints.max_free);
						queue_state(next_state, &state, next_cost);
					}

					{ //move to top:
						next_state.r_next_bottom = next_state.l_prev_bottom = false;
						Cost next_cost = cost;
						NeedleRollGoal nrg = top[r];
						nrg.needle += ofs;
						next_cost.penalty += nrg.penalty(constraints.min_free, constraints.max_free);
						queue_state(next_state, &state, next_cost);
					}
				}

			}

	}

	auto expand_state = [&](State const &state, Cost const &cost) {
		assert(state.l <= state.r);
		assert(state.r < top.size());
		if (bottom.empty() && l == 0 && r + 1 == top.size()) {
			special_expand_state(state, cost);
			return;
		}
		if (bottom.empty()) {
			assert(l > 0 && r + 1 < top.size()); //moved at least two stitches already.
		}

		//General case collapse moves:

		//handy to know if it's okay to stack up l/r with their previous stitches:
		//does the stitch at top[l] have the same [eventual] goal as the stitch to its left?
		bool can_stack_l;
		if (state.l > 0) {
			can_stack_l = top[state.l].has_same_real_goal_as(top[state.l-1]);
		} else { assert(state.l == 0);
			//the first stitch needs special treatment, because the thing to its left is maybe on a different bed:
			if (bottom.empty()) {
				//well, if no bottom stitches, still easy to compare:
				can_stack_l = top[0].has_same_real_goal_as(top.back());
			} else {
				//roll bottom stitch around left to compare:
				NeedleRollGoal temp = bottom[0];
				temp.roll = -(temp.roll + 1);
				can_stack_l = top[0].has_same_real_goal_as(temp);
			}
		}

		bool can_stack_r;
		if (state.r + 1 < top.size()) {
			can_stack_r = top[state.r].has_same_real_goal_as(top[state.r+1]);
		} else { assert(state.r + 1 == top.size());
			//the last stitch needs special treatment, because the thing to its left is maybe on a different bed:
			if (bottom.empty()) {
				//well, if no bottom stitches, still easy to compare:
				can_stack_r = top.back().has_same_real_goal_as(top[0]);
			} else {
				//roll bottom stitch around right to compare:
				NeedleRollGoal temp = bottom.back();
				temp.roll = -(temp.roll - 1);
				can_stack_r = top.back().has_same_real_goal_as(temp);
			}
		}

			for (int32_t ofs = -int32_t(constraints.max_racking); ofs <= int32_t(constraints.max_racking); ++ofs) {
				//what offsets one xfer at?
				if (!no_limits && std::abs(top[state.l].needle + ofs - state.l_prev_needle) > top[state.l].left_slack) continue;
				if (!no_limits && std::abs(top[state.r].needle + ofs - state.r_next_needle) > top[state.r].right_slack) continue;

				//collapse l:
				if (constraints.min_free <= top[l].needle + ofs && top[l].needle + ofs <= constraints.max_free) {
					State next_state;
					next_state.l = state.l + 1;
					next_state.r = state.r;
					next_state.l_prev_needle = top[l].needle + ofs;

					//roll to bottom:
					if (state.l_prev_bottom && ( //must have rolled previous stitch as well, and...
						state.l_prev_needle > next_state.l_prev_needle //continue moving left
						|| (state.l_prev_needle == next_state.l_prev_needle && can_stack_l) //or stack
					)) {
						next_state.l_prev_bottom = true;
						Cost next_cost = cost;
						NeedleRollGoal nrg = top[l];
						nrg.needle = next_state.l_prev_needle;
						nrg.roll = -(nrg.roll + 1);
						next_cost.penalty += nrg.penalty(constraints.min_free, constraints.max_free);
						queue_state(next_state, &state, next_cost);
					}
					//move to top:
					if ((state.l_prev_bottom //if last was on bottom, okay to put anywhere
					 || state.l_prev_needle < next_state.l_prev_needle //if last was on top must go right
					 || (state.l_prev_needle == next_state.l_prev_needle && can_stack_l) //or stack
					) && ( //make sure not moving past right edge
						state.r_prev_bottom //..if right is still on the bottom, no problem
						|| state.r + 1 == top.size() //..if right still hasn't moved first, no problem
						|| next_state.l_prev_needle < next_state.r_prev_needle
					) { //move to top:
						next_state.l_prev_bottom = false;
						Cost next_cost = cost;
						NeedleRollGoal nrg = top[l];
						nrg.needle = next_state.l_prev_needle;
						next_cost.penalty += nrg.penalty(constraints.min_free, constraints.max_free);
						queue_state(next_state, &state, next_cost);
					}
				}

				//move r:
				if (constraints.min_free <= top[l].needle + ofs && top[l].needle + ofs <= constraints.max_free) {
					State next_state;
					next_state.l = state.l + 1;
					next_state.r = state.r;
					next_state.l_prev_needle = top[l].needle + ofs;

					if (state.l_prev_bottom) { //roll to bottom:
						next_state.l_prev_bottom = true;
						Cost next_cost = cost;
						NeedleRollGoal nrg = top[l];
						nrg.needle = next_state.l_prev_needle;
						nrg.roll = -(nrg.roll + 1);
						next_cost.penalty += nrg.penalty(constraints.min_free, constraints.max_free);
						queue_state(next_state, &state, next_cost);
					}

					{ //move to top:
						next_state.l_prev_bottom = false;
						Cost next_cost = cost;
						NeedleRollGoal nrg = top[l];
						nrg.needle = next_state.l_prev_needle;
						next_cost.penalty += nrg.penalty(constraints.min_free, constraints.max_free);
						queue_state(next_state, &state, next_cost);
					}
				}

			}
		}
	};

	//Initial state(s):
	{
		//NOTE: if bottom is empty, the prev_needle / next_needle fields don't matter until after first xfer
		State init;
		init.l = 0;
		init.l_prev_needle = (bottom.empty() ? top[0].needle : bottom[0].needle);
		init.l_prev_bottom = true;
		init.r = top.size()-1;
		init.r_next_needle = (bottom.empty() ? top.back().needle : bottom.back().needle);
		init.r_next_bottom = true;
		Cost cost;
		cost.penalty = 0;
		queue_state(init, nullptr, cost);
	}

	//Actual search:
	const State *best = nullptr;
	while (!todo.empty()) {
		Cost cost = todo.begin()->first;
		const State *state = todo.begin()->second;
		todo.erase(todo.begin());
		{ //see if this is the first time the state has been expanded:
			auto f = best_source.find(*state);
			assert(f != best_source.end());
			assert(&(f->first) == state);
			if (f->second.second < cost) continue;
			assert(f->second.second == cost);
		}
		//if this is an ending state, end:
		if (state->l > state->r) {
			best = state;
			break;
		}
		//otherwise, expand:
		expand_state(*state, cost);
	}
	assert(best && "Must have gotten to some ending state.");

	//read back operations from best:
	std::vector< Transfer > ops;
	while (best) {
		const State *from;
		{
			auto f = best_source.find(*best);
			assert(f != best_source.end());
			from = f->second.first;
		}
		if (from == nullptr) break;

		//reconstruct operation:
		ops.emplace_back();
		if (best->l != from->l) {
			assert(best->l == from->l + 1);
			ops.back().from.bed = top_bed;
			ops.back().from.needle = top[from->l].needle;
			ops.back().to.needle = best->l_prev_needle;
			ops.back().to.bed = (best->l_prev_bottom ? to_bottom_bed : to_top_bed);
		} else { assert(best->r != from->r);
			assert(best->r == from->r - 1);
			ops.back().from.bed = top_bed;
			ops.back().from.needle = top[from->r].needle;
			ops.back().to.needle = best->r_next_needle;
			ops.back().to.bed = (best->r_next_bottom ? to_bottom_bed : to_top_bed);
		}

		best = from;
	}

	plan.insert(plan.end(), ops.rbegin(), ops.rend()); //reverse ops + add to plan

}

