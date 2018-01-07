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
		enum : uint8_t {
			LRollInvalid = -10,
			LRoll2 = -2,
			LRoll1 = -1,
			Roll0 = 0
			RRoll1 = 1,
			RRoll2 = 2,
			RRollInvalid = 10
		};
		int8_t l_prev_roll; //{-10, -2, -1, 0} adjacent stitch is on the bottom bed to the left?
		int8_t r_next_roll; //{10, 2, 1, 0} adjacent stitch is on the bottom bed to the right?
		bool operator==(State const &o) const {
			return l == o.l
			    && r == o.r
			    && l_prev_needle == o.l_prev_needle
			    && r_next_needle == o.r_next_needle
			    && l_prev_roll == o.l_prev_roll
			    && r_next_roll == o.r_next_roll;
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


	//Let's do this in terms of the actions that can be applied:
	struct Action {
		enum Type : uint8_t {
			MoveLeft,
			MoveRight,
			RollLeft,
			RollRight,
			Roll2Left,
			Roll2Right
		} type;
		int32_t needle;

		Action(Type type_, int32_t needle_) : type(type_), needle(needle_) { }
	};

	auto apply_action = [&queue_state,&top,&bottom](Action const &action, State const &state, Cost const &cost) {
		State next_state = state;
		Cost next_cost = cost;
		if        (action.type == MoveLeft) {
			assert(state.l < top.size());

			next_cost.penalty += top[state.l].after_offset_and_roll(action.needle - top[state.l].needle, 0).penalty(constraints.min_free, constraints.max_free);

			next_state.l += 1;
			next_state.l_prev_needle = action.needle;
			next_state.l_prev_roll = State::Roll0;

			//if this is the first stitch, track it:
			if (state.r_next_roll == State::RRollInvalid) {
				assert(bottom.empty() && state.r + 1 == top.size());
				state.r_next_needle = action.needle;
				state.r_next_roll = State::RRoll2;
			}
		} else if (action.type == MoveRight) {
			assert(state.r >= 0);

			next_cost.penalty += top[state.r].after_offset_and_roll(action.needle - top[state.r].needle, 0).penalty(constraints.min_free, constraints.max_free);

			next_state.r -= 1;
			next_state.r_next_needle = action.needle;
			next_state.r_next_roll = State::Roll0;

			//if this is the first stitch, track it:
			if (state.l_prev_roll == State::LRollInvalid) {
				assert(bottom.empty() && state.l == 0);
				state.l_prev_needle = action.needle;
				state.l_prev_roll = State::RRoll2;
			}

		} else if (action.type == RollLeft) {
			assert(state.l < top.size());
			
			next_cost.penalty += top[state.l].after_offset_and_roll(action.needle - top[state.l].needle, -1).penalty(constraints.min_free, constraints.max_free);

			next_state.l += 1;
			next_state.l_prev_needle = action.needle;
			next_state.l_prev_roll = State::LRoll1;

			//if this is the first stitch, track it on the other side:
			if (state.r_next_roll == State::RRollInvalid) {
				assert(bottom.empty() && state.r + 1 == top.size());
				state.r_next_needle = action.needle;
				state.r_next_roll = State::RRoll1;
			}
		} else if (action.type == RollRight) {
			assert(state.r >= 0);
			
			next_cost.penalty += top[state.r].after_offset_and_roll(action.needle - top[state.r].needle, +1).penalty(constraints.min_free, constraints.max_free);

			next_state.r -= 1;
			next_state.r_next_needle = action.needle;
			next_state.r_next_roll = State::RRoll1;

			//if this is the first stitch, track it on the other side:
			if (state.l_prev_roll == State::LRollInvalid) {
				assert(bottom.empty() && state.l == 0);
				state.l_prev_needle = action.needle;
				state.l_prev_roll = State::LRoll1;
			}
		} else if (action.type == Roll2Left) {
			assert(state.l < top.size());
			assert(state.l_prev_roll == State::LRoll2);
			assert(top[state.l].can_stack_left); //Roll2's are always stacking
			assert(bottom.empty()); //... on already-xferred right stitches


			assert(state.r + 1 < top.size());
			assert(state.r_next_roll == State::Roll0);

			//do we add to penalty when stacking? I guess so.
			next_cost.penalty += top[state.l].after_offset_and_roll(action.needle - top[state.l].needle, -2).penalty(constraints.min_free, constraints.max_free);

			next_state.l += 1;
			next_state.l_prev_needle = action.needle;
			next_state.l_prev_roll = State::LRoll2;

			//shouldn't ever need to inform other side.

		} else if (action.type == Roll2Right) {
			assert(state.r >= 0);
			assert(state.r_next_roll == State::RRoll2);
			assert(top[state.r].can_stack_right); //Roll2's are always stacking
			assert(bottom.empty()); //... on already-xferred right stitches


			assert(state.l > 0);
			assert(state.l_prev_roll == State::Roll0);


			//do we add to penalty when stacking? I guess so.
			next_cost.penalty += top[state.r].after_offset_and_roll(action.needle - top[state.r].needle, 2).penalty(constraints.min_free, constraints.max_free);

			next_state.r -= 1;
			next_state.r_next_needle = action.needle;
			next_state.r_next_roll = State::RRoll2;

			//shouldn't ever need to inform other side.
		} else {
			assert(0 && "Unhandled action type.");
		}

		queue_state(next_state, &state, next_cost);
	};

	auto expand_state = [&](State const &state, Cost const &cost) {
		assert(state.l <= state.r);
		assert(state.r < top.size());

		assert(state.l_prev_roll <= State::Roll0);
		assert(state.r_next_roll >= State::Roll0);

		//First, and most important range: what do the current bridges, constraints, and slack allow in terms of racking?
		int32_t min_ofs = -int32_t(constraints.max_racking);
		int32_t max_ofs = int32_t(constraints.max_racking);

		if (bottom.empty() && state.l == 0 && state.r + 1 == top.size()) {
			//no bridges to worry about!
		} else {
			//can't have | ofs + top[l].needle - state.l_prev_needle | > top[l].left_slack
			//want -top[l].left_slack <= ofs + top[l].needle - state.l_prev_needle <= top[l].left_slack
			// -top[l].left_slack - (top[l].needle - state.l_prev_needle) <= ofs <= top[l].left_slack - (top[l].needle - state.l_prev_needle)
			min_ofs = std::max(min_ofs, -top[l].left_slack - (top[l].needle - state.l_prev_needle));
			max_ofs = std::min(max_ofs,  top[l].left_slack - (top[l].needle - state.l_prev_needle));
			min_ofs = std::max(min_ofs, -top[r].right_slack - (top[r].needle - state.r_next_needle));
			max_ofs = std::min(max_ofs,  top[r].right_slack - (top[r].needle - state.r_next_needle));
		}

		{ //"roll" moves for left stitch:
			int32_t l_roll_min = min_ofs;
			int32_t l_roll_max = max_ofs;

			//limit based on left stitches:
			if (state.l_prev_roll == State::LRollInvalid) {
				//nothing on the other bed, do whatever!
			} else if (state.l_prev_roll == State::LRoll2) {
				//must arrive to the left of l_prev_needle, as it's on the front:
				l_roll_max = std::min(l_roll_max, state.l_prev_needle - 1);
			} else if (state.l_prev_roll == State::LRoll1) {
				//can arrive to the left of l_prev_needle or stack:
				l_roll_max = std::min(l_roll_max, state.l_prev_needle + (top[state.l].can_stack_left ? 0 : -1));
			} else { assert(state.l_prev_roll == State::Roll0);
				//have already moved a stitch, so can't continue to roll:
				l_roll_min = std::numeric_limits< int32_t >::max();
				l_roll_max = std::numeric_limits< int32_t >::min();
			}
			//limit based on right stitches:
			if (state.r_prev_roll == State::Roll0) {
				//must be to the left of the top-bed r_prev_needle:
				l_roll_max = std::min(l_roll_max, state.r_prev_needle - 1);
			} else if (state.r_prev_roll == State::RRoll2) {
				assert(state.l_prev_roll == State::Roll0); //would need to limit, but already on front bed
			}

			for (int32_t needle = l_roll_min; needle <= l_roll_max; ++needle) {
				apply_action(Action(RollLeft, needle), state, cost);
			}

		}

		//TODO: the other move types

		//-------- OLD STUFF BELOW HERE -----------

		//TODO: clamp ofs for existing stitch locations as well (can't stack unless ... ?)

		//Next, what stitches is it okay for l to roll onto?
		int32_t l_roll_min, l_roll_max;

		if (bottom.empty() && l == 0) {
			if (state.r + 1 == top.size()) {
				//nothing has been moved, so l can roll to anywhere at all.
				l_roll_min = min_ofs;
				l_roll_max = max_ofs;
			} else {
				//nothing has moved from this side yet
				bool can_stack_l;
				if (state.l_prev_bottom) {
					can_stack_l = top[0].has_same_real_goal_as(top.back());
				} else {
					can_stack_l = false; //TODO: handle "double-roll" case where first stitch stacks with move but subsequent stitches are rolls
				}
				l_roll_min = min_ofs;
				l_roll_max = std::min(max_ofs, state.l_prev_needle + (can_stack_l ? 0 : -1));
			}
		} else if (!state.l_prev_bottom) {
			//TODO: there is a "double roll" case when bottom is empty and first r stitch was moved that is not yet handled here
			//done rolling, so use invalid range:
			l_roll_min = std::numeric_limits< int32_t >::max();
			l_roll_max = std::numeric_limits< int32_t >::min();
		} else { assert(state.l_prev_bottom);
			l_roll_min = min_ofs;
			//must roll to the left of the last stitch
			//...unless it has the same goal:
			bool can_stack_l;
			if (l == 0) {
				assert(!bottom.empty() && "handled this case above");
				NeedleRollGoal temp = bottom[0];
				temp.roll = -(temp.roll + 1);
				can_stack_l = top[0].has_same_real_goal_as(temp);
			} else {
				can_stack_l = top[l-1].has_same_real_goal_as(top[l]);
			}
			l_roll_max = std::min(max_ofs, state.l_prev_needle + (can_stack_l ? 0 : -1));
		}

		int32_t l_move_min = std::numeric_limits< int32_t >::max();
		int32_t l_move_max = std::numeric_limits< int32_t >::min();

		if (l == 0) {
			if (bottom.empty() && r + 1 == top.size()) {
			}
		}


		// <--- I WAS HERE

		/// ----- END OF NEW CODE -----

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

