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
		to_top.clear();
		to_bottom = bottom;
		return;
	}

	//Otherwise, Dijkstra-style search for optimal sequence of moves:

	//Search state:
	#pragma pack(push,1)
	struct State {
		int32_t l; //index of left stitch on top
		int32_t r; //index of right stitch on top
		int32_t l_prev_needle; //needle of stitch to the left of 'l'
		int32_t r_next_needle; //needle of stitch to the right of 'r'
		enum : int8_t {
			LRollInvalid = -10,
			LRoll2 = -2,
			LRoll1 = -1,
			Roll0 = 0,
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
		std::string to_string() const {
			return std::to_string(l_prev_needle) + "r" + std::to_string(int32_t(l_prev_roll))
				+ " [" + std::to_string(l) + "," + std::to_string(r) + "] "
				+ std::to_string(r_next_needle) + "r" + std::to_string(int32_t(r_next_roll))
			;
		}
	};
	#pragma pack(pop)
	static_assert(sizeof(State) == 4*4+2*1, "collapse's State is packed");

	struct HashState {
		size_t operator()(State const &state) const {
			static std::hash< std::string > hash;
			return hash(std::string(reinterpret_cast< char const * >(&state), reinterpret_cast< char const * >(&state) + sizeof(state)));
		}
	};

	//Let's do this in terms of the actions that can be applied:
	struct Action {
		enum Type : uint8_t {
			None,
			MoveLeft,
			MoveRight,
			RollLeft,
			RollRight,
			Roll2Left,
			Roll2Right,
		} type;
		int32_t needle;

		Action(Type type_, int32_t needle_) : type(type_), needle(needle_) { }

		std::string to_string() const {
			if (type == None) return "CNone";
			else if (type == MoveLeft) return "CMoveLeft to " + std::to_string(needle);
			else if (type == MoveRight) return "CMoveRight to " + std::to_string(needle);
			else if (type == RollLeft) return "CRollLeft to " + std::to_string(needle);
			else if (type == RollRight) return "CRollRight to " + std::to_string(needle);
			else if (type == Roll2Left) return "CRoll2Left to " + std::to_string(needle);
			else if (type == Roll2Right) return "CRoll2Right to " + std::to_string(needle);
			else {
				assert(0 && "invalid move type");
				return "!"; //never reached
			}
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

	struct StateInfo {
		Cost cost;
		State const *source;
		Action action;

		StateInfo(Cost const &cost_, State const *source_, Action const &action_) : cost(cost_), source(source_), action(action_) { }
	};

	std::multimap< Cost, const State * > todo;
	std::unordered_map< State, StateInfo, HashState > best_source;

	auto queue_state = [&](State const &state, Cost const &cost, State const *from, Action const &action) {
		auto ret = best_source.insert(std::make_pair(state, StateInfo(cost, from, action)));
		if (cost < ret.first->second.cost) {
			ret.first->second = StateInfo(cost, from, action);
			ret.second = true;
		}
		if (ret.second) {
			todo.insert(std::make_pair( cost, &ret.first->first ));
		}
	};


	auto apply_action = [&queue_state,&top,&bottom,&constraints](Action const &action, State const &state, Cost const &cost) {
		//std::cout << "  doing '" << action.to_string() << "'" << std::endl; //DEBUG
		State next_state = state;
		Cost next_cost = cost;
		if        (action.type == Action::MoveLeft) {
			assert(state.l >= 0 && state.l < int32_t(top.size()));

			next_cost.penalty += top[state.l].after_offset_and_roll(action.needle - top[state.l].needle, 0).penalty(constraints.min_free, constraints.max_free);

			next_state.l += 1;
			next_state.l_prev_needle = action.needle;
			next_state.l_prev_roll = State::Roll0;

			//if this is the first stitch, track it:
			if (state.r_next_roll == State::RRollInvalid) {
				assert(bottom.empty() && state.r + 1 == int32_t(top.size()));
				next_state.r_next_needle = action.needle;
				next_state.r_next_roll = State::RRoll2;
			}
		} else if (action.type == Action::MoveRight) {
			assert(state.r >= 0 && state.r < int32_t(top.size()));

			next_cost.penalty += top[state.r].after_offset_and_roll(action.needle - top[state.r].needle, 0).penalty(constraints.min_free, constraints.max_free);

			next_state.r -= 1;
			next_state.r_next_needle = action.needle;
			next_state.r_next_roll = State::Roll0;

			//if this is the first stitch, track it:
			if (state.l_prev_roll == State::LRollInvalid) {
				assert(bottom.empty() && state.l == 0);
				next_state.l_prev_needle = action.needle;
				next_state.l_prev_roll = State::LRoll2;
			}

		} else if (action.type == Action::RollLeft) {
			assert(state.l >= 0 && state.l < int32_t(top.size()));
			
			next_cost.penalty += top[state.l].after_offset_and_roll(action.needle - top[state.l].needle, -1).penalty(constraints.min_free, constraints.max_free);

			next_state.l += 1;
			next_state.l_prev_needle = action.needle;
			next_state.l_prev_roll = State::LRoll1;

			//if this is the first stitch, track it on the other side:
			if (state.r_next_roll == State::RRollInvalid) {
				assert(bottom.empty() && state.r + 1 == int32_t(top.size()));
				next_state.r_next_needle = action.needle;
				next_state.r_next_roll = State::RRoll1;
			}
		} else if (action.type == Action::RollRight) {
			assert(state.r >= 0 && state.r < int32_t(top.size()));
			
			next_cost.penalty += top[state.r].after_offset_and_roll(action.needle - top[state.r].needle, +1).penalty(constraints.min_free, constraints.max_free);

			next_state.r -= 1;
			next_state.r_next_needle = action.needle;
			next_state.r_next_roll = State::RRoll1;

			//if this is the first stitch, track it on the other side:
			if (state.l_prev_roll == State::LRollInvalid) {
				assert(bottom.empty() && state.l == 0);
				next_state.l_prev_needle = action.needle;
				next_state.l_prev_roll = State::LRoll1;
			}
		} else if (action.type == Action::Roll2Left) {
			assert(state.l >= 0 && state.l < int32_t(top.size()));
			assert(bottom.empty());
			assert(state.l_prev_roll == State::LRoll2 || state.l_prev_roll == State::LRollInvalid);

			//do we add to penalty when stacking? I guess so.
			next_cost.penalty += top[state.l].after_offset_and_roll(action.needle - top[state.l].needle, -2).penalty(constraints.min_free, constraints.max_free);

			next_state.l += 1;
			next_state.l_prev_needle = action.needle;
			next_state.l_prev_roll = State::LRoll2;

			if (state.r_next_roll == State::RRollInvalid) {
				next_state.r_next_needle = action.needle;
				next_state.r_next_roll = State::Roll0;
			}

		} else if (action.type == Action::Roll2Right) {
			assert(state.r >= 0 && state.r < int32_t(top.size()));
			assert(bottom.empty());
			assert(state.r_next_roll == State::RRoll2 || state.r_next_roll == State::RRollInvalid);

			//do we add to penalty when stacking? I guess so.
			next_cost.penalty += top[state.r].after_offset_and_roll(action.needle - top[state.r].needle, 2).penalty(constraints.min_free, constraints.max_free);

			next_state.r -= 1;
			next_state.r_next_needle = action.needle;
			next_state.r_next_roll = State::RRoll2;

			if (state.l_prev_roll == State::RRollInvalid) {
				next_state.l_prev_needle = action.needle;
				next_state.l_prev_roll = State::Roll0;
			}


			//shouldn't ever need to inform other side.
		} else {
			assert(0 && "Unhandled action type.");
		}

		queue_state(next_state, next_cost, &state, action);
	};

	auto expand_state = [&](State const &state, Cost const &cost) {
		assert(state.l <= state.r);
		assert(state.r < int32_t(top.size()));

		assert(state.l_prev_roll <= State::Roll0);
		assert(state.r_next_roll >= State::Roll0);

		std::vector< std::pair< int32_t, Action > > offset_action;

		//First, and most important range: what do the current bridges, constraints, and slack allow in terms of racking?
		int32_t min_ofs = -int32_t(constraints.max_racking);
		int32_t max_ofs = int32_t(constraints.max_racking);

		if (bottom.empty() && state.l == 0 && state.r + 1 == int32_t(top.size())) {
			//no bridges to worry about!
			assert(state.l_prev_roll == State::LRollInvalid && state.r_next_roll == State::RRollInvalid);
		} else {
			//can't have | ofs + top[l].needle - state.l_prev_needle | > top[l].left_slack
			//want -top[l].left_slack <= ofs + top[l].needle - state.l_prev_needle <= top[l].left_slack
			// -top[l].left_slack - (top[l].needle - state.l_prev_needle) <= ofs <= top[l].left_slack - (top[l].needle - state.l_prev_needle)
			if (top[state.l].left_slack != SlackForNoYarn) {
				min_ofs = std::max(min_ofs, -top[state.l].left_slack - (top[state.l].needle - state.l_prev_needle));
				max_ofs = std::min(max_ofs,  top[state.l].left_slack - (top[state.l].needle - state.l_prev_needle));
			}
			if (top[state.r].right_slack != SlackForNoYarn) {
				min_ofs = std::max(min_ofs, -top[state.r].right_slack - (top[state.r].needle - state.r_next_needle));
				max_ofs = std::min(max_ofs,  top[state.r].right_slack - (top[state.r].needle - state.r_next_needle));
			}
		}

		{ //"roll2" moves for left stitch:
			int32_t min = min_ofs + top[state.l].needle;
			int32_t max = max_ofs + top[state.l].needle;
			if (state.l_prev_roll == State::LRollInvalid) {
				//can roll2 to ~anywhere~
			} else if (state.l_prev_roll == State::LRoll2) {
				min = std::max(min, state.l_prev_needle + (top[state.l].can_stack_left ? 0 : +1));
				max = std::min(max, state.r_next_needle);
			} else {
				min = std::numeric_limits< int32_t >::max();
				max = std::numeric_limits< int32_t >::min();
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				offset_action.emplace_back(needle - top[state.l].needle, Action(Action::Roll2Left, needle));
				//apply_action(Action(Action::Roll2Left, needle), state, cost);
			}
		}

		{ //"roll2" moves for right stitch:
			int32_t min = min_ofs + top[state.r].needle;
			int32_t max = max_ofs + top[state.r].needle;
			if (state.r_next_roll == State::LRollInvalid) {
				//can roll2 to ~anywhere~
			} else if (state.r_next_roll == State::RRoll2) {
				min = std::max(min, state.l_prev_needle);
				max = std::min(max, state.r_next_needle + (top[state.r].can_stack_right ? 0 : -1));
			} else {
				min = std::numeric_limits< int32_t >::max();
				max = std::numeric_limits< int32_t >::min();
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				offset_action.emplace_back(needle - top[state.r].needle, Action(Action::Roll2Right, needle));
				//apply_action(Action(Action::Roll2Right, needle), state, cost);
			}
		}

		{ //"roll" moves for left stitch:
			int32_t min = min_ofs + top[state.l].needle;
			int32_t max = max_ofs + top[state.l].needle;

			//limit based on left stitches:
			if (state.l_prev_roll == State::LRollInvalid) {
				//nothing on the other bed, do whatever!
			} else if (state.l_prev_roll == State::LRoll2) {
				//must arrive to the left of l_prev_needle, as it's on the front:
				max = std::min(max, state.l_prev_needle - 1);
			} else if (state.l_prev_roll == State::LRoll1) {
				//can arrive to the left of l_prev_needle or stack:
				max = std::min(max, state.l_prev_needle + (top[state.l].can_stack_left ? 0 : -1));
			} else { assert(state.l_prev_roll == State::Roll0);
				//have already moved a stitch, so can't continue to roll:
				min = std::numeric_limits< int32_t >::max();
				max = std::numeric_limits< int32_t >::min();
			}
			//limit based on right stitches:
			if (state.r_next_roll == State::Roll0) {
				//must be to the left of the top-bed r_prev_needle:
				max = std::min(max, state.r_next_needle - 1);
			} else if (state.r_next_roll == State::RRoll2) {
				assert(state.l_prev_roll == State::Roll0); //would need to limit, but already on front bed
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				offset_action.emplace_back(needle - top[state.l].needle, Action(Action::RollLeft, needle));
				//apply_action(Action(Action::RollLeft, needle), state, cost);
			}
		}

		{ //"roll" moves for right stitch:
			int32_t min = min_ofs + top[state.r].needle;
			int32_t max = max_ofs + top[state.r].needle;

			//limit based on right stitches:
			if (state.r_next_roll == State::RRollInvalid) {
				//nothing on the other bed, do whatever!
			} else if (state.r_next_roll == State::RRoll2) {
				//must arrive to the right of r_next_needle, as it's on the front:
				min = std::max(min, state.r_next_needle + 1);
			} else if (state.r_next_roll == State::RRoll1) {
				//can arrive to the right of r_prev_needle or stack:
				min = std::max(min, state.r_next_needle + (top[state.r].can_stack_right ? 0 : 1));
			} else { assert(state.r_next_roll == State::Roll0);
				//have already moved a stitch, so can't continue to roll:
				min = std::numeric_limits< int32_t >::max();
				max = std::numeric_limits< int32_t >::min();
			}

			//limit based on left stitches:
			if (state.l_prev_roll == State::Roll0) {
				//must be to the right of the top-bed l_prev_needle:
				min = std::max(min, state.l_prev_needle + 1);
			} else if (state.l_prev_roll == State::LRoll2) {
				assert(state.r_next_roll == State::Roll0); //would need to limit, but already on front bed
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				offset_action.emplace_back(needle - top[state.r].needle, Action(Action::RollRight, needle));
				//apply_action(Action(Action::RollRight, needle), state, cost);
			}
		}

		{ //"move" moves for left stitch:
			int32_t min = min_ofs + top[state.l].needle;
			int32_t max = max_ofs + top[state.l].needle;

			//limit based on left stitches:
			if (state.l_prev_roll == State::LRollInvalid) {
				//nothing on the other bed, do whatever!
			} else if (state.l_prev_roll == State::LRoll2) {
				//must arrive to the left of l_prev_needle, as it's on the front; but r_next_needle case below will deal with this
			} else if (state.l_prev_roll == State::LRoll1) {
				//nothing to interfere on the left, do whatever!
			} else { assert(state.l_prev_roll == State::Roll0);
				//must place to the right of existing front bed stuff:
				min = std::max(min, state.l_prev_needle + (top[state.l].can_stack_left ? 0 : 1));
			}
			//limit based on right stitches:
			if (state.r_next_roll == State::Roll0) {
				//must be to the left of the top-bed r_prev_needle:
				max = std::min(max, state.r_next_needle - 1);
			} else if (state.r_next_roll == State::RRoll2) {
				assert(state.l_prev_roll == State::Roll0); //would need to limit, but already on front bed
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				offset_action.emplace_back(needle - top[state.l].needle, Action(Action::MoveLeft, needle));
				//apply_action(Action(Action::MoveLeft, needle), state, cost);
			}
		}

		{ //"move" moves for right stitch:
			int32_t min = min_ofs + top[state.r].needle;
			int32_t max = max_ofs + top[state.r].needle;

			//limit based on left stitches:
			if (state.r_next_roll == State::RRollInvalid) {
				//nothing on the other bed, do whatever!
			} else if (state.r_next_roll == State::RRoll2) {
				//must arrive to the right of r_next_needle, as it's on the front; but l_prev_needle case below will deal with this
			} else if (state.r_next_roll == State::RRoll1) {
				//nothing to interfere on the left, do whatever!
			} else { assert(state.r_next_roll == State::Roll0);
				//must place to the left of existing front bed stuff:
				max = std::min(max, state.r_next_needle + (top[state.r].can_stack_right ? 0 : -1));
			}
			//limit based on left stitches:
			if (state.l_prev_roll == State::Roll0) {
				//must be to the right of the top-bed l_prev_needle:
				min = std::max(min, state.l_prev_needle + 1);
			} else if (state.l_prev_roll == State::LRoll2) {
				assert(state.r_next_roll == State::Roll0); //would need to limit, but already on front bed
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				offset_action.emplace_back(needle - top[state.r].needle, Action(Action::MoveRight, needle));
				//apply_action(Action(Action::MoveRight, needle), state, cost);
			}
		}

		std::stable_sort(offset_action.begin(), offset_action.end(), [](std::pair< int32_t, Action > const &a, std::pair< int32_t, Action > const &b) -> bool {
			return std::abs(a.first) < std::abs(b.first);
		});

		for (auto const &oa : offset_action) {
			apply_action(oa.second, state, cost);
		}
	};

	//Initial state(s):
	{
		//NOTE: if bottom is empty, the prev_needle / next_needle fields don't matter until after first xfer
		State init;
		init.l = 0;
		init.l_prev_needle = (bottom.empty() ? top[0].needle : bottom[0].needle);
		init.l_prev_roll = (bottom.empty() ? State::LRollInvalid : State::LRoll1);
		init.r = top.size()-1;
		init.r_next_needle = (bottom.empty() ? top.back().needle : bottom.back().needle);
		init.r_next_roll = (bottom.empty() ? State::RRollInvalid : State::RRoll1);
		Cost cost;
		cost.penalty = 0;
		queue_state(init, cost, nullptr, Action(Action::None, 0));
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
			if (f->second.cost < cost) continue;
			assert(f->second.cost == cost);
		}
		//std::cout << "Considering " << state->to_string() << std::endl; //DEBUG

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
		auto f = best_source.find(*best);
		assert(f != best_source.end());
		if (f->second.source == nullptr) break;
		State const &state = *f->second.source;
		Action const &action = f->second.action;

		if (action.type == Action::MoveLeft) {
			assert(state.l >= 0 && state.l < int32_t(top.size()));
			ops.emplace_back(BedNeedle(top_bed, top[state.l].needle), BedNeedle(to_top_bed, action.needle));
		} else if (action.type == Action::MoveRight) {
			assert(state.r >= 0 && state.r < int32_t(top.size()));
			ops.emplace_back(BedNeedle(top_bed, top[state.r].needle), BedNeedle(to_top_bed, action.needle));
		} else if (action.type == Action::RollLeft) {
			assert(state.l >= 0 && state.l < int32_t(top.size()));
			ops.emplace_back(BedNeedle(top_bed, top[state.l].needle), BedNeedle(to_bottom_bed, action.needle));
		} else if (action.type == Action::RollRight) {
			assert(state.r >= 0 && state.r < int32_t(top.size()));
			ops.emplace_back(BedNeedle(top_bed, top[state.r].needle), BedNeedle(to_bottom_bed, action.needle));
		} else if (action.type == Action::Roll2Left) {
			assert(state.l >= 0 && state.l < int32_t(top.size()));
			ops.emplace_back(BedNeedle(top_bed, top[state.l].needle), BedNeedle(to_top_bed, action.needle));
		} else if (action.type == Action::Roll2Right) {
			assert(state.r >= 0 && state.r < int32_t(top.size()));
			ops.emplace_back(BedNeedle(top_bed, top[state.r].needle), BedNeedle(to_top_bed, action.needle));
		} else {
			assert(0 && "Invalid action type.");
		}
		ops.back().why = state.to_string() + "; " + action.to_string();

		best = f->second.source;
	}

	std::reverse(ops.begin(), ops.end());

/*
	std::cout << "  Final plan:\n"; //DEBUG
	for (auto const &op : ops) {
		std::cout << "    " << op.to_string() << '\n';
	}
	std::cout.flush(); //DEBUG
	*/

/*
	std::cout << "Before Collapse:\n"; //DEBUG
	draw_beds(top_bed, top, bottom_bed, bottom); //DEBUG
*/
	run_transfers(constraints,
		top_bed, top,
		bottom_bed, bottom,
		ops,
		to_top_bed, &to_top,
		to_bottom_bed, &to_bottom);

/*
	std::cout << "After Collapse:\n"; //DEBUG
	draw_beds(to_top_bed, to_top, to_bottom_bed, to_bottom); //DEBUG
*/

	plan.insert(plan.end(), ops.begin(), ops.end());

}
