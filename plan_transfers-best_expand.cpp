#include "plan_transfers-helpers.hpp"

#include <iostream>
#include <map>
#include <unordered_map>


void best_expand(
	Constraints const &constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const &top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const &bottom,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal > *to_top_,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal > *to_bottom_,
	std::vector< Transfer > *plan_
) {
	//Expand won't change the top bed's location, but will change the bottom bed's:
	assert(top_bed == to_top_bed);
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

	//if no stitches on bottom, nothing to move, so done:
	if (bottom.empty()) {
		to_top = top;
		to_bottom.clear();
		return;
	}

	//Otherwise, Dijkstra-style search for optimal sequence of moves:

	//Moves in 'expand' move a stitch from the bottom bed (or roll a stitch from the top bed) to the to_bottom bed

	//Search state:
	#pragma pack(push,1)
	struct State {
		//if 'l' or 'r' is < 0 or >= bottom.size(), this means that some top stitches are getting rolled
		int32_t l; //index of left stitch
		int32_t r; //index of right stitch
		int32_t l_next_needle; //needle that the stitch at index l+1 'l' was moved to
		int32_t r_prev_needle; //needle that the stitch at index r-1 was moved to
		bool operator==(State const &o) const {
			return l == o.l
			    && r == o.r
			    && l_next_needle == o.l_next_needle
			    && r_prev_needle == o.r_prev_needle
			;
		};
		std::string to_string() const {
			return std::to_string(l) + " " + std::to_string(l_next_needle) + "." + std::to_string(r_prev_needle) + " " + std::to_string(r);
		}
	};
	#pragma pack(pop)
	static_assert(sizeof(State) == 4*4, "expand's State is packed");

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
			Finish //charge for all the remaining top stitches
		} type;
		int32_t needle;

		Action(Type type_, int32_t needle_) : type(type_), needle(needle_) { }

		std::string to_string() const {
			if (type == None) return "ENone";
			else if (type == MoveLeft) return "EMoveLeft to " + std::to_string(needle);
			else if (type == MoveRight) return "EMoveRight to " + std::to_string(needle);
			else if (type == Finish) return "EFinish";
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

		int32_t top_l = -(state.l + 1);
		int32_t top_r = int32_t(top.size()) - 1 - (state.r - int32_t(bottom.size()));

		assert(state.l <= state.r);
		assert(top_l <= top_r);

		if        (action.type == Action::MoveLeft) {
			if (state.l >= 0) {
				//l is on bottom bed
				next_cost.penalty += bottom[state.l].after_offset_and_roll(action.needle - bottom[state.l].needle, 0).penalty(constraints.min_free, constraints.max_free);

				next_state.l -= 1;
				next_state.l_next_needle = action.needle;

				if (state.r == state.l) {
					next_state.r += 1;
					next_state.r_prev_needle = action.needle;
				}
			} else { assert(state.l < 0);
				//l has wrapped to top bed
				assert(top_l >= 0 && uint32_t(top_l) < top.size());

				next_cost.penalty += top[top_l].after_offset_and_roll(action.needle - top[top_l].needle, -1).penalty(constraints.min_free, constraints.max_free);

				next_state.l -= 1;
				next_state.l_next_needle = action.needle;

				//Ignore top_l == top_r case.
			}
		} else if (action.type == Action::MoveRight) {
			if (state.r < int32_t(bottom.size())) {
				//r is on bottom bed
				assert(state.r >= 0);
				next_cost.penalty += bottom[state.r].after_offset_and_roll(action.needle - bottom[state.r].needle, 0).penalty(constraints.min_free, constraints.max_free);

				next_state.r += 1;
				next_state.r_prev_needle = action.needle;

				if (state.l == state.r) {
					next_state.l -= 1;
					next_state.l_next_needle = action.needle;
				}
			} else { assert(state.r >= int32_t(bottom.size()));
				//r has wrapped to top bed
				assert(top_r >= 0 && uint32_t(top_r) < top.size());

				next_cost.penalty += top[top_r].after_offset_and_roll(action.needle - top[top_r].needle, 1).penalty(constraints.min_free, constraints.max_free);

				next_state.r += 1;
				next_state.r_prev_needle = action.needle;

				//Ignore top_l == top_r case; as it's already a final state.
			}
		} else if (action.type == Action::Finish) {
			//got to have done all bottom stitches:
			assert(top_l >= 0 && top_r < int32_t(top.size()));

			//charge for all the not-done top stitches:
			for (int32_t i = top_l; i <= top_r; ++i) {
				next_cost.penalty += top[i].penalty(constraints.min_free, constraints.max_free);
			}

			//cross the indices:
			next_state.l = -1 - int32_t(top.size());
			next_state.r = bottom.size() + top.size();

			assert( -(next_state.l + 1) == int32_t(top.size()) );
			assert( int32_t(top.size()) - 1 - (next_state.r - int32_t(bottom.size())) == -1 );
			assert(next_state.r - next_state.l >= int32_t(top.size() + bottom.size()) );

		} else {
			assert(0 && "Invalid action type");
		}
		//std::cout << "     penalty " << next_cost.penalty  << " from " << cost.penalty << std::endl; //DEBUG
		assert(next_cost.penalty >= cost.penalty);

		queue_state(next_state, next_cost, &state, action);
	};

	auto expand_state = [&](State const &state, Cost const &cost) {
		assert(state.l <= state.r);

		int32_t top_l = -(state.l + 1);
		int32_t top_r = int32_t(top.size()) - 1 - (state.r - int32_t(bottom.size()));

		assert(top_l <= top_r); //would have been flagged a final state otherwise

		//First, and most important range: what do the current bridges, constraints, and slack allow in terms of racking?
		int32_t min_ofs = -int32_t(constraints.max_racking);
		int32_t max_ofs = int32_t(constraints.max_racking);

		if (state.l == state.r) {
			//no bridges to worry about!
		} else {
			//can't have | ofs + top[l].needle - state.l_prev_needle | > top[l].left_slack
			//want -top[l].left_slack <= ofs + top[l].needle - state.l_prev_needle <= top[l].left_slack
			// -top[l].left_slack - (top[l].needle - state.l_prev_needle) <= ofs <= top[l].left_slack - (top[l].needle - state.l_prev_needle)

			int32_t l_needle = 0;
			Slack l_slack = SlackForNoYarn;
			if (state.l >= 0) {
				assert(state.l >= 0 && state.l < int32_t(bottom.size()));
				l_needle = bottom[state.l].needle;
				l_slack = bottom[state.l].right_slack;
			} else if (top_l < int32_t(top.size())) {
				assert(top_l >= 0 && top_l < int32_t(top.size()));
				l_needle = top[top_l].needle;
				l_slack = top[top_l].left_slack;
			} else if (top_l == int32_t(top.size())) {
				assert(!bottom.empty());
				l_needle = bottom.back().needle;
				l_slack = bottom.back().right_slack;
			}

			int32_t r_needle = 0;
			Slack r_slack = SlackForNoYarn;
			if (state.r < int32_t(bottom.size())) {
				assert(state.r >= 0 && state.r < int32_t(bottom.size()));
				r_needle = bottom[state.r].needle;
				r_slack = bottom[state.r].left_slack;
			} else if (top_r >= 0) {
				assert(top_r >= 0 && top_r < int32_t(top.size()));
				r_needle = top[top_r].needle;
				r_slack = top[top_r].right_slack;
			} else if (top_r == -1) {
				assert(!bottom.empty());
				r_needle = bottom[0].needle;
				r_slack = bottom[0].left_slack;
			}

			if (l_slack != SlackForNoYarn) {
				min_ofs = std::max(min_ofs, -l_slack - (l_needle - state.l_next_needle));
				max_ofs = std::min(max_ofs,  l_slack - (l_needle - state.l_next_needle));
			}
			if (r_slack != SlackForNoYarn) {
				min_ofs = std::max(min_ofs, -r_slack - (r_needle - state.r_prev_needle));
				max_ofs = std::min(max_ofs,  r_slack - (r_needle - state.r_prev_needle));
			}
		}

		//if it is possible to finish, try the finishing move:
		if (state.l < 0 && state.r >= int32_t(bottom.size())) {
			assert(top_l >= 0 && top_r < int32_t(top.size()));
			//only allow finish if zero offset is valid:
			if (min_ofs <= 0 && 0 <= max_ofs) {
				apply_action(Action(Action::Finish, 0), state, cost);
			}
		}

		{ //left moves:
			int32_t min, max;
			if (state.l >= 0) {
				assert(state.l < int32_t(bottom.size()));
				//where can bottom stitch be dropped?

				min = bottom[state.l].needle + min_ofs;
				max = bottom[state.l].needle + max_ofs;
				if (state.l == state.r) {
					//first move, no limit.
				} else {
					//NOTE: mis-use of 'can_stack_' in asymmetric case; probably will use stacking priorities for this anyway,however.
					//must place to the left of last placed stitch:
					max = std::min(max, state.l_next_needle + (bottom[state.l].can_stack_right ? 0 : -1));
				}
			} else if (top_l >= int32_t(top.size())) {
				//if we've somehow rolled *all* of the top stitches, not much to be done:
				max = std::numeric_limits< int32_t >::min();
				min = std::numeric_limits< int32_t >::max();
			} else {
				assert(top_l >= 0 && top_l < int32_t(top.size()));

				min = top[top_l].needle + min_ofs;
				max = top[top_l].needle + max_ofs;

				//must place left of last-placed stitch:
				//NOTE: mis-use of 'can_stack_' in asymmetric case; probably will use stacking priorities for this anyway,however.
				max = std::min(max, state.l_next_needle + (top[top_l].can_stack_left ? 0 : -1));

				if (state.r < int32_t(bottom.size())) {
					//can't move at all if pinned by bottom stitch:
					bool pinned = false;
					//    l --- o
					//  r ...
					if (bottom[state.r].needle <= top[top_l].needle) {
						pinned = true;
					}

					//  l --- o
					//     r ...
					if (top_l + 1 < int32_t(top.size()) && top[top_l+1].needle > bottom[state.r].needle) {
						pinned = true;
					}

					//  l ------.
					//     r -- o
					if (top_l + 1 == int32_t(top.size()) && state.r + 1 < int32_t(bottom.size())) {
						pinned = true;
					}

					if (pinned) {
						max = std::numeric_limits< int32_t >::min();
						min = std::numeric_limits< int32_t >::max();
					}
				}
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				apply_action(Action(Action::MoveLeft, needle), state, cost);
			}

		}

		{ //right moves:
			int32_t min, max;
			if (state.r < int32_t(bottom.size())) {
				assert(state.r >= 0);
				//where can bottom stitch be dropped?

				min = bottom[state.r].needle + min_ofs;
				max = bottom[state.r].needle + max_ofs;
				if (state.l == state.r) {
					//first move, no limit.
				} else {
					//NOTE: mis-use of 'can_stack_' in asymmetric case; probably will use stacking priorities for this anyway,however.
					//must place to the right of last placed stitch:
					min = std::max(min, state.r_prev_needle + (bottom[state.r].can_stack_left ? 0 : +1));
				}
			} else if (top_r < 0) {
				//if we've somehow rolled *all* of the top stitches, not much to be done
				max = std::numeric_limits< int32_t >::min();
				min = std::numeric_limits< int32_t >::max();
			} else {
				assert(top_r >= 0 && top_r < int32_t(top.size()));

				min = top[top_r].needle + min_ofs;
				max = top[top_r].needle + max_ofs;

				//must place right of last-placed stitch:
				//NOTE: mis-use of 'can_stack_' in asymmetric case; probably will use stacking priorities for this anyway,however.
				min = std::max(min, state.r_prev_needle + (top[top_r].can_stack_right ? 0 : +1));

				if (state.l >= 0) {
					//can't move at all if pinned by bottom stitch:
					bool pinned = false;
					//  o --- r
					//     .... l
					if (bottom[state.l].needle >= top[top_r].needle) {
						pinned = true;
					}

					//  o --- r
					// ... l
					if (top_r - 1 >= 0 && top[top_r-1].needle < bottom[state.l].needle) {
						pinned = true;
					}
					// /--- r
					// o ... l
					if (top_r == 0 && state.l > 0) {
						pinned = true;
					}

					if (pinned) {
						max = std::numeric_limits< int32_t >::min();
						min = std::numeric_limits< int32_t >::max();
					}
				}
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				apply_action(Action(Action::MoveRight, needle), state, cost);
			}

		}

	};


	//Initial states:
	for (int32_t i = 0; uint32_t(i) < bottom.size(); ++i) {
		State init;
		init.l = i;
		init.l_next_needle = 0;
		init.r = i;
		init.r_prev_needle = 0;
		Cost cost;
		cost.penalty = 0;
		queue_state(init, cost, nullptr, Action(Action::None, 0));
	}
	//TODO: consider states that start on *top* needles (seems like a rare and unneeded case)

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
		//std::cout << "Considering " << state->to_string() << "  [penalty: " << cost.penalty << "]" << std::endl; //DEBUG

		//if this is an ending state, end:
		if (state->l < 0 && state->r >= int32_t(bottom.size()) && state->r - state->l > int32_t(top.size() + bottom.size())) {
			best = state;
			break;
		}
		//otherwise, expand:
		expand_state(*state, cost);
	}
	assert(best && "Must have gotten to some ending state.");

	//read back operations from best:
	std::vector< Transfer > ops;

	bool is_first = true;
	while (best) {
		auto f = best_source.find(*best);
		assert(f != best_source.end());
		if (f->second.source == nullptr) break;
		State const &state = *f->second.source;
		Action const &action = f->second.action;

		int32_t top_l = -(state.l + 1);
		int32_t top_r = int32_t(top.size()) - 1 - (state.r - int32_t(bottom.size()));

		assert(state.l <= state.r);
		assert(top_l <= top_r);

		if (action.type == Action::MoveLeft) {
			if (state.l >= 0) {
				ops.emplace_back(BedNeedle(bottom_bed, bottom[state.l].needle), BedNeedle(to_bottom_bed, action.needle));
			} else { assert(state.l < 0);
				assert(top_l >= 0 && uint32_t(top_l) < top.size());
				ops.emplace_back(BedNeedle(top_bed, top[top_l].needle), BedNeedle(to_bottom_bed, action.needle));
			}
		} else if (action.type == Action::MoveRight) {
			if (state.r < int32_t(bottom.size())) {
				ops.emplace_back(BedNeedle(bottom_bed, bottom[state.r].needle), BedNeedle(to_bottom_bed, action.needle));
			} else { assert(state.r >= int32_t(bottom.size()));
				assert(top_r >= 0 && uint32_t(top_r) < top.size());
				ops.emplace_back(BedNeedle(top_bed, top[top_r].needle), BedNeedle(to_bottom_bed, action.needle));
			}
		} else if (action.type == Action::Finish) {
			assert(is_first);
		} else {
			assert(0 && "Invalid action type.");
		}
		if (action.type != Action::Finish) {
			assert(!ops.empty() && ops.back().why == "");
			ops.back().why = state.to_string() + "; " + action.to_string();
		}

		best = f->second.source;
		is_first = false;
	}
	std::cout.flush(); //DEBUG

	std::reverse(ops.begin(), ops.end());

/*
	std::cout << "  Final plan:\n"; //DEBUG
	for (auto const &op : ops) {
		std::cout << "    " << op.to_string() << '\n';
	}
	std::cout.flush(); //DEBUG
*/

/*
	std::cout << "Before Expand:\n"; //DEBUG
	draw_beds(top_bed, top, bottom_bed, bottom); //DEBUG
*/
	run_transfers(constraints,
		top_bed, top,
		bottom_bed, bottom,
		ops,
		to_top_bed, &to_top,
		to_bottom_bed, &to_bottom);
/*
	std::cout << "After Expand:\n"; //DEBUG
	draw_beds(to_top_bed, to_top, to_bottom_bed, to_bottom); //DEBUG
*/

	plan.insert(plan.end(), ops.begin(), ops.end());


}

