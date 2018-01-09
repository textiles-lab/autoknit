#include "plan_transfers-helpers.hpp"

//!!NOTE: transfers don't know if they are rolling or not. This might be an issue!

void run_transfers(
	Constraints const &constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const &top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const &bottom,
	std::vector< Transfer > const &plan,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal > *to_top_,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal > *to_bottom_
) {
	//make sure all output arrays exist:
	assert(to_top_);
	auto &to_top = *to_top_;
	assert(to_bottom_);
	auto &to_bottom = *to_bottom_;

	//clear output arrays:
	to_top.clear();
	to_bottom.clear();

	//no stitches -> nothing to do!
	if (top.size() + bottom.size() == 0) {
		assert(plan.empty());
		return;
	}


	//transform stitches into ccw list:
	//NOTE: this is ccw when viewed with 'top' bed above 'bottom' bed; this can actually be cw in practice, but what matters is that this computation is consistent with how stitches are extracted later.
	std::vector< std::pair< BedNeedle::Bed, NeedleRollGoal > > ccw;
	ccw.reserve(bottom.size() + top.size());
	for (auto si = bottom.begin(); si != bottom.end(); ++si) {
		ccw.emplace_back(std::make_pair(bottom_bed, *si));
	}
	for (auto si = top.rbegin(); si != top.rend(); ++si) {
		ccw.emplace_back(std::make_pair(top_bed, *si));
	}
	assert(ccw.size() == bottom.size() + top.size());

	//DEBUG:
	std::cout << "before:";
	for (auto const &bs : ccw) {
		std::cout << ' ' << BedNeedle(bs.first, bs.second.needle).to_string();
	}
	std::cout << "\n";

	//everything must be on the 'from' beds:
	for (auto const &bs : ccw) {
		assert(bs.first == bottom_bed || bs.first == top_bed);
	}

	//run transfers on the ccw list:
	for (auto const &t : plan) {
		assert(t.from.bed == top_bed || t.from.bed == bottom_bed);
		assert(t.to.bed == to_top_bed || t.to.bed == to_bottom_bed);

		assert(constraints.min_free <= t.to.needle && t.to.needle <= constraints.max_free); //make sure the transfer goes to a valid needle

		//POTENTIAL OPTIMIZATION: use a hash table to avoid visiting every stitch here:
		uint32_t source = -1U;
		uint32_t target = -1U;
		for (auto &bs : ccw) {
			if (bs.first == t.to.bed && bs.second.needle == t.to.needle) {
				assert(target == -1U);
				target = &bs - &ccw[0];
			}
			if (bs.first == t.from.bed && bs.second.needle == t.from.needle) {
				assert(source == -1U);
				source = &bs - &ccw[0];
			}
		}
		assert(source != -1U);

		auto &bs = ccw[source];
		{ //transform source:
			bs.first = t.to.bed;
			bool from_is_top = (t.from.bed == top_bed);
			bool to_is_top = (t.to.bed == to_top_bed);

			//NOTE: will recompute roll, so only parity matters:
			int32_t roll = (from_is_top == to_is_top ? 0 : 1);
			bs.second = bs.second.after_offset_and_roll(t.to.needle - t.from.needle, roll);
			assert(bs.second.needle == t.to.needle);
		}

		if (target != -1U) {
			auto &tbs = ccw[target];
			assert(tbs.first == bs.first && tbs.second.needle == bs.second.needle);

			//can only stack things with the same eventual goal:
			assert(bs.second.has_same_real_goal_as(tbs.second));

			//handle stacking (merge slack / can_stack):
			if ((source + 1) % ccw.size() == target) {
				if (bs.first == to_top_bed) {
					//target is left of source, and source is stacking atop
					assert(bs.second.can_stack_left);

					bs.second.can_stack_left = tbs.second.can_stack_left;
					bs.second.left_slack = tbs.second.left_slack;
					tbs.second.can_stack_right = tbs.second.can_stack_right;
					tbs.second.right_slack = tbs.second.right_slack;
				} else { assert(bs.first == to_bottom_bed);
					//target is right of source, and source is stacking under
					assert(tbs.second.can_stack_left);

					tbs.second.can_stack_left = bs.second.can_stack_left;
					tbs.second.left_slack = bs.second.left_slack;
					bs.second.can_stack_right = bs.second.can_stack_right;
					bs.second.right_slack = bs.second.right_slack;
				}
			} else { assert((target + 1) % ccw.size() == source);
				if (bs.first == to_top_bed) {
					//target is right of source, and source is stacking atop
					assert(bs.second.can_stack_right);

					tbs.second.can_stack_left = bs.second.can_stack_left;
					tbs.second.left_slack = bs.second.left_slack;
					bs.second.can_stack_right = bs.second.can_stack_right;
					bs.second.right_slack = bs.second.right_slack;

				} else { assert(bs.first == to_bottom_bed);
					//target is left of source, and source is stacking under
					assert(tbs.second.can_stack_right);

					bs.second.can_stack_left = tbs.second.can_stack_left;
					bs.second.left_slack = tbs.second.left_slack;
					tbs.second.can_stack_right = tbs.second.can_stack_right;
					tbs.second.right_slack = tbs.second.right_slack;
				}
			}

			//TODO: eventually track how many stitches are on this location with 'weight':
			//bs.weight = tbs.weight = bs.weight + tbs.weight;

			//delete source or target stitch:
			ccw.erase(ccw.begin() + std::max(source, target));
		}
	}

	//DEBUG:
	std::cout << "after:";
	for (auto const &bs : ccw) {
		std::cout << ' ' << BedNeedle(bs.first, bs.second.needle).to_string();
	}
	std::cout << "\n";

	//must have placed everything on the 'to' beds:
	for (auto const &bs : ccw) {
		assert(bs.first == to_bottom_bed || bs.first == to_top_bed);
	}

	{ //update roll:
		auto swaps = [&to_bottom_bed,&to_top_bed](BedNeedle const &p, BedNeedle const &n) -> int8_t {
			if (p.bed == n.bed) {
				if (p.bed == to_bottom_bed) {
					return (p.needle <= n.needle ? 0 : 2);
				} else { assert(p.bed == to_top_bed);
					return (p.needle >= n.needle ? 0 : 2);
				}
			} else {
				return 1;
			}
		};

		auto from_bn = [&to_bottom_bed,&to_top_bed](std::pair< BedNeedle::Bed, NeedleRollGoal > const &bs) -> BedNeedle {
			return BedNeedle(bs.first, bs.second.needle);
		};
		auto to_bn = [&to_bottom_bed,&to_top_bed](std::pair< BedNeedle::Bed, NeedleRollGoal > const &bs) -> BedNeedle {
			if ((bs.first == to_bottom_bed) == (bs.second.roll % 2 == 0)) {
				//on bottom bed and not rolling, or on top bed and rolling
				return BedNeedle(to_bottom_bed, bs.second.goal);
			} else {
				//on top bed and not rolling, or on bottom bed and rolling
				return BedNeedle(to_top_bed, bs.second.goal);
			}
		};

		std::vector< int32_t > winding;
		winding.reserve(ccw.size());
		winding.emplace_back(from_bn(ccw[0]).bed == to_bn(ccw[0]).bed ? 0 : 1);
		for (uint32_t i = 1; i < ccw.size(); ++i) {
			winding.emplace_back(
				winding.back()
				- swaps(from_bn(ccw[i-1]), from_bn(ccw[i]))
				+ swaps(to_bn(ccw[i-1]), to_bn(ccw[i]))
			);
		}
		//make sure winding is consistent:
		assert((winding.back()
			- swaps(from_bn(ccw.back()), from_bn(ccw[0]))
			+ swaps(to_bn(ccw.back()), to_bn(ccw[0]))
			- winding[0]) % 2 == 0);

		auto minimize_winding = [](std::vector< int32_t > &winding) {
			if (winding.empty()) return;
			//pull out closest multiple of two to the median:
			std::vector< int32_t > temp = winding;
			std::sort(temp.begin(), temp.end());
			int32_t twice_median = temp[temp.size()/2] + temp[(temp.size()+1)/2];
			int32_t close_multiple = (twice_median / 4) * 2;
			int32_t DEBUG_before = 0;
			int32_t DEBUG_after = 0;
			int32_t DEBUG_after_minus = 0;
			int32_t DEBUG_after_plus = 0;
			for (auto &w : winding) {
				DEBUG_before += std::abs(w);
				w -= close_multiple;
				DEBUG_after += std::abs(w);
				DEBUG_after_minus += std::abs(w - 2);
				DEBUG_after_plus += std::abs(w + 2);
			}
			assert(DEBUG_after <= DEBUG_before);
			assert(DEBUG_after <= DEBUG_after_minus);
			assert(DEBUG_after <= DEBUG_after_plus);
		};

		minimize_winding(winding);

		assert(winding.size() == ccw.size());
		for (uint32_t i = 0; i < ccw.size(); ++i) {
			//make sure winding minimization hasn't changed meaning of roll:
			int32_t new_roll = (ccw[i].first == to_bottom_bed ? winding[i] : -winding[i]);
			assert((ccw[i].second.roll  - new_roll) % 2 == 0);
			ccw[i].second.roll = new_roll;
		}
	}


	//transform ccw list back into beds:
	{
		uint32_t first_stitch = 0; //bottom-most, ccw-most stitch
		for (uint32_t i = 1; i < ccw.size(); ++i) {
			if (ccw[i].first == to_bottom_bed) {
				if (ccw[first_stitch].first == to_top_bed || ccw[first_stitch].second.needle > ccw[i].second.needle) {
					first_stitch = i;
				}
			} else { assert(ccw[i].first == to_top_bed);
				if (ccw[first_stitch].first == to_top_bed && ccw[first_stitch].second.needle < ccw[i].second.needle) {
					first_stitch = i;
				}
			}
		}
		std::rotate(ccw.begin(), ccw.begin() + first_stitch, ccw.end());

		//ccw now looks like either bbbbbbtttttt [or just ttttt]
		auto ci = ccw.begin();

		to_bottom.clear();
		while (ci != ccw.end() && ci->first == to_bottom_bed) {
			to_bottom.emplace_back(ci->second);
			++ci;
		}

		to_top.clear();
		while (ci != ccw.end() && ci->first == to_top_bed) {
			to_top.emplace_back(ci->second);
			++ci;
		}

		assert(ci == ccw.end());

		std::reverse(to_top.begin(), to_top.end());
	}

}
