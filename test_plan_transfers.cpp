#include "plan_transfers.hpp"

#include <cassert>
#include <list>
#include <vector>
#include <random>
#include <algorithm>
#include <iostream>

//simulate_transfers applies the transfers in 'transfers' to the stitches in 'from_ccw'.
// it will return 'true' if the transfers succeed in moving the stitches to 'to_ccw' without
// violating the constraints in 'constraints' or the given 'slack'.
//
bool simulate_transfers(
	Constraints const &constraints,
	std::vector< BedNeedle > const &from_ccw,
	std::vector< BedNeedle > const &to_ccw,
	std::vector< Slack > const &slack,
	std::vector< Transfer > const &transfers,
	std::string *error = nullptr) {

	//somewhat hack-y error reporting macro; probably 'throw' would be more elegant:
	#define ERROR_UNLESS( COND, MESSAGE ) \
		do { \
			if (!( COND )) { \
				if (error) *error = (MESSAGE); \
				return false; \
			} \
		} while (0)


	// ---------------------------------
	//Check problem setup: (valid problem assumed, so these are asserts)

	assert(from_ccw.size() == to_ccw.size());
	assert(from_ccw.size() == slack.size());
	auto check_stitches = [&constraints,&slack](std::vector< BedNeedle > const &ccw) {
		if (ccw.empty()) return; //empty stitches list trivially passes all tests

		//stitches sit in the free needle range (and not on the sliders):
		for (auto s : ccw) {
			assert(constraints.min_free <= s.needle && s.needle <= constraints.max_free);
			assert(s.bed == BedNeedle::Back || s.bed == BedNeedle::Front);
		}
		//stitches obey slack:
		for (uint32_t i = 0; i < ccw.size(); ++i) {
			uint32_t n = (i + 1 == ccw.size() ? 0 : i + 1);
			int32_t diff = std::abs(int32_t(ccw[i].needle) - int32_t(ccw[n].needle));
			assert(diff <= slack[i]);
		}
		//stitches are oriented ccw:

		//count the various sorts of edges:
		uint32_t back_left = 0;
		uint32_t back_right = 0;
		uint32_t back_front = 0;
		uint32_t front_left = 0;
		uint32_t front_right = 0;
		uint32_t front_back = 0;

		auto is_front = [&ccw](uint32_t i) {
			if (ccw[i].bed == BedNeedle::Front) {
				return true;
			} else { assert(ccw[i].bed == BedNeedle::Back);
				return false;
			}
		};

		for (uint32_t i = 0; i < ccw.size(); ++i) {
			uint32_t n = (i + 1 == ccw.size() ? 0 : i + 1);
			if (is_front(i)) {
				if (is_front(n)) {
					if (ccw[n].needle < ccw[i].needle) {
						front_left += 1;
					} else if (ccw[i].needle < ccw[n].needle) {
						front_right += 1;
					}
				} else {
					front_back += 1;
				}
			} else {
				if (is_front(n)) {
					back_front += 1;
				} else {
					if (ccw[n].needle < ccw[i].needle) {
						back_left += 1;
					} else if (ccw[i].needle < ccw[n].needle) {
						back_right += 1;
					}
				}
			}
		}

		if (front_back == 0 && back_front == 0) { //on one bed
			//exactly one "seam" edge:
			assert(front_left + back_right == 1);
		} else { //on two beds
			//must cross only once:
			assert(front_back == 1 && back_front == 1);
			//no "seam" edge:
			assert(front_left == 0 && back_right == 0);
		}
	};
	check_stitches(from_ccw);
	check_stitches(to_ccw);

	// ---------------------------------
	//Check that transfers are valid:
	// (valid transfers not assumed, so these are ERROR_UNLESS)
	for (auto const &t : transfers) {
		ERROR_UNLESS(
			(t.from.bed == BedNeedle::Back && t.to.bed == BedNeedle::Front)
			|| (t.from.bed == BedNeedle::Back && t.to.bed == BedNeedle::FrontSliders)
			|| (t.from.bed == BedNeedle::Front && t.to.bed == BedNeedle::Back)
			|| (t.from.bed == BedNeedle::Front && t.to.bed == BedNeedle::BackSliders),
			"Transfer does not have valid source and destination beds."
		);
		ERROR_UNLESS(constraints.min_free <= t.from.needle && t.from.needle <= constraints.max_free, "Transfer originates outside of free needle range.");
		ERROR_UNLESS(constraints.min_free <= t.to.needle && t.to.needle <= constraints.max_free, "Transfer ends outside of free needle range.");
		ERROR_UNLESS(uint32_t(std::abs(t.to.needle - t.from.needle)) <= constraints.max_racking, "Transfer does not obey racking limit.");
	}

	// ---------------------------------
	//Check that transfers list solves the problem:

	//stitches track where each stitch is currently located:
	struct Stitch : public BedNeedle {
		Stitch(BedNeedle const &bn) : BedNeedle(bn) { }
	};

	std::vector< Stitch > stitches;
	stitches.reserve(from_ccw.size());
	for (auto f : from_ccw) {
		stitches.emplace_back(f);
	}

	//yarns track yarn constraints between stitches:
	struct Yarn {
		Yarn(Stitch *from_, Stitch *to_, int32_t slack_) : from(from_), to(to_), slack(slack_) {
			assert(from);
			assert(to);
		}
		Stitch *from;
		Stitch *to;
		int32_t slack;
	};

	std::vector< Yarn > yarns;
	yarns.reserve(slack.size());
	for (uint32_t i = 0; i < slack.size(); ++i) {
		if (slack[i] == SlackForNoYarn) continue; //don't create yarn when no yarn exists
		uint32_t n = (i + 1 == slack.size() ? 0 : i + 1);
		yarns.emplace_back( &stitches[i], &stitches[n], slack[i] );
	}

	//needles provide a spatial look-up table for stitches and pins:
	struct Needle {
		std::vector< Stitch * > back;
		std::vector< Stitch * > back_sliders;
		std::vector< Stitch * > front_sliders;
		std::vector< Stitch * > front;
		std::vector< Stitch * > &on(BedNeedle::Bed bed) {
			if (bed == BedNeedle::Back) return back;
			else if (bed == BedNeedle::BackSliders) return back_sliders;
			else if (bed == BedNeedle::FrontSliders) return front_sliders;
			else if (bed == BedNeedle::Front) return front;
			else assert(0);
		}
	};

	//figure out the min and max needle the code needs to represent:
	int32_t min_needle = std::numeric_limits< int32_t >::max();
	int32_t max_needle = std::numeric_limits< int32_t >::min();
	for (auto const &s : from_ccw) {
		min_needle = std::min(min_needle, s.needle);
		max_needle = std::max(max_needle, s.needle);
	}
	for (auto const &s : to_ccw) {
		min_needle = std::min(min_needle, s.needle);
		max_needle = std::max(max_needle, s.needle);
	}
	for (auto const &t : transfers) {
		min_needle = std::min(min_needle, t.from.needle);
		min_needle = std::min(min_needle, t.to.needle);
		max_needle = std::max(max_needle, t.from.needle);
		max_needle = std::max(max_needle, t.to.needle);
	}

	assert(min_needle <= max_needle);

	//needles[0] corresponds to min_needle:
	std::vector< Needle > needles(max_needle - min_needle + 1);

	//helpers to avoid that index math:
	auto stitches_on = [&needles,&min_needle,&max_needle](BedNeedle const &bn) -> std::vector< Stitch * > & {
		assert(bn.needle >= min_needle && bn.needle <= max_needle);
		return needles[bn.needle - min_needle].on(bn.bed);
	};

	//walk through transfers one by one and check state:
	for (auto const &t : transfers) {
		//find all stitches that are in the 'from' location:
		std::vector< Stitch * > &from = stitches_on(t.from);

		ERROR_UNLESS(!from.empty(), "Transferring from an empty needle."); //NOTE: could be a WARNING because it's not wrong just inefficient


		//make sure that transfer isn't going through any stitches on sliders:
		if (t.from.bed == BedNeedle::Back) {
			ERROR_UNLESS(stitches_on(BedNeedle(BedNeedle::BackSliders, t.from.needle)).empty(), "Transfer from hook must have empty slider.");
		} else if (t.from.bed == BedNeedle::Front) {
			ERROR_UNLESS(stitches_on(BedNeedle(BedNeedle::FrontSliders, t.from.needle)).empty(), "Transfer from hook must have empty slider.");
		}
		if (t.to.bed == BedNeedle::Back) {
			ERROR_UNLESS(stitches_on(BedNeedle(BedNeedle::BackSliders, t.to.needle)).empty(), "Transfer to hook must have empty slider.");
		} else if (t.to.bed == BedNeedle::Front) {
			ERROR_UNLESS(stitches_on(BedNeedle(BedNeedle::FrontSliders, t.to.needle)).empty(), "Transfer to hook must have empty slider.");
		}

		//check that racking is valid (assuming that no yarn is pinned)
		int32_t racking; //measured as back - front
		if (t.from.bed == BedNeedle::Back || t.from.bed == BedNeedle::BackSliders) {
			racking = t.from.needle - t.to.needle;
		} else {
			racking = t.to.needle - t.from.needle;
		}

		for (auto const &y : yarns) {
			bool from_is_front = (y.from->bed == BedNeedle::Front || y.from->bed == BedNeedle::FrontSliders);
			bool to_is_front = (y.to->bed == BedNeedle::Front || y.to->bed == BedNeedle::FrontSliders);
			if (from_is_front && !to_is_front) {
				int32_t dis = std::abs((y.to->needle + racking) - y.from->needle);
				ERROR_UNLESS(dis < y.slack, "Transfer requires racking that stretches yarn too much.");
			}
			if (!from_is_front && to_is_front) {
				int32_t dis = std::abs(y.to->needle - (y.from->needle + racking));
				ERROR_UNLESS(dis < y.slack, "Transfer requires racking that stretches yarn too much.");
			}
		}

		//move stitches from from to to:
		std::vector< Stitch * > &to = stitches_on(t.to);

		for (auto s = from.rbegin(); s != from.rend(); ++s) {
			(*s)->bed = t.to.bed;
			(*s)->needle = t.to.needle;
			to.emplace_back(*s);
		}
		from.clear();
	}

	//check that all stitches arrived at their targets:
	for (uint32_t i = 0; i < to_ccw.size(); ++i) {
		ERROR_UNLESS(stitches[i].bed == to_ccw[i].bed && stitches[i].needle == to_ccw[i].needle, "Stitches must actually arrive at their destinations.");
	}

	return true;

}

void dump_layout(std::vector< BedNeedle > const &ccw, std::vector< Slack > *slack = nullptr) {
	if (slack) {
		assert(slack->size() == ccw.size());
	}

	for (uint32_t i = 0; i < ccw.size(); ++i) {
		if (i != 0) std::cout << ' ';
		if (ccw[i].bed == BedNeedle::Front) std::cout << 'f';
		else if (ccw[i].bed == BedNeedle::Back) std::cout << 'b';
		else std::cout << '?';
		std::cout << ccw[i].needle;
		if (slack) {
			std::cout << ' ';
			if ((*slack)[i] == SlackForNoYarn) std::cout << '*';
			else std::cout << (*slack)[i];
		}
	}
	std::cout << '\n';

	int32_t min_needle = std::numeric_limits< int32_t >::max();
	int32_t max_needle = std::numeric_limits< int32_t >::min();
	for (auto const &bn : ccw) {
		min_needle = std::min(min_needle, bn.needle);
		max_needle = std::max(max_needle, bn.needle);
	}
	uint32_t needle_size = 1;
	uint32_t edge_size = 0;
	if (slack) {
		for (auto s : *slack) {
			if (s == SlackForNoYarn) continue;
			edge_size = std::max< uint32_t >(edge_size, std::to_string(s).size());
		}
	}
	edge_size += 2;

	std::string back = std::string(edge_size + (max_needle - min_needle + 1) * (needle_size + edge_size), ' ');
	std::string front = std::string(edge_size + (max_needle - min_needle + 1) * (needle_size + edge_size), ' ');

	for (auto const &bn : ccw) {
		std::string label = "o";
		if (&bn == &ccw[0]) label = "*";
		assert(label.size() == needle_size);
		if (bn.bed == BedNeedle::Front) {
			for (uint32_t i = 0; i < label.size(); ++i) {
				uint32_t bi = edge_size + (bn.needle - min_needle) * (needle_size + edge_size) + i;
				assert(bi < front.size());
				front[bi] = label[i];
			}
		} else { assert(bn.bed == BedNeedle::Back);
			for (uint32_t i = 0; i < label.size(); ++i) {
				uint32_t bi = edge_size + (bn.needle - min_needle) * (needle_size + edge_size) + i;
				assert(bi < back.size());
				back[bi] = label[i];
			}
		}
	}

	if (slack) {
		for (uint32_t s = 0; s < ccw.size(); ++s) {
			if ((*slack)[s] == SlackForNoYarn) continue;
			BedNeedle const &bn = ccw[s];
			BedNeedle const &next = ccw[s+1 < ccw.size() ? s + 1 : 0];
			if (bn == next) continue;
			if (bn.bed == BedNeedle::Front && next.bed == BedNeedle::Front) {
				if (bn.needle < next.needle) {
					uint32_t begin = edge_size + (bn.needle - min_needle) * (needle_size + edge_size) + needle_size;
					uint32_t end = edge_size + (next.needle - min_needle) * (needle_size + edge_size);
					std::string label = std::to_string((*slack)[s]);
					while (label.size() < end - begin) {
						label += '-';
						if (label.size() < end - begin) label = '-' + label;
					}
					for (uint32_t i = begin; i < end; ++i) {
						assert(i < front.size() && front[i] == ' ');
						front[i] = label[i - begin];
					}
				} else {
					//TODO: wrap around the entire back
				}
			} else if (bn.bed == BedNeedle::Back && next.bed == BedNeedle::Back) {
				if (next.needle < bn.needle) {
					uint32_t begin = edge_size + (next.needle - min_needle) * (needle_size + edge_size) + needle_size;
					uint32_t end = edge_size + (bn.needle - min_needle) * (needle_size + edge_size);
					std::string label = std::to_string((*slack)[s]);
					while (label.size() < end - begin) {
						label += '-';
						if (label.size() < end - begin) label = '-' + label;
					}
					for (uint32_t i = begin; i < end; ++i) {
						assert(i < back.size() && back[i] == ' ');
						back[i] = label[i - begin];
					}
				} else {
					//TODO: wrap around the entire front
				}
			} else if (bn.bed == BedNeedle::Front && next.bed == BedNeedle::Back) {
				//right edge
				uint32_t front_begin = edge_size + (bn.needle - min_needle) * (needle_size + edge_size) + needle_size;
				uint32_t back_begin = edge_size + (next.needle - min_needle) * (needle_size + edge_size) + needle_size;
				uint32_t end = edge_size + (max_needle - min_needle + 1) * (needle_size + edge_size);
				std::string front_label = "";
				std::string back_label = "";
				if (back_begin <= front_begin) {
					back_label = std::to_string((*slack)[s]);
				} else {
					front_label = std::to_string((*slack)[s]);
				}
				while (front_label.size() < end - front_begin) {
					front_label += '-';
					if (front_label.size() < end - front_begin) front_label = '-' + front_label;
				}
				front_label.back() = '/';
				while (back_label.size() < end - back_begin) {
					back_label += '-';
					if (back_label.size() < end - back_begin) back_label = '-' + back_label;
				}
				back_label.back() = '\\';

				for (uint32_t i = front_begin; i < end; ++i) {
					assert(i < front.size() && front[i] == ' ');
					front[i] = front_label[i - front_begin];
				}
				for (uint32_t i = back_begin; i < end; ++i) {
					assert(i < back.size() && back[i] == ' ');
					back[i] = back_label[i - back_begin];
				}
			} else if (bn.bed == BedNeedle::Back && next.bed == BedNeedle::Front) {
				//left edge
				uint32_t begin = 0;
				uint32_t front_end = edge_size + (next.needle - min_needle) * (needle_size + edge_size);
				uint32_t back_end = edge_size + (bn.needle - min_needle) * (needle_size + edge_size);

				std::string front_label = "";
				std::string back_label = "";
				if (back_end >= front_end) {
					back_label = std::to_string((*slack)[s]);
				} else {
					front_label = std::to_string((*slack)[s]);
				}
				while (front_label.size() < front_end - begin) {
					front_label += '-';
					if (front_label.size() < front_end - begin) front_label = '-' + front_label;
				}
				front_label[0] = '\\';
				while (back_label.size() < back_end - begin) {
					back_label += '-';
					if (back_label.size() < back_end - begin) back_label = '-' + back_label;
				}
				back_label[0] = '/';

				for (uint32_t i = begin; i < front_end; ++i) {
					assert(i < front.size() && front[i] == ' ');
					front[i] = front_label[i - begin];
				}
				for (uint32_t i = begin; i < back_end; ++i) {
					assert(i < back.size() && back[i] == ' ');
					back[i] = back_label[i - begin];
				}
			} else {
				assert(0 && "I think we handle all cases");
			}
		}
	}

	std::cout << back << '\n' << front << '\n';
	std::cout.flush();
};


bool test_plan_transfers() {
	static std::mt19937 mt(0xdeadbeef);

	Constraints constraints;
	//pick max racking in range [1,20]; generally pick 4 which is a realistic value.
	if (mt() < mt.max() / 4) {
		constraints.max_racking = 1 + (mt() % 20);
	} else {
		constraints.max_racking = 4;
	}

	uint32_t count;
	{ //select count:
		uint32_t val = mt();
		//let's be overly fancy in our size distribution and pick sizes up to 10 with uniform probability:
		//if (val < mt.max() / 2) {
			count = 1 + (val % 10);
		/*} else {
			//...and larger sizes with some sort of decaying probability:
			float amt = (val - mt.max()/2) / float(mt.max()/2);
			count = 10 + std::floor(100.0f * (1.0f - std::pow(1.0f / 100.0f, amt)));
		}*/
	}


	std::vector< Slack > slack;
	{ //build random slacks:
		slack.resize(count);
		for (auto &s : slack) {
			uint32_t val = mt();
			if (val < mt.max() / 32) {
				s = val % 5 + 1;
			} else if (val < mt.max() / 8) {
				s = 3;
			} else if (val < mt.max() / 4) {
				s = 2;
			} else {
				s = 1;
			}
		}
	}

	std::vector< BedNeedle > from_ccw;
	std::vector< BedNeedle > to_ccw;
	{ //build layout:
		//decrease some slack randomly:
		auto decrease_some_slack = [](std::vector< Slack > &layout_slack) {
			for (auto &s : layout_slack) {
				uint32_t val = mt();
				if (val < mt.max() / 32) {
					s = 0;
				} else if (val < mt.max() / 8 && s >= 2) {
					s -= 2;
				} else if (val < mt.max() / 4 && s >= 1) {
					s -= 1;
				}
				assert(s >= 0);
			}
		};

		//decrease slack -- stitches can be arranged with less-than-max slack between them:
		std::vector< Slack > from_slack = slack;
		decrease_some_slack(from_slack);
		std::vector< Slack > to_slack = slack;
		decrease_some_slack(to_slack);

		//though, if stitches are stuck together in 'from' they need to be stuck together in 'to':
		for (uint32_t i = 0; i < to_slack.size(); ++i) {
			if (from_slack[i] == 0) to_slack[i] = 0;
		}

		//add a break with some probability:
		if ((mt() & 3) == 0) {
			from_slack.back() = to_slack.back() = SlackForNoYarn;
		}

		//make layouts given layout slack:
		auto make_layout = [](std::vector< Slack > const &layout_slack) {
			std::vector< BedNeedle > ccw;
			if (layout_slack.back() == SlackForNoYarn) {
				//single-bed layout:
				int32_t n = 0;
				for (auto const &s : layout_slack) {
					ccw.emplace_back(BedNeedle::Front, n);
					n += s;
				}
			} else {
				//do a two-bed layout:
				int32_t total_slack = 0;
				for (auto s : layout_slack) total_slack += s;

				int32_t front_slack = total_slack / 2;

				int32_t index = mt() % (total_slack + 1);
				for (auto s : layout_slack) {
					if (index > total_slack) index -= total_slack;

					if (index < front_slack) {
						ccw.emplace_back(BedNeedle::Front, index);
					} else {
						ccw.emplace_back(BedNeedle::Back, 2 * front_slack - index);
					}

					index += s;
				}

			}
			return ccw;
		};
		from_ccw = make_layout(from_slack);
		to_ccw = make_layout(to_slack);

		//swap front and back beds:
		if (mt() % 2) {
			for (auto &s : from_ccw) {
				s.bed = (s.bed == BedNeedle::Front ? BedNeedle::Back : BedNeedle::Front);
				s.needle = -s.needle;
			}
		}
		if (mt() % 2) {
			for (auto &s : to_ccw) {
				s.bed = (s.bed == BedNeedle::Front ? BedNeedle::Back : BedNeedle::Front);
				s.needle = -s.needle;
			}
		}

		//pick random offsets for beds:
		auto assign_random_offsets = [&constraints,&slack](std::vector< BedNeedle > &ccw) {
			assert(ccw.size() == slack.size());

			int32_t front_min = std::numeric_limits< int32_t >::max();
			int32_t front_max = std::numeric_limits< int32_t >::min();
			int32_t back_min = std::numeric_limits< int32_t >::max();
			int32_t back_max = std::numeric_limits< int32_t >::min();
			int32_t right_slack = -1; //slack between back and front stitch
			int32_t left_slack = -1;
			int32_t right_offset = 0; //needle index of back - front stitch
			int32_t left_offset = 0;
			for (auto const &bn : ccw) {
				uint32_t i = &bn - &ccw[0];
				uint32_t n = (i + 1 < ccw.size() ? i + 1 : 0); //next index
				if (bn.bed == BedNeedle::Front) {
					front_min = std::min(front_min, bn.needle);
					front_max = std::max(front_max, bn.needle);
					if (ccw[n].bed == BedNeedle::Back) {
						assert(right_slack == -1);
						right_slack = slack[i];
						right_offset = ccw[n].needle - bn.needle;
					}
				} else { assert(bn.bed == BedNeedle::Back);
					back_min = std::min(back_min, bn.needle);
					back_max = std::max(back_max, bn.needle);
					if (ccw[n].bed == BedNeedle::Front) {
						assert(left_slack == -1);
						left_slack = slack[i];
						left_offset = bn.needle - ccw[n].needle;
					}
				}
			}
			int32_t new_front_min;
			int32_t new_back_min;
			if (front_min > front_max) {
				assert(back_min <= back_max);
				assert(left_slack == -1 && right_slack == -1);
				//back bed only
				new_back_min = new_front_min = (mt() % 21) - 10;
			} else if (back_min > back_max) {
				assert(front_min <= front_max);
				assert(left_slack == -1 && right_slack == -1);
				//front bed only
				new_back_min = new_front_min = (mt() % 21) - 10;
			} else { assert(front_min <= front_max && back_min <= back_max);
				assert(left_slack > 0 && right_slack > 0);
				//both beds

				std::cout << "left ofs/slack: " << left_offset << "/" << left_slack << " right ofs/slack: " << right_offset << "/" << right_slack << std::endl; //DEBUG

				//currently, works at offset zero, right?
				assert(std::abs(left_offset) <= left_slack);
				assert(std::abs(right_offset) <= right_slack);

				//see how far we can slide left and still have this work:
				int32_t min_offset = 0;
				while ( std::abs(min_offset - 1 + left_offset) < left_slack
				     && std::abs(min_offset - 1 + right_offset) < right_slack ) {
					min_offset -= 1;
				}

				//see how far we can slide right and still have this work:
				int32_t max_offset = 0;
				while ( std::abs(max_offset + 1 + left_offset) < left_slack
				     && std::abs(max_offset + 1 + right_offset) < right_slack ) {
					max_offset += 1;
				}

				//remember that one can also force the machine to be racked at the start of the exercise:
				// (can one?)
				min_offset -= int32_t(constraints.max_racking);
				max_offset += int32_t(constraints.max_racking);

				new_front_min = (mt() % 21) - 10;
				new_back_min = new_front_min + (mt() % (max_offset - min_offset + 1) + min_offset);
			}

			for (auto &bn : ccw) {
				if (bn.bed == BedNeedle::Front) {
					bn.needle = (bn.needle - front_min) + new_front_min;
				} else { assert(bn.bed == BedNeedle::Back);
					bn.needle = (bn.needle - back_min) + new_back_min;
				}
			}
		};
		assign_random_offsets(from_ccw);
		assign_random_offsets(to_ccw);

		{ //roll stitch arrays:
			uint32_t roll = mt() % from_ccw.size();
			std::rotate(from_ccw.begin(), from_ccw.begin() + roll, from_ccw.end());
			std::rotate(to_ccw.begin(), to_ccw.begin() + roll, to_ccw.end());
		}

	}

	{ //pick free range:
		int32_t min_front_needle = std::numeric_limits< int32_t >::max();
		int32_t max_front_needle = std::numeric_limits< int32_t >::min();
		int32_t min_back_needle = std::numeric_limits< int32_t >::max();
		int32_t max_back_needle = std::numeric_limits< int32_t >::min();
		for (auto const &bn : from_ccw) {
			if (bn.bed == BedNeedle::Front) {
				min_front_needle = std::min(min_front_needle, bn.needle);
				max_front_needle = std::max(max_front_needle, bn.needle);
			} else { assert(bn.bed == BedNeedle::Back);
				min_back_needle = std::min(min_back_needle, bn.needle);
				max_back_needle = std::max(max_back_needle, bn.needle);
			}
		}
		for (auto const &bn : to_ccw) {
			if (bn.bed == BedNeedle::Front) {
				min_front_needle = std::min(min_front_needle, bn.needle);
				max_front_needle = std::max(max_front_needle, bn.needle);
			} else { assert(bn.bed == BedNeedle::Back);
				min_back_needle = std::min(min_back_needle, bn.needle);
				max_back_needle = std::max(max_back_needle, bn.needle);
			}
		}
		{
			uint32_t val = mt();
			if (val < mt.max() / 4) {
				//leave constraints min free as infinite
				constraints.min_free = std::numeric_limits< int32_t >::min();
			} else if (val < mt.max() / 2) {
				//constrain to only one needle
				constraints.min_free = std::min(min_front_needle, min_back_needle) - 1;
			} else {
				//give several needles
				constraints.min_free = std::min(min_front_needle, min_back_needle) - 1 - (mt() % 10);
			}
		}
		{
			uint32_t val = mt();
			if (val < mt.max() / 4) {
				//leave constraints max free as infinite
				constraints.max_free = std::numeric_limits< int32_t >::max();
			} else if (val < mt.max() / 2) {
				//constrain to only one needle
				constraints.max_free = std::max(max_front_needle, max_back_needle) + 1;
			} else {
				//give several needles
				constraints.max_free = std::max(max_front_needle, max_back_needle) + 1 + (mt() % 10);
			}
		}

	}


	std::cout << "From:\n";
	dump_layout(from_ccw, &slack); //DEBUG
	std::cout << "To:\n";
	dump_layout(to_ccw, &slack); //DEBUG


	{ //actually run transfer planning:
		std::vector< Transfer > transfers;
		std::string error = "";

		bool result = plan_transfers(
			constraints,
			from_ccw,
			to_ccw,
			slack,
			&transfers,
			&error);

		if (!result) {
			std::cerr << "ERROR: plan_transfers failed:\n" << error << std::endl;
			return false;
		}

		if (!simulate_transfers(
			constraints,
			from_ccw,
			to_ccw,
			slack,
			transfers,
			&error)) {
			std::cerr << "ERROR: result of plan_transfers failed checks:\n" << error << std::endl;
			return false;
		}
	}

	return true;
}

int main(int argc, char **argv) {
	for (uint32_t i = 0; i < 1000; ++i) {
		if (!test_plan_transfers()) return 1;
	}
	return 0;
}
