#include "plan_transfers-helpers.hpp"

#include <iostream>
#include <sstream>

void draw_beds(
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const &top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const &bottom
) {
	//Should make something like:
	//         -10       -9        -8          -7     <-- needle #
	// bs:    [9+2>-2-[9+2]-------1------[            <-- blob per nrg, line per edge w/ slack
	// f :     ...


	int32_t min_needle = std::numeric_limits< int32_t >::max();
	int32_t max_needle = std::numeric_limits< int32_t >::min();
	for (auto const &nrg : top) {
		min_needle = std::min(min_needle, nrg.needle);
		max_needle = std::max(max_needle, nrg.needle);
	}
	for (auto const &nrg : bottom) {
		min_needle = std::min(min_needle, nrg.needle);
		max_needle = std::max(max_needle, nrg.needle);
	}

	auto make_label = [](NeedleRollGoal const &nrg) -> std::string {
		std::string ret = "";
		if (nrg.left_slack != SlackForNoYarn) ret += std::to_string(nrg.left_slack);
		ret += (nrg.can_stack_left ? "<" : "[");
		ret += std::to_string(nrg.goal);
		if (nrg.roll > 0) ret += "+" + std::to_string(nrg.roll);
		if (nrg.roll < 0) ret += std::to_string(nrg.roll);
		ret += (nrg.can_stack_right ? ">" : "]");
		if (nrg.right_slack != SlackForNoYarn) ret += std::to_string(nrg.right_slack);
		return ret;
	};

	std::vector< std::string > needle_labels(max_needle - min_needle + 1, "");
	std::vector< std::string > top_labels(max_needle - min_needle + 1, "");
	std::vector< std::string > bottom_labels(max_needle - min_needle + 1, "");

	for (int32_t n = min_needle; n <= max_needle; ++n) {
		needle_labels[n - min_needle] = std::to_string(n);
	}

	for (auto const &s : top) {
		assert(top_labels[s.needle - min_needle] == "");
		top_labels[s.needle - min_needle] = make_label(s);
	}
	for (auto const &s : bottom) {
		assert(bottom_labels[s.needle - min_needle] == "");
		bottom_labels[s.needle - min_needle] = make_label(s);
	}

	int32_t column_width = 1;
	for (auto const &l : needle_labels) column_width = std::max< int32_t >(column_width, l.size());
	for (auto const &l : top_labels) column_width = std::max< int32_t >(column_width, l.size());
	for (auto const &l : bottom_labels) column_width = std::max< int32_t >(column_width, l.size());
	column_width += 1;

	std::ostringstream needle_str, top_str, bottom_str;

	{ //row labels:
		auto make_row_header = [](BedNeedle::Bed bed) {
			if      (bed == BedNeedle::Front) return "f :";
			else if (bed == BedNeedle::FrontSliders) return "fs:";
			else if (bed == BedNeedle::BackSliders) return "bs:";
			else if (bed == BedNeedle::Back) return "b :";
			else {
				assert(0 && "Doesn't work.");
				return "!"; //never executed
			}
		};
		needle_str << "   ";
		top_str << make_row_header(top_bed);
		bottom_str << make_row_header(bottom_bed);
	}

	auto pad = [](std::string str, char p, uint32_t width) -> std::string {
		while (str.size() < width) {
			str += p;
			if (str.size() < width) str = p + str;
		}
		return str;
	};

	for (auto &l : needle_labels) l = pad(l, ' ', column_width);
	for (auto &l : top_labels) l = pad(l, ' ', column_width);
	for (auto &l : bottom_labels) l = pad(l, ' ', column_width);

/*
	//Someday, nice edge labels
	for (uint32_t i = 0; i + 1 < top.size(); ++i) {
		if (top[i].needle == top[i+1].needle) continue;
		std::string label = std::to_string(top[i].right_slack);
		uint32_t space = 0;
		for (uint32_t n = top[i].needle; n <= top[i+1].needle; ++n) {
			if (n == top[i].needle) {
			} else if (n == top[i+1].needle) {
			} else {
				assert(top_labels[i]
			}
		}
		assert(top_labels[top[i].needle - min_needle] != "");
		uint32_t space = 
		if (top[i].
	}
*/

	for (int32_t n = min_needle; n <= max_needle; ++n) {
		needle_str << needle_labels[n - min_needle];
		top_str << top_labels[n - min_needle];
		bottom_str << bottom_labels[n - min_needle];
	}

	std::cout << needle_str.str() << '\n';
	std::cout << top_str.str() << '\n';
	std::cout << bottom_str.str() << '\n';
	std::cout.flush();
}
