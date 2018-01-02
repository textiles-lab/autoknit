#include "plan_transfers.hpp"

#include <cassert>
#include <list>
#include <vector>

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
			assert(diff < slack[i]);
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

	//pins track where yarn is trapped under stitches:
	struct Pin {
		BedNeedle::Bed under = BedNeedle::Front;
		int32_t needle = 0;
	};

	std::list< Pin > pins;

	//yarns track yarn constraints between stitches:
	struct Yarn {
		Yarn(Stitch *from_, Stitch *to_, int32_t slack_) : from(from_), to(to_), slack(slack_) {
			assert(from);
			assert(to);
		}
		Stitch *from;
		std::vector< Pin * > pins; //any places the yarn is pinned between from and to
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
		std::vector< Pin * > under_back;
		std::vector< Stitch * > back;
		std::vector< Pin * > under_back_sliders;
		std::vector< Stitch * > back_sliders;
		std::vector< Stitch * > front_sliders;
		std::vector< Pin * > under_front_sliders;
		std::vector< Stitch * > front;
		std::vector< Pin * > under_front;
		std::vector< Pin * > &under(BedNeedle::Bed bed) {
			if (bed == BedNeedle::Back) return under_back;
			else if (bed == BedNeedle::BackSliders) return under_back_sliders;
			else if (bed == BedNeedle::FrontSliders) return under_front_sliders;
			else if (bed == BedNeedle::Front) return under_front;
			else assert(0);
		}
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
	auto pins_under = [&needles,&min_needle,&max_needle](BedNeedle const &bn) -> std::vector< Pin * > & {
		assert(bn.needle >= min_needle && bn.needle <= max_needle);
		return needles[bn.needle - min_needle].under(bn.bed);
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

		//figure out what yarns are getting pinned by this move:
		std::vector< Pin * > new_pins;
		for (auto &yarn : yarns) {
			for (int32_t i = -1; i < int32_t(yarn.pins.size()); ++i) {
				int32_t from;
				int32_t to;
				if (i == -1) from = yarn.from->needle;
				else from = yarn.pins[i]->needle;
				if (i + 1 < int32_t(yarn.pins.size())) to = yarn.pins[i+1]->needle;
				else to = yarn.to->needle;

				//question: does the line t.from.needle -> t.to.needle cross the line from -> to?

				//TODO: seems like one needs to think about left/right on pins/stitches to answer this.
				//<--- I WAS HERE
			}
		}

		//TODO: move stitches from from to to

	}

	//TODO: check that all stitches arrived at their targets

	return true;

}


int main(int argc, char **argv) {

	

	return 0;
}
