#include "Stitch.hpp"

#include "TaggedArguments.hpp"

#include <deque>
#include <map>
#include <set>

//Loop held on a needle:
struct Loop {
	constexpr Loop(uint32_t stitch_, uint32_t idx_) : stitch(stitch_), idx(idx_) {
	}
	uint32_t stitch;
	uint32_t idx;
	bool operator==(Loop const &o) const {
		return (stitch == o.stitch && idx == o.idx);
	}
};
constexpr const Loop GAP = Loop(-1U, -1U);


int main(int argc, char **argv) {
	std::string in_st = "";
	std::string out_js = "";
	{ //parse arguments:
		TaggedArguments args;
		args.emplace_back("st", &in_st, "input stitches file (required)");
		args.emplace_back("js", &out_js, "output knitting file");
		bool usage = !args.parse(argc, argv);
		if (!usage && in_st == "") {
			std::cerr << "ERROR: 'st:' argument is required." << std::endl;
			usage = true;
		}
		if (usage) {
			std::cerr << "Usage:\n\t./schedule [tag:value] [...]\n" << args.help_string() << std::endl;
			return 1;
		}
	}

	//------------------------------

	std::vector< Stitch > stitches;
	if (!load_stitches(in_st, &stitches)) {
		std::cerr << "ERROR: failed to load stitches from '" << in_st << "'." << std::endl;
		return 1;
	}
	std::cout << "Read " << stitches.size() << " stitches from '" << in_st << "'." << std::endl;

	//------------------------------

	//New scheduling workflow:
	// (1) split into "steps" ==> tube-supported bits of knitting that will be done at once
	//    - each step eventually needs a shape + roll + offset for its output loops
	//      (implies a shape for input loops)
	// (2) pick a consistent shape + roll for "interesting" steps
	//    - these are steps that take loops from more than one output or input
	//    - effectively, this finds an upward planar embedding, where the edges are chains of construction steps and the vertices occur when steps have more than one tube as a parent or child.
	// (2) figure out a layout (shape + roll) for each step


	//Cycle is a tube-supported chunk of loops on the bed.
	struct Cycle : public std::deque< Loop > {
		//cycle is stored in CCW order.
		//Loops created by corresponding stitch will either be at the .back() [ccw stitch] or .front() [cw stitch]
		// one or more 'GAP' stitches may exist in open cycles.

		//TODO: shape!
		//TODO: info about bridges! (== past merges)
	};

	//Cycles are 1-1 with stitches -- each stitch creates a new cycle from an old cycle

	std::vector< Cycle > cycles;

	#define REPORT_ERROR( X ) do { std::cerr << (X) << std::endl; exit(1); } while(0)

	{ //build cycles:
		//current yarn position w.r.t. loops:
		struct YarnInfo {
			Loop loop = GAP;
			char direction = '\0';
		};
		std::map< uint32_t, YarnInfo > active_yarns;

		//current version of all cycles:
		std::set< uint32_t > active_cycles;

		/*
		//helper: get the last-constructed loop from a stitch index:
		auto last_loop = [&stitches](uint32_t idx) -> Loop {
			assert(idx < stitches.size());
			auto const &s = stitches[idx];
			if (s.out[1] < stitches.size()) {
				return Loop(idx, 1);
			} else {
				assert(s.out[0] < stitches.size());
				return Loop(idx, 0);
			}
		};
		*/

		//helper: make a new active cycle by:
		// (1) finding the 'find' loops
		// (2) merging cycles so that they are in the same cycle
		// (3) replacing them with the 'replace' loops
		auto make_cycle = [&cycles, &active_cycles](std::vector< Loop > const &find, std::vector< Loop > const &replace) {
			//(Not as general a function as it appears.)

			//Find each non-GAP loop:
			struct CycleIndex {
				uint32_t cycle = -1U;
				uint32_t idx = -1U;
			};

			//build an index of active loops:
			std::map< Loop, CycleIndex > index;
			std::vector< CycleIndex > gaps;
			for (auto ac : active_cycles) {
				auto const &cycle = cycles[ac];
				for (auto const &loop : cycle) {
					CycleIndex ci(ac, &loop - &cycle[0]);
					if (loop == GAP) {
						gaps.emplace_back(ci);
					} else {
						auto res = index.insert(std::make_pair(loop, ci));
						assert(res.second);
					}
				}
			}

			//helper: find loop using index
			auto find_loop = [&index](Loop const &loop) -> CycleIndex const & {
				auto f = index.find(loop);
				assert(f != index.end() && "input should always be in an active cycle");
				return f->second;
			};

			//helper: remove cycle from active cycles and return it
			auto deactivate_cycle = [&active_cycles, &cycles](uint32_t c) -> Cycle const & {
				assert(c < cycles.size());
				//deactivate and return a cycle:
				auto f = active_cycles.find(c);
				assert(f != active_cycles.end());
				active_cycles.erase(f);
				return cycles[c];
			};

			auto roll_to_back = [](Cycle &cycle, uint32_t index) {
				assert(index < cycle.size());
				std::rotate(cycle.begin(), cycle.begin() + index + 1, cycle.end());
			};

			//always will make a new cycle:
			active_cycles.insert(cycles.size());

			//Special case the easy things:
			if (find.size() == 1 && find[0] == GAP) {
				//finding *just* a GAP always creates a new cycle.
				cycles.emplace_back();
				cycles.back().assign(replace);
			} else if (find.size() == 1 && find[0] != GAP) {
				//finding a non-gap is straightforward:
				CycleIndex ci = find_loop(find[0]);
				//deactivate old cycle; use as basis for new cycle:
				cycles.emplace_back(deactivate_cycle(ci.cyle));
				Cycle &cycle = cycles.back();
				//perform replacement:
				cycle.erase(cycle.begin() + ci.index);
				cycle.insert(cycle.begin() + ci.index, replace.begin(), replace.end());
			} else if (find.size() == 2 && find[0] != GAP) {
				CycleIndex ci0 = find_loop(find[0]);

				//start with the ci0 cycle:
				cycles.emplace_back(deactivate_cycle(ci0.cycle));
				Cycle &cycle = cycles.back();

				assert(cycle.size() >= 2); //all cycles are at least size 2

				//arrange find[0] at the back:
				roll_to_back(cycle, ci0.index);
				assert(cycle.back() == find[0]);

				if (cycle.front() == find[1]) {
					//all in one cycle; great!
					cycle.pop_back(); //remove find[0] from the back
					cycle.pop_front(); //remove find[1] from the front
					//insert replacement:
					cycle.insert(cycle.end(), replace.begin(), replace.end());
				} else {
					//not all in one cycle, so...
					if (find[1] == GAP) { //...add gap:
						cycle.add_bridge(); //always adds bridge at the back/front interface
						assert(cycle.back() == GAP);
						cycle.pop_back();
						assert(cycle.back() == find[0]);
						cycle.pop_back();
						//insert replacement:
						cycle.insert(cycle.end(), replace.begin(), replace.end());
					} else {
						CycleIndex ci1 = find_loop(find[1]);
						if (ci1.cycle == ci0.cycle) { //...split cycle:
							auto f = std::find(cycle.begin(), cycle.end(), find[1]);
							assert(f != cycle.end());
							//<--- I WAS HERE
							// NOTES:
							//   This needs to "split" a cycle.. sort of.
							//   One way to do this:
							//     track active *loops* instead of cycles,
							//     and have the split portion of the cycle co-opt those
							//     loops it contains.
							//   Maintains 1-1 property (good!)
							//   Do we need to add a bridge somewhere?
							//   Thought: maybe bridges are tracked independently of cycles?
							//   Thought: perhaps he way to go is track active constraints and implied cycles? (dump 1-1 idea)
						} else { //...merge cycles:
						assert(ci1.cycle != ci0.cycle);
						Cycle to_merge = deactivate_cycle(ci1.cycle);
						roll_to_back(to_merge, ci1.index);
						assert(to_merge.back() == find[1]);

						cycles.emplace_back(merge_cycles(
					}
				}
				
			} else {
				assert(0 && "make_cycle handles very few cases for find/replace");
			}

			{ //(1) Look for 'find' as a subset of a single cycle:
				
			for (auto const &f : find) {
				if (f === GAP)
				if (f == GAP) {
					//handle later:
					found.emplace_back();
					continue;
				}
				for (auto ac : active_cycles) {
					auto const &cycle = cycles[ac];
					for (auto const &loop : cycle) {
						if (f == loop) {
							found.
						}
					}
				}
			}
			std::vector< uint32_t > inds;


		};

		for (uint32_t si = 0; si < stitches.size(); ++si) {
			Stitch const &s = stitches[si];
			YarnInfo &yarn = active_yarns[s.yarn];

			if (s.type == Stitch::Start) {
			/* someday...
				Loop loop(si, 0);
				if (pi == -1U) {
					//new yarn -> new cycle.
					Cycle cycle;
					cycle.emplace_back(Loop::GAP);
					if (s.direction == Stitch::CCW) {
						cycle.emplace_back(loop);
					} else {
						cycle.emplace_front(loop);
					}
					active_cycles.insert(cycles.size());
					cycles.emplace_back(cycle);
				} else {
					Loop prev = last_loop(pi);

					Stitch const &p = stitches[pi];
					if (p.direction != s.directon) {
						REPORT_ERROR("Can't switch direction and cast on.");
					}

					//find the cycle with the previous stitch's loop:
					auto f = lookup_loop(prev);

					assert(pi + 1 == cycles.size());
					Cycle cycle = cycles[pi]; //copy previous cycle
					assert(cycle.size() >= 1 && "no stitch-associated cycle will be empty"); //must be -- can't create a zero-size cycle
					if (s.direction == Stitch::CCW) {
						assert(cycle.back().stitch == pi && "previous cycle should respect stitch position invariant");
						if (cycle.front() != Loop::GAP) {
							//TODO: build a bridge
							REPORT_ERROR("Can't cast on without a gap.");
						}
						cycle.emplace_back(loop);
					} else {
						assert(cycle.front().stitch == pi && "previous cycle should respect stitch position invariant");
						if (cycle.back() != Loop::GAP) {
							//TODO: build a bridge
							REPORT_ERROR("Can't cast on without a gap.");
						}
						cycle.emplace_front(loop);
					}
				}
			*/
			} else if (s.type == Stitch::Tuck || s.type == Stitch::Miss || s.type == Stitch::Knit) {
				assert(s.in[0] < si && s.in[1] == -1U && s.out[0] != -1U && s.out[0] > si && s.out[1] == -1U && "valid 1-1 stitch");
				Loop in0(s.in[0], stitches[s.in[0]].find_out(si));
				Loop out0(si, 0);

				std::vector< Loop > find, replace;
				if (yarn.loop == GAP) {
					//bringing in yarn, I suppose:
					assert("TODO: bring in yarn on 1-1 stitch(!?!)"); //kinda weird case
					find = {in0};
					replace = {out0};
				} else if (yarn.loop == in0 && yarn.direction != s.direction) {
					//turning around:
					find = {in0};
					replace = {out0};
				} else {
					//adjacent stitch:
					find = {yarn.loop, in0};
					replace = {yarn.loop, out0};
					if (s.direction != Stitch::CCW) {
						std::swap(find[0], find[1]);
						std::swap(replace[0], replace[1]);
					}
				}

				//update active cycles:
				make_cycle(find, replace);
				assert(cycles.size() == si + 1); //make sure cycle was made.

				//update all yarns at in0:
				for (auto &y : active_yarns) {
					if (y.second.loop == in0) {
						y.second.loop = out0;
					}
				}

				//update yarn for stitch:
				yarn.loop = out0;
				yarn.direction = s.direction;

			}

		}
		
		
	}

	return 0;
}
