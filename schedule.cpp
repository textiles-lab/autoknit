#include "Stitch.hpp"

#include "TaggedArguments.hpp"

#include <deque>
#include <map>
#include <set>
#include <algorithm>

//Loop held on a needle:
struct Loop {
	constexpr Loop(uint32_t stitch_, uint32_t idx_) : stitch(stitch_), idx(idx_) {
	}
	uint32_t stitch;
	uint32_t idx;
	bool operator==(Loop const &o) const {
		return (stitch == o.stitch && idx == o.idx);
	}
	bool operator!=(Loop const &o) const {
		return (stitch != o.stitch || idx != o.idx);
	}
	bool operator<(Loop const &o) const {
		if (stitch != o.stitch) return stitch < o.stitch;
		else return idx < o.idx;
	}
	std::string to_string() const {
		if (stitch == -1U && idx == -1U) return "GAP";
		else return std::to_string(stitch) + "_" + std::to_string(idx);
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

	//Cycles are 1-1 with stitches -- each stitch creates a new cycle from [part of] an old cycle

	std::vector< Cycle > cycles;

	#define REPORT_ERROR( X ) do { std::cerr << (X) << std::endl; exit(1); } while(0)

	{ //build cycles:
		//current yarn position w.r.t. loops:
		struct YarnInfo {
			Loop loop = GAP;
			char direction = '\0';
		};
		std::map< uint32_t, YarnInfo > active_yarns;

		//Cycle locations of all active loops:
		struct CycleIndex {
			CycleIndex() = default;
			CycleIndex(uint32_t cycle_, uint32_t index_) : cycle(cycle_), index(index_) { }
			uint32_t cycle = -1U;
			uint32_t index = -1U;
		};
		std::map< Loop, CycleIndex > active_loops;

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
		// (2a) merging cycles so that they are in the same cycle
		// (2b) splitting cycle if crossing (note: stitch_direction used to figure out which part to take)
		// (3) replacing them with the 'replace' loops
		auto make_cycle = [&cycles, &active_loops](std::vector< Loop > const &find, std::vector< Loop > const &replace /*, char stitch_direction */) {
			//(Not as general a function as it appears.)

			//helper: find loop using index
			auto find_loop = [&active_loops](Loop const &loop) -> CycleIndex const & {
				auto f = active_loops.find(loop);
				assert(f != active_loops.end() && "expecting an active loop");
				return f->second;
			};

			auto roll_to_back = [](Cycle &cycle, uint32_t index) {
				assert(index < cycle.size());
				std::rotate(cycle.begin(), cycle.begin() + index + 1, cycle.end());
			};

			auto roll_to_front = [](Cycle &cycle, uint32_t index) {
				assert(index < cycle.size());
				std::rotate(cycle.begin(), cycle.begin() + index, cycle.end());
			};

			//always will make a new cycle:
			cycles.emplace_back();
			Cycle &cycle = cycles.back();

			//Special case the easy things:
			if (find.size() == 1 && find[0] == GAP) {
				//finding *just* a GAP always creates a new cycle.
				cycle.assign(replace.begin(), replace.end());
			} else if (find.size() == 1 && find[0] != GAP) {
				//finding a non-gap is straightforward:
				CycleIndex ci = find_loop(find[0]);
				//grab cycle holding non-gap:
				cycle = cycles[ci.cycle];
				//perform replacement:
				cycle.erase(cycle.begin() + ci.index);
				cycle.insert(cycle.begin() + ci.index, replace.begin(), replace.end());
			} else if (find.size() == 2 && find[0] != GAP) {
				CycleIndex ci0 = find_loop(find[0]);

				//start with the ci0 cycle:
				cycle = cycles[ci0.cycle];

				assert(cycle.size() >= 2); //all cycles are at least size 2

				//arrange find[0] at the back:
				roll_to_back(cycle, ci0.index);
				assert(cycle.back() == find[0]);

				//first non-gap stitch at the front of cycle:
				auto non_gap = cycle.begin();
				while (non_gap != cycle.end() && *non_gap == GAP) {
					++non_gap;
				}
				assert(non_gap != cycle.end());

				if (cycle.front() == find[1]) {
					//all in one cycle; great!
					assert(cycle.back() == find[0] && cycle.front() == find[1]);
				} else if (*non_gap == find[1]) {
					//all in one cycle + closing a gap.
					cycle.erase(cycle.begin(), non_gap);
					assert(cycle.back() == find[0] && cycle.front() == find[1]);
				} else {
					//not all in one cycle, so...
					if (find[1] == GAP) { //...add gap:
						//TODO: add_bridge(cycle.back(), cycle.front())
						cycle.emplace_front(GAP);
						assert(cycle.back() == find[0] && cycle.front() == find[1]);
					} else {
						CycleIndex ci1 = find_loop(find[1]);
						if (ci1.cycle == ci0.cycle) { //...split cycle:
							//want find[0] find[1] to be CCW-ordered in new cycle, so erase the proper bits:
							cycle.erase(cycle.begin(), cycle.begin() + ci1.index);

							assert(cycle.back() == find[0] && cycle.front() == find[1]);
						} else { //...merge cycles:
							Cycle cycle1 = cycles[ci1.cycle];
							roll_to_front(cycle1, ci1.index);
							assert(cycle1.size() >= 2);
							assert(cycle1.front() == find[1]);
							//NOTE: might end up with doubled GAP, but that doesn't matter(?)

							//TODO: add_bridge(cycle.back(), cycle.front())
							//TODO: add_bridge(cycle1.back(), cycle1.front())
							cycle.insert(cycle.begin(), cycle1.begin(), cycle1.end());
							assert(cycle.back() == find[0] && cycle.front() == find[1]);
						}
					}
				}
				//remove pattern:
				assert(cycle.back() == find[0] && cycle.front() == find[1]);
				cycle.pop_back(); //remove find[0] from the back
				cycle.pop_front(); //remove find[1] from the front
				//insert replacement:
				cycle.insert(cycle.end(), replace.begin(), replace.end());
			} else {
				assert(0 && "make_cycle handles very few cases for find/replace");
			}

			assert(cycle.size() >= 2); //all cycles are at least size 2

			//remove everything in 'find' from active_loop:
			for (auto const &l : find) {
				if (l != GAP) {
					auto f = active_loops.find(l);
					assert(f != active_loops.end());
					active_loops.erase(f);
				}
			}
			//update everything in 'cycle' in active_loops:
			for (uint32_t index = 0; index < cycle.size(); ++index) {
				active_loops[cycle[index]] = CycleIndex(cycles.size()-1, index);
			}


		};

		for (uint32_t si = 0; si < stitches.size(); ++si) {
			Stitch const &s = stitches[si];
			YarnInfo &yarn = active_yarns[s.yarn];

			if (s.type == Stitch::Start) {
				assert(s.in[0] == -1U && s.in[1] == -1U && s.out[0] != -1U && s.out[0] > si && s.out[1] == -1U && "valid 0-1 stitch");
				Loop out0(si, 0);

				std::vector< Loop > find, replace;
				if (yarn.loop == GAP) {
					//bringing in yarn, I suppose:
					find = {GAP};
					replace = {out0, GAP};
				} else if (yarn.direction != s.direction) {
					//turning around:
					REPORT_ERROR("Can't turn around on start stitch.");
				} else {
					//adjacent stitch:
					find = {yarn.loop, GAP};
					replace = {yarn.loop, out0};
					if (s.direction != Stitch::CCW) {
						std::swap(find[0], find[1]);
						std::swap(replace[0], replace[1]);
					}
				}

				//DEBUG: dump find/replace:
				std::cout << "  find:";
				for (auto const &l : find) std::cout << ' ' << l.to_string();
				std::cout << '\n';
				std::cout << "  replace:";
				for (auto const &l : replace) std::cout << ' ' << l.to_string();
				std::cout << '\n';

				//update active cycles:
				make_cycle(find, replace);
				assert(cycles.size() == si + 1); //make sure cycle was made.

				//update yarn for stitch:
				yarn.loop = out0;
				yarn.direction = s.direction;

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

			
			//DEBUG: dump new cycle:
			std::cout << "cycles[" << (cycles.size()-1) << "]:";
			for (auto const &l : cycles.back()) {
				std::cout << ' ' << l.to_string();
			}
			std::cout << '\n';

		}
		
		
	}

	return 0;
}
