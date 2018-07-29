#include "Stitch.hpp"
#include "Shape.hpp"
#include "typeset.hpp"
#include "ScheduleCost.hpp"
#include "embed_DAG.hpp"

#include "TaggedArguments.hpp"

#include <deque>
#include <map>
#include <set>
#include <algorithm>
#include <list>

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
constexpr const Loop INVALID_LOOP = Loop(-1U, -1U);


struct CycleIndex {
	CycleIndex(uint32_t cycle_, uint32_t index_) : cycle(cycle_), index(index_) { }
	uint32_t cycle = -1U;
	uint32_t index = -1U;
	bool operator!=(CycleIndex const &o) const {
		return cycle != o.cycle && index != o.index;
	}
	std::string to_string() const {
		if (cycle == -1U && index == -1U) return ".";
		return std::to_string(int32_t(cycle)) + "-" + std::to_string(int32_t(index));
	}
	std::string to_string_simple() const {
		if (cycle == -1U && index == -1U) return ".";
		if (cycle < 26) {
			if (index == 0) {
				return std::string() + char('A' + cycle);
			} else if (index == 1) {
				return std::string() + '+';
			} else {
				return std::string() + char('a' + cycle);
			}
		} else {
			if (index == 0) {
				return std::string() + '*';
			} else {
				return std::string() + 'x';
			}
		}
	}
};


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
	// (1) split into "steps" ==> bits of knitting that will be done at once
	//    - each step eventually needs a shape + roll + offset for its output loops
	//      (implies a shape for input loops)
	// (2) pick a consistent shape + roll for "interesting" steps
	//    - these are steps that take loops from more than one output or input
	//    - effectively, this finds an upward planar embedding, where the edges are chains of construction steps and the vertices occur when steps have more than one tube as a parent or child.
	// (2) figure out a layout (shape + roll) for each step

	//Construction will proceed in "Steps".
	// Each step will:
	//  a) shift loops stored on the bed, if needed.
	//  b) split/merge cycles if needed.
	//  c) reshape stored loops (roll/increase/decrease)
	//  d) create some new stitches

	// In other words, we can build a graph where the edges are loops stored on the bed, and the vertices are steps.
	// Scheduling means assigning compatible shapes to each edge (equiv: each shape's output)

	//Storage connects steps:
	typedef uint32_t StorageIdx;
	struct Storage : std::deque< Loop > {
		//Shape shape;
	};

	struct ScheduleOption {
		ScheduleCost cost;
		std::vector< uint32_t > in_order;
		std::vector< PackedShape > in_shapes;
		PackedShape inter_shape;
		std::vector< uint32_t > out_order;
		std::vector< PackedShape > out_shapes;
	};

	struct Step {
		uint32_t begin = 0, end = 0; //stitch range constructed in this step
		std::vector< StorageIdx > in, out;
		Storage inter; //configuration of loops after input arrangement but before performing stitches
		std::vector< uint32_t > inter_to_out; //how stitches in inter are transformed to get to out

		//Step reinterprets storages[in[*]] as inter + storages[out[1...]]
		// stretches and rolls inter,
		// and performs knitting on inter to produce storages[out[0]]

		std::vector< ScheduleOption > options;
	};

	#define REPORT_ERROR( X ) do { std::cerr << (X) << std::endl; exit(1); } while(0)

	std::vector< Storage > storages;
	std::vector< Step > steps;

	storages.reserve(2 * stitches.size()); //just how much storage can there actually be?
	steps.reserve(stitches.size());

	//get the out loop for a given 'in':
	auto source_loop = [&stitches](uint32_t si, uint32_t ii) {
		assert(stitches[si].in[ii] < stitches.size());
		Loop ret(stitches[si].in[ii], 0);
		if (stitches[ret.stitch].out[ret.idx] != si) ret.idx = 1;
		assert(stitches[ret.stitch].out[ret.idx] == si);
		return ret;
	};

	{ //build steps:
		//current yarn position w.r.t. active loops:
		struct YarnInfo {
			Loop loop = INVALID_LOOP;
			char direction = '\0';
		};
		std::map< uint32_t, YarnInfo > active_yarns;

		std::map< Loop, StorageIdx > active_loops;

		struct Bridge {
			Loop a,b;
			uint32_t distance = 0;
		};
		std::vector< Bridge > bridges;

		uint32_t next_begin = 0;
		while (next_begin < stitches.size()) {
			//std::cout << "steps[" << steps.size() << "]:\n"; //DEBUG

			steps.emplace_back();
			Step &step = steps.back();

			{ //figure out begin/end range for step:
				step.begin = next_begin;
				step.end = step.begin + 1;
				//steps take as many stitches as they can, given some constraints:
				uint32_t increases = 0, decreases = 0;
				uint32_t last_shaping = step.begin;
				if (stitches[step.begin].type == Stitch::Increase) {
					++increases;
				}
				if (stitches[step.begin].type == Stitch::Decrease) {
					++decreases;
				}
				auto basic_type = [](Stitch const &s) {
					if (s.type == Stitch::Start) return '+';
					else if (s.type == Stitch::End) return '-';
					else return '|';
				};
				while (step.end < stitches.size()) {
					//same basic type:
					if (basic_type(stitches[step.end]) != basic_type(stitches[step.begin])) break;
					//same direction:
					if (stitches[step.end].direction != stitches[step.begin].direction) break;
					//same yarn:
					if (stitches[step.end].yarn != stitches[step.begin].yarn) break;
					//same course: (not sure this is really needed, but it does mean that every loop ends up in a Storage at some point)
					if (stitches[step.end].in[0] != -1U && stitches[step.end].in[0] >= step.begin) break;
					if (stitches[step.end].in[1] != -1U && stitches[step.end].in[1] >= step.begin) break;
					//not too many decreases/increases in a segment (...and, if there might be, cut things back a bit)
					if (stitches[step.end].type == Stitch::Increase) {
						++increases;
						if (increases >= 4) {
							step.end = (step.end - last_shaping + 1) / 2 + last_shaping;
							break;
						}
						last_shaping = step.end;
					}
					if (stitches[step.end].type == Stitch::Decrease) {
						++decreases;
						if (decreases >= 4) {
							step.end = (step.end - last_shaping + 1) / 2 + last_shaping;
							break;
						}
						last_shaping = step.end;
					}
					++step.end;
				}
				next_begin = step.end;
			}

			std::vector< Loop > in_chain;
			std::vector< Loop > out_chain;

			//stitches take in_chain -> out_chain:
			for (uint32_t si = step.begin; si != step.end; ++si) {
				for (uint32_t i = 0; i < 2; ++i) {
					if (stitches[si].in[i] == -1U) continue;
					Loop l = source_loop(si, i);
					in_chain.emplace_back(l);
				}
				for (uint32_t o = 0; o < 2; ++o) {
					if (stitches[si].out[o] == -1U) continue;
					out_chain.emplace_back(si, o);
				}
			}

			if (stitches[step.begin].direction == Stitch::CW) {
				std::reverse(in_chain.begin(), in_chain.end());
				std::reverse(out_chain.begin(), out_chain.end());
			}

			/*{ //DEBUG
				std::cout << "  in:";
				for (auto const &l : in_chain) {
					std::cout << " " << l.to_string();
				}
				std::cout << '\n';
				std::cout << "  out:";
				for (auto const &l : out_chain) {
					std::cout << " " << l.to_string();
				}
				std::cout << '\n';
				std::cout.flush();
			}*/

			YarnInfo &yarn = active_yarns[stitches[step.begin].yarn];

			{ //fill in step's 'in' array based on storages holding in_chain:
				std::set< StorageIdx > used;
				for (auto const &l : in_chain) {
					auto f = active_loops.find(l);
					assert(f != active_loops.end() && "stitches should depend only on active stuff");
					used.insert(f->second);
				}
				if (yarn.loop != INVALID_LOOP) {
					auto f = active_loops.find(yarn.loop);
					assert(f != active_loops.end() && "active yarn should reference active stuff");
					used.insert(f->second);
				}

				step.in.assign(used.begin(), used.end());
			}

			{ //do some loop surgery to figure out 'out' storages:
				std::list< Storage > outs;
				for (auto storage_idx : step.in) {
					outs.emplace_back(storages[storage_idx]);
				}

				//flip and re-jigger yarn storages:
				auto link_ccw = [&outs](Loop const &a, Loop const &b) {
					std::list< Storage >::iterator sa = outs.end();
					uint32_t la = -1U;
					std::list< Storage >::iterator sb = outs.end();
					uint32_t lb = -1U;

					for (auto oi = outs.begin(); oi != outs.end(); ++oi) {
						for (uint32_t li = 0; li < oi->size(); ++li) {
							if ((*oi)[li] == a) {
								assert(sa == outs.end());
								sa = oi;
								assert(la == -1U);
								la = li;
							}
							if ((*oi)[li] == b) {
								assert(sb == outs.end());
								sb = oi;
								assert(lb == -1U);
								lb = li;
							}
						}
					}
					assert(sa != outs.end() && la < sa->size());
					assert(sb != outs.end() && lb < sb->size());

					if (sa == sb) {
						if ((la + 1) % sa->size() == lb) {
							//already ccw! great.
						} else {
							//must split loop.
							Storage non_ab;
							Storage ab;
							if (la < lb) {
								ab.insert(ab.end(), sa->begin(), sa->begin() + la + 1);
								non_ab.insert(non_ab.end(), sa->begin() + la + 1, sa->begin() + lb);
								ab.insert(ab.end(), sa->begin() + lb, sa->end());
							} else { assert(lb < la);
								non_ab.insert(non_ab.end(), sa->begin(), sa->begin() + lb);
								ab.insert(ab.end(), sa->begin() + lb, sa->begin() + la + 1);
								non_ab.insert(non_ab.end(), sa->begin() + la + 1, sa->end());
							}
							//PARANOIA:
							assert(non_ab.size() + ab.size() == sa->size());
							la = -1U;
							lb = -1U;
							for (uint32_t li = 0; li < ab.size(); ++li) {
								if (ab[li] == a) {
									assert(la == -1U);
									la = li;
								}
								if (ab[li] == b) {
									assert(lb == -1U);
									lb = li;
								}
							}
							assert(la != -1U && lb != -1U);
							assert((la + 1) % ab.size() == lb);
							for (auto const &l : non_ab) {
								assert(l != a && l != b);
							}
							//end PARANOIA
							*sa = std::move(ab);
							outs.emplace_back(non_ab);
						}
					} else {
						//must merge loops
						std::rotate(sa->begin(), sa->begin() + (la + 1), sa->end());
						assert(sa->back() == a);
						std::rotate(sb->begin(), sb->begin() + lb, sb->end());
						assert(sb->front() == b);
						//TODO: add_bridge(sa->front(), sa->back())
						//TODO: add_bridge(sb->front(), sb->back())
						sa->insert(sa->end(), sb->begin(), sb->end());
						outs.erase(sb);
					}

				};

				for (uint32_t i = 0; i + 1 < in_chain.size(); ++i) {
					link_ccw(in_chain[i], in_chain[i+1]);
				}

				if (yarn.loop != INVALID_LOOP) {
					if (yarn.direction != stitches[step.begin].direction) {
						//reversing direction should happen on the same stitch:
						assert(!in_chain.empty() && yarn.loop == (stitches[step.begin].direction == Stitch::CCW ? in_chain[0] : in_chain.back()));
					} else {
						assert(!in_chain.empty() && "Don't handle linking yarn into starts, though that might be useful.");
						if (yarn.direction == Stitch::CCW) { //formed ccw, so yarn link is before first elt of chain:
							link_ccw(yarn.loop, in_chain[0]);
						} else {
							link_ccw(in_chain.back(), yarn.loop);
						}
					}
				}

				//Now 'outs' should have in_chain in one ccw list.

				// --> replace with out_chain.
				if (in_chain.empty()) {
					//empty in chain -> create a new cycle~

					step.inter = Storage(); //inter is an empty cycle

					Storage result;
					result.insert(result.end(), out_chain.begin(), out_chain.end());
					outs.emplace_front(result);
				} else {
					//check that in_chain is indeed in order:
					std::list< Storage >::iterator s = outs.end();
					uint32_t l = -1U;

					for (auto oi = outs.begin(); oi != outs.end(); ++oi) {
						for (uint32_t li = 0; li < oi->size(); ++li) {
							if ((*oi)[li] == in_chain[0]) {
								assert(s == outs.end());
								s = oi;
								assert(l == -1U);
								l = li;
							}
						}
					}
					assert(s != outs.end());
					assert(l < s->size());
					std::rotate(s->begin(), s->begin() + l, s->end());

					//move 's' to outs[0]:
					std::swap(*outs.begin(), *s);
					s = outs.begin();

					//store the pre-transform version of outs[0] as the step's "inter":
					step.inter = *s;

					assert(s->size() >= in_chain.size());
					for (uint32_t i = 0; i < in_chain.size(); ++i) {
						assert((*s)[i] == in_chain[i]);
					}

					s->erase(s->begin(), s->begin() + in_chain.size());
					s->insert(s->begin(), out_chain.begin(), out_chain.end());
					if (s->empty()) {
						outs.erase(s);
						//TODO: special flag for steps whose inter doesn't actually turn into outs[0] ?
					}
				}
				uint32_t base = storages.size();
				storages.insert(storages.end(), outs.begin(), outs.end());
				for (uint32_t i = base; i < storages.size(); ++i) {
					step.out.push_back(i);
				}
			}

			if (step.in.size() > 0 && step.out.size() > 0) {
				//make a map from step.inter -> step.out[0]
				std::vector< uint32_t > inter_to_out(step.inter.size(), -1U);
		
				std::map< Loop, CycleIndex > inter_loops;
				std::map< Loop, CycleIndex > out_loops;

				for (uint32_t oi = 0; oi < step.out.size(); ++oi) {
					Storage const &inter = (oi == 0 ? step.inter : storages[step.out[oi]]);
					for (uint32_t i = 0; i < inter.size(); ++i) {
						auto const &l = inter[i];
						auto ret = inter_loops.insert(std::make_pair(l, CycleIndex(oi, i)));
						assert(ret.second);
					}
					Storage const &out = storages[step.out[oi]];
					for (uint32_t i = 0; i < out.size(); ++i) {
						auto const &l = out[i];
						auto ret = out_loops.insert(std::make_pair(l, CycleIndex(oi, i)));
						assert(ret.second);
					}
				}
		
				//First, the loops that actually change because of stitches:
				for (uint32_t si = step.begin; si != step.end; ++si) {
					Stitch const &stitch = stitches[si];
					if (stitch.type == Stitch::Tuck || stitch.type == Stitch::Miss || stitch.type == Stitch::Knit
					 || stitch.type == Stitch::Increase) {
						auto f = inter_loops.find(source_loop(si,0));
						assert(f != inter_loops.end());
						assert(f->second.cycle == 0 && f->second.index < inter_to_out.size());
	
						auto t = out_loops.find(Loop(si, 0));
						assert(t != out_loops.end());
						assert(t->second.cycle == 0 && t->second.index < storages[step.out[0]].size());
	
						assert(inter_to_out[f->second.index] == -1U);
						inter_to_out[f->second.index] = t->second.index;

					} else if (stitch.type == Stitch::Decrease) {
						auto t = out_loops.find(Loop(si, 0));
						assert(t != out_loops.end());
						assert(t->second.cycle == 0 && t->second.index < storages[step.out[0]].size());
						for (uint32_t i = 0; i < 2; ++i) {
							auto f = inter_loops.find(source_loop(si,i));
							assert(f != inter_loops.end());
							assert(f->second.cycle == 0 && f->second.index < inter_to_out.size());

							assert(inter_to_out[f->second.index] == -1U);
							inter_to_out[f->second.index] = t->second.index;
						}
					} else {
						std::cout << "What is '" << stitch.type << "'?" << std::endl; //DEBUG
						assert(0 && "Unsupported stitch type.");
					}
				}

				//add things in the cycle that aren't touched by stitches:
				for (uint32_t ii = 0; ii < step.inter.size(); ++ii) {
					if (inter_to_out[ii] != -1U) continue;
					auto t = out_loops.find(step.inter[ii]);
					assert(t != out_loops.end());
					assert(t->second.cycle == 0 && t->second.index < storages[step.out[0]].size());
					inter_to_out[ii] = t->second.index;
				}

				step.inter_to_out = inter_to_out;
			}

			/*{ //DEBUG
				for (StorageIdx s : step.in) {
					std::cout << "  uses storages[" << s << "]:";
					for (auto const &l : storages[s]) std::cout << " " << l.to_string();
					std::cout << "\n";
				}
				for (StorageIdx s : step.out) {
					std::cout << "  makes storages[" << s << "]:";
					for (auto const &l : storages[s]) std::cout << " " << l.to_string();
					std::cout << "\n";
				}
				std::cout.flush();
			}*/

			{ //update active yarn:
				if (out_chain.empty()) {
					//no outs -> take yarn out.
					auto f = active_yarns.find(stitches[step.begin].yarn);
					assert(f != active_yarns.end());
					active_yarns.erase(f);
				} else {
					yarn.loop = Loop(
						step.end-1,
						stitches[step.end-1].out[1] != -1U ? 1 : 0
					);
					assert(stitches[step.end-1].out[yarn.loop.idx] != -1U);
					yarn.direction = stitches[step.end-1].direction;
				}
			}

			{ //update ~other ~ yarns:
				std::map< Loop, Loop > cw_out;
				std::map< Loop, Loop > ccw_out;
				for (uint32_t si = step.begin; si < step.end; ++si) {
					Loop cw = Loop(si, -1U);
					Loop ccw = Loop(si, -1U);
					if (stitches[si].out[0] == -1U) {
						//leave with '-1U' for out idx
					} else if (stitches[si].out[1] == -1U) {
						cw.idx = ccw.idx = 0;
					} else {
						if (stitches[si].direction == Stitch::CCW) {
							cw.idx = 0;
							ccw.idx = 1;
						} else {
							cw.idx = 1;
							ccw.idx = 0;
						}
					}
					for (uint32_t i = 0; i < 2; ++i) {
						if (stitches[si].in[i] == -1U) continue;
						Loop from(stitches[si].in[i], stitches[stitches[si].in[0]].find_out(si));
						auto res = cw_out.insert(std::make_pair(from, cw));
						assert(res.second);
						res = ccw_out.insert(std::make_pair(from, ccw));
						assert(res.second);
					}

					//TODO: (this would be a good place to update bridges as well)
				}
				for (auto &i_ay : active_yarns) {
					auto &ay = i_ay.second;
					if (ay.direction == Stitch::CCW) {
						auto f = ccw_out.find(ay.loop);
						if (f != ccw_out.end()) {
							ay.loop = f->second;
						}
					} else {
						auto f = cw_out.find(ay.loop);
						if (f != cw_out.end()) {
							ay.loop = f->second;
						}
					}
					if (ay.loop.idx == -1U) {
						//TODO: need to take yarn out, probably; does this ever come up though?
						assert(false && "yarn-out case that probably never comes up");
					}
				}
			}

			{ //update active_loops:
				for (StorageIdx i : step.in) {
					for (auto const &l : storages[i]) {
						auto f = active_loops.find(l);
						assert(f != active_loops.end());
						active_loops.erase(f);
					}
				}
				for (StorageIdx i : step.out) {
					for (auto const &l : storages[i]) {
						assert(!active_loops.count(l));
						active_loops[l] = i;
					}
				}
			}

		} //while (more stitches)

	} //end of build steps


	//Figure out possible shapes for storages near *interesting* steps:
	std::cout << "Figuring out shapes for interesting steps:" << std::endl; //DEBUG
	for (auto &step : steps) {
		//interesting steps have more than one out/in:
		if (step.in.size() <= 1 && step.out.size() <= 1) continue;

		std::cout << "\nsteps[" << (&step - &steps[0]) << "]:\n"; //DEBUG

		{ //DEBUG
			for (StorageIdx s : step.in) {
				std::cout << "  uses storages[" << s << "]:";
				for (auto const &l : storages[s]) std::cout << " " << l.to_string();
				std::cout << "\n";
			}
			for (StorageIdx s : step.out) {
				std::cout << "  makes storages[" << s << "]:";
				for (auto const &l : storages[s]) std::cout << " " << l.to_string();
				std::cout << "\n";
			}
			std::cout.flush();
		}

		//Construction model:
		// input cycles come in some order and shape
		// (input cycles are perhaps reshaped) <-- "late reshape"
		// middle (intermediate) cycle is created
		// intermediate cycle is reshaped
		// knitting happens

	
		//get local copies of storages:
		std::vector< Storage > ins, outs, intermediates;
		for (auto si : step.in) {
			ins.emplace_back(storages[si]);
		}
		for (auto si : step.out) {
			outs.emplace_back(storages[si]);
		}
		intermediates.emplace_back(step.inter);
		for (uint32_t i = 1; i < step.out.size(); ++i) {
			intermediates.emplace_back(storages[step.out[i]]);
		}

		//The contents of 'inter + outs[1...]' and 'ins' should be the same. Construct a mapping:
		std::vector< std::vector< CycleIndex > > in_to_inter;
		{
			std::map< Loop, CycleIndex > loops;
			for (auto const &inter : intermediates) {
				for (uint32_t i = 0; i < inter.size(); ++i) {
					auto const &l = inter[i];
					auto ret = loops.insert(std::make_pair(l, CycleIndex(&inter - &intermediates[0], i)));
					assert(ret.second);
				}
			}
			uint32_t used = 0;
			in_to_inter.reserve(ins.size());
			for (auto const &in : ins) {
				in_to_inter.emplace_back();
				in_to_inter.back().reserve(in.size());
				for (auto const &l : in) {
					auto f = loops.find(l);
					assert(f != loops.end());
					in_to_inter.back().emplace_back(f->second);
					++used;
				}
				assert(in_to_inter.back().size() == in.size());
			}
			assert(in_to_inter.size() == ins.size());
			assert(used == loops.size());
		}

		//intermediates[0] (== step.inter) gets knit on and changes contents before its output:
		assert(step.inter_to_out.size() == step.inter.size());

		//TODO: ~cache of transfer costs~


		//enumerate input orders and shapes:
		//  -- for each order, test if the intermediate storages are all cycles
		//     then figure out cost to reshape for knitting (for all possible reshapes)
		//     add [cost, input order, input shapes, output shapes, output order] to possible shapes list
		std::vector< uint32_t > remaining;
		remaining.reserve(ins.size());
		for (uint32_t i = 0; i < ins.size(); ++i) {
			remaining.emplace_back(i);
		}
		std::vector< uint32_t > in_order;
		std::vector< PackedShape > in_shapes(ins.size(), -1U);
		std::vector< CycleIndex > front;
		std::vector< CycleIndex > back;

		//Given a current order, test to see if it leaves intermediates in a valid configuration:
		auto test_order = [&]() -> const char * {
			struct Info {
				uint32_t back_min = -1U;
				uint32_t back_max = 0;
				uint32_t front_min = -1U;
				uint32_t front_max = 0;
				uint32_t zero_index = -1U;
				bool zero_on_back = false;
			};
			std::vector< Info > inter_info(intermediates.size());
			for (uint32_t i = 0; i < back.size(); ++i) {
				auto const &ci = back[i];
				if (ci.cycle == -1U) continue; //gap
				assert(ci.cycle < inter_info.size());
				auto &info = inter_info[ci.cycle];
				info.back_min = std::min(info.back_min, i);
				info.back_max = std::max(info.back_max, i);
				if (ci.index == 0) {
					assert(info.zero_index == -1U);
					info.zero_index = i;
					info.zero_on_back = true;
				}
			}
			for (uint32_t i = 0; i < front.size(); ++i) {
				auto const &ci = front[i];
				if (ci.cycle == -1U) continue; //gap
				assert(ci.cycle < inter_info.size());
				auto &info = inter_info[ci.cycle];
				info.front_min = std::min(info.front_min, i);
				info.front_max = std::max(info.front_max, i);
				if (ci.index == 0) {
					assert(info.zero_index == -1U);
					info.zero_index = i;
					info.zero_on_back = false;
				}
			}

			std::vector< PackedShape > inter_shapes;
			inter_shapes.reserve(intermediates.size());
			for (auto const &info : inter_info) {
				uint32_t cycle = &info - &inter_info[0];
				//check for gaps inside cycle:
				uint32_t count = 0;
				if (info.back_min <= info.back_max) count += (info.back_max - info.back_min + 1);
				if (info.front_min <= info.front_max) count += (info.front_max - info.front_min + 1);
				if (count > intermediates[cycle].size()) {
					//gap found; bail out
					return "gap";
				}
				assert(count == intermediates[cycle].size()); //must have seen all loops
				assert(info.zero_index != -1U); //including the zero element, of course

				Shape shape(0,0);

				//check that things are aligned:
				if (info.back_min > info.back_max) {
					//only on front bed
					assert(info.front_min <= info.front_max);
					if (info.front_min != info.front_max) {
						return "flat (front)";
					}
					shape.nibbles = Shape::BackLeft;
					shape.roll = 0;
				} else if (info.front_min > info.front_max) {
					//only on back bed
					assert(info.back_min <= info.back_max);
					if (info.back_min != info.back_max) {
						return "flat (back)";
					}
					shape.nibbles = Shape::FrontLeft;
					shape.roll = 0;
				} else {
					//on both beds, so check alignment between beds:
					if (std::abs(int32_t(info.back_min) - int32_t(info.front_min)) > 1) {
						return "mis-aligned (left)";
					}
					if (std::abs(int32_t(info.back_max) - int32_t(info.front_max)) > 1) {
						return "mis-aligned (right)";
					}
					shape.nibbles = 
						  (info.back_min > info.front_min ? Shape::BackLeft : 0)
						| (info.back_min < info.front_min ? Shape::FrontLeft : 0)
						| (info.back_max < info.front_max ? Shape::BackRight : 0)
						| (info.back_max > info.front_max ? Shape::FrontRight : 0)
					;
					if (info.zero_on_back) {
						shape.roll = (info.front_max - info.front_min + 1) + (info.back_max - info.zero_index);
					} else {
						shape.roll = info.zero_index - info.front_min;
					}
				}

				std::vector< CycleIndex > data;
				data.reserve(intermediates[cycle].size());
				for (uint32_t i = 0; i < intermediates[cycle].size(); ++i) {
					data.emplace_back(cycle, i);
				}
				std::vector< CycleIndex > expected_front, expected_back;
				shape.append_to_beds(data, CycleIndex(-1U, -1U), &expected_front, &expected_back);

				/*//DEBUG
				std::cout << "    Testing:\n";
				std::string front_string, back_string;
				typeset_beds< CycleIndex >(expected_front, expected_back, [](CycleIndex const &ci){
					return ci.to_string_simple();
				}, "", &front_string, &back_string);
				std::cout << "      " << back_string << std::endl;
				std::cout << "      " << front_string << std::endl;
				*/

				for (uint32_t i = info.back_min; i <= info.back_max; ++i) {
					uint32_t r = i - std::min(info.back_min, info.front_min);
					assert(r < expected_back.size());
					if (back[i] != expected_back[r]) {
						return "bad shape (back)";
					}
				}
				for (uint32_t i = info.front_min; i <= info.front_max; ++i) {
					uint32_t r = i - std::min(info.back_min, info.front_min);
					assert(r < expected_front.size());
					if (front[i] != expected_front[r]) {
						return "bad shape (front)";
					}
				}

				inter_shapes.emplace_back(shape.pack());
			}

			std::vector< uint32_t > out_order;
			out_order.reserve(outs.size());
			{ //figure out order of outs:
				std::vector< bool > seen(outs.size(), false);
				auto do_out = [&seen, &out_order](uint32_t cycle) {
					if (cycle == -1U) return;
					assert(cycle < seen.size());
					if (!seen[cycle]) {
						seen[cycle] = true;
						out_order.emplace_back(cycle);
					}
				};
				for (uint32_t i = 0; i < std::max(front.size(), back.size()); ++i) {
					if (i < front.size()) do_out(front[i].cycle);
					if (i < back.size()) do_out(back[i].cycle);
				}
				//make sure we saw all out cycles:
				for (auto s : seen) {
					assert(s);
				}
				assert(out_order.size() == outs.size());
			}


			//Order passed the test, now add possible transforms with their costs:

			/*//DEBUG:
			std::cout << "Have a working order:" << std::endl;
			std::string front_string, back_string;
			typeset_beds< CycleIndex >(front, back, [](CycleIndex const &ci){
				return ci.to_string();
			}, " ", &front_string, &back_string);
			std::cout << "  " << back_string << std::endl;
			std::cout << "  " << front_string << std::endl;
			*/


			std::vector< Shape > out0_shapes = Shape::make_shapes_for(outs[0].size());
			for (auto const &out0_shape : out0_shapes) {
				ScheduleCost cost;
				for (auto const &inter_shape : inter_shapes) {
					cost += ScheduleCost::shape_cost(Shape::unpack(inter_shape));
				}
				cost += ScheduleCost::shape_cost(out0_shape);
				cost += ScheduleCost::transfer_cost(
					intermediates[0].size(), Shape::unpack(inter_shapes[0]),
					outs[0].size(), out0_shape,
					step.inter_to_out);

				ScheduleOption option;
				option.cost = cost;

				option.in_order = in_order;
				option.in_shapes = in_shapes;

				option.inter_shape = inter_shapes[0];

				option.out_order = out_order;
				option.out_shapes = inter_shapes;
				option.out_shapes[0] = out0_shape.pack();

				step.options.emplace_back(option);
			}

			return nullptr;

		};

		std::function< void() > enumerate_orders = [&]() {
			if (remaining.empty()) {

				/*//DEBUG:
				std::cout << "  Trying:\n";
				std::string front_string, back_string;
				typeset_beds< CycleIndex >(front, back, [](CycleIndex const &ci){
					return ci.to_string_simple();
				}, "", &front_string, &back_string);
				std::cout << "    " << back_string << std::endl;
				std::cout << "    " << front_string << std::endl;
				*/


				test_order();
				/*
				const char *reason =
				if (reason != nullptr) { //DEBUG
					std::cout << "   FAILED: " << reason << "." << std::endl;
				}*/
				return;
			}
			//order isn't done yet, try all shapes:
			for (uint32_t r = 0; r < remaining.size(); ++r) {
				std::swap(remaining[r], remaining.back());
				uint32_t idx = remaining.back();
				remaining.pop_back();
				in_order.emplace_back(idx);

				std::vector< Shape > shapes = Shape::make_shapes_for(ins[idx].size());
				for (auto const &shape : shapes) {
					in_shapes[idx] = shape.pack();

					uint32_t front_size_before = front.size();
					uint32_t back_size_before = back.size();

					shape.append_to_beds(in_to_inter[idx], CycleIndex(-1U, -1U), &front, &back);

					enumerate_orders();

					front.erase(front.begin() + front_size_before, front.end());
					back.erase(back.begin() + back_size_before, back.end());

					in_shapes[idx] = -1U;
				}

				assert(!in_order.empty() && in_order.back() == idx);
				in_order.pop_back();
				remaining.push_back(idx);
				std::swap(remaining[r], remaining.back());
			}
		};

		enumerate_orders();

		std::cout << "In total, step had " << step.options.size() << " scheduling options." << std::endl;

	} //interesting steps
	std::cout << "(done)" << std::endl;


	//Figure out possible shapes for storages near *boring* steps:
	std::cout << "Figuring out shapes for boring steps:" << std::endl; //DEBUG
	for (auto &step : steps) {
		if (!(step.in.size() <= 1 && step.out.size() <= 1)) continue; //skip exciting steps

		if (step.in.size() == 0 && step.out.size() == 1) {
			//Starts a cycle -- can do so in any order.
			assert(stitches[step.begin].type == Stitch::Start);

			std::vector< Shape > out_shapes = Shape::make_shapes_for(storages[step.out[0]].size());
			ScheduleOption option;
			option.in_order.clear();
			option.in_shapes.clear();
			option.inter_shape = 0;
			option.out_order = {0};
			for (auto const &out_shape : out_shapes) {
				option.cost = ScheduleCost::shape_cost(out_shape);
				option.out_shapes = {out_shape.pack()};
				step.options.emplace_back(option);
			}
			
		} else if (step.in.size() == 1 && step.out.size() == 0) {
			//Ends a cycle -- can do so in any order.
			assert(stitches[step.begin].type == Stitch::End);

			std::vector< Shape > in_shapes = Shape::make_shapes_for(storages[step.in[0]].size());
			ScheduleOption option;
			option.out_order.clear();
			option.out_shapes.clear();
			option.in_order = {0};
			for (auto const &in_shape : in_shapes) {
				option.cost = ScheduleCost::shape_cost(in_shape);
				option.inter_shape = in_shape.pack();
				option.in_shapes = {in_shape.pack()};
				step.options.emplace_back(option);
			}

		} else { assert(step.in.size() == 1 && step.out.size() == 1);

			std::vector< Shape > in_shapes = Shape::make_shapes_for(storages[step.in[0]].size());
			std::vector< Shape > out_shapes = Shape::make_shapes_for(storages[step.out[0]].size());
			assert(step.inter.size() == storages[step.in[0]].size());
			assert(step.inter.size() == step.inter_to_out.size());

			ScheduleOption option;
			option.in_order = {0};
			option.out_order = {0};
			for (auto const &in_shape : in_shapes) {
				option.in_shapes = {in_shape.pack()};
				option.inter_shape = in_shape.pack();
				for (auto const &out_shape : out_shapes) {
					option.out_shapes = {out_shape.pack()};

					option.cost = ScheduleCost::shape_cost(in_shape);
					option.cost += ScheduleCost::shape_cost(out_shape);
					option.cost += ScheduleCost::transfer_cost(
						step.inter.size(), in_shape,
						storages[step.out[0]].size(), out_shape,
						step.inter_to_out
					);

					step.options.emplace_back(option);
				}
			}

		}
		std::cout << "steps[" << (&step - &steps[0]) << "] had " << step.options.size() << " scheduling options." << std::endl;
	}

	//---------------------
	{ //Next up: do upward-planar embedding.

		//For every 'interesting' step, build a DAGNode
		//For every chain of 1-in / 1-out steps, build a DAGEdge

		std::vector< DAGNode > nodes;
		std::vector< DAGEdge > edges;

		std::map< uint32_t, DAGEdgeIndex > active_edges; //maps from storages to edges currently in progress

		std::vector< uint32_t > node_steps;

		//dag edges corresponding to each step's ins/outs:
		std::vector< std::vector< DAGEdgeIndex > > step_in_edges;
		std::vector< std::vector< DAGEdgeIndex > > step_out_edges;
		step_in_edges.reserve(steps.size());
		step_out_edges.reserve(steps.size());

		std::cout << "Building DAG "; std::cout.flush();
		for (auto const &step : steps) {
			std::cout << "."; std::cout.flush();
			step_in_edges.emplace_back();
			step_out_edges.emplace_back();
			if (step.in.size() == 1 && step.out.size() == 1) {
				//1-1 step; accumulate DAGEdge cost matrix:
				auto f = active_edges.find(step.in[0]);
				assert(f != active_edges.end());
				assert(f->second < edges.size());
				DAGEdge &edge = edges[f->second];
				assert(edge.to == step.in[0]);
				active_edges.erase(f);

				std::vector< DAGCost > old_costs = std::move(edge.costs);
				uint32_t old_to_shapes = edge.to_shapes;
				assert(edge.to_shapes == Shape::count_shapes_for(storages[step.in[0]].size()));

				edge.to = step.out[0];
				edge.to_shapes = Shape::count_shapes_for(storages[step.out[0]].size());
				edge.costs.resize(edge.from_shapes * edge.to_shapes, ScheduleCost::max());

				//accumulate step costs into edge costs:
				for (auto const &option : step.options) {
					assert(option.in_order.size() == 1);
					assert(option.in_shapes.size() == 1);
					assert(option.out_order.size() == 1);
					assert(option.out_shapes.size() == 1);

					uint32_t in_shape = Shape::unpack(option.in_shapes[0]).index_for(storages[step.in[0]].size());
					uint32_t out_shape = Shape::unpack(option.out_shapes[0]).index_for(storages[step.out[0]].size());

					for (uint32_t f = 0; f < edge.from_shapes; ++f) {
						DAGCost &cost = edge.costs[f * edge.to_shapes + out_shape];
						DAGCost test = old_costs[f * old_to_shapes + in_shape];
						if (!(test == DAGCost::max())) {
							test += option.cost;
							cost = std::min(cost, test);
						}
					}
				}

				//step is <within> this edge:
				step_in_edges.back().emplace_back(&edge - &edges[0]);
				step_out_edges.back().emplace_back(&edge - &edges[0]);

				auto ret = active_edges.insert(std::make_pair(step.out[0], &edge - &edges[0]));
				assert(ret.second);
			} else { assert(step.in.size() != 1 || step.out.size() != 1);
				//Exciting step; build a DAGNode, consume some DAGEdges, start some DAGEdges.

				node_steps.emplace_back(&step - &steps[0]);
				nodes.emplace_back();
				DAGNode &node = nodes.back();

				std::vector< DAGEdgeIndex > in_edges;
				std::vector< DAGEdgeIndex > out_edges;

				//look up indices of in edges:
				for (uint32_t i = 0; i < step.in.size(); ++i) {
					auto f = active_edges.find(step.in[i]);
					assert(f != active_edges.end());
					DAGEdge &edge = edges[f->second];
					assert(edge.to == step.in[i]);
					edge.to = &node - &nodes[0];
					in_edges.emplace_back(f->second);
					active_edges.erase(f);
				}

				//create out edges:
				for (uint32_t i = 0; i < step.out.size(); ++i) {
					auto ret = active_edges.insert(std::make_pair(step.out[i], edges.size()));
					assert(ret.second);
					out_edges.emplace_back(edges.size());

					edges.emplace_back();
					DAGEdge &edge = edges.back();
					edge.from = &node - &nodes[0];
					edge.from_shapes = Shape::count_shapes_for(storages[step.out[i]].size());
					edge.to = step.out[i];
					edge.to_shapes = Shape::count_shapes_for(storages[step.out[i]].size());

					//NOTE: this ends up costing a bit of extra time, because multiple-source-shape edges always need a square cost matrix, but we don't actually have any costs to put into it, so we just make up a matrix that does not permit rotation:
					edge.costs.resize(edge.from_shapes * edge.to_shapes, ScheduleCost::max());
					assert(edge.from_shapes == edge.to_shapes);
					for (uint32_t i = 0; i < edge.from_shapes; ++i) {
						edge.costs[i * edge.to_shapes + i] = ScheduleCost::zero();
					}
				}

				step_in_edges.back() = in_edges;
				step_out_edges.back() = out_edges;

				//PARANOIA: no duplicate in_edges or out_edges, right?
				assert(std::set< uint32_t >(in_edges.begin(), in_edges.end()).size() == in_edges.size());
				assert(std::set< uint32_t >(out_edges.begin(), out_edges.end()).size() == out_edges.size());

				node.options.reserve(step.options.size());
				for (auto const &option : step.options) {
					node.options.emplace_back();
					DAGOption &dag_option = node.options.back();
					dag_option.cost = option.cost;
					dag_option.in_order.reserve(step.in.size());
					dag_option.in_shapes.reserve(step.in.size());
					dag_option.out_order.reserve(step.out.size());
					dag_option.out_shapes.reserve(step.out.size());
					for (uint32_t i = 0; i < step.in.size(); ++i) {
						uint32_t in = option.in_order[i];
						assert(in < in_edges.size());
						dag_option.in_order.emplace_back(in_edges[in]);
						assert(in < step.in.size());
						dag_option.in_shapes.emplace_back(Shape::unpack(option.in_shapes[in]).index_for(storages[step.in[in]].size()));
						assert(dag_option.in_shapes.back() < edges[in_edges[in]].to_shapes); //DEBUG
					}
					assert(dag_option.in_order.size() == step.in.size());
					assert(dag_option.in_shapes.size() == step.in.size());
					for (uint32_t o = 0; o < step.out.size(); ++o) {
						uint32_t out = option.out_order[o];
						assert(out < out_edges.size());
						dag_option.out_order.emplace_back(out_edges[out]);
						assert(out < step.out.size());
						dag_option.out_shapes.emplace_back(Shape::unpack(option.out_shapes[out]).index_for(storages[step.out[out]].size()));
						assert(dag_option.out_shapes.back() < edges[out_edges[out]].from_shapes); //DEBUG
					}
					assert(dag_option.out_order.size() == step.out.size());
					assert(dag_option.out_shapes.size() == step.out.size());
					//PARANOIA: no duplicates in order, right?
					assert(std::set< uint32_t >(dag_option.in_order.begin(), dag_option.in_order.end()).size() == dag_option.in_order.size());
					assert(std::set< uint32_t >(dag_option.out_order.begin(), dag_option.out_order.end()).size() == dag_option.out_order.size());
				}
				assert(node.options.size() == step.options.size());
			}
		} //for(steps)
		std::cout << " done." << std::endl;

		assert(node_steps.size() == nodes.size());

		std::cout << "Have " << edges.size() << " edges and " << nodes.size() << " nodes." << std::endl;

		std::vector< uint32_t > node_options;
		std::vector< int32_t > node_positions, edge_positions;

		if (!embed_DAG(nodes, edges, &node_options, &node_positions, &edge_positions)) {
			std::cerr << "ERROR: failed to find an upward-planar embedding." << std::endl;
			return 1;
		}

		assert(node_options.size() == nodes.size());
		assert(node_positions.size() == nodes.size());
		assert(edge_positions.size() == edges.size());

		//Go through steps and assign selected option:
		std::vector< uint32_t > step_options(steps.size(), -1U);
		assert(node_options.size() == nodes.size());
		for (uint32_t n = 0; n < nodes.size(); ++n) {
			assert(node_steps[n] < steps.size());
			std::cout << "Step " << node_steps[n] << " gets option " << node_options[n] << " of " << nodes[n].options.size() << std::endl;
			assert(node_options[n] < steps[node_steps[n]].options.size());
			assert(nodes[n].options.size() == steps[node_steps[n]].options.size());
			step_options[node_steps[n]] = node_options[n];
		}

		//assign storage positions:
		std::vector< int32_t > storage_positions(storages.size(), std::numeric_limits< int32_t >::max());
		for (auto const &step : steps) {
			uint32_t si = &step - &steps[0];
			if (step_in_edges[si].size() == 1 && step_out_edges[si].size() == 1) {
				assert(step_options[si] == -1U); //wasn't part of the dag nodes.
				assert(step_in_edges[si][0] == step_out_edges[si][0]);
			}
			assert(step_in_edges[si].size() == step.in.size());
			for (uint32_t in = 0; in < step.in.size(); ++in) {
				assert(step_in_edges[si][in] < edge_positions.size());
				assert(storage_positions[step.in[in]] == edge_positions[step_in_edges[si][in]]);
			}
			assert(step_out_edges[si].size() == step.out.size());
			for (uint32_t out = 0; out < step.out.size(); ++out) {
				assert(step_out_edges[si][out] < edge_positions.size());
				assert(storage_positions[step.out[out]] == std::numeric_limits< int32_t >::max());
				storage_positions[step.out[out]] = edge_positions[step_out_edges[si][out]];
			}
		}

		//build up cost matrices to figure out all unassigned options:
		struct ShapeCost {
			ScheduleCost cost = ScheduleCost::max();
			uint32_t option = -1U; //option that cost comes via
		};
		std::vector< std::vector< ShapeCost > > storage_costs;
		storage_costs.reserve(storages.size());
		for (auto const &storage : storages) {
			storage_costs.emplace_back(Shape::count_shapes_for(storage.size()));
		}
		std::vector< uint32_t > storage_steps(storages.size(), -1U);
		for (auto const &step : steps) {
			for (auto o : step.out) {
				assert(o < storage_steps.size());
				assert(storage_steps[o] == -1U);
				storage_steps[o] = &step - &steps[0];
			}
		}
		for (auto const &step : steps) {
			uint32_t si = &step - &steps[0];
			if (step_options[si] != -1U) {
				//already-assigned option:
				auto const &option = step.options[step_options[si]];
				//chase back storage_costs chains from ins:
				for (uint32_t in = 0; in < step.in.size(); ++in) {
					uint32_t storage = step.in[in];
					uint32_t idx = Shape::unpack(option.in_shapes[in]).index_for(storages[step.in[in]].size());
					while (true) {
						uint32_t si = storage_steps[storage];
						assert(si < steps.size());
						if (step_options[si] != -1U) break; //finished
						assert(steps[si].in.size() == 1 && steps[si].out.size() == 1);
						assert(steps[si].out[0] == storage);

						assert(!(storage_costs[storage][idx].cost == ScheduleCost::max()));
						uint32_t opt = storage_costs[storage][idx].option;
						assert(opt < steps[si].options.size());

						uint32_t in_idx = Shape::unpack(steps[si].options[opt].in_shapes[0]).index_for(storages[steps[si].in[0]].size());
						uint32_t out_idx = Shape::unpack(steps[si].options[opt].out_shapes[0]).index_for(storages[steps[si].out[0]].size());
						assert(out_idx == idx);

						step_options[si] = opt;

						idx = in_idx;
						storage = steps[si].in[0];
					}
				}
				//start storage_costs chains from outs:
				for (uint32_t out = 0; out < step.out.size(); ++out) {
					uint32_t idx = Shape::unpack(option.out_shapes[out]).index_for(storages[step.out[out]].size());
					assert(idx < storage_costs[step.out[out]].size());
					storage_costs[step.out[out]][idx].cost = ScheduleCost::zero();
					storage_costs[step.out[out]][idx].option = step_options[si];
				}
			} else {
				//was part of an edge, need to propagate:
				assert(step.in.size() == 1 && step.out.size() == 1);

				for (auto const &option : step.options) {
					assert(option.in_order.size() == 1);
					assert(option.in_shapes.size() == 1);
					assert(option.out_order.size() == 1);
					assert(option.out_shapes.size() == 1);

					uint32_t in_shape = Shape::unpack(option.in_shapes[0]).index_for(storages[step.in[0]].size());
					uint32_t out_shape = Shape::unpack(option.out_shapes[0]).index_for(storages[step.out[0]].size());

					assert(in_shape < storage_costs[step.in[0]].size());
					ScheduleCost test = storage_costs[step.in[0]][in_shape].cost;
					if (!(test == ScheduleCost::max())) {
						test += option.cost;
						assert(out_shape < storage_costs[step.out[0]].size());
						if (test < storage_costs[step.out[0]][out_shape].cost) {
							storage_costs[step.out[0]][out_shape].cost = test;
							storage_costs[step.out[0]][out_shape].option = &option - &step.options[0];
						}
					}
				}
			}
		}

	
		//record shapes:
		std::vector< PackedShape > storage_shapes(storages.size(), -1U);
		for (uint32_t si = 0; si < steps.size(); ++si) {
			assert(step_options[si] < steps[si].options.size());
			auto const &option = steps[si].options[step_options[si]];
			for (uint32_t in = 0; in < steps[si].in.size(); ++in) {
				assert(storage_shapes[steps[si].in[in]] == option.in_shapes[in]);
			}
			for (uint32_t out = 0; out < steps[si].out.size(); ++out) {
				assert(storage_shapes[steps[si].out[out]] == -1U);
				storage_shapes[steps[si].out[out]] = option.out_shapes[out];
			}
			//PARANOIA: check order vs storage positions:
			for (uint32_t i = 1; i < option.in_order.size(); ++i) {
				assert(option.in_order[i-1] < steps[si].in.size());
				assert(option.in_order[i] < steps[si].in.size());
				assert(
					storage_positions[steps[si].in[option.in_order[i-1]]]
					< storage_positions[steps[si].in[option.in_order[i]]]
				);
			}
			for (uint32_t i = 1; i < steps[si].options[step_options[si]].out_order.size(); ++i) {
				assert(option.out_order[i-1] < steps[si].out.size());
				assert(option.out_order[i] < steps[si].out.size());
				assert(
					storage_positions[steps[si].out[option.out_order[i-1]]]
					< storage_positions[steps[si].out[option.out_order[i]]]
				);
			}
		}

		//TODO: PARANOIA: check that step option orders match recorded storage positions.

		auto summarize = [](Shape const &shape, uint32_t size, uint32_t mark) -> std::string {
			//shape as braille dots
			// e.g.  ⠲⠶⠷⠶
			std::vector< char > data(size, '.');
			if (mark < size) data[mark] = '*';
			std::vector< char > front, back;
			shape.append_to_beds(data, ' ', &front, &back);

			if (front.size() % 2 == 1) front.emplace_back(' ');
			if (back.size() % 2 == 1) back.emplace_back(' ');

			std::vector< uint8_t > dots(std::max(front.size(), back.size())/2, 0);
			for (uint32_t i = 0; i + 1 < front.size(); i += 2) {
				if (front[i]   == '.') dots[i/2] |= 0x4;
				if (front[i+1] == '.') dots[i/2] |= 0x20;
				if (front[i]   == '*') dots[i/2] |= 0x4 | 0x40;
				if (front[i+1] == '*') dots[i/2] |= 0x20 | 0x80;
			}
			for (uint32_t i = 0; i + 1 < back.size(); i += 2) {
				if (back[i]   == '.') dots[i/2] |= 0x2;
				if (back[i+1] == '.') dots[i/2] |= 0x10;
				if (back[i]   == '*') dots[i/2] |= 0x2 | 0x1;
				if (back[i+1] == '*') dots[i/2] |= 0x10 | 0x8;
			}

			std::string ret;
			for (auto d : dots) {
				uint16_t c = 0x2800 | d;
				ret += char(0xe0 | ((c >> 12) & 0x1f));
				ret += char(0x80 | ((c >> 6) & 0x3f));
				ret += char(0x80 | ((c) & 0x3f));
			}
			return ret;
		};

		//DEBUG: dump selections:
		for (auto const &step : steps) {
			uint32_t si = &step - &steps[0];
			std::cout << "Step " << si << " gets option " << step_options[si] << " of " << steps[si].options.size() << "\n";
			auto const &option = steps[si].options[step_options[si]];
			std::cout << "  takes";
			for (auto i : option.in_order) {
				uint32_t in = step.in[i];
				uint32_t mark = 0;
				for (uint32_t s = 1; s < storages[in].size(); ++s) {
					if (storages[in][mark] < storages[in][s]) mark = s;
				}
				std::cout << " s" << in << " in shape " << storage_shapes[in] << " " << summarize(Shape::unpack(storage_shapes[in]), storages[in].size(), mark) << " in lane " << storage_positions[in] << ",";
			}
			std::cout << "\n";
			std::cout << "  makes";
			for (auto o : option.out_order) {
				uint32_t out = step.out[o];
				uint32_t mark = 0;
				for (uint32_t s = 1; s < storages[out].size(); ++s) {
					if (storages[out][mark] < storages[out][s]) mark = s;
				}
				std::cout << " s" << out << " in shape " << storage_shapes[out] << " " << summarize(Shape::unpack(storage_shapes[out]), storages[out].size(), mark) << " in lane " << storage_positions[out] << ",";
			}
			std::cout << "\n";
			std::cout.flush();
		}

	}

	return 0;
}
