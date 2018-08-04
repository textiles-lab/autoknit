#include "pipeline.hpp"

#include <glm/gtx/norm.hpp>
#include <glm/gtx/hash.hpp>

#include <iostream>
#include <map>
#include <set>
#include <functional>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

//helper to re-label closest points so they have a possible flattening on the bed, given the supplied costs for relabelling:
void flatten(std::vector< uint32_t > &closest, std::vector< float > const &costs, bool is_loop);

void ak::link_chains(
	Parameters const &parameters,
	Model const &slice, //in: slice on which the chains reside
	std::vector< float > const &slice_times, //in: time field (times @ vertices), for slice
	std::vector< std::vector< uint32_t > > const &active_chains, //in: current active chains (slice vertex #'s)
	std::vector< std::vector< Stitch > > const &active_stitches, //in: current active stitches
	std::vector< std::vector< uint32_t > > const &next_chains, //in: current next chains (slice vertex #'s)
	//need this or slice_times (above) std::vector< std::vector< bool > > const &discard_segments,
	std::vector< std::vector< Stitch > > *next_stitches_, //out: next active stitches
	std::vector< Link > *links_ //out: active_chains[from_chain][from_vertex] -> linked_next_chains[to_chain][to_vertex] links
) {
	assert(slice_times.size() == slice.vertices.size());

	for (auto const &chain : active_chains) {
		assert(chain.size() >= 2);
		assert(chain[0] != chain.back() || chain.size() >= 3);
		for (auto v : chain) {
			assert(v < slice.vertices.size());
		}
	}
	assert(active_stitches.size() == active_chains.size());

	for (auto const &stitches : active_stitches) {
		for (auto const &s : stitches) {
			assert(s.t >= 0.0f);
			assert(s.t < 1.0f);
			if (&s > &stitches[0]) {
				assert(s.t > (&s-1)->t);
			}
		}
	}

	for (auto const &chain : next_chains) {
		assert(chain.size() >= 2);
		assert(chain[0] != chain.back() || chain.size() >= 3);
		for (auto v : chain) {
			assert(v < slice.vertices.size());
		}
	}

	assert(next_stitches_);
	auto &next_stitches = *next_stitches_;
	next_stitches.clear();

	assert(links_);
	auto &links = *links_;
	links.clear();

	//figure out the time to trim after:
	float active_max_time = -std::numeric_limits< float >::infinity();
	for (auto const &chain : active_chains) {
		for (auto v : chain) {
			active_max_time = std::max(active_max_time, slice_times[v]);
		}
	}

	//first, write down length to each vertex (makes it easy to move to/from parameter space):
	auto make_lengths = [&slice](std::vector< std::vector< uint32_t > > const &chains) {
		std::vector< std::vector< float > > all_lengths;
		all_lengths.reserve(chains.size());
		for (auto const &chain : chains) {
			all_lengths.emplace_back();
			std::vector< float > &lengths = all_lengths.back();
			lengths.reserve(chain.size());
			float total_length = 0.0f;
			lengths.emplace_back(total_length);
			for (uint32_t vi = 1; vi < chain.size(); ++vi) {
				glm::vec3 const &a = slice.vertices[chain[vi-1]];
				glm::vec3 const &b = slice.vertices[chain[vi]];
				total_length += glm::length(b-a);
				lengths.emplace_back(total_length);
			}
			assert(lengths.size() == chain.size());
		}
		return all_lengths;
	};
	std::vector< std::vector< float > > active_lengths = make_lengths(active_chains);
	std::vector< std::vector< float > > next_lengths = make_lengths(next_chains);

	//figure out ranges (in parameter space) that are to be discarded:
	std::vector< std::vector< std::pair< float, bool > > > next_discard_after;
	next_discard_after.reserve(next_chains.size());
	for (auto const &chain : next_chains) {
		uint32_t ci = &chain - &next_chains[0];
		std::vector< float > const &lengths = next_lengths[ci];
		assert(lengths.size() == chain.size());
		next_discard_after.emplace_back();
		auto &discard_after = next_discard_after.back();
		discard_after.emplace_back(0.0f, (slice_times[chain[0]] > active_max_time));
		for (uint32_t vi = 0; vi + 1 < chain.size(); ++vi) {
			float ta = slice_times[chain[vi]];
			float tb = slice_times[chain[vi+1]];
			float la = lengths[vi];
			float lb = lengths[vi+1];
			//note: for the purposes of assigning ranges, treat time t == active_max_time as t - epsilon.
			if ((ta <= active_max_time && tb > active_max_time)
			 || (ta > active_max_time && tb <= active_max_time) ) {
				float m = (active_max_time - ta) / (tb - ta);
				float l = m * (lb - la) + la;
				discard_after.emplace_back(l / lengths.back(), (tb > active_max_time));
			} else {
				assert((ta > active_max_time) == (tb > active_max_time));
				//do nothing!
			}
		}

		bool is_loop = (chain[0] == chain.back());

		//remove any non-discard segments shorter than 1.5 stitches:
		if (discard_after.size() > 1) {
			float const MinLength = 1.5f * parameters.stitch_width_mm / parameters.model_units_mm;

			float first_len = discard_after[1].first - 0.0f;
			float last_len = 1.0f - discard_after.back().first;
			if (is_loop) {
				first_len = last_len = first_len + last_len;
			}

			for (uint32_t s = 0; s < discard_after.size(); ++s) {
				if (discard_after[s].second != false) continue;
				float len;
				if (s == 0) len = first_len;
				else if (s + 1 == discard_after.size()) len = last_len;
				else len = discard_after[s+1].first - discard_after[s].first;
				if (len * lengths.back() < MinLength) {
					discard_after[s].second = true;
				}
			}
			if (is_loop) {
				assert(discard_after[0].second == discard_after.back().second);
			}
			for (uint32_t s = 0; s + 1 < discard_after.size(); /* later */) {
				if (discard_after[s].second == discard_after[s+1].second) {
					discard_after.erase(discard_after.begin() + (s+1));
				} else {
					++s;
				}
			}
		}

		//remove any discard segments shorter than 0.5 stitches:
		if (discard_after.size() > 1) {
			float const MinLength = 0.5f * parameters.stitch_width_mm / parameters.model_units_mm;

			float first_len = discard_after[1].first - 0.0f;
			float last_len = 1.0f - discard_after.back().first;
			if (is_loop) {
				first_len = last_len = first_len + last_len;
			}

			for (uint32_t s = 0; s < discard_after.size(); ++s) {
				if (discard_after[s].second != true) continue;
				float len;
				if (s == 0) len = first_len;
				else if (s + 1 == discard_after.size()) len = last_len;
				else len = discard_after[s+1].first - discard_after[s].first;
				if (len * lengths.back() < MinLength) {
					discard_after[s].second = false;
				}
			}
			if (is_loop) {
				assert(discard_after[0].second == discard_after.back().second);
			}
			for (uint32_t s = 0; s + 1 < discard_after.size(); /* later */) {
				if (discard_after[s].second == discard_after[s+1].second) {
					discard_after.erase(discard_after.begin() + (s+1));
				} else {
					++s;
				}
			}
		}

	}


	{ //if all segments are marked 'discard', then mark everything 'accept':
		bool only_discard = true;
		for (auto const &discard_after : next_discard_after) {
			for (auto const &d : discard_after) {
				if (!d.second) {
					only_discard = false;
					break;
				}
			}
			if (only_discard == false) break;
		}

		if (only_discard) {
			std::cout << "Marking everything accept because it was all marked discard." << std::endl;
			for (auto &discard_after : next_discard_after) {
				assert(discard_after.size() == 1);
				assert(discard_after[0] == std::make_pair(0.0f, true));
				discard_after[0].second = false;
			}
		} else {
			std::cout << "Have a mix of discard and accept." << std::endl;
		}
	}


	//------------------------------------------------------------------------------------
	//find segments of active and next chains that are mutual nearest neighbors:

	std::vector< std::vector< uint32_t > > active_closest;
	std::vector< std::vector< uint32_t > > next_closest;

	{ //find closest pairs:

		std::vector< std::vector< uint32_t > > adj(slice.vertices.size());
		{ //adjacency matrix -- always handy:
			std::unordered_set< glm::uvec2 > edges;
			for (auto const &tri : slice.triangles) {
				auto do_edge = [&edges](uint32_t a, uint32_t b) {
					if (b > a) std::swap(a,b);
					edges.insert(glm::uvec2(a,b));
				};
				do_edge(tri.x, tri.y);
				do_edge(tri.y, tri.z);
				do_edge(tri.z, tri.x);
			}
			for (auto const &e : edges) {
				adj[e.x].emplace_back(e.y);
				adj[e.y].emplace_back(e.x);
			}
		}

		auto closest_source_chain = [&slice,&adj](
			std::vector< std::vector< uint32_t > > const &sources,
			std::vector< std::vector< uint32_t > > const &targets) {

			std::vector< float > dis(slice.vertices.size(), std::numeric_limits< float >::infinity());
			std::vector< uint32_t > from(slice.vertices.size(), -1U);
			std::vector< std::pair< float, uint32_t > > todo;

			auto queue = [&dis, &todo, &from](uint32_t n, float d, uint32_t f) {
				assert(d < dis[n]);
				dis[n] = d;
				from[n] = f;
				todo.emplace_back(std::make_pair(-d, n));
				std::push_heap(todo.begin(), todo.end());
			};

			for (auto const &chain : sources) {
				uint32_t ci = &chain - &sources[0];
				for (auto const &v : chain) {
					if (v == -1U) continue;
					if (0.0f < dis[v]) queue(v, 0.0f, ci); //some verts appear twice
				}
			}

			while (!todo.empty()) {
				std::pop_heap(todo.begin(), todo.end());
				uint32_t at = todo.back().second;
				float d = -todo.back().first;
				todo.pop_back();
				if (d > dis[at]) continue;
				assert(d == dis[at]);
				for (auto n : adj[at]) {
					float nd = d + glm::length(slice.vertices[n] - slice.vertices[at]);
					if (nd < dis[n]) queue(n, nd, from[at]);
				}
			}

			std::vector< std::vector< uint32_t > > closest;
			closest.reserve(targets.size());
			for (auto const &chain : targets) {
				closest.emplace_back();
				closest.back().reserve(targets.size());
				for (auto v : chain) {
					if (v == -1U) {
						closest.back().emplace_back(-1U);
					} else {
						closest.back().emplace_back(from[v]);
					}
				}
			}
			return closest;
		};

		active_closest = closest_source_chain(next_chains, active_chains);
		next_closest = closest_source_chain(active_chains, next_chains);

		//HACK: treat these as 'segment' values instead of 'vertex' (but really every segment just gets first vertex's value):
		for (auto &c : active_closest) {
			if (!c.empty()) c.pop_back();
		}
		for (auto &c : next_closest) {
			if (!c.empty()) c.pop_back();
		}
#if 0
		//DEBUG:
		for (auto const &ac : active_closest) {
			std::cout << "active_closest[" << (&ac - &active_closest[0]) << "]:";
			for (auto c : ac) {
				std::cout << " " << int32_t(c);
			}
			std::cout << "\n";
		}
		std::cout.flush();
		for (auto const &nc : next_closest) {
			std::cout << "next_closest[" << (&nc - &next_closest[0]) << "]:";
			for (auto c : nc) {
				std::cout << " " << int32_t(c);
			}
			std::cout << "\n";
		}
		std::cout.flush();
#endif
	}

	struct BeginEnd {
		BeginEnd(float begin_, float end_) : begin(begin_), end(end_) { }
		float begin;
		float end;
	};

	//sort active and next into matching parametric ranges:
	struct Match {
		std::vector< BeginEnd > active; //vertex ranges
		std::vector< BeginEnd > next; //segment ranges
	};

	std::map< std::pair< uint32_t, uint32_t >, Match > matches;


	for (auto &closest : active_closest) {
		uint32_t ai = &closest - &active_closest[0];
		auto const &lengths = active_lengths[ai];
		assert(lengths.size() == closest.size() + 1);

		std::vector< float > weights; weights.reserve(closest.size());
		for (uint32_t i = 1; i < lengths.size(); ++i) {
			weights.emplace_back(lengths[i] - lengths[i-1]);
		}
		flatten(closest, weights, active_chains[ai].empty() || active_chains[ai][0] == active_chains[ai].back());

		for (uint32_t begin = 0; begin < closest.size(); /* later */) {
			uint32_t end = begin + 1;
			while (end < closest.size() && closest[end] == closest[begin]) ++end;
			assert(end < lengths.size());
			matches[std::make_pair(ai, closest[begin])].active.emplace_back(lengths[begin] / lengths.back(), lengths[end] / lengths.back());
			begin = end;
		}
	}

	for (auto &closest : next_closest) {
		uint32_t ni = &closest - &next_closest[0];

		assert(ni < next_lengths.size());
		auto const &lengths = next_lengths[ni];
		assert(lengths.size() == closest.size() + 1);

		std::vector< float > weights; weights.reserve(closest.size());
		for (uint32_t i = 1; i < lengths.size(); ++i) {
			weights.emplace_back(lengths[i] - lengths[i-1]);
		}
		flatten(closest, weights, next_chains[ni].empty() || next_chains[ni][0] == next_chains[ni].back());


		for (uint32_t begin = 0; begin < closest.size(); /* later */) {
			uint32_t end = begin + 1;
			while (end < closest.size() && closest[end] == closest[begin]) ++end;
			assert(end < lengths.size());
			matches[std::make_pair(closest[begin], ni)].next.emplace_back(lengths[begin] / lengths.back(), lengths[end] / lengths.back());
			begin = end;
		}
	}

	std::vector< uint32_t > active_matches(active_chains.size(), 0);
	std::vector< uint32_t > next_matches(next_chains.size(), 0);
	uint32_t empty_matches = 0;
	for (auto const &anm : matches) {
		if (anm.first.first != -1U) active_matches[anm.first.first] += 1;
		if (anm.first.second != -1U) next_matches[anm.first.second] += 1;
		if (anm.first.first == -1U || anm.first.second == -1U) ++empty_matches;
	}
	std::cout << "NOTE: have " << empty_matches << " segments that match with nothing." << std::endl;

	{ //If there are any merges or splits, all participating next cycles are marked 'accept':
		std::set< uint32_t > to_mark;
		for (auto const &anm : matches) {
			if (anm.first.first == -1U || anm.first.second == -1U) continue;
			bool is_split_or_merge = (active_matches[anm.first.first] > 1 || next_matches[anm.first.second] > 1);
			if (is_split_or_merge) {
				to_mark.insert(anm.first.second);
			}
		}
		uint32_t were_marked = 0;
		for (auto ni : to_mark) {
			for (auto const &td : next_discard_after[ni]) {
				if (td.second) ++were_marked;
			}
			next_discard_after[ni].assign(1, std::make_pair(0.0f, false));
		}
		if (!to_mark.empty() && were_marked != 0) {
			std::cerr << "Marked " << were_marked << " segments on " << to_mark.size() << " next cycles as 'accept' based on participating in a merge/split." << std::endl;
		}
	}


	next_stitches.assign(next_chains.size(), std::vector< ak::Stitch >());

	//allocate stitch counts based on segment lengths + source stitches:
	//(and limit based on active stitch counts)
	//TODO: account for balance also, somehow!
	//then make next stitches
	for (auto const &anm : matches) {
		Match const &match = anm.second;

		if (match.active.empty()) {
			std::cout << "Ignoring match with empty active chain." << std::endl;
			continue;
		} else if (match.next.empty()) {
			if (active_matches[anm.first.first] == 1) {
				std::cout << "WARNING: active chain matches nothing at all; will not be linked and will thus be discarded." << std::endl;
			} else {
				std::cout << "Ignoring match with empty next chain." << std::endl;
			}

			continue;
		}
		assert(!match.active.empty());
		assert(!match.next.empty());

		bool is_split_or_merge = (active_matches[anm.first.first] > 1 || next_matches[anm.first.second] > 1);
		//if is_split_or_merge will only link 1-1

		//compute min/max totals from active stitches:
		uint32_t active_ones = 0;
		uint32_t active_anys = 0;
		{
			assert(anm.first.first < active_stitches.size());
			auto const &stitches = active_stitches[anm.first.first];
			auto si = stitches.begin();
			for (auto const &be : match.active) {
				while (si != stitches.end() && (si->t < be.begin)) ++si;
				while (si != stitches.end() && (si->t < be.end)) {
					if      (si->flag == Stitch::FlagLinkOne) ++active_ones;
					else if (si->flag == Stitch::FlagLinkAny) ++active_anys;
					else assert(si->flag == Stitch::FlagLinkOne || si->flag == Stitch::FlagLinkAny);
					++si;
				}
			}
		}

		//compute number of ones based on discard segments:
		//Ideally we'd like this to be true (saves on discarded stitches):
		//  discarded segments contain at least one stitch marked 'LinkOne'
		//  while kept segments contain at least two stitches marked 'LinkOne'
		//What we're doing below actually puts two stitches in discarded segments as well.
		uint32_t next_ones = 0;

		bool next_is_loop = (next_chains[anm.first.second][0] == next_chains[anm.first.second].back());
		auto const &discard_after = next_discard_after[anm.first.second];
		{ //compute next_ones:
			if (next_is_loop) {
				assert(discard_after[0].second == discard_after.back().second);
			}
			//every discard/non-discard edge needs at least one stitch next to it on each side:
			for (auto const &be : match.next) {
				for (auto tdi = discard_after.begin(); tdi != discard_after.end(); ++tdi) {
					if (next_is_loop && tdi == discard_after.begin()) continue;
					if (be.begin <= tdi->first && tdi->first < be.end) ++next_ones;
					if (be.begin < tdi->first && tdi->first <= be.end) ++next_ones;
				}
			}

			if (next_ones > 2 * active_anys + active_ones) {
				std::cerr << "ERROR: more discard/non-discard ends are required (" << next_ones << ") than are permitted by the current active flags (" << active_anys << "*2 + " << active_ones << "); code to fix this (by removing shortest same-discard segment) not yet implemented." << std::endl;
				assert(next_ones <= 2 * active_anys + active_ones);
			}
		}

		//compute desired stitch count based on segment lengths:
		assert(anm.first.second < next_lengths.size());
		std::vector< float > const &lengths = next_lengths[anm.first.second];
		float total_length = 0.0f;
		std::cout << "Matching to"; //DEBUG
		for (auto const &be : match.next) {
			std::cout << " [" << be.begin << ", " << be.end << ")"; //DEBUG
			assert(be.begin < be.end);
			total_length += be.end - be.begin;
		}
		std::cout << std::endl; //DEBUG
		total_length *= lengths.back();

		float stitch_width = parameters.stitch_width_mm / parameters.model_units_mm;
		uint32_t stitches = std::max(1, int32_t(std::round(total_length / stitch_width)));

		{ //adjust for possible links:
			//least is to link 1-1 for every next_ones and then link everything else 2-1:
			uint32_t lower = next_ones //next ones to link to one stitch
				+ (std::max(0, int32_t(active_ones + active_anys) - int32_t(next_ones)) + 1) / 2 //other stitches link 2-1
			;
			uint32_t upper = active_ones + 2 * active_anys; //most is to increase from all anys
			assert(lower <= upper);

			if (is_split_or_merge) {
				assert(lower <= active_ones + active_anys && active_ones + active_anys <= upper);
				std::cout << "NOTE: setting stitches from " << stitches << " to ";
				stitches = active_ones + active_anys;
				std::cout << stitches << " to make split/merge 1-1." << std::endl;
				//stitches = active_ones + active_anys;
			}

			if (stitches < lower || stitches > upper) {
				std::cout << "NOTE: stitches (" << stitches << ") will be clamped to possible range [" << lower << ", " << upper << "], which might cause some shape distortion." << std::endl;
				stitches = std::max(lower, std::min(upper, stitches));
			}
			std::cout << "Will make " << stitches << " stitches, given active with " << active_ones << " ones, " << active_anys << " anys; next with " << next_ones << " ones." << std::endl; //DEBUG
		}

		std::vector< ak::Stitch > new_stitches;
		if (stitches > 0) {
			//spread stitches among "allocation ranges" (same discard status)
			struct Alloc {
				Alloc(float begin_, float end_, bool first_one_, bool last_one_) : begin(begin_), end(end_), first_one(first_one_), last_one(last_one_) { }
				float begin, end;
				bool first_one, last_one;
				uint32_t stitches = 0;
				float length = std::numeric_limits< float >::quiet_NaN();
			};
			std::vector< Alloc > alloc;
			for (auto const &be : match.next) {
				alloc.emplace_back(be.begin, be.end, false, false);
				//split allocation range on discards:
				for (auto tdi = discard_after.begin(); tdi != discard_after.end(); ++tdi) {
					if (next_is_loop && tdi == discard_after.begin()) continue;
					if (tdi->first < alloc.back().begin) {
					} else if (tdi->first == alloc.back().begin) {
						alloc.back().first_one = true;
					} else if (tdi->first < alloc.back().end) {
						alloc.back().last_one = true;
						alloc.back().end = tdi->first;
						alloc.emplace_back(tdi->first, be.end, true, false);
					} else if (tdi->first == alloc.back().end) {
						alloc.back().last_one = true;
					} else { assert(tdi->first > alloc.back().end);
					}
				}
			}
			uint32_t total_ones = 0;
			for (auto &a : alloc) {
				a.stitches = (a.first_one ? 1 : 0) + (a.last_one ? 1 : 0);
				total_ones += a.stitches;
				a.length = (a.end - a.begin) * lengths.back();
			}

			assert(total_ones == next_ones); //better have the same number of ones as we accounted for previously

			//add remaining stitches to alloc ranges based on which range has the least-dense stitches:
			for (uint32_t s = total_ones; s < stitches; ++s) {
				uint32_t best = -1U;
				float best_density = 0.0f;
				for (auto const &a : alloc) {
					float d = a.length / float(a.stitches + 1);
					if (d > best_density) {
						best = &a - &alloc[0];
						best_density = d;
					}
				}
				assert(best < alloc.size());
				alloc[best].stitches += 1;
			}

			//actually create stitches from allocation ranges:
			for (auto const &a : alloc) {
				if (a.stitches == 0) continue;
				for (uint32_t s = 0; s < a.stitches; ++s) {
					float t = (s + 0.5f) / float(a.stitches) * (a.end - a.begin) + a.begin;
					ak::Stitch::Flag flag = ak::Stitch::FlagLinkAny;
					if (s == 0 && a.first_one) flag = ak::Stitch::FlagLinkOne;
					if (s + 1 == a.stitches && a.last_one) flag = ak::Stitch::FlagLinkOne;
					new_stitches.emplace_back( t, flag );

				}
			}
		}
		assert(new_stitches.size() == stitches);

		next_stitches[anm.first.second].insert(next_stitches[anm.first.second].end(), new_stitches.begin(), new_stitches.end());
	} //end stitch allocation

	for (auto &stitches : next_stitches) {
		std::stable_sort(stitches.begin(), stitches.end(), [](ak::Stitch const &a, ak::Stitch const &b){
			return a.t < b.t;
		});
	}

	for (auto const &anm : matches) {
		Match const &match = anm.second;

		if (match.active.empty()) {
			std::cout << "Ignoring match with empty active chain." << std::endl;
			continue;
		} else if (match.next.empty()) {
			std::cout << "Ignoring match with empty next chain." << std::endl;
			continue;
		}
		assert(!match.active.empty());
		assert(!match.next.empty());

		//read back all stitches in a particular set of ranges:
		auto extract_stitch_info = [&slice](
			std::vector< uint32_t > const &chain,
			std::vector< float > const &lengths,
			std::vector< ak::Stitch > const &stitches,
			std::vector< BeginEnd > const &ranges,
			std::vector< uint32_t > *stitch_indices_,
			std::vector< glm::vec3 > *stitch_locations_,
			std::vector< bool > *stitch_linkones_ ) {

			assert(chain.size() == lengths.size());

			assert(stitch_indices_);
			auto &stitch_indices = *stitch_indices_;
			stitch_indices.clear();

			assert(stitch_locations_);
			auto &stitch_locations = *stitch_locations_;
			stitch_locations.clear();

			assert(stitch_linkones_);
			auto &stitch_linkones = *stitch_linkones_;
			stitch_linkones.clear();

			auto li = lengths.begin();
			auto si = stitches.begin();
			for (auto const &be : ranges) {
				while (si != stitches.end() && si->t < be.begin) ++si;
				while (si != stitches.end() && si->t <= be.end) {
					float l = si->t * lengths.back();
					while (li != lengths.end() && *li <= l) ++li;
					assert(li != lengths.end());
					assert(li != lengths.begin());
					uint32_t i = li - lengths.begin();

					stitch_indices.emplace_back(si - stitches.begin());
					stitch_locations.emplace_back(
						glm::mix(
							slice.vertices[chain[i-1]],
							slice.vertices[chain[i]],
							(l - *(li-1)) / (*li - *(li-1))
						)
					);
					stitch_linkones.emplace_back( si->flag == ak::Stitch::FlagLinkOne );
					++si;
				}
			}
		};

		std::vector< uint32_t > next_stitch_indices;
		std::vector< glm::vec3 > next_stitch_locations;
		std::vector< bool > next_stitch_linkones;
		extract_stitch_info(
			next_chains[anm.first.second],
			next_lengths[anm.first.second],
			next_stitches[anm.first.second],
			match.next,
			&next_stitch_indices,
			&next_stitch_locations,
			&next_stitch_linkones );

		std::vector< uint32_t > active_stitch_indices;
		std::vector< glm::vec3 > active_stitch_locations;
		std::vector< bool > active_stitch_linkones;
		extract_stitch_info(
			active_chains[anm.first.first],
			active_lengths[anm.first.first],
			active_stitches[anm.first.first],
			match.active,
			&active_stitch_indices,
			&active_stitch_locations,
			&active_stitch_linkones );

		//figure out how to link to next stitches:
		{ //DEBUG:
			uint32_t new_ones = 0;
			for (auto o : next_stitch_linkones) {
				if (o) ++new_ones;
			}
			uint32_t active_ones = 0;
			for (auto o : active_stitch_linkones) {
				if (o) ++active_ones;
			}
			std::cout << " About to connect " << active_stitch_locations.size() << " active stitches (" << active_ones << " linkones) to " << next_stitch_locations.size() << " new stitches (" << new_ones << " linkones)." << std::endl;
		}

		//actually build links:
		{ //least-clever linking solution: FlagLinkOne's link 1-1, others link to keep arrays mostly in sync

			std::vector< std::pair< uint32_t, uint32_t > > best_links;
			float best_cost = std::numeric_limits< float >::infinity();

			auto try_links_even = [&](uint32_t roll_active, uint32_t roll_new) {
				std::vector< std::pair< uint32_t, uint32_t > > possible_links;

				if (active_stitch_locations.size() <= next_stitch_locations.size()) {
					//evenly distribute increases among the non-linkone stitches:
					uint32_t total = 0;
					for (auto l : active_stitch_linkones) {
						if (!l) ++total;
					}
					uint32_t increases = next_stitch_locations.size() - active_stitch_locations.size();
					std::vector< bool > inc(total, false);
					for (uint32_t i = 0; i < increases; ++i) {
						assert(inc[i * total / increases] == false);
						inc[i * total / increases] = true;
					}
					uint32_t n = 0;
					uint32_t i = 0;
					for (uint32_t a = 0; a < active_stitch_locations.size(); ++a) {
						uint32_t ra = (a + roll_active) % active_stitch_locations.size();
						uint32_t rn = (n + roll_new) % next_stitch_locations.size();
						possible_links.emplace_back(ra, rn);
						++n;
						if (!active_stitch_linkones[ra]) {
							assert(i < inc.size());
							if (inc[i]) {
								rn = (n + roll_new) % next_stitch_locations.size();
								possible_links.emplace_back(ra, rn);
								++n;
							}
							++i;
						}
					}
					assert(i == total);
					assert(n == next_stitch_locations.size());
				} else if (active_stitch_locations.size() > next_stitch_locations.size()) {
					//evenly distribute decreases among the non-linkone stitches:
					uint32_t total = 0;
					for (auto l : next_stitch_linkones) {
						if (!l) ++total;
					}
					uint32_t decreases = active_stitch_locations.size() - next_stitch_locations.size();
					std::vector< bool > dec(total, false);
					for (uint32_t i = 0; i < decreases; ++i) {
						assert(dec[i * total / decreases] == false);
						dec[i * total / decreases] = true;
					}
					uint32_t a = 0;
					uint32_t i = 0;
					for (uint32_t n = 0; n < next_stitch_locations.size(); ++n) {
						uint32_t ra = (a + roll_active) % active_stitch_locations.size();
						uint32_t rn = (n + roll_new) % next_stitch_locations.size();
						possible_links.emplace_back(ra, rn);
						++a;
						if (!next_stitch_linkones[rn]) {
							assert(i < dec.size());
							if (dec[i]) {
								ra = (a + roll_active) % active_stitch_locations.size();
								possible_links.emplace_back(ra, rn);
								++a;
							}
							++i;
						}
					}
					assert(i == total);
					assert(a == active_stitch_locations.size());
				}
				float const row_height = 2.0f * parameters.stitch_height_mm / parameters.model_units_mm;
				float cost = 0.0f;
				for (auto const &p : possible_links) {
					float len = glm::length(next_stitch_locations[p.second] - active_stitch_locations[p.first]);
					cost += (len - row_height) * (len - row_height);
				}

				if (cost < best_cost) {
					best_cost = cost;
					best_links = possible_links;
				}
			};


			if (active_stitch_locations.size() >= next_stitch_locations.size()) {
				for (uint32_t roll_active = 0; roll_active < active_stitch_locations.size(); ++roll_active) {
					try_links_even(roll_active,0);
				}
			} else {
				for (uint32_t roll_new = 0; roll_new < next_stitch_locations.size(); ++roll_new) {
					try_links_even(0,roll_new);
				}
			}


			for (auto const &p : best_links) {
				Link link;
				link.from_chain = anm.first.first;
				assert(p.first < active_stitch_indices.size());
				link.from_stitch = active_stitch_indices[p.first];

				link.to_chain = anm.first.second;
				assert(p.second < next_stitch_indices.size());
				link.to_stitch = next_stitch_indices[p.second];
				links.emplace_back(link);
			}
		}
	}

	//mark next stitches in discard range as 'discard':
	assert(next_discard_after.size() == next_stitches.size());
	uint32_t marked = 0;
	uint32_t total = 0;
	for (auto &stitches : next_stitches) {
		auto const &discard_after = next_discard_after[&stitches - &next_stitches[0]];

		auto di = discard_after.begin();
		for (auto &s : stitches) {
			assert(di != discard_after.end());
			while (di + 1 != discard_after.end() && (di + 1)->first <= s.t) ++di;
			assert(di->first <= s.t);
			if (di->second) {
				s.flag = ak::Stitch::FlagDiscard;
				++marked;
			}
			++total;
		}
	}
	std::cout << "Marked " << marked << " of " << total << " newly created stitches as 'discard'." << std::endl;


}


//-------------------------------------------------------------------

void flatten(std::vector< uint32_t > &closest, std::vector< float > const &weights, bool is_loop) {
	assert(closest.size() == weights.size());
	if (closest.empty()) return;

	//make sure that 'closest' looks like:
	//  a a a b b b c c c
	//   a a b b b b c c
	// that is, can be flattened to the knitting machine bed while preserving constituent cycles
	// One view of this: if you start at some stitch A, then the left side should be a mirror of the right side
	//  (with the possible exception that some symbols may be skipped on the left or right)

	//(a) condense closest into short list of symbols:
	std::vector< std::pair< uint32_t, float > > symbols;
	symbols.reserve(closest.size()); //certainly no more symbols than closest
	for (uint32_t i = 0; i < closest.size(); ++i) {
		uint32_t symb = closest[i];
		if (symbols.empty() || symbols.back().first != symb) {
			symbols.emplace_back(std::make_pair(symb, 0.0f));
		}
		assert(symbols.back().first == symb);
		symbols.back().second += weights[i];
	}
	assert(!symbols.empty());

	//(b) early-out in certain easy-to-check conditions:
	if (symbols.size() == 1) return; //single symbol
	//DEBUG: don't do these checks; exercise the code a bit more instead:
	//if (symbols.size() == 2) return; //two symbols without alternation
	//if (is_loop && symbols.size() == 3 && symbols[0] == symbols.back()) return; //two symbols without alternation (loop version)

	//symbols -> bits
	std::vector< std::pair< uint16_t, float > > bit_symbols;
	uint32_t bits;
	{
		bit_symbols.reserve(symbols.size());
		std::map< uint32_t, uint16_t > symbol_bit;
		for (auto &sw : symbols) {
			auto ret = symbol_bit.insert(std::make_pair(sw.first, uint16_t(1 << symbol_bit.size())));
			assert(ret.first->second && "Only have enough bits for 16-symbol flattening");
			bit_symbols.emplace_back(ret.first->second, sw.second);
		}
		assert(bit_symbols.size() == symbols.size());
		bits = symbol_bit.size();
	}

	struct State {
		uint16_t used = 0; //have seen all symbols with bits in used
		uint8_t min = 0; //symbols that are strictly between min and max have been processed
		uint8_t max = 0;
		uint16_t current = 0; //most recently kept symbol
		uint16_t padding = 0; //padding to make state 64 bits long

		typedef uint64_t Packed;
		Packed pack() const {
			return *reinterpret_cast< Packed const * >(this);
		}
		static State unpack(Packed packed) {
			return *reinterpret_cast< State const * >(&packed);
		}

		std::string to_string(uint32_t bits = 16) const {
			std::string ret = "(" + std::to_string(min) + "," + std::to_string(max) + ")";
			for (uint32_t i = bits-1; i < bits; --i) {
				if (used & (1 << i)) {
					if ((1 << i) == current) {
						ret += "*";
					} else {
						ret += "x";
					}
				} else {
					ret += ".";
				}
			}
			return ret;
		}
	};
	static_assert(sizeof(State) == 8, "packed state");


	struct {
		State::Packed state;
		float cost = std::numeric_limits< float >::infinity();
		State::Packed from;
	} finished;
	std::unordered_map< State::Packed, std::pair< float, State::Packed > > visited;
	std::vector< std::pair< float, State::Packed > > todo;
	static std::greater< std::pair< float, State::Packed > > const TODOCompare;

	auto queue_state = [&visited, &finished, &todo, bits](State const state, float const cost, State const from) {
		assert(state.min != from.min || state.max != from.max); //must have done *something*

		std::cout << state.to_string(bits) << " from " << from.to_string(bits) << " cost " << cost << std::endl; //DEBUG

		if ((state.min != from.min && state.min == from.max)
		 || (state.max != from.max && state.max == from.min)) {
			//pointers crossed or met -> state is finished!
			assert(state.min == state.max || (state.min == from.max && state.max == from.min));
			if (cost < finished.cost) {
				finished.state = state.pack();
				finished.cost = cost;
				finished.from = from.pack();
			}
			return;
		}

		//queue/indicate regular
		auto ret = visited.insert(std::make_pair(state.pack(), std::make_pair(cost, from.pack())));
		if (ret.second || ret.first->second.first > cost) {
			ret.first->second = std::make_pair(cost, from.pack());
			todo.emplace_back(std::make_pair(cost, state.pack()));
			std::push_heap(todo.begin(), todo.end(), TODOCompare);
		}

	};

	auto expand_state = [&bit_symbols, &queue_state](State const state, float const cost) {
		auto const used = state.used;
		auto const min = state.min;
		auto const max = state.max;
		auto const current = state.current;
		assert(state.padding == 0);

		//*_next_symbol is the symbol that is advanced over when moving min/max,
		// leading to some asymmetry in indexing:
		// a(bc)d -> (abc)d <-- min_next_symbol is 'a' (at index of min_next)
		// a(bc)d -> a(bcd) <-- max_next_symbol is 'd' (at index of max)

		auto min_next = (min == 0 ? bit_symbols.size() - 1 : min - 1);
		auto min_next_symbol = bit_symbols[min_next];
		auto max_next = (max + 1U < bit_symbols.size() ? max + 1 : 0);
		auto max_next_symbol = bit_symbols[max];

		//actions:
		//no reason not to keep if symbol is current:
		if (min_next_symbol.first == current || max_next_symbol.first == current) {
			assert(used & current);
			State next;
			next.used = used;
			next.min = (min_next_symbol.first == current ? min_next : min);
			next.max = (max_next_symbol.first == current ? max_next : max);
			next.current = current;

			float next_cost = cost;

			queue_state(next, next_cost, state);

			return; //no other actions worth taking; this one was free!
		}

		//keep min (symbol must be unused):
		if (!(used & min_next_symbol.first)) {
			State next;
			next.used = used | min_next_symbol.first;
			next.min = min_next;
			next.max = max;
			next.current = min_next_symbol.first;

			float next_cost = cost;

			queue_state(next, next_cost, state);
		}

		//keep max (symbol must be unused):
		if (!(used & max_next_symbol.first)) {
			State next;
			next.used = used | max_next_symbol.first;
			next.min = min;
			next.max = max_next;
			next.current = max_next_symbol.first;

			float next_cost = cost;

			queue_state(next, next_cost, state);
		}

		//discard min:
		{
			State next;
			next.used = used;
			next.min = min_next;
			next.max = max;
			next.current = current;

			float next_cost = cost + min_next_symbol.second;

			queue_state(next, next_cost, state);
		}

		//discard max:
		{
			State next;
			next.used = used;
			next.min = min;
			next.max = max_next;
			next.current = current;

			float next_cost = cost + max_next_symbol.second;

			queue_state(next, next_cost, state);
		}
	};

	//queue starting states:
	for (uint32_t s = 0; s < bit_symbols.size(); ++s) {
		State init;
		init.used = 0;
		init.min = s;
		init.max = s;
		init.current = 0;
		expand_state(init, 0.0f);
		if (!is_loop) break;
	}

	while (!todo.empty()) {
		std::pop_heap(todo.begin(), todo.end(), TODOCompare);
		auto state = State::unpack(todo.back().second);
		float cost = todo.back().first;
		todo.pop_back();
		//if the cheapest thing is more expensive than the finish, we're done:
		if (cost >= finished.cost) break;

		{ //cost should either be stale or what is stored in 'visited':
			auto f = visited.find(state.pack());
			assert(f != visited.end());
			if (cost > f->second.first) continue;
			assert(cost == f->second.first);
		}
		expand_state(state, cost);
	}
	assert(finished.cost != std::numeric_limits< float >::infinity()); //found ~some~ path

	//read back states:
	std::vector< State::Packed > path;
	path.emplace_back(finished.state);
	path.emplace_back(finished.from);
	while (true) {
		auto f = visited.find(path.back());
		if (f == visited.end()) break;
		path.emplace_back(f->second.second);
	}
	std::reverse(path.begin(), path.end());
	std::cout << "----" << std::endl; //DEBUG

	std::vector< int8_t > keep(bit_symbols.size(), -1);
	for (uint32_t i = 1; i < path.size(); ++i) {
		State state = State::unpack(path[i-1]);
		State next = State::unpack(path[i]);

		std::cout << state.to_string(bits) << " -> " << next.to_string(bits) << ": "; std::cout.flush(); //DEBUG
		if (state.min != next.min && state.max != next.max) {
			//a(bc)d -> (abcd), keep 'a' (next.min), 'd' (state.max)
			assert(bit_symbols[next.min].first == bit_symbols[state.max].first);
			assert(next.current = bit_symbols[state.min].first);
			std::cout << "keep " << int32_t(next.min) << " (\"" << symbols[next.min].first << "\")" << ", " << int32_t(state.max) << " (\"" << symbols[state.max].first << "\")" << std::endl; //DEBUG
			assert(keep[next.min] == -1);
			assert(keep[state.max] == -1);
			keep[next.min] = keep[state.max] = 1;
		} else if (state.min != next.min) {
			//a(bc)d -> (abc)d, keep/discard next.min
			if (bit_symbols[next.min].first == next.current) {
				std::cout << "keep " << int32_t(next.min) << " (\"" << symbols[next.min].first << "\")" << std::endl; //DEBUG
				assert(keep[next.min] == -1);
				keep[next.min] = 1;
			} else {
				std::cout << "discard " << int32_t(next.min) << " (\"" << symbols[next.min].first << "\")" << std::endl; //DEBUG
				assert(keep[next.min] == -1);
				keep[next.min] = 0;
			}
		} else { assert(state.max != next.max);
			//a(bc)d -> a(bcd), keep/discard state.max
			if (bit_symbols[state.max].first == next.current) {
				std::cout << "keep " << int32_t(state.max) << " (\"" << symbols[state.max].first << "\")" << std::endl; //DEBUG
				assert(keep[state.max] == -1);
				keep[state.max] = 1;
			} else {
				std::cout << "discard " << int32_t(state.max) << " (\"" << symbols[state.max].first << "\")" << std::endl; //DEBUG
				assert(keep[state.max] == -1);
				keep[state.max] = 0;
			}
		}
	}

	//DEBUG: was at least one thing kept?
	bool kept_at_least_one = false;
	for (auto k : keep) {
		assert(k != -1);
		if (k == 1) kept_at_least_one = true;
	}
	assert(kept_at_least_one);

	//use keep to figure out which elements of closest to re-label.
	std::vector< bool > relabel; relabel.reserve(closest.size());
	{
		std::cout << "relabel:"; //DEBUG
		auto si = symbols.begin();
		for (auto c : closest) {
			assert(si != symbols.end());
			if (si->first != c) ++si;
			assert(si != symbols.end());
			assert(si->first == c);
			relabel.emplace_back(keep[si - symbols.begin()] == 0);
			if (relabel.back()) {
				std::cout << ' ' << 'x' << int32_t(c) << 'x'; //DEBUG
			} else {
				std::cout << ' ' << ' ' << int32_t(c) << ' '; //DEBUG
			}
		}
		std::cout << std::endl; //DEBUG
		assert(relabel.size() == closest.size());
		assert(si != symbols.end());
		++si;
		assert(si == symbols.end());

		bool have_keep = false;
		for (uint32_t seed = 0; seed < closest.size(); ++seed) {
			if (relabel[seed] == false) have_keep = true;
		}
		assert(have_keep);
	}

	{ //do relabelling:
		auto relabel_range = [&closest,&weights,&relabel](uint32_t first, uint32_t last) {
			uint32_t before = (first == 0 ? closest.back() : closest[first-1]);
			uint32_t after  = (last + 1 == closest.size() ? closest[0] : closest[last+1]);
			std::cout << "Relabelling [" << first << ", " << last << "] using " << before << "/" << after << std::endl; //DEBUG

			assert(!relabel[(first == 0 ? closest.size() : first) - 1]);
			assert(!relabel[(last + 1 == closest.size() ? 0 : last) + 1]);
			assert(relabel[first]);
			assert(relabel[last]);

			float total = 0.0f;
			{ //first pass: figure out total weight sum:
				uint32_t i = first;
				while (true) {
					assert(relabel[i]);
					total += weights[i];
					if (i == last) break;
					++i;
				}
			}
			float sum = 0.0f;
			{ //second pass: assign before/after based on weight sum:
				uint32_t i = first;
				while (true) {
					if (sum + 0.5f * weights[i] < 0.5f * total) {
						closest[i] = before;
					} else {
						closest[i] = after;
					}
					relabel[i] = false;
					sum += weights[i];
					if (i == last) break;
					++i;
				}
			}
		};
		for (uint32_t seed = 0; seed < closest.size(); ++seed) {
			if (!relabel[seed]) continue;
			uint32_t first = seed;
			uint32_t last = seed;
			while (relabel[last + 1 < closest.size() ? last + 1 : 0]) {
				last = (last + 1 < closest.size() ? last + 1 : 0);
			}
			relabel_range(first, last);
		}
	}

}

