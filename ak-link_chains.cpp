#include "pipeline.hpp"

#include <glm/gtx/norm.hpp>
#include <glm/gtx/hash.hpp>

#include <iostream>
#include <map>
#include <set>
#include <functional>
#include <algorithm>
#include <unordered_set>

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

	for (auto const &closest : active_closest) {
		//TODO: make sure nothing appears more than twice in closest
		uint32_t ai = &closest - &active_closest[0];

		auto const &lengths = active_lengths[ai];
		assert(lengths.size() == closest.size() + 1);

		for (uint32_t begin = 0; begin < closest.size(); /* later */) {
			uint32_t end = begin + 1;
			while (end < closest.size() && closest[end] == closest[begin]) ++end;
			assert(end < lengths.size());
			matches[std::make_pair(ai, closest[begin])].active.emplace_back(lengths[begin] / lengths.back(), lengths[end] / lengths.back());
			begin = end;
		}
	}

	for (auto const &closest : next_closest) {
		//TODO: make sure nothing appears more than twice in closest
		uint32_t ni = &closest - &next_closest[0];

		assert(ni < next_lengths.size());
		auto const &lengths = next_lengths[ni];
		assert(lengths.size() == closest.size() + 1);

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
			bool is_split_or_merge = (active_matches[anm.first.first] > 1 || active_matches[anm.first.second] > 1);
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
			std::cout << "Ignoring match with empty next chain." << std::endl;
			continue;
		}
		assert(!match.active.empty());
		assert(!match.next.empty());

		//bool is_split_or_merge = (active_matches[anm.first.first] > 1 || active_matches[anm.first.second] > 1);

		//if (is_split_or_merge) continue; //TODO: handle splits/merges! somehow!

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
			active_chains[anm.first.second],
			active_lengths[anm.first.second],
			active_stitches[anm.first.second],
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

