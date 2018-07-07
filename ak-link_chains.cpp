#include "pipeline.hpp"

#include <glm/gtx/norm.hpp>

#include <iostream>
#include <map>
#include <set>
#include <functional>
#include <algorithm>

void ak::link_chains(
	ak::Parameters const &parameters,
	ak::Model const &model, //in: model
	std::vector< float > const &times,          //in: time field (times @ vertices)
	std::vector< std::vector< ak::EmbeddedVertex > > const &active_chains, //in: current active chains
	std::vector< std::vector< ak::Flag > > const &active_flags, //in: stitches
	std::vector< std::vector< ak::EmbeddedVertex > > const &next_chains_in, //in: next chains
	std::vector< std::vector< ak::EmbeddedVertex > > *linked_next_chains_, //out: next chains
	std::vector< std::vector< ak::Flag > > *linked_next_flags_, //out: flags indicating status of vertices on next chains
	std::vector< ak::Link > *links_ //out: active_chains[from_chain][from_vertex] -> linked_next_chains[to_chain][to_vertex] links
) {
	assert(times.size() == model.vertices.size());

	for (auto const &chain : active_chains) {
		assert(chain.size() >= 2);
		assert(chain[0] != chain.back() || chain.size() >= 3);
	}

	for (auto const &chain : next_chains_in) {
		assert(chain.size() >= 2);
		assert(chain[0] != chain.back() || chain.size() >= 3);
	}

	assert(linked_next_chains_);
	auto &linked_next_chains = *linked_next_chains_;
	linked_next_chains.clear();

	assert(linked_next_flags_);
	auto &linked_next_flags = *linked_next_flags_;
	linked_next_flags.clear();

	assert(links_);
	auto &links = *links_;
	links.clear();

	//figure out the time to trim after:
	float active_max_time = -std::numeric_limits< float >::infinity();
	for (auto const &chain : active_chains) {
		for (auto const &ev : chain) {
			active_max_time = std::max(active_max_time, ev.interpolate(times));
		}
	}

	//make a copy of next_chains_in with vertices inserted at time field crossings:
	std::vector< std::vector< ak::EmbeddedVertex > > next_chains;
	std::vector< std::vector< bool > > next_segment_discards;
	for (auto const &chain_in : next_chains_in) {
		next_chains.emplace_back();
		std::vector< ak::EmbeddedVertex > &chain = next_chains.back();
		next_segment_discards.emplace_back();
		std::vector< bool > &discards = next_segment_discards.back();

		chain.emplace_back(chain_in[0]);
		for (uint32_t ci = 0; ci + 1 < chain_in.size(); ++ci) {
			float a = chain_in[ci].interpolate(times);
			float b = chain_in[ci+1].interpolate(times);
			//NOTE: treat exactly active_max_time as active_max_time + epsilon
			if (a < active_max_time && b >= active_max_time) {
				float m = (active_max_time - a) / (b - a);
				discards.emplace_back(false);
				chain.emplace_back(EmbeddedVertex::mix(chain_in[ci], chain_in[ci+1], m));
				discards.emplace_back(true);
			} else if (b < active_max_time && a >= active_max_time) {
				float m = (active_max_time - a) / (b - a);
				discards.emplace_back(true);
				chain.emplace_back(EmbeddedVertex::mix(chain_in[ci], chain_in[ci+1], m));
				discards.emplace_back(false);
			} else {
				assert((a >= active_max_time) == (b >= active_max_time));
				discards.emplace_back(a >= active_max_time);
			}
			chain.emplace_back(chain_in[ci+1]);
		}
		assert(discards.size() + 1 == chain.size());
	}

	//PARANOIA: subsequent chain verts need common simplex.
	for (auto const &chain : next_chains) {
		for (uint32_t ci = 0; ci + 1 < chain.size(); ++ci) {
			EmbeddedVertex::common_simplex(chain[ci].simplex, chain[ci+1].simplex);
		}
	}
	//end PARANOIA

	//compute lengths for each segment:
	std::vector< std::vector< float > > next_segment_lengths;
	for (auto const &chain : next_chains) {
		next_segment_lengths.emplace_back();
		std::vector< float > &lengths = next_segment_lengths.back();
		lengths.reserve(chain.size()-1);
		for (uint32_t ci = 0; ci + 1 < chain.size(); ++ci) {
			glm::vec3 a = chain[ci].interpolate(model.vertices);
			glm::vec3 b = chain[ci+1].interpolate(model.vertices);
			lengths.emplace_back(glm::length(b-a));
		}
	}

	//now flip discards flags based on lengths (make sure segments aren't too short):
	auto if_mixed_then_flatten_and_call = [](bool is_loop,
		std::vector< float > &lengths, std::vector< bool > &discards,
		std::function< void(std::vector< float >&, std::vector< bool > &)> const &call){
		//check that discards is mixed:
		bool have_false = false;
		bool have_true = false;
		for (auto d : discards) {
			if (d) have_true = true;
			else have_false = true;
		}
		if (!(have_false && have_true)) return;
		
		//flatten: if it's a loop, roll so that first element and last element have different discard flags
		uint32_t old_first = 0;
		auto DEBUG_old_lengths = lengths;
		if (is_loop) {
			uint32_t new_first = 0;
			if (discards.back() == discards[0]) {
				new_first = 1;
				while (new_first < discards.size() && discards[new_first-1] == discards[new_first]) ++new_first;
			}
			assert(new_first < discards.size());

			if (new_first == 0) {
				old_first = 0;
			} else {
				old_first = discards.size() - new_first;
			}

			std::rotate(discards.begin(), discards.begin() + new_first, discards.end());
			std::rotate(lengths.begin(), lengths.begin() + new_first, lengths.end());
			assert(discards[0] != discards.back());
		}

		call(lengths, discards);

		if (is_loop) {
			std::rotate(discards.begin(), discards.begin() + old_first, discards.end());
			std::rotate(lengths.begin(), lengths.begin() + old_first, lengths.end());
			assert(lengths == DEBUG_old_lengths);
		}
	};

	for (uint32_t idx = 0; idx < next_chains.size(); ++idx) {
		assert(next_chains[idx].size() == next_segment_lengths[idx].size() + 1);
		assert(next_chains[idx].size() == next_segment_discards[idx].size() + 1);
		bool is_loop = next_chains[idx][0] == next_chains[idx].back();

		//first, remove any non-discard segment shorter than 1.5 stitches:
		if_mixed_then_flatten_and_call(is_loop, next_segment_lengths[idx], next_segment_discards[idx], [&parameters](std::vector< float > &lengths, std::vector< bool > &discards){

			float MinSegmentLength = 1.5f * parameters.stitch_width_mm / parameters.model_units_mm;

			for (uint32_t begin = 0; begin < discards.size(); /* later */) {
				if (discards[begin]) {
					++begin;
					continue;
				}

				uint32_t end = begin + 1;
				float length = lengths[begin];
				while (end < discards.size() && !discards[end]) {
					length += lengths[end];
					++end;
				}

				if (length < MinSegmentLength) {
					std::cout << "setting discard on [" << begin << ", " << end << ") of [0, " << discards.size() << ")" << std::endl; //DEBUG
					//discard too-short segment:
					for (uint32_t i = begin; i < end; ++i) {
						discards[i] = true;
					}
				}

				begin = end;
			}
		});

		//then, remove any discard segment shorter than 0.5 stitches:
		if_mixed_then_flatten_and_call(is_loop, next_segment_lengths[idx], next_segment_discards[idx], [&parameters](std::vector< float > &lengths, std::vector< bool > &discards){

			float MinSegmentLength = 0.5f * parameters.stitch_width_mm / parameters.model_units_mm;

			for (uint32_t begin = 0; begin < discards.size(); /* later */) {
				if (!discards[begin]) {
					++begin;
					continue;
				}

				uint32_t end = begin + 1;
				float length = lengths[begin];
				while (end < discards.size() && discards[end]) {
					length += lengths[end];
					++end;
				}

				if (length < MinSegmentLength) {
					//mark too-short segment:
					std::cout << "setting keep on [" << begin << ", " << end << ") of [0, " << discards.size() << ")" << std::endl; //DEBUG
					for (uint32_t i = begin; i < end; ++i) {
						discards[i] = false;
					}
				}

				begin = end;
			}
		});

	}

	{ //if all segments are marked 'discard', then mark everything 'accept':
		bool only_discard = true;
		for (auto const &discards : next_segment_discards) {
			for (auto d : discards) {
				if (!d) {
					only_discard = false;
					break;
				}
			}
			if (only_discard == false) break;
		}

		if (only_discard) {
			std::cout << "Marking everything accept because it was all marked discard." << std::endl;
			for (auto &discards : next_segment_discards) {
				discards.assign(discards.size(), false);
			}
		} else {
			std::cout << "Have a mix of discard and accept." << std::endl;
		}
	}


	//find segments of active and next chains that are mutual nearest neighbors:

	auto make_locations = [&model](std::vector< std::vector< ak::EmbeddedVertex > > const &chains) {
		std::vector< std::vector< glm::vec3 > > locations;
		locations.reserve(chains.size());
		for (auto const &chain : chains) {
			locations.emplace_back();
			locations.back().reserve(chain.size());
			for (auto const &ev : chain) {
				locations.back().emplace_back(ev.interpolate(model.vertices));
			}
		}
		return locations;
	};

	std::vector< std::vector< glm::vec3 > > active_locations = make_locations(active_chains);
	std::vector< std::vector< glm::vec3 > > next_locations = make_locations(next_chains);

	auto make_vertex_closest = [](std::vector< std::vector< glm::vec3 > > const &locations, std::vector< std::vector< glm::vec3 > > const &targets) {
		std::vector< std::vector< uint32_t > > closest;
		closest.reserve(locations.size());
		for (auto const &location : locations) {
			closest.emplace_back();
			closest.back().reserve(location.size());
			for (auto const &l : location) {
				uint32_t close = -1U;
				float best2 = std::numeric_limits< float >::infinity();
				for (auto const &target : targets) {
					float dis2 = std::numeric_limits< float >::infinity();
					for (auto const &t : target) {
						dis2 = std::min(dis2, glm::length2(t - l));
					}
					if (dis2 < best2) {
						best2 = dis2;
						close = &target - &targets[0];
					}
				}
				closest.back().emplace_back(close);
			}
		}
		return closest;
	};

	auto make_segment_closest = [](std::vector< std::vector< glm::vec3 > > const &locations, std::vector< std::vector< glm::vec3 > > const &targets) {
		std::vector< std::vector< uint32_t > > closest;
		closest.reserve(locations.size()-1);
		for (auto const &location : locations) {
			closest.emplace_back();
			closest.back().reserve(location.size());
			for (uint32_t li = 0; li + 1 < location.size(); ++li) {
				glm::vec3 l = 0.5f * (location[li] + location[li+1]);
				uint32_t close = -1U;
				float best2 = std::numeric_limits< float >::infinity();
				for (auto const &target : targets) {
					float dis2 = std::numeric_limits< float >::infinity();
					for (auto const &t : target) {
						dis2 = std::min(dis2, glm::length2(t - l));
					}
					if (dis2 < best2) {
						best2 = dis2;
						close = &target - &targets[0];
					}
				}
				closest.back().emplace_back(close);
			}
			assert(closest.back().size() + 1 == location.size());
		}
		return closest;
	};



	std::vector< std::vector< uint32_t > > active_closest = make_vertex_closest(active_locations, next_locations);
	std::vector< std::vector< uint32_t > > next_segment_closest = make_segment_closest(next_locations, active_locations);

	//sort active and next into matching sub-chains:
	struct BeginEnd {
		BeginEnd(uint32_t begin_, uint32_t end_) : begin(begin_), end(end_) { }
		uint32_t begin;
		uint32_t end;
	};

	struct Match {
		std::vector< BeginEnd > active; //vertex ranges
		std::vector< BeginEnd > next; //segment ranges
	};

	std::map< std::pair< uint32_t, uint32_t >, Match > matches;

	for (auto const &closest : active_closest) {
		//TODO: make sure nothing appears more than twice in closest
		uint32_t ai = &closest - &active_closest[0];

		for (uint32_t begin = 0; begin < closest.size(); /* later */) {
			uint32_t end = begin + 1;
			while (end < closest.size() && closest[end] == closest[begin]) ++end;
			matches[std::make_pair(ai, closest[begin])].active.emplace_back(begin, end);
			begin = end;
		}
	}

	for (auto const &closest : next_segment_closest) {
		//TODO: make sure nothing appears more than twice in closest
		uint32_t ni = &closest - &next_segment_closest[0];

		for (uint32_t begin = 0; begin < closest.size(); /* later */) {
			uint32_t end = begin + 1;
			while (end < closest.size() && closest[end] == closest[begin]) ++end;
			matches[std::make_pair(closest[begin], ni)].next.emplace_back(begin, end);
			begin = end;
		}
	}

	std::vector< uint32_t > active_matches(active_chains.size(), 0);
	std::vector< uint32_t > next_matches(next_chains.size(), 0);
	for (auto const &anm : matches) {
		active_matches[anm.first.first] += 1;
		next_matches[anm.first.second] += 1;
	}

	{ //If there are any merges or splits, all participating next cycles are marked 'accept':
		std::set< uint32_t > to_mark;
		for (auto const &anm : matches) {
			bool is_split_or_merge = (active_matches[anm.first.first] > 1 || active_matches[anm.first.second] > 1);
			if (is_split_or_merge) {
				to_mark.insert(anm.first.second);
			}
		}
		uint32_t were_marked = 0;
		for (auto ni : to_mark) {
			for (auto d : next_segment_discards[ni]) {
				if (d) ++were_marked;
			}
			next_segment_discards[ni].assign(next_segment_discards[ni].size(), false);
		}
		if (!to_mark.empty() && were_marked != 0) {
			std::cerr << "Marked " << were_marked << " segments on " << to_mark.size() << " next cycles as 'accept' based on participating in a merge/split." << std::endl;
		}
	}

	struct NewStitch {
		NewStitch(uint32_t chain_, uint32_t segment_, float along_, ak::Flag flag_) : chain(chain_), segment(segment_), along(along_), flag(flag_) { }
		uint32_t chain;
		uint32_t segment;
		float along;
		ak::Flag flag;
	};

	std::vector< NewStitch > all_new_stitches;

	//allocate stitch counts based on segment lengths + source stitches:
	//(and limit based on active stitch counts)
	//TODO: account for balance also
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

		//compute min/max totals from active chain flags:
		uint32_t active_ones = 0;
		uint32_t active_anys = 0;
		uint32_t DEBUG_stitch_count = 0;
		{
			assert(anm.first.first < active_flags.size());
			bool is_loop = (active_chains[anm.first.first][0] == active_chains[anm.first.first].back());
			for (auto const &be : match.active) {
				for (uint32_t i = be.begin; i < be.end; ++i) {
					if (is_loop && i + 1 == active_flags[anm.first.first].size()) continue;
					auto f = active_flags[anm.first.first][i];
					if (f == ak::FlagLinkOne) ++active_ones;
					else if (f == ak::FlagLinkAny) ++active_anys;
				}
			}
			DEBUG_stitch_count = active_ones + active_anys;
		}

		//compute min for next segments based on short row ends constraints:

		std::vector< bool > const &discards = next_segment_discards[anm.first.second];
		bool next_is_cycle = (next_chains[anm.first.second][0] == next_chains[anm.first.second].back());
		bool discards_before_front = (next_is_cycle ? discards.back() : discards[0]);
		bool discards_after_back = (next_is_cycle ? discards[0] : discards.back());

		uint32_t next_ones = 0;

		{ //compute next_ones:
			//every discard/non-discard edge needs at least one stitch next to it on each side:
			for (auto const &be : match.next) {
				for (uint32_t i = be.begin; i < be.end; ++i) {
					if (discards[i] != (i > 0 ? discards[i-1] : discards_before_front)) next_ones += 1;
					if (discards[i] != (i + 1 < discards.size() ? discards[i+1] : discards_after_back)) next_ones += 1;
				}
			}

			if (next_ones > 2 * active_anys + active_ones) {
				std::cerr << "ERROR: more discard/non-discard ends are required (" << next_ones << ") than are permitted by the current active flags (" << active_anys << "*2 + " << active_ones << "); code to fix this (by removing shortest same-discard segment) not yet implemented." << std::endl;
				assert(next_ones <= 2 * active_anys + active_ones);
			}
		}

		//compute desired stitch count based on segment lengths:
		assert(anm.first.second < next_segment_lengths.size());
		std::vector< float > const &lengths = next_segment_lengths[anm.first.second];
		float total_length = 0.0f;
		for (auto const &be : match.next) {
			assert(be.end <= lengths.size());
			for (uint32_t i = be.begin; i < be.end; ++i) {
				total_length += lengths[i];
			}
		}
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

		std::vector< NewStitch > new_stitches;
		if (stitches > 0) {
			//spread stitches among "allocation ranges" (same discard status)
			struct Alloc {
				Alloc(uint32_t begin_, uint32_t end_, bool first_one_, bool last_one_, float length_) : begin(begin_), end(end_), first_one(first_one_), last_one(last_one_), stitches((first_one_ ? 1 : 0) + (last_one_ ? 1 : 0)), length(length_) { }
				uint32_t begin;
				uint32_t end;
				bool first_one;
				bool last_one;
				uint32_t stitches;
				float length;
			};
			std::vector< Alloc > alloc;
			uint32_t total_ones = 0;
			for (auto const &be : match.next) {
				for (uint32_t begin = be.begin; begin < be.end; /* later */) {
					uint32_t end = begin;
					while (end < be.end && discards[begin] == discards[end]) ++end;
					bool first_one = (discards[begin] != (begin > 0 ? discards[begin-1] : discards_before_front));
					bool last_one = (discards[end-1] != (end < discards.size() ? discards[end] : discards_after_back));

					float length = 0.0f;
					for (uint32_t i = begin; i < end; ++i) {
						length += lengths[i];
					}

					alloc.emplace_back(begin, end, first_one, last_one, length);
					total_ones += alloc.back().stitches;

					begin = end;
				}
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
				float stitch_step = a.length / float(a.stitches);
				float stitch_acc = 0.5f * stitch_step; //TODO: try different offsets?
				uint32_t s = 0;
				for (uint32_t i = a.begin; i < a.end; ++i) {
					float length = lengths[i];
					float remain = length;
					while (stitch_acc < remain) {
						remain -= stitch_acc;
						ak::Flag flag = ak::FlagLinkAny;
						if (s == 0 && a.first_one) flag = ak::FlagLinkOne;
						if (s + 1 == a.stitches && a.last_one) flag = ak::FlagLinkOne;
						new_stitches.emplace_back( anm.first.second, i, (length - remain) / length, flag);
						++s;
						stitch_acc = stitch_step;
					}
					stitch_acc -= remain;
				}
				assert(s == a.stitches);
			}
		}
		assert(new_stitches.size() == stitches);


		std::vector< bool > new_stitch_linkone;
		std::vector< glm::vec3 > new_stitch_locations;
		new_stitch_locations.reserve(new_stitches.size());
		for (auto const &n : new_stitches) {
			assert(n.chain < next_chains.size());
			assert(n.segment < next_chains[n.chain].size());
			assert(n.segment + 1 < next_chains[n.chain].size());
			new_stitch_locations.emplace_back(
				glm::mix(
					next_chains[n.chain][n.segment].interpolate(model.vertices),
					next_chains[n.chain][n.segment+1].interpolate(model.vertices),
					n.along
				)
			);
			new_stitch_linkone.emplace_back(n.flag == ak::FlagLinkOne);
		}

		//figure out how to link to next stitches:
		std::vector< uint32_t > active_stitch_indices;
		std::vector< bool > active_stitch_linkone;
		std::vector< glm::vec3 > active_stitch_locations;
		{
			assert(anm.first.first < active_chains.size());
			auto const &active_chain = active_chains[anm.first.first];
			auto const &active_flag = active_flags[anm.first.first];
			assert(active_chain.size() == active_flag.size());
			bool is_loop = (active_chain[0] == active_chain.back());
			for (auto const &be : match.active) {
				for (uint32_t i = be.begin; i < be.end; ++i) {
					if (is_loop && (i + 1 == active_flag.size())) continue; //skip last index in loops
					if (active_flag[i] == ak::FlagLinkAny || active_flag[i] == ak::FlagLinkOne ) {
						active_stitch_indices.emplace_back( i );
						active_stitch_locations.emplace_back( active_chain[i].interpolate(model.vertices) );
						active_stitch_linkone.emplace_back( active_flag[i] == ak::FlagLinkOne );
					}
				}
			}
		}
		assert(active_stitch_indices.size() == DEBUG_stitch_count);

		{ //DEBUG:
			uint32_t new_ones = 0;
			for (auto o : new_stitch_linkone) {
				if (o) ++new_ones;
			}
			uint32_t active_ones = 0;
			for (auto o : active_stitch_linkone) {
				if (o) ++active_ones;
			}
			std::cout << " About to connect " << active_stitch_locations.size() << " active stitches (" << active_ones << " linkones) to " << new_stitch_locations.size() << " new stitches (" << new_ones << " linkones)." << std::endl;
		}

		//actually build links:
		{ //least-clever linking solution: FlagLinkOne's link 1-1, others link to keep arrays mostly in sync

			std::vector< std::pair< uint32_t, uint32_t > > best_links;
			float best_cost = std::numeric_limits< float >::infinity();

			auto try_links_even = [&](uint32_t roll_active, uint32_t roll_new) {
				std::vector< std::pair< uint32_t, uint32_t > > possible_links;

				if (active_stitch_locations.size() <= new_stitch_locations.size()) {
					//evenly distribute increases among the non-linkone stitches:
					uint32_t total = 0;
					for (auto l : active_stitch_linkone) {
						if (!l) ++total;
					}
					uint32_t increases = new_stitch_locations.size() - active_stitch_locations.size();
					std::vector< bool > inc(total, false);
					for (uint32_t i = 0; i < increases; ++i) {
						assert(inc[i * total / increases] == false);
						inc[i * total / increases] = true;
					}
					uint32_t n = 0;
					uint32_t i = 0;
					for (uint32_t a = 0; a < active_stitch_locations.size(); ++a) {
						uint32_t ra = (a + roll_active) % active_stitch_locations.size();
						uint32_t rn = (n + roll_new) % new_stitch_locations.size();
						possible_links.emplace_back(ra, rn);
						++n;
						if (!active_stitch_linkone[ra]) {
							assert(i < inc.size());
							if (inc[i]) {
								rn = (n + roll_new) % new_stitch_locations.size();
								possible_links.emplace_back(ra, rn);
								++n;
							}
							++i;
						}
					}
					assert(i == total);
					assert(n == new_stitch_locations.size());
				} else if (active_stitch_locations.size() > new_stitch_locations.size()) {
					//evenly distribute decreases among the non-linkone stitches:
					uint32_t total = 0;
					for (auto l : new_stitch_linkone) {
						if (!l) ++total;
					}
					uint32_t decreases = active_stitch_locations.size() - new_stitch_locations.size();
					std::vector< bool > dec(total, false);
					for (uint32_t i = 0; i < decreases; ++i) {
						assert(dec[i * total / decreases] == false);
						dec[i * total / decreases] = true;
					}
					uint32_t a = 0;
					uint32_t i = 0;
					for (uint32_t n = 0; n < new_stitch_locations.size(); ++n) {
						uint32_t ra = (a + roll_active) % active_stitch_locations.size();
						uint32_t rn = (n + roll_new) % new_stitch_locations.size();
						possible_links.emplace_back(ra, rn);
						++a;
						if (!new_stitch_linkone[rn]) {
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
					float len = glm::length(new_stitch_locations[p.second] - active_stitch_locations[p.first]);
					cost += (len - row_height) * (len - row_height);
				}

				if (cost < best_cost) {
					best_cost = cost;
					best_links = possible_links;
				}
			};


			if (active_stitch_locations.size() >= new_stitch_locations.size()) {
				for (uint32_t roll_active = 0; roll_active < active_stitch_locations.size(); ++roll_active) {
					try_links_even(roll_active,0);
				}
			} else {
				for (uint32_t roll_new = 0; roll_new < new_stitch_locations.size(); ++roll_new) {
					try_links_even(0,roll_new);
				}
			}


			for (auto const &p : best_links) {
				Link link;
				link.from_chain = anm.first.first;
				assert(p.first < active_stitch_indices.size());
				link.from_vertex = active_stitch_indices[p.first];

				//will get translated to real index in divide-at-stitches code:
				link.to_chain = -1U;
				assert(p.second < new_stitches.size());
				link.to_vertex = all_new_stitches.size() + p.second;
				links.emplace_back(link);
			}
		}


		all_new_stitches.insert(all_new_stitches.end(), new_stitches.begin(), new_stitches.end());
	}

	//insert new stitches into next chains to make linked_next_chains.
	std::map< uint32_t, std::map< std::pair< uint32_t, float >, NewStitch * > > splits;
	for (auto &ns : all_new_stitches) {
		auto &m = splits[ns.chain];
		auto res = m.insert(std::make_pair( std::make_pair(ns.segment, ns.along), &ns ));
		assert(res.second); //shouldn't have two overlapping stitches!
	}

	std::vector< std::pair< uint32_t, uint32_t > > ns_index(all_new_stitches.size(), std::make_pair(-1U, -1U));

	for (auto const &chain : next_chains) {
		auto const &discards = next_segment_discards[&chain - &next_chains[0]];

		bool next_is_cycle = (chain[0] == chain.back());
		bool discards_before_front = (next_is_cycle ? discards.back() : discards[0]);
		bool discards_after_back = (next_is_cycle ? discards[0] : discards.back());

		auto const &m = splits[&chain - &next_chains[0]];

		linked_next_chains.emplace_back();
		auto &linked_chain = linked_next_chains.back();
		linked_next_flags.emplace_back();
		auto &flags = linked_next_flags.back();

		uint32_t s = 0;
		auto add_until = [&](uint32_t until) {
			while (s <= until) {
				assert(s < chain.size());
				linked_chain.emplace_back(chain[s]);
				bool discard = 
					(s > 0 ? discards[s-1] : discards_before_front)
					&& (s < discards.size() ? discards[s] : discards_after_back);
				if (discard) flags.emplace_back(ak::FlagDiscard);
				else flags.emplace_back(ak::FlagLinkNone);
				++s;
			}
		};
		for (auto const &san : m) {
			add_until(san.first.first);
			assert(s == san.first.first+1);
			assert(s < chain.size());
			assert(ns_index[san.second - &all_new_stitches[0]] == std::make_pair(-1U, -1U));
			ns_index[san.second - &all_new_stitches[0]] = std::make_pair(
				&linked_chain - &linked_next_chains[0],
				linked_chain.size()
			);
			linked_chain.emplace_back( EmbeddedVertex::mix( chain[s-1], chain[s], san.first.second ) );
			if (discards[s-1]) flags.emplace_back(ak::FlagDiscard);
			else flags.emplace_back(san.second->flag);
		}
		add_until(chain.size()-1);
	}

	//PARANOIA: subsequent chain verts need common simplex.
	for (auto const &chain : linked_next_chains) {
		for (uint32_t ci = 0; ci + 1 < chain.size(); ++ci) {
			EmbeddedVertex::common_simplex(chain[ci].simplex, chain[ci+1].simplex);
		}
	}
	//end PARANOIA

	//translate links:
	for (auto &link : links) {
		assert(link.to_chain == -1U);
		assert(link.to_vertex < ns_index.size());
		link.to_chain = ns_index[link.to_vertex].first;
		link.to_vertex = ns_index[link.to_vertex].second;
		assert(link.to_chain < linked_next_chains.size());
		assert(link.to_vertex < linked_next_chains[link.to_chain].size());
	}

}

