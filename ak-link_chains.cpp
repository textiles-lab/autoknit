#include "pipeline.hpp"

void ak::link_chains(
	ak::Parameters const &parameters,
	ak::Model const &model, //in: model
	std::vector< float > const &times,          //in: time field (times @ vertices)
	std::vector< std::vector< ak::EmbeddedVertex > > const &active_chains, //in: current active chains
	std::vector< std::vector< ak::Flag > > const &active_flags, //in: stitches
	std::vector< std::vector< ak::EmbeddedVertex > > const &next_chains, //in: next chains
	std::vector< std::vector< ak::EmbeddedVertex > > *linked_next_chains_, //out: next chains
	std::vector< std::vector< ak::Flag > > *linked_next_flags_, //out: flags indicating status of vertices on next chains
	std::vector< ak::Link > *links_ //out: active_chains[from_chain][from_vertex] -> linked_next_chains[to_chain][to_vertex] links
) {
	assert(times.size() == model.vertices.size());

	assert(linked_next_chains_);
	auto &linked_next_chains = *linked_next_chains_;
	linked_next_chains.clear();

	assert(linked_next_flags_);
	auto &linked_next_flags = *linked_next_flags_;
	linked_next_flags.clear();


	//mark next chains for discard based on time values from active chains:
	std::vector< std::vector< ak::Flag > > next_flags;
	{
		float active_max_time = -std::numeric_limits< float >::infinity();
		for (auto const &chain : active_chains) {
			for (auto const &ev : chain) {
				active_max_time = std::max(active_max_time, ev.interpolate(times));
			}
		}

		next_flags.reserve(next_chains.size());
		for (auto const &chain : next_chains) {
			next_flags.emplace_back();
			next_flags.back().reserve(chain.size());
			auto &flags = next_flags.back();

			uint32_t discard_count = 0;
			uint32_t any_count = 0;
			//mark any location on next with time > active_max_time as 'discard':
			for (auto const &ev : chain) {
				float time = ev.interpolate(times);
				if (time > active_max_time) {
					flags.emplace_back(ak::FlagDiscard);
					++discard_count;
				} else {
					next_flags.back().emplace_back(ak::FlagLinkAny);
					++any_count;
				}
			}

			//associate a bit of length with each sample on the chain:
			std::vector< float > length(chain.size(), 0.0f);
			for (uint32_t ci = 0; ci + 1 < chain.size(); ++ci) {
				glm::vec3 a = chain[ci].interpolate(model.vertices);
				glm::vec3 b = chain[ci+1].interpolate(model.vertices);
				float l = glm::length(a,b);
				length[ci] += 0.5f * l;
				length[ci+1] += 0.5f * l;
			}
			if (chain[0] == chain.back()) {
				//it's a loop, so both front and back get sample (summed) length:
				length[0] += length.back();
				length.back() = length[0];
			}

			//Now there are two transforms that are easier to do when things are flattened, so here's a helper for that:

			auto if_mixed_when_flat = [&flags,&lengths,&chain](std::function< void(std::vector< ak::Flag > &, std;:vector< float > &) > const &do_this){
				assert(chain.size() == flags.size());
				assert(chain.size() == lengths.size());
				//if chain is a loop, then so must be flags and lengths:
				bool is_loop = (chain[0] == chain.back());
				assert(!is_loop || flags[0] == flags.back());
				assert(!is_loop || lengths[0] == lengths.back());

				//only proceed if there is a mix of stitch types in 'flags':
				bool have_discard = false;
				bool have_any = false;
				for (auto f : flags) {
					if (f == ak::FlagDiscard) have_discard = true;
					else have_any = true;
				}
				if (!(have_discard && have_any)) return;

				//reorganize flags and lengths so that they start on a flag-segment boundary:
				uint32_t old_first = 0;
				auto DEBUG_old_lengths = lengths;
				if (is_loop) {
					uint32_t new_first = 0;
					while (new_first + 1 < flags.size() && flags[new_first] == flags[new_first+1]) ++new_first;
					assert(new_first + 1 < flags.size());

					assert(flags.back() == flags[0]);
					flags.pop_back();
					assert(lengths.back() == lengths[0]);
					lengths.pop_back();

					old_first = (new_first == 0 ? 0 : flags.size() - new_first);

					std::rotate(flags.begin(), flags.begin() + new_first, flags.end());
					std::rotate(lengths.begin(), lengths.begin() + new_first, lengths.end());

					assert(flags[0] != flags.back());
				}

				do_this(flags, lengths);

				//now un-reorganize:
				if (is_loop) {
					std::rotate(flags.begin(), flags.begin() + old_first, flags.end());
					std::rotate(lengths.begin(), lengths.begin() + old_first, lengths.end());

					flags.emplace_back(flags[0]);
					lengths.emplace_back(lengths[0]);

					assert(lengths == DEUBG_old_lengths);
				}
			};

			//if any non-discard segment is shorter than ~two stitches, mark it discard:
			if_mixed_when_flat([](std::vector< ak::Flag > &flags, std::vector< float > &lengths){
				float MinSegmentLength = 1.5f * parameters.stitch_width_mm / parameters.model_units_mm;

				for (uint32_t begin = 0; begin < flags.size(); /* later */) {
					if (flags[begin] == ak::FlagDiscard) {
						++begin;
						continue;
					}
					assert(flags[begin] == ak::FlagLinkAny);

					uint32_t end = begin + 1;
					float length = lengths[begin];
					while (end < flags.size() && flags[end] == flags[begin]) {
						length += lengths[end];
						++end;
					}

					if (length < MinSegmentLength) {
						//discard too-short segment:
						for (uint32_t i = begin; i < end; ++i) {
							flags[i] = ak::FlagDiscard;
						}
					}

					begin = end;
				}
			});

			//if any discard segment is shorter than ~one stitch, mark it non-discard:
			if_mixed_when_flat([](std::vector< ak::Flag > &flags, std::vector< float > &lengths){
				float MinSegmentLength = 0.5f * parameters.stitch_width_mm / parameters.model_units_mm;

				for (uint32_t begin = 0; begin < flags.size(); /* later */) {
					if (flags[begin] != ak::FlagDiscard) {
						++begin;
						continue;
					}
					assert(flags[begin] == ak::FlagDiscard);

					uint32_t end = begin + 1;
					float length = lengths[begin];
					while (end < flags.size() && flags[end] == flags[begin]) {
						length += lengths[end];
						++end;
					}

					if (length < MinSegmentLength) {
						//discard too-short segment:
						for (uint32_t i = begin; i < end; ++i) {
							flags[i] = ak::FlagLinkAny;
						}
					}

					begin = end;
				}
			});

		}

		//if all flags are marked 'discard', then mark everything 'accept':
		bool only_discard = true;
		for (auto const &flags : next_flags) {
			for (auto f : flags) {
				if (f != ak::FlagDiscard) {
					only_discard = false;
					break;
				}
			}
			if (only_discard == false) break;
		}

		if (only_discard) {
			std::cout << "All flags were discard, so marking all to any." << std::endl;

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
				locations.back().emplace_back(ev.interpolate(model));
			}
		}
		return locations;
	};

	std::vector< std::vector< glm::vec3 > > active_locations = make_locations(active_chains);
	std::vector< std::vector< glm::vec3 > > next_locations = make_locations(next_chains);

	auto make_closest = [](std::vector< std::vector< glm::vec3 > > const &locations, std::vector< std::vector< glm::vec3 > > const &targets) {
		std::vector< std::vector< uint32_t > > closest;
		closest.reserve(locations.size());
		for (auto const &location : locations) {
			closest.emplace_back();
			closest.back().reserve(location.size());
			for (auto const &l : location) {
				uint32_t close = -1U;
				float best2 = std::numeric_limits< float >::infinity();
				for (auto const &target : targets) {
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


	std::vector< std::vector< uint32_t > > active_closest = make_closest(active_locations, next_locations);
	std::vector< std::vector< uint32_t > > next_closest = make_closest(next_locations, active_locations);

	//sort active and next into matching segments:
	struct Segment {
		Segment(uint32_t begin_, uint32_t end_) : begin(begin_), end(end_) { }
		uint32_t begin;
		uint32_t end;
	};

	struct Match {
		std::vector< Segment > active;
		std::vector< Segment > next;
		std::vector< uint32_t > allocation; //stitch counts for next segments
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

	for (auto const &closest : next_closest) {
		//TODO: make sure nothing appears more than twice in closest
		uint32_t ni = &closest - &next_closest[0];

		for (uint32_t begin = 0; begin < closest.size(); /* later */) {
			uint32_t end = begin + 1;
			while (end < closest.size() && closest[end] == closest[begin]) ++end;
			matches[std::make_pair(closest[begin], ni)].next.emplace_back(begin, end);
			begin = end;
		}
	}

	//allocate stitch counts based on segment lengths:
	//(and limit based on active stitch counts)
	//TODO: account for balance also
	for (auto const &anm : matches) {
		assert(anm.first.second < next_chains.size());
		std::vector< ak::EmbeddedVertex > const &next = next_chains[anm.first.second];
		Match const &match = anm.second;
		//compute min/max totals from active chain flags:
		uint32_t min = 0;
		uint32_t max = 0;
		{
			uint32_t link_one = 0;
			uint32_t link_any = 0;
			assert(anm.first.first < active_flags.size());
			for (auto f : active_flags[anm.first.first]) {
				if (f == ak::FlagLinkOne) ++link_one;
				else if (f == ak::FlagLinkAny) ++link_any;
			}
			//link_any stitches can be decreases:
			min = link_one + (link_any + 1) / 2;
			//link_any stitches can be increases:
			max = link_one + link_any * 2;
		}

	}


	//TODO: allocate stitch counts to matches in a way that respects balance.

	for (auto const &anm : matches) {
		uint32_t active = anm.first.first;
		uint32_t next = anm.first.first;
		Match const &match = anm.second;

		if (match.active.empty()) {
			std::cerr << "WARNING: active stitches with no next chain." << std::endl;
			continue;
		}
		if (match.next.empty()) {
			std::cerr << "WARNING: next stitches with no active chain. (Unlike the active-with-no-next case, this should *NEVER* happen.)" << std::endl;
			continue;
		}
		//Steps
	}

}

