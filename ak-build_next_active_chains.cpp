#include "pipeline.hpp"

#include <iostream>
#include <map>
#include <unordered_set>

void ak::build_next_active_chains(
	ak::Parameters const &parameters,
	ak::Model const &model,
	std::vector< std::vector< ak::EmbeddedVertex > > const &active_chains, //in: current active chains
	std::vector< std::vector< ak::Flag > > const &active_flags_in, //in: flags for current active
	std::vector< std::vector< ak::EmbeddedVertex > > const &next_chains, //in: next chains
	std::vector< std::vector< ak::Flag > > const &next_flags_in, //in: flags for next active
	std::vector< ak::Link > const &links, //in: links between active and next
	std::vector< std::vector< ak::EmbeddedVertex > > *next_active_chains_, //out: next active chains
	std::vector< std::vector< ak::Flag > > *next_active_flags_ //out: next stitch flags
) {
	assert(active_chains.size() == active_flags_in.size());
	for (uint32_t a = 0; a < active_chains.size(); ++a) {
		assert(active_chains[a].size() == active_flags_in[a].size());
	}

	assert(next_chains.size() == next_flags_in.size());
	for (uint32_t a = 0; a < next_chains.size(); ++a) {
		assert(next_chains[a].size() == next_flags_in[a].size());
	}

	for (auto const &l : links) {
		assert(l.from_chain < active_chains.size());
		assert(l.from_vertex < active_chains[l.from_chain].size());
		assert(l.to_chain < next_chains.size());
		assert(l.to_vertex < next_chains[l.to_chain].size());
	}

	assert(next_active_chains_);
	auto &next_active_chains = *next_active_chains_;
	next_active_chains.clear();

	assert(next_active_flags_);
	auto &next_active_flags = *next_active_flags_;
	next_active_flags.clear();


	//build a lookup structure for links:
	struct ChainVertex {
		ChainVertex(uint32_t chain_, uint32_t vertex_) : chain(chain_), vertex(vertex_) { }
		uint32_t chain;
		uint32_t vertex;
		bool operator<(ChainVertex const &o) const {
			if (chain != o.chain) return chain < o.chain;
			else return vertex < o.vertex;
		}
	};

	std::multimap< ChainVertex, ChainVertex > active_next;
	std::multimap< ChainVertex, ChainVertex > next_active;

	for (auto const &l : links) {
		active_next.insert(std::make_pair(
			ChainVertex(l.from_chain, l.from_vertex),
			ChainVertex(l.to_chain, l.to_vertex)
		));
		next_active.insert(std::make_pair(
			ChainVertex(l.to_chain, l.to_vertex),
			ChainVertex(l.from_chain, l.from_vertex)
		));
	}


	//flood fill discards right up to stitches:
	auto fill_discards = [](std::vector< std::vector< ak::EmbeddedVertex > > const &in_chains,
		std::vector< std::vector< ak::Flag > > &in_flags) {
		for (auto &flags : in_flags) {
			auto const &chain = in_chains[&flags - &in_flags[0]];
			bool is_loop = (chain[0] == chain.back());
			for (uint32_t i = 0; i + 1 < flags.size(); ++i) {
				if (flags[i] == ak::FlagDiscard && flags[i+1] == ak::FlagLinkNone) {
					flags[i+1] = ak::FlagDiscard;
				}
			}
			if (is_loop && flags[0] != flags.back()) {
				flags[0] = flags.back();
				for (uint32_t i = 0; i + 1 < flags.size(); ++i) {
					if (flags[i] == ak::FlagDiscard && flags[i+1] == ak::FlagLinkNone) {
						flags[i+1] = ak::FlagDiscard;
					}
				}
			}
			for (uint32_t i = flags.size() - 1; i > 0; --i) {
				if (flags[i] == ak::FlagDiscard && flags[i-1] == ak::FlagLinkNone) {
					flags[i-1] = ak::FlagDiscard;
				}
			}
			if (is_loop && flags[0] != flags.back()) {
				flags.back() = flags[0];
				for (uint32_t i = flags.size() - 1; i > 0; --i) {
					if (flags[i] == ak::FlagDiscard && flags[i-1] == ak::FlagLinkNone) {
						flags[i-1] = ak::FlagDiscard;
					}
				}
			}
		}
		return in_flags;
	};

	//for next_flags, flood fill existing discards:
	std::vector< std::vector< ak::Flag > > next_flags = next_flags_in;
	fill_discards(next_chains, next_flags);

	//next_edges tells you if a given next index is at the edge of a discard range:
	std::vector< std::unordered_set< uint32_t > > next_edges(next_chains.size());
	for (uint32_t c = 0; c < next_chains.size(); ++c) {
		auto const &chain = next_chains[c];
		auto const &flags = next_flags[c];
		bool is_loop = (chain[0] == chain.back());
		ak::Flag before_front = (is_loop ? flags[flags.size()-2] : flags[0]);
		for (uint32_t i = 0; i < flags.size(); ++i) {
			//ignore last flag on a loop (handled as first):
			if (i + 1 == flags.size() && is_loop) continue;
			//ignore things that aren't stitches:
			if (!(flags[i] == ak::FlagLinkOne || flags[i] == ak::FlagLinkAny)) continue;
			if ( ((i > 0 ? flags[i-1] : before_front) == ak::FlagDiscard)
			  ||  (i + 1 < flags.size() && flags[i+1] == ak::FlagDiscard) ) {
				assert(flags[i] == ak::FlagLinkOne); //edges are marked LinkOne always
				next_edges[c].insert(i);
			}
		}
	}


	//for active_flags, mark any stitch with links to non-edge stitches as 'discard', then fill:
	std::vector< std::vector< ak::Flag > > active_flags = active_flags_in;
	std::vector< std::vector< uint32_t > > active_to_none(active_flags.size());
	for (uint32_t c = 0; c < active_chains.size(); ++c) {
		auto &flags = active_flags[c];
		bool is_loop = (active_chains[c][0] == active_chains[c].back());
		for (uint32_t i = 0; i < flags.size(); ++i) {
			//handle last flag same as first in loop:
			if (i + 1 == flags.size() && is_loop) {
				flags[i] = flags[0];
				continue;
			}
			//ignore non-stitches:
			if (!(flags[i] == ak::FlagLinkOne || flags[i] == ak::FlagLinkAny)) continue;
			//mark any stitches linked to non-edges for discard:
			auto r = active_next.equal_range(ChainVertex(c, i));
			bool have_edge = false;
			bool have_link = false;
			for (auto ri = r.first; ri != r.second; ++ri) {
				ak::Flag f = next_flags[ri->second.chain][ri->second.vertex];
				if (f == ak::FlagDiscard) {
					//ignore this link, target is discarded
				} else if (next_edges[ri->second.chain].count(ri->second.vertex)) {
					assert(f == ak::FlagLinkOne); //edges are always marked LinkOne
					have_edge = true;
				} else {
					assert(f == ak::FlagLinkOne || f == ak::FlagLinkAny); //this should be a link to a valid-flagged stitch
					have_link = true;
				}
			}
			if (have_edge) {
				active_to_none[c].emplace_back(i);
				if (i == 0 && is_loop) {
					active_to_none[c].emplace_back(flags.size()-1);
				}
			} else if (have_link) {
				flags[i] = ak::FlagDiscard;
			}
		}
	}
	fill_discards(active_chains, active_flags);
	//finally, mark the stitches that were kept because they were linked to an edge as 'LinkNone':
	for (uint32_t c = 0; c < active_to_none.size(); ++c) {
		for (auto i : active_to_none[c]) {
			assert( active_flags[c][i] == ak::FlagLinkOne || active_flags[c][i] == ak::FlagLinkAny );
			active_flags[c][i] = ak::FlagLinkNone;
		}
		//make sure loop-ness was preserved:
		assert(active_chains[c][0] != active_chains[c].back() || active_flags[c][0] == active_flags[c].back());
	}


	//divide chains into a series of 'fragments' between discarded segments:
	struct Fragment {
		enum {
			OnActive,
			OnNext
		} on;
		uint32_t chain;
		std::vector< uint32_t > inds;
	};

	std::vector< Fragment > fragments;

	//fragments are between discarded segments:
	auto add_fragments = [&fragments,&next_active_chains,&next_active_flags](
		std::vector< std::vector< ak::EmbeddedVertex > > const &src_chains,
		std::vector< std::vector< ak::Flag > > const &src_flags,
		decltype(Fragment().on) on) {

		for (auto const &chain : src_chains) {
			uint32_t chain_index = &chain - &src_chains[0];
			auto const &flags = src_flags[chain_index];

			bool has_discard = false;
			bool has_stitches = false;
			for (auto f : flags) {
				if (f == ak::FlagDiscard) has_discard = true;
				if (f == ak::FlagLinkOne || f == ak::FlagLinkAny) has_stitches = true;
			}

			//skip entirely if no stitches exist (e.g., it's all-discard):
			if (!has_stitches) continue;

			//just copy to output if no discards exist:
			if (!has_discard) {
				next_active_chains.emplace_back(chain);
				next_active_flags.emplace_back(flags);
				continue;
			}

			//have some stitches and some discards, so make fragments:
			std::vector< std::vector< uint32_t > > frags;
			for (uint32_t i = 0; i < flags.size(); ++i) {
				if (flags[i] == ak::FlagDiscard) continue;
				frags.emplace_back();
				std::vector< uint32_t > &inds = frags.back();
				inds.emplace_back(i);
				while (i+1 < flags.size() && flags[i+1] != ak::FlagDiscard) {
					inds.emplace_back(i+1);
					++i;
				}
			}
			assert(!frags.empty());
			if (chain[0] == chain.back()) {
				//for loops, might have a fragment that bridges first/last:
				if (frags[0][0] == 0) {
					//... if so, merge first/last fragments:
					assert(frags.size() >= 2); //must have more than one fragment because there was some discard in there
					assert(frags.back().back() + 1 == flags.size()); //last frag should have grabbed last ind (same as first)
					std::vector< uint32_t > inds = std::move(frags.back());
					inds.insert(inds.end(), frags[0].begin(), frags[0].end());
					frags[0] = std::move(inds);
					frags.pop_back();
				}
			}

			//transform recorded index chains into actual fragments:
			for (auto &inds : frags) {
				fragments.emplace_back();
				fragments.back().on = on;
				fragments.back().chain = chain_index;
				fragments.back().inds = std::move(inds);
			}
		}
	};

	add_fragments(next_chains, next_flags, Fragment::OnNext);
	add_fragments(active_chains, active_flags, Fragment::OnActive);


	//DEBUG: dump fragments as individual output chains
	for (auto const &frag : fragments) {
		next_active_chains.emplace_back();
		next_active_flags.emplace_back();
		std::vector< ak::EmbeddedVertex > const &src_chain = (frag.on == Fragment::OnNext ? next_chains : active_chains)[frag.chain];
		std::vector< ak::Flag > const &src_flags = (frag.on == Fragment::OnNext ? next_flags : active_flags)[frag.chain];
		for (auto i : frag.inds) {
			next_active_chains.back().emplace_back(src_chain[i]);
			next_active_flags.back().emplace_back(src_flags[i]);
		}
	}


#if 0 //I'm thinking this walking strategy gets messy with, say, non-linked chains and figuring out good seeding

	{ //now walk around active and next chains, reading out new chains.
		//basic step:
		// if on next, is next stitch a discard?
		//   yes -> follow link to active, mark active as LinkNone
		//   no  -> continue to next, preserve flag
		// if on active, is there a link to non-discard?
		//   yes -> follow link to non-discard, mark active as LinkNone
		//   no  -> continue to next, preserve flag

		struct ChainPoint {
			enum {
				OnActive,
				OnNext
			} on;
			uint32_t chain;
			uint32_t vertex;
		};

		std::unordered_set< glm::uvec2 > next_visited;

		auto walk_chains = [&](ChainPoint const &seed) {
			//check if this seed is already used:
			assert(seed.on == ChainPoint::OnNext);
			if (next_visited.count(glm::uvec2(seed.chain, seed.vertex))) return;
			//if not, start walking!
			next_active_chains.emplace_back();
			next_active_flags.emplace_back();
			std::vector< std::vector< ak::EmbeddedVertex > > *chain = &next_active_chains.back();
			std::vector< std::vector< ak::Flag > > *flags = &next_active_flags.back();
			ChainPoint cp = seed;
			while (true) {
				if (cp.on == ChainPoint::OnNext) {
					assert(cp
					chain->emplace_back(next_chains[cp.chain][cp.vertex]);
					flags->emplace_back(
				} else { assert(cp.on == ChainPoint::OnActive);
				}
			}

		};

		for (uint32_t seed_chain = 0; seed_chain < next_chains.size(); ++seed_chain) {
			//a bit of extra handling to avoid seeding on last stitch of non-loop chain:
			uint32_t end = next_chains[seed_chain].size();
			if (next_chains[seed_chain][0] == next_chains[seed_chain].back()) end = next_chains[seed_chain].size()-1;
			for (uint32_t seed_vertex = 0; seed_vertex < end; ++seed_vertex) {
				ak::Flag f = next_flags[seed_chain][seed_vertex];
				if (f == ak::FlagLinkOne || f == ak::FlagLinkAny) {
					ChainPoint cp;
					cp.on = ChainPoint::OnNext;
					cp.chain = seed_chain;
					cp.vertex = seed_vertex;
					walk_chains(cp);
				}
			}
		}


	struct VertexFlag {
		ak::EmbeddedVertex vertex;
		ak::Flag flag;
	};

	for (auto const &chain : next_chains) {
		auto const &flags = next_flags[&chain - &next_chains[0]];

		bool has_discard = false;
		bool has_stitches = false;
		for (auto f : flags) {
			if (f == ak::FlagDiscard) has_discard = true;
			if (f == ak::FlagLinkOne || f == ak::FlagLinkAny) has_stitches = true;
		}
		//don't create new chains from chains that don't contain stitches:
		if (!has_stitches) continue;

		if (chain[0] == chain.back()) { //chain is loop
			//start/end output chain on this stitch:
			uint32_t first = 0;
			while (first < chain.size() && !(flags[first] == ak::FlagLinkOne || flags[first] == ak::FlagLinkAny)) {
				++first;
			}
			//roll so that first 
		} else { //chain is not loop
			std::cerr << "ERROR: code for non-loop chains not yet written." << std::endl;
			continue;
		}
	}

#endif //kill old code

}
