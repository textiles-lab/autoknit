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
		if (active_chains[l.from_chain][0] == active_chains[l.from_chain].back()) {
			assert(l.from_vertex + 1 < active_chains[l.from_chain].size());
		}
		assert(l.to_chain < next_chains.size());
		assert(l.to_vertex < next_chains[l.to_chain].size());
		if (next_chains[l.to_chain][0] == next_chains[l.to_chain].back()) {
			assert(l.to_vertex + 1 < next_chains[l.to_chain].size());
		}
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

	std::map< ChainVertex, std::vector< ChainVertex > > active_next;
	std::map< ChainVertex, std::vector< ChainVertex > > next_active;

	for (auto const &l : links) {
		active_next[ ChainVertex(l.from_chain, l.from_vertex) ]
			.emplace_back(l.to_chain, l.to_vertex);
		next_active[ ChainVertex(l.to_chain, l.to_vertex) ]
			.emplace_back(l.from_chain, l.from_vertex);
	}

	//now sort the links so that they run in increasing chain order:
	{ //sort active->next links:
		//need to know where stitches are to check:
		std::vector< std::unordered_set< uint32_t > > next_stitches(next_chains.size());
		//stitches are flagged stitches:
		for (auto const &flags : next_flags_in) {
			uint32_t ci = &flags - &next_flags_in[0];
			bool is_loop = next_chains[ci][0] == next_chains[ci].back();
			for (uint32_t i = 0; i < flags.size(); ++i) {
				if (i + 1 == flags.size() && is_loop) break;
				if (flags[i] == ak::FlagLinkOne || flags[i] == ak::FlagLinkAny) {
					next_stitches[ci].insert(i);
				}
			}
		}
		//stitches are the destination of links (even if flagged discard):
		for (auto const &l : links) {
			next_stitches[l.to_chain].insert(l.to_vertex);
		}

		for (auto &fv : active_next) {
			std::vector< ChainVertex > &v = fv.second;
			assert(!v.empty());
			if (v.size() == 1) continue;
			assert(v.size() == 2);
			assert(v[0].chain == v[1].chain); //no splits/merges
			assert(v[0].vertex != v[1].vertex);
			uint32_t ci = v[0].chain;
			auto stitch_in = [&](uint32_t begin, uint32_t end) {
				for (uint32_t i = begin; i < end; ++i) {
					if (next_stitches[ci].count(i)) return true;
				}
				return false;
			};
			bool is_loop = (next_chains[ci][0] == next_chains[ci].back());
			if (is_loop) {
				if (v[0].vertex < v[1].vertex) {
					if (stitch_in(v[0].vertex+1, v[1].vertex)) std::swap(v[0], v[1]);
				}

				if (v[0].vertex < v[1].vertex) {
					assert(!stitch_in(v[0].vertex+1,v[1].vertex));
					assert(stitch_in(0, v[0].vertex) || stitch_in(v[1].vertex+1, next_chains[ci].size()));
				} else { assert(v[0].vertex > v[1].vertex);
					assert(!stitch_in(v[0].vertex+1, next_chains[ci].size()) && !stitch_in(0, v[1].vertex));
					assert(stitch_in(v[1].vertex+1,v[0].vertex));
				}
			} else {
				if (v[0].vertex > v[1].vertex) std::swap(v[0], v[1]);
				assert(!stitch_in(v[0].vertex+1,v[1].vertex));
			}
		}
	}

	//now sort the links so that they run in increasing chain order:
	{ //sort next->active links:
		//need to know where stitches are to check:
		std::vector< std::unordered_set< uint32_t > > active_stitches(active_chains.size());
		//stitches are flagged stitches:
		for (auto const &flags : active_flags_in) {
			uint32_t ci = &flags - &active_flags_in[0];
			bool is_loop = (active_chains[ci][0] == active_chains[ci].back());
			for (uint32_t i = 0; i < flags.size(); ++i) {
				if (i + 1 == flags.size() && is_loop) break;
				if (flags[i] == ak::FlagLinkOne || flags[i] == ak::FlagLinkAny) {
					active_stitches[ci].insert(i);
				}
			}
		}
		//stitches are the source of links [but should always be flagged anyway]:
		for (auto const &l : links) {
			ak::Flag f = active_flags_in[l.from_chain][l.from_vertex];
			assert(f == ak::FlagLinkOne || f == ak::FlagLinkAny);
			assert(active_stitches[l.from_chain].count(l.from_vertex));
		}

		for (auto &fv : next_active) {
			std::vector< ChainVertex > &v = fv.second;
			assert(!v.empty());
			if (v.size() == 1) continue;
			assert(v.size() == 2);
			assert(v[0].chain == v[1].chain); //no splits/merges
			assert(v[0].vertex != v[1].vertex);
			uint32_t ci = v[0].chain;
			auto stitch_in = [&](uint32_t begin, uint32_t end) {
				for (uint32_t i = begin; i < end; ++i) {
					if (active_stitches[ci].count(i)) return true;
				}
				return false;
			};
			bool is_loop = (active_chains[ci][0] == active_chains[ci].back());
			if (is_loop) {
				if (v[0].vertex < v[1].vertex) {
					if (stitch_in(v[0].vertex+1, v[1].vertex)) std::swap(v[0], v[1]);
				}

				if (v[0].vertex < v[1].vertex) {
					assert(!stitch_in(v[0].vertex+1,v[1].vertex));
					assert(stitch_in(0, v[0].vertex) || stitch_in(v[1].vertex+1, next_chains[ci].size()));
				} else { assert(v[0].vertex > v[1].vertex);
					assert(!stitch_in(v[0].vertex+1, next_chains[ci].size()) && !stitch_in(0, v[1].vertex));
					assert(stitch_in(v[1].vertex+1,v[0].vertex));
				}
			} else {
				if (v[0].vertex > v[1].vertex) std::swap(v[0], v[1]);
				assert(!stitch_in(v[0].vertex+1,v[1].vertex));
			}
		}
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
	std::vector< std::unordered_set< uint32_t > > next_discard_after(next_chains.size());
	std::vector< std::unordered_set< uint32_t > > next_discard_before(next_chains.size());
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
			bool discard_before = ((i > 0 ? flags[i-1] : before_front) == ak::FlagDiscard);
			bool discard_after = (i + 1 < flags.size() && flags[i+1] == ak::FlagDiscard);
			if (discard_before) next_discard_before[c].insert(i);
			if (discard_after) next_discard_after[c].insert(i);
			if (discard_before || discard_after) {
				assert(flags[i] == ak::FlagLinkOne); //edges are marked LinkOne always
				next_edges[c].insert(i);
			}
			assert(!(discard_before && discard_after));
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

			ak::Flag ignore = ak::FlagDiscard;
			ak::Flag *prev = (i > 0 ? &flags[i-1] : (is_loop ? &flags[flags.size()-2] : &ignore));
			ak::Flag *next = (is_loop ? &flags[(i+1)%(flags.size()-1)] : (i + 1 < flags.size() ? &flags[i+1] : &ignore));

			//mark any stitches linked to non-edges for discard:
			bool have_edge = false;
			bool have_link = false;
			auto f = active_next.find(ChainVertex(c,i));
			if (f != active_next.end()) {
				auto &v = f->second;
				assert(!v.empty());
				for (auto &cv : v) {
					ak::Flag f = next_flags[cv.chain][cv.vertex];
					if (f == ak::FlagDiscard) {
						//ignore.
					} else if (next_edges[cv.chain].count(cv.vertex)) {
						assert(f == ak::FlagLinkOne); //edges are always marked LinkOne
						have_edge = true;
					} else {
						assert(f == ak::FlagLinkOne || f == ak::FlagLinkAny); //this should be a link to a valid-flagged stitch
						have_link = true;
					}
				}

				//add more discard marks if this is connected to a fragment-ending stitch:
				if (next_discard_before[v[0].chain].count(v[0].vertex)) {
					if (*next != ak::FlagDiscard) {
						assert(*next == ak::FlagLinkNone); //requiring at least one location in each chain *between* stitches.
						*next = ak::FlagDiscard;
					}
				}

				if (next_discard_after[v.back().chain].count(v.back().vertex)) {
					if (*prev != ak::FlagDiscard) {
						assert(*prev == ak::FlagLinkNone); //requiring at least one location in each chain *between* stitches.
						*prev = ak::FlagDiscard;
					}
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


	/*

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

	//DEBUG: dump paths as individual chains as well
	for (auto const &frag : fragments) {
		if (frag.on == Fragment::OnNext) {
			auto r = next_active.equal_range(ChainVertex(frag.chain, frag.inds.back()));
			//expecting exactly one link, given that this is the start/end of a range:
			assert(r.first != r.second);
			auto t = r.first;
			++t;
			assert(t == r.second);

			ak::EmbeddedVertex from = next_chains[frag.chain][frag.inds.back()];
			ak::EmbeddedVertex to = active_chains[r.first->second.chain][r.first->second.vertex];

			std::vector< ak::EmbeddedVertex > path;

			ak::embedded_path(parameters, model, from, to, &path);

			next_active_chains.emplace_back(path);
			next_active_flags.emplace_back(next_active_chains.back().size(), ak::FlagLinkNone);
			
		} else {
			assert(frag.on == Fragment::OnActive);
		}
	}
	*/

	std::vector< bool > used(fragments.size(), false);
	for (uint32_t seed = 0; seed < fragments.size(); ++seed) {
		if (used[seed]) continue;
		uint32_t at = seed;
		next_active_chains.emplace_back();
		next_active_flags.emplace_back();
		auto &chain = next_active_chains.back();
		auto &flags = next_active_flags.back();
		while (true) {
			assert(!used[at]);
			used[at] = true;
			//do fragment:
			auto const &frag = fragments[at];
			auto const &src_chain = (frag.on == Fragment::OnNext ? next_chains : active_chains)[frag.chain];
			auto const &src_flags = (frag.on == Fragment::OnNext ? next_flags : active_flags)[frag.chain];
			for (auto i : frag.inds) {
				chain.emplace_back(src_chain[i]);
				flags.emplace_back(src_flags[i]);
			}
			//do connection to next fragment:
			std::vector< ak::EmbeddedVertex > path;
			if (frag.on == Fragment::OnNext) {
				auto f = next_active.find(ChainVertex(frag.chain, frag.inds.back()));
				//expecting exactly one link, given that this is the start/end of a range:
				assert(f != next_active.end());
				assert(f->second.size() == 1);

				ChainVertex to = f->second[0];

				if (active_flags[to.chain][to.vertex] == ak::FlagDiscard) {
					//this is the end of a strip
					std::cerr << "End of a strip; but we don't really have strip-backward-extending code written." << std::endl;
					assert(false);
				}

				ak::embedded_path(parameters, model,
					next_chains[frag.chain][frag.inds.back()],
					active_chains[to.chain][to.vertex],
					&path);

				//now find fragment that starts with to.chain / to.vertex:
				uint32_t found = -1U;
				for (uint32_t f = 0; f < fragments.size(); ++f) {
					if (fragments[f].on == Fragment::OnActive && fragments[f].chain == to.chain && fragments[f].inds[0] == to.vertex) {
						assert(found == -1U);
						found = f;
					}
				}
				assert(found != -1U);
				at = found;

			} else { assert(frag.on == Fragment::OnActive);

				auto f = active_next.find(ChainVertex(frag.chain, frag.inds.back()));

				if (f == active_next.end()) {
					std::cerr << "End of a strip; but we don't really have strip-backward-extending code written." << std::endl;
					assert(false);
				}
				assert(!f->second.empty());

				//figure out which exactly one edge this vertex links to (it shouldn't link to two because that would imply two edges without a discard between, which implies a split/merge, which means no increases are allowed)
				ChainVertex to(f->second[0]);
				assert(f->second.size() == 1 || f->second.size() == 2);
				if (f->second.size() == 2) {
					if (!next_edges[to.chain].count(to.vertex)) {
						to = f->second[1];
						assert(next_edges[to.chain].count(to.vertex));
					}
				}
				assert(to.vertex != -1U);
				assert(next_edges[to.chain].count(to.vertex));

				ak::embedded_path(parameters, model,
					active_chains[frag.chain][frag.inds.back()],
					next_chains[to.chain][to.vertex],
					&path);

				//now find fragment that starts with to.chain / to.vertex:
				uint32_t found = -1U;
				for (uint32_t f = 0; f < fragments.size(); ++f) {
					if (fragments[f].on == Fragment::OnNext && fragments[f].chain == to.chain && fragments[f].inds[0] == to.vertex) {
						assert(found == -1U);
						found = f;
					}
				}
				assert(found != -1U);
				at = found;
			}

			for (uint32_t i = 1; i + 1 < path.size(); ++i) {
				chain.emplace_back(path[i]);
				flags.emplace_back(ak::FlagLinkNone);
			}

			if (at == seed) {
				//loop closed.
				assert(chain[0] == path.back());
				chain.emplace_back(chain[0]);
				flags.emplace_back(flags[0]);
				break;
			}
		}

	}

	uint32_t lines = 0;
	uint32_t loops = 0;
	for (auto const &chain : next_active_chains) {
		auto const &flags = next_active_flags[&chain - &next_active_chains[0]];
		if (chain[0] == chain.back()) {
			assert(flags[0] == flags.back());
			++loops;
		} else {
			++lines;
		}
	}

	std::cout << "Have " << next_active_chains.size() << " next active chains; " << lines << " lines and " << loops << " loops." << std::endl;


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
