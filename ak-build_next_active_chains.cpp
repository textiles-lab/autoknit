#include "pipeline.hpp"

#include <iostream>
#include <map>
#include <unordered_set>
#include <set>

struct OnChainVertex {
	enum On : uint8_t { OnNone, OnActive, OnNext } on;
	uint32_t chain;
	uint32_t vertex;
	OnChainVertex(On on_ = OnNone, uint32_t chain_ = -1U, uint32_t vertex_ = -1U) : on(on_), chain(chain_), vertex(vertex_) { }

	bool operator<(OnChainVertex const &o) const {
		if (on != o.on) return on < o.on;
		else if (chain != o.chain) return chain < o.chain;
		else return vertex < o.vertex;
	}
	bool operator==(OnChainVertex const &o) const {
		return on == o.on && chain == o.chain && vertex == o.vertex;
	}
	bool operator!=(OnChainVertex const &o) const {
		return !(*this == o);
	}
};

std::ostream &operator<<(std::ostream &out, OnChainVertex const &ocv) {
	if (ocv.on == OnChainVertex::OnActive) out << 'a';
	else if (ocv.on == OnChainVertex::OnNext) out << 'n';
	else if (ocv.on == OnChainVertex::OnNone) out << 'x';
	else out << '?';
	if (ocv.chain == -1U) out << '.';
	else out << ocv.chain;
	out << ':';
	if (ocv.vertex == -1U) out << '.';
	else out << ocv.vertex;
	return out;
};

void ak::build_next_active_chains(
	ak::Parameters const &parameters,
	ak::Model const &model,
	std::vector< std::vector< ak::EmbeddedVertex > > const &active_chains, //in: current active chains
	std::vector< std::vector< ak::Flag > > const &active_flags, //in: flags for current active
	std::vector< std::vector< ak::EmbeddedVertex > > const &next_chains, //in: next chains
	std::vector< std::vector< ak::Flag > > const &next_flags, //in: flags for next active
	std::vector< ak::Link > const &links_in, //in: links between active and next
	std::vector< std::vector< ak::EmbeddedVertex > > *next_active_chains_, //out: next active chains
	std::vector< std::vector< ak::Flag > > *next_active_flags_ //out: next stitch flags
) {
	assert(active_chains.size() == active_flags.size());
	for (uint32_t a = 0; a < active_chains.size(); ++a) {
		assert(active_chains[a].size() == active_flags[a].size());
	}

	assert(next_chains.size() == next_flags.size());
	for (uint32_t a = 0; a < next_chains.size(); ++a) {
		assert(next_chains[a].size() == next_flags[a].size());
	}

	for (auto const &l : links_in) {
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

	//filter to links that target non-discarded stitches only:
	std::vector< ak::Link > links;
	for (auto const &l : links_in) {
		if (next_flags[l.to_chain][l.to_vertex] == ak::FlagDiscard) continue;
		links.emplace_back(l);
	}

	//build a lookup structure for links:
	struct ChainVertex {
		ChainVertex(uint32_t chain_, uint32_t vertex_) : chain(chain_), vertex(vertex_) { }
		uint32_t chain;
		uint32_t vertex;
		bool operator<(ChainVertex const &o) const {
			if (chain != o.chain) return chain < o.chain;
			else return vertex < o.vertex;
		}
		bool operator==(ChainVertex const &o) const {
			return chain == o.chain && vertex == o.vertex;
		}
		bool operator!=(ChainVertex const &o) const {
			return !(*this == o);
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
		for (auto const &flags : next_flags) {
			uint32_t ci = &flags - &next_flags[0];
			bool is_loop = next_chains[ci][0] == next_chains[ci].back();
			for (uint32_t i = 0; i < flags.size(); ++i) {
				if (i + 1 == flags.size() && is_loop) break;
				if (flags[i] == ak::FlagLinkOne || flags[i] == ak::FlagLinkAny) {
					next_stitches[ci].insert(i);
				}
			}
		}
		//stitches are the destination of links (even if flagged discard):
		for (auto const &l : links_in) {
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
		for (auto const &flags : active_flags) {
			uint32_t ci = &flags - &active_flags[0];
			bool is_loop = (active_chains[ci][0] == active_chains[ci].back());
			for (uint32_t i = 0; i < flags.size(); ++i) {
				if (i + 1 == flags.size() && is_loop) break;
				if (flags[i] == ak::FlagLinkOne || flags[i] == ak::FlagLinkAny) {
					active_stitches[ci].insert(i);
				}
			}
		}
		//stitches are the source of links [but should always be flagged anyway]:
		for (auto const &l : links_in) {
			ak::Flag f = active_flags[l.from_chain][l.from_vertex];
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

	//record whether the segments adjacent to every next stitch is are marked as "discard" or "keep":
	std::vector< std::map< uint32_t, std::pair< bool, bool > > > keep_adj(next_chains.size());
	for (uint32_t nc = 0; nc < next_chains.size(); ++nc) {
		auto const &chain = next_chains[nc];
		auto const &flags = next_flags[nc];
		bool is_loop = (chain[0] == chain.back());

		auto &ka = keep_adj[nc];

		uint32_t end = (is_loop ? flags.size() - 1 : flags.size());
		bool have_discard = false;
		uint32_t prev_nv = -1U;

		//record if there is a discard before and after every stitch:
		for (uint32_t nv = 0; nv < end; ++nv) {
			if (flags[nv] == ak::FlagLinkOne || flags[nv] == ak::FlagLinkAny) {
				if (prev_nv != -1U) ka[prev_nv].second = !have_discard;
				ka[nv].first = !have_discard;
				prev_nv = nv;
				have_discard = false;
			} else if (flags[nv] == ak::FlagDiscard) {
				have_discard = true;
			}
		}
		if (prev_nv != -1U) ka[prev_nv].second = !have_discard;

		if (ka.empty()) continue;

		if (ka.size() == 1) {
			std::cout << "WARNING: hit the 'only one retained stitch on a next chain' case, which is rare (impossible?) and probably not tested." << std::endl;
			ka.begin()->second.first = false;
			ka.begin()->second.second = false;
		}

		if (is_loop) {
			//fuse info in first and last stitch:
			ka.begin()->second.first = ka.rbegin()->second.second = (ka.begin()->second.first && ka.rbegin()->second.second);
		} else {
			//always discard edges in chain:
			ka.begin()->second.first = ka.rbegin()->second.second = false;
		}
	}

	//any active stitch without a link marks the segment over it as "discard":
	for (uint32_t ac = 0; ac < active_chains.size(); ++ac) {

		auto discard_after = [&next_chains,&next_flags,&keep_adj](ChainVertex const &cv) {
			assert(cv.chain < next_flags.size());
			assert(cv.vertex < next_flags[cv.chain].size());
			assert(next_flags[cv.chain][cv.vertex] == ak::FlagLinkOne || next_flags[cv.chain][cv.vertex] == ak::FlagLinkAny);
			bool is_loop = (next_chains[cv.chain][0] == next_chains[cv.chain].back());

			auto &ka = keep_adj[cv.chain];
			auto f = ka.find(cv.vertex);
			assert(f != ka.end());
			//mark after cv as not keep:
			f->second.second = false;
			//mark before next stitch as not keep:
			++f;
			if (f == ka.end() && is_loop) f = ka.begin();
			if (f != ka.end()) f->second.first = false;
		};

		auto discard_before = [&next_chains,&next_flags,&keep_adj](ChainVertex const &cv) {
			assert(cv.chain < next_flags.size());
			assert(cv.vertex < next_flags[cv.chain].size());
			assert(next_flags[cv.chain][cv.vertex] == ak::FlagLinkOne || next_flags[cv.chain][cv.vertex] == ak::FlagLinkAny);
			bool is_loop = (next_chains[cv.chain][0] == next_chains[cv.chain].back());

			auto &ka = keep_adj[cv.chain];
			auto f = ka.find(cv.vertex);
			assert(f != ka.end());
			//mark before cv as not keep:
			f->second.first = false;
			//mark after previous stitch as not keep:
			if (f == ka.begin() && is_loop) f = ka.end();
			if (f != ka.begin()) {
				--f;
				f->second.second = false;
			}
		};

		auto do_segment = [&](uint32_t av1, uint32_t av2) {
			auto f1 = active_next.find(ChainVertex(ac, av1));
			auto f2 = active_next.find(ChainVertex(ac, av2));
			if (f1 != active_next.end() && f2 == active_next.end()) {
				// n0 n1
				//  \ /    x
				//   a1 -- a2
				discard_after(f1->second.back());
			} else if (f1 == active_next.end() && f2 != active_next.end()) {
				//      n0 n1
				//  x    \ /   
				// a1 --- a2
				discard_before(f2->second[0]);
			}
		};

		auto const &chain = active_chains[ac];
		auto const &flags = active_flags[ac];
		bool is_loop = (chain[0] == chain.back());

		uint32_t prev_av = -1U;

		if (is_loop) {
			for (uint32_t av = flags.size()-2; av < flags.size(); --av) {
				if (flags[av] == ak::FlagLinkOne || flags[av] == ak::FlagLinkAny) {
					prev_av = av;
					break;
				}
			}
		}
		uint32_t end = (is_loop ? flags.size() - 1 : flags.size());
		for (uint32_t av = 0; av < end; ++av) {
			if (flags[av] == ak::FlagLinkOne || flags[av] == ak::FlagLinkAny) {
				if (prev_av != -1U) do_segment(prev_av, av);
				prev_av = av;
			}
		}
	}


	//Now that discard information is known,
	//segment-stitch links are made as follows:

	// subsequent non-discard segments: [n1,n2] -> n3
	// n1 --> n2 --> n3
	//  |     |
	// a1 --- a2

	// non-discard followed by discard go to active: [n1, n2] -> a3
	// n1 --> n2 -X> n3
	//  |     / \       ...
	// a1 - a2  a3

	//this shouldn't happen, but if it does the logical things is [a2,n2] -> a3
	// n1 -x> n2 -X> n3
	//        / \         ....
	//      a2  a3

	//  non-linked active: [a1, a2] -> a3
	// (?)     x     (?)
	// a1 --> a2 --> a3

	// non-linked to linked active: [a1, a2] -> n1
	//      n1  n2
	//       \  /
	// a1 --> a2 --> a3

	//linked to non-linked active: [n2, a2] -> a3
	//      n1  n2
	//       \  /
	// a1 --> a2 --> a3


	//build a lookup structure for stitches:


	std::map< std::pair< OnChainVertex, OnChainVertex >, OnChainVertex > next_vertex;

	for (uint32_t nc = 0; nc < next_chains.size(); ++nc) {
		auto const &chain = next_chains[nc];
		bool is_loop = (chain[0] == chain.back());
		auto const &ak = keep_adj[nc];

		if (ak.empty()) continue;

		auto prev = ak.end();
		auto cur = ak.end();
		if (is_loop) {
			--cur;
			prev = cur;
			if (prev == ak.begin()) prev = ak.end();
			--prev;
		}
		for (auto next = ak.begin(); next != ak.end(); ++next) {
			//std::cout << "Prev: " << (prev == ak.end() ? 'x' : prev->first)
			//	<< " Cur: " << (cur == ak.end() ? 'x' : cur->first)
			//	<< " Next: " << next->first << std::endl; //DEBUG
			if (cur != ak.end()) {
				OnChainVertex cur_ocv(OnChainVertex::OnNext, nc, cur->first);
				OnChainVertex prev_ocv;
				OnChainVertex next_ocv;
				if (cur->second.first == true) {
					assert(prev != ak.end());
					assert(prev->second.second == true);
					//two stitches with a keep range between them.
					prev_ocv.on = OnChainVertex::OnNext;
					prev_ocv.chain = nc;
					prev_ocv.vertex = prev->first;
				} else { assert(cur->second.first == false);
					//range to prev stitch is not marked 'keep':
					assert(prev == ak.end() || prev->second.second == false);
					//traverse link down to active, if exists:
					auto link = next_active.find(ChainVertex(nc, cur->first));
					if (!link->second.empty()) {
						prev_ocv.on = OnChainVertex::OnActive;
						prev_ocv.chain = link->second[0].chain;
						prev_ocv.vertex = link->second[0].vertex;
					}
				}
				if (cur->second.second) {
					assert(next != ak.end());
					assert(next->second.first == true);
					//have a keep range to the next stitch:
					next_ocv.on = OnChainVertex::OnNext;
					next_ocv.chain = nc;
					next_ocv.vertex = next->first;
				} else { assert(cur->second.second == false);
					//range to next stitch is not marked 'keep':
					assert(next == ak.end() || next->second.first == false);
					//traverse link down to active, if exists:
					auto link = next_active.find(ChainVertex(nc, cur->first));
					if (!link->second.empty()) {
						next_ocv.on = OnChainVertex::OnActive;
						next_ocv.chain = link->second.back().chain;
						next_ocv.vertex = link->second.back().vertex;
					}
				}

				//std::cout << "   " << prev_ocv << " | " << cur_ocv << " | " << next_ocv << std::endl; //DEBUG

				if (prev_ocv.on != OnChainVertex::OnNone && next_ocv.on != OnChainVertex::OnNone) {
					auto ret = next_vertex.insert(std::make_pair(
						std::make_pair(prev_ocv, cur_ocv),
						next_ocv
					));
					//std::cout << "Inserted [" << ret.first->first.first << ", " << ret.first->first.second << "] -> " << ret.first->second << std::endl; //DEBUG
					assert(ret.second);
				}
			}
			
			prev = cur;
			cur = next;
		}

	}


	//now edges from active chains:
	for (uint32_t ac = 0; ac < active_chains.size(); ++ac) {
		auto const &chain = active_chains[ac];
		bool is_loop = (chain[0] == chain.back());

		std::vector< uint32_t > stitches;
		{
			auto const &flags = active_flags[ac];
			uint32_t end = (is_loop ? flags.size() - 1 : flags.size());
			for (uint32_t av = 0; av < end; ++av) {
				if (flags[av] == ak::FlagLinkOne || flags[av] == ak::FlagLinkAny) {
					stitches.emplace_back(av);
				}
			}
		}

		if (stitches.empty()) continue;

		auto prev = stitches.end();
		auto cur = stitches.end();
		if (is_loop) {
			--cur;
			prev = cur;
			if (prev == stitches.begin()) prev = stitches.end();
			--prev;
		}

		for (auto next = stitches.begin(); next != stitches.end(); ++next) {
			if (cur != stitches.end()) {
				OnChainVertex cur_ocv(OnChainVertex::OnActive, ac, *cur);
				auto cur_link = active_next.find(ChainVertex(ac, *cur));
				if (cur_link == active_next.end()) {
					//no link, so edge is previous to next:
					//      x    
					// p -> c -> n
					if (prev != stitches.end() && next != stitches.end()) {
						auto ret = next_vertex.insert(std::make_pair(
							std::make_pair(OnChainVertex(OnChainVertex::OnActive, ac, *prev), cur_ocv),
							OnChainVertex(OnChainVertex::OnActive, ac, *next)
						));
						assert(ret.second);
					}
				} else {
					//have a link.
					// -x> n0 n1 -x>
					//      \ /
					// p --> c --> n

					//link previous up (if not under a keep):
					if (prev != stitches.end()) {
						auto n = cur_link->second[0];
						//don't link up if this is right of another stitch that is also connected:
						// -x> n0
						//    /  \    x
						//   p -> c
						auto n_link = next_active.find(n);
						assert(n_link != next_active.end());
						if (n_link->second[0] == ChainVertex(ac, *cur)) {
							auto f = keep_adj[n.chain].find(n.vertex);
							assert(f != keep_adj[n.chain].end());
							if (f->second.first == false) {
								auto ret = next_vertex.insert(std::make_pair(
									std::make_pair(OnChainVertex(OnChainVertex::OnActive, ac, *prev), cur_ocv),
									OnChainVertex(OnChainVertex::OnNext, n.chain, n.vertex)
								));
								assert(ret.second);
							}
						}
					}

					//link down to next (if not under a keep):
					if (next != stitches.end()) {
						auto n = cur_link->second.back();
						//don't link down if this is left of another stitch that is also connected:
						//     n0 -x>
						//    /  \    x
						//   c -> n
						auto n_link = next_active.find(n);
						assert(n_link != next_active.end());
						if (n_link->second.back() == ChainVertex(ac, *cur)) {
							auto f = keep_adj[n.chain].find(n.vertex);
							assert(f != keep_adj[n.chain].end());
							if (f->second.second == false) {
								auto ret = next_vertex.insert(std::make_pair(
									std::make_pair(OnChainVertex(OnChainVertex::OnNext, n.chain, n.vertex), cur_ocv),
									OnChainVertex(OnChainVertex::OnActive, ac, *next)
								));
								assert(ret.second);
							}
						}
					}

					//deal with the (rare? impossible?) case where have edges to a discard segment:
					if (cur_link->second.size() == 2) {
						auto n0 = cur_link->second[0];
						auto n1 = cur_link->second.back();
						auto f0 = keep_adj[n0.chain].find(n0.vertex);
						assert(f0 != keep_adj[n0.chain].end());
						auto f1 = keep_adj[n1.chain].find(n1.vertex);
						assert(f1 != keep_adj[n1.chain].end());
						if (f0->second.second == false || f1->second.first == false) {
							assert(f0->second.second == false && f1->second.first == false);
							std::cerr << "WARNING: encountered very odd (impossible?) case where chain dips to active for one [increase] stitch." << std::endl;
							auto ret = next_vertex.insert(std::make_pair(
								std::make_pair(OnChainVertex(OnChainVertex::OnNext, n0.chain, n0.vertex), cur_ocv),
								OnChainVertex(OnChainVertex::OnNext, n1.chain, n1.vertex)
							));
							assert(ret.second);
						}
					}
				}
			}
			
			prev = cur;
			cur = next;
		}

	}



#if 0
	//DEBUG:
	for (uint32_t c = 0; c < next_chains.size(); ++c) {
		std::cout << "next[" << c << "] ";
		for (auto f : next_flags[c]) {
			if      (f == ak::FlagDiscard)  std::cout << '-';
			else if (f == ak::FlagLinkNone) std::cout << '.';
			else if (f == ak::FlagLinkOne)  std::cout << '1';
			else if (f == ak::FlagLinkAny)  std::cout << '2';
			else std::cout << '?';
		}
		std::cout << std::endl;
	}

	//DEBUG:
	for (uint32_t c = 0; c < active_chains.size(); ++c) {
		std::cout << "active[" << c << "] ";
		for (auto f : active_flags[c]) {
			if      (f == ak::FlagDiscard)  std::cout << '-';
			else if (f == ak::FlagLinkNone) std::cout << '.';
			else if (f == ak::FlagLinkOne)  std::cout << '1';
			else if (f == ak::FlagLinkAny)  std::cout << '2';
			else std::cout << '?';
		}
		std::cout << std::endl;
	}

	//DEBUG:
	for (auto const &cvv : active_next) {
		std::cout << "active[" << cvv.first.chain << "][" << cvv.first.vertex << "] ->";
		for (auto const &s : cvv.second) {
			std::cout << " next[" << s.chain << "][" << s.vertex << "]";
		}
		std::cout << "\n";
	}
	std::cout.flush();
#endif


	//Walk through created edges array, creating chains therefrom:

	std::vector< std::vector< OnChainVertex > > loops;
	std::map< std::pair< OnChainVertex, OnChainVertex >, std::vector< OnChainVertex > > partials;

	while (!next_vertex.empty()) {
		std::vector< OnChainVertex > chain;
		chain.emplace_back(next_vertex.begin()->first.first);
		chain.emplace_back(next_vertex.begin()->first.second);
		chain.emplace_back(next_vertex.begin()->second);
		assert(chain[0] != chain[1] && chain[0] != chain[2] && chain[1] != chain[2]);
		next_vertex.erase(next_vertex.begin());
		while (true) {
			auto f = next_vertex.find(std::make_pair(chain[chain.size()-2], chain[chain.size()-1]));
			if (f == next_vertex.end()) break;
			chain.emplace_back(f->second);
			next_vertex.erase(f);
		}

		{ //check if a partial chain comes after this one; if so, append it:
			auto f = partials.find(std::make_pair(chain[chain.size()-2], chain[chain.size()-1]));
			if (f != partials.end()) {
				chain.pop_back();
				chain.pop_back();
				chain.insert(chain.end(), f->second.begin(), f->second.end());
				partials.erase(f);
			}
		}

		//loops should look like abcdab
		//because abcd -> ab-c bc-d cd-a da-b
		if (chain[0] == chain[chain.size()-2] && chain[1] == chain[chain.size()-1]) {
			//great -- full loop.
			chain.pop_back();
			loops.emplace_back(chain);
		} else {
			//partial loop -- save for later
			auto ret = partials.insert(std::make_pair(std::make_pair(chain[0], chain[1]), chain));
			assert(ret.second);
		}
	}

	auto output = [&](std::vector< OnChainVertex > const &ocvs) {
		std::vector< ak::EmbeddedVertex > chain;
		std::vector< ak::Flag > flags;
		for (uint32_t i = 0; i + 1 < ocvs.size(); ++i) {
			auto const &ocv0 = ocvs[i];
			auto const &ocv1 = ocvs[i+1];

			std::vector< ak::EmbeddedVertex > const &src_chain = (ocv0.on == OnChainVertex::OnActive ? active_chains : next_chains)[ocv0.chain];
			std::vector< ak::Flag > const &src_flags = (ocv0.on == OnChainVertex::OnActive ? active_flags : next_flags)[ocv0.chain];
			if (i == 0) {
				chain.emplace_back(src_chain[ocv0.vertex]);
				flags.emplace_back(src_flags[ocv0.vertex]);
			}

			if (ocv0.on == ocv1.on) {
				assert(ocv0.chain == ocv1.chain);
				if (ocv0.vertex < ocv1.vertex) {
					for (uint32_t j = ocv0.vertex + 1; j <= ocv1.vertex; ++j) {
						chain.emplace_back(src_chain[j]);
						flags.emplace_back(src_flags[j]);
					}
				} else { assert(ocv1.vertex < ocv0.vertex);
					assert(src_chain[0] == src_chain.back());
					for (uint32_t j = ocv0.vertex + 1; j + 1 < src_chain.size(); ++j) {
						chain.emplace_back(src_chain[j]);
						flags.emplace_back(src_flags[j]);
					}
					for (uint32_t j = 0; j <= ocv1.vertex; ++j) {
						chain.emplace_back(src_chain[j]);
						flags.emplace_back(src_flags[j]);
					}
				}
			} else {
				std::vector< ak::EmbeddedVertex > const &dst_chain = (ocv1.on == OnChainVertex::OnActive ? active_chains : next_chains)[ocv1.chain];
				std::vector< ak::Flag > const &dst_flags = (ocv1.on == OnChainVertex::OnActive ? active_flags : next_flags)[ocv1.chain];

				std::vector< ak::EmbeddedVertex > path;

				ak::embedded_path(parameters, model,
					src_chain[ocv0.vertex],
					dst_chain[ocv1.vertex],
					&path);

				assert(path.size() >= 2);
				assert(path[0] == src_chain[ocv0.vertex]);
				assert(path.back() == dst_chain[ocv1.vertex]);

				assert(flags.size() == chain.size());
				if (ocv0.on == OnChainVertex::OnActive) {
					assert(chain.back() == src_chain[ocv0.vertex]);
					flags.back() = ak::FlagLinkNone;
				}

				chain.insert(chain.end(), path.begin() + 1, path.end());
				while (flags.size() + 1 < chain.size()) {
					flags.emplace_back(ak::FlagLinkNone);
				}
				flags.emplace_back(dst_flags[ocv1.vertex]);

				assert(flags.size() == chain.size());

				if (ocv1.on == OnChainVertex::OnActive) {
					assert(chain.back() == dst_chain[ocv1.vertex]);
					flags.back() = ak::FlagLinkNone;
				}

			}
		}
		if (chain[0] == chain.back()) {
			if (flags[0] == ak::FlagLinkNone || flags.back() == ak::FlagLinkNone) {
				flags[0] = flags.back() = ak::FlagLinkNone;
			}
		}

		next_active_chains.emplace_back(chain);
		next_active_flags.emplace_back(flags);
	};

	std::cout << "Found " << loops.size() << " loops and " << partials.size() << " chains." << std::endl;
	for (auto const &loop : loops) {
		output(loop);
	}
	for (auto const &pp : partials) {
		output(pp.second);
	}



/* for later:
	ak::embedded_path(parameters, model,
		next_chains[frag.chain][frag.inds.back()],
		active_chains[to.chain][to.vertex],
		&path);
*/

}
