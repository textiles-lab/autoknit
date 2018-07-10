#include "pipeline.hpp"

#include <iostream>
#include <map>
#include <unordered_set>
#include <set>
#include <algorithm>

struct OnChainStitch {
	enum On : uint8_t { OnNone, OnActive, OnNext } on;
	uint32_t chain;
	uint32_t stitch;
	enum Type : uint8_t { TypeBegin, TypeStitch, TypeEnd } type = TypeStitch;
	OnChainStitch(On on_ = OnNone, uint32_t chain_ = -1U, uint32_t stitch_ = -1U) : on(on_), chain(chain_), stitch(stitch_) { }

	bool operator<(OnChainStitch const &o) const {
		if (on != o.on) return on < o.on;
		else if (chain != o.chain) return chain < o.chain;
		else return stitch < o.stitch;
	}
	bool operator==(OnChainStitch const &o) const {
		return on == o.on && chain == o.chain && stitch == o.stitch;
	}
	bool operator!=(OnChainStitch const &o) const {
		return !(*this == o);
	}
};

std::ostream &operator<<(std::ostream &out, OnChainStitch const &ocv) {
	if (ocv.on == OnChainStitch::OnActive) out << 'a';
	else if (ocv.on == OnChainStitch::OnNext) out << 'n';
	else if (ocv.on == OnChainStitch::OnNone) out << 'x';
	else out << '?';
	if (ocv.chain == -1U) out << '.';
	else out << ocv.chain;
	out << ':';
	if (ocv.stitch == -1U) {
		if (ocv.type == OnChainStitch::TypeStitch) out << '!';
		else if (ocv.type == OnChainStitch::TypeBegin) out << 'B';
		else if (ocv.type == OnChainStitch::TypeEnd) out << 'E';
		else out << 'q';
	} else out << ocv.stitch;
	return out;
};


void ak::build_next_active_chains(
	ak::Parameters const &parameters,
	ak::Model const &slice,
	std::vector< ak::EmbeddedVertex > const &slice_on_model, //in: vertices of slice (on model)
	std::vector< std::vector< uint32_t > > const &active_chains,  //in: current active chains (on slice)
	std::vector< std::vector< ak::Stitch > > const &active_stitches, //in: current active stitches
	std::vector< std::vector< uint32_t > > const &next_chains, //in: next chains (on slice)
	std::vector< std::vector< ak::Stitch > > const &next_stitches, //in: next stitches
	std::vector< ak::Link > const &links_in, //in: links between active and next
	std::vector< std::vector< ak::EmbeddedVertex > > *next_active_chains_, //out: next active chains (on model)
	std::vector< std::vector< ak::Stitch > > *next_active_stitches_ //out: next active stitches
) {
	for (auto const &chain : active_chains) {
		for (auto v : chain) {
			assert(v < slice.vertices.size());
		}
		for (uint32_t i = 1; i < chain.size(); ++i) {
			assert(chain[i-1] != chain[i]);
		}
	}

	assert(active_stitches.size() == active_chains.size());

	for (auto const &chain : next_chains) {
		for (auto v : chain) {
			assert(v < slice.vertices.size());
		}
		for (uint32_t i = 1; i < chain.size(); ++i) {
			assert(chain[i-1] != chain[i]);
		}
	}

	assert(next_stitches.size() == next_chains.size());

	for (auto const &l : links_in) {
		assert(l.from_chain < active_chains.size());
		assert(l.from_stitch < active_stitches[l.from_chain].size());
		assert(l.to_chain < next_chains.size());
		assert(l.to_stitch < next_stitches[l.to_chain].size());
	}

	assert(next_active_chains_);
	auto &next_active_chains = *next_active_chains_;
	next_active_chains.clear();

	assert(next_active_stitches_);
	auto &next_active_stitches = *next_active_stitches_;
	next_active_stitches.clear();

	//filter to links that target non-discarded stitches only:
	std::vector< ak::Link > links;
	for (auto const &l : links_in) {
		if (next_stitches[l.to_chain][l.to_stitch].flag == ak::Stitch::FlagDiscard) continue;
		links.emplace_back(l);
	}

	//build a lookup structure for links:
	struct ChainStitch {
		ChainStitch(uint32_t chain_, uint32_t stitch_) : chain(chain_), stitch(stitch_) { }
		uint32_t chain;
		uint32_t stitch;
		bool operator<(ChainStitch const &o) const {
			if (chain != o.chain) return chain < o.chain;
			else return stitch < o.stitch;
		}
		bool operator==(ChainStitch const &o) const {
			return chain == o.chain && stitch == o.stitch;
		}
		bool operator!=(ChainStitch const &o) const {
			return !(*this == o);
		}
	};

	std::map< ChainStitch, std::vector< ChainStitch > > active_next;
	std::map< ChainStitch, std::vector< ChainStitch > > next_active;

	//NOTE: link_chains guarantees that links are in "direction of chain" order, so old sorting code removed:
	for (auto const &l : links) {
		active_next[ ChainStitch(l.from_chain, l.from_stitch) ]
			.emplace_back(l.to_chain, l.to_stitch);
		next_active[ ChainStitch(l.to_chain, l.to_stitch) ]
			.emplace_back(l.from_chain, l.from_stitch);
	}


	//record whether the segments adjacent to every next stitch is are marked as "discard" or "keep":
	std::vector< std::vector< std::pair< bool, bool > > > keep_adj(next_chains.size());
	for (uint32_t nc = 0; nc < next_chains.size(); ++nc) {
		auto const &chain = next_chains[nc];
		bool is_loop = (chain[0] == chain.back());
		auto const &stitches = next_stitches[nc];
		auto &ka = keep_adj[nc];
		ka.assign(stitches.size(), std::make_pair(true, true));
		for (uint32_t ns = 0; ns < stitches.size(); ++ns) {
			bool discard_before = false;
			bool discard_after = false;
			//This first ("easy") case is needed because of the case:
			// n0 --- x --- x --- n1
			//  \    /       \   /
			//    a0  -------  a1
			// which should get marked as
			// n0 xxx  xxx  xxx n1
			//   \             /
			//    a0  ------ a1
			// (otherwise, could not prune links_in -> links and probably skip this case.)
			if (stitches[ns].flag == ak::Stitch::FlagDiscard) {
				discard_before = true;
				discard_after = true;
			}
			//stitches with a link to an active stitch whose next stitch doesn't have a link are marked discard-adj:
			discard_before = discard_before || [&]() -> bool {
				//okay, which active stitch is linked to this?
				auto fa = next_active.find(ChainStitch(nc, ns));
				if (fa == next_active.end()) return false; //nothing? no reason to discard before (though weird, I guess)
				//something? take earlier something:
				ChainStitch a = fa->second[0];
				auto const &a_chain = active_chains[a.chain];
				bool a_is_loop = (a_chain[0] == a_chain.back());
				auto const &a_stitches = active_stitches[a.chain];
				assert(a.stitch < a_stitches.size());
				{ //check for the case where a has an earlier non-discard link:
					auto fn = active_next.find(a);
					assert(fn != active_next.end());
					assert(fn->second.size() == 1 || fn->second.size() == 2);
					if (fn->second.size() == 2 && fn->second.back() == ChainStitch(nc, ns)) {
						// p -- n
						//  \  /
						//   a
						return false;
					} else {
						assert(fn->second[0] == ChainStitch(nc, ns));
					}
				}
				//check previous stitch:
				if (a.stitch == 0 && !a_is_loop) return false; //no previous stitch
				{
					ChainStitch pa(a.chain, (a.stitch > 0 ? a.stitch - 1 : a_stitches.size() - 1));
					auto fn = active_next.find(pa);
					if (fn == active_next.end()) {
						//        n  
						//   x    | /
						//  pa -- a
						return true; //aha! unlinked stitch; mark segment as not-keep.
					}
				}
				return false; //didn't have unlinked previous stitch
			}();
			discard_after = discard_after || [&]() -> bool {
				//okay, which active stitch is linked to this?
				auto fa = next_active.find(ChainStitch(nc, ns));
				if (fa == next_active.end()) return false; //nothing? no reason to discard after (though weird, I guess)
				//something? take later something:
				ChainStitch a = fa->second.back();
				auto const &a_chain = active_chains[a.chain];
				bool a_is_loop = (a_chain[0] == a_chain.back());
				auto const &a_stitches = active_stitches[a.chain];
				assert(a.stitch < a_stitches.size());
				{ //check for the case where a has a later non-discard link:
					auto fn = active_next.find(a);
					assert(fn != active_next.end());
					assert(fn->second.size() == 1 || fn->second.size() == 2);
					if (fn->second.size() == 2 && fn->second[0] == ChainStitch(nc, ns)) {
						// n -- nn
						//  \  /
						//   a
						return false;
					} else {
						assert((fn->second.size() == 1 && fn->second[0] == ChainStitch(nc, ns))
							|| (fn->second.size() == 2 && fn->second.back() == ChainStitch(nc, ns)) );
					}
				}
				//check next stitch:
				if (a.stitch + 1 == a_stitches.size() && !a_is_loop) return false; //no next stitch
				{
					ChainStitch na(a.chain, (a.stitch + 1 < a_stitches.size() ? a.stitch + 1 : 0));
					auto fn = active_next.find(na);
					if (fn == active_next.end()) {
						//   n  
						// \ |     x
						//   a --- na
						return true; //aha! unlinked stitch; mark segment as not-keep.
					}
				}
				return false; //didn't have an unlinked next stitch
			}();

			if (discard_before) {
				if (ns > 0) ka[ns-1].second = false;
				else if (is_loop) ka.back().second = false;
				ka[ns].first = false;
			}

			if (discard_after) {
				ka[ns].second = false;
				if (ns + 1 < stitches.size()) ka[ns+1].first = false;
				else if (is_loop) ka[0].first = false;
			}
		}
	}

	//TODO: dump all links and keeps in glorious ascii-art for "debugging purposes"

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

	//this shouldn't happen, but if it does the logical thing is [a2,n2] -> a3
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

	std::map< std::pair< OnChainStitch, OnChainStitch >, OnChainStitch > next_vertex;

	//edges where middle vertex is on a next chain:
	for (uint32_t nc = 0; nc < next_chains.size(); ++nc) {
		auto const &chain = next_chains[nc];
		bool is_loop = (chain[0] == chain.back());
		auto const &ak = keep_adj[nc];
		auto const &stitches = next_stitches[nc];

		assert(ak.size() == stitches.size());

		for (uint32_t ns = 0; ns < stitches.size(); ++ns) {
			if (stitches[ns].flag == ak::Stitch::FlagDiscard) continue;
			OnChainStitch cur_ocs(OnChainStitch::OnNext, nc, ns);
			OnChainStitch prev_ocs;
			if (ak[ns].first) {
				//keep before segment, so prev is previous stitch (or to the begin of a non-loop):
				if (ns > 0 || is_loop) {
					prev_ocs.on = OnChainStitch::OnNext;
					prev_ocs.chain = nc;
					prev_ocs.stitch = (ns > 0 ? ns - 1 : stitches.size()-1);
				} else {
					prev_ocs.on = OnChainStitch::OnNext;
					prev_ocs.chain = nc;
					prev_ocs.stitch = -1U;
					prev_ocs.type = OnChainStitch::TypeBegin;
				}
			} else {
				//don't keep before segment, so prev involves walking down a link:
				auto fa = next_active.find(ChainStitch(nc, ns));
				assert(fa != next_active.end()); //there should always be a link to walk down, right?
				prev_ocs.on = OnChainStitch::OnActive;
				prev_ocs.chain = fa->second[0].chain;
				prev_ocs.stitch = fa->second[0].stitch;
			}
			OnChainStitch next_ocs;
			if (ak[ns].second) {
				//keep after segment, so next is next stitch (or to the end of a non-loop):
				if (ns + 1 < stitches.size() || is_loop) {
					next_ocs.on = OnChainStitch::OnNext;
					next_ocs.chain = nc;
					next_ocs.stitch = (ns + 1 < stitches.size() ? ns + 1 : 0);
				} else {
					next_ocs.on = OnChainStitch::OnNext;
					next_ocs.chain = nc;
					next_ocs.stitch = -1U;
					next_ocs.type = OnChainStitch::TypeEnd;
				}
			} else {
				//don't keep next segment, so next involves walking down a link:
				auto fa = next_active.find(ChainStitch(nc, ns));
				assert(fa != next_active.end()); //there should always be a link to walk down, right?
				next_ocs.on = OnChainStitch::OnActive;
				next_ocs.chain = fa->second.back().chain;
				next_ocs.stitch = fa->second.back().stitch;
			}

			if (prev_ocs.on != OnChainStitch::OnNone && next_ocs.on != OnChainStitch::OnNone) {
				auto ret = next_vertex.insert(std::make_pair(
					std::make_pair(prev_ocs, cur_ocs),
					next_ocs
				));
				std::cout << "Inserted [" << ret.first->first.first << ", " << ret.first->first.second << "] -> " << ret.first->second << std::endl; //DEBUG
				assert(ret.second);
			}
		}
	}


	//now edges from active chains:
	for (uint32_t ac = 0; ac < active_chains.size(); ++ac) {
		auto const &chain = active_chains[ac];
		bool is_loop = (chain[0] == chain.back());
		auto const &stitches = active_stitches[ac];

		for (uint32_t as = 0; as < stitches.size(); ++as) {
			OnChainStitch cur_ocs(OnChainStitch::OnActive, ac, as);

			OnChainStitch prev_ocs;
			prev_ocs.on = OnChainStitch::OnActive;
			prev_ocs.chain = ac;
			if (as > 0 || is_loop) {
				prev_ocs.stitch = (as > 0 ? as - 1 : stitches.size()-1);
			} else {
				prev_ocs.type = OnChainStitch::TypeBegin;
			}
			OnChainStitch next_ocs;
			next_ocs.on = OnChainStitch::OnActive;
			next_ocs.chain = ac;
			if (as + 1 < stitches.size() || is_loop) {
				next_ocs.stitch = (as + 1 < stitches.size() ? as + 1 : 0);
			} else {
				next_ocs.type = OnChainStitch::TypeEnd;
			}

			auto fn = active_next.find(ChainStitch(ac, as));
			if (fn == active_next.end()) {
				//no link, so edge is previous to next:
				//      x    
				// p -> c -> n
				auto ret = next_vertex.insert(std::make_pair(
					std::make_pair(prev_ocs, cur_ocs),
					next_ocs
				));
				std::cout << "Inserted [" << ret.first->first.first << ", " << ret.first->first.second << "] -> " << ret.first->second << std::endl; //DEBUG
				assert(ret.second);
			} else {
				//have a link. check for discarded segments:

				//previous segment is discarded, so link prev up:
				if (!keep_adj[fn->second[0].chain][fn->second[0].stitch].first) {
					auto n = fn->second[0];
					auto fa = next_active.find(n);
					assert(fa != next_active.end());
					if (fa->second[0] == ChainStitch(ac, as)) {
						//   xxx n0
						//       | /
						// p --- c
						auto ret = next_vertex.insert(std::make_pair(
							std::make_pair(prev_ocs, cur_ocs),
								OnChainStitch(OnChainStitch::OnNext, n.chain, n.stitch)
							));
						std::cout << "Inserted [" << ret.first->first.first << ", " << ret.first->first.second << "] -> " << ret.first->second << std::endl; //DEBUG
						assert(ret.second);
					} else {
						assert(fa->second.size() == 2 && fa->second.back() == ChainStitch(ac, as));
						//don't link (this sort of case):
						//   xxx n0
						//     /  | 
						//   p -- c
					}
				}

				//next segment is discarded, so link down to next:
				if (!keep_adj[fn->second.back().chain][fn->second.back().stitch].second) {
					auto n = fn->second.back();
					auto fa = next_active.find(n);
					assert(fa != next_active.end());
					if (fa->second.back() == ChainStitch(ac, as)) {
						//   n0 xxx
						// \ |
						//   c --- n
						auto ret = next_vertex.insert(std::make_pair(
							std::make_pair(OnChainStitch(OnChainStitch::OnNext, n.chain, n.stitch), cur_ocs),
								next_ocs
							));
						std::cout << "Inserted [" << ret.first->first.first << ", " << ret.first->first.second << "] -> " << ret.first->second << std::endl; //DEBUG
						assert(ret.second);
					} else {
						assert(fa->second.size() == 2 && fa->second[0] == ChainStitch(ac, as));
						//don't link (this sort of case):
						//    n0 xxx
						//   / |
						// c - n
					}
				}
			}
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

	//DEBUG

	auto check_ocs = [&](OnChainStitch const &ocs) {
		assert(ocs.on == OnChainStitch::OnActive || ocs.on == OnChainStitch::OnNext);
		auto const &chains = (ocs.on == OnChainStitch::OnActive ? active_chains : next_chains);
		assert(ocs.chain < chains.size());
		auto const &stitches = (ocs.on == OnChainStitch::OnActive ? active_stitches : next_stitches);
		assert(stitches.size() == chains.size());
		auto const &stitche = stitches[ocs.chain];
		if (ocs.type == OnChainStitch::TypeStitch) {
			assert(ocs.stitch < stitche.size());
		} else {
			assert(ocs.type == OnChainStitch::TypeBegin || ocs.type == OnChainStitch::TypeEnd);
			assert(ocs.stitch == -1U);
		}
	};

	for (auto const &pv : next_vertex) {
		std::cout << "(" << pv.first.first << ", " << pv.first.second << ") -> " << pv.second << "\n";
		check_ocs(pv.first.first);
		check_ocs(pv.first.second);
		check_ocs(pv.second);
	}
	std::cout.flush();


	//Walk through created edges array, creating chains therefrom:

	std::vector< std::vector< OnChainStitch > > loops;
	std::map< std::pair< OnChainStitch, OnChainStitch >, std::vector< OnChainStitch > > partials;

	while (!next_vertex.empty()) {
		std::vector< OnChainStitch > chain;
		chain.emplace_back(next_vertex.begin()->first.first);
		chain.emplace_back(next_vertex.begin()->first.second);
		chain.emplace_back(next_vertex.begin()->second);
		//assert(chain[0] != chain[1] && chain[0] != chain[2] && chain[1] != chain[2]); //<-- not always true in two-stitch-loop cases (do we want two-stitch loops? Probably not.)
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



	//need lengths to figure out where stitches are on chains:
	//(duplicated, inelegantly, from link-chains)
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
		}
		return all_lengths;
	};
	std::vector< std::vector< float > > active_lengths = make_lengths(active_chains);
	std::vector< std::vector< float > > next_lengths = make_lengths(next_chains);

	auto output = [&](std::vector< OnChainStitch > const &path) {
		assert(path.size() >= 2);

		//build embedded vertices for all stitches:
		std::vector< ak::EmbeddedVertex > path_evs;
		std::vector< uint32_t > path_lefts; //also loop up i such that the stitch is in [ chain[i], chain[l+1] )
		path_evs.reserve(path.size());
		path_lefts.reserve(path.size());
		for (auto const &ocs : path) {
			std::vector< uint32_t > const &src_chain = (ocs.on == OnChainStitch::OnActive ? active_chains : next_chains)[ocs.chain];
			std::vector< float > const &src_lengths = (ocs.on == OnChainStitch::OnActive ? active_lengths : next_lengths)[ocs.chain];
			std::vector< ak::Stitch > const &src_stitches = (ocs.on == OnChainStitch::OnActive ? active_stitches : next_stitches)[ocs.chain];
			assert(!src_chain.empty());
			assert(src_lengths.size() == src_chain.size());
			std::cout << (path_evs.size()-1) << " " << ocs << " ";
			if (ocs.type == OnChainStitch::TypeBegin) {
				assert(src_chain[0] != src_chain.back());
				path_evs.emplace_back(ak::EmbeddedVertex::on_vertex(src_chain[0]));
				path_lefts.emplace_back(0);
				std::cout << "Begin: " << src_chain[0] << std::endl; //DEBUG
			} else if (ocs.type == OnChainStitch::TypeEnd) {
				assert(src_chain[0] != src_chain.back());
				path_evs.emplace_back(ak::EmbeddedVertex::on_vertex(src_chain.back()));
				path_lefts.emplace_back(src_chain.size()-2);
				std::cout << "End: " << src_chain.back() << std::endl; //DEBUG
			} else { assert(ocs.type == OnChainStitch::TypeStitch);
				assert(ocs.stitch < src_stitches.size());
				float l = src_lengths.back() * src_stitches[ocs.stitch].t;
				auto li = std::upper_bound(src_lengths.begin(), src_lengths.end(), l);
				assert(li != src_lengths.end());
				assert(li != src_lengths.begin());
				float m = (l - *(li-1)) / (*li - *(li-1));
				uint32_t i = li - src_lengths.begin();
				assert(i > 0);
				path_evs.emplace_back(ak::EmbeddedVertex::on_edge(src_chain[i-1], src_chain[i], m));
				path_lefts.emplace_back(i-1);
				std::cout << "Stitch: " << src_chain[i-1] << "-" << src_chain[i] << " at " << m << std::endl; //DEBUG
			}
		}

		std::vector< ak::EmbeddedVertex > chain;
		std::vector< ak::Stitch > stitches;
		float length = 0.0f;

		auto append_ev = [&chain,&slice,&length](ak::EmbeddedVertex const &ev) {
			if (!chain.empty()) {
				assert(ev != chain.back());
				ak::EmbeddedVertex::common_simplex(chain.back().simplex, ev.simplex); //make sure this works
				length += glm::length(
					chain.back().interpolate(slice.vertices) - ev.interpolate(slice.vertices)
				);
			}
			chain.emplace_back(ev);
		};

		for (uint32_t pi = 0; pi + 1 < path.size(); ++pi) {
			auto const &a = path[pi];
			auto const &a_ev = path_evs[pi];
			auto const &a_left = path_lefts[pi];
			std::vector< uint32_t > const &a_chain = (a.on == OnChainStitch::OnActive ? active_chains : next_chains).at(a.chain);
			assert(a_left + 1 < a_chain.size());
			auto const &b = path[pi+1];
			auto const &b_ev = path_evs[pi+1];
			auto const &b_left = path_lefts[pi+1];
			std::vector< uint32_t > const &b_chain = (b.on == OnChainStitch::OnActive ? active_chains : next_chains).at(b.chain);
			assert(b_left + 1 < b_chain.size());

			check_ocs(a); //DEBUG
			check_ocs(b); //DEBUG
			std::cout << "From " << a << " to " << b << std::endl; //DEBUG

			if (pi == 0) append_ev(a_ev);
			else assert(!chain.empty() && chain.back() == a_ev);
			if (a.type == OnChainStitch::TypeBegin) {
				assert(pi == 0);
				assert(a_left == 0);
			} else { assert(a.type == OnChainStitch::TypeStitch);
				std::vector< ak::Stitch > const &a_stitches = (a.on == OnChainStitch::OnActive ? active_stitches : next_stitches).at(a.chain);
				assert(a.stitch < a_stitches.size());
				ak::Stitch::Flag a_flag = a_stitches[a.stitch].flag;
				if (pi == 0) stitches.emplace_back(length, a_flag);
				else assert(!stitches.empty() && stitches.back().t == length && stitches.back().flag == a_flag);
			}

			if (a.on == b.on) {
				assert(a.chain == b.chain);
				bool is_loop = (a_chain[0] == a_chain.back());
				if (a_left != b_left) {
					assert(b_left + 1 < a_chain.size()); //RIIIIIGHT?
					uint32_t v = a_left;
					do {
						//advance v
						assert(v + 1 < a_chain.size());
						v += 1;
						if (v + 1 == a_chain.size()) {
							assert(is_loop);
							v = 0;
						}
						append_ev(ak::EmbeddedVertex::on_vertex(a_chain[v]));
					} while (v != b_left);
				}
			} else {
				//std::cout << "Building embedded path." << std::endl; //DEBUG
				//find an embedded path between a and b:
				std::vector< ak::EmbeddedVertex > ab;
				ak::embedded_path( parameters, slice, a_ev, b_ev, &ab);
				assert(ab[0] == a_ev);
				assert(ab.back() == b_ev);
				for (uint32_t i = 1; i + 1 < ab.size(); ++i) {
					append_ev(ab[i]);
				}
			}

			append_ev(b_ev);
			if (b.type == OnChainStitch::TypeEnd) {
				assert(pi + 2 == path.size());
				assert(b_left == b_chain.size() - 2);
			} else { assert(b.type == OnChainStitch::TypeStitch);
				std::vector< ak::Stitch > const &b_stitches = (b.on == OnChainStitch::OnActive ? active_stitches : next_stitches).at(b.chain);
				assert(b.stitch < b_stitches.size());
				ak::Stitch::Flag b_flag = b_stitches[b.stitch].flag;
				stitches.emplace_back(length, b_flag);
			}
		}

		//should turn loops into loops:
		assert((chain[0] == chain.back()) == (path[0] == path.back()));

		//remove last stitch from loops:
		if (path[0] == path.back()) {
			assert(!stitches.empty());
			assert(stitches[0].flag == stitches.back().flag);
			assert(stitches[0].t == 0.0f && stitches.back().t == length);
			stitches.pop_back();
		}

		//convert stitch t values from lengths to [0,1) range:
		for (auto &s : stitches) {
			s.t /= length;
		}

		//convert chain from being embedded on slice to being embedded on model:
		for (auto &ev : chain) {
			glm::uvec3 simplex = slice_on_model[ev.simplex.x].simplex;
			if (ev.simplex.y != -1U) simplex = ak::EmbeddedVertex::common_simplex(simplex, slice_on_model[ev.simplex.y].simplex);
			if (ev.simplex.z != -1U) simplex = ak::EmbeddedVertex::common_simplex(simplex, slice_on_model[ev.simplex.z].simplex);

			glm::vec3 weights = ev.weights.x * slice_on_model[ev.simplex.x].weights_on(simplex);
			if (ev.simplex.y != -1U) weights += ev.weights.y * slice_on_model[ev.simplex.y].weights_on(simplex);
			if (ev.simplex.z != -1U) weights += ev.weights.z * slice_on_model[ev.simplex.z].weights_on(simplex);

			ev = ak::EmbeddedVertex::canonicalize(simplex, weights);
		}

		next_active_chains.emplace_back(chain);
		next_active_stitches.emplace_back(stitches);
	};

	std::cout << "Found " << loops.size() << " loops and " << partials.size() << " chains." << std::endl;
	for (auto const &loop : loops) {
		output(loop);
	}
	for (auto const &pp : partials) {
		output(pp.second);
	}

	//HACK: sometimes duplicate vertices after splatting back to model somehow:
	uint32_t trimmed = 0;
	for (auto &chain : next_active_chains) {
		for (uint32_t i = 1; i < chain.size(); /* later */) {
			if (chain[i-1] == chain[i]) {
				chain.erase(chain.begin() + i);
				++trimmed;
			} else {
				++i;
			}
		}
	}

	if (trimmed) {
		std::cout << "Trimmed " << trimmed << " identical-after-moving-to-model vertices from next active chains." << std::endl;
	}


	//PARANOIA:
	assert(next_active_stitches.size() == next_active_chains.size());
	for (auto const &chain : next_active_chains) {
		for (uint32_t i = 1; i < chain.size(); ++i) {
			assert(chain[i-1] != chain[i]);
		}
	}
}
