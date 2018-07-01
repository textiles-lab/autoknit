#include "pipeline.hpp"

#include <glm/gtx/norm.hpp>
#include <glm/gtx/hash.hpp>

#include <iostream>
#include <unordered_map>
#include <unordered_set>


void ak::peel_chains(
	ak::Parameters const &parameters,
	ak::Model const &model, //in: model
	std::vector< float > const &times,          //in: time field (times @ vertices)
	std::vector< std::vector< ak::EmbeddedVertex > > const &active_chains, //in: current active chains
	std::vector< std::vector< ak::EmbeddedVertex > > *next_chains_ //out: next chains (may be different size than active_chains)
) {
	assert(next_chains_);
	auto &next_chains = *next_chains_;
	next_chains.clear();

	//This version of the code just uses the 3D distance to the curve.
	//might have problems with models that get really close to themselves.

	std::vector< float > values(model.vertices.size(), std::numeric_limits< float >::infinity());

	auto do_seg = [&values,&model](glm::vec3 const &a, glm::vec3 const &b) {
		if (a == b) return;
		glm::vec3 ab = b-a;
		float limit = glm::dot(ab,ab);
		float inv_limit = 1.0f / limit;
		for (auto const &v : model.vertices) {
			float amt = glm::dot(v-a, ab);
			amt = std::max(0.0f, std::min(limit, amt));
			glm::vec3 pt = (amt * inv_limit) * (b-a) + a;
			float dis2 = glm::length2(v - pt);
			float &best2 = values[&v - &model.vertices[0]];
			best2 = std::min(best2, dis2);
		}
	};

	for (auto const &chain : active_chains) {
		for (uint32_t i = 0; i + 1 < chain.size(); ++i) {
			do_seg(chain[i].interpolate(model.vertices), chain[i+1].interpolate(model.vertices));
		}
	}

	for (auto &v : values) {
		v = std::sqrt(v);
	}

	//TODO: blur/smooth distance fn somehow?

	//extract candidate chains:
	std::vector< std::vector< EmbeddedVertex > > possible_chains;

	//for (int step = 1; step < 10; ++step) { //DEBUG: show several distances
	{	int step = 1;

		float level = step * (parameters.stitch_height_mm / parameters.model_units_mm);

		std::vector< std::vector< EmbeddedVertex > > temp_chains;

		ak::extract_level_chains(model, values, level, &temp_chains);
		
		possible_chains.insert(possible_chains.end(), temp_chains.begin(), temp_chains.end());
	}

	//Trim chains that are not immediately left-of the active chains (should deal with problems of the surface looping back on itself?)

	//IDEA: start at active chain vertices and move to left-of until other chains are hit. Mark these chains as "okay".

	//these are (chain, index) pairs on sorted-index-order edges:
	struct EdgeVertex {
		EdgeVertex(float pos_, uint32_t chain_, uint32_t index_) : pos(pos_), chain(chain_), index(index_) { }
		float pos;
		uint32_t chain;
		uint32_t index;
	};
	std::unordered_multimap< glm::uvec2, EdgeVertex > edge_vertices;
	std::unordered_set< uint32_t > marked_vertices; //traversal will stop at marked vertices

	//add vertices for active chains to avoid crossing:
	for (auto const &ac : active_chains) {
		for (auto const &ev : ac) {
			if (ev.simplex.z != -1U) {
				//inside a triangle; ignore
			} else if (ev.simplex.y != -1U) {
				//edge vertex, easy!
				glm::uvec2 s = glm::uvec2(ev.simplex);
				assert(s.x < s.y);
				edge_vertices.insert(std::make_pair(s, EdgeVertex(ev.weights.y, -1U, -1U)));
			} else {
				//somewhat awkward; but at least stop traversal here:
				marked_vertices.insert(ev.simplex.x);
			}
		}
	}

	//add vertices for potential next chains because that's what we do:
	for (auto const &nc : possible_chains) {
		uint32_t nci = &nc - &possible_chains[0];
		for (auto const &ev : nc) {
			uint32_t evi = &ev - &nc[0];
			assert(ev.simplex.z == -1U);
			assert(ev.simplex.y != -1U);
			assert(ev.simplex.x < ev.simplex.y);
			edge_vertices.insert(std::make_pair(glm::uvec2(ev.simplex), EdgeVertex(ev.weights.y, nci, evi)));
		}
	}

	//adjacency list: seems useful
	std::vector< std::vector< uint32_t > > adj(model.vertices.size());
	{
		std::unordered_set< glm::uvec2 > sorted_edges;
		for (auto const &tri : model.triangles) {
			if (tri.x < tri.y) sorted_edges.insert(glm::uvec2(tri.x, tri.y));
			else               sorted_edges.insert(glm::uvec2(tri.y, tri.x));
			if (tri.y < tri.z) sorted_edges.insert(glm::uvec2(tri.y, tri.z));
			else               sorted_edges.insert(glm::uvec2(tri.z, tri.y));
			if (tri.z < tri.x) sorted_edges.insert(glm::uvec2(tri.z, tri.x));
			else               sorted_edges.insert(glm::uvec2(tri.x, tri.z));
		}
		for (auto const &e : sorted_edges) {
			adj[e.x].emplace_back(e.y);
			adj[e.y].emplace_back(e.x);
		}
	}

	//list of verts opposite directed edges: also seems useful
	std::unordered_map< glm::uvec2, uint32_t > opposite;
	{
		for (auto const &tri : model.triangles) {
			auto ret = opposite.insert(std::make_pair(glm::uvec2(tri.x, tri.y), tri.z));
			assert(ret.second);
			ret = opposite.insert(std::make_pair(glm::uvec2(tri.y, tri.z), tri.x));
			assert(ret.second);
			ret = opposite.insert(std::make_pair(glm::uvec2(tri.z, tri.x), tri.y));
			assert(ret.second);
		}
	}

	struct WalkAt {
		WalkAt(uint32_t from_, uint32_t to_, float pos_) : from(from_), to(to_), pos(pos_) {
		}
		uint32_t from, to;
		float pos;
	};
	std::vector< WalkAt > todo;

	//start left-of the active chain:
	for (auto const &ac : active_chains) {
		for (uint32_t ci = 0; ci + 1 < ac.size(); ++ci) {
			EmbeddedVertex const &a = ac[ci];
			EmbeddedVertex const &b = ac[ci+1];
			if (a.simplex.z != -1U && b.simplex.z != -1U) {
				//both inside a triangle, no edges to start from
				continue;
			}
			glm::uvec3 common = EmbeddedVertex::common_simplex(a.simplex, b.simplex);
			if (common.y == -1U) {
				//both on the same vertex -- weird!
				std::cerr << "WEIRD CASE: two chain verts on same vertex." << std::endl;
			} else if (common.z == -1U) {
				//both on the same edge; seed if one is at vertex:
				if (a.simplex.y == -1U) {
					//a is at vertex
					uint32_t to = (a.simplex.x == b.simplex.x ? b.simplex.y : b.simplex.x);
					auto f = opposite.find(glm::uvec2(a.simplex.x, to));
					if (f != opposite.end()) {
						todo.emplace_back(a.simplex.x, f->second, -1.0f);
					}

				}
				if (b.simplex.y == -1U) {
					//b is at vertex
					uint32_t from = (b.simplex.x == a.simplex.x ? a.simplex.y : a.simplex.x);
					auto f = opposite.find(glm::uvec2(from, b.simplex.x));
					if (f != opposite.end()) {
						todo.emplace_back(b.simplex.x, f->second, -1.0f);
					}
				}
			} else {
				assert(common.x < common.y && common.y < common.z && common.z != -1U);
				//both on the same triangle; seed if one is at edge (or vertex?):
				if (a.simplex.y == -1U) {
					//TODO: don't skip this case
					// rare(?) case: a is on vertex.
				} else if (a.simplex.z == -1U) {
					//a is on edge, figure out edge orientation:
					auto f1 = opposite.find(glm::uvec2(a.simplex.x, a.simplex.y));
					auto f2 = opposite.find(glm::uvec2(a.simplex.y, a.simplex.x));
					bool f1_matches = (f1 != opposite.end() && (f1->second == common.x || f1->second == common.y || f1->second == common.z));
					bool f2_matches = (f2 != opposite.end() && (f2->second == common.x || f2->second == common.y || f2->second == common.z));
					if (f1_matches) {
						assert(!f2_matches);
						todo.emplace_back(a.simplex.y, a.simplex.x, a.weights.x);
					} else { assert(f2_matches);
						todo.emplace_back(a.simplex.x, a.simplex.y, a.weights.y);
					}
				}

				if (b.simplex.y == -1U) {
					//TODO: don't skip this case
					// rare(?) case: b is on vertex.
				} else if (b.simplex.z == -1U) {
					//b is on edge, figure out edge orientation:
					auto f1 = opposite.find(glm::uvec2(b.simplex.x, b.simplex.y));
					auto f2 = opposite.find(glm::uvec2(b.simplex.y, b.simplex.x));
					bool f1_matches = (f1 != opposite.end() && (f1->second == common.x || f1->second == common.y || f1->second == common.z));
					bool f2_matches = (f2 != opposite.end() && (f2->second == common.x || f2->second == common.y || f2->second == common.z));
					if (f1_matches) {
						assert(!f2_matches);
						todo.emplace_back(b.simplex.x, b.simplex.y, b.weights.y);
					} else { assert(f2_matches);
						todo.emplace_back(b.simplex.y, b.simplex.x, b.weights.x);
					}
				}

			}
		}
	}

	std::unordered_set< uint32_t > found;

	//okay, now keep walking until marked vertices are encountered:
	while (!todo.empty()) {
		WalkAt at = todo.back();
		todo.pop_back();
		//is there a next planar map vertex along this edge?
		if (at.from < at.to) {
			auto r = edge_vertices.equal_range(glm::uvec2(at.from, at.to));
			auto next = r.second;
			float pos = at.pos;
			float next_at = std::numeric_limits< float >::infinity();
			for (auto ni = r.first; ni != r.second; ++ni) {
				if (pos < ni->second.pos && ni->second.pos < next_at) {
					next = ni;
					next_at = ni->second.pos;
				}
			}
			if (next != r.second) {
				//mark whatever we hit!
				if (next->second.chain != -1U) {
					found.insert(next->second.chain);
				}
				continue;
			}
		} else { assert(at.to < at.from);
			auto r = edge_vertices.equal_range(glm::uvec2(at.to, at.from));
			auto next = r.second;
			float pos = 1.0f - at.pos;
			float next_at = -std::numeric_limits< float >::infinity();
			for (auto ni = r.first; ni != r.second; ++ni) {
				if (next_at < ni->second.pos && ni->second.pos < pos) {
					next = ni;
					next_at = ni->second.pos;
				}
			}
			if (next != r.second) {
				//mark whatever we hit!
				if (next->second.chain != -1U) {
					found.insert(next->second.chain);
				}
				continue;
			}
		}
		//well, we got to 'to':
		if (marked_vertices.insert(at.to).second) {
			//if it was the first visit, got toward other adjacent things:
			for (auto a : adj[at.to]) {
				if (a == at.from) continue;
				todo.emplace_back(at.to, a, -1.0f);
			}
		}
	}


	std::cout << "Found " << found.size() << " of " << possible_chains.size() << " chains." << std::endl;

	next_chains.reserve(found.size());
	for (auto f : found) {
		next_chains.emplace_back();
		//next_chains.back() = possible_chains[f];
		sample_chain(parameters.get_chain_sample_spacing(), model, possible_chains[f], &next_chains.back());
	}

}

