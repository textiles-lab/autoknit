#include "ak-pipeline.hpp"

#include "EmbeddedPlanarMap.hpp"

#include <glm/gtx/hash.hpp>

#include <iostream>

inline std::ostream &operator<<(std::ostream &out, ak::EmbeddedVertex const &ev) {
	out << "(" << ev.weights.x << ", " << ev.weights.y << ", " << ev.weights.z << ")@[" << int32_t(ev.simplex.x) << ", " << int32_t(ev.simplex.y) << ", " << int32_t(ev.simplex.z) << "]";
	return out;
}

//EPM value that can track edge splits:
struct Edge {
	uint32_t a;
	uint32_t b;
	enum Type : uint8_t {
		Initial, //a,b are vertex indices
		Reverse, //a,b are the same edge index
		Combine, //a,b are edge indices
		SplitFirst, //a,b are same edge index
		SplitSecond, //a,b are edge edge index
	} type;
	Edge(Type type_, uint32_t a_, uint32_t b_) : a(a_), b(b_), type(type_) { }
};
std::vector< Edge > edges;

struct Value {
	int32_t sum;
	uint32_t edge;
	struct Reverse { static inline void reverse(Value *v) {
		v->sum = - v->sum;
		edges.emplace_back(Edge::Reverse, v->edge, v->edge);
		v->edge = edges.size() - 1;
	} };
	struct Combine { static inline void combine(Value *v, Value const &b) {
		v->sum += b.sum;
		edges.emplace_back(Edge::Combine, v->edge, b.edge);
		v->edge = edges.size() - 1;
	} };
	struct Split { static inline void split(Value const &v, Value *first, Value *second) {
		first->sum = v.sum;
		second->sum = v.sum;

		edges.emplace_back(Edge::SplitFirst, v.edge, v.edge);
		first->edge = edges.size() - 1;
		edges.emplace_back(Edge::SplitSecond, v.edge, v.edge);
		second->edge = edges.size() - 1;
	} };
};

void ak::trim_model(
	ak::Model const &model, //in: model
	std::vector< std::vector< ak::EmbeddedVertex > > const &left_of,
	std::vector< std::vector< ak::EmbeddedVertex > > const &right_of,
	ak::Model *clipped_, //out: portion of model's surface that is left_of the left_of chains and right_of the right_of chains
	std::vector< ak::EmbeddedVertex > *clipped_vertices_, //out: map from clipped vertices to source mesh
	std::vector< std::vector< uint32_t > > *left_of_vertices_, //out (optional): indices of vertices corresponding to left_of chains [may be some rounding]
	std::vector< std::vector< uint32_t > > *right_of_vertices_ //out (optional): indices of vertices corresponding to right_of chains [may be some rounding]
) {
	assert(clipped_);
	auto &clipped = *clipped_;
	clipped.clear();

	assert(clipped_vertices_);
	auto &clipped_vertices = *clipped_vertices_;
	clipped_vertices.clear();

	{ //PARANOIA: make sure all chains are loops or edge-to-edge:
		std::unordered_set< glm::uvec2 > edges;
		for (auto const &tri : model.triangles) {
			auto do_edge = [&edges](uint32_t a, uint32_t b) {
				if (a > b) std::swap(a,b);
				auto ret = edges.insert(glm::uvec2(a,b));
				if (!ret.second) edges.erase(ret.first);
			};
			do_edge(tri.x, tri.y);
			do_edge(tri.y, tri.z);
			do_edge(tri.z, tri.x);
		}
		std::unordered_set< uint32_t > edge_verts;
		for (auto const &e : edges) {
			edge_verts.insert(e.x);
			edge_verts.insert(e.y);
		}

		auto on_edge = [&edges, &edge_verts](ak::EmbeddedVertex const &ev) -> bool {
			if (ev.simplex.z != -1U) {
				return false;
			} else if (ev.simplex.y != -1U) {
				return edges.count(glm::uvec2(ev.simplex.x, ev.simplex.y)) != 0;
			} else {
				return edge_verts.count(ev.simplex.x) != 0;
			}
		};

		for (auto const &chain : left_of) {
			assert(chain.size() >= 2);
			if (chain[0] == chain.back()) continue;
			assert(on_edge(chain[0]));
			assert(on_edge(chain.back()));
		}

		for (auto const &chain : right_of) {
			assert(chain.size() >= 2);
			if (chain[0] == chain.back()) continue;
			assert(on_edge(chain[0]));
			assert(on_edge(chain.back()));
		}
	}

	//embed chains using planar map:

	edges.clear(); //global list used for edge tracking in epm; awkward but should work.
	EmbeddedPlanarMap< Value, Value::Reverse, Value::Combine, Value::Split > epm;
	std::unordered_set< glm::uvec2 > empty_edges; //when it's the same vertex after rounding
	uint32_t total_chain_edges = 0;
	uint32_t fresh_id = 0;
	for (auto const &chain : left_of) {
		uint32_t prev = epm.add_vertex(chain[0]);
		uint32_t prev_id = fresh_id++;
		for (uint32_t i = 1; i < chain.size(); ++i) {
			uint32_t cur = epm.add_vertex(chain[i]);
			uint32_t cur_id = fresh_id++;
			if (prev == cur) {
				std::cout << "NOTE: vertex " << chain[i-1] << " and " << chain[i] << " (in a left_of chain) round to the same value." << std::endl; //DEBUG
				empty_edges.insert(glm::uvec2(prev_id, cur_id));
			}
			Value value;
			value.sum = 1;
			edges.emplace_back(Edge::Initial, prev_id, cur_id);
			value.edge = edges.size() - 1;
			epm.add_edge(prev, cur, value);
			prev = cur;
			prev_id = cur_id;
			++total_chain_edges;
		}
	}

	for (auto const &chain : right_of) {
		uint32_t prev = epm.add_vertex(chain[0]);
		uint32_t prev_id = fresh_id++;
		for (uint32_t i = 1; i < chain.size(); ++i) {
			uint32_t cur = epm.add_vertex(chain[i]);
			uint32_t cur_id = fresh_id++;
			if (prev == cur) {
				std::cout << "NOTE: vertex " << chain[i-1] << " and " << chain[i] << " (in a right_of chain) round to the same value." << std::endl; //DEBUG
				empty_edges.insert(glm::uvec2(prev_id, cur_id));
			}
			Value value;
			value.sum = (1 << 8);
			edges.emplace_back(Edge::Initial, cur_id, prev_id);
			value.edge = edges.size() - 1;
			epm.add_edge(cur, prev, value);
			prev = cur;
			prev_id = cur_id;
			++total_chain_edges;
		}
	}

	uint32_t total_simplex_edges = 0;
	for (const auto &edges : epm.simplex_edges) {
		total_simplex_edges += edges.second.size();
	}
	std::cout << "EPM has " << epm.vertices.size() << " vertices." << std::endl;
	std::cout << "EPM has " << epm.simplex_vertices.size() << " simplices with vertices." << std::endl;
	std::cout << "EPM has " << epm.simplex_edges.size() << " simplices with edges (" << total_simplex_edges << " edges from " << total_chain_edges << " chain edges)." << std::endl;


	//clean up any small loops that may exist in chains:

	//read out left_of_vertices and right_of_vertices from edge information:
	std::vector< std::vector< uint32_t > > left_of_epm, right_of_epm;
	{
		//for each original id->id edge, extract the chain of vertices that it expands to:
		struct Source {
			Source(uint32_t a_, uint32_t b_, uint32_t num_, uint32_t den_) : a(a_), b(b_), num(num_), den(den_) { }
			uint32_t a, b; //a / b vertices
			uint32_t num, den; //position along (subdivided?) edge
		};
		std::vector< std::vector< Source > > sources;
		sources.reserve(edges.size());
		for (uint32_t e = 0; e < edges.size(); ++e) {
			assert(e == sources.size());
			sources.emplace_back();
			if (edges[e].type == Edge::Initial) {
				assert(edges[e].a != edges[e].b);
				sources.back().emplace_back(edges[e].a, edges[e].b, 1, 2);
			} else if (edges[e].type == Edge::Reverse) {
				assert(edges[e].a == edges[e].b);
				assert(edges[e].a < e);
				for (auto const &s : sources[edges[e].a]) {
					sources.back().emplace_back(s.b, s.a, s.den - s.num, s.den);
				}
			} else if (edges[e].type == Edge::Combine) {
				assert(edges[e].a < e);
				assert(edges[e].b < e);
				for (auto const &s : sources[edges[e].a]) {
					sources.back().emplace_back(s.a, s.b, s.num, s.den);
				}
				for (auto const &s : sources[edges[e].b]) {
					sources.back().emplace_back(s.a, s.b, s.num, s.den);
				}
			} else if (edges[e].type == Edge::SplitFirst) {
				assert(edges[e].a == edges[e].b);
				assert(edges[e].a < e);
				for (auto const &s : sources[edges[e].a]) {
					sources.back().emplace_back(s.a, s.b, s.num * 2 - 1, s.den * 2);
				}
			} else if (edges[e].type == Edge::SplitSecond) {
				assert(edges[e].a == edges[e].b);
				assert(edges[e].a < e);
				for (auto const &s : sources[edges[e].a]) {
					assert(uint32_t(s.den * 2) > s.den); //make sure we're not overflowing
					sources.back().emplace_back(s.a, s.b, s.num * 2 + 1, s.den * 2);
				}
			} else {
				assert(false);
			}
		}

		struct SubEdge {
			SubEdge(uint32_t a_, uint32_t b_, uint32_t num_, uint32_t den_) : a(a_), b(b_), num(num_), den(den_) { }
			uint32_t a, b; //a / b vertices (epm indices)
			uint32_t num, den; //position of center
		};

		//now iterate actual epm edges and check where they end up mapping.
		std::unordered_map< glm::uvec2, std::vector< SubEdge > > edge_subedges; //<-- indexed by 'vertex_id' values, contains epm.vertices indices
		for (auto const &se : epm.simplex_edges) {
			for (auto const &ee : se.second) {
				assert(ee.value.edge < sources.size());
				assert(!sources[ee.value.edge].empty());
				//copy subedge(s) corresponding to this edge to their source edges:
				for (auto const &s : sources[ee.value.edge]) {
					assert(s.a != s.b);
					if (s.a < s.b) {
						//if (s.den > 2) std::cout << s.a << "/" << s.b << " -> [" << ee.first << "-" << ee.second << "] (non-flipped)" << std::endl; //DEBUG
						edge_subedges[glm::uvec2(s.a, s.b)].emplace_back(ee.first, ee.second, s.num, s.den);
					} else {
						//if (s.den > 2) std::cout << s.b << "/" << s.a << " -> [" << ee.second << "-" << ee.first << "] (flipped)" << std::endl; //DEBUG
						edge_subedges[glm::uvec2(s.b, s.a)].emplace_back(ee.second, ee.first, s.den - s.num, s.den);
					}
				}
			}
		}
		//make sure subedge chains are logically sorted:
		for (auto &ese : edge_subedges) {
			std::vector< SubEdge > &subedges = ese.second;
			assert(!subedges.empty());
			std::stable_sort(subedges.begin(), subedges.end(), [](SubEdge const &a, SubEdge const &b){
				//    a.num / a.den < b.num / b.den
				// => a.num * b.den < b.num * a.den
				return uint64_t(a.num) * uint64_t(b.den) < uint64_t(b.num) * uint64_t(a.den);
			});
			/*if (subedges.size() > 1) {
				//DEBUG:
				for (auto &s : subedges) {
					std::cout << " [" << s.a << "-" << s.b << "]@(" << s.num << "/" << s.den << ")";
				}
				std::cout << std::endl;
			}*/
			for (uint32_t i = 0; i + 1 < subedges.size(); ++i) {
				assert(subedges[i].b == subedges[i+1].a);
			}
		}

		//okay, now read back vertex chains by querying by ID values:

		uint32_t fresh_id = 0;
		left_of_epm.reserve(left_of.size());
		for (auto const &chain : left_of) {
			left_of_epm.emplace_back();
			std::vector< uint32_t > &epm_chain = left_of_epm.back();
			uint32_t prev_id = fresh_id++;
			for (uint32_t i = 1; i < chain.size(); ++i) {
				uint32_t cur_id = fresh_id++;
				auto f = edge_subedges.find(glm::uvec2(prev_id, cur_id));
				if (empty_edges.count(glm::uvec2(prev_id, cur_id))) {
					assert(f == edge_subedges.end());
				} else {
					assert(f != edge_subedges.end());
					for (auto const &se : f->second) {
						if (epm_chain.empty()) epm_chain.emplace_back(se.a);
						else assert(epm_chain.back() == se.a);
						epm_chain.emplace_back(se.b);
					}
				}
				prev_id = cur_id;
			}
			assert((chain[0] == chain.back()) == (epm_chain[0] == epm_chain.back()));
		}
		right_of_epm.reserve(right_of.size());
		for (auto const &chain : right_of) {
			right_of_epm.emplace_back();
			std::vector< uint32_t > &epm_chain = right_of_epm.back();
			uint32_t prev_id = fresh_id++;
			for (uint32_t i = 1; i < chain.size(); ++i) {
				uint32_t cur_id = fresh_id++;
				auto f = edge_subedges.find(glm::uvec2(prev_id, cur_id));
				if (empty_edges.count(glm::uvec2(prev_id, cur_id))) {
					assert(f == edge_subedges.end());
				} else {
					assert(f != edge_subedges.end());
					for (auto const &se : f->second) {
						if (epm_chain.empty()) epm_chain.emplace_back(se.a);
						else assert(epm_chain.back() == se.a);
						epm_chain.emplace_back(se.b);
					}
				}
				prev_id = cur_id;
			}
			assert((chain[0] == chain.back()) == (epm_chain[0] == epm_chain.back()));
		}
	}

	auto cleanup_chain = [&](std::vector< uint32_t > &epm_chain, int32_t value) {
		//want chain to be loop-free --> look for loops!
		float length_removed = 0.0f;
		uint32_t verts_removed = 0;

		//useful:
		auto remove_from_epm = [&epm_chain, &epm, value](uint32_t first, uint32_t last) {
			assert(first <= last);
			assert(last < epm_chain.size());
			//remove value of all segments [first,last] from planar map:
			for (uint32_t i = first; i + 1 <= last; ++i) {
				auto &a = epm.vertices[epm_chain[i]];
				auto &b = epm.vertices[epm_chain[i+1]];
				glm::uvec3 common = IntegerEmbeddedVertex::common_simplex(a.simplex, b.simplex);
				auto f = epm.simplex_edges.find(common);
				assert(f != epm.simplex_edges.end());
				bool found = false;
				for (auto &e : f->second) {
					if (e.first == epm_chain[i] && e.second == epm_chain[i+1]) {
						e.value.sum -= value;
						found = true;
						break;
					} else if (e.second == epm_chain[i] && e.first == epm_chain[i+1]) {
						e.value.sum += value;
						found = true;
						break;
					}
				}
				assert(found);
			}
		};

		float initial_length = std::numeric_limits< float >::quiet_NaN();

		bool again = true;
		while (again) {
			again = false;

			std::vector< float > lengths;
			lengths.reserve(epm_chain.size());
			lengths.emplace_back(0.0f);
			for (uint32_t i = 1; i < epm_chain.size(); ++i) {
				glm::vec3 a = epm.vertices[epm_chain[i-1]].interpolate(model.vertices);
				glm::vec3 b = epm.vertices[epm_chain[i]].interpolate(model.vertices);
				lengths.emplace_back(lengths.back() + glm::length(b-a));
			}
			if (!(initial_length == initial_length)) initial_length = lengths.back();

			std::unordered_map< uint32_t, uint32_t > visited;
			visited.reserve(epm.vertices.size());

			uint32_t first_i = (epm_chain[0] == epm_chain.back() ? 1 : 0);
			for (uint32_t i = first_i; i < epm_chain.size(); ++i) {
				auto ret = visited.insert(std::make_pair(epm_chain[i], i));
				if (ret.second) continue;

				uint32_t first = ret.first->second;
				uint32_t last = i;
				assert(first < last);
				assert(epm_chain[first] == epm_chain[last]);

				float length = lengths[last] - lengths[first];
				float length_outer = (lengths[first] - lengths[0]) + (lengths.back() - lengths[last]);

				
				if (epm_chain[0] == epm_chain.back() && length_outer < length) {
					//remove the 'outer' loop -- (last,back] + [0,first)
					verts_removed += first + (epm_chain.size() - (last + 1));
					length_removed += length_outer;

					remove_from_epm(0, first);
					remove_from_epm(last, epm_chain.size()-1);
					epm_chain.erase(epm_chain.begin() + last + 1, epm_chain.end());
					epm_chain.erase(epm_chain.begin(), epm_chain.begin() + first);
					assert(epm_chain[0] == epm_chain.back()); //preserve circularity, right?
				} else {
					//remove the inner loop -- (first, last)
					verts_removed += (last + 1) - (first + 1);
					length_removed += length;

					remove_from_epm(first, last);
					epm_chain.erase(epm_chain.begin() + first + 1, epm_chain.begin() + last + 1);
				}
				again = true;
				break;
			}
		} //while (again)
		
		if (verts_removed) {
			std::cout << "Removed " << verts_removed << " vertices (that's " << length_removed << " units; " << length_removed / initial_length * 100.0 << "% of the initial length of " << initial_length << " units)." << std::endl;
		}

	};

	for (auto &epm_chain : left_of_epm) {
		cleanup_chain(epm_chain, 1);
	}

	for (auto &epm_chain : right_of_epm) {
		cleanup_chain(epm_chain, -(1 << 8));
	}

	//----- end loop cleanup -----

	
	//build split mesh:
	std::vector< ak::EmbeddedVertex > split_verts;
	std::vector< glm::uvec3 > split_tris;
	std::vector< uint32_t > epm_to_split;
	epm.split_triangles(model.vertices, model.triangles, &split_verts, &split_tris, &epm_to_split);


	//transfer edge values to split mesh:
	std::unordered_map< glm::uvec2, int32_t > edge_values;
	for (auto const &se : epm.simplex_edges) {
		for (auto const &ee : se.second) {
			uint32_t a = epm_to_split[ee.first];
			uint32_t b = epm_to_split[ee.second];
			int32_t value = ee.value.sum;
			auto ret = edge_values.insert(std::make_pair(glm::uvec2(a,b), value));
			assert(ret.second);

			ret = edge_values.insert(std::make_pair(glm::uvec2(b,a), -value));
			assert(ret.second);
		}
	}
	
	//tag triangles with values:
	std::unordered_map< glm::uvec2, uint32_t > edge_to_tri;
	edge_to_tri.reserve(split_tris.size() * 3);
	for (auto const &tri : split_tris) {
		uint32_t ti = &tri - &split_tris[0];
		auto ret = edge_to_tri.insert(std::make_pair(glm::uvec2(tri.x, tri.y), ti)); assert(ret.second);
		ret = edge_to_tri.insert(std::make_pair(glm::uvec2(tri.y, tri.z), ti)); assert(ret.second);
		ret = edge_to_tri.insert(std::make_pair(glm::uvec2(tri.z, tri.x), ti)); assert(ret.second);
	}

	constexpr int32_t const Unvisited = std::numeric_limits< int32_t >::max();
	std::vector< int32_t > values(split_tris.size(), Unvisited);

	std::vector< bool > keep(split_tris.size(), false);

	for (uint32_t seed = 0; seed < split_tris.size(); ++seed) {
		if (values[seed] != Unvisited) continue;

		std::vector< uint32_t > component;
		component.reserve(split_tris.size());

		values[seed] = (128 << 8) | (128);
		component.emplace_back(seed);

		for (uint32_t ci = 0; ci < component.size(); ++ci) {
			uint32_t ti = component[ci];
			uint32_t value = values[ti];
			assert(value != Unvisited);
			glm::uvec3 tri = split_tris[ti];
			auto over = [&](uint32_t a, uint32_t b) {
				auto f = edge_to_tri.find(glm::uvec2(b,a));
				if (f == edge_to_tri.end()) return;
				int32_t nv = value;
				auto f2 = edge_values.find(glm::uvec2(b,a));
				if (f2 != edge_values.end()) nv += f2->second;

				if (values[f->second] == Unvisited) {
					values[f->second] = nv;
					component.emplace_back(f->second);
				} else {
					assert(values[f->second] == nv);
				}
			};
			over(tri.x, tri.y);
			over(tri.y, tri.z);
			over(tri.z, tri.x);
		}

		int32_t max_right = 128;
		int32_t max_left = 128;
		for (auto ti : component) {
			int32_t right = values[ti] >> 8;
			int32_t left = values[ti] & 0xff;
			max_right = std::max(max_right, right);
			max_left = std::max(max_left, left);
		}

		int32_t keep_value = (max_right << 8) | max_left;
		for (auto ti : component) {
			if (values[ti] == keep_value) {
				keep[ti] = true;
			}
		}
	}

	std::vector< uint32_t > split_vert_to_clipped_vertex(split_verts.size(), -1U);
	auto use_vertex = [&](uint32_t v) {
		if (split_vert_to_clipped_vertex[v] == -1U) {
			assert(v < split_verts.size());
			split_vert_to_clipped_vertex[v] = clipped_vertices.size();
			clipped_vertices.emplace_back(split_verts[v]);
		}
		return split_vert_to_clipped_vertex[v];
	};

	for (auto const &tri : split_tris) {
		if (!keep[&tri - &split_tris[0]]) continue;
		clipped.triangles.emplace_back(
			use_vertex(tri.x),
			use_vertex(tri.y),
			use_vertex(tri.z)
		);
	}
	clipped.vertices.reserve(clipped_vertices.size());
	for (auto const &v : clipped_vertices) {
		clipped.vertices.emplace_back(v.interpolate(model.vertices));
	}

	std::cout << "Trimmed model from " << model.triangles.size() << " triangles on " << model.vertices.size() << " vertices to " << clipped.triangles.size() << " triangles on " << clipped.vertices.size() << " vertices." << std::endl;

	//transform vertex indices for left_of and right_of vertices -> clipped model:
	auto transform_chain = [&](std::vector< uint32_t > const &epm_chain) {
		assert(!epm_chain.empty());
		std::vector< uint32_t > split_chain;
		split_chain.reserve(epm_chain.size());
		for (auto v : epm_chain) {
			assert(v < epm_to_split.size());
			v = epm_to_split[v];
			assert(v != -1U);
			assert(v < split_vert_to_clipped_vertex.size());
			v = split_vert_to_clipped_vertex[v];
			//assert(v != -1U); //<-- sometimes chains don't include the edge loops (like when they cross)
			split_chain.emplace_back(v);
		}
		assert(!split_chain.empty());
		assert((epm_chain[0] == epm_chain.back()) == (split_chain[0] == split_chain.back()));
		return split_chain;
	};
	if (left_of_vertices_) {
		left_of_vertices_->clear();
		left_of_vertices_->reserve(left_of_epm.size());
		for (auto const &chain : left_of_epm) {
			//std::cout << "left_of[" << (&chain - &left_of_epm[0]) << "]:"; //DEBUG
			left_of_vertices_->emplace_back(transform_chain(chain));
		}
	}
	if (right_of_vertices_) {
		right_of_vertices_->clear();
		right_of_vertices_->reserve(right_of_epm.size());
		for (auto const &chain : right_of_epm) {
			//std::cout << "right_of[" << (&chain - &right_of_epm[0]) << "]:"; //DEBUG
			right_of_vertices_->emplace_back(transform_chain(chain));
		}
	}


}

