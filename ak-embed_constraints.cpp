#include "pipeline.hpp"

#include <glm/gtx/norm.hpp>
#include <glm/gtx/hash.hpp>

#include <algorithm>
#include <deque>
#include <iostream>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>


#define WEIGHT_SUM  (1024U)
struct IntegerEmbeddedVertex {
	glm::uvec3 simplex;
	glm::ivec3 weights;

	IntegerEmbeddedVertex(const glm::uvec3 &simplex_, const glm::ivec3 &weights_) : simplex(simplex_), weights(weights_) {
		assert(weights.x + weights.y + weights.z == WEIGHT_SUM);
	}

	static IntegerEmbeddedVertex simplify(const IntegerEmbeddedVertex &in) {
		assert(in.simplex.x < in.simplex.y && in.simplex.y <= in.simplex.z);
		assert(in.weights.x + in.weights.y + in.weights.z == WEIGHT_SUM);
		IntegerEmbeddedVertex ret(glm::uvec3(0, -1U, -1U), glm::ivec3(1024, 0,0));
		uint32_t o = 0;
		for (uint32_t i = 0; i < 3; ++i) {
			if (in.weights[i] != 0) {
				ret.simplex[o] = in.simplex[i];
				ret.weights[o] = in.weights[i];
				++o;
			}
		}
		return ret;
	}

	IntegerEmbeddedVertex(const ak::EmbeddedVertex &src) : simplex(src.simplex) {
		assert(simplex.x != -1U);
		assert(simplex.x < simplex.y);
		assert((simplex.y == -1U && simplex.z == -1U) || simplex.y < simplex.z);

		float x = src.weights.x;
		float xy = x + src.weights.y;
		float xyz = xy + src.weights.z;

		int32_t ix = std::round(WEIGHT_SUM * (x / xyz));
		int32_t ixy = std::round(WEIGHT_SUM * (xy / xyz));
		int32_t ixyz = std::round(WEIGHT_SUM * (xyz / xyz));

		weights = glm::ivec3(ix, ixy - ix, ixyz - ixy);
		assert(weights.x + weights.y + weights.z == WEIGHT_SUM);
	};

	bool operator==(const IntegerEmbeddedVertex &o) const {
		return simplex == o.simplex && weights == o.weights;
	}

	glm::ivec3 weights_on(glm::uvec3 simplex2) const {
		glm::ivec3 ret(0,0,0);
		uint32_t o = 0;
		for (uint32_t i = 0; i < 3; ++i) {
			if (simplex[i] == -1U) break;
			while (simplex2[o] < simplex[i]) {
				++o;
				assert(o < 3);
			}
			assert(simplex2[o] == simplex[i]);
			ret[o] = weights[i];
		}
		assert(ret.x + ret.y + ret.z == WEIGHT_SUM);

		return ret;
	}

	static glm::uvec3 common_simplex(const glm::uvec3 &a, const glm::uvec3 &b) {
		return ak::EmbeddedVertex::common_simplex(a,b);
	}
};

struct EmbeddedEdge {
	EmbeddedEdge(uint32_t first_, uint32_t second_, float value_) : first(first_), second(second_), value(value_) { }
	uint32_t first;
	uint32_t second;
	float value;
};

struct EmbeddedPlanarMap {
	std::vector< IntegerEmbeddedVertex > vertices;

	std::unordered_map< glm::uvec3, std::vector< uint32_t > > simplex_vertices;
	std::unordered_map< glm::uvec3, std::vector< EmbeddedEdge > > simplex_edges;

	enum LineSide : int8_t {
		Left = 1,
		On = 0,
		Right = -1,
	};

	//is there a point in the interior of both a-b and a2-b2?
	bool segments_intersect(const glm::ivec2 &a, const glm::ivec2 &b, const glm::ivec2 &a2, const glm::ivec2 &b2) {
		if (a == b) return false;
		if (a2 == b2) return false;

		int32_t perp_a2 = (a2 - a).x * -(b - a).y + (a2 - a).y * (b - a).x;
		int32_t perp_b2 = (b2 - a).x * -(b - a).y + (b2 - a).y * (b - a).x;

		if (perp_a2 == 0 && perp_b2 == 0) {
			//annoying (colinear) special case
			int32_t along_a2 = (a2 - a).x * (b - a).x + (a2 - a).y * (b - a).y;
			int32_t along_b2 = (b2 - a).x * (b - a).x + (b2 - a).y * (b - a).y;
			int32_t limit = (b - a).x * (b - a).x + (b - a).y * (b - a).y;

			if (along_a2 <= 0 && along_b2 <= 0) return false;
			if (along_a2 >= limit && along_b2 >= limit) return false;
			return true;
		}

		if (perp_a2 <= 0 && perp_b2 <= 0) return false;
		if (perp_a2 >= 0 && perp_b2 >= 0) return false;

		int32_t perp_a = (a - a2).x * -(b2 - a2).y + (a - a2).y * (b2 - a2).x;
		int32_t perp_b = (b - a2).x * -(b2 - a2).y + (b - a2).y * (b2 - a2).x;

		assert(!(perp_a == 0 && perp_b == 0)); //should have been handled above

		if (perp_a <= 0 && perp_b <= 0) return false;
		if (perp_a >= 0 && perp_b >= 0) return false;

		return true;
	}

	glm::ivec2 rounded_intersection(const glm::ivec2 &a, const glm::ivec2 &b, const glm::ivec2 &a2, const glm::ivec2 &b2) {
		//NOTE: should call only when actually intersecting.
		//NOTE2: should ~probably~ use the 2x2 matrix inverse formulation here
		int32_t perp_a2 = (a2 - a).x * -(b - a).y + (a2 - a).y * (b - a).x;
		int32_t perp_b2 = (b2 - a).x * -(b - a).y + (b2 - a).y * (b - a).x;
		assert(perp_a2 != perp_b2);

		double t = double(0 - perp_a2) / double(perp_b2 - perp_a2);
		return glm::ivec2(
			std::round((b.x - a.x) * t + double(a.x)),
			std::round((b.y - a.y) * t + double(a.y))
		);
	}


	//is point 'pt' in the interior of line segment a-b?
	bool point_in_segment(const glm::ivec2 &pt, const glm::ivec2 &a, const glm::ivec2 &b) {
		glm::ivec2 ab = b - a;
		glm::ivec2 ap = pt - a;

		int32_t perp = ap.x * -ab.y + ap.y * ab.x;
		if (perp != 0) return false;
		
		int32_t along = ap.x * ab.x + ap.y * ab.y;
		if (along <= 0) return false;

		int32_t limit = ab.x * ab.x + ab.y * ab.y;
		if (along >= limit) return false;

		return true;
	}

	bool point_in_segment(const IntegerEmbeddedVertex &pt_, const IntegerEmbeddedVertex &a_, const IntegerEmbeddedVertex &b_) {
		//work in barycentric coordinates:
		glm::uvec3 common = IntegerEmbeddedVertex::common_simplex(pt_.simplex, IntegerEmbeddedVertex::common_simplex(a_.simplex, b_.simplex));

		glm::ivec2 pt = glm::ivec2(pt_.weights_on(common));
		glm::ivec2 a = glm::ivec2(a_.weights_on(common));
		glm::ivec2 b = glm::ivec2(b_.weights_on(common));

		return point_in_segment(pt, a, b);
	}

	
	uint32_t add_vertex(const IntegerEmbeddedVertex &v_) {
		IntegerEmbeddedVertex v = IntegerEmbeddedVertex::simplify(v_);

		auto &verts = simplex_vertices[v.simplex];
		auto &edges = simplex_edges[v.simplex];

		for (auto i : verts) {
			if (vertices[i] == v) return i;
		}
		/* //PARANOIA:
		assert(v.simplex.x < v.simplex.y && v.simplex.y <= v.simplex.z);
		assert(v.weights.x + v.weights.y + v.weights.z == WEIGHT_SUM);
		for (const auto &v2 : vertices) {
			assert(!(v2 == v));
		}
		//end PARANOIA*/

		uint32_t idx = vertices.size();
		vertices.emplace_back(v);
		verts.emplace_back(idx);

		uint32_t old_size = edges.size();

		for (uint32_t e = 0; e < old_size; ++e) {
			if (point_in_segment(v, vertices[edges[e].first], vertices[edges[e].second])) {
				auto second_half = edges[e];
				second_half.first = idx;
				edges[e].second = idx;
				edges.emplace_back(second_half);
			}
		}
		return idx;
	}

	void add_edge(uint32_t ai, uint32_t bi, float value) {
		assert(ai < vertices.size() && bi < vertices.size());

		if (ai == bi) return;

		const auto &a_ = vertices[ai];
		const auto &b_ = vertices[bi];
		glm::uvec3 common = IntegerEmbeddedVertex::common_simplex(a_.simplex, b_.simplex);

		glm::ivec2 a = glm::ivec2(a_.weights_on(common));
		glm::ivec2 b = glm::ivec2(b_.weights_on(common));

		//split edge at any vertices in simplex:
		auto &verts = simplex_vertices[common];
		for (auto vi : verts) {
			assert(vi < vertices.size());
			glm::ivec2 v = glm::vec2(vertices[vi].weights_on(common));
			if (point_in_segment(v, a, b)) {
				add_edge(ai, vi, value);
				add_edge(vi, bi, value);
				return;
			}
		}

		//split edge (and add new vertex) if there is an intersection:
		auto &edges = simplex_edges[common];
		for (uint32_t e = 0; e < edges.size(); ++e) {

			//if it matches the edge, over-write value & done!
			if ((edges[e].first == ai && edges[e].second == bi) || (edges[e].first == bi && edges[e].second == ai)) {
				edges[e].value = value;
				return;
			}

			glm::ivec2 a2 = glm::ivec2(vertices[edges[e].first].weights_on(common));
			glm::ivec2 b2 = glm::ivec2(vertices[edges[e].second].weights_on(common));

			//if endpoints are interior to an existing edge, split existing edge:
			if (point_in_segment(a, a2, b2)) {
				auto second_half = edges[e];
				second_half.first = ai;
				edges.emplace_back(second_half);
				edges[e].second = ai;
				b2 = a;
			}
			if (point_in_segment(b, a2, b2)) {
				auto second_half = edges[e];
				second_half.first = bi;
				edges.emplace_back(second_half);
				edges[e].second = bi;
				b2 = b;
			}

			//if edges cross, remove, add intersection, and re-insert:
			if (segments_intersect(a,b, a2,b2)) {
				glm::ivec3 pt = glm::ivec3(rounded_intersection(a,b,a2,b2), 0);
				pt.z = WEIGHT_SUM - pt.x - pt.y;
				uint32_t pti = add_vertex(IntegerEmbeddedVertex(common, pt));

				float ai2 = edges[e].first;
				float bi2 = edges[e].second;
				float value2 = edges[e].value;
				edges.erase(edges.begin() + e);

				add_edge(ai2, pti, value2);
				add_edge(pti, bi2, value2);
				add_edge(ai, pti, value);
				add_edge(pti, bi, value);

				return;
			}
			
		}

		//if got to this point, no intersections:
		edges.emplace_back(ai, bi, value);


	}
	void add_edge(const ak::EmbeddedVertex &a, const ak::EmbeddedVertex &b, float value) {
		uint32_t ai = add_vertex(a);
		uint32_t bi = add_vertex(b);
		add_edge(ai, bi, value);
	}
};

void ak::embed_constraints(
	ak::Parameters const &parameters,
	ak::Model const &model,
	std::vector< ak::Constraint > const &constraints,
	ak::Model *constrained_model_,
	std::vector< float > *constrained_values_, //same size as out_model's vertices
	std::vector< std::vector< glm::vec3 > > *DEBUG_chain_paths,
	std::vector< std::vector< glm::vec3 > > *DEBUG_chain_loops
) {
	assert(constrained_model_);
	auto &constrained_model = *constrained_model_;
	constrained_model = ak::Model();

	assert(constrained_values_);
	auto &constrained_values = *constrained_values_;
	constrained_values = std::vector< float >();

	if (DEBUG_chain_paths) {
		*DEBUG_chain_paths = std::vector< std::vector< glm::vec3 > >(constraints.size());
	}

	if (DEBUG_chain_loops) {
		*DEBUG_chain_loops = std::vector< std::vector< glm::vec3 > >(constraints.size());
	}

	//No constraints => return input model:
	if (constraints.empty()) {
		constrained_model = model;
		constrained_values.assign(constrained_model.vertices.size(), std::numeric_limits< float >::quiet_NaN());
		return;
	}

	std::vector< std::vector< std::pair< uint32_t, float > > > adj(model.vertices.size());
	{ //extract edges from model:
		std::set< std::pair< uint32_t, uint32_t > > edges;
		for (auto const &tri : model.triangles) {
			edges.insert(std::minmax(tri.x, tri.y));
			edges.insert(std::minmax(tri.y, tri.z));
			edges.insert(std::minmax(tri.z, tri.x));
		}
		for (auto const &e : edges) {
			float len = glm::length(model.vertices[e.second] - model.vertices[e.first]);
			adj[e.first].emplace_back(e.second, len);
			adj[e.second].emplace_back(e.first, len);
		}
	}

	//find chain paths on original model:
	std::vector< std::vector< uint32_t > > paths;
	for (auto const &cons : constraints) {
		if (cons.chain.empty()) continue;
		std::vector< uint32_t > path;
		for (uint32_t goal : cons.chain) {
			if (path.empty()) {
				path.emplace_back(goal);
				continue;
			}
			std::vector< std::pair< float, uint32_t > > todo;
			std::vector< std::pair< float, uint32_t > > visited(model.vertices.size(), std::make_pair(std::numeric_limits< float >::infinity(), -1U));
			auto visit = [&todo, &visited](uint32_t vertex, float distance, uint32_t from) {
				if (distance < visited[vertex].first) {
					visited[vertex] = std::make_pair(distance, from);
					todo.emplace_back(distance, vertex);
					std::push_heap(todo.begin(), todo.end(), std::greater< std::pair< float, uint32_t > >());
				}
			};
			visit(goal, 0.0f, -1U);
			while (!todo.empty()) {
				std::pop_heap(todo.begin(), todo.end(), std::greater< std::pair< float, uint32_t > >());
				auto at = todo.back();
				todo.pop_back();
				if (at.first > visited[at.second].first) continue;
				if (at.second == path.back()) break;
				for (auto const &a : adj[at.second]) {
					visit(a.first, at.first + a.second, at.second);
				}
			}
			while (path.back() != goal) {
				if (visited[path.back()].second == -1U) {
					std::cerr << "ERROR: constraint chain moves between connected components." << std::endl;
					break;
				}
				path.emplace_back(visited[path.back()].second);
			}
		}
			paths.emplace_back(path);
	}

	//Now create a higher-resolution mesh for trimming / eventually interpolation:

	const float MaxEdgeLength = parameters.get_max_edge_length(); //largest allowed edge length
	constexpr const float MinEdgeRatio = 0.3f; //smallest allowed smallest-to-largest edge ratio in a triangle

	const float MaxEdgeLength2 = MaxEdgeLength * MaxEdgeLength;
	constexpr const float MinEdgeRatio2 = MinEdgeRatio * MinEdgeRatio;

	std::cout << "Max edge length: " << MaxEdgeLength << " model units." << std::endl;

	std::vector< glm::vec3 > verts = model.vertices;
	std::vector< glm::uvec3 > tris = model.triangles;

	//PARANOIA: no degenerate triangles, please
	for (auto const &tri : tris) {
		glm::vec3 const &x = verts[tri.x];
		glm::vec3 const &y = verts[tri.y];
		glm::vec3 const &z = verts[tri.z];
		assert(tri.x != tri.y && tri.x != tri.z && tri.y != tri.z);
		assert(x != y && x != z && y != z);
	}

	/*
	std::vector< EmbeddedVertex > everts;
	everts.reserve(verts.size());
	for (uint32_t i = 0; i < verts.size(); ++i) {
		everts.emplace_back(EmbeddedVertex::on_vertex(i));
	}
	*/

	auto divide = [&verts, &tris, &paths](std::unordered_set< glm::uvec2 > const &marked) {
		assert(!marked.empty());
		std::unordered_map< glm::uvec2, uint32_t > marked_verts;
		marked_verts.reserve(marked.size());

		{ //create new verts in the middle of edges:
			std::vector< glm::ivec2 > edges(marked.begin(), marked.end());
			//sort to avoid any system-specific hash ordering:
			std::sort(edges.begin(), edges.end(), [](glm::uvec2 const &a, glm::uvec2 const &b){
				if (a.x != b.x) return a.x < b.x;
				else return a.y < b.y;
			});
			for (auto const &e : edges) {
				marked_verts.insert(std::make_pair(e, verts.size()));
				verts.emplace_back((verts[e.x] + verts[e.y]) / 2.0f);
			}
		}

		auto lookup = [&marked_verts](uint32_t a, uint32_t b) {
			auto f = marked_verts.find((a < b ? glm::uvec2(a,b) : glm::uvec2(b,a)));
			if (f != marked_verts.end()) return f->second;
			else return -1U;
		};

		//subdivide all paths:
		for (auto &path : paths) {
			std::vector< uint32_t > new_path;
			new_path.emplace_back(path[0]);
			for (uint32_t i = 1; i < path.size(); ++i) {
				uint32_t v = lookup(path[i-1], path[i]);
				if (v != -1U) new_path.emplace_back(v);
				new_path.emplace_back(path[i]);
			}
			path = std::move(new_path);
		}

		//subdivide all tris:
		std::vector< glm::uvec3 > new_tris;

		auto quad = [&new_tris, &verts](uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
			float ac = glm::length2(verts[c] - verts[a]);
			float bd = glm::length2(verts[d] - verts[b]);
			if (ac < bd) {
				new_tris.emplace_back(a,b,c);
				new_tris.emplace_back(c,d,a);
			} else {
				new_tris.emplace_back(a,b,d);
				new_tris.emplace_back(b,c,d);
			}
		};
		for (auto const &tri : tris) {
			uint32_t a = tri.x;
			uint32_t b = tri.y;
			uint32_t c = tri.z;
			uint32_t ab = lookup(a,b);
			uint32_t bc = lookup(b,c);
			uint32_t ca = lookup(c,a);

			if (ab != -1U && bc != -1U && ca != -1U) {
				//1 -> 4 subdiv!
				new_tris.emplace_back(a, ab, ca);
				new_tris.emplace_back(b, bc, ab);
				new_tris.emplace_back(c, ca, bc);
				new_tris.emplace_back(ab, bc, ca);
			} else if (ab != -1U && bc != -1U && ca == -1U) {
				//1 -> 3 subdiv!
				//NOTE: should consider recursively subdividing to avoid this case
				quad(a, ab, bc, c);
				new_tris.emplace_back(ab, b, bc);
			} else if (ab != -1U && bc == -1U && ca != -1U) {
				new_tris.emplace_back(a, ab, ca);
				quad(ab, b, c, ca);
			} else if (ab == -1U && bc != -1U && ca != -1U) {
				quad(a, b, bc, ca);
				new_tris.emplace_back(bc, c, ca);
			} else if (ab != -1U && bc == -1U && ca == -1U) {
				//1 -> 2 subdiv!
				new_tris.emplace_back(a, ab, c);
				new_tris.emplace_back(b, c, ab);
			} else if (ab == -1U && bc != -1U && ca == -1U) {
				new_tris.emplace_back(a, b, bc);
				new_tris.emplace_back(bc, c, a);
			} else if (ab == -1U && bc == -1U && ca != -1U) {
				new_tris.emplace_back(a, b, ca);
				new_tris.emplace_back(b, c, ca);
			} else { assert(ab == -1U && bc == -1U && ca == -1U);
				//no subdiv!
				new_tris.emplace_back(a, b, c);
			}
		}
		tris = std::move(new_tris);
	};

	//edge length subdivision:
	while (true) {
		//mark edges for subdivision:
		std::unordered_set< glm::uvec2 > marked;
		auto mark = [&marked](uint32_t a, uint32_t b) {
			if (b < a) std::swap(a,b);
			marked.insert(glm::uvec2(a,b));
		};
		auto is_marked = [&marked](uint32_t a, uint32_t b) {
			if (b < a) std::swap(a,b);
			return marked.find(glm::uvec2(a,b)) != marked.end();
		};
		(void)is_marked;
		(void)MinEdgeRatio2;

		//mark for length:
		for (auto const &tri : tris) {
			float len_ab2 = glm::length2(verts[tri.y] - verts[tri.x]);
			float len_bc2 = glm::length2(verts[tri.z] - verts[tri.y]);
			float len_ca2 = glm::length2(verts[tri.x] - verts[tri.z]);
			if (len_ab2 > MaxEdgeLength2) mark(tri.x, tri.y);
			if (len_bc2 > MaxEdgeLength2) mark(tri.y, tri.z);
			if (len_ca2 > MaxEdgeLength2) mark(tri.z, tri.x);
		}
		/*//avoid 1->3 subdivisions:
		while (true) {
			uint32_t old_size = marked.size();
			for (auto const &tri : tris) {
				uint32_t count =
					  (is_marked(tri.x, tri.y) ? 1 : 0)
					+ (is_marked(tri.y, tri.z) ? 1 : 0)
					+ (is_marked(tri.z, tri.x) ? 1 : 0);
				if (count == 2) {
					mark(tri.x, tri.y);
					mark(tri.y, tri.z);
					mark(tri.z, tri.x);
				}
			}
			if (marked.size() == old_size) break;
		}*/

		//std::cout << "  marked " << marked.size() << " for length." << std::endl;
		/* This seems broken [makes way too many triangles]:
		if (marked.empty()) {
			//mark for ratio:
			while (true) {
				uint32_t old_size = marked.size();
				for (auto const &tri : tris) {
					float len_ab2 = glm::length2(verts[tri.y] - verts[tri.x]);
					float len_bc2 = glm::length2(verts[tri.z] - verts[tri.y]);
					float len_ca2 = glm::length2(verts[tri.x] - verts[tri.z]);
					if (is_marked(tri.x, tri.y)) len_ab2 /= 4.0f;
					if (is_marked(tri.y, tri.z)) len_bc2 /= 4.0f;
					if (is_marked(tri.z, tri.x)) len_ca2 /= 4.0f;

					if (std::min(len_bc2, len_ca2) / len_ab2 < MinEdgeRatio2) mark(tri.x, tri.y);
					if (std::min(len_ab2, len_ca2) / len_bc2 < MinEdgeRatio2) mark(tri.y, tri.z);
					if (std::min(len_ab2, len_bc2) / len_ca2 < MinEdgeRatio2) mark(tri.z, tri.x);
				}
				if (marked.size() == old_size) break;
			}
			std::cout << "  marked " << marked.size() << " for ratio." << std::endl;
		}
		*/

		if (marked.empty()) {
			break;
		}
		divide(marked);
	}
	//std::cout << "After division, have " << tris.size() << " triangles on " << verts.size() << " vertices." << std::endl;

	//PARANOIA: no degenerate triangles, please?
	for (auto const &tri : tris) {
		glm::vec3 const &x = verts[tri.x];
		glm::vec3 const &y = verts[tri.y];
		glm::vec3 const &z = verts[tri.z];
		assert(tri.x != tri.y && tri.x != tri.z && tri.y != tri.z);
		assert(x != y && x != z && y != z);
	}

	if (DEBUG_chain_paths) {
		for (auto const &path : paths) {
			auto &DEBUG_chain_path = (*DEBUG_chain_paths)[&path - &paths[0]];
			for (uint32_t v : path) {
				DEBUG_chain_path.emplace_back(verts[v]);
			}
		}
	}

	adj.assign(verts.size(), std::vector< std::pair< uint32_t, float > >());
	{ //extract edges from subdivided model:
		std::set< std::pair< uint32_t, uint32_t > > edges;
		for (auto const &tri : tris) {
			edges.insert(std::minmax(tri.x, tri.y));
			edges.insert(std::minmax(tri.y, tri.z));
			edges.insert(std::minmax(tri.z, tri.x));
		}
		for (auto const &e : edges) {
			float len = glm::length(verts[e.second] - verts[e.first]);
			adj[e.first].emplace_back(e.second, len);
			adj[e.second].emplace_back(e.first, len);
		}
	}

	std::unordered_map< glm::uvec2, uint32_t > opposite; //vertex opposite each [oriented] triangle edge
	opposite.reserve(tris.size() * 3);
	for (auto const &tri : tris) {
		auto ret_xy = opposite.insert(std::make_pair(glm::uvec2(tri.x, tri.y), tri.z));
		assert(ret_xy.second);
		auto ret_yz = opposite.insert(std::make_pair(glm::uvec2(tri.y, tri.z), tri.x));
		assert(ret_yz.second);
		auto ret_zx = opposite.insert(std::make_pair(glm::uvec2(tri.z, tri.x), tri.y));
		assert(ret_zx.second);
	}

	{ //build (+ add to adj) extra "shortcut" edges by unwrapping triangle neighborhoods:
		std::unordered_map< glm::uvec2, float > min_dis;
		auto get_dis = [&](uint32_t a, uint32_t b) -> float & {
			if (a > b) std::swap(a,b);
			return min_dis.insert(std::make_pair(glm::uvec2(a,b), std::numeric_limits< float >::infinity())).first->second;
		};
		for (auto const &tri : tris) {
			glm::vec2 flat_x, flat_y, flat_z; //original verts
			{
				glm::vec3 const &x = verts[tri.x];
				glm::vec3 const &y = verts[tri.y];
				glm::vec3 const &z = verts[tri.z];
				flat_x = glm::vec2(0.0f, 0.0f);
				flat_y = glm::vec2(glm::length(y-x), 0.0f);

				glm::vec3 xy = glm::normalize(y-x);
				glm::vec3 perp_xy = glm::normalize(glm::cross(glm::cross(y-x, z-x), y-x));
				float along = glm::dot(z-x, xy);
				float perp = glm::dot(z-x, perp_xy);

				flat_z = glm::vec2(along, perp);

				//std::cout << "x: (" << x.x << ", " << x.y << ", " << x.z << ") -> (" << flat_x.x << ", " << flat_x.y << ")" << std::endl; //DEBUG
				//std::cout << "y: (" << y.x << ", " << y.y << ", " << y.z << ") -> (" << flat_y.x << ", " << flat_y.y << ")" << std::endl; //DEBUG
				//std::cout << "z: (" << z.x << ", " << z.y << ", " << z.z << ") -> (" << flat_z.x << ", " << flat_z.y << ")" << std::endl; //DEBUG

			}

			//look through edge [ai,bi] from point [root], where edge [ai,bi] is ccw oriented.

			auto is_ccw = [](glm::vec2 const &a, glm::vec2 const &b, glm::vec2 const &c) {
				return glm::dot(glm::vec2(-(b.y-a.y),(b.x-a.x)), c-a) > 0.0f;
			};

			std::function< void(uint32_t, uint32_t, glm::vec2 const &, uint32_t, glm::vec2 const &, uint32_t, glm::vec2 const &, glm::vec2 const &, glm::vec2 const &) > unfold = [&](uint32_t depth, uint32_t root, glm::vec2 const &flat_root, uint32_t ai, glm::vec2 const &flat_a, uint32_t bi, glm::vec2 const &flat_b, glm::vec2 const &limit_a, glm::vec2 const &limit_b) {
				//std::cout << "r: " << root << ": (" << flat_root.x << ", " << flat_root.y << ")" << std::endl; //DEBUG
				//std::cout << "a: " << ai << ": (" << flat_a.x << ", " << flat_a.y << ")" << std::endl; //DEBUG
				//std::cout << "b: " << bi << ": (" << flat_b.x << ", " << flat_b.y << ")" << std::endl; //DEBUG
				assert(is_ccw(flat_root, flat_a, flat_b));
				//should go 'a - limit_a - limit_b - b':
				//assert(flat_a == limit_a || is_ccw(flat_root, flat_a, limit_a));
				assert(is_ccw(flat_root, limit_a, limit_b));
				//assert(flat_b == limit_b || is_ccw(flat_root, limit_b, flat_b));

				uint32_t ci;
				glm::vec2 flat_c;
				{ //if there is a triangle over the ai->bi edge, find other vertex and flatten it:
					auto f = opposite.find(glm::uvec2(bi, ai));
					if (f == opposite.end()) return;
					ci = f->second;
					//figure out c's position along ab and distance from ab:
					glm::vec3 const &a = verts[ai];
					glm::vec3 const &b = verts[bi];
					glm::vec3 const &c = verts[ci];

					glm::vec3 ab = glm::normalize(b-a);
					float along = glm::dot(c-a, ab);
					float perp = -glm::length(c-a - ab*along);

					glm::vec2 flat_ab = glm::normalize(flat_b - flat_a);
					glm::vec2 flat_perp_ab = glm::vec2(-flat_ab.y, flat_ab.x);

					flat_c = flat_a + flat_ab * along + flat_perp_ab * perp;
				}

				//std::cout << "c: " << ci << ": (" << flat_c.x << ", " << flat_c.y << ")" << std::endl; //DEBUG

				//flat_a and flat_b should always be outside limit, it seems like we need to test anyway (thanks, numerics)

				bool ccw_rac = is_ccw(flat_root, limit_a, flat_c) && is_ccw(flat_root, flat_a, flat_c);
				bool ccw_rcb = is_ccw(flat_root, flat_c, limit_b) && is_ccw(flat_root, flat_c, flat_b);

				if (ccw_rac && ccw_rcb) {
					float &dis = get_dis(root, ci);
					dis = std::min(dis, glm::length(flat_root - flat_c));

					//PARANOIA:
					float dis3 = glm::length(verts[root] - verts[ci]);
					if (dis3 > dis + 1e-6) {
						std::cerr << "dis3: " << dis3 << " vs flat dis " << dis << " seems bad!" << std::endl;
						std::cerr << "  ra3: " << glm::length(verts[root] - verts[ai]) << " vs ra: " << glm::length(flat_root - flat_a) << std::endl;
						std::cerr << "  rb3: " << glm::length(verts[root] - verts[bi]) << " vs rb: " << glm::length(flat_root - flat_b) << std::endl;
						std::cerr << "  ab3: " << glm::length(verts[ai] - verts[bi]) << " vs ab: " << glm::length(flat_a - flat_b) << std::endl;
						std::cerr << "  ac3: " << glm::length(verts[ai] - verts[ci]) << " vs ac: " << glm::length(flat_a - flat_c) << std::endl;
						std::cerr << "  bc3: " << glm::length(verts[bi] - verts[ci]) << " vs bc: " << glm::length(flat_b - flat_c) << std::endl;
						assert(dis3 < dis + 1e-6);
					}

					if (depth > 1) {
						assert(is_ccw(flat_root, flat_a, flat_c));
						unfold(depth - 1, root, flat_root, ai, flat_a, ci, flat_c, limit_a, flat_c);
						assert(is_ccw(flat_root, flat_c, flat_b));
						unfold(depth - 1, root, flat_root, ci, flat_c, bi, flat_b, flat_c, limit_b);
					}
				} else if (ccw_rac && !ccw_rcb) {
					if (depth > 1) {
						//assert(!is_ccw(flat_root, flat_c, limit_b)); //DEBUG
						//assert(is_ccw(flat_root, limit_b, flat_c)); //DEBUG -- fails sometimes [thanks, numerics]
						assert(is_ccw(flat_root, flat_a, flat_c));
						unfold(depth - 1, root, flat_root, ai, flat_a, ci, flat_c, limit_a, limit_b);
					}
				} else if (!ccw_rac && ccw_rcb) {
					if (depth > 1) {
						assert(is_ccw(flat_root, flat_c, flat_b));
						unfold(depth - 1, root, flat_root, ci, flat_c, bi, flat_b, limit_a, limit_b);
					}
				}
			};

			const constexpr uint32_t D = 3; //depth to unfold triangles to for more adjacency information; makes slightly nicer geodesics at the expense of increased compute time.

			if (D > 0) {
				unfold(D, tri.x, flat_x, tri.y, flat_y, tri.z, flat_z, flat_y, flat_z);
				unfold(D, tri.y, flat_y, tri.z, flat_z, tri.x, flat_x, flat_z, flat_x);
				unfold(D, tri.z, flat_z, tri.x, flat_x, tri.y, flat_y, flat_x, flat_y);
			}
		}
		for (uint32_t x = 0; x < verts.size(); ++x) {
			for (auto const &yd : adj[x]) {
				float &dis = get_dis(x, yd.first);
				dis = std::min(dis, yd.second);
			}
		}

		//clear adj + re-create from min_dis:
		uint32_t old_adj = 0;
		for (auto const &a : adj) {
			old_adj += a.size();
		}

		adj.assign(verts.size(), std::vector< std::pair< uint32_t, float > >());

		for (auto const &xyd : min_dis) {
			assert(xyd.first.x != xyd.first.y);
			adj[xyd.first.x].emplace_back(xyd.first.y, xyd.second);
			adj[xyd.first.y].emplace_back(xyd.first.x, xyd.second);
		}

		uint32_t new_adj = 0;
		for (auto const &a : adj) {
			new_adj += a.size();
		}

		//std::cout << "Went from " << old_adj << " to " << new_adj << " by unfolding triangles." << std::endl;

		//for consistency:
		for (auto &a : adj) {
			std::sort(a.begin(), a.end());
		}
	}

	//uint32_t used_edges = 0;

	std::vector< std::vector< EmbeddedVertex > > embedded_chains;

	for (auto const &cons : constraints) {
		embedded_chains.emplace_back();

		auto const &path = paths[&cons - &constraints[0]];
		if (cons.radius == 0.0f) {
			//add directly to embedded constrained edges.
			for (auto v : path) {
				assert(v < verts.size());
				embedded_chains.back().emplace_back(EmbeddedVertex::on_vertex(v));
			}
			continue;
		}
		//generate distance field from constraint:
		std::vector< std::pair< float, uint32_t > > todo;
		std::vector< float > distances(verts.size(), std::numeric_limits< float >::infinity());
		auto visit = [&todo, &distances](uint32_t vertex, float distance) {
			if (distance < distances[vertex]) {
				distances[vertex] = distance;
				todo.emplace_back(distance, vertex);
				std::push_heap(todo.begin(), todo.end(), std::greater< std::pair< float, uint32_t > >());
			}
		};
		for (uint32_t i = 0; i < path.size(); ++i) {
			visit(path[i], -cons.radius);
		}
		/*auto do_edge = [&](uint32_t ai, uint32_t bi) {
			auto f = opposite.find(glm::uvec2(ai, bi));
			if (f == opposite.end()) return;
			uint32_t ci = f->second;
			glm::vec3 const &a = verts[ai];
			glm::vec3 const &b = verts[bi];
			glm::vec3 const &c = verts[ci];
			float along = glm::dot(c - a, b - a);
			if (along <= 0.0f) return;
			float lim = glm::dot(b - a, b - a);
			if (along >= lim) return;
			//++used_edges;
			glm::vec3 close = glm::mix(a, b, along / lim);
			visit(ci, glm::length(c - close) - cons.radius);
		};
		for (uint32_t i = 1; i < path.size(); ++i) {
			do_edge(path[i-1], path[i]);
			do_edge(path[i], path[i-1]);
		}*/

		while (!todo.empty()) {
			std::pop_heap(todo.begin(), todo.end(), std::greater< std::pair< float, uint32_t > >());
			auto at = todo.back();
			todo.pop_back();
			if (at.first > distances[at.second]) continue;
			if (at.first > 0.0f) break; //once we start expanding things that are past the contour, no need to continue (TODO: consider blur radius)
			for (auto const &a : adj[at.second]) {
				visit(a.first, at.first + a.second);
			}
		}

		//read back embedded path.

		std::unordered_map< glm::uvec2, EmbeddedVertex > embedded_pts;
		std::unordered_map< glm::uvec2, glm::vec3 > pts;
		auto add = [&distances,&verts,&pts,&embedded_pts](uint32_t a, uint32_t b) {
			assert(distances[a] < 0.0f && distances[b] >= 0.0f);
			float mix = (0.0f - distances[a]) / (distances[b] - distances[a]);
			pts[glm::uvec2(a,b)] = glm::mix(verts[a], verts[b], mix);
			embedded_pts[glm::uvec2(a,b)] = EmbeddedVertex::on_edge(a,b,mix);
			return glm::uvec2(a,b);
		};
		std::unordered_map< glm::uvec2, glm::uvec2 > links;
		std::unordered_map< glm::uvec2, glm::uvec2 > back_links;
		auto link = [&links,&back_links](glm::uvec2 f, glm::uvec2 t) {
			auto res = links.insert(std::make_pair(f, t));
			assert(res.second);
			auto res2 = back_links.insert(std::make_pair(t, f));
			assert(res2.second);
		};
		for (auto const &tri : tris) {
			uint32_t a = tri.x;
			uint32_t b = tri.y;
			uint32_t c = tri.z;
			//spin triangle until 'a' is the minimum distance value:
			for (uint32_t i = 0; i < 3; ++i) {
				if (distances[a] <= distances[b] && distances[a] <= distances[c]) break;
				uint32_t t = a; a = b; b = c; c = t;
			}
			//NOTE: we treat 0.0f as "0.0f + epsilon"
			if (distances[a] >= 0.0f) continue; //all above border
			assert(distances[a] < 0.0f);

			if (distances[b] >= 0.0f && distances[c] >= 0.0f) {
				//edge is from ab to ca
				link(add(a,b), add(a,c));
			} else if (distances[b] >= 0.0f && distances[c] < 0.0f) {
				//edge is from ab to bc
				link(add(a,b), add(c,b));
			} else if (distances[b] < 0.0f && distances[c] >= 0.0f) {
				//edge is from bc to ca
				link(add(b,c), add(a,c));
			} else {
				assert(distances[b] < 0.0f && distances[c] < 0.0f);
				//all below border, nothing to do.
			}
		}

		//read back path from links:
		if (!links.empty()) {
			std::deque< glm::uvec2 > loop;
			loop.emplace_back(links.begin()->first);
			while (true) {
				auto f = links.find(loop.back());
				if (f == links.end()) break;
				loop.emplace_back(f->second);
				if (f->second == loop[0]) break;
			}
			if (loop[0] != loop.back()) {
				while (true) {
					auto f = back_links.find(loop[0]);
					if (f == back_links.end()) break;
					if (f->second == loop.back()) break;
					loop.emplace_front(f->second);
				}
			}

			for (glm::uvec2 e : loop) {
				auto f = embedded_pts.find(e);
				assert(f != embedded_pts.end());
				embedded_chains.back().emplace_back(f->second);
			}

			if (DEBUG_chain_loops) {
				auto &DEBUG_chain_loop = (*DEBUG_chain_loops)[&cons - &constraints[0]];
				for (glm::uvec2 e : loop) {
					auto f = pts.find(e);
					assert(f != pts.end());
					DEBUG_chain_loop.emplace_back(f->second);
				}
			}

		}
	}

	//should have a chain per constraint:
	assert(embedded_chains.size() == constraints.size());

	//embed chains using planar map:
	EmbeddedPlanarMap epm;
	uint32_t total_chain_edges = 0;
	for (uint32_t c = 0; c < constraints.size(); ++c) {
		uint32_t first = 0;
		uint32_t last = 0;
		for (uint32_t i = 0; i + 1 < embedded_chains[c].size(); ++i) {
			uint32_t a = epm.add_vertex(embedded_chains[c][i]);
			uint32_t b = epm.add_vertex(embedded_chains[c][i+1]);
			epm.add_edge(a,b,constraints[c].value);
			++total_chain_edges;
			if (i == 0) first = a;
			if (i + 2 == embedded_chains[c].size()) last = b;
		}
		if (first != last) std::cout << "NOTE: have open chain." << std::endl;
	}
	uint32_t total_simplex_edges = 0;
	for (const auto &edges : epm.simplex_edges) {
		total_simplex_edges += edges.second.size();
	}
	std::cout << "EPM has " << epm.vertices.size() << " vertices." << std::endl;
	std::cout << "EPM has " << epm.simplex_vertices.size() << " simplices with vertices." << std::endl;
	std::cout << "EPM has " << epm.simplex_edges.size() << " simplices with edges (" << total_simplex_edges << " edges from " << total_chain_edges << " chain edges)." << std::endl;

	/*//DEBUG:
	for (const auto &se : epm.simplex_edges) {
		assert(se.first.x <= se.first.y && se.first.y <= se.first.z);
		if (se.first.z != -1U) {
			std::cout << se.first.x << ", " << se.first.y << ", " << se.first.z << std::endl;
		}
	}*/

	{ //Build a mesh that is split at the embedded edges:
		std::vector< glm::vec3 > split_verts = verts;
		std::vector< glm::uvec3 > split_tris;
		std::vector< uint32_t > epm_to_split;

		epm_to_split.reserve(epm.vertices.size());
		for (uint32_t i = 0; i < epm.vertices.size(); ++i) {
			const auto &v = epm.vertices[i];
			assert(v.simplex.x < verts.size());
			assert(v.weights.x + v.weights.y + v.weights.z == WEIGHT_SUM);
			if (v.simplex.y == -1U) {
				assert(v.simplex.z == -1U);
				epm_to_split.emplace_back(v.simplex.x);
			} else if (v.simplex.z == -1U) {
				assert(v.simplex.y < verts.size());
				assert(v.weights.z == 0);
				epm_to_split.emplace_back(split_verts.size());
				split_verts.emplace_back((
					  verts[v.simplex.x] * float(v.weights.x)
					+ verts[v.simplex.y] * float(v.weights.y)
					) / float(WEIGHT_SUM));
			} else {
				assert(v.simplex.y < verts.size());
				assert(v.simplex.z < verts.size());
				epm_to_split.emplace_back(split_verts.size());
				split_verts.emplace_back((
					  verts[v.simplex.x] * float(v.weights.x)
					+ verts[v.simplex.y] * float(v.weights.y)
					+ verts[v.simplex.z] * float(v.weights.z)
					) / float(WEIGHT_SUM));
			}
		}


		for (const auto &tri : tris) {
			glm::uvec3 simplex = tri;
			bool need_flip = false;
			if (simplex.x > simplex.y) {
				std::swap(simplex.x, simplex.y);
				need_flip = !need_flip;
			}
			if (simplex.y > simplex.z) {
				std::swap(simplex.y, simplex.z);
				need_flip = !need_flip;
			}
			if (simplex.x > simplex.y) {
				std::swap(simplex.x, simplex.y);
				need_flip = !need_flip;
			}
			assert(simplex.x < simplex.y && simplex.y < simplex.z);
			/*if (!epm.simplex_edges.count(simplex)) {
				//DEBUG, was: split_tris.emplace_back(tri);
				continue;
			}
			assert(!epm.simplex_edges[simplex].empty());
			*/
			std::vector< uint32_t > source_verts; //in split_verts
			std::vector< glm::ivec2 > coords; //[x,y] in weight space
			std::unordered_multimap< uint32_t, uint32_t > half_edges; //relative to coords

			std::unordered_map< uint32_t, uint32_t > split_verts_to_coords;
			auto ref_vert = [&](uint32_t split_vert, const glm::ivec3 &weights_on) -> uint32_t {
				auto res = split_verts_to_coords.insert(std::make_pair(split_vert, coords.size()));
				if (res.second) {
					//std::cout << "Inserting " << split_vert << " / " << weights_on.x << ", " << weights_on.y << ", " << weights_on.z << std::endl; //DEBUG
					assert(res.first->second == coords.size());
					coords.emplace_back(weights_on);
					source_verts.emplace_back(split_vert);
				}
				return res.first->second;
			};

			//add two half-edges for all internal edges:
			//std::cout << " ---- internal edges ---- " << std::endl; //DEBUG
			auto f = epm.simplex_edges.find(simplex);
			if (f != epm.simplex_edges.end()) {
				for (const auto &edge : f->second) {
					uint32_t a = ref_vert(epm_to_split[edge.first], epm.vertices[edge.first].weights_on(simplex));
					uint32_t b = ref_vert(epm_to_split[edge.second], epm.vertices[edge.second].weights_on(simplex));
					half_edges.insert(std::make_pair(a,b));
					half_edges.insert(std::make_pair(b,a));
				}
			}
			//std::cout << " ---- sides ---- " << std::endl; //DEBUG
			//add one half-edge along all sides:
			auto do_side = [&](uint32_t ifrom, uint32_t ito) {
				glm::uvec3 edge_simplex(simplex[ifrom], simplex[ito], -1U);
				if (edge_simplex.x > edge_simplex.y) std::swap(edge_simplex.x, edge_simplex.y);
				uint32_t from;
				uint32_t to;
				{
					glm::ivec3 from_weights = glm::ivec3(0);
					from_weights[ifrom] = WEIGHT_SUM;
					glm::ivec3 to_weights = glm::ivec3(0);
					to_weights[ito] = WEIGHT_SUM;
					from = ref_vert(simplex[ifrom], from_weights);
					to = ref_vert(simplex[ito], to_weights);
				}
				
				std::map< int32_t, uint32_t > weight_to_vert;
				auto f = epm.simplex_vertices.find(edge_simplex);
				if (f != epm.simplex_vertices.end()) {
					for (auto sv : f->second) {
						glm::ivec3 weights = epm.vertices[sv].weights_on(simplex);
						auto res = weight_to_vert.insert(std::make_pair(weights[ito], ref_vert(epm_to_split[sv], weights)));
						assert(res.second);
					}
				}
				weight_to_vert.insert(std::make_pair(0, from));
				weight_to_vert.insert(std::make_pair(WEIGHT_SUM, to));

				assert(weight_to_vert[0] == from);
				assert(weight_to_vert[WEIGHT_SUM] == to);
				auto at = weight_to_vert.begin();
				uint32_t prev = at->second;
				++at;
				for (; at != weight_to_vert.end(); ++at) {
					half_edges.insert(std::make_pair(prev, at->second));
					prev = at->second;
				}
			};
			do_side(0,1);
			do_side(1,2);
			do_side(2,0);

			/*//DEBUG:
			std::cout << "---------------\n";
			for (uint32_t c = 0; c < coords.size(); ++c) {
				std::cout << "coords[" << c << "] = (" << coords[c].x << ", " << coords[c].y << ")\n";
			}
			for (auto e : half_edges) {
				std::cout << "  " << e.first << " -> " << e.second << "\n";
			}
			std::cout.flush();*/

			while (!half_edges.empty()) {
				std::vector< uint32_t > loop;
				loop.emplace_back(half_edges.begin()->first);
				loop.emplace_back(half_edges.begin()->second);

				//std::cout << "Seed: " << loop[0] << " -> " << loop[1] << std::endl; //DEBUG
				bool had_reflex = false;

				while (true) {
					const glm::ivec2 &from = coords[loop[loop.size()-2]];
					const glm::ivec2 &at = coords[loop.back()];
					//want to find he ccw-most exit, relative to from->at
					auto r = half_edges.equal_range(loop.back());
					assert(r.first != r.second); //ran out of edges to make a loop(?!?!)
					uint32_t best_to = -1U;
					glm::ivec2 best_d = glm::ivec2(0,0); //best == largest y/x
					int best_quad = 4;
					auto best_iter = half_edges.end();
					for (auto ri = r.first; ri != r.second; ++ri) {
						uint32_t i = ri->second;
						const glm::ivec2 &next = coords[i];
						glm::ivec2 d;
						d.x = (next - at).x * (at - from).x + (next - at).y * (at - from).y;
						d.y = (next - at).x *-(at - from).y + (next - at).y * (at - from).x;
						assert(!(d.x == 0 && d.y == 0));
						int quad;
						if (d.x <= 0 && d.y > 0) {
							quad = 0;
							d.x = -d.x;
							std::swap(d.x, d.y);
						} else if (d.x > 0 && d.y >= 0) {
							quad = 1;
						} else if (d.x >= 0 && d.y < 0) {
							d.y = -d.y;
							std::swap(d.x, d.y);
							quad = 2;
						} else if (d.x < 0 && d.y <= 0) {
							d.x = -d.x;
							d.y = -d.y;
							quad = 3;
						} else {
							assert(false);
						}

						//std::cout << "  " << ri->first << " -> " << ri->second << " quad " << quad << " y/x " << d.y << "/" << d.x << std::endl; //DEBUG

						
						if (quad < best_quad || (quad == best_quad && uint64_t(d.y) * uint64_t(best_d.x) > uint64_t(best_d.y) * uint64_t(d.x))) {
							best_to = i;
							best_quad = quad;
							best_d = d;
							best_iter = ri;
						}
					}
					assert(best_to != -1U);
					if (best_quad >= 2) had_reflex = true;
					//std::cout << "Best is " << best_iter->first << " to " << best_iter->second << std::endl; //DEBUG
					if (best_iter->first == loop[0] && best_iter->second == loop[1]) {
						half_edges.erase(best_iter);
						break;
					}
					half_edges.erase(best_iter);
					loop.emplace_back(best_to);
				}
				//std::cout << "Peeled a loop of " << loop.size() << " verts." << std::endl; //DEBUG
				assert(loop[0] == loop.back());
				loop.pop_back();
				assert(loop.size() >= 3);

				if (had_reflex) {
					std::cerr << "ERROR: no code for reflex verts in planar map yet." << std::endl;
				}
				for (uint32_t i = 1; i + 1 < loop.size(); ++i) {
					split_tris.emplace_back(source_verts[loop[0]], source_verts[loop[i]], source_verts[loop[i+1]]);
					if (need_flip) {
						std::swap(split_tris.back().y, split_tris.back().z);
					}
				}
			}
		}

		//record constrained edges in terms of split_verts:
		std::unordered_map< glm::uvec2, float > constrained_edges;
		std::vector< float > split_values(split_verts.size(), std::numeric_limits< float >::quiet_NaN());
		for (const auto &se : epm.simplex_edges) {
			for (auto const &e : se.second) {
				glm::uvec2 ab = glm::uvec2(epm_to_split[e.first], epm_to_split[e.second]);
				if (ab.x > ab.y) std::swap(ab.x, ab.y);
				constrained_edges.insert(std::make_pair(ab, e.value));
				//also grab vertex values:
				split_values[epm_to_split[e.first]] = e.value;
				split_values[epm_to_split[e.second]] = e.value;
			}
		}
		std::cout << constrained_edges.size() << " constrained edges." << std::endl;

		
		std::vector< uint32_t > tri_component(split_tris.size(), -1U);
		std::vector< bool > component_keep;
		{ //mark connected components + delete the "wrong" ones
			std::unordered_map< glm::uvec2, uint32_t > over;
			for (const auto &tri : split_tris) {
				uint32_t ti = &tri - &split_tris[0];
				auto res = over.insert(std::make_pair(glm::uvec2(tri.x, tri.y), ti));
				assert(res.second);
				res = over.insert(std::make_pair(glm::uvec2(tri.y, tri.z), ti));
				assert(res.second);
				res = over.insert(std::make_pair(glm::uvec2(tri.z, tri.x), ti));
				assert(res.second);
			}
			for (uint32_t seed = 0; seed < split_tris.size(); ++seed) {
				if (tri_component[seed] != -1U) continue;
				//std::cout << "Doing CC with seed " << seed << std::endl; //DEBUG
				uint32_t component = component_keep.size();
				tri_component[seed] = component;
				std::set< float > values;
				std::vector< uint32_t > todo;
				todo.emplace_back(seed);
				auto do_edge = [&](uint32_t a, uint32_t b) {
					{ //if edge is constrained, don't traverse over:
						glm::uvec2 e(a,b);
						if (e.x > e.y) std::swap(e.x,e.y);
						auto v = constrained_edges.find(e);
						if (v != constrained_edges.end()) {
							values.insert(v->second);
							return;
						}
					}
					//otherwise, traverse over:
					auto f = over.find(glm::uvec2(b,a));
					if (f != over.end()) {
						if (tri_component[f->second] != component) {
							assert(tri_component[f->second] == -1U);
							tri_component[f->second] = component;
							todo.emplace_back(f->second);
						}
					}
				};
				while (!todo.empty()) {
					uint32_t at = todo.back();
					todo.pop_back();
					assert(tri_component[at] == component);
					do_edge(split_tris[at].x, split_tris[at].y);
					do_edge(split_tris[at].y, split_tris[at].z);
					do_edge(split_tris[at].z, split_tris[at].x);
				}
				component_keep.emplace_back(values.size() > 1);
			}
			std::cout << "Have " << component_keep.size() << " connected components." << std::endl;
		}


		//remove any split_verts that aren't used:
		std::vector< glm::vec3 > compressed_verts;
		std::vector< float > compressed_values;
		std::vector< glm::uvec3 > compressed_tris;
		for (uint32_t ti = 0; ti < split_tris.size(); ++ti) {
			if (component_keep[tri_component[ti]]) {
				compressed_tris.emplace_back(split_tris[ti]);
			}
		}
		compressed_verts.reserve(split_verts.size());
		std::vector< uint32_t > to_compressed(split_verts.size(), -1U);
		auto add_vert = [&](uint32_t vi) {
			if (to_compressed[vi] == -1U) {
				to_compressed[vi] = compressed_verts.size();
				compressed_verts.emplace_back(split_verts[vi]);
				compressed_values.emplace_back(split_values[vi]);
			}
			return to_compressed[vi];
		};
		for (auto &tri : compressed_tris) {
			tri.x = add_vert(tri.x);
			tri.y = add_vert(tri.y);
			tri.z = add_vert(tri.z);
		}

		std::cout << "Went from " << tris.size() << " to (via split) " << split_tris.size() << " to (via discard) " << compressed_tris.size() << " triangles." << std::endl; //DEBUG

		constrained_model.vertices = compressed_verts;
		constrained_model.triangles = compressed_tris;

		constrained_values = compressed_values;

	}


	//std::cout << "Used " << used_edges << " edges." << std::endl; //DEBUG

	//TODO: split at distance field level set
	//TODO: constrain values at distance field border

/*
	constrained_model.vertices = verts;
	constrained_model.triangles = tris;
	*/

/*

	constrained_values.assign(constrained_model.vertices.size(), std::numeric_limits< float >::quiet_NaN());
	//Quick hack to test interpolation:
	uint32_t lowest = 0;
	uint32_t highest = 0;
	for (uint32_t i = 0; i < constrained_model.vertices.size(); ++i) {
		if (constrained_model.vertices[i].z < constrained_model.vertices[lowest].z) lowest = i;
		if (constrained_model.vertices[i].z > constrained_model.vertices[highest].z) highest = i;
	}
	if (!constrained_values.empty()) {
		constrained_values[lowest] =-1.0f;
		constrained_values[highest] = 1.0f;
	}
*/

}
