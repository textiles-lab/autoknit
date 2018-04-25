#include "pipeline.hpp"

#include <glm/gtx/norm.hpp>
#include <glm/gtx/hash.hpp>

#include <algorithm>
#include <deque>
#include <iostream>
#include <set>
#include <unordered_map>
#include <unordered_set>

void ak::embed_constraints(
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

	constexpr const float MaxEdgeLength = 0.1f; //largest allowed edge length
	constexpr const float MinEdgeRatio = 0.3f; //smallest allowed smallest-to-largest edge ratio in a triangle

	constexpr const float MaxEdgeLength2 = MaxEdgeLength * MaxEdgeLength;
	constexpr const float MinEdgeRatio2 = MinEdgeRatio * MinEdgeRatio;

	std::vector< glm::vec3 > verts = model.vertices;
	std::vector< glm::uvec3 > tris = model.triangles;
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


	for (auto const &cons : constraints) {
		auto const &path = paths[&cons - &constraints[0]];
		if (cons.radius == 0.0f) {
			//add directly to embedded constrained edges.
			break;
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

		std::unordered_map< glm::uvec2, glm::vec3 > pts;
		auto add = [&distances,&verts,&pts](uint32_t a, uint32_t b) {
			assert(distances[a] < 0.0f && distances[b] >= 0.0f);
			float mix = (0.0f - distances[a]) / (distances[b] - distances[a]);
			pts[glm::uvec2(a,b)] = glm::mix(verts[a], verts[b], mix);
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
				if (f->second == loop[0]) {
					loop.emplace_back(f->second);
					break;
				}
				loop.emplace_back(f->second);
			}
			if (loop[0] != loop.back()) {
				while (true) {
					auto f = back_links.find(loop[0]);
					if (f == back_links.end()) break;
					if (f->second == loop.back()) break;
					loop.emplace_front(f->second);
				}
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

	//TODO: split at distance field level set
	//TODO: constrain values at distance field border

	constrained_model.vertices = verts;
	constrained_model.triangles = tris;



}
