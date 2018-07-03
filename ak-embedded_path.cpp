#include "pipeline.hpp"

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <glm/gtx/hash.hpp>

void ak::embedded_path(
	ak::Parameters const &parameters,
	ak::Model const &model,
	ak::EmbeddedVertex const &source,
	ak::EmbeddedVertex const &target,
	std::vector< ak::EmbeddedVertex > *path_ //out: path; path[0] will be source and path.back() will be target
) {
	assert(source != target);

	assert(path_);
	auto &path = *path_;
	path.clear();


	//idea: place distance storage along each edge and on each corner.
	std::vector< ak::EmbeddedVertex > loc_ev;
	std::vector< glm::vec3 > loc_pos;

	//TODO: eventually, incrementally build locs starting from source / target verts

	//add locs for each vertex in the mesh:
	std::vector< uint32_t > vertex_locs;
	vertex_locs.reserve(model.vertices.size());
	for (uint32_t vi = 0; vi < model.vertices.size(); ++vi) {
		vertex_locs.emplace_back(loc_ev.size());
		loc_ev.emplace_back(EmbeddedVertex::on_vertex(vi));
		loc_pos.emplace_back(loc_ev.back().interpolate(model.vertices));
	}

	//make an edge-to-triangle look-up structure:
	std::unordered_multimap< glm::uvec2, uint32_t > edge_triangles;
	std::unordered_map< glm::uvec3, uint32_t > simplex_triangle;
	std::unordered_set< glm::uvec2 > edges;
	for (auto const &tri : model.triangles) {
		uint32_t ti = &tri - &model.triangles[0];
		auto do_edge = [&](uint32_t a, uint32_t b) {
			if (a > b) std::swap(a,b);
			edge_triangles.insert(std::make_pair(glm::uvec2(a,b), ti));
			edges.insert(glm::uvec2(a,b));
		};
		do_edge(tri.x, tri.y);
		do_edge(tri.y, tri.z);
		do_edge(tri.z, tri.x);

		glm::uvec3 simplex = tri;
		if (simplex.x > simplex.y) std::swap(simplex.x, simplex.y);
		if (simplex.y > simplex.z) std::swap(simplex.y, simplex.z);
		if (simplex.x > simplex.y) std::swap(simplex.x, simplex.y);
		auto ret = simplex_triangle.insert(std::make_pair(simplex, ti));
		assert(ret.second);
	}

	//add (several?) locs along each edge:
	std::unordered_map< glm::uvec2, std::pair< uint32_t, uint32_t > > edge_locs;
	edge_locs.reserve(edges.size());
	float const max_spacing = parameters.get_max_path_sample_spacing();
	for (auto const &e : edges) {
		uint32_t count = std::max(0, int32_t(std::floor(glm::length(model.vertices[e.y] - model.vertices[e.x]) / max_spacing)));
		uint32_t begin = loc_ev.size();
		uint32_t end = begin + count;
		edge_locs.insert(std::make_pair(e, std::make_pair(begin, end)));
		for (uint32_t i = 0; i < count; ++i) {
			loc_ev.emplace_back(EmbeddedVertex::on_edge(e.x, e.y, (i + 0.5f) / float(count)));
			loc_pos.emplace_back(loc_ev.back().interpolate(model.vertices));
		}
		assert(loc_ev.size() == end);
	}

	//build adjacency lists for each triangle:
	std::vector< std::vector< uint32_t > > loc_tris(loc_ev.size());
	std::vector< std::vector< uint32_t > > tri_adj(model.triangles.size());
	for (auto const &tri : model.triangles) {
		uint32_t ti = &tri - &model.triangles[0];
		auto do_edge = [&](uint32_t a, uint32_t b) {
			if (a > b) std::swap(b,a);
			auto f = edge_locs.find(glm::uvec2(a,b));
			assert(f != edge_locs.end());
			for (uint32_t i = f->second.first; i < f->second.second; ++i) {
				loc_tris[i].emplace_back(ti);
				tri_adj[ti].emplace_back(i);
			}
		};

		auto do_vertex = [&](uint32_t a) {
			uint32_t i = vertex_locs[a];
			loc_tris[i].emplace_back(ti);
			tri_adj[ti].emplace_back(i);
		};

		do_edge(tri.x, tri.y);
		do_edge(tri.y, tri.z);
		do_edge(tri.z, tri.x);
		do_vertex(tri.x);
		do_vertex(tri.y);
		do_vertex(tri.z);
	}


	//add source and target to the locs lists:
	auto add_embedded = [&](ak::EmbeddedVertex const &ev) {
		assert(ev.simplex.x != -1U);
		if (ev.simplex.y == -1U) {
			//at a vertex
			return vertex_locs[ev.simplex.x];
		} else if (ev.simplex.z == -1U) {
			//on an edge
			uint32_t idx = loc_ev.size();
			loc_ev.emplace_back(ev);
			loc_pos.emplace_back(loc_ev.back().interpolate(model.vertices));

			loc_tris.emplace_back();
			auto r = edge_triangles.equal_range(glm::uvec2(ev.simplex));
			assert(r.first != r.second);
			for (auto ri = r.first; ri != r.second; ++ri) {
				uint32_t ti = ri->second;
				loc_tris.back().emplace_back(ti);
				tri_adj[ti].emplace_back(idx);
			}

			return idx;
		} else {
			//on a triangle
			uint32_t idx = loc_ev.size();
			loc_ev.emplace_back(ev);
			loc_pos.emplace_back(loc_ev.back().interpolate(model.vertices));

			auto f = simplex_triangle.find(ev.simplex);
			assert(f != simplex_triangle.end());
			uint32_t ti = f->second;

			loc_tris.emplace_back();
			loc_tris.back().emplace_back(ti);
			tri_adj[ti].emplace_back(idx);

			return idx;
		}
	};

	uint32_t source_idx = add_embedded(source);
	uint32_t target_idx = add_embedded(target);

	/*//DEBUG:
	std::cout << "Source is loc " << source_idx << " with adj\n";
	for (auto t : loc_tris[source_idx]) {
		std::cout << " " << t << ":";
		for (auto a : tri_adj[t]) {
			std::cout << " " << a;
		}
		std::cout << "\n";
	}
	std::cout.flush();*/

	//now do actual search:

	std::vector< float > loc_dis(loc_pos.size(), std::numeric_limits< float >::infinity());
	std::vector< uint32_t > loc_from(loc_pos.size(), -1U);

	glm::vec3 target_pos = target.interpolate(model.vertices);

	std::vector< std::pair< float, std::pair< uint32_t, float > > > todo;

	auto queue = [&](uint32_t at, float distance, uint32_t from) {
		assert(distance < loc_dis[at]);
		loc_dis[at] = distance;
		loc_from[at] = from;

		float heuristic = glm::length(target_pos - loc_pos[at]);
		todo.emplace_back(std::make_pair(-(heuristic + distance), std::make_pair(at, distance)));
		std::push_heap(todo.begin(), todo.end());
	};

	queue(source_idx, 0.0f, -1U);
	while (!todo.empty()) {
		std::pop_heap(todo.begin(), todo.end());
		uint32_t at = todo.back().second.first;
		float distance = todo.back().second.second;
		todo.pop_back();

		if (distance > loc_dis[at]) continue;
		if (at == target_idx) break; //bail out early -- don't need distances to everything.

		assert(distance == loc_dis[at]);
		for (auto t : loc_tris[at]) {
			for (auto n : tri_adj[t]) {
				if (n == at) continue;
				float d = distance + glm::length(loc_pos[n] - loc_pos[at]);
				if (d < loc_dis[n]) queue(n, d, at);
			}
		}
	}

	//read back path:
	if (loc_from[target_idx] == -1U) {
		throw std::runtime_error("embedded_path requested between disconnected vertices");
	}

	uint32_t at = target_idx;
	do {
		path.emplace_back(loc_ev[at]);
		at = loc_from[at];
	} while (at != -1U);
	assert(path.size() >= 2);
	std::reverse(path.begin(), path.end());

	assert(path[0] == source);
	assert(path.back() == target);
}

