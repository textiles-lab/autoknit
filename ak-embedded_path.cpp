#include "pipeline.hpp"

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <glm/gtx/hash.hpp>
#include <glm/gtx/norm.hpp>

void embedded_path_simple(
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
		loc_ev.emplace_back(ak::EmbeddedVertex::on_vertex(vi));
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
			loc_ev.emplace_back(ak::EmbeddedVertex::on_edge(e.x, e.y, (i + 0.5f) / float(count)));
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


void ak::embedded_path(
	ak::Parameters const &parameters,
	ak::Model const &model,
	ak::EmbeddedVertex const &source,
	ak::EmbeddedVertex const &target,
	std::vector< ak::EmbeddedVertex > *path_ //out: path; path[0] will be source and path.back() will be target
) {

	assert(path_);
	auto &path = *path_;
	path.clear();

	//first do a vertex-to-vertex distance computation to bound the computation:

	std::vector< std::vector< uint32_t > > adj(model.vertices.size());

	std::unordered_set< glm::uvec2 > edges;
	for (auto const &tri : model.triangles) {
		auto do_edge = [&](uint32_t a, uint32_t b) {
			if (a > b) std::swap(a,b);
			edges.insert(glm::uvec2(a,b));
		};
		do_edge(tri.x, tri.y);
		do_edge(tri.y, tri.z);
		do_edge(tri.z, tri.x);
	}

	for (auto const &e : edges) {
		adj[e.x].emplace_back(e.y);
		adj[e.y].emplace_back(e.x);
	}

	uint32_t target_idx = target.simplex.x;

	std::vector< float > dis(model.vertices.size(), std::numeric_limits< float >::infinity());

	std::vector< std::pair< float, std::pair< uint32_t, float > > > todo;

	auto queue = [&](uint32_t at, float distance) {
		assert(distance < dis[at]);
		dis[at] = distance;

		float heuristic = glm::length(model.vertices[target_idx] - model.vertices[at]);
		todo.emplace_back(std::make_pair(-(heuristic + distance), std::make_pair(at, distance)));
		std::push_heap(todo.begin(), todo.end());
	};

	queue(source.simplex.x, glm::length(source.interpolate(model.vertices) - model.vertices[source.simplex.x]));

	while (!todo.empty()) {
		std::pop_heap(todo.begin(), todo.end());
		uint32_t at = todo.back().second.first;
		float distance = todo.back().second.second;
		todo.pop_back();

		if (distance > dis[at]) continue;
		assert(distance == dis[at]);

		if (at == target_idx) break; //bail out early -- don't need distances to everything.

		for (auto n : adj[at]) {
			float d = distance + glm::length(model.vertices[n] - model.vertices[at]);
			if (d < dis[n]) queue(n, d);
		}
	}

	//okay, so this is a conservative (long) estimate of path length:
	float dis2 = dis[target_idx] + glm::length(target.interpolate(model.vertices) - model.vertices[target_idx]);
	dis2 = dis2*dis2;

	//come up with a model containing only triangles that might be used in the path:

	Model trimmed;
	trimmed.vertices.reserve(model.vertices.size());
	trimmed.triangles.reserve(model.triangles.size());

	std::vector< uint32_t > to_trimmed(model.vertices.size(), -1U);
	std::vector< uint32_t > from_trimmed;
	from_trimmed.reserve(model.vertices.size());
	auto vertex_to_trimmed = [&to_trimmed,&from_trimmed,&trimmed,&model](uint32_t v) {
		if (to_trimmed[v] == -1U) {
			to_trimmed[v] = trimmed.vertices.size();
			from_trimmed.emplace_back(v);
			trimmed.vertices.emplace_back(model.vertices[v]);
		}
		return to_trimmed[v];
	};

	{ //keep triangles that are close enough to source and target that the path could possible pass through them:
		glm::vec3 src = source.interpolate(model.vertices);
		glm::vec3 tgt = target.interpolate(model.vertices);
		for (auto const &tri : model.triangles) {
			glm::vec3 min = glm::min(model.vertices[tri.x], glm::min(model.vertices[tri.y], model.vertices[tri.z]));
			glm::vec3 max = glm::max(model.vertices[tri.x], glm::max(model.vertices[tri.y], model.vertices[tri.z]));

			float len2_src = glm::length2(glm::max(min, glm::min(max, src)) - src);
			float len2_tgt = glm::length2(glm::max(min, glm::min(max, tgt)) - tgt);
			if (len2_src + len2_tgt < dis2) {
				trimmed.triangles.emplace_back(glm::uvec3(
					vertex_to_trimmed(tri.x), vertex_to_trimmed(tri.y), vertex_to_trimmed(tri.z)
				));
			}
		}
	}

	ak::EmbeddedVertex trimmed_source = source;
		
	trimmed_source.simplex.x = vertex_to_trimmed(trimmed_source.simplex.x);
	if (trimmed_source.simplex.y != -1U) trimmed_source.simplex.y = vertex_to_trimmed(trimmed_source.simplex.y);
	if (trimmed_source.simplex.z != -1U) trimmed_source.simplex.z = vertex_to_trimmed(trimmed_source.simplex.z);
	trimmed_source = ak::EmbeddedVertex::canonicalize(trimmed_source.simplex, trimmed_source.weights);

	ak::EmbeddedVertex trimmed_target = target;
	trimmed_target.simplex.x = vertex_to_trimmed(trimmed_target.simplex.x);
	if (trimmed_target.simplex.y != -1U) trimmed_target.simplex.y = vertex_to_trimmed(trimmed_target.simplex.y);
	if (trimmed_target.simplex.z != -1U) trimmed_target.simplex.z = vertex_to_trimmed(trimmed_target.simplex.z);
	trimmed_target = ak::EmbeddedVertex::canonicalize(trimmed_target.simplex, trimmed_target.weights);

	assert(from_trimmed.size() == trimmed.vertices.size());

	/*//DEBUG:
	std::cout << "Trimmed has " << trimmed.vertices.size() << " verts and " << trimmed.triangles.size() << " tris." << std::endl;

	std::cout << "source on: " << (int)source.simplex.x << ", " << (int)source.simplex.y << ", " << (int)source.simplex.z << std::endl;
	std::cout << "trimmed_source on: " << (int)trimmed_source.simplex.x << ", " << (int)trimmed_source.simplex.y << ", " << (int)trimmed_source.simplex.z << std::endl;
	std::cout << "target on: " << (int)target.simplex.x << ", " << (int)target.simplex.y << ", " << (int)target.simplex.z << std::endl;
	std::cout << "trimmed_target on: " << (int)trimmed_target.simplex.x << ", " << (int)trimmed_target.simplex.y << ", " << (int)trimmed_target.simplex.z << std::endl;

	//PARANOIA:
	bool found_source = false;
	bool found_target = false;
	for (auto simplex : trimmed.triangles) {
		if (simplex.x > simplex.y) std::swap(simplex.x, simplex.y);
		if (simplex.y > simplex.z) std::swap(simplex.y, simplex.z);
		if (simplex.x > simplex.y) std::swap(simplex.x, simplex.y);

		if (simplex == trimmed_source.simplex) found_source = true;
		if (simplex == trimmed_target.simplex) found_target = true;
		simplex.z = -1U;
		if (simplex == trimmed_source.simplex) found_source = true;
		if (simplex == trimmed_target.simplex) found_target = true;
		simplex.y = -1U;
		if (simplex == trimmed_source.simplex) found_source = true;
		if (simplex == trimmed_target.simplex) found_target = true;
	}
	assert(found_source);
	assert(found_target);
	*/


	embedded_path_simple(
		parameters,
		trimmed,
		trimmed_source,
		trimmed_target,
		&path);

	for (auto &v : path) {
		v.simplex.x = from_trimmed[v.simplex.x];
		if (v.simplex.y != -1U) v.simplex.y = from_trimmed[v.simplex.y];
		if (v.simplex.z != -1U) v.simplex.z = from_trimmed[v.simplex.z];
		v = ak::EmbeddedVertex::canonicalize(v.simplex, v.weights);
	}

}
