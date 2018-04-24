#include "pipeline.hpp"

#include <set>
#include <algorithm>
#include <iostream>

void ak::embed_constraints(
	ak::Model const &model,
	std::vector< ak::Constraint > const &constraints,
	ak::Model *_constrained_model,
	std::vector< float > *_constrained_values, //same size as out_model's vertices
	std::vector< std::vector< glm::vec3 > > *DEBUG_chain_paths
) {
	assert(_constrained_model);
	auto &constrained_model = *_constrained_model;
	constrained_model = ak::Model();

	assert(_constrained_values);
	auto &constrained_values = *_constrained_values;
	constrained_values = std::vector< float >();

	if (DEBUG_chain_paths) {
		*DEBUG_chain_paths = std::vector< std::vector< glm::vec3 > >();
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
		if (DEBUG_chain_paths) {
			DEBUG_chain_paths->emplace_back();
			for (uint32_t v : path) {
				DEBUG_chain_paths->back().emplace_back(model.vertices[v]);
			}
		}
	}

}
