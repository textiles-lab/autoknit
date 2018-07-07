#include "pipeline.hpp"

#include "EmbeddedPlanarMap.hpp"

#include <glm/gtx/hash.hpp>

#include <iostream>

void ak::trim_model(
	ak::Model const &model, //in: model
	std::vector< std::vector< ak::EmbeddedVertex > > const &left_of,
	std::vector< std::vector< ak::EmbeddedVertex > > const &right_of,
	ak::Model *clipped_, //out: portion of model's surface that is left_of the left_of chains and right_of the right_of chains
	std::vector< ak::EmbeddedVertex > *clipped_vertices_ //out: map from clipped vertices to source mesh
) {
	assert(clipped_);
	auto &clipped = *clipped_;
	clipped.clear();

	assert(clipped_vertices_);
	auto &clipped_vertices = *clipped_vertices_;
	clipped_vertices.clear();

	//embed chains using planar map:
	EmbeddedPlanarMap< int32_t, NegativeValue< int32_t >, SumValues< int32_t > > epm;
	uint32_t total_chain_edges = 0;
	for (auto const &chain : left_of) {
		for (uint32_t i = 0; i + 1 < chain.size(); ++i) {
			epm.add_edge(chain[i], chain[i+1], 1);
			++total_chain_edges;
		}
	}
	for (auto const &chain : right_of) {
		for (uint32_t i = 0; i + 1 < chain.size(); ++i) {
			epm.add_edge(chain[i+1], chain[i], (1 << 8) );
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
			int32_t value = ee.value;
			auto ret = edge_values.insert(std::make_pair(glm::uvec2(a,b), value));
			assert(ret.second);

			epm.reverse_value(&value);
			ret = edge_values.insert(std::make_pair(glm::uvec2(b,a), value));
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

}

