#include "pipeline.hpp"

#include <unordered_map>
#include <iostream>

#include <glm/gtx/hash.hpp>

void ak::find_first_active_chains(
	ak::Parameters const &parameters,
	ak::Model const &model,
	std::vector< float > const &times,
	std::vector< std::vector< ak::EmbeddedVertex > > *active_chains_,
	std::vector< std::vector< Stitch > > *active_stitches_,
	ak::RowColGraph *graph_
) {

	assert(active_chains_);
	auto &active_chains = *active_chains_;
	active_chains.clear();

	assert(active_stitches_);
	auto &active_stitches = *active_stitches_;
	active_stitches.clear();

	//PARANOIA: triangles must reference valid time values:
	for (glm::uvec3 const &tri : model.triangles) {
		assert(tri.x < times.size());
		assert(tri.y < times.size());
		assert(tri.z < times.size());
	}
	//end PARANOIA

	//find boundary loops:

	//build half-edge -> other vertex map:
	std::unordered_map< glm::uvec2, uint32_t > next_vertex;
	{
		auto do_edge = [&next_vertex](uint32_t a, uint32_t b, uint32_t c) {
			assert(a != b && a != c && b != c);
			auto ret = next_vertex.insert(std::make_pair(glm::uvec2(a,b), c));
			if (!ret.second) {
				throw std::runtime_error("ERROR: non-manifold mesh [or inconsistent orientation] -- directed edge appears twice.");
			}
		};
		for (glm::uvec3 const &tri : model.triangles) {
			do_edge(tri.x, tri.y, tri.z);
			do_edge(tri.y, tri.z, tri.x);
			do_edge(tri.z, tri.x, tri.y);
		}
	}

	//extract boundary (== unpaired half-edges):
	std::unordered_map< uint32_t, uint32_t > boundary;
	{
		for (const auto &ev : next_vertex) {
			if (next_vertex.count(glm::uvec2(ev.first.y, ev.first.x))) continue; //skip paired edges
			auto ret = boundary.insert(std::make_pair(ev.first.x, ev.first.y)); //half-edge is on the boundary
			if (!ret.second) {
				throw std::runtime_error("ERROR: non-manifold mesh [or inconsistent orientation] -- vertex is source of multiple boundary edges.");
			}
		}
	}

	//peel chains (well, loops, actually) from boundary:
	while (!boundary.empty()) {
		std::vector< uint32_t > chain;
		chain.emplace_back(boundary.begin()->first);
		chain.emplace_back(boundary.begin()->second);
		boundary.erase(boundary.begin());
		do {
			auto f = boundary.find(chain.back());
			assert(f != boundary.end());
			chain.emplace_back(f->second);
			boundary.erase(f);
		} while (chain.back() != chain[0]);

		//to be marked active, chain must be constant value and < adjacent time values.
		float chain_min = std::numeric_limits< float >::infinity();
		float chain_max = -std::numeric_limits< float >::infinity();
		float adj_min = std::numeric_limits< float >::infinity();
		float adj_max = -std::numeric_limits< float >::infinity();
		for (uint32_t i = 0; i + 1 < chain.size(); ++i) {
			chain_min = std::min(chain_min, times[chain[i]]);
			chain_max = std::max(chain_max, times[chain[i]]);
			auto f = next_vertex.find(glm::uvec2(chain[i], chain[i+1]));
			assert(f != next_vertex.end());
			adj_min = std::min(adj_min, times[f->second]);
			adj_max = std::max(adj_max, times[f->second]);
		}

		//std::cout << "Considering chain with value range [" << chain_min << ", " << chain_max << "] and neighbor value range [" << adj_min << ", " << adj_max << "]." << std::endl;

		if (chain_min != chain_max) {
			std::cerr << "WARNING: discarding chain with non-constant value range [" << chain_min << ", " << chain_max << "]." << std::endl;
			continue;
		}
		//the 1e-3 is to add some slop in case of somewhat noisy interpolation
		if (!(adj_min > chain_max - 1e-3)) {
			if (adj_max < chain_min) {
				//this is a maximum chain
			} else { //this is a mixed min/max chain; weird
				std::cerr << "WARNING: discarding chain with value range [" << chain_min << ", " << chain_max << "] because neighbors have value range [" << adj_min << ", " << adj_max << "]." << std::endl;
			}
			continue;
		}

		//sample chain to get output active chain:
		std::vector< EmbeddedVertex > embedded_chain;
		embedded_chain.reserve(chain.size());

		for (auto c : chain) {
			embedded_chain.emplace_back(EmbeddedVertex::on_vertex(c));
		}

		//subdivide the chain to make sure there are enough samples for later per-sample computations:
		std::vector< EmbeddedVertex > divided_chain;
		sample_chain(parameters.get_chain_sample_spacing(), model, embedded_chain, &divided_chain);

		//further subdivide and place stitches:
		float total_length = 0.0f;
		for (uint32_t ci = 1; ci < divided_chain.size(); ++ci) {
			glm::vec3 a = divided_chain[ci-1].interpolate(model.vertices);
			glm::vec3 b = divided_chain[ci].interpolate(model.vertices);
			total_length += glm::length(b-a);
		}

		float stitch_width = parameters.stitch_width_mm / parameters.model_units_mm;
		uint32_t stitches = std::max(3, int32_t(std::round(total_length / stitch_width)));

		active_chains.emplace_back(divided_chain);
		active_stitches.emplace_back();
		active_stitches.back().reserve(stitches);
		for (uint32_t s = 0; s < stitches; ++s) {
			active_stitches.back().emplace_back((s + 0.5f) / float(stitches), Stitch::FlagLinkAny);
		}
	}

	assert(active_chains.size() == active_stitches.size());

	std::cout << "Found " << active_chains.size() << " first active chains." << std::endl;

	if (graph_) {
		for (uint32_t ci = 0; ci < active_chains.size(); ++ci) {
			auto const &chain = active_chains[ci];

			std::vector< float > lengths;
			lengths.reserve(chain.size());
			lengths.emplace_back(0.0f);
			for (uint32_t i = 1; i < chain.size(); ++i) {
				glm::vec3 a = chain[i-1].interpolate(model.vertices);
				glm::vec3 b = chain[i].interpolate(model.vertices);
				lengths.emplace_back(lengths.back() + glm::length(b-a));
			}
			assert(lengths.size() == chain.size());

			auto li = lengths.begin();
			for (auto &s : active_stitches[ci]) {
				float l = lengths.back() * s.t;

				while (li != lengths.end() && *li <= l) ++li;
				assert(li != lengths.begin());
				assert(li != lengths.end());

				float m = (l - *(li-1)) / (*li - *(li -1));
				uint32_t i = li - lengths.begin();

				assert(s.vertex == -1U);
				s.vertex = graph_->vertices.size();
				graph_->vertices.emplace_back();
				graph_->vertices.back().at = ak::EmbeddedVertex::mix(
					chain[i-1], chain[i], m
				);
			}

			uint32_t prev = (chain[0] == chain.back() ? active_stitches[ci].back().vertex : -1U);
			for (auto &s : active_stitches[ci]) {
				if (prev != -1U) {
					assert(prev < graph_->vertices.size());
					assert(s.vertex < graph_->vertices.size());
					assert(graph_->vertices[prev].row_out == -1U);
					graph_->vertices[prev].row_out = s.vertex;
					assert(graph_->vertices[s.vertex].row_in == -1U);
					graph_->vertices[s.vertex].row_in = prev;
				}
				prev = s.vertex;
			}
		}
	}

}
