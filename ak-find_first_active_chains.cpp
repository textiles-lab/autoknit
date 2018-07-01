#include "pipeline.hpp"

#include <unordered_map>
#include <iostream>

#include <glm/gtx/hash.hpp>

void ak::find_first_active_chains(
	ak::Parameters const &parameters,
	ak::Model const &model, //in: model (vertices & triangles)
	std::vector< float > const &times,          //in: time field (times @ vertices)
	std::vector< std::vector< ak::EmbeddedVertex > > *active_chains_, //out: all mesh boundaries that contain a minimum
	std::vector< std::vector< ak::Flag > > *active_flags_ //out: all mesh boundaries that contain a minimum
) {

	assert(active_chains_);
	auto &active_chains = *active_chains_;
	active_chains.clear();

	assert(active_flags_);
	auto &active_flags = *active_flags_;
	active_flags.clear();

	//PARANOIA: triangles must reference valid time values:
	for (glm::uvec3 const &tri : model.triangles) {
		assert(tri.x < times.size());
		assert(tri.y < times.size());
		assert(tri.z < times.size());
	}

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
		std::vector< float > lengths;
		lengths.reserve(divided_chain.size()-1);
		float total_length = 0.0f;
		for (uint32_t ci = 0; ci + 1 < divided_chain.size(); ++ci) {
			glm::vec3 a = divided_chain[ci].interpolate(model.vertices);
			glm::vec3 b = divided_chain[ci+1].interpolate(model.vertices);
			float l = glm::length(b-a);
			lengths.emplace_back(l);
			total_length += l;
		}
		assert(lengths.size() == divided_chain.size()-1);

		float stitch_width = parameters.stitch_width_mm / parameters.model_units_mm;
		uint32_t stitches = std::max(3, int32_t(std::round(total_length / stitch_width)));

		active_chains.emplace_back();
		active_flags.emplace_back();

		{ //subdivide chain further
			float stitch_step = total_length / stitches;
			float stitch_acc = 0.5f * stitch_step;
			for (uint32_t di = 0; di + 1 < divided_chain.size(); ++di) {

				active_chains.back().emplace_back(divided_chain[di]);
				//stitches close enough to an existing vertex get snapped:
				if (stitch_acc < 0.01f * stitch_width) {
					active_flags.back().emplace_back(ak::FlagLinkAny);
					stitch_acc += stitch_step;
				} else {
					active_flags.back().emplace_back(ak::FlagLinkNone);
				}

				if (di + 1 == divided_chain.size()) break;


				float length = lengths[di];

				float remain = length;
				while (stitch_acc < remain - 0.01f * stitch_width) {
					remain -= stitch_acc;
					active_chains.back().emplace_back(EmbeddedVertex::mix(
						divided_chain[di], divided_chain[di+1],
						(length - remain) / length
					));
					active_flags.back().emplace_back(ak::FlagLinkAny);
					stitch_acc = stitch_step;
				}
				stitch_acc -= remain;
			}
			active_chains.back().emplace_back(divided_chain.back());

			assert(active_flags.back().size() == active_chains.back().size());

			assert((divided_chain[0] == divided_chain.back()) == (active_chains.back()[0] == active_chains.back().back()));

			//make sure this worked:
			uint32_t any_count = 0;
			for (auto f : active_flags.back()) {
				if (f == ak::FlagLinkAny) ++any_count;
			}
			assert(any_count == stitches);
		}
	}

	assert(active_chains.size() == active_flags.size());

	std::cout << "Found " << active_chains.size() << " first active chains." << std::endl;

}
