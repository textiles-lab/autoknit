#include "pipeline.hpp"

#include <glm/gtx/norm.hpp>

#include <iostream>

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

	//for (int step = 1; step < 10; ++step) { //DEBUG: show several distances
	{	int step = 1;

		float level = step * (parameters.stitch_height_mm / parameters.model_units_mm);

		std::vector< std::vector< EmbeddedVertex > > temp_chains;

		ak::extract_level_chains(model, values, level, &temp_chains);
		
		next_chains.insert(next_chains.end(), temp_chains.begin(), temp_chains.end());
	}

	//Trim chains that are not immediately left-of the active chains (should deal with problems of the surface looping back on itself?)

	//IDEA: start at active chain vertices and move to left-of- 
	//std::unordered_multimap< glm::uvec2, std::pair< uint32_t, uint32_t > > edge_vertices;
	


}

