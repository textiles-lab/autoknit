#include "pipeline.hpp"

void ak::sample_chain(
	float spacing,
	ak::Model const &model, //in: model over which chain is defined
	std::vector< ak::EmbeddedVertex > const &chain, //in: chain to be sampled
	std::vector< ak::EmbeddedVertex > *sampled_chain_ //out: sub-sampled chain
	//std::vector< ak::Flag > *sampled_flags_ //out: flags (possibly with linkNone on points needed in chain for consistency but not sampled?)
) {
	assert(sampled_chain_);
	auto &sampled_chain = *sampled_chain_;
	sampled_chain.clear();

	//assert(sampled_flags_);
	//auto &sampled_flags = *sampled_flags_;
	//sampled_flags.clear();

	for (uint32_t ci = 0; ci + 1 < chain.size(); ++ci) {
		sampled_chain.emplace_back(chain[ci]);

		glm::vec3 a = chain[ci].interpolate(model.vertices);
		glm::vec3 b = chain[ci+1].interpolate(model.vertices);
		glm::uvec3 common = EmbeddedVertex::common_simplex(chain[ci].simplex, chain[ci+1].simplex);
		glm::vec3 wa = chain[ci].weights_on(common);
		glm::vec3 wb = chain[ci+1].weights_on(common);

		float length = glm::length(b-a);
		int32_t insert = std::floor(length / spacing);
		for (int32_t i = 0; i < insert; ++i) {
			float m = float(i+1) / float(insert+1);
			sampled_chain.emplace_back(common, glm::mix(wa, wb, m));
		}
	}
	sampled_chain.emplace_back(chain.back());

	//sampled_flags.assign(sampled_chain.size(), FlagLinkAny);
}
