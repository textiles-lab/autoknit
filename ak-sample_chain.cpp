#include "pipeline.hpp"

void ak::sample_chain(
	ak::Model const &model, //in: model over which chain is defined
	std::vector< ak::EmbeddedVertex > const &chain, //in: chain to be sampled
	std::vector< ak::EmbeddedVertex > *sampled_chain_, //out: sub-sampled chain
	std::vector< ak::Flag > *sampled_flags_ //out: flags (possibly with linkNone on points needed in chain for consistency but not sampled?)
) {
	assert(sampled_chain_);
	auto &sampled_chain = *sampled_chain_;
	sampled_chain.clear();

	assert(sampled_flags_);
	auto &sampled_flags = *sampled_flags_;
	sampled_flags.clear();

	//TODO: do something with nice (jittered?) uniform spacing

	sampled_chain = chain;
	sampled_flags.assign(sampled_chain.size(), FlagLinkAny);
}
