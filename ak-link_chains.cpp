#include "pipeline.hpp"

void ak::link_chains(
	ak::Parameters const &parameters,
	ak::Model const &model, //in: model
	std::vector< float > const &times,          //in: time field (times @ vertices)
	std::vector< std::vector< ak::EmbeddedVertex > > const &active_chains, //in: current active chains
	std::vector< std::vector< ak::Flag > > const &active_flags, //in: stitches
	std::vector< std::vector< ak::EmbeddedVertex > > const &next_chains, //in: next chains
	std::vector< std::vector< ak::EmbeddedVertex > > *linked_next_chains_, //out: next chains
	std::vector< std::vector< ak::Flag > > *linked_next_flags_, //out: flags indicating status of vertices on next chains
	std::vector< ak::Link > *links_ //out: active_chains[from_chain][from_vertex] -> linked_next_chains[to_chain][to_vertex] links
) {


}

