#pragma once

#include "ScheduleCost.hpp"

#include <vector>

typedef ScheduleCost DAGCost;

typedef uint32_t DAGNodeIndex;
typedef uint32_t DAGEdgeIndex;

struct DAGEdge {
	DAGNodeIndex from = -1U;
	DAGNodeIndex to = -1U;
	uint32_t from_shapes = 0;
	uint32_t to_shapes = 0;
	std::vector< DAGCost > costs;
	//If from != -1U && to != -1U: cost[f * to_shapes + t] --> cost given from/to shape
	//If from == -1U: cost[t] is min-cost given to shape
	//If to == -1U: cost[f] is min-cost given from shape
};

struct DAGOption {
	DAGCost cost;
	std::vector< DAGEdgeIndex > in_order;
	std::vector< DAGEdgeIndex > out_order;
	std::vector< uint32_t > in_shapes;
	std::vector< uint32_t > out_shapes;
};

struct DAGNode {
	std::vector< DAGOption > options;
};

//NOTE: edges must have from < to (i.e., nodes should be topologically sorted)
bool embed_DAG(
	std::vector< DAGNode > const &nodes,
	std::vector< DAGEdge > const &edges,
	std::vector< uint32_t > *node_options,
	std::vector< int32_t > *node_positions, //positions give total left-to-right order of edges/nodes
	std::vector< int32_t > *edge_positions
);
