#include "embed_DAG.hpp"

#include <set>
#include <map>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <string>

bool embed_DAG(
	std::vector< DAGNode > const &nodes,
	std::vector< DAGEdge > const &edges,
	std::vector< uint32_t > *node_options,
	std::vector< int32_t > *node_positions, //positions give total left-to-right order of edges/nodes
	std::vector< int32_t > *edge_positions
) {

	for (auto const &node : nodes) {
		if (node.options.empty()) {
			std::cerr << "WARNING: embed_DAG will fail because a node has no options." << std::endl;
			return false;
		}
	}
	
	{ //PARANOIA:
		for (auto const &node : nodes) {
			assert(!node.options.empty());
			std::set< DAGEdgeIndex > ins(node.options[0].in_order.begin(), node.options[0].in_order.end());
			std::set< DAGEdgeIndex > outs(node.options[0].out_order.begin(), node.options[0].out_order.end());
			//ins and outs reference node properly:
			for (auto i : ins) {
				assert(i < edges.size());
				assert(edges[i].to == &node - &nodes[0]);
			}
			for (auto o : outs) {
				assert(o < edges.size());
				assert(edges[o].from == &node - &nodes[0]);
			}

			for (auto const &option : node.options) {
				assert(option.in_order.size() == ins.size());
				assert(option.in_shapes.size() == ins.size());
				assert(option.out_order.size() == outs.size());
				assert(option.out_shapes.size() == outs.size());
				//node orders all talk about the same edges:
				for (auto i : option.in_order) assert(ins.count(i));
				for (auto o : option.out_order) assert(outs.count(o));
				//node shapes don't index outside edge cost matrices:
				for (uint32_t i = 0; i < option.in_shapes.size(); ++i) {
					assert(option.in_shapes[i] < edges[option.in_order[i]].to_shapes);
				}
				for (uint32_t o = 0; o < option.out_shapes.size(); ++o) {
					assert(option.out_shapes[o] < edges[option.out_order[o]].from_shapes);
				}
			}
		}

		for (auto const &edge : edges) {
			//edges have properly-sized cost matrix:
			assert(edge.costs.size() == edge.from_shapes * edge.to_shapes);
			//edges reference nodes that reference them:
			{
				assert(edge.from < nodes.size());
				DAGNode const &node = nodes[edge.from];
				assert(!node.options.empty());
				auto f = std::find(node.options[0].out_order.begin(), node.options[0].out_order.end(), &edge - &edges[0]);
				assert(f != node.options[0].out_order.end());
			}
			{
				assert(edge.to < nodes.size());
				DAGNode const &node = nodes[edge.to];
				assert(!node.options.empty());
				auto f = std::find(node.options[0].in_order.begin(), node.options[0].in_order.end(), &edge - &edges[0]);
				assert(f != node.options[0].in_order.end());
			}
			//nodes are sorted:
			assert(edge.from < edge.to);
		}
	} //end PARANOIA

	std::vector< uint32_t > select_order;

	//TODO: ~smart ordering~
	for (uint32_t i = 0; i < nodes.size(); ++i) {
		select_order.emplace_back(i);
	}
	std::stable_sort(select_order.begin(), select_order.end(), [&nodes](uint32_t a, uint32_t b) -> bool{
		return nodes[a].options.size() < nodes[b].options.size();
	});

	struct State {
		std::vector< uint32_t > selected; //selected options for each node
		//derived:
		uint32_t step = 0;
		std::vector< bool > left_of; //partial order on the edges
		std::vector< uint32_t > from_shape, to_shape; //selected shapes on edges

		bool operator==(State const &o) const {
			return selected == o.selected;
		}
	};

	struct HashState {
		size_t operator()(State const &state) const {
			static std::hash< std::string > hash;
			return hash(std::string(
				reinterpret_cast< const char * >(&state.selected[0]),
				state.selected.size() * sizeof(uint32_t)));
		}
	};
	

	std::vector< std::pair< DAGCost, const State * > > to_expand;
	std::unordered_map< State, DAGCost, HashState > visited;
	//std::unordered_set< State, HashState > expanded;

	const uint32_t Initial = 8000000;
	to_expand.reserve(Initial);
	visited.reserve(Initial);
	//expanded.reserve(Initial);
	
	std::greater< std::pair< DAGCost, const State * > > CompareCost;

	auto queue_state = [&to_expand, &visited, &CompareCost](State const &state, DAGCost const &cost) {
		auto ret = visited.insert(std::make_pair(state, DAGCost::max()));
		if (cost < ret.first->second) {
			ret.first->second = cost;
			to_expand.emplace_back(cost, &(ret.first->first));
			std::push_heap(to_expand.begin(), to_expand.end(), CompareCost);
		}
	};


	//Any edge that crosses a node must be left_of *all* in/out edges of that node:
	std::vector< std::vector< std::vector< DAGEdgeIndex > > > edge_excludes;
	edge_excludes.reserve(edges.size());
	for (auto const &edge : edges) {
		edge_excludes.emplace_back();
		std::vector< std::vector< DAGEdgeIndex > > &excl = edge_excludes.back();
		for (uint32_t n = edge.from + 1; n + 1 < edge.to; ++n) {
			DAGNode const &node = nodes[n];
			excl.emplace_back();
			excl.back().insert(excl.back().begin(), node.options[0].in_order.begin(), node.options[0].in_order.end());
			excl.back().insert(excl.back().begin(), node.options[0].out_order.begin(), node.options[0].out_order.end());
			std::sort(excl.back().begin(), excl.back().end());
		}
	}

	std::function< bool(std::vector< bool > &, uint32_t, uint32_t) > add_left_of;

	auto run_excludes = [&add_left_of,&edges,&edge_excludes](std::vector< bool > &left_of, uint32_t a) -> bool {
		for (std::vector< uint32_t > const &excl : edge_excludes[a]) {
			bool is_left = false;
			bool is_right = false;
			for (uint32_t b : excl) {
				if (left_of[a * edges.size() + b]) is_left = true;
				if (left_of[b * edges.size() + a]) is_right = true;
			}
			if (is_left && is_right) return false;
			if (is_left) {
				for (uint32_t b : excl) {
					if (!add_left_of(left_of, a, b)) return false;
				}
			}
			if (is_right) {
				for (uint32_t b : excl) {
					if (!add_left_of(left_of, b, a)) return false;
				}
			}
		}
		return true;
	};

	add_left_of = [&edges,&add_left_of,&run_excludes](std::vector< bool > &left_of, uint32_t a, uint32_t b) -> bool {
		if (left_of[a * edges.size() + b]) return true;
		if (left_of[b * edges.size() + a]) return false;
		left_of[a * edges.size() + b] = true;
		for (uint32_t c = 0; c < edges.size(); ++c) {
			if (left_of[c * edges.size() + a] && !add_left_of(left_of, c, b)) return false;
			if (left_of[b * edges.size() + c] && !add_left_of(left_of, a, c)) return false;
		}
		if (!run_excludes(left_of, a)) return false;
		if (!run_excludes(left_of, b)) return false;
		return true;
	};

	auto expand_state = [&queue_state, &select_order, &add_left_of, &nodes, &edges](State const &state, DAGCost const &cost) {
		assert(state.step < select_order.size());
		uint32_t n = select_order[state.step];
		DAGNode const &node = nodes[n];
		for (auto const &option : node.options) {
			State next = state;
			next.step += 1;
			next.selected[n] = &option - &node.options[0];

			bool bad = false;
			for (uint32_t i = 1; i < option.in_order.size(); ++i) {
				if (!add_left_of(next.left_of, option.in_order[i-1], option.in_order[i])) {
					bad = true;
					break;
				}
			}
			if (bad) continue;
			for (uint32_t o = 1; o < option.out_order.size(); ++o) {
				if (!add_left_of(next.left_of, option.out_order[o-1], option.out_order[o])) {
					bad = true;
					break;
				}
			}
			if (bad) continue;

			DAGCost next_cost = cost;
			next_cost += option.cost;
			for (uint32_t i = 0; i < option.in_order.size(); ++i) {
				uint32_t e = option.in_order[i];
				assert(next.to_shape[e] == -1U);
				next.to_shape[e] = option.in_shapes[i];
				if (next.from_shape[e] != -1U) {
					//add cost from edge if edge got completed:
					auto const &edge_cost = edges[e].costs[next.from_shape[e] * edges[e].to_shapes + next.to_shape[e]];
					if (edge_cost == DAGCost::max()) {
						bad = true;
						break;
					}
					next_cost += edge_cost;
				}
			}
			if (bad) continue;

			for (uint32_t o = 0; o < option.out_order.size(); ++o) {
				uint32_t e = option.out_order[o];
				assert(next.from_shape[e] == -1U);
				next.from_shape[e] = option.out_shapes[o];
				if (next.to_shape[e] != -1U) {
					//add cost from edge if edge got completed:
					auto const &edge_cost = edges[e].costs[next.from_shape[e] * edges[e].to_shapes + next.to_shape[e]];
					if (edge_cost == DAGCost::max()) {
						bad = true;
						break;
					}
					next_cost += edge_cost;
				}
			}
			if (bad) continue;


			queue_state(next, next_cost);

		}
	};

	auto set_output = [&](State const &state) {
		assert(state.step == select_order.size());
		assert(state.selected.size() == nodes.size());
		if (node_options) {
			*node_options = state.selected;
		}
		//come up with a total ordering:
		std::vector< bool > left_of = state.left_of;
		for (uint32_t a = 0; a < edges.size(); ++a) {
			for (uint32_t b = 0; b < edges.size(); ++b) {
				//break ordering ties in an arbitrary way:
				if (a != b && (!left_of[a * edges.size() + b] && !left_of[b * edges.size() + a])) {
					bool ret = add_left_of(left_of, a, b);
					assert(ret);
				}
			}
		}
		//store edge and node positions (easy to compute from total order):
		if (node_positions) {
			node_positions->assign(nodes.size(), std::numeric_limits< int32_t >::max());
			node_positions->reserve(nodes.size());
		}
		if (edge_positions) {
			edge_positions->clear();
			edge_positions->reserve(edges.size());
		}
		for (uint32_t a = 0; a < edges.size(); ++a) {
			int32_t pos = edges.size();
			for (uint32_t b = 0; b < edges.size(); ++b) {
				if (left_of[a * edges.size() + b]) {
					pos -= 1;
				}
			}
			if (edge_positions) {
				edge_positions->emplace_back(pos);
			}
			if (node_positions) {
				if (edges[a].from != -1U) (*node_positions)[edges[a].from] = std::min((*node_positions)[edges[a].from], pos);
				if (edges[a].to != -1U) (*node_positions)[edges[a].to] = std::min((*node_positions)[edges[a].to], pos);
			}
		}

	};

	{
		State start;
		start.selected.resize(nodes.size(), -1U);
		start.step = 0;
		start.left_of.resize(edges.size() * edges.size(), false);
		start.from_shape.resize(edges.size(), -1U);
		start.to_shape.resize(edges.size(), -1U);
		queue_state(start, DAGCost::zero());
	}

	uint32_t step = 0;

	while (!to_expand.empty()) {
		std::pop_heap(to_expand.begin(), to_expand.end(), CompareCost);
		DAGCost cost = to_expand.back().first;
		const State &state = *to_expand.back().second;
		to_expand.pop_back();
		auto f = visited.find(state);
		assert(f != visited.end());
		assert(!(cost < f->second));
		if (cost == f->second) {
			//auto res = expanded.insert(state);
			//assert(res.second);
			if ((++step) % 10000 == 0) {
				//DEBUG:
				std::cout << /*expanded.size() << "/" <<*/ visited.size() << "/" << to_expand.size() << "   ";
				std::cout << "[" << state.step << "]";
				for (auto s : state.selected) std::cout << ' ' << (s == -1U ? std::string(".") : std::to_string(s));
				std::cout << std::endl;
			}

			if (state.step < select_order.size()) {
				expand_state(state, cost);
			} else {
				//Found the cheapest selected state!
				set_output(state);

				return true;
			}
		}
	}


	return false;
}
