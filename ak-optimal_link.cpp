#include "pipeline.hpp"

#include <unordered_map>
#include <stdexcept>
#include <cstring>

void ak::optimal_link(
	float target_distance, bool do_roll,
	std::vector< glm::vec3 > const &source,
	std::vector< bool > const &source_linkone,
	std::vector< glm::vec3 > const &target,
	std::vector< bool > const &target_linkone,
	std::vector< std::pair< uint32_t, uint32_t > > *links_) {

	assert(source.size() == source_linkone.size());
	assert(target.size() == target_linkone.size());

	assert(links_);
	auto &links = *links_;
	links.clear();

	//options:
	// s[i]  s[i] s[i+1]   s[i]
	//  |       \ /       /  |  
	// t[i]     t[i]    t[i] t[i+1]

	struct State {
		uint16_t source_idx = 0xffff;
		uint16_t target_idx = 0xffff;
		uint16_t source_remain = 0x0; //unlinked source
		uint16_t target_remain = 0x0; //unlinked target
		uint64_t pack() {
			uint64_t ret;
			memcpy(reinterpret_cast< char * >(&ret), this, sizeof(uint64_t));
			return ret;
		}
		static State unpack(uint64_t packed) {
			State ret;
			memcpy(reinterpret_cast< char * >(&ret), &packed, sizeof(uint64_t));
			return ret;
		}
	};
	static_assert(sizeof(State) == sizeof(uint64_t), "State should be packed");

	struct Action {
		Action(uint8_t take_source_, uint8_t take_target_) : take_source(take_source_), take_target(take_target_) { }
		Action() = default;
		uint8_t take_source = 0;
		uint8_t take_target = 0;
	};

	std::unordered_map< uint64_t, std::pair< float, Action > > distance_via;
	std::vector< std::pair< float, uint64_t > > heap;

	auto visit = [&distance_via, &heap](uint64_t state, float distance, Action const &action) {
		auto f = distance_via.insert(std::make_pair(state, std::make_pair(std::numeric_limits< float >::infinity(), Action()))).first;
		if (distance < f->second.first) {
			f->second = std::make_pair(distance, action);
			heap.emplace_back(-distance, state);
			std::push_heap(heap.begin(), heap.end());
		}
	};

	if (do_roll) {
		for (uint32_t t = 0; t < target.size(); ++t) {
			State s;
			s.source_idx = 0;
			s.target_idx = t;
			s.source_remain = uint32_t(source.size());
			s.target_remain = uint32_t(target.size());
			visit(s.pack(), 0.0f, Action());

			s.source_idx = 1;
			visit(s.pack(), 0.0f, Action());
		}
	} else {
		State s;
		s.source_idx = 0;
		s.target_idx = 0;
		s.source_remain = uint32_t(source.size());
		s.target_remain = uint32_t(target.size());
		visit(s.pack(), 0.0f, Action());
	}

	auto penalty = [&source, &target, &target_distance](uint32_t si, uint32_t ti) {
		assert(si < source.size());
		assert(ti < target.size());
		float dis = glm::length(source[si] - target[ti]) - target_distance;
		return dis*dis;
	};

	//TODO: could refine this a fair bit & build a table with real link counts
	auto is_possible = [](State const &state) {
		if (state.source_remain * 2 < state.target_remain) return false;
		if (state.target_remain * 2 < state.source_remain) return false;
		return true;
	};

	State best;
	while (!heap.empty()) {
		std::pop_heap(heap.begin(), heap.end());
		State state = State::unpack(heap.back().second);
		float distance = -heap.back().first;
		heap.pop_back();
		auto f = distance_via.find(state.pack());
		assert(f != distance_via.end());
		assert(f->second.first <= distance);
		if (f->second.first < distance) continue;
		if (state.source_remain == 0 && state.target_remain == 0) {
			best = state;
			break;
		}
		assert(is_possible(state));
		//try the actions:
		{ // 1-1
			float next_distance = distance + penalty(state.source_idx, state.target_idx);
			State next = state;
			next.source_idx = (next.source_idx + 1) % source.size();
			next.target_idx = (next.target_idx + 1) % target.size();
			next.source_remain -= 1;
			next.target_remain -= 1;
			if (is_possible(next)) visit(next.pack(), next_distance, Action(1,1));
		}
		// 1-2
		if ( state.target_remain >= 2
			&& !source_linkone[state.source_idx]
			&& !target_linkone[state.target_idx]
			&& !target_linkone[(state.target_idx+1) % target.size()] ) {
			float next_distance = distance
				+ penalty(state.source_idx, state.target_idx)
				+ penalty(state.source_idx, (state.target_idx+1) % target.size());
			State next = state;
			next.source_idx = (next.source_idx + 1) % source.size();
			next.target_idx = (next.target_idx + 2) % target.size();
			next.source_remain -= 1;
			next.target_remain -= 2;
			if (is_possible(next)) visit(next.pack(), next_distance, Action(1,2));
		}
		// 2-1
		if ( state.source_remain >= 2
			&& !source_linkone[state.source_idx]
			&& !source_linkone[(state.source_idx+1) % source.size()]
			&& !target_linkone[state.target_idx] ) {
			float next_distance = distance
				+ penalty(state.source_idx, state.target_idx)
				+ penalty((state.source_idx+1)%source.size(), state.target_idx);
			State next = state;
			next.source_idx = (next.source_idx + 2) % source.size();
			next.target_idx = (next.target_idx + 1) % target.size();
			next.source_remain -= 2;
			next.target_remain -= 1;
			if (is_possible(next)) visit(next.pack(), next_distance, Action(2,1));
		}
	}

	if (best.source_idx == 0xffff) {
		//Failed!
		throw std::runtime_error("Failed to link");
	}

	//Read back:
	State at = best;
	while (true) {
		auto f = distance_via.find(at.pack());
		assert(f != distance_via.end());
		Action action = f->second.second;
		if (action.take_source == 0 && action.take_target == 0) break;
		if (action.take_source == 1) {
			at.source_idx = (at.source_idx + uint32_t(source.size()) - 1) % uint32_t(source.size());
			at.source_remain += 1;
			while (action.take_target > 0) {
				--action.take_target;
				at.target_idx = (at.target_idx + uint32_t(target.size()) - 1) % uint32_t(target.size());
				at.target_remain += 1;
				links.emplace_back(at.source_idx, at.target_idx);
			}
		} else if (action.take_target == 1) {
			at.target_idx = (at.target_idx + uint32_t(target.size()) - 1) % uint32_t(target.size());
			at.target_remain += 1;
			while (action.take_source > 0) {
				--action.take_source;
				at.source_idx = (at.source_idx + uint32_t(source.size()) - 1) % uint32_t(source.size());
				at.source_remain += 1;
				links.emplace_back(at.source_idx, at.target_idx);
			}
		} else {
			assert(action.take_source == 1 || action.take_target == 1);
		}
	}
	assert(at.source_remain == source.size() && at.target_remain == target.size());
	assert(at.source_idx == 0 || at.source_idx == 1);

	std::reverse(links.begin(), links.end());
}

