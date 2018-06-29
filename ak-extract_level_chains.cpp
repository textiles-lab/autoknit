#include "pipeline.hpp"

#include <unordered_map>
#include <iostream>
#include <deque>

#include <glm/gtx/hash.hpp>

void ak::extract_level_chains(
	ak::Model const &model, //in: model on which to embed vertices
	std::vector< float > const &values, //in: values at vertices
	float const level, //in: level at which to extract chains
	std::vector< std::vector< ak::EmbeddedVertex > > *chains_ //chains of edges at given level
) {
	assert(chains_);
	auto &chains = *chains_;
	chains.clear();

	//embed points along all edges that start below level and end at or above it:

	std::vector< glm::vec3 > const &verts = model.vertices;
	std::vector< glm::uvec3 > const &tris = model.triangles;

	std::unordered_map< glm::uvec2, EmbeddedVertex > embedded_pts;
	std::unordered_map< glm::uvec2, glm::vec3 > pts;
	auto add = [&](uint32_t a, uint32_t b) {
		assert(values[a] < level && values[b] >= level);
		float mix = (level - values[a]) / (values[b] - values[a]);
		pts[glm::uvec2(a,b)] = glm::mix(verts[a], verts[b], mix);
		embedded_pts[glm::uvec2(a,b)] = EmbeddedVertex::on_edge(a,b,mix);
		return glm::uvec2(a,b);
	};
	std::unordered_map< glm::uvec2, glm::uvec2 > links;
	std::unordered_map< glm::uvec2, glm::uvec2 > back_links;
	auto link = [&links,&back_links](glm::uvec2 f, glm::uvec2 t) {
		auto res = links.insert(std::make_pair(f, t));
		assert(res.second);
		auto res2 = back_links.insert(std::make_pair(t, f));
		assert(res2.second);
	};
	for (auto const &tri : tris) {
		uint32_t a = tri.x;
		uint32_t b = tri.y;
		uint32_t c = tri.z;
		//spin triangle until 'a' is the minimum distance value:
		for (uint32_t i = 0; i < 3; ++i) {
			if (values[a] <= values[b] && values[a] <= values[c]) break;
			uint32_t t = a; a = b; b = c; c = t;
		}
		//NOTE: we treat level as "level + epsilon"
		if (values[a] >= level) continue; //all above border
		assert(values[a] < level);
		//NOTE: if values increase along +y, chains should be oriented in the +x direction
		//assuming ccw oriented triangles, this means:

		if (values[b] >= level && values[c] >= level) {
			//edge is from ca to ab
			link(add(a,c), add(a,b));
		} else if (values[b] >= level && values[c] < level) {
			//edge is from bc to ab
			link(add(c,b), add(a,b));
		} else if (values[b] < level && values[c] >= level) {
			//edge is from ca to bc
			link(add(a,c), add(b,c));
		} else {
			assert(values[b] < level && values[c] < level);
			//all below border, nothing to do.
		}
	}

	uint32_t found_chains = 0;
	uint32_t found_loops = 0;

	//read back path from links:
	while (!links.empty()) {
		std::deque< glm::uvec2 > loop;
		loop.emplace_back(links.begin()->first);
		loop.emplace_back(links.begin()->second);

		//remove seed link:
		links.erase(links.begin());
		{
			auto b = back_links.find(loop.back());
			assert(b != back_links.end());
			assert(b->second == loop[0]);
			back_links.erase(b);
		}

		//extend forward:
		while (true) {
			auto f = links.find(loop.back());
			if (f == links.end()) break;
			loop.emplace_back(f->second);
			//remove link:
			auto b = back_links.find(loop.back());
			assert(b != back_links.end());
			assert(b->second == loop[loop.size()-2]);
			links.erase(f);
			back_links.erase(b);
		}
		//extend backward:
		while (true) {
			auto b = back_links.find(loop[0]);
			if (b == back_links.end()) break;
			loop.emplace_front(b->second);
			//remove link:
			auto f = links.find(loop.front());
			assert(f != links.end());
			assert(f->second == loop[1]);
			back_links.erase(b);
			links.erase(f);
		}

		if (loop.front() == loop.back()) ++found_loops;
		else ++found_chains;

		chains.emplace_back();
		chains.back().reserve(loop.size());
		for (glm::uvec2 e : loop) {
			auto f = embedded_pts.find(e);
			assert(f != embedded_pts.end());
			chains.back().emplace_back(f->second);
		}
	}

	std::cout << "extract_level_chains found " << found_loops << " loops and " << found_chains << " chains." << std::endl;
}
