#include "pipeline.hpp"

#include <glm/gtx/norm.hpp>
#include <glm/gtx/hash.hpp>

#include <iostream>
#include <unordered_map>
#include <unordered_set>


void ak::peel_slice(
	Parameters const &parameters,
	Model const &model,
	std::vector< std::vector< EmbeddedVertex > > const &active_chains,
	Model *slice_,
	std::vector< EmbeddedVertex > *slice_on_model_,
	std::vector< std::vector< uint32_t > > *slice_active_chains_,
	std::vector< std::vector< uint32_t > > *slice_next_chains_,
	std::vector< bool > *used_boundary_
) {
	assert(slice_);
	auto &slice = *slice_;
	slice.clear();

	assert(slice_on_model_);
	auto &slice_on_model = *slice_on_model_;
	slice_on_model.clear();

	assert(slice_active_chains_);
	auto &slice_active_chains = *slice_active_chains_;
	slice_active_chains.clear();

	assert(slice_next_chains_);
	auto &slice_next_chains = *slice_next_chains_;
	slice_next_chains.clear();

	{
		//DEBUG:
		uint32_t loops = 0;
		uint32_t lines = 0;
		for (auto const &chain : active_chains) {
			if (chain[0] == chain.back()) ++loops;
			else ++lines;
		}
		std::cout << "---- peel slice on [" << loops << " loops and " << lines << " lines] ----" << std::endl;
	}

	Model clipped;
	std::vector< ak::EmbeddedVertex > clipped_on_model;
	ak::trim_model(model, active_chains, std::vector< std::vector< ak::EmbeddedVertex > >(), &clipped, &clipped_on_model);

	//This version of the code just uses the 3D distance to the curve.
	//might have problems with models that get really close to themselves.

	std::vector< float > values(clipped.vertices.size(), std::numeric_limits< float >::infinity());

	auto do_seg = [&values,&clipped](glm::vec3 const &a, glm::vec3 const &b) {
		if (a == b) return;
		glm::vec3 ab = b-a;
		float limit = glm::dot(ab,ab);
		float inv_limit = 1.0f / limit;
		for (auto const &v : clipped.vertices) {
			float amt = glm::dot(v-a, ab);
			amt = std::max(0.0f, std::min(limit, amt));
			glm::vec3 pt = (amt * inv_limit) * (b-a) + a;
			float dis2 = glm::length2(v - pt);
			float &best2 = values[&v - &clipped.vertices[0]];
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

	//extract chains:

	std::vector< std::vector< ak::EmbeddedVertex > > next_chains;
	{
		float level = 2.0f * parameters.stitch_height_mm / parameters.model_units_mm;

		std::vector< std::vector< EmbeddedVertex > > level_chains;

		ak::extract_level_chains(clipped, values, level, &level_chains);

		{ //(sort-of) hack: make all chains into loops by including portions of the boundary if needed:
			std::unordered_map< glm::uvec2, uint32_t > next;
			auto do_edge = [&next](uint32_t a, uint32_t b, uint32_t c) {
				auto ret = next.insert(std::make_pair(glm::uvec2(a,b), c));
				assert(ret.second);
			};
			for (auto const &tri : clipped.triangles) {
				do_edge(tri.x, tri.y, tri.z);
				do_edge(tri.y, tri.z, tri.x);
				do_edge(tri.z, tri.x, tri.y);
			}
			struct ChainEnd {
				ChainEnd(float along_, uint32_t chain_, bool is_start_) : along(along_), chain(chain_), is_start(is_start_) { }
				float along;
				uint32_t chain;
				bool is_start;
			};
			std::unordered_map< glm::uvec2 , std::vector< ChainEnd > > on_edge;
			auto do_end = [&](EmbeddedVertex const &ev, uint32_t chain, bool is_start) {
				assert(ev.simplex.x != -1U);
				assert(ev.simplex.y != -1U);
				assert(ev.simplex.z == -1U);

				glm::uvec2 e = glm::uvec2(ev.simplex.x, ev.simplex.y);
				float amt = ev.weights.y;
				if (next.count(e)) {
					e = glm::uvec2(ev.simplex.y, ev.simplex.x);
					amt = ev.weights.x;
				}

				assert(next.count(e) + next.count(glm::uvec2(e.y,e.x)) == 1); //e should be a boundary edge
				assert(next.count(e) == 0);

				on_edge[e].emplace_back(amt, chain, is_start);
			};
			for (auto &chain : level_chains) {
				if (chain[0] == chain.back()) continue; //ignore loops
				do_end(chain[0], &chain - &level_chains[0], true);
				do_end(chain.back(), &chain - &level_chains[0], false);
			}
			std::vector< uint32_t > append(level_chains.size(), -1U);
			std::vector< bool > used_boundary(level_chains.size(), true);

			//loops marked as such before connection-making:
			for (uint32_t c = 0; c < level_chains.size(); ++c) {
				if (level_chains[c][0] == level_chains[c].back()) {
					append[c] = c;
					used_boundary[c] = false;
				}
			}

			auto chase_path = [&](glm::uvec2 begin_e, ChainEnd begin_ce) {
				assert(!begin_ce.is_start); //start at the end of a path, chase to the start of another
				std::vector< EmbeddedVertex > path;
				//path.emplace_back(EmbeddedVertex::on_edge(begin_e, begin_ce.along));

				ChainEnd const *end_ce = nullptr;

				glm::uvec2 e = begin_e;
				ChainEnd ce = begin_ce;
				while (true) {
					{ //find next point along edge, if it exists:
						ChainEnd const *found_ce = nullptr;
						auto f = on_edge.find(e);
						if (f != on_edge.end()) {
							for (auto const &nce : f->second) {
								if (nce.along <= ce.along) continue;
								if (found_ce == nullptr || nce.along < found_ce->along) {
									found_ce = &nce;
								}
							}
						}
						if (found_ce) {
							assert(found_ce->is_start);
							end_ce = found_ce;
							break;
						}
					}
					//next point is end of edge.
					//add vertex:
					path.emplace_back(EmbeddedVertex::on_vertex(e.y));
					//circulate to next edge:
					glm::uvec2 old_e = e;
					while (true) {
						auto f = next.find(glm::uvec2(e.y, e.x));
						if (f == next.end()) break;
						e = glm::uvec2(f->second, e.y);
					}
					assert(e.y == old_e.y);
					assert(e.x != old_e.x);
					e = glm::uvec2(e.y, e.x);
					ce.chain = -1U;
					ce.along = 0.0f;
				}
				assert(end_ce);
				assert(end_ce->is_start); //start of (another?) path.
				path.emplace_back(level_chains[end_ce->chain][0]);
				level_chains[begin_ce.chain].insert(level_chains[begin_ce.chain].end(), path.begin(), path.end());
				//flag that append code should append other chain:
				assert(append[begin_ce.chain] == -1U);
				append[begin_ce.chain] = end_ce->chain;
			};

			for (auto const &seed_ece : on_edge) {
				for (auto const &seed : seed_ece.second) {
					if (!seed.is_start) {
						chase_path(seed_ece.first, seed);
					}
				}
			}

			for (uint32_t c = 0; c < level_chains.size(); ++c) {
				if (append[c] == -1U) continue; //marked for discard
				while (append[c] != c) {
					uint32_t a = append[c];
					assert(a < level_chains.size());
					assert(!level_chains[a].empty());
					assert(level_chains[c].back() == level_chains[a][0]); //already have first vertex
					level_chains[c].insert(level_chains[c].end(), level_chains[a].begin()+1, level_chains[a].end());
					append[c] = append[a];

					append[a] = -1U; //mark for discard
					level_chains[a].clear();
				}
			}
			for (uint32_t c = 0; c < level_chains.size(); /* later */) {
				if (append[c] == -1U) {
					assert(level_chains[c].empty());
					std::swap(used_boundary[c], used_boundary.back());
					used_boundary.pop_back();
					std::swap(level_chains[c], level_chains.back());
					level_chains.pop_back();
				} else {
					assert(!level_chains[c].empty());
					assert(level_chains[c][0] == level_chains[c].back());
					++c;
				}
			}

			if (used_boundary_) *used_boundary_ = used_boundary;

		}

		uint32_t loops = 0;
		uint32_t lines = 0;

		next_chains.reserve(level_chains.size());
		for (auto &chain : level_chains) {
			//chain is embedded on 'clipped' which is embedded on 'model'; re-embed on just 'model':
			for (auto &v : chain) {
				glm::uvec3 simplex = clipped_on_model[v.simplex.x].simplex;
				if (v.simplex.y != -1U) simplex = ak::EmbeddedVertex::common_simplex(simplex, clipped_on_model[v.simplex.y].simplex);
				if (v.simplex.z != -1U) simplex = ak::EmbeddedVertex::common_simplex(simplex, clipped_on_model[v.simplex.z].simplex);
				glm::vec3 weights = v.weights.x * clipped_on_model[v.simplex.x].weights_on(simplex);
				if (v.simplex.y != -1U) weights += v.weights.y * clipped_on_model[v.simplex.y].weights_on(simplex);
				if (v.simplex.z != -1U) weights += v.weights.z * clipped_on_model[v.simplex.z].weights_on(simplex);
				v.simplex = simplex;
				v.weights = weights;
			}

			//subdivide chain and add to outputs:
			if (chain[0] == chain.back()) ++loops;
			else ++lines;
			next_chains.emplace_back();
			sample_chain(parameters.get_chain_sample_spacing(), model, chain, &next_chains.back());
		}
		std::cout << "  extracted " << loops << " loops and " << lines << " lines." << std::endl;
	}

	//PARANOIA:
	for (auto const &chain : active_chains) {
		for (auto const &v : chain) {
			assert(v.simplex.x < model.vertices.size());
			assert(v.simplex.y == -1U || v.simplex.y < model.vertices.size());
			assert(v.simplex.z == -1U || v.simplex.z < model.vertices.size());
		}
		for (uint32_t i = 1; i < chain.size(); ++i) {
			assert(chain[i-1] != chain[i]);
		}
	}

	for (auto const &chain : next_chains) {
		for (auto const &v : chain) {
			assert(v.simplex.x < model.vertices.size());
			assert(v.simplex.y == -1U || v.simplex.y < model.vertices.size());
			assert(v.simplex.z == -1U || v.simplex.z < model.vertices.size());
		}
		for (uint32_t i = 1; i < chain.size(); ++i) {
			assert(chain[i-1] != chain[i]);
		}
	}
	//end PARANOIA

	//now actually pull out the proper slice:
	ak::trim_model(model, active_chains, next_chains, &slice, &slice_on_model, &slice_active_chains, &slice_next_chains);

	//sometimes this can combine vertices, in which case the output chains should be trimmed:
	uint32_t trimmed = 0;
	for (auto &chain : slice_active_chains) {
		for (auto v : chain) {
			assert(v < slice.vertices.size());
		}
		for (uint32_t i = 1; i < chain.size(); /* later */) {
			if (chain[i-1] == chain[i]) {
				chain.erase(chain.begin() + i);
				++trimmed;
			} else {
				++i;
			}
		}
	}

	for (auto &chain : slice_next_chains) {
		for (auto v : chain) {
			assert(v < slice.vertices.size());
		}
		for (uint32_t i = 1; i < chain.size(); ++i) {
			if (chain[i-1] == chain[i]) {
				chain.erase(chain.begin() + i);
				++trimmed;
			} else {
				++i;
			}
		}
	}

	if (trimmed) {
		std::cout << "Trimmed " << trimmed << " too-close-for-epm vertices from slice chains." << std::endl;
	}
}

