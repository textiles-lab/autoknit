#pragma once

#include "pipeline.hpp"

#include <glm/gtx/hash.hpp>

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <iostream>
#include <algorithm>

#define WEIGHT_SUM  (1024U)
struct IntegerEmbeddedVertex {
	glm::uvec3 simplex;
	glm::ivec3 weights;

	IntegerEmbeddedVertex(const glm::uvec3 &simplex_, const glm::ivec3 &weights_) : simplex(simplex_), weights(weights_) {
		assert(weights.x + weights.y + weights.z == WEIGHT_SUM);
	}

	static IntegerEmbeddedVertex simplify(const IntegerEmbeddedVertex &in) {
		assert(in.simplex.x < in.simplex.y && in.simplex.y <= in.simplex.z);
		assert(in.weights.x + in.weights.y + in.weights.z == WEIGHT_SUM);
		IntegerEmbeddedVertex ret(glm::uvec3(0, -1U, -1U), glm::ivec3(1024, 0,0));
		uint32_t o = 0;
		for (uint32_t i = 0; i < 3; ++i) {
			if (in.weights[i] != 0) {
				ret.simplex[o] = in.simplex[i];
				ret.weights[o] = in.weights[i];
				++o;
			}
		}
		return ret;
	}

	IntegerEmbeddedVertex(const ak::EmbeddedVertex &src) : simplex(src.simplex) {
		assert(simplex.x != -1U);
		assert(simplex.x < simplex.y);
		assert((simplex.y == -1U && simplex.z == -1U) || simplex.y < simplex.z);

		float x = src.weights.x;
		float xy = x + src.weights.y;
		float xyz = xy + src.weights.z;

		int32_t ix = std::round(WEIGHT_SUM * (x / xyz));
		int32_t ixy = std::round(WEIGHT_SUM * (xy / xyz));
		int32_t ixyz = std::round(WEIGHT_SUM * (xyz / xyz));

		weights = glm::ivec3(ix, ixy - ix, ixyz - ixy);
		assert(weights.x + weights.y + weights.z == WEIGHT_SUM);
	};

	bool operator==(const IntegerEmbeddedVertex &o) const {
		return simplex == o.simplex && weights == o.weights;
	}

	template< typename T >
	T interpolate(std::vector< T > const &values) const {
		T ret = values[simplex.x] * (float(weights.x) / float(WEIGHT_SUM));
		if (simplex.y != -1U) ret += values[simplex.y] * (float(weights.y) / float(WEIGHT_SUM));
		if (simplex.z != -1U) ret += values[simplex.z] * (float(weights.z) / float(WEIGHT_SUM));
		return ret;
	}

	glm::ivec3 weights_on(glm::uvec3 simplex2) const {
		glm::ivec3 ret(0,0,0);
		uint32_t o = 0;
		for (uint32_t i = 0; i < 3; ++i) {
			if (simplex[i] == -1U) break;
			while (simplex2[o] < simplex[i]) {
				++o;
				assert(o < 3);
			}
			assert(simplex2[o] == simplex[i]);
			ret[o] = weights[i];
		}
		assert(ret.x + ret.y + ret.z == WEIGHT_SUM);

		return ret;
	}

	static glm::uvec3 common_simplex(const glm::uvec3 &a, const glm::uvec3 &b) {
		return ak::EmbeddedVertex::common_simplex(a,b);
	}


};


template< typename VALUE >
struct SameValue {
	static void reverse(VALUE *v) { }
};

template< typename VALUE >
struct NegativeValue {
	static void reverse(VALUE *v) { *v = -*v; }
};

template< typename VALUE >
struct ReplaceValue {
	static void combine(VALUE *tgt, VALUE const &src) { *tgt = src; }
};

template< typename VALUE >
struct SumValues {
	static void combine(VALUE *tgt, VALUE const &src) { *tgt += src; }
};

template< typename VALUE >
struct CopyValue {
	static void split(VALUE const &src, VALUE *first, VALUE *second) { *first = *second = src; }
};


template< typename VALUE >
struct EmbeddedEdge {
	EmbeddedEdge(uint32_t first_, uint32_t second_, VALUE const &value_) : first(first_), second(second_), value(value_) { }
	uint32_t first;
	uint32_t second;
	VALUE value;
};

template< typename VALUE, class REVERSE_VALUE = SameValue< VALUE >, class COMBINE_VALUES = ReplaceValue< VALUE >, class SPLIT_VALUE = CopyValue< VALUE > >
struct EmbeddedPlanarMap {
	std::vector< IntegerEmbeddedVertex > vertices;

	std::unordered_map< glm::uvec3, std::vector< uint32_t > > simplex_vertices;
	std::unordered_map< glm::uvec3, std::vector< EmbeddedEdge< VALUE > > > simplex_edges;

	static inline void reverse_value(VALUE *value) { REVERSE_VALUE::reverse(value); }
	static inline void combine_values(VALUE *value, VALUE const &incoming) { COMBINE_VALUES::combine(value, incoming); }
	static inline void split_value(VALUE const &value, VALUE *first, VALUE *second) { SPLIT_VALUE::split(value, first, second); }

	enum LineSide : int8_t {
		Left = 1,
		On = 0,
		Right = -1,
	};

	//is there a point in the interior of both a-b and a2-b2?
	bool segments_intersect(const glm::ivec2 &a, const glm::ivec2 &b, const glm::ivec2 &a2, const glm::ivec2 &b2) {
		if (a == b) return false;
		if (a2 == b2) return false;

		int32_t perp_a2 = (a2 - a).x * -(b - a).y + (a2 - a).y * (b - a).x;
		int32_t perp_b2 = (b2 - a).x * -(b - a).y + (b2 - a).y * (b - a).x;

		if (perp_a2 == 0 && perp_b2 == 0) {
			//annoying (colinear) special case
			int32_t along_a2 = (a2 - a).x * (b - a).x + (a2 - a).y * (b - a).y;
			int32_t along_b2 = (b2 - a).x * (b - a).x + (b2 - a).y * (b - a).y;
			int32_t limit = (b - a).x * (b - a).x + (b - a).y * (b - a).y;

			if (along_a2 <= 0 && along_b2 <= 0) return false;
			if (along_a2 >= limit && along_b2 >= limit) return false;
			return true;
		}

		if (perp_a2 <= 0 && perp_b2 <= 0) return false;
		if (perp_a2 >= 0 && perp_b2 >= 0) return false;

		int32_t perp_a = (a - a2).x * -(b2 - a2).y + (a - a2).y * (b2 - a2).x;
		int32_t perp_b = (b - a2).x * -(b2 - a2).y + (b - a2).y * (b2 - a2).x;

		assert(!(perp_a == 0 && perp_b == 0)); //should have been handled above

		if (perp_a <= 0 && perp_b <= 0) return false;
		if (perp_a >= 0 && perp_b >= 0) return false;

		return true;
	}

	glm::ivec2 rounded_intersection(const glm::ivec2 &a, const glm::ivec2 &b, const glm::ivec2 &a2, const glm::ivec2 &b2) {
		//NOTE: should call only when actually intersecting.
		//NOTE2: should ~probably~ use the 2x2 matrix inverse formulation here
		int32_t perp_a2 = (a2 - a).x * -(b - a).y + (a2 - a).y * (b - a).x;
		int32_t perp_b2 = (b2 - a).x * -(b - a).y + (b2 - a).y * (b - a).x;
		if (perp_a2 == perp_b2) {
			//DEBUG:
			std::cout << "a  = (" << a.x << ", " << a.y << ") b  = (" << b.x << ", " << b.y << ")\n";
			std::cout << "a2 = (" << a2.x << ", " << a2.y << ") b2 = (" << b2.x << ", " << b2.y << ")\n";
			std::cout << "perp_a2 is " << perp_a2 << "\n";
			std::cout << "perp_b2 is " << perp_b2 << "\n";
			std::cout.flush();
		}
		assert(perp_a2 != perp_b2);

		double t = double(0 - perp_a2) / double(perp_b2 - perp_a2);
		return glm::ivec2(
			std::round((b.x - a.x) * t + double(a.x)),
			std::round((b.y - a.y) * t + double(a.y))
		);
	}


	//is point 'pt' in the interior of line segment a-b?
	bool point_in_segment(const glm::ivec2 &pt, const glm::ivec2 &a, const glm::ivec2 &b) {
		glm::ivec2 ab = b - a;
		glm::ivec2 ap = pt - a;

		int32_t perp = ap.x * -ab.y + ap.y * ab.x;
		if (perp != 0) return false;
		
		int32_t along = ap.x * ab.x + ap.y * ab.y;
		if (along <= 0) return false;

		int32_t limit = ab.x * ab.x + ab.y * ab.y;
		if (along >= limit) return false;

		return true;
	}

	bool point_in_segment(const IntegerEmbeddedVertex &pt_, const IntegerEmbeddedVertex &a_, const IntegerEmbeddedVertex &b_) {
		//work in barycentric coordinates:
		glm::uvec3 common = IntegerEmbeddedVertex::common_simplex(pt_.simplex, IntegerEmbeddedVertex::common_simplex(a_.simplex, b_.simplex));

		glm::ivec2 pt = glm::ivec2(pt_.weights_on(common));
		glm::ivec2 a = glm::ivec2(a_.weights_on(common));
		glm::ivec2 b = glm::ivec2(b_.weights_on(common));

		return point_in_segment(pt, a, b);
	}

	
	uint32_t add_vertex(const IntegerEmbeddedVertex &v_) {
		IntegerEmbeddedVertex v = IntegerEmbeddedVertex::simplify(v_);

		auto &verts = simplex_vertices[v.simplex];
		auto &sedges = simplex_edges[v.simplex];

		for (auto i : verts) {
			if (vertices[i] == v) return i;
		}
		/* //PARANOIA:
		assert(v.simplex.x < v.simplex.y && v.simplex.y <= v.simplex.z);
		assert(v.weights.x + v.weights.y + v.weights.z == WEIGHT_SUM);
		for (const auto &v2 : vertices) {
			assert(!(v2 == v));
		}
		//end PARANOIA*/

		uint32_t idx = vertices.size();
		vertices.emplace_back(v);
		verts.emplace_back(idx);

		uint32_t old_size = sedges.size();

		for (uint32_t e = 0; e < old_size; ++e) {
			if (point_in_segment(v, vertices[sedges[e].first], vertices[sedges[e].second])) {
				auto second_half = sedges[e];
				second_half.first = idx;
				sedges[e].second = idx;
				sedges.emplace_back(second_half);
			}
		}
		return idx;
	}

	void add_edge(uint32_t ai, uint32_t bi, VALUE const &value) {
		assert(ai < vertices.size() && bi < vertices.size());

		if (ai == bi) return;

		const auto &a_ = vertices[ai];
		const auto &b_ = vertices[bi];
		glm::uvec3 common = IntegerEmbeddedVertex::common_simplex(a_.simplex, b_.simplex);

		glm::ivec2 a = glm::ivec2(a_.weights_on(common));
		glm::ivec2 b = glm::ivec2(b_.weights_on(common));

		//split edge at any vertices in simplex:
		auto &verts = simplex_vertices[common];
		for (auto vi : verts) {
			assert(vi < vertices.size());
			glm::ivec2 v = glm::vec2(vertices[vi].weights_on(common));
			if (point_in_segment(v, a, b)) {
				VALUE value_first, value_second;
				split_value(value, &value_first, &value_second);
				add_edge(ai, vi, value_first);
				add_edge(vi, bi, value_second);
				return;
			}
		}

		//split edge (and add new vertex) if there is an intersection:
		auto &sedges = simplex_edges[common];
		for (uint32_t e = 0; e < sedges.size(); ++e) {

			//if it matches the edge, over-write value & done!
			if (sedges[e].first == ai && sedges[e].second == bi) {
				combine_values(&sedges[e].value, value);
				return;
			}
			if (sedges[e].first == bi && sedges[e].second == ai) {
				VALUE temp = value;
				reverse_value(&temp);
				combine_values(&sedges[e].value, temp);
				return;
			}

			glm::ivec2 a2 = glm::ivec2(vertices[sedges[e].first].weights_on(common));
			glm::ivec2 b2 = glm::ivec2(vertices[sedges[e].second].weights_on(common));

			//if endpoints are interior to an existing edge, split existing edge:
			if (point_in_segment(a, a2, b2)) {
				assert(false); //THIS SHOULD NEVER HAPPEN (should have been avoided by vertex insertion?)
				auto second_half = sedges[e];
				second_half.first = ai;
				sedges.emplace_back(second_half);
				sedges[e].second = ai;
				b2 = a;
			}
			if (point_in_segment(b, a2, b2)) {
				assert(false); //THIS SHOULD NEVER HAPPEN (should have been avoided by vertex insertion?)
				auto second_half = sedges[e];
				second_half.first = bi;
				sedges.emplace_back(second_half);
				sedges[e].second = bi;
				b2 = b;
			}

			//if edges cross, remove, add intersection, and re-insert:
			if (segments_intersect(a,b, a2,b2)) {
				glm::ivec3 pt = glm::ivec3(rounded_intersection(a,b,a2,b2), 0);
				pt.z = WEIGHT_SUM - pt.x - pt.y;
				uint32_t pti = add_vertex(IntegerEmbeddedVertex(common, pt));

				float ai2 = sedges[e].first;
				float bi2 = sedges[e].second;
				VALUE value2 = sedges[e].value;
				sedges.erase(sedges.begin() + e);

				VALUE value2_first, value2_second;
				split_value(value2, &value2_first, &value2_second);
				add_edge(ai2, pti, value2_first);
				add_edge(pti, bi2, value2_second);

				VALUE value_first, value_second;
				split_value(value, &value_first, &value_second);
				add_edge(ai, pti, value_first);
				add_edge(pti, bi, value_second);

				return;
			}
			
		}

		//if got to this point, no intersections:
		sedges.emplace_back(ai, bi, value);
	}
	void add_edge(const ak::EmbeddedVertex &a, const ak::EmbeddedVertex &b, VALUE const &value) {
		uint32_t ai = add_vertex(a);
		uint32_t bi = add_vertex(b);
		add_edge(ai, bi, value);
	}

	void split_triangles(
		std::vector< glm::vec3 > const &verts, //in: mesh vertices
		std::vector< glm::uvec3 > const &tris, //in: mesh triangles
		std::vector< ak::EmbeddedVertex > *split_verts_, //out: split vertices
		std::vector< glm::uvec3 > *split_tris_, //out: split triangles
		std::vector< uint32_t > *epm_to_split_ //out: epm vertices -> split vertices
		) {

		assert(split_verts_);
		auto &split_verts = *split_verts_;
		split_verts.clear();

		assert(split_tris_);
		auto &split_tris = *split_tris_;
		split_tris.clear();

		assert(epm_to_split_);
		auto &epm_to_split = *epm_to_split_;
		epm_to_split.clear();

		split_verts.reserve(verts.size());
		for (uint32_t v = 0; v < verts.size(); ++v) {
			split_verts.emplace_back(ak::EmbeddedVertex::on_vertex(v));
		}

		epm_to_split.reserve(vertices.size());
		for (uint32_t i = 0; i < vertices.size(); ++i) {
			const auto &v = vertices[i];
			assert(v.simplex.x < verts.size());
			assert(v.weights.x + v.weights.y + v.weights.z == WEIGHT_SUM);
			if (v.simplex.y == -1U) {
				assert(v.simplex.z == -1U);
				epm_to_split.emplace_back(v.simplex.x);
			} else if (v.simplex.z == -1U) {
				assert(v.simplex.y < verts.size());
				assert(v.weights.z == 0);
				epm_to_split.emplace_back(split_verts.size());
				split_verts.emplace_back(ak::EmbeddedVertex(v.simplex, glm::vec3(v.weights) / float(WEIGHT_SUM)));
			} else {
				assert(v.simplex.y < verts.size());
				assert(v.simplex.z < verts.size());
				epm_to_split.emplace_back(split_verts.size());
				split_verts.emplace_back(ak::EmbeddedVertex(v.simplex, glm::vec3(v.weights) / float(WEIGHT_SUM)));
			}
		}

		uint32_t did_reflex = 0;
		uint32_t did_simple = 0;

		std::vector< glm::uvec3 > DEBUG_split_tris;

		//now, for each triangle, locally triangulate the planar map:
		for (const auto &tri : tris) {
			glm::uvec3 simplex = tri;
			bool need_flip = false;
			if (simplex.x > simplex.y) {
				std::swap(simplex.x, simplex.y);
				need_flip = !need_flip;
			}
			if (simplex.y > simplex.z) {
				std::swap(simplex.y, simplex.z);
				need_flip = !need_flip;
			}
			if (simplex.x > simplex.y) {
				std::swap(simplex.x, simplex.y);
				need_flip = !need_flip;
			}
			assert(simplex.x < simplex.y && simplex.y < simplex.z);

			std::vector< uint32_t > source_verts; //in split_verts
			std::vector< glm::ivec2 > coords; //[x,y] in weight space
			std::unordered_multimap< uint32_t, uint32_t > half_edges; //relative to coords

			std::unordered_map< uint32_t, uint32_t > split_verts_to_coords;
			auto ref_vert = [&](uint32_t split_vert, const glm::ivec3 &weights_on) -> uint32_t {
				auto res = split_verts_to_coords.insert(std::make_pair(split_vert, coords.size()));
				if (res.second) {
					//std::cout << "Inserting " << split_vert << " / " << weights_on.x << ", " << weights_on.y << ", " << weights_on.z << std::endl; //DEBUG
					assert(res.first->second == coords.size());
					coords.emplace_back(weights_on);
					source_verts.emplace_back(split_vert);
				}
				return res.first->second;
			};

			//add two half-edges for all internal edges:
			//std::cout << " ---- internal edges ---- " << std::endl; //DEBUG
			auto f = simplex_edges.find(simplex);
			if (f != simplex_edges.end()) {
				for (const auto &edge : f->second) {
					uint32_t a = ref_vert(epm_to_split[edge.first], vertices[edge.first].weights_on(simplex));
					uint32_t b = ref_vert(epm_to_split[edge.second], vertices[edge.second].weights_on(simplex));
					half_edges.insert(std::make_pair(a,b));
					half_edges.insert(std::make_pair(b,a));
				}
			}
			//std::cout << " ---- sides ---- " << std::endl; //DEBUG
			//add one half-edge along all sides:
			auto do_side = [&](uint32_t ifrom, uint32_t ito) {
				glm::uvec3 edge_simplex(simplex[ifrom], simplex[ito], -1U);
				if (edge_simplex.x > edge_simplex.y) std::swap(edge_simplex.x, edge_simplex.y);
				uint32_t from;
				uint32_t to;
				{
					glm::ivec3 from_weights = glm::ivec3(0);
					from_weights[ifrom] = WEIGHT_SUM;
					glm::ivec3 to_weights = glm::ivec3(0);
					to_weights[ito] = WEIGHT_SUM;
					from = ref_vert(simplex[ifrom], from_weights);
					to = ref_vert(simplex[ito], to_weights);
				}
				
				std::map< int32_t, uint32_t > weight_to_vert;
				auto f = simplex_vertices.find(edge_simplex);
				if (f != simplex_vertices.end()) {
					for (auto sv : f->second) {
						glm::ivec3 weights = vertices[sv].weights_on(simplex);
						auto res = weight_to_vert.insert(std::make_pair(weights[ito], ref_vert(epm_to_split[sv], weights)));
						assert(res.second);
					}
				}
				weight_to_vert.insert(std::make_pair(0, from));
				weight_to_vert.insert(std::make_pair(WEIGHT_SUM, to));

				assert(weight_to_vert[0] == from);
				assert(weight_to_vert[WEIGHT_SUM] == to);
				auto at = weight_to_vert.begin();
				uint32_t prev = at->second;
				++at;
				for (; at != weight_to_vert.end(); ++at) {
					half_edges.insert(std::make_pair(prev, at->second));
					prev = at->second;
				}
			};
			do_side(0,1);
			do_side(1,2);
			do_side(2,0);

			/*//DEBUG:
			std::cout << "---------------\n";
			for (uint32_t c = 0; c < coords.size(); ++c) {
				std::cout << "coords[" << c << "] = (" << coords[c].x << ", " << coords[c].y << ")\n";
			}
			for (auto e : half_edges) {
				std::cout << "  " << e.first << " -> " << e.second << "\n";
			}
			std::cout.flush();*/

			while (!half_edges.empty()) {
				std::vector< uint32_t > loop;
				loop.emplace_back(half_edges.begin()->first);
				loop.emplace_back(half_edges.begin()->second);

				//std::cout << "Seed: " << loop[0] << " -> " << loop[1] << std::endl; //DEBUG
				bool had_reflex = false;

				while (true) {
					const glm::ivec2 &from = coords[loop[loop.size()-2]];
					const glm::ivec2 &at = coords[loop.back()];
					//want to find he ccw-most exit, relative to from->at
					auto r = half_edges.equal_range(loop.back());
					assert(r.first != r.second); //ran out of edges to make a loop(?!?!)
					uint32_t best_to = -1U;
					glm::ivec2 best_d = glm::ivec2(0,0); //best == largest y/x
					int best_quad = 4;
					auto best_iter = half_edges.end();
					for (auto ri = r.first; ri != r.second; ++ri) {
						uint32_t i = ri->second;
						const glm::ivec2 &next = coords[i];
						glm::ivec2 d;
						d.x = (next - at).x * (at - from).x + (next - at).y * (at - from).y;
						d.y = (next - at).x *-(at - from).y + (next - at).y * (at - from).x;
						assert(!(d.x == 0 && d.y == 0));
						int quad;
						if (d.x <= 0 && d.y > 0) {
							quad = 0;
							d.x = -d.x;
							std::swap(d.x, d.y);
						} else if (d.x > 0 && d.y >= 0) {
							quad = 1;
						} else if (d.x >= 0 && d.y < 0) {
							d.y = -d.y;
							std::swap(d.x, d.y);
							quad = 2;
						} else if (d.x < 0 && d.y <= 0) {
							d.x = -d.x;
							d.y = -d.y;
							quad = 3;
						} else {
							assert(false);
							quad = 0; //"potentially unused local variable" warning on cl.exe
						}

						//std::cout << "  " << ri->first << " -> " << ri->second << " quad " << quad << " y/x " << d.y << "/" << d.x << std::endl; //DEBUG

						
						if (quad < best_quad || (quad == best_quad && uint64_t(d.y) * uint64_t(best_d.x) > uint64_t(best_d.y) * uint64_t(d.x))) {
							best_to = i;
							best_quad = quad;
							best_d = d;
							best_iter = ri;
						}
					}
					assert(best_to != -1U);
					if (best_quad >= 2) had_reflex = true;
					//std::cout << "Best is " << best_iter->first << " to " << best_iter->second << std::endl; //DEBUG
					if (best_iter->first == loop[0] && best_iter->second == loop[1]) {
						half_edges.erase(best_iter);
						break;
					}
					half_edges.erase(best_iter);
					loop.emplace_back(best_to);
				}
				//std::cout << "Peeled a loop of " << loop.size() << " verts." << std::endl; //DEBUG
				assert(loop[0] == loop.back());
				loop.pop_back();
				assert(loop.size() >= 3);

				if (had_reflex || true) { //DEBUG -- looks like 'had_reflex' doesn't always get set on reflex verts (specifically on doubling-back loops like a b c b d a); for now, always using the slower special-case code.
					++did_reflex;
					//std::cerr << "ERROR: no code for reflex verts in planar map yet." << std::endl;

					std::vector< uint32_t > remain = loop;
				
					while (remain.size() >= 3) {
						//std::cout << "loop:";
						//for (auto r : remain) {
						//	std::cout << " (" << coords[r].x << ", " << coords[r].y << ")";
						//}
						//std::cout << std::endl;

						bool found = false;
						for (uint32_t i = 0; i < remain.size(); ++i) {
							uint32_t previ = (i > 0 ? i - 1 : remain.size() - 1);
							uint32_t nexti = (i + 1 < remain.size() ? i + 1 : 0);
							glm::ivec2 const &prev = coords[remain[previ]];
							glm::ivec2 const &at = coords[remain[i]];
							glm::ivec2 const &next = coords[remain[nexti]];

							//check if 'at' is reflex; if it is, not an ear:
							int32_t perp_dot = -(next.y - at.y) * (prev.x - at.x) + (next.x - at.x) * (prev.y - at.y);
							if (perp_dot <= 0) continue;
							//make sure link doesn't intersect:
							bool inside = false;
							if (remain.size() > 3) {
								for (uint32_t j = 0; j < remain.size(); ++j) {
									if (remain[j] == remain[previ] || remain[j] == remain[i] || remain[j] == remain[nexti]) continue;
									glm::ivec2 const &pt = coords[remain[j]];
									int32_t right_dot = -(next.y - at.y) * (pt.x - at.x) + (next.x - at.x) * (pt.y - at.y);
									if (right_dot < 0) continue;
									int32_t left_dot  = -(at.y - prev.y) * (pt.x - prev.x) + (at.x - prev.x) * (pt.y - prev.y);
									if (left_dot < 0) continue;
									int32_t top_dot  = -(prev.y - next.y) * (pt.x - next.x) + (prev.x - next.x) * (pt.y - next.y);
									if (top_dot < 0) continue;

									//DEBUG: std::cout << "[" << previ << " " << i << " " << nexti << "] contains " << j << std::endl;

									inside = true;
									break;
								}
							}
							if (!inside) {
								found = true;
								split_tris.emplace_back(source_verts[remain[previ]], source_verts[remain[i]], source_verts[remain[nexti]]);
								if (need_flip) {
									std::swap(split_tris.back().y, split_tris.back().z);
								}
								remain.erase(remain.begin() + i);
								break;
							}
						}
						assert(found);
					}
				} else {
					++did_simple;
					//simple fan triangulation:
					for (uint32_t i = 1; i + 1 < loop.size(); ++i) {
						split_tris.emplace_back(source_verts[loop[0]], source_verts[loop[i]], source_verts[loop[i+1]]);
						if (need_flip) {
							std::swap(split_tris.back().y, split_tris.back().z);
						}
					}
				}
			}
		}

		if (did_reflex) std::cout << "  Note: used reflex-vertex special-case code in " << did_reflex << " of " << (did_reflex + did_simple) << " cases." << std::endl;

	}
};

