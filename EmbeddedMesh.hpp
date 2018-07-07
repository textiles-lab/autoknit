#pragma once

namespace ak {

template< typename T, int SUM >
struct EmbeddedVertexT {
	typedef glm::vec< 3, T, glm::highp > weights_type;
	glm::uvec3 simplex;
	weights_type weights;

	EmbeddedVertex() = default;
	EmbeddedVertex(glm::uvec3 const &simplex_, weights_type const &weights_) : simplex(simplex_), weights(weights_) { }

	bool operator==(EmbeddedVertex const &o) const {
		return simplex == o.simplex && weights == o.weights;
	}

	static EmbeddedVertex on_vertex(uint32_t a) {
		return EmbeddedVertex(glm::uvec3(a, -1U, -1U), weights_type(T(SUM), T(0), T(0)));
	}
	static EmbeddedVertex on_edge(uint32_t a, uint32_t b, T wb) {
		T wa = T(SUM) - wb;
		if (a > b) {
			std::swap(a,b);
			std::swap(wa,wb);
		}
		return EmbeddedVertex(glm::uvec3(a, b, -1U), weights_type(wa, wb, T(0)));
	}

	static EmbeddedVertex canonicalize(glm::uvec3 s, weights_type w) {
		if (s.x > s.y) { std::swap(s.x,s.y); std::swap(w.x,w.y); }
		if (s.y > s.z) { std::swap(s.y,s.z); std::swap(w.y,w.z); }
		if (s.x > s.y) { std::swap(s.x,s.y); std::swap(w.x,w.y); }

		return EmbeddedVertex(s, w);
	}

	static EmbeddedVertex mix(EmbeddedVertex const &a, EmbeddedVertex const &b, float m) {
		glm::uvec3 common = common_simplex(a.simplex, b.simplex);
		return EmbeddedVertex(
			common,
			glm::mix(a.weights_on(common), b.weights_on(common), m)
		);
	}

	template< typename T >
	T interpolate(std::vector< T > const &values) const {
		T ret = values[simplex.x] * weights.x;
		if (simplex.y != -1U) ret += values[simplex.y] * weights.y;
		if (simplex.z != -1U) ret += values[simplex.z] * weights.z;
		return ret;
	}

	glm::vec3 weights_on(glm::uvec3 simplex2) const {
		glm::vec3 ret(0,0,0);
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
		return ret;
	}


	static glm::uvec3 common_simplex(const glm::uvec3 &a, const glm::uvec3 &b) {
		glm::ivec3 ret;
		uint32_t ia = 0;
		uint32_t ib = 0;
		for (uint32_t o = 0; o < 3; ++o) {
			if (a[ia] == b[ib]) {
				ret[o] = a[ia];
				++ia; ++ib;
			} else if (a[ia] < b[ib]) {
				ret[o] = a[ia];
				++ia;
			} else { assert(a[ia] > b[ib]);
				ret[o] = b[ib];
				++ib;
			}
		}
		assert(ia == 3 || a[ia] == -1U);
		assert(ib == 3 || b[ib] == -1U);
		return ret;
	}
};

typedef EmbeddedVertexT< float > EmbeddedVertex;

}
