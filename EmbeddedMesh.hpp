#pragma once

namespace ak {

template< typename T >
struct EmbeddedVertex {
	typedef glm::vec< 3, T, glm::highp > weights_type;
	glm::uvec3 simplex;
	weights_type weights;

	EmbeddedVertex() = default;
	EmbeddedVertex(glm::uvec3 const &simplex_, glm::vec3 const &weights_) : simplex(simplex_), weights(weights_) { }

	bool operator==(EmbeddedVertex const &o) const {
		return simplex == o.simplex && weights == o.weights;
	}

	static EmbeddedVertex on_vertex(uint32_t a) {
		return EmbeddedVertex(glm::uvec3(a, -1U, -1U), glm::vec3(1.0f, 0.0f, 0.0f));
	}
	static EmbeddedVertex on_edge(uint32_t a, uint32_t b, float mix) {
		if (a > b) {
			std::swap(a,b);
			mix = 1.0f - mix;
		}
		return EmbeddedVertex(glm::uvec3(a, b, -1U), glm::vec3(1.0f - mix, mix, 0.0f));
	}

	template< typename T >
	T interpolate(std::vector< T > const &values) const {
		T ret = values[simplex.x] * weights.x;
		if (simplex.y != -1U) ret += values[simplex.y] * weights.y;
		if (simplex.z != -1U) ret += values[simplex.z] * weights.z;
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

}
