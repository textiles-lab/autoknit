#pragma once

#include <glm/glm.hpp>

#include <vector>
#include <string>

// The autoknit pipeline in data formats and transformation functions.

namespace ak {

// Input model: vertices and triangles, loaded from an .obj file:
struct Model {
	std::vector< glm::vec3 > vertices;
	std::vector< glm::uvec3 > triangles;
	void clear() { triangles.clear(); vertices.clear(); }
};

//Load an object file into a Model structure.
//NOTE: throws on error
void load_obj(
	std::string const &file, //in: file to load
	Model *model //out: model to fill with loaded data
);

// Constraint, stored as [chain of] points on the model's surface.
struct Constraint {
	std::vector< uint32_t > chain;
	float value = 0.0f;
	float radius = 0.0f;
};

// Parameters: used to influence various steps
struct Parameters {
	//stitch size in millimeters:
	float stitch_width_mm = 3.66f;
	float stitch_height_mm = 1.73f;

	//model unit size in millimeters:
	float model_units_mm = 1.0f;

	//maximum edge length for embed_constraints:
	float get_max_edge_length() const {
		return 0.5f * std::min(stitch_width_mm, 2.0f * stitch_height_mm) / model_units_mm;
	}

	//potential stitch spacing for sample_chain:
	float get_chain_sample_spacing() const {
		return 0.25f * stitch_width_mm / model_units_mm;
	}
};

//Load list of constraints from a file
//NOTE: throws on error
void load_constraints(
	Model const &model, //in: model for vertex lookup
	std::string const &file, //in: file to load
	std::vector< Constraint > *constraints //out: list of constraints
);

//Save list of constraints to a file:
void save_constraints(
	Model const &model, //in: model for vertex lookup
	std::vector< Constraint > const &constraints, //in: list of constraints
	std::string const &file //in: file name to save to
);

//Given list of constraints, properly trim (maybe remesh?) and constrain a model:
//(in the case of no constraints, returns the input model with all-NaN constrained values)
void embed_constraints(
	Parameters const &parameters,
	Model const &model,
	std::vector< Constraint > const &constraints,
	Model *constrained_model,
	std::vector< float > *constrained_values, //same size as out_model's vertices
	std::vector< std::vector< glm::vec3 > > *DEBUG_chain_paths = nullptr, //out, optional: paths computed for each constraint chain
	std::vector< std::vector< glm::vec3 > > *DEBUG_chain_loops = nullptr //out, optional: loops around the paths computed for radius-based trimming
);

//Given list of values, fill missing with as-smooth-as-possible interpolation:
void interpolate_values(
	Model const &model, //in: model to embed constraints on
	std::vector< float > const &constraints, //same size as model.vertices; if non-NaN, fixes value
	std::vector< float > *values //smooth interpolation of (non-NaN) constraints
);

//peel/link to create embedded row/column meshes:

struct EmbeddedVertex {
	glm::uvec3 simplex;
	glm::vec3 weights;
	EmbeddedVertex() = default;
	EmbeddedVertex(glm::uvec3 const &simplex_, glm::vec3 const &weights_) : simplex(simplex_), weights(weights_) { }

	bool operator==(EmbeddedVertex const &o) const {
		return simplex == o.simplex && weights == o.weights;
	}
	bool operator!=(EmbeddedVertex const &o) const {
		return !(*this == o);
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

//helper: extract embedded level sets given values at vertices:
//NOTE: chain orientation is along +x (if values increase along +y)
void extract_level_chains(
	Model const &model, //in: model on which to embed vertices
	std::vector< float > const &values, //in: values at vertices
	float const level, //in: level at which to extract chains
	std::vector< std::vector< EmbeddedVertex > > *chains //chains of edges at given level
);

//NOTE: chains represented as [a,b,c,d,a] if a loop, [a,b,c,d] if a chain. Always represented in CCW order.

//an active chain is a list of embedded vertices on the mesh. If it is a loop, the first and last vertex are the same.
// active chains are accompanied by flags for each vertex, indicating whether that vertex has been selected to be a stitch.

enum Flag : int8_t {
	FlagDiscard = -1,
	FlagLinkNone = 0, //for non-stitches or already-linked stitches
	FlagLinkOne  = 1, //for short-row ends -- must link to and from one stitch only
	FlagLinkAny  = 2,
};

//helper used by the find_first and peel methods to sample chain at sub-stitch length scale:
void sample_chain(
	float vertex_spacing, //in: maximum space between chain vertices
	Model const &model, //in: model over which chain is defined
	std::vector< EmbeddedVertex > const &chain, //in: chain to be sampled
	std::vector< EmbeddedVertex > *sampled_chain //out: sub-sampled chain
	//std::vector< Flag > *sampled_flags //out: flags (possibly with linkNone on points needed in chain for consistency but not sampled?)
);

//The first active chains are along boundaries that are at minimums:
void find_first_active_chains(
	Parameters const &parameters,
	Model const &model, //in: model
	std::vector< float > const &times,          //in: time field (times @ vertices)
	std::vector< std::vector< EmbeddedVertex > > *active_chains, //out: all mesh boundaries that contain a minimum
	std::vector< std::vector< Flag > > *active_flags //out: all mesh boundaries that contain a minimum
);

void peel_chains(
	Parameters const &parameters,
	Model const &model, //in: model
	std::vector< float > const &times,          //in: time field (times @ vertices)
	std::vector< std::vector< EmbeddedVertex > > const &active_chains, //in: current active chains
	std::vector< std::vector< EmbeddedVertex > > *next_chains //out: next chains (may be different size than active_chains)
);

struct Link {
	uint32_t from_chain, from_vertex;
	uint32_t to_chain, to_vertex;
};

void link_chains(
	Parameters const &parameters,
	Model const &model, //in: model
	std::vector< float > const &times,          //in: time field (times @ vertices)
	std::vector< std::vector< EmbeddedVertex > > const &active_chains, //in: current active chains
	std::vector< std::vector< Flag > > const &active_flags, //in: stitches
	std::vector< std::vector< EmbeddedVertex > > const &next_chains, //in: next chains
	std::vector< std::vector< EmbeddedVertex > > *linked_next_chains, //out: next chains
	std::vector< std::vector< Flag > > *linked_next_flags, //out: flags indicating status of vertices on next chains
	std::vector< Link > *links //out: active_chains[from_chain][from_vertex] -> linked_next_chains[to_chain][to_vertex] links
);

//maybe:
//void reposition_course -- somehow update next course's position based on linking result.

void build_next_active_chains(
	Parameters const &parameters,
	Model const &model,
	std::vector< std::vector< EmbeddedVertex > > const &active_chains, //in: current active chains
	std::vector< std::vector< Flag > > const &active_flags, //in: flags for current active
	std::vector< std::vector< EmbeddedVertex > > const &next_chains, //in: next chains
	std::vector< std::vector< Flag > > const &next_flags, //in: flags for next active
	std::vector< Link > const &links, //in: links between active and next
	std::vector< std::vector< EmbeddedVertex > > *next_active_chains, //out: next active chains
	std::vector< std::vector< Flag > > *next_active_flags //out: next stitch flags
);

//probably is in driver code (different cases for build_first and build_next):
//void extract_stitches(

struct Stitch {
	EmbeddedVertex at;
	uint32_t yarn_id = -1U;
	uint32_t ins[2] = {-1U, -1U};
	uint32_t outs[2] = {-1U, -1U};
	enum Type : uint8_t {
		None = 0,
	} type = None;
};

//cycles -> stitches

void trace_stitches(
	std::vector< EmbeddedVertex > const &vertices, //in: copied to stitches
	std::vector< glm::uvec2 > const &course_links, //in: ccw-oriented yarn edges
	std::vector< glm::uvec2 > const &wale_links, //in: early-to-late oriented loop edges
	std::vector< Stitch > *traced
	//out: list of stitches
);

void schedule_stitches(
	std::vector< Stitch > const &stitches
	//in: list of stitches
	//out: knitout-ish
);

} //namespace ak
