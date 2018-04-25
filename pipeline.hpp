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

//Given list of constraints, properly trim and constrain a model:
void embed_constraints(
	Model const &model,
	std::vector< Constraint > const &constraints,
	Model *constrained_model,
	std::vector< float > *constrained_values, //same size as out_model's vertices
	std::vector< std::vector< glm::vec3 > > *DEBUG_chain_paths = nullptr //out, optional: paths computed for each constraint chain
);

//Given list of values, fill missing with as-smooth-as-possible interpolation:
void interpolate_values(
	Model const &model, //in: model to embed constraints on
	std::vector< float > const &constraints, //same size as model.vertices; if non-NaN, fixes value
	std::vector< float > *values //smooth interpolation of (non-NaN) constraints
);

//make time field:
void interpolate_constraints(
	//in: mesh, embedded constraints(?)
	//out: time field
);

//peel/link to create embedded row/column meshes:

struct EmbeddedVertex {
	glm::uvec3 vertices;
	glm::vec3 weights;
	EmbeddedVertex() = default;
	EmbeddedVertex(glm::uvec3 const &vertices_, glm::vec3 const &weights_) : vertices(vertices_), weights(weights_) { }

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
};

//NOTE: chains represented as [a,b,c,d,a] if a loop, [a,b,c,d] if a chain. Always represented in CCW order.

enum Flag : int8_t {
	FlagDiscard = -1,
	FlagLinkNone = 0, //for non-stitches or already-linked stitches
	FlagLinkOne  = 1, //for short-row ends
	FlagLinkAny  = 2,
};

void find_first_active_chains(
	std::vector< glm::uvec3 > const &triangles, //in: list of model triangles
	std::vector< float > const &times,          //in: time field (times @ vertices)
	std::vector< std::vector< EmbeddedVertex > > *active_chains, //out: all mesh boundaries that contain a minimum
	std::vector< std::vector< Flag > > *active_flags //out: all mesh boundaries that contain a minimum
);

void peel_chains(
	std::vector< glm::vec3 > const &vertices,   //in: list of model vertices
	std::vector< glm::uvec3 > const &triangles, //in: list of model triangles
	std::vector< float > const &times,          //in: time field (times @ vertices)
	std::vector< std::vector< EmbeddedVertex > > const &active_chains, //in: current active chains
	std::vector< std::vector< EmbeddedVertex > > *next_chains //out: next chains (may be different size than active_chains)
);

struct Link {
	uint32_t from_chain, from_vertex;
	uint32_t to_chain, to_vertex;
};

void link_chains(
	std::vector< glm::vec3 > const &vertices,   //in: list of model vertices
	std::vector< glm::uvec3 > const &triangles, //in: list of model triangles
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
