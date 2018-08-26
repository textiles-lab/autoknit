#pragma once

#include <glm/glm.hpp>

#include <vector>
#include <string>
#include <algorithm>

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

	//sample spacing for sample_chain:
	float get_chain_sample_spacing() const {
		return 0.25f * stitch_width_mm / model_units_mm;
	}

	//edge sample spacing for embedded_path:
	float get_max_path_sample_spacing() const {
		return 0.02f * std::min(stitch_width_mm, 2.0f * stitch_height_mm) / model_units_mm;
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

	static EmbeddedVertex canonicalize(glm::uvec3 s, glm::vec3 w) {
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

struct Stitch {
	float t; //position along chain [0,1)
	uint32_t vertex; //used to track stitches when building RowColGraph
	enum Flag : char {
		FlagDiscard = -1,
		FlagLinkOne =  1,
		FlagLinkAny =  2,
	};
	Flag flag; //what sort of links are allowed to this stitch (or if this stitch is just marked for removal)
	Stitch(float t_, Flag flag_, uint32_t vertex_ = -1U) : t(t_), vertex(vertex_), flag(flag_) { }
};


struct RowColGraph {
	struct Vertex {
		EmbeddedVertex at;
		uint32_t row_in = -1U;
		uint32_t row_out = -1U;
		uint32_t col_in[2] = {-1U, -1U};
		uint32_t col_out[2] = {-1U, -1U};
		void add_col_in(uint32_t i) {
			if (col_in[0] == -1U) col_in[0] = i;
			else if (col_in[1] == -1U) col_in[1] = i;
			else assert(col_in[0] == -1U || col_in[1] == -1U); //no room!
		}
		void add_col_out(uint32_t i) {
			if (col_out[0] == -1U) col_out[0] = i;
			else if (col_out[1] == -1U) col_out[1] = i;
			else assert(col_out[0] == -1U || col_out[1] == -1U); //no room!
		}
	};
	std::vector< Vertex > vertices;
	void clear() {
		vertices.clear();
	}
};


//helper used by the find_first and peel methods to sample chain at sub-stitch length scale:
void sample_chain(
	float vertex_spacing, //in: maximum space between chain vertices
	Model const &model, //in: model over which chain is defined
	std::vector< EmbeddedVertex > const &chain, //in: chain to be sampled
	std::vector< EmbeddedVertex > *sampled_chain //out: sub-sampled chain
	//std::vector< Flag > *sampled_flags //out: flags (possibly with linkNone on points needed in chain for consistency but not sampled?)
);

//helper:
// remove all of the model that is:
//  - right of anything in left_of
//  - left of anything in right_of.
// so this leaves things that are left of left_of and right of right_of
// (the double-negative definition above is used because it makes more sense in the case of empty components)
void trim_model(
	Model const &model, //in: model
	std::vector< std::vector< EmbeddedVertex > > const &left_of,
	std::vector< std::vector< EmbeddedVertex > > const &right_of,
	Model *clipped, //out: portion of model's surface that is left_of the left_of chains and right_of the right_of chains
	std::vector< EmbeddedVertex > *clipped_vertices,//out: map from clipped vertices to source mesh
	std::vector< std::vector< uint32_t > > *left_of_vertices = nullptr, //out (optional): indices of vertices corresponding to left_of chains [may be some rounding]
	std::vector< std::vector< uint32_t > > *right_of_vertices = nullptr //out (optional): indices of vertices corresponding to right_of chains [may be some rounding]
);

//The first active chains are along boundaries that are at minimums:
void find_first_active_chains(
	Parameters const &parameters,
	Model const &model, //in: model
	std::vector< float > const &times,          //in: time field (times @ vertices)
	std::vector< std::vector< EmbeddedVertex > > *active_chains, //out: all mesh boundaries that contain a minimum
	std::vector< std::vector< Stitch > > *active_stitches, //out: evenly-spaced stitch locations along boundaries.
	RowColGraph *graph = nullptr //in/out: graph to update [optional]
);

void peel_slice(
	Parameters const &parameters,
	Model const &model, //in: model
	std::vector< std::vector< EmbeddedVertex > > const &active_chains, //in: current active chains
	Model *slice, //out: slice of model from active chains to next chains
	std::vector< EmbeddedVertex > *slice_on_model, //out: map from slice vertices to model vertices
	std::vector< std::vector< uint32_t > > *slice_active_chains, //out: active chains on slice
	std::vector< std::vector< uint32_t > > *slice_next_chains, //out: next chains on slice
	std::vector< bool > *used_boundary = nullptr //out:does slice_next_chains[i] include part of a boundary?
);

struct Link {
	uint32_t from_chain, from_stitch;
	uint32_t to_chain, to_stitch;
};

//NOTE:
// in the output of link_chains, links are ordered so that if a vertex appears
// in more than one link, the links in which it appears are adjacent, and the
// order is increasing along the opposite chain.
// e.g.
//  -------a--->
//       / |
//  ----b--c--->
//  the links array will contain
//        ... (a,b) (a,c) ...
//     or ... (b,a) (c,a) ...
//  (depending on which chain is 'from' and which is 'to')
void link_chains(
	Parameters const &parameters,
	Model const &slice, //in: slice on which the chains reside
	std::vector< float > const &slice_times, //in: time field (times @ vertices), for slice
	std::vector< std::vector< uint32_t > > const &active_chains, //in: current active chains (slice vertex #'s)
	std::vector< std::vector< Stitch > > const &active_stitches, //in: current active stitches, sorted by time
	std::vector< std::vector< uint32_t > > const &next_chains, //in: current next chains (slice vertex #'s)
	std::vector< bool > const &next_used_boundary, //in: did next chain use boundary? (forces no discard)
	//need this or slice_times (above) std::vector< std::vector< bool > > const &discard_segments,
	std::vector< std::vector< Stitch > > *next_stitches, //out: next active stitches
	std::vector< Link > *links //out: active_chains[from_chain][from_vertex] -> linked_next_chains[to_chain][to_vertex] links
);

//helper:
// find a shortest path between two embedded vertices
//  note: will throw runtime_error if no path exists
void embedded_path(
	Parameters const &parameters,
	Model const &model,
	EmbeddedVertex const &source,
	EmbeddedVertex const &target,
	std::vector< EmbeddedVertex > *path //out: path; path[0] will be source and path.back() will be target
);

//maybe:
//void reposition_course -- somehow update next course's position based on linking result.


void build_next_active_chains(
	Parameters const &parameters,
	Model const &slice,
	std::vector< EmbeddedVertex > const &slice_on_model, //in: vertices of slice (on model)
	std::vector< std::vector< uint32_t > > const &active_chains,  //in: current active chains (on slice)
	std::vector< std::vector< Stitch > > const &active_stitches, //in: current active stitches
	std::vector< std::vector< uint32_t > > const &next_chains, //in: next chains (on slice)
	std::vector< std::vector< Stitch > > const &next_stitches, //in: next stitches
	std::vector< bool > const &next_used_boundary, //in: next chains used boundary?
	std::vector< Link > const &links, //in: links between active and next
	std::vector< std::vector< EmbeddedVertex > > *next_active_chains, //out: next active chains (on model)
	std::vector< std::vector< Stitch > > *next_active_stitches, //out: next active stitches
	RowColGraph *graph = nullptr //in/out: graph to update [optional]
);


struct TracedStitch {
	uint32_t yarn = -1U; //yarn ID (why is this on a yarn_in? I guess the schedule.cpp code will tell me someday.
	//ins and outs are in construction order (OLD was: CW direction):
	uint32_t ins[2] = {-1U, -1U};
	uint32_t outs[2] = {-1U, -1U};
	enum Type : char {
		None = '\0',
		Start = 's',
		End = 'e',
		Tuck = 't',
		Miss = 'm',
		Knit = 'k',
		//I am going to add these because the scheduling re-write cares about them, though I'm not sure if they are a good idea to have in a general sense:
		Increase = 'i',
		Decrease = 'd',
	} type = None;
	enum Dir : char {
		CW = 'c', Clockwise = CW,
		AC = 'a', Anticlockwise = AC, CCW = AC, Counterclocwise = AC,
	} dir = AC;

	//useful for debugging and visualization:
	uint32_t vertex = -1U; //vertex of rowcolgraph where created
	glm::vec3 at = glm::vec3(std::numeric_limits< float >::quiet_NaN());
};

//cycles -> stitches

void trace_graph(
	RowColGraph const &graph, //in: row-column graph
	std::vector< TracedStitch > *traced, //out:traced list of stitches
	Model *DEBUG_model = nullptr //in (optional): model; stitches' .at will be set using its vertices
);

void schedule_stitches(
	std::vector< TracedStitch > const &stitches
	//in: list of stitches
	//out: knitout-ish
);

} //namespace ak
