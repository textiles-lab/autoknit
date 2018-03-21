#pragma once

#include <glm/glm.hpp>

#include <vector>
#include <string>

// The autoknit pipeline in data formats and transformation functions.

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

// Constraints: [chains of] points on the model's surface.

struct Constraints {
	std::vector< glm::vec3 > points;
	std::vector< float > point_values;
	std::vector< glm::uvec2 > geodesics; //geodesics between points may also be constrained
};

void embed_constraints(
	Model const &model, //in: model to embed constraints on
	Constraints const &constraints, //in: constraints to embed
	Model *constrained_model, //out: model, re-triangulated with vertices at each of the constraint points
	std::vector< uint32_t > *constrainted_indices,
	std::vector< float > *constrainted_values
);


