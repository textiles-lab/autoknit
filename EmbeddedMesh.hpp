#pragma once

struct EmbeddedVertex {
	//EmbeddedVertex is either:
	//inds = (a,-1U,-1U) //on a vertex
	//inds = (a,b,-1U) //on an edge
	//inds = (a,b,c) //on a triangle
	glm::uvec3 inds = glm::uvec3(-1U);
	glm::vec3 weights = glm::vec3(0.0f); //weights always sum to 1.0f
};

struct EmbeddedMesh {
	std::vector< EmbeddedVertex > vertices;
	std::vector< glm::uvec2 > edges;
};
