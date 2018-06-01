#include "pipeline.hpp"

#include <Eigen/SparseQR>

void ak::interpolate_values(
	Model const &model,
	std::vector< float > const &constraints,
	std::vector< float > *values_
) {
	assert(constraints.size() == model.vertices.size());
	assert(values_);
	auto &values = *values_;
	values = constraints;

	std::vector< uint32_t > dofs;
	dofs.reserve(constraints.size());
	uint32_t total_dofs = 0;
	for (auto c : constraints) {
		if (c == c) dofs.emplace_back(-1U);
		else dofs.emplace_back(total_dofs++);
	}

	std::cout << "Have " << total_dofs << " degrees of freedom and " << (constraints.size() - total_dofs) << " constraints." << std::endl;

	if (total_dofs == constraints.size()) {
		throw std::runtime_error("Cannot interpolate from no constraints.");
	}

	std::unordered_map< glm::uvec2, float > edge_weights;

	for (const auto &tri : model.triangles) {
		const glm::vec3 &a = model.vertices[tri.x];
		const glm::vec3 &b = model.vertices[tri.y];
		const glm::vec3 &c = model.vertices[tri.z];

		float weight_ab = glm::dot(a-c, b-c) / glm::length(glm::cross(a-c, b-c));
		float weight_bc = glm::dot(b-a, c-a) / glm::length(glm::cross(b-a, c-a));
		float weight_ca = glm::dot(c-b, a-b) / glm::length(glm::cross(c-b, a-b));

		//<---- I WAS HERE
	}

	std::vector< Eigen:Triplet< double > > coefficients;
	coefficients.reserve(model.triangles.size()*3); //more than is needed



}
