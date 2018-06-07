#include "pipeline.hpp"

//#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
//#include <Eigen/IterativeLinearSolvers>
//#pragma GCC diagnostic pop

#include <iostream>
#include <map>

void ak::interpolate_values(
	Model const &model,
	std::vector< float > const &constraints,
	std::vector< float > *values_
) {
	assert(constraints.size() == model.vertices.size());
	assert(values_);
	auto &values = *values_;

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

	std::map< std::pair< uint32_t, uint32_t >, float > edge_weights;

	for (const auto &tri : model.triangles) {
		const glm::vec3 &a = model.vertices[tri.x];
		const glm::vec3 &b = model.vertices[tri.y];
		const glm::vec3 &c = model.vertices[tri.z];

		float weight_ab = glm::dot(a-c, b-c) / glm::length(glm::cross(a-c, b-c));
		float weight_bc = glm::dot(b-a, c-a) / glm::length(glm::cross(b-a, c-a));
		float weight_ca = glm::dot(c-b, a-b) / glm::length(glm::cross(c-b, a-b));

		edge_weights.insert(std::make_pair(std::minmax(tri.x, tri.y), 0.0f)).first->second += weight_ab;
		edge_weights.insert(std::make_pair(std::minmax(tri.y, tri.z), 0.0f)).first->second += weight_bc;
		edge_weights.insert(std::make_pair(std::minmax(tri.z, tri.x), 0.0f)).first->second += weight_ca;
	}

	//turn edge weights vector into adjacency lists:
	std::vector< std::vector< std::pair< uint32_t, float > > > adj(model.vertices.size());
	for (const auto &ew : edge_weights) {
		adj[ew.first.first].emplace_back(ew.first.second, ew.second);
		adj[ew.first.second].emplace_back(ew.first.first, ew.second);
	}


	std::vector< Eigen::Triplet< double > > coefficients;
	coefficients.reserve(model.triangles.size()*3); //more than is needed
	Eigen::VectorXd rhs(total_dofs);

	for (uint32_t i = 0; i < dofs.size(); ++i) {
		if (dofs[i] == -1U) continue;
		//sum adj[x] + one * 1 - c * x = 0.0f
		float sum = 0.0f;
		float one = 0.0f;
		for (auto a : adj[i]) {
			if (dofs[a.first] == -1U) {
				one += a.second * constraints[a.first];
			} else {
				coefficients.emplace_back(dofs[i], dofs[a.first], a.second);
			}
			sum += a.second;
		}
		coefficients.emplace_back(dofs[i], dofs[i], -sum);
		rhs[dofs[i]] = -one;
	}

	Eigen::SparseMatrix< double > A(total_dofs, total_dofs);
	A.setFromTriplets(coefficients.begin(), coefficients.end());
	//A = A * A.transpose();
	A.makeCompressed(); //redundant?

	//Eigen::SparseLU< Eigen::SparseMatrix< double > > solver;
	//Eigen::SparseQR< Eigen::SparseMatrix< double >, Eigen::COLAMDOrdering< int > > solver;
	Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > solver;
	//Eigen::ConjugateGradient< Eigen::SparseMatrix< double > > solver;
	solver.compute(A);
	if (solver.info() != Eigen::Success) {
		std::cerr << "Decomposition failed." << std::endl;
		exit(1);
	}
	Eigen::VectorXd x = solver.solve(rhs);
	if (solver.info() != Eigen::Success) {
		std::cerr << "Solving failed." << std::endl;
		exit(1);
	}
	//std::cout << solver.iterations() << " interations later..." << std::endl; //DEBUG
	//std::cout << solver.error() << " (estimated error)..." << std::endl; //DEBUG

	values = constraints;
	for (uint32_t i = 0; i < dofs.size(); ++i) {
		if (dofs[i] != -1U) values[i] = x[dofs[i]];
	}

}
