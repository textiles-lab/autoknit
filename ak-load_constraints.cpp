#include "pipeline.hpp"

#include <glm/gtx/norm.hpp>

#include <iostream>
#include <fstream>

template< typename S >
void read_scalar(std::istream &in, S *_out, std::string const &name) {
	assert(_out);
	if (!in.read(reinterpret_cast< char * >(_out), sizeof(S))) {
		throw std::runtime_error("Failed to read scalar " + name);
	}
}


template< typename S >
void read_vector(std::istream &in, std::vector< S > *_out, std::string const &name) {
	assert(_out);
	auto &out = *_out;
	uint32_t count;
	read_scalar(in, &count, name + " count");
	out.assign(count, S());
	if (!in.read(reinterpret_cast< char * >(out.data()), sizeof(S) * out.size())) {
		throw std::runtime_error("Failed to read vector data for " + name);
	}
}

void read_eof(std::istream &in, std::string const &name) {
	if (std::istream::traits_type::not_eof( in.get() )) {
		throw std::runtime_error("Trailing data reading " + name);
	}
}

void ak::load_constraints(
	ak::Model const &model, //in: model for vertex lookup
	std::string const &filename, //in: file to load
	std::vector< ak::Constraint > *_constraints //out: list of constraints
) {
	assert(_constraints);
	auto &constraints = *_constraints;
	constraints.clear();

	std::vector< glm::vec3 > verts;
	struct StoredConstraint {
		uint32_t verts_count;
		float value;
		float radius;
	};
	std::vector< StoredConstraint > stored_constraints;

	{ //read file:
		std::ifstream in(filename, std::ios::binary);
		read_vector(in, &verts, "verts");
		read_vector(in, &stored_constraints, "constraints");
		read_eof(in, "constraints " + filename);
	}

	uint32_t missing_verts = 0;

	{ //interpret constraints:
		uint32_t begin_vert = 0;
		for (auto const &sc : stored_constraints) {
			uint32_t end_vert = begin_vert + sc.verts_count;
			if (end_vert < begin_vert + 1 || end_vert > verts.size()) {
				throw std::runtime_error("Stored vertex size is invalid.");
			}
			constraints.emplace_back();
			constraints.back().value = sc.value;
			constraints.back().radius = std::max(0.0f, sc.radius);
			//look up verts:
			for (uint32_t vi = begin_vert; vi < end_vert; ++vi) {
				glm::vec3 const &v = verts[vi];
				uint32_t closest = -1U;
				float closest_dis = 0.01f * 0.01f;
				for (auto const &mv : model.vertices) {
					float dis = glm::length2(mv - v);
					if (dis < closest_dis) {
						closest_dis = dis;
						closest = &mv - &model.vertices[0];
					}
				}
				if (closest == -1U) {
					++missing_verts;
				} else {
					constraints.back().chain.emplace_back(closest);
				}
			}
			if (constraints.back().chain.empty()) constraints.pop_back();
		}
	}

	if (missing_verts) {
		std::cerr << "WARNING: had " << missing_verts << " missing verts loading constraints from '" << filename << "'" << std::endl;
	}
}
