#include "pipeline.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

void load_obj(
	std::string const &file,
	Model *model_
) {
	assert(model_);
	auto &model = *model_;

	model.vertices.clear();
	model.triangles.clear();

	std::ifstream in(file);

	uint32_t tri_faces = 0;

	std::string line;
	while (std::getline(in, line)) {
		{ //strip comments
			auto idx = line.find('#');
			if (idx != std::string::npos) {
				line = line.substr(0, idx);
			}
		}
		std::istringstream iss(line);
		std::vector< std::string > tokens;
		{
			std::string tok;
			while (iss >> tok) {
				tokens.emplace_back(tok);
			}
		}
		if (tokens.empty()) continue; //blank line

		if (tokens[0] == "v") {
			//vertex position (x,y,z, [w])
			if (!(tokens.size() == 3 || tokens.size() == 4)) {
				throw std::runtime_error("'v' command should have 3 or 4 arguments");
			}
			glm::vec3 pos;
			pos.x = std::stof(tokens[1]);
			pos.y = std::stof(tokens[2]);
			pos.z = std::stof(tokens[3]);
			model.vertices.emplace_back(pos);
		} else if (tokens[0] == "f") {
			//face
			if (tokens.size() < 4) {
				throw std::runtime_error("'f' command should have at least 3 arguments");
			}
			if (tokens.size() > 4) {
				++tri_faces;
			}

			glm::ivec3 tri;
			tri.x = std::stoul(tokens[1]);
			if (tri.x < 0) tri.x += model.vertices.size();
			tri.y = std::stoul(tokens[2]);
			if (tri.y < 0) tri.y += model.vertices.size();
			//turn face into a triangle fan:
			for (uint32_t i = 3; i < tokens.size(); ++i) {
				tri.z = std::stoul(tokens[3]);
				if (tri.z < 0) tri.z += model.vertices.size();
				model.triangles.emplace_back(tri);
				tri.y = tri.z;
			}
		} else if (tokens[0] == "vt") {
			//"texture coordinate" -- ignored
		} else if (tokens[0] == "vn") {
			//"vertex normal" -- ignored
		} else if (tokens[0] == "vp") {
			//"parameter-space vertex" -- ignored
		} else if (tokens[0] == "l") {
			//"line" -- ignored
		} else if (tokens[0] == "mtllib") {
			//"material template library" -- ignored
		} else if (tokens[0] == "usemtl") {
			//"material name" -- ignored
		} else if (tokens[0] == "o") {
			//"object" -- ignored
		} else if (tokens[0] == "g") {
			//"group" -- ignored
		} else {
			std::cerr << "WARNING: unknown obj command '" << tokens[0] << "'" << std::endl;
		}
	}

	std::cout << "Read " << model.vertices.size() << " vertices and " << model.triangles.size() << " triangles from '" << file << "'." << std::endl;
	if (tri_faces) {
		std::cerr << "WARNING: had to triangulate " << tri_faces << " faces." << std::endl;
	}

	//validate + properly index triangles:
	for (auto &t : model.triangles) {
		for (uint32_t i = 0; i < 3; ++i) {
			if (t[i] < 1 || t[i] > model.vertices.size()) {
				throw std::runtime_error("Invalid triangle index.");
			}
			t[i] -= 1;
		}
	}
}
