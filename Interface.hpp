#pragma once

#include <kit/kit.hpp>

#include <kit/GLBuffer.hpp>
#include <kit/GLVertexArray.hpp>

#include "pipeline.hpp"

struct Interface : public kit::Mode {
	Interface();
	virtual ~Interface();
	virtual void update(float elapsed) override;
	virtual void draw() override;
	virtual void pointer_action(kit::PointerID pointer, kit::PointerAction action, kit::Pointer const &old_state, kit::Pointer const &new_state) override;

	//camera:
	struct Camera {
		glm::vec3 pos = glm::vec3(0.0f, 0.0f, 0.0f);
		glm::quat rot = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);
		float fovy = 60.0f;
	} camera;

	//mouse tracking:
	//TODO

	//original model:
	Model model;
	//model buffer: (vertices, normals)
	GLAttribBuffer< glm::vec3, glm::vec3 > model_triangles;
	void set_model(Model const &new_model);

};
