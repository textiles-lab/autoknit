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
		glm::vec3 center = glm::vec3(0.0f);
		float radius = 5.0f;
		//all in radians:
		float azimuth = glm::radians(30.0f);
		float elevation = glm::radians(45.0f);
		float fovy = glm::radians(60.0f);

		//matrix that takes positions to camera space:
		glm::mat4 mv() const {
			float ca = std::cos(azimuth);
			float sa = std::sin(azimuth);
			float ce = std::cos(elevation);
			float se = std::sin(elevation);
			glm::vec3 right = glm::vec3(     -sa,      ca, 0.0f);
			glm::vec3 up    = glm::vec3(-se * ca,-se * sa, ce);
			glm::vec3 out   = glm::vec3( ce * ca, ce * sa, se);
			glm::vec3 at = out * radius;
			return glm::mat4(
				right.x, up.x, out.x, 0.0f,
				right.y, up.y, out.y, 0.0f,
				right.z, up.z, out.z, 0.0f,
				glm::dot(-at, right), glm::dot(-at, up), glm::dot(-at, out), 1.0f
			);
		}
	} camera;

	//mouse tracking:
	enum Drag {
		DragNone = 0,
		DragCamera,
		DragCameraFlipX, //for dragging when upside-down
	};
	Drag drag = DragNone;

	//original model:
	Model model;
	//model buffer: (vertices, normals)
	GLAttribBuffer< glm::vec3, glm::vec3 > model_triangles;
	GLVertexArray model_triangles_for_model_draw;
	void set_model(Model const &new_model);

	//place camera to view whole model:
	void reset_camera();

};
