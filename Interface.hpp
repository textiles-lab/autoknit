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
	virtual void handle_event(SDL_Event const &) override;

	//framebuffers + textures
	glm::uvec2 fb_size = glm::uvec2(0);

	GLuint color_id_fb = 0; //(color_tex, id_tex) + depth_tex
	GLuint color_tex = 0;
	GLuint id_tex = 0;
	GLuint depth_tex = 0; //depth values (for soft edges)

	void alloc_fbs(); //(re-)allocate framebuffers given current display size

	//camera:
	struct Camera {
		glm::vec3 center = glm::vec3(0.0f);
		float radius = 5.0f;
		//all in radians:
		float azimuth = glm::radians(30.0f);
		float elevation = glm::radians(45.0f);
		float fovy = glm::radians(60.0f);

		glm::vec3 at() const {
			float ca = std::cos(azimuth);
			float sa = std::sin(azimuth);
			float ce = std::cos(elevation);
			float se = std::sin(elevation);
			glm::vec3 out   = glm::vec3( ce * ca, ce * sa, se);
			glm::vec3 at = center + out * radius;
			return at;
		}
		//matrix that takes positions to camera space:
		glm::mat4 mv() const {
			float ca = std::cos(azimuth);
			float sa = std::sin(azimuth);
			float ce = std::cos(elevation);
			float se = std::sin(elevation);
			glm::vec3 right = glm::vec3(     -sa,      ca, 0.0f);
			glm::vec3 up    = glm::vec3(-se * ca,-se * sa, ce);
			glm::vec3 out   = glm::vec3( ce * ca, ce * sa, se);
			glm::vec3 at = center + out * radius;
			return glm::mat4(
				right.x, up.x, out.x, 0.0f,
				right.y, up.y, out.y, 0.0f,
				right.z, up.z, out.z, 0.0f,
				glm::dot(-at, right), glm::dot(-at, up), glm::dot(-at, out), 1.0f
			);
		}
		glm::mat4 mvp() const {
			return glm::infinitePerspective(fovy, kit::display.aspect, 0.1f) * mv();
		}
	} camera;

	//mouse tracking:
	enum Drag {
		DragNone = 0,
		DragCamera,
		DragCameraFlipX, //for dragging when upside-down
		DragCameraPan,
	};
	Drag drag = DragNone;

	struct {
		glm::vec2 at = glm::vec2(std::numeric_limits< float >::quiet_NaN());
		bool moved = false;
	} mouse;

	struct {
		glm::vec3 point = glm::vec3(std::numeric_limits< float >::quiet_NaN());
		uint32_t tri = -1U;
		glm::vec3 coords = glm::vec3(std::numeric_limits< float >::quiet_NaN());
		uint32_t vert = -1U;

		uint32_t cons = -1U;

		void clear() {
			point = glm::vec3(std::numeric_limits< float >::quiet_NaN());
			tri = -1U;
			coords = glm::vec3(std::numeric_limits< float >::quiet_NaN());
			vert = -1U;
			cons = -1U;
		}
	} hovered;

	void update_hovered();

	//original model:
	ak::Model model;
	void set_model(ak::Model const &model);
	//model buffer: (vertices, normals, ids)
	void update_model_triangles();
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > model_triangles;
	GLVertexArray model_triangles_for_model_draw;

	//constraints:
	//ak::Constraints constraints;
	void set_constraints(ak::Constraints const &constraints);
	std::vector< std::pair< uint32_t, float > > constrained_vertices;
	//std::vector< std::vector< std::pair< uint32_t, float > > > constrained_paths;
	/*
	//constrained model:
	void update_constrained_model();
	ak::Model constrained_model;
	std::vector< float > constrained_values;
	//constrained model buffer: (vertices [w is value], normals, ids)
	void update_constrained_model_triangles();
	GLAttribBuffer< glm::vec4, glm::vec3, glm::u8vec4 > constrained_model_triangles;
	GLVertexArray constrained_model_triangles_for_model_draw;
	*/


	//place camera to view whole model:
	void reset_camera();

};
