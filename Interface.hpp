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
		DragConsPt,
		DragConsRadius,
	};
	Drag drag = DragNone;
	struct {
		uint32_t cons = -1U;
		uint32_t cons_pt = -1U;
		void clear() {
			cons = -1U;
			cons_pt = -1U;
		}
	} dragging;

	enum {
		ShowModel,
		ShowConstrainedModel,
	} show = ShowModel;

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
		uint32_t cons_pt = -1U;

		void clear() {
			point = glm::vec3(std::numeric_limits< float >::quiet_NaN());
			tri = -1U;
			coords = glm::vec3(std::numeric_limits< float >::quiet_NaN());
			vert = -1U;
			cons = -1U;
			cons_pt = -1U;
		}
	} hovered;

	void update_hovered();

	//parameters:
	ak::Parameters parameters;

	//original model:
	ak::Model model;
	void set_model(ak::Model const &model);
	//model buffer: (vertices, normals, ids)
	void update_model_triangles();
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > model_triangles;
	GLVertexArray model_triangles_for_model_draw;

	//constraints:
	std::vector< ak::Constraint > constraints;
	void set_constraints(std::vector< ak::Constraint > const &constraints);
	void update_constraints();
	bool constraints_dirty = true;
	std::vector< std::vector< glm::vec3 > > DEBUG_constraint_paths;
	std::vector< std::vector< glm::vec3 > > DEBUG_constraint_loops;
	ak::Model constrained_model;
	std::vector< float > constrained_values;

	//interpolation (also computed by update_constraints):
	std::vector< float > interpolated_values;

	//peeling information:
	void start_peeling();
	std::vector< std::vector< ak::EmbeddedVertex > > active_chains;
	std::vector< std::vector< ak::Flag > > active_flags;

	void step_peeling();
	std::vector< std::vector< ak::EmbeddedVertex > > next_chains;


	std::string save_constraints_file = ""; //if not "", will save constraints to this file after every change
	void save_constraints();

	void update_DEBUG_constraint_paths_tristrip();
	//position, normal, color:
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > DEBUG_constraint_paths_tristrip;
	GLVertexArray DEBUG_constraint_paths_tristrip_for_path_draw;

	void update_DEBUG_constraint_loops_tristrip();
	//position, normal, color:
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > DEBUG_constraint_loops_tristrip;
	GLVertexArray DEBUG_constraint_loops_tristrip_for_path_draw;


	void update_constrained_model_triangles();
	//constrained model buffer: position, normal, id, texcoord
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4, glm::vec2 > constrained_model_triangles;
	GLVertexArray constrained_model_triangles_for_textured_draw;

	void update_active_chains_tristrip();
	//position, normal, color:
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > active_chains_tristrip;
	GLVertexArray active_chains_tristrip_for_path_draw;

	void update_next_chains_tristrip();
	//position, normal, color:
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > next_chains_tristrip;
	GLVertexArray next_chains_tristrip_for_path_draw;


	//place camera to view whole model:
	void reset_camera();

};
