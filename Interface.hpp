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

	enum : uint32_t {
		ShowModel            = (1 << 0),
		ShowConstraints      = (1 << 1),
		ShowTimesModel       = (1 << 2),
		ShowSlice            = (1 << 3),
		ShowActiveChains     = (1 << 4),
		ShowSliceChains      = (1 << 5),
		ShowLinks            = (1 << 6),
		ShowNextActiveChains = (1 << 7),
		ShowRowColGraph      = (1 << 8),
		ShowTraced           = (1 << 9),

		ShowModelBits = ShowModel | ShowTimesModel | ShowSlice,
		ShowStepBits = ShowActiveChains | ShowSliceChains | ShowLinks | ShowNextActiveChains,
		ShowGraphBits = ShowRowColGraph | ShowTraced,
	};
	uint32_t show = ShowModel | ShowConstraints | ShowRowColGraph;


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

	//-------------------------------
	//original model:
	ak::Model model;
	void set_model(ak::Model const &model);

	//visualization:
	//model buffer: (vertices, normals, ids)
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > model_triangles;
	GLVertexArray model_triangles_for_model_draw;
	bool model_triangles_dirty = true;
	void update_model_triangles();

	//-------------------------------
	//constraints:
	std::vector< ak::Constraint > constraints;
	ak::Model constrained_model;
	std::vector< float > constrained_values;
	std::vector< std::vector< glm::vec3 > > DEBUG_constraint_paths;
	std::vector< std::vector< glm::vec3 > > DEBUG_constraint_loops;
	void clear_constraints();

	void set_constraints(std::vector< ak::Constraint > const &constraints);
	bool constraints_dirty = true;
	void update_constraints();

	//data wrangling:
	std::string save_constraints_file = ""; //if not "", will save constraints to this file after every change
	void save_constraints();

	//constraints paths/loops; position, normal, color:
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > constraints_tristrip;
	GLVertexArray constraints_tristrip_for_path_draw;
	bool constraints_tristrip_dirty = true;
	void update_constraints_tristrip();

	//-------------------------------
	//interpolation:
	std::vector< float > times;
	void clear_times();
	bool times_dirty = true;
	void update_times();

	//visualization: (constrained model colored with times)
	//constrained model buffer: position, normal, id, texcoord
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4, glm::vec2 > times_model_triangles;
	GLVertexArray times_model_triangles_for_textured_draw;
	bool times_model_triangles_dirty = true;
	void update_times_model_triangles();


	//-------------------------------
	//peeling:
	uint32_t peel_step = 0;
	enum {
		PeelBegin = 0,
		PeelSlice = 1,
		PeelLink = 2,
		PeelBuild = 3,
		PeelRepeat = 4,
	} peel_action = PeelBegin;

	ak::RowColGraph rowcol_graph;
	
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > rowcol_graph_tristrip;
	GLVertexArray rowcol_graph_tristrip_for_path_draw;
	bool rowcol_graph_tristrip_dirty = true;
	void update_rowcol_graph_tristrip();


	// - - - - - - - - - - - - - - - - 
	//peeling - begin / step:
	std::vector< std::vector< ak::EmbeddedVertex > > active_chains;
	std::vector< std::vector< ak::Stitch > > active_stitches;

	//active chains + stitches; position, normal, color:
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > active_chains_tristrip;
	GLVertexArray active_chains_tristrip_for_path_draw;
	bool active_chains_tristrip_dirty = true;
	void update_active_chains_tristrip();

	// - - - - - - - - - - - - - - - - 
	//peeling - slice:
	ak::Model slice;
	std::vector< ak::EmbeddedVertex > slice_on_model;
	std::vector< std::vector< uint32_t > > slice_active_chains;
	std::vector< std::vector< uint32_t > > slice_next_chains;
	std::vector< float > slice_times;

	//sliced model: position, normal, color
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > slice_triangles;
	GLVertexArray slice_triangles_for_path_draw;
	bool slice_triangles_dirty = true;
	void update_slice_triangles();

	//chains on slice; position, normal, color:
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > slice_chains_tristrip;
	GLVertexArray slice_chains_tristrip_for_path_draw;
	bool slice_chains_tristrip_dirty = true;
	void update_slice_chains_tristrip();

	// - - - - - - - - - - - - - - - - 
	//peeling - link:
	std::vector< std::vector< ak::Stitch > > next_stitches;
	std::vector< ak::Link > links;

	//links (+ stitches); position, normal, color:
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > links_tristrip;
	GLVertexArray links_tristrip_for_path_draw;
	bool links_tristrip_dirty = true;
	void update_links_tristrip();

	// - - - - - - - - - - - - - - - - 
	//peeling - build:
	std::vector< std::vector< ak::EmbeddedVertex > > next_active_chains;
	std::vector< std::vector< ak::Stitch > > next_active_stitches;

	//active chains + stitches; position, normal, color:
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > next_active_chains_tristrip;
	GLVertexArray next_active_chains_tristrip_for_path_draw;
	bool next_active_chains_tristrip_dirty = true;
	void update_next_active_chains_tristrip();

	// - - - - - - - - - - - - - - - - 
	//driver functions that step through the above:
	void clear_peeling();
	bool step_peeling();

	//-------------------------------
	//tracing:

	std::vector< ak::TracedStitch > traced;
	void clear_traced();
	bool traced_dirty = true;
	void update_traced();


	std::string save_traced_file = ""; //if not "", will save traced stitches to this file after every change
	void save_traced();


	bool traced_tristrip_dirty = true;
	//traced yarns + stitches: position, normal, color:
	GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 > traced_tristrip;
	GLVertexArray traced_tristrip_for_path_draw;
	void update_traced_tristrip();


	//link bottom constraints directly to top constraints (to test linking function)
	void DEBUG_test_linking(bool flip = false);


	//place camera to view whole model:
	void reset_camera();

};
