#include "Interface.hpp"

#include <kit/GLProgram.hpp>
#include <kit/Load.hpp>
#include <kit/gl_errors.hpp>
#include <kit/check_fb.hpp>

GLuint model_draw_p2c = -1U;
GLuint model_draw_p2l = -1U;
GLuint model_draw_n2l = -1U;

kit::Load< GLProgram > model_draw(kit::LoadTagDefault, [](){
	GLProgram *ret = new GLProgram(
		"#version 330\n"
		"#line " STR(__LINE__) "\n"
		"uniform mat4 p2c;\n" //position to clip space
		"uniform mat4x3 p2l;\n" //position to light space
		"uniform mat3 n2l;\n" //normal to light space
		"in vec4 Position;\n"
		"in vec3 Normal;\n"
		"in vec4 ID;\n"
		"out vec3 normal;\n"
		"out vec3 position;\n"
		"out vec4 id;\n"
		"void main() {\n"
		"	gl_Position = p2c * Position;\n"
		"	position = p2l * Position;\n"
		"	normal = n2l * Normal;\n"
		"	id = ID;\n"
		"}\n"
		,
		"#version 330\n"
		"#line " STR(__LINE__) "\n"
		"in vec3 normal;\n"
		"in vec3 position;\n"
		"in vec4 id;\n"
		"layout(location = 0) out vec4 fragColor;\n"
		"layout(location = 1) out vec4 fragID;\n"
		"void main() {\n"
		"	vec3 n = normalize(normal);\n"
		"	vec3 l = vec3(0.0, 0.0, 1.0);\n"
		"	float nl = dot(n, l) * 0.5 + 0.5;\n"
		"	fragColor = vec4(nl * vec3(1.0), 1.0);\n"
		"	fragID = id;\n"
		"}\n"
	);

	model_draw_p2c = ret->getUniformLocation("p2c", GLProgram::MissingIsError);
	model_draw_p2l = ret->getUniformLocation("p2l", GLProgram::MissingIsWarning);
	model_draw_n2l = ret->getUniformLocation("n2l", GLProgram::MissingIsWarning);

	return ret;
});

//------------------------------------

GLuint marker_draw_p2c = -1U;
GLuint marker_draw_p2l = -1U;
GLuint marker_draw_n2l = -1U;

GLuint marker_draw_id = -1U;
GLuint marker_draw_color = -1U;

kit::Load< GLProgram > marker_draw(kit::LoadTagDefault, [](){
	GLProgram *ret = new GLProgram(
		"#version 330\n"
		"#line " STR(__LINE__) "\n"
		"uniform mat4 p2c;\n" //position to clip space
		"uniform mat4x3 p2l;\n" //position to light space
		"uniform mat3 n2l;\n" //normal to light space
		"in vec4 Position;\n"
		"in vec3 Normal;\n"
		"out vec3 normal;\n"
		"out vec3 position;\n"
		"void main() {\n"
		"	gl_Position = p2c * Position;\n"
		"	position = p2l * Position;\n"
		"	normal = n2l * Normal;\n"
		"}\n"
		,
		"#version 330\n"
		"#line " STR(__LINE__) "\n"
		"uniform vec3 color;\n"
		"uniform vec4 id;\n"
		"in vec3 normal;\n"
		"in vec3 position;\n"
		"layout(location = 0) out vec4 fragColor;\n"
		"layout(location = 1) out vec4 fragID;\n"
		"void main() {\n"
		"	vec3 n = normalize(normal);\n"
		"	vec3 l = vec3(0.0, 0.0, 1.0);\n"
		"	float nl = dot(n, l) * 0.5 + 0.5;\n"
		"	fragColor = vec4(nl * color, 1.0);\n"
		"	fragID = id;\n"
		"}\n"
	);

	marker_draw_p2c = ret->getUniformLocation("p2c", GLProgram::MissingIsError);
	marker_draw_p2l = ret->getUniformLocation("p2l", GLProgram::MissingIsWarning);
	marker_draw_n2l = ret->getUniformLocation("n2l", GLProgram::MissingIsWarning);

	marker_draw_color = ret->getUniformLocation("color", GLProgram::MissingIsWarning);
	marker_draw_id = ret->getUniformLocation("id", GLProgram::MissingIsWarning);

	return ret;
});

//------------------------------------


kit::Load< GLProgram > copy_fb(kit::LoadTagDefault, [](){
	GLProgram *ret = new GLProgram(
		//ye olde "fullscreen triangle":
		"#version 330\n"
		"#line " STR(__LINE__) "\n"
		"void main() {\n"
		"	gl_Position = vec4(4 * (gl_VertexID & 1) - 1,  2 * (gl_VertexID & 2) - 1, 0.0, 1.0);\n"
		"}\n"
		,
		//infinite perspective matrix takes z to (-z -2n) / -z = 1 + 2n/z (so z ranges from -n to -inf go to -1 to 1
		//depth-component texture takes [-1,1] to [0,1]
		"#version 330\n"
		"#line " STR(__LINE__) "\n"
		"uniform sampler2D color_tex;\n"
		"uniform sampler2D depth_tex;\n"
		"layout(location = 0) out vec4 fragColor;\n"
		"void main() {\n"
		"	ivec2 px = ivec2(gl_FragCoord.xy);\n"
		"	float NEAR = 0.1;\n"
		"	float ref = texelFetchOffset(depth_tex, px, 0, ivec2(0,0)).r;\n"
		"	float z = min(0.9999, ref);\n"
		//over from depth component space -> projected z value:
		"	z = 2.0 * (z - 0.5);\n"
		//from projected z value to world z value:
		"	z = (2.0 * NEAR) / (z - 1.0);\n"
		//z is now in range [-inf, -NEAR]
		"	z += 0.5;\n"
		"	z = 1.0 + (2.0 * NEAR) / z;\n"
		"	z = 0.5 * z + 0.5;\n"
		"	float amt = 1.0 / (z - ref);\n"
		"	#define DO(X,Y) { \\\n"
		"		float val = texelFetchOffset(depth_tex, px, 0, ivec2(X,Y)).r; \\\n"
		"		tint += clamp((val - ref) * amt, 0.0, 1.0);\\\n"
		"	}\n"
		"	float tint = 0.0;\n"
		"	           DO(-1,  2) DO( 0,  2) DO( 1,  2)          \n"
		"	DO(-2,  1) DO(-1,  1) DO( 0,  1) DO( 1,  1) DO(2,  1)\n"
		"	DO(-2,  0) DO(-1,  0)            DO( 1,  0) DO(2,  0)\n"
		"	DO(-2, -1) DO(-1, -1) DO( 0, -1) DO( 1, -1) DO(2, -1)\n"
		"	           DO(-1, -2) DO( 0, -2) DO( 1, -2)          \n"
		"	tint = 1.0 - (tint / 20.0);\n"
		"	fragColor = vec4(tint, tint, tint, 1.0) * texelFetch(color_tex, px, 0);\n"
		"}\n"
	);

	glUseProgram(ret->program);
	glUniform1i(ret->getUniformLocation("color_tex", GLProgram::MissingIsWarning), 0);
	glUniform1i(ret->getUniformLocation("depth_tex", GLProgram::MissingIsWarning), 1);
	glUseProgram(0);

	return ret;
});

//------------------------------------

//layout: (position)
kit::Load< GLAttribBuffer< glm::vec3 > > sphere_tristrip(kit::LoadTagInit, [](){
	GLAttribBuffer< glm::vec3 > *ret = new GLAttribBuffer< glm::vec3 >();
	const uint32_t Rings = 10;
	const uint32_t Slices = 16;

	std::vector< glm::vec2 > circle;
	circle.reserve(Slices + 1);
	for (uint32_t a = 0; a < Slices; ++a) {
		float ang = a / float(Slices) * M_PI * 2.0f;
		circle.emplace_back(std::cos(ang), std::sin(ang));
	}
	circle.emplace_back(circle[0]);

	std::vector< GLAttribBuffer< glm::vec3 >::Vertex > attribs;
	attribs.reserve((Rings) * 2 * (Slices+1));

	glm::vec2 prev(0.0f, -1.0f);
	for (uint32_t a = 1; a <= Rings; ++a) {
		glm::vec2 next;
		if (a == Rings) {
			next = glm::vec2(0.0f, 1.0f);
		} else {
			float ang = (a / float(Rings) - 0.5f) * M_PI;
			next = glm::vec2(std::cos(ang), std::sin(ang));
		}

		for (auto const &c : circle) {
			attribs.emplace_back(glm::vec3(prev.x * c.x, prev.x * c.y, prev.y));
			attribs.emplace_back(glm::vec3(next.x * c.x, next.x * c.y, next.y));
		}

		prev = next;
	}

	ret->set(attribs, GL_STATIC_DRAW);

	return ret;
});

kit::Load< GLVertexArray > sphere_tristrip_for_marker_draw(kit::LoadTagDefault, [](){
	return new GLVertexArray(GLVertexArray::make_binding(marker_draw->program, {
		{marker_draw->getAttribLocation("Position", GLProgram::MissingIsError), (*sphere_tristrip)[0]},
		{marker_draw->getAttribLocation("Normal", GLProgram::MissingIsError), (*sphere_tristrip)[0]}
	}));
});

kit::Load< GLVertexArray > empty_vertex_array(kit::LoadTagDefault);


Interface::Interface() {
	model_triangles_for_model_draw = GLVertexArray::make_binding(model_draw->program, {
		{model_draw->getAttribLocation("Position", GLProgram::MissingIsError), model_triangles[0]},
		{model_draw->getAttribLocation("Normal", GLProgram::MissingIsWarning), model_triangles[1]},
		{model_draw->getAttribLocation("ID", GLProgram::MissingIsWarning), model_triangles[2]}
	});
}

Interface::~Interface() {
	if (color_id_fb) {
		glDeleteFramebuffers(1, &color_id_fb);
		color_id_fb = 0;
	}
	if (color_tex) {
		glDeleteTextures(1, &color_tex);
	}
	if (id_tex) {
		glDeleteTextures(1, &id_tex);
	}
	if (depth_tex) {
		glDeleteTextures(1, &depth_tex);
	}
}

void Interface::update(float elapsed) {
	if (mouse.moved) {
		update_hovered();
		mouse.moved = false;
	}
}

void Interface::draw() {
	if (fb_size != kit::display.size) alloc_fbs();

	glViewport(0, 0, kit::display.size.x, kit::display.size.y);

	glBindFramebuffer(GL_FRAMEBUFFER, color_id_fb);

	{
		glClearColor(0.9, 0.9, 0.9, 0.0);
		glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

		GLfloat zeros[4] = {0.0f, 0.0f, 0.0f, 0.0f};
		glClearBufferfv(GL_COLOR, 1, zeros);
	}

	glEnable(GL_DEPTH_TEST);


	{ //draw the model:
		glUseProgram(model_draw->program);

		//Position-to-clip matrix:
		glm::mat4 p2c = camera.mvp();
		//Position-to-light matrix:
		glm::mat4x3 p2l = glm::mat4(1.0f);
		//Normal-to-light matrix:
		glm::mat3 n2l = glm::inverse(glm::transpose(glm::mat3(p2l)));

		glUniformMatrix4fv(model_draw_p2c, 1, GL_FALSE, glm::value_ptr(p2c));
		glUniformMatrix4x3fv(model_draw_p2l, 1, GL_FALSE, glm::value_ptr(p2l));
		glUniformMatrix3fv(model_draw_n2l, 1, GL_FALSE, glm::value_ptr(n2l));

		glBindVertexArray(model_triangles_for_model_draw.array);

		glDrawArrays(GL_TRIANGLES, 0, model_triangles.count);

		glBindVertexArray(0);

		glUseProgram(0);
	}

	{ //draw marker spheres:
		glUseProgram(marker_draw->program);
		glBindVertexArray(sphere_tristrip_for_marker_draw->array);

		auto sphere = [&](glm::vec3 const &at, float r, glm::vec3 const &color, glm::u8vec4 const &id) {
			//Position-to-light matrix:
			glm::mat4x3 p2l = glm::mat4x3(
				glm::vec3(r, 0.0f, 0.0f),
				glm::vec3(0.0f, r, 0.0f),
				glm::vec3(0.0f, 0.0f, r),
				glm::vec3(at)
			);
			//Normal-to-light matrix:
			glm::mat3 n2l = glm::inverse(glm::transpose(glm::mat3(p2l)));
			//Position-to-clip matrix:
			glm::mat4 p2c = camera.mvp() * glm::mat4(p2l);

			glUniformMatrix4fv(marker_draw_p2c, 1, GL_FALSE, glm::value_ptr(p2c));
			glUniformMatrix4x3fv(marker_draw_p2l, 1, GL_FALSE, glm::value_ptr(p2l));
			glUniformMatrix3fv(marker_draw_n2l, 1, GL_FALSE, glm::value_ptr(n2l));

			glUniform4f(marker_draw_id, id.r / 255.0f, id.g / 255.0f, id.b / 255.0f, id.a / 255.0f);
			glUniform3fv(marker_draw_color, 1, glm::value_ptr(color));

			glDrawArrays(GL_TRIANGLE_STRIP, 0, sphere_tristrip->count);
		};

		//mouse cursor sphere (don't draw into ID buffer):
		glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
		if (hovered.tri < model.triangles.size()) {
			glm::vec3 const &a = model.vertices[model.triangles[hovered.tri].x];
			glm::vec3 const &b = model.vertices[model.triangles[hovered.tri].y];
			glm::vec3 const &c = model.vertices[model.triangles[hovered.tri].z];
			glm::vec3 at = (a + b + c) / 3.0f;
			sphere(at, 0.2f, glm::vec3(1.0f, 0.0f, 0.0f), glm::u8vec4(0x00));
		}
		glColorMaski(1, GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

		glBindVertexArray(0);
		glUseProgram(0);
	}

	GL_ERRORS();

	//---------------------------------------------------
	//copy rendered color to screen:

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	glClearColor(1.0f, 0.0f, 1.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	glUseProgram(copy_fb->program);

	glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, color_tex);
	glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, depth_tex);

	glBindVertexArray(empty_vertex_array->array);
	glDrawArrays(GL_TRIANGLES, 0, 3);
	glBindVertexArray(0);

	glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, 0);

	glUseProgram(0);

	GL_ERRORS();

}

void Interface::pointer_action(kit::PointerID pointer, kit::PointerAction action, kit::Pointer const &old_state, kit::Pointer const &new_state) {
	if (mouse.at != new_state.at) {
		mouse.at = new_state.at;
		mouse.moved = true;
	}
	
	if (action == kit::PointerDown) {
		uint8_t pressed = ~old_state.buttons & new_state.buttons;
		if (pressed == kit::ButtonRight) {
			if (drag == DragNone) {
				if (std::cos(camera.elevation) > 0.0f) {
					drag = DragCamera;
				} else {
					drag = DragCameraFlipX;
				}
			}
		}
	} else if (action == kit::PointerUp) {
		uint8_t released = old_state.buttons & ~new_state.buttons;
		if (released == kit::ButtonRight) {
			if (drag == DragCamera || drag == DragCameraFlipX) {
				drag = DragNone;
			}
		}
	} else if (action == kit::PointerMove) {
		if (drag == DragCamera || drag == DragCameraFlipX) {
			glm::vec2 d = new_state.at - old_state.at;
			if (drag == DragCameraFlipX) d.x = -d.x;

			camera.azimuth -= d.x;
			camera.elevation -= d.y;
		}
	}
}

void Interface::update_hovered() {
	hovered.clear();

	glm::vec2 px(0.5f * (mouse.at.x + 1.0f) * fb_size.x, 0.5f * (mouse.at.y + 1.0f) * fb_size.y);

	if (px.x < 0 || px.x >= fb_size.x || px.y < 0 || px.y >= fb_size.y) return;

	glBindFramebuffer(GL_FRAMEBUFFER, color_id_fb);
	glReadBuffer(GL_COLOR_ATTACHMENT1);
	glm::u8vec4 col;
	glReadPixels(px.x, px.y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, &col);

	std::cout << int(col.r) << " " << int(col.g) << " " << int(col.b) << " " << int(col.a) << std::endl; //DEBUG

	if (col.r == 1) {
		uint32_t idx = uint32_t(col.a) | (uint32_t(col.b) << 8) | (uint32_t(col.g) << 16);
		if (idx < model.triangles.size()) {
			hovered.tri = idx;
		}
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void Interface::alloc_fbs() {
	fb_size = kit::display.size;

	//allocate textures:

	if (color_tex == 0) glGenTextures(1, &color_tex);
	glBindTexture(GL_TEXTURE_2D, color_tex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, fb_size.x, fb_size.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_2D, 0);

	if (id_tex == 0) glGenTextures(1, &id_tex);
	glBindTexture(GL_TEXTURE_2D, id_tex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, fb_size.x, fb_size.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_2D, 0);

	if (depth_tex == 0) glGenTextures(1, &depth_tex);
	glBindTexture(GL_TEXTURE_2D, depth_tex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, fb_size.x, fb_size.y, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_2D, 0);

	GL_ERRORS();

	//set up framebuffer:

	if (color_id_fb == 0) glGenFramebuffers(1, &color_id_fb);
	glBindFramebuffer(GL_FRAMEBUFFER, color_id_fb);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, color_tex, 0);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, id_tex, 0);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depth_tex, 0);
	GLenum bufs[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
	glDrawBuffers(2, bufs);

	check_fb();
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	GL_ERRORS();
}

void Interface::set_model(ak::Model const &new_model) {
	model = new_model;

	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > attribs;
	attribs.reserve(3 * model.triangles.size());

	for (auto const &t : model.triangles) {
		glm::vec3 const &a = model.vertices[t.x];
		glm::vec3 const &b = model.vertices[t.y];
		glm::vec3 const &c = model.vertices[t.z];

		glm::vec3 n = glm::normalize(glm::cross(b-a, c-a));

		uint32_t idx = &t - &model.triangles[0];
		glm::u8vec4 id = glm::u8vec4(
			0x1,
			(idx >> 16) & 0xff,
			(idx >> 8) & 0xff,
			idx & 0xff
		);

		attribs.emplace_back(a, n, id);
		attribs.emplace_back(b, n, id);
		attribs.emplace_back(c, n, id);
	}

	model_triangles.set(attribs, GL_STATIC_DRAW);

	reset_camera();
}

void Interface::reset_camera() {
	glm::vec3 min = glm::vec3(std::numeric_limits< float >::infinity());
	glm::vec3 max = glm::vec3(-std::numeric_limits< float >::infinity());
	for (auto const &t : model.triangles) {
		for (uint32_t i = 0; i < 3; ++i) {
			min = glm::min(min, model.vertices[t[i]]);
			max = glm::max(max, model.vertices[t[i]]);
		}
	}
	camera.center = 0.5f * (max + min);
	camera.radius = glm::length(max - min);
	camera.azimuth = glm::radians(30.0f);
	camera.elevation = glm::radians(45.0f);
}
