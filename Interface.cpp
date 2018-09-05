#include "Interface.hpp"

#include "Stitch.hpp"

#include <kit/GLProgram.hpp>
#include <kit/GLTexture.hpp>
#include <kit/Load.hpp>
#include <kit/gl_errors.hpp>
#include <kit/check_fb.hpp>

#include <glm/gtx/norm.hpp>
#include <glm/gtx/hash.hpp>

#include <SDL.h>

#include <algorithm>
#include <unordered_set>

//given normalized 0..1 time:
glm::vec3 time_color(float time) {
	const constexpr size_t Size = 3;
	static glm::vec3 grad[Size] = {
		glm::vec3(0.2f, 0.2f, 1.0f),
		glm::vec3(0.8f, 0.8f, 0.8f),
		glm::vec3(0.8f, 0.2f, 0.2f)
	};
	time *= (Size-1);
	int32_t i = std::max(0, std::min(int32_t(Size)-2, int32_t(std::floor(time))));
	float f = std::max(0.0f, std::min(1.0f, time - i));
	return glm::mix(grad[i], grad[i+1], f);
}

//-------------------------------------------------------------

constexpr const uint32_t TimeTexSize = 16;

kit::Load< GLTexture > time_tex(kit::LoadTagDefault, [](){
	GLTexture *ret = new GLTexture();
	std::vector< glm::vec3 > data(TimeTexSize);
	for (uint32_t i = 0; i < data.size(); ++i) {
		data[i] = time_color(i / float(TimeTexSize-1));
	}

	glBindTexture(GL_TEXTURE_2D, ret->texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, data.size(), 1, 0, GL_RGB, GL_FLOAT, data.data());
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_2D, 0);

	return ret;
});

//-------------------------------------------------------------

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

//-------------------------------------------------------------

GLuint textured_draw_p2c = -1U;
GLuint textured_draw_p2l = -1U;
GLuint textured_draw_n2l = -1U;

kit::Load< GLProgram > textured_draw(kit::LoadTagDefault, [](){
	GLProgram *ret = new GLProgram(
		"#version 330\n"
		"#line " STR(__LINE__) "\n"
		"uniform mat4 p2c;\n" //position to clip space
		"uniform mat4x3 p2l;\n" //position to light space
		"uniform mat3 n2l;\n" //normal to light space
		"in vec4 Position;\n"
		"in vec3 Normal;\n"
		"in vec4 ID;\n"
		"in vec2 TexCoord;\n"
		"out vec3 normal;\n"
		"out vec3 position;\n"
		"out vec4 id;\n"
		"out vec2 texCoord;\n"
		"void main() {\n"
		"	gl_Position = p2c * Position;\n"
		"	position = p2l * Position;\n"
		"	normal = n2l * Normal;\n"
		"	id = ID;\n"
		"	texCoord = TexCoord;\n"
		"}\n"
		,
		"#version 330\n"
		"#line " STR(__LINE__) "\n"
		"uniform sampler2D tex;\n"
		"in vec3 normal;\n"
		"in vec3 position;\n"
		"in vec4 id;\n"
		"in vec2 texCoord;\n"
		"layout(location = 0) out vec4 fragColor;\n"
		"layout(location = 1) out vec4 fragID;\n"
		"void main() {\n"
		"	vec3 n = normalize(normal);\n"
		"	vec3 l = vec3(0.0, 0.0, 1.0);\n"
		"	float nl = dot(n, l) * 0.5 + 0.5;\n"
		"	fragColor = vec4(nl * texture(tex, texCoord).rgb, 1.0);\n"
		"	fragID = id;\n"
		"}\n"
	);

	textured_draw_p2c = ret->getUniformLocation("p2c", GLProgram::MissingIsError);
	textured_draw_p2l = ret->getUniformLocation("p2l", GLProgram::MissingIsWarning);
	textured_draw_n2l = ret->getUniformLocation("n2l", GLProgram::MissingIsWarning);

	glUseProgram(ret->program);
	glUniform1i(ret->getUniformLocation("tex", GLProgram::MissingIsWarning), 0);
	glUseProgram(0);

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

GLuint path_draw_p2c = -1U;
GLuint path_draw_p2l = -1U;
GLuint path_draw_n2l = -1U;

GLuint path_draw_id = -1U;

kit::Load< GLProgram > path_draw(kit::LoadTagDefault, [](){
	GLProgram *ret = new GLProgram(
		"#version 330\n"
		"#line " STR(__LINE__) "\n"
		"uniform mat4 p2c;\n" //position to clip space
		"uniform mat4x3 p2l;\n" //position to light space
		"uniform mat3 n2l;\n" //normal to light space
		"in vec4 Position;\n"
		"in vec3 Normal;\n"
		"in vec3 Color;\n"
		"out vec3 normal;\n"
		"out vec3 position;\n"
		"out vec3 color;\n"
		"void main() {\n"
		"	gl_Position = p2c * Position;\n"
		"	position = p2l * Position;\n"
		"	normal = n2l * Normal;\n"
		"	color = Color;\n"
		"}\n"
		,
		"#version 330\n"
		"#line " STR(__LINE__) "\n"
		"uniform vec4 id;\n"
		"in vec3 normal;\n"
		"in vec3 position;\n"
		"in vec3 color;\n"
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

	path_draw_p2c = ret->getUniformLocation("p2c", GLProgram::MissingIsError);
	path_draw_p2l = ret->getUniformLocation("p2l", GLProgram::MissingIsWarning);
	path_draw_n2l = ret->getUniformLocation("n2l", GLProgram::MissingIsWarning);

	path_draw_id = ret->getUniformLocation("id", GLProgram::MissingIsWarning);

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
	std::cout << "Setting up various buffer bindings." << std::endl; //DEBUG

	model_triangles_for_model_draw = GLVertexArray::make_binding(model_draw->program, {
		{model_draw->getAttribLocation("Position", GLProgram::MissingIsError), model_triangles[0]},
		{model_draw->getAttribLocation("Normal", GLProgram::MissingIsWarning), model_triangles[1]},
		{model_draw->getAttribLocation("ID", GLProgram::MissingIsWarning), model_triangles[2]}
	});

	constraints_tristrip_for_path_draw = GLVertexArray::make_binding(path_draw->program, {
		{path_draw->getAttribLocation("Position", GLProgram::MissingIsError), constraints_tristrip[0]},
		{path_draw->getAttribLocation("Normal", GLProgram::MissingIsWarning), constraints_tristrip[1]},
		{path_draw->getAttribLocation("Color", GLProgram::MissingIsWarning), constraints_tristrip[2]}
	});

	times_model_triangles_for_textured_draw = GLVertexArray::make_binding(textured_draw->program, {
		{textured_draw->getAttribLocation("Position", GLProgram::MissingIsError), times_model_triangles[0]},
		{textured_draw->getAttribLocation("Normal", GLProgram::MissingIsWarning), times_model_triangles[1]},
		{textured_draw->getAttribLocation("ID", GLProgram::MissingIsWarning), times_model_triangles[2]},
		{textured_draw->getAttribLocation("TexCoord", GLProgram::MissingIsWarning), times_model_triangles[3]}
	});

	rowcol_graph_tristrip_for_path_draw = GLVertexArray::make_binding(path_draw->program, {
		{path_draw->getAttribLocation("Position", GLProgram::MissingIsError), rowcol_graph_tristrip[0]},
		{path_draw->getAttribLocation("Normal", GLProgram::MissingIsWarning), rowcol_graph_tristrip[1]},
		{path_draw->getAttribLocation("Color", GLProgram::MissingIsWarning), rowcol_graph_tristrip[2]}
	});

	active_chains_tristrip_for_path_draw = GLVertexArray::make_binding(path_draw->program, {
		{path_draw->getAttribLocation("Position", GLProgram::MissingIsError), active_chains_tristrip[0]},
		{path_draw->getAttribLocation("Normal", GLProgram::MissingIsWarning), active_chains_tristrip[1]},
		{path_draw->getAttribLocation("Color", GLProgram::MissingIsWarning), active_chains_tristrip[2]}
	});

	slice_triangles_for_path_draw = GLVertexArray::make_binding(path_draw->program, {
		{path_draw->getAttribLocation("Position", GLProgram::MissingIsError), slice_triangles[0]},
		{path_draw->getAttribLocation("Normal", GLProgram::MissingIsWarning), slice_triangles[1]},
		{path_draw->getAttribLocation("Color", GLProgram::MissingIsWarning), slice_triangles[2]},
	});

	slice_chains_tristrip_for_path_draw = GLVertexArray::make_binding(path_draw->program, {
		{path_draw->getAttribLocation("Position", GLProgram::MissingIsError), slice_chains_tristrip[0]},
		{path_draw->getAttribLocation("Normal", GLProgram::MissingIsWarning), slice_chains_tristrip[1]},
		{path_draw->getAttribLocation("Color", GLProgram::MissingIsWarning), slice_chains_tristrip[2]}
	});

	links_tristrip_for_path_draw = GLVertexArray::make_binding(path_draw->program, {
		{path_draw->getAttribLocation("Position", GLProgram::MissingIsError), links_tristrip[0]},
		{path_draw->getAttribLocation("Normal", GLProgram::MissingIsWarning), links_tristrip[1]},
		{path_draw->getAttribLocation("Color", GLProgram::MissingIsWarning), links_tristrip[2]}
	});

	next_active_chains_tristrip_for_path_draw = GLVertexArray::make_binding(path_draw->program, {
		{path_draw->getAttribLocation("Position", GLProgram::MissingIsError), next_active_chains_tristrip[0]},
		{path_draw->getAttribLocation("Normal", GLProgram::MissingIsWarning), next_active_chains_tristrip[1]},
		{path_draw->getAttribLocation("Color", GLProgram::MissingIsWarning), next_active_chains_tristrip[2]}
	});

	traced_tristrip_for_path_draw = GLVertexArray::make_binding(path_draw->program, {
		{path_draw->getAttribLocation("Position", GLProgram::MissingIsError), traced_tristrip[0]},
		{path_draw->getAttribLocation("Normal", GLProgram::MissingIsWarning), traced_tristrip[1]},
		{path_draw->getAttribLocation("Color", GLProgram::MissingIsWarning), traced_tristrip[2]}
	});



	GL_ERRORS();

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
	if (constraints_dirty) {
		update_constraints();
	}
}

void Interface::draw() {
	if (fb_size != kit::display.size) alloc_fbs();

	glViewport(0, 0, kit::display.size.x, kit::display.size.y);

	glBindFramebuffer(GL_FRAMEBUFFER, color_id_fb);

	{
		glClearColor(0.9f, 0.9f, 0.9f, 0.0f);
		glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

		GLfloat zeros[4] = {0.0f, 0.0f, 0.0f, 0.0f};
		glClearBufferfv(GL_COLOR, 1, zeros);
	}

	glEnable(GL_DEPTH_TEST);

	if (show & ShowModel) { //draw the model:
		if (model_triangles_dirty) update_model_triangles();
		glUseProgram(model_draw->program);

		//Position-to-clip matrix:
		glm::mat4 p2c = camera.mvp();
		//Position-to-light matrix:
		glm::mat4x3 p2l = camera.mv(); //glm::mat4(1.0f);
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

	if (show & ShowTimesModel) { //draw the constrained model:
		if (times_model_triangles_dirty) update_times_model_triangles();
		glUseProgram(textured_draw->program);

		//Position-to-clip matrix:
		glm::mat4 p2c = camera.mvp();
		//Position-to-light matrix:
		glm::mat4x3 p2l = camera.mv(); //glm::mat4(1.0f);
		//Normal-to-light matrix:
		glm::mat3 n2l = glm::inverse(glm::transpose(glm::mat3(p2l)));

		glUniformMatrix4fv(textured_draw_p2c, 1, GL_FALSE, glm::value_ptr(p2c));
		glUniformMatrix4x3fv(textured_draw_p2l, 1, GL_FALSE, glm::value_ptr(p2l));
		glUniformMatrix3fv(textured_draw_n2l, 1, GL_FALSE, glm::value_ptr(n2l));

		glBindVertexArray(times_model_triangles_for_textured_draw.array);
		glBindTexture(GL_TEXTURE_2D, time_tex->texture);

		glDrawArrays(GL_TRIANGLES, 0, times_model_triangles.count);

		glBindTexture(GL_TEXTURE_2D, 0);
		glBindVertexArray(0);

		glUseProgram(0);
	}

	if (show & ShowSlice) { //draw the clipped model debug out:
		if (slice_triangles_dirty) update_slice_triangles();
		glUseProgram(path_draw->program);

		//Position-to-clip matrix:
		glm::mat4 p2c = camera.mvp();
		//Position-to-light matrix:
		glm::mat4x3 p2l = camera.mv(); //glm::mat4(1.0f);
		//Normal-to-light matrix:
		glm::mat3 n2l = glm::inverse(glm::transpose(glm::mat3(p2l)));

		glUniformMatrix4fv(path_draw_p2c, 1, GL_FALSE, glm::value_ptr(p2c));
		glUniformMatrix4x3fv(path_draw_p2l, 1, GL_FALSE, glm::value_ptr(p2l));
		glUniformMatrix3fv(path_draw_n2l, 1, GL_FALSE, glm::value_ptr(n2l));

		glBindVertexArray(slice_triangles_for_path_draw.array);

		glDrawArrays(GL_TRIANGLES, 0, slice_triangles.count);

		glBindVertexArray(0);

		glUseProgram(0);
	}


	if (show & ShowModel) { //draw marker spheres:
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

		//constrained vertices (do draw into ID buffer):
		float min_time = -1.0f;
		float max_time = 1.0f;
		for (auto const &c : constraints) {
			min_time = glm::min(min_time, c.value);
			max_time = glm::max(max_time, c.value);
		}
		for (auto const &c : constraints) {
			uint32_t idx = &c - &constraints[0];
			glm::vec3 color = time_color((c.value - min_time) / (max_time - min_time));
			for (uint32_t i = 0; i < c.chain.size(); ++i) {
				glm::u8vec4 id(2, i, (idx >> 8) & 0xff, idx & 0xff);
				sphere(model.vertices[c.chain[i]], 0.04f, color, id);
			}
		}

		//mouse cursor sphere (don't draw into ID buffer):
		glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
		if (hovered.tri < model.triangles.size()) {
			glm::vec3 const &a = model.vertices[model.triangles[hovered.tri].x];
			glm::vec3 const &b = model.vertices[model.triangles[hovered.tri].y];
			glm::vec3 const &c = model.vertices[model.triangles[hovered.tri].z];
			sphere(a, 0.02f, glm::vec3(1.0f, 0.0f, 0.0f), glm::u8vec4(0x00));
			sphere(b, 0.02f, glm::vec3(0.0f, 1.0f, 0.0f), glm::u8vec4(0x00));
			sphere(c, 0.02f, glm::vec3(0.0f, 0.0f, 1.0f), glm::u8vec4(0x00));
			sphere(hovered.point, 0.02f, glm::vec3(0.4f, 0.4f, 0.4f), glm::u8vec4(0x00));
		}

		glColorMaski(1, GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

		glBindVertexArray(0);
		glUseProgram(0);
	}

	//draw constraint paths:
	if (show & ShowConstraints) {
		if (constraints_tristrip_dirty) update_constraints_tristrip();

		glm::mat4 p2c = camera.mvp();
		glm::mat4x3 p2l = camera.mv();
		glm::mat3 n2l = glm::inverse(glm::transpose(glm::mat3(p2l)));

		glUseProgram(path_draw->program);
		glBindVertexArray(constraints_tristrip_for_path_draw.array);

		glUniformMatrix4fv(path_draw_p2c, 1, GL_FALSE, glm::value_ptr(p2c));
		glUniformMatrix4x3fv(path_draw_p2l, 1, GL_FALSE, glm::value_ptr(p2l));
		glUniformMatrix3fv(path_draw_n2l, 1, GL_FALSE, glm::value_ptr(n2l));

		glUniform4f(path_draw_id, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f);

		//don't draw into ID array:
		glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

		glDrawArrays(GL_TRIANGLE_STRIP, 0, constraints_tristrip.count);

		glColorMaski(1, GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

		glBindVertexArray(0);
		glUseProgram(0);
	}

	//draw current row-column graph (peeling):
	if (show & ShowRowColGraph) {
		if (rowcol_graph_tristrip_dirty) update_rowcol_graph_tristrip();

		glm::mat4 p2c = camera.mvp();
		glm::mat4x3 p2l = camera.mv();
		glm::mat3 n2l = glm::inverse(glm::transpose(glm::mat3(p2l)));

		glUseProgram(path_draw->program);
		glBindVertexArray(rowcol_graph_tristrip_for_path_draw.array);

		glUniformMatrix4fv(path_draw_p2c, 1, GL_FALSE, glm::value_ptr(p2c));
		glUniformMatrix4x3fv(path_draw_p2l, 1, GL_FALSE, glm::value_ptr(p2l));
		glUniformMatrix3fv(path_draw_n2l, 1, GL_FALSE, glm::value_ptr(n2l));

		glUniform4f(path_draw_id, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f);

		//don't draw into ID array:
		glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

		glDrawArrays(GL_TRIANGLE_STRIP, 0, rowcol_graph_tristrip.count);

		glColorMaski(1, GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

		glBindVertexArray(0);
		glUseProgram(0);
	}


	//draw active chains (peeling):
	if (show & ShowActiveChains) {
		if (active_chains_tristrip_dirty) update_active_chains_tristrip();

		glm::mat4 p2c = camera.mvp();
		glm::mat4x3 p2l = camera.mv();
		glm::mat3 n2l = glm::inverse(glm::transpose(glm::mat3(p2l)));

		glUseProgram(path_draw->program);
		glBindVertexArray(active_chains_tristrip_for_path_draw.array);

		glUniformMatrix4fv(path_draw_p2c, 1, GL_FALSE, glm::value_ptr(p2c));
		glUniformMatrix4x3fv(path_draw_p2l, 1, GL_FALSE, glm::value_ptr(p2l));
		glUniformMatrix3fv(path_draw_n2l, 1, GL_FALSE, glm::value_ptr(n2l));

		glUniform4f(path_draw_id, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f);

		//don't draw into ID array:
		glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

		glDrawArrays(GL_TRIANGLE_STRIP, 0, active_chains_tristrip.count);

		glColorMaski(1, GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

		glBindVertexArray(0);
		glUseProgram(0);
	}

	//draw slice-relative chains (peeling):
	if (show & ShowSliceChains) {
		if (slice_chains_tristrip_dirty) update_slice_chains_tristrip();

		glm::mat4 p2c = camera.mvp();
		glm::mat4x3 p2l = camera.mv();
		glm::mat3 n2l = glm::inverse(glm::transpose(glm::mat3(p2l)));

		glUseProgram(path_draw->program);
		glBindVertexArray(slice_chains_tristrip_for_path_draw.array);

		glUniformMatrix4fv(path_draw_p2c, 1, GL_FALSE, glm::value_ptr(p2c));
		glUniformMatrix4x3fv(path_draw_p2l, 1, GL_FALSE, glm::value_ptr(p2l));
		glUniformMatrix3fv(path_draw_n2l, 1, GL_FALSE, glm::value_ptr(n2l));

		glUniform4f(path_draw_id, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f);

		//don't draw into ID array:
		glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

		glDrawArrays(GL_TRIANGLE_STRIP, 0, slice_chains_tristrip.count);

		glColorMaski(1, GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

		glBindVertexArray(0);
		glUseProgram(0);
	}


	//draw links + next stitches:
	if (show & ShowLinks) {
		if (links_tristrip_dirty) update_links_tristrip();

		glm::mat4 p2c = camera.mvp();
		glm::mat4x3 p2l = camera.mv();
		glm::mat3 n2l = glm::inverse(glm::transpose(glm::mat3(p2l)));

		glUseProgram(path_draw->program);
		glBindVertexArray(links_tristrip_for_path_draw.array);

		glUniformMatrix4fv(path_draw_p2c, 1, GL_FALSE, glm::value_ptr(p2c));
		glUniformMatrix4x3fv(path_draw_p2l, 1, GL_FALSE, glm::value_ptr(p2l));
		glUniformMatrix3fv(path_draw_n2l, 1, GL_FALSE, glm::value_ptr(n2l));

		glUniform4f(path_draw_id, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f);

		//don't draw into ID array:
		glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

		glDrawArrays(GL_TRIANGLE_STRIP, 0, links_tristrip.count);

		glColorMaski(1, GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

		glBindVertexArray(0);
		glUseProgram(0);
	}

	//draw next active chains:
	if (show & ShowNextActiveChains) {
		if (next_active_chains_tristrip_dirty) update_next_active_chains_tristrip();

		glm::mat4 p2c = camera.mvp();
		glm::mat4x3 p2l = camera.mv();
		glm::mat3 n2l = glm::inverse(glm::transpose(glm::mat3(p2l)));

		glUseProgram(path_draw->program);
		glBindVertexArray(next_active_chains_tristrip_for_path_draw.array);

		glUniformMatrix4fv(path_draw_p2c, 1, GL_FALSE, glm::value_ptr(p2c));
		glUniformMatrix4x3fv(path_draw_p2l, 1, GL_FALSE, glm::value_ptr(p2l));
		glUniformMatrix3fv(path_draw_n2l, 1, GL_FALSE, glm::value_ptr(n2l));

		glUniform4f(path_draw_id, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f);

		//don't draw into ID array:
		glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

		glDrawArrays(GL_TRIANGLE_STRIP, 0, next_active_chains_tristrip.count);

		glColorMaski(1, GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

		glBindVertexArray(0);
		glUseProgram(0);
	}

	//draw current traced stitches (tracing):
	if (show & ShowTraced) {
		if (traced_tristrip_dirty) update_traced_tristrip();

		glm::mat4 p2c = camera.mvp();
		glm::mat4x3 p2l = camera.mv();
		glm::mat3 n2l = glm::inverse(glm::transpose(glm::mat3(p2l)));

		glUseProgram(path_draw->program);
		glBindVertexArray(traced_tristrip_for_path_draw.array);

		glUniformMatrix4fv(path_draw_p2c, 1, GL_FALSE, glm::value_ptr(p2c));
		glUniformMatrix4x3fv(path_draw_p2l, 1, GL_FALSE, glm::value_ptr(p2l));
		glUniformMatrix3fv(path_draw_n2l, 1, GL_FALSE, glm::value_ptr(n2l));

		glUniform4f(path_draw_id, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f, 0 / 255.0f);

		//don't draw into ID array:
		glColorMaski(1, GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

		glDrawArrays(GL_TRIANGLE_STRIP, 0, traced_tristrip.count);

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

void Interface::handle_event(SDL_Event const &evt) {
	#define MAPX( X ) ((((X) + 0.5f) / kit::display.window_size.x) * 2.0f - 1.0f)
	#define MAPY( Y ) ((((Y) + 0.5f) / kit::display.window_size.y) *-2.0f + 1.0f)

	if (evt.type == SDL_MOUSEMOTION) {
		glm::vec2 old_at = mouse.at;
		mouse.at = glm::vec2( MAPX(evt.motion.x), MAPY(evt.motion.y) );
		mouse.moved = true;
		if (drag == DragCamera || drag == DragCameraFlipX) {
			glm::vec2 d = mouse.at - old_at;
			if (drag == DragCameraFlipX) d.x = -d.x;

			camera.azimuth -= d.x;
			camera.elevation -= d.y;
		} else if (drag == DragCameraPan) {
			glm::vec2 d = mouse.at - old_at;
			glm::mat3 frame = glm::transpose(glm::mat3(camera.mv()));
			camera.center -= 0.5f * (d.x * frame[0] + d.y * frame[1]) * camera.radius;
		} else if (drag == DragConsPt) {
			if (dragging.cons < constraints.size() && dragging.cons_pt < constraints[dragging.cons].chain.size()) {
				if (hovered.vert < model.vertices.size()) {
					constraints[dragging.cons].chain[dragging.cons_pt] = hovered.vert;
					constraints_dirty = true; //TODO: only set when vert has changed.
				}
			} else {
				drag = DragNone;
				dragging.clear();
			}
		} else if (drag == DragConsRadius) {
			if (dragging.cons < constraints.size() && dragging.cons_pt < constraints[dragging.cons].chain.size()) {
				if (hovered.tri < model.triangles.size()) {
					constraints[dragging.cons].radius = glm::length(hovered.point - model.vertices[constraints[dragging.cons].chain[dragging.cons_pt]]);
					constraints_dirty = true;
				}
			} else {
				drag = DragNone;
				dragging.clear();
			}
		}
	} else if (evt.type == SDL_MOUSEBUTTONDOWN) {
		mouse.at = glm::vec2( MAPX(evt.button.x), MAPY(evt.button.y) );
		mouse.moved = true;
		if (evt.button.button == SDL_BUTTON_LEFT) {
			if (drag == DragNone) {
				if (hovered.cons < constraints.size() && hovered.cons_pt < constraints[hovered.cons].chain.size()) {
					drag = DragConsPt;
					dragging.cons = hovered.cons;
					dragging.cons_pt = hovered.cons_pt;
				}
			}
		} else if (evt.button.button == SDL_BUTTON_RIGHT) {
			if (drag == DragNone) {
				if (SDL_GetModState() & KMOD_SHIFT) {
					drag = DragCameraPan;
				} else {
					if (std::cos(camera.elevation) > 0.0f) {
						drag = DragCamera;
					} else {
						drag = DragCameraFlipX;
					}
				}
			}
		}
	} else if (evt.type == SDL_MOUSEBUTTONUP) {
		mouse.at = glm::vec2( MAPX(evt.button.x), MAPY(evt.button.y) );
		mouse.moved = true;
		if (evt.button.button == SDL_BUTTON_LEFT) {
			if (drag == DragConsPt || drag == DragConsRadius) {
				drag = DragNone;
				dragging.clear();
			}
		} else if (evt.button.button == SDL_BUTTON_RIGHT) {
			if (drag == DragCamera || drag == DragCameraFlipX || drag == DragCameraPan) {
				drag = DragNone;
			}
		}
	} else if (evt.type == SDL_MOUSEWHEEL) {
		camera.radius *= std::pow(0.5f, evt.wheel.y / 10.0f);
		if (camera.radius < 0.1f) camera.radius = 0.1f;
	} else if (evt.type == SDL_KEYDOWN) {
		if (evt.key.keysym.scancode == SDL_SCANCODE_S) {
			if (show & ShowModel) show = (show & ~ShowModelBits) | ShowTimesModel;
			else if (show & ShowTimesModel) show = (show & ~ShowModelBits) | ShowSlice;
			else /*if (show & ShowSlice)*/ show = (show & ~ShowModelBits) | ShowModel;
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_G) {
			show = show ^ ShowRowColGraph;
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_L) {
			//show = ShowConstrainedModel;
			DEBUG_test_linking();
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_P) {
			if (evt.key.keysym.mod & KMOD_SHIFT) {
				clear_peeling();
			} else {
				step_peeling();
			}
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_T) {
			update_traced();
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_C) {
			if (hovered.cons < constraints.size()) {
				if (drag == DragNone) {
					auto &cons = constraints[hovered.cons];
					cons.chain.insert(cons.chain.begin() + hovered.cons_pt, cons.chain[hovered.cons_pt]);
					drag = DragConsPt;
					dragging.cons = hovered.cons;
					dragging.cons_pt = hovered.cons_pt + 1;
				}
			} else if (hovered.vert < model.vertices.size()) {
				constraints.emplace_back();
				constraints.back().chain.emplace_back(hovered.vert);
				constraints_dirty = true;
			}
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_R) {
			if (hovered.cons < constraints.size()) {
				if (drag == DragNone) {
					if (constraints[hovered.cons].radius == 0.0f) {
						drag = DragConsRadius;
						dragging.cons = hovered.cons;
						dragging.cons_pt = hovered.cons_pt;
					} else {
						constraints[hovered.cons].radius = 0.0f;
						constraints_dirty = true;
					}
				}
			}
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_X) {
			if (hovered.cons < constraints.size() && hovered.cons_pt < constraints[hovered.cons].chain.size()) {
				auto &cons = constraints[hovered.cons];
				cons.chain.erase(
					cons.chain.begin() + hovered.cons_pt
				);
				if (cons.chain.empty()) {
					constraints.erase(constraints.begin() + hovered.cons);
				}
				constraints_dirty = true;
			}
			
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_EQUALS || evt.key.keysym.scancode == SDL_SCANCODE_KP_PLUS) {
			if (hovered.cons < constraints.size()) {
				constraints[hovered.cons].value += 0.1f;
				constraints_dirty = true;
			}
		} else if (evt.key.keysym.scancode == SDL_SCANCODE_MINUS || evt.key.keysym.scancode == SDL_SCANCODE_KP_MINUS) {
			if (hovered.cons < constraints.size()) {
				constraints[hovered.cons].value -= 0.1f;
				constraints_dirty = true;
			}
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

	//std::cout << int(col.r) << " " << int(col.g) << " " << int(col.b) << " " << int(col.a) << std::endl; //DEBUG

	if (col.r == 1) {
		uint32_t idx = uint32_t(col.a) | (uint32_t(col.b) << 8) | (uint32_t(col.g) << 16);
		if (idx < model.triangles.size()) {
			hovered.tri = idx;
			//find barycentric coords that project to mouse ray.
			glm::vec3 const &a = model.vertices[model.triangles[hovered.tri].x];
			glm::vec3 const &b = model.vertices[model.triangles[hovered.tri].y];
			glm::vec3 const &c = model.vertices[model.triangles[hovered.tri].z];
			glm::vec3 norm = glm::normalize(glm::cross(b-a, c-a));

			glm::vec3 from = camera.at();
			glm::vec3 dir = glm::normalize(glm::inverse(glm::mat3(camera.mvp())) * glm::vec3(mouse.at, 1.0f));

			//want t such that t * (dir * norm) = (a - from) * norm
			float dn = glm::dot(dir, norm);
			glm::vec3 pt;
			if (std::abs(dn) < 1e-6) {
				pt = (a + b + c) / 3.0f;
			} else {
				pt = (glm::dot(a - from, norm) / dn) * dir + from;
			}
			//hovered.point = pt;

			float cc = glm::dot(pt-a, glm::cross(b-a,norm));
			float ca = glm::dot(pt-b, glm::cross(c-b,norm));
			float cb = glm::dot(pt-c, glm::cross(a-c,norm));

			float sum = ca + cb + cc;
			if (std::abs(sum) < 1e-6) {
				ca = cb = cc = 1.0f / 3.0f;
			} else {
				ca /= sum;
				cb /= sum;
				cc /= sum;
				//clamp to triangle:
				ca = glm::clamp(ca, 0.0f, 1.0f);
				cb = glm::clamp(cb, 0.0f, 1.0f);
				cc = glm::clamp(cc, 0.0f, 1.0f);
				float r = 1.0f - ca - cb - cc;
				ca += r / 3.0f;
				cb += r / 3.0f;
				cc += r / 3.0f;
			}
			hovered.point = ca * a + cb * b + cc * c;
			hovered.coords = glm::vec3(ca, cb, cc);

			if (ca >= cb && ca >= cc) {
				hovered.vert = model.triangles[hovered.tri].x;
			} else if (cb >= cc) {
				hovered.vert = model.triangles[hovered.tri].y;
			} else {
				hovered.vert = model.triangles[hovered.tri].z;
			}
		}
	} else if (col.r == 2) {
		uint32_t idx = uint32_t(col.a) | (uint32_t(col.b) << 8);
		uint32_t pt = uint32_t(col.g);
		if (idx < constraints.size() && pt < constraints[idx].chain.size()) {
			hovered.cons = idx;
			hovered.cons_pt = pt;
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
	model_triangles_dirty = true;
	set_constraints(std::vector< ak::Constraint >());

	reset_camera();
}


void Interface::clear_constraints() {
	constraints.clear();
	constrained_model.clear();
	constrained_values.clear();
	DEBUG_constraint_paths.clear();
	DEBUG_constraint_loops.clear();

	constraints_dirty = true;
	constraints_tristrip_dirty = true;

	DEBUG_constraint_paths.clear();
	DEBUG_constraint_loops.clear();

	clear_times();
}

void Interface::set_constraints(std::vector< ak::Constraint > const &constraints_) {
	clear_constraints();
	constraints = constraints_;
}

void Interface::update_constraints() {
	constraints_dirty = false;

	save_constraints();

	constrained_model.clear();
	constrained_values.clear();
	DEBUG_constraint_paths.clear();
	DEBUG_constraint_loops.clear();

	ak::embed_constraints(parameters, model, constraints, &constrained_model, &constrained_values, &DEBUG_constraint_paths, &DEBUG_constraint_loops);

	constraints_tristrip_dirty = true;

	clear_times();
}

void Interface::save_constraints() {
	if (save_constraints_file == "") return;
	ak::save_constraints(model, constraints, save_constraints_file);
}

void Interface::clear_times() {
	times.clear();

	times_dirty = true;
	times_model_triangles_dirty = true;

	clear_peeling();
}

void Interface::update_times() {
	if (constraints_dirty) update_constraints();
	times_dirty = false;

	try {
		ak::interpolate_values(constrained_model, constrained_values, &times);
	} catch (std::exception &e) {
		std::cout << "ERROR during interpoation: " << e.what() << std::endl;
		times.clear();
	}

	times_model_triangles_dirty = true;

	clear_peeling();
}

void Interface::clear_peeling() {
	peel_step = 0;
	peel_action = PeelBegin;

	rowcol_graph.clear();
	rowcol_graph_tristrip_dirty = true;

	active_chains.clear();
	active_stitches.clear();
	active_chains_tristrip_dirty = true;

	slice.clear();
	slice_on_model.clear();
	slice_active_chains.clear();
	slice_next_chains.clear();
	slice_next_used_boundary.clear();
	slice_times.clear();
	slice_triangles_dirty = true;
	slice_chains_tristrip_dirty = true;

	next_stitches.clear();
	links.clear();
	links_tristrip_dirty = true;

	next_active_chains.clear();
	next_active_stitches.clear();
	next_active_chains_tristrip_dirty = true;

}

bool Interface::step_peeling() {
	if (times_dirty) update_times();
	if (times.empty()) return false; //can't step if no time info

	if (peel_action == PeelBegin || peel_action == PeelRepeat) {
		auto old_peel_step = peel_step;
		auto old_peel_action = peel_action;
		auto old_next_active_chains = next_active_chains;
		auto old_next_active_stitches = next_active_stitches;
		auto old_rowcol_graph = rowcol_graph;
		auto old_rowcol_graph_tristrip_dirty = rowcol_graph_tristrip_dirty;
		clear_peeling();
		peel_step = old_peel_step;
		peel_action = old_peel_action;
		rowcol_graph = old_rowcol_graph;
		rowcol_graph_tristrip_dirty = old_rowcol_graph_tristrip_dirty;

		if (peel_action == PeelBegin) {
			std::cout << " -- peel begin [step " << peel_step << "]--" << std::endl;
			//read lower boundary:
			ak::find_first_active_chains(parameters, constrained_model, times, &active_chains, &active_stitches, &rowcol_graph);
			rowcol_graph_tristrip_dirty = true;

			assert(peel_step == 0);
		} else { assert(peel_action == PeelRepeat);
			std::cout << " -- repeat [step " << peel_step << "]--" << std::endl;
			//copy active chains from next_active arrays:
			active_chains = old_next_active_chains;
			active_stitches = old_next_active_stitches;
		}
		show = ShowTimesModel | ShowActiveChains;
		if (active_chains.empty()) return false;
		peel_action = PeelSlice;
		peel_step += 1;

	} else if (peel_action == PeelSlice) {
		std::cout << " -- slice [step " << peel_step << "]--" << std::endl;
		ak::peel_slice(parameters, constrained_model, active_chains, &slice, &slice_on_model, &slice_active_chains, &slice_next_chains, &slice_next_used_boundary);
		slice_times.clear();
		slice_times.reserve(slice_on_model.size());
		for (auto &ev : slice_on_model) {
			slice_times.emplace_back(ev.interpolate(times));
		}

		slice_triangles_dirty = true;
		slice_chains_tristrip_dirty = true;
		show = ShowSlice | ShowSliceChains;

		peel_action = PeelLink;
		peel_step += 1;
	} else if (peel_action == PeelLink) {
		std::cout << " -- link [step " << peel_step << "]--" << std::endl;
		ak::link_chains(parameters, slice, slice_times, slice_active_chains, active_stitches, slice_next_chains, slice_next_used_boundary, &next_stitches, &links);

		links_tristrip_dirty = true;
		show = ShowSlice | ShowSliceChains | ShowLinks;

		peel_action = PeelBuild;
		peel_step += 1;
	} else if (peel_action == PeelBuild) {
		std::cout << " -- build [step " << peel_step << "]--" << std::endl;
		ak::build_next_active_chains(parameters, slice, slice_on_model, slice_active_chains, active_stitches, slice_next_chains, next_stitches, slice_next_used_boundary, links, &next_active_chains, &next_active_stitches, &rowcol_graph);

		rowcol_graph_tristrip_dirty = true;
		next_active_chains_tristrip_dirty = true;
		show = ShowSlice | ShowNextActiveChains;

		peel_action = PeelRepeat;
		peel_step += 1;
	}

	return true;
}

void Interface::clear_traced() {
	traced.clear();

	traced_dirty = true;
	traced_tristrip_dirty = true;
}

void Interface::update_traced() {
	//I guess just update from current rowcol graph, whatever that may be
	traced_dirty = false;

	ak::trace_graph(rowcol_graph, &traced, &constrained_model);

	save_traced();

	traced_tristrip_dirty = true;

	show |= ShowTraced; //<-- slightly hack-y; should really have UI for this sort of stuff
}


void Interface::save_traced() {
	if (save_traced_file == "") return;
	std::cout << "Saving traced stitches to '" << save_traced_file << "'." << std::endl;
	std::vector< Stitch > stitches;
	stitches.reserve(traced.size());
	for (auto const &ts : traced) {
		stitches.emplace_back();
		stitches.back().yarn = ts.yarn;
		stitches.back().type = ts.type;
		stitches.back().direction = ts.dir;
		stitches.back().in[0] = ts.ins[0];
		stitches.back().in[1] = ts.ins[1];
		stitches.back().out[0] = ts.outs[0];
		stitches.back().out[1] = ts.outs[1];
		stitches.back().at = ts.at;
	}

	save_stitches(save_traced_file, stitches);
}


void Interface::update_model_triangles() {
	model_triangles_dirty = false;

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
	std::cout << "Set model_triangles to have " << model_triangles.count << " vertices." << std::endl; //DEBUG
	GL_ERRORS();
}


static void make_sphere(
	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > *attribs,
	glm::vec3 const &c, float r, glm::u8vec4 const &color ) {
	assert(attribs);

	const constexpr uint32_t Slices = 14;
	const constexpr uint32_t Rings = 8;
	static glm::vec2 Ring[Slices];
	static glm::vec2 Slice[Rings+1];

	static bool inited = [](){
		for (uint32_t a = 0; a < Slices; ++a) {
			float ang = a / float(Slices) * 2.0f * float(M_PI);
			Ring[a] = glm::vec2(std::cos(ang), std::sin(ang));
		}

		for (uint32_t a = 0; a <= Rings; ++a) {
			float ang = (a / float(Rings) - 0.5f) * float(M_PI);
			Slice[a] = glm::vec2(std::cos(ang), std::sin(ang));
		}
		return true;
	}();
	(void)inited;

	auto dir = [&](uint32_t ri, uint32_t si) -> glm::vec3 {
		return glm::vec3(
			Slice[ri].x * Ring[si].x,
			Slice[ri].x * Ring[si].y,
			Slice[ri].y
		);
	};

	for (uint32_t ri = 0; ri < Rings; ++ri) {
		attribs->emplace_back(c + r * dir(ri, Slices-1), dir(ri, Slices-1), color);
		attribs->emplace_back(attribs->back());
		attribs->emplace_back(c + r * dir(ri+1, Slices-1), dir(ri+1, Slices-1), color);
		for (uint32_t si = 0; si < Slices; ++si) {
			attribs->emplace_back(c + r * dir(ri, si), dir(ri, si), color);
			attribs->emplace_back(c + r * dir(ri+1, si), dir(ri+1, si), color);
		}
		attribs->emplace_back(attribs->back());
	}
}



static void make_tube(
	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > *attribs,
	glm::vec3 const &a, glm::vec3 const &b, float r, glm::u8vec4 const &color ) {
	assert(attribs);

	static std::vector< glm::vec2 > circle = [](){
		const constexpr uint32_t Angles = 16;
		std::vector< glm::vec2 > ret;
		ret.reserve(Angles);
		for (uint32_t a = 0; a < Angles; ++a) {
			float ang = a / float(Angles) * 2.0f * float(M_PI);
			ret.emplace_back(std::cos(ang), std::sin(ang));
		}
		return ret;
	}();

	glm::vec3 along = glm::normalize(b-a);
	glm::vec3 p1,p2;
	if (std::abs(along.x) <= std::abs(along.y) && std::abs(along.x) <= std::abs(along.z)) {
		p1 = glm::vec3(1.0f, 0.0f, 0.0f);
	} else if (std::abs(along.y) <= std::abs(along.z)) {
		p1 = glm::vec3(0.0f, 1.0f, 0.0f);
	} else {
		p1 = glm::vec3(0.0f, 0.0f, 1.0f);
	}
	p1 = glm::normalize(p1 - glm::dot(along, p1) * along);
	p2 = glm::cross(along, p1);
	glm::mat2x3 xy = glm::mat2x3(p1, p2);

	attribs->emplace_back(a + xy * (r * circle.back()), xy * circle.back(), color);
	attribs->emplace_back(attribs->back());
	attribs->emplace_back(b + xy * (r * circle.back()), xy * circle.back(), color);
	for (auto const &c : circle) {
		attribs->emplace_back(a + xy * (r * c), xy * c, color);
		attribs->emplace_back(b + xy * (r * c), xy * c, color);
	}
	attribs->emplace_back(attribs->back());
}

void Interface::update_constraints_tristrip() {
	constraints_tristrip_dirty = false;

	assert(DEBUG_constraint_paths.size() == constraints.size());

	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > attribs;

	float min_time = -1.0f;
	float max_time = 1.0f;
	for (auto const &c : constraints) {
		min_time = glm::min(min_time, c.value);
		max_time = glm::max(max_time, c.value);
	}

	for (auto const &path : DEBUG_constraint_paths) {
		glm::vec3 color = time_color((constraints[&path - &DEBUG_constraint_paths[0]].value - min_time) / (max_time - min_time));
		glm::u8vec4 color8 = glm::u8vec4(255 * color.r, 255 * color.g, 255 * color.b, 255);
		//generate some sort of path thing (from uncapped tubes, for now):
		for (uint32_t pi = 0; pi + 1 < path.size(); ++pi) {
			glm::vec3 a = path[pi];
			glm::vec3 b = path[pi+1];
			constexpr const float r = 0.02f;
			make_tube(&attribs, a, b, r, color8);
		}
	}

	for (auto const &path : DEBUG_constraint_loops) {
		glm::vec3 color = time_color((constraints[&path - &DEBUG_constraint_loops[0]].value - min_time) / (max_time - min_time));
		glm::u8vec4 color8 = glm::u8vec4(255 * color.r, 255 * color.g, 255 * color.b, 255);
		//generate some sort of path thing (from uncapped tubes, for now):
		for (uint32_t pi = 0; pi + 1 < path.size(); ++pi) {
			glm::vec3 a = path[pi];
			glm::vec3 b = path[pi+1];
			constexpr const float r = 0.01f;
			make_tube(&attribs, a, b, r, color8);
		}
	}

	constraints_tristrip.set(attribs, GL_STATIC_DRAW);
}

void Interface::update_times_model_triangles() {
	times_model_triangles_dirty = false;

	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4, glm::vec2 >::Vertex > attribs;
	attribs.reserve(3 * constrained_model.triangles.size());

	float min =-1.0f;
	float max = 1.0f;

	for (auto v : times) {
		min = std::min(min, v);
		max = std::max(max, v);
	}

	std::vector< glm::vec2 > texcoords;
	if (times.empty()) {
		texcoords.assign(constrained_model.vertices.size(), glm::vec2(0.5f, 0.5f));
	} else {
		assert(times.size() == constrained_model.vertices.size());
		texcoords.reserve(times.size());
		for (auto v : times) {
			texcoords.emplace_back(
				(((v - min) / (max - min)) * (TimeTexSize-1) + 0.5f) / float(TimeTexSize),
				0.5f);
		}
	}

	for (auto const &t : constrained_model.triangles) {
		glm::vec3 const &a = constrained_model.vertices[t.x];
		glm::vec3 const &b = constrained_model.vertices[t.y];
		glm::vec3 const &c = constrained_model.vertices[t.z];

		glm::vec3 n = glm::normalize(glm::cross(b-a, c-a));

		uint32_t idx = &t - &constrained_model.triangles[0];
		glm::u8vec4 id = glm::u8vec4(
			0x3,
			(idx >> 16) & 0xff,
			(idx >> 8) & 0xff,
			idx & 0xff
		);

		glm::vec3 m = (a + b + c) / 3.0f;
		attribs.emplace_back(glm::mix(a, m, 0.1f), n, id, texcoords[t.x]);
		attribs.emplace_back(glm::mix(b, m, 0.1f), n, id, texcoords[t.y]);
		attribs.emplace_back(glm::mix(c, m, 0.1f), n, id, texcoords[t.z]);
	}

	times_model_triangles.set(attribs, GL_STATIC_DRAW);
}


std::vector< std::vector< glm::vec3 > > interpolate_stitch_locations(std::vector< std::vector< glm::vec3 > > const &chains, std::vector< std::vector< ak::Stitch > > const &stitches) {
	assert(stitches.size() == chains.size());

	std::vector< std::vector< glm::vec3 > > locations;
	locations.reserve(locations.size());

	for (auto const &chain : chains) {
		locations.emplace_back();

		uint32_t ci = &chain - &chains[0];
		if (stitches[ci].empty()) continue;

		std::vector< float > lengths;
		lengths.reserve(chain.size());
		lengths.emplace_back(0.0f);
		for (uint32_t pi = 0; pi + 1 < chain.size(); ++pi) {
			glm::vec3 a = chain[pi];
			glm::vec3 b = chain[pi+1];
			lengths.emplace_back(lengths.back() + glm::length(b-a));
		}
		assert(lengths.size() == chain.size());

		locations.back().reserve(stitches[ci].size());
		auto li = lengths.begin();
		for (auto const &s : stitches[ci]) {
			float l = s.t * lengths.back();
			while (li != lengths.end() && *li <= l) ++li;
			assert(li != lengths.end());
			assert(li != lengths.begin());
			float m = (l - *(li-1)) / (*li - *(li-1));
			uint32_t i = li - lengths.begin();
			locations.back().emplace_back( glm::mix(chain[i-1], chain[i], m) );
		}
	}

	return locations;
}


void make_chains_tristrip(
	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > *attribs_,
	std::vector< std::vector< glm::vec3 > > const &chains,
	std::vector< std::vector< ak::Stitch > > const &stitches,
	glm::u8vec4 color8,
	float r) {

	assert(attribs_);
	auto &attribs = *attribs_;

	for (auto const &chain : chains) {
		//the chain itself:
		for (uint32_t pi = 0; pi + 1 < chain.size(); ++pi) {
			glm::vec3 a = chain[pi];
			glm::vec3 b = chain[pi+1];
			//TODO: some sort of direction indication?

			make_sphere(&attribs, a, r, color8);
			make_tube(&attribs, a, b, r, color8);
		}
		if (chain[0] != chain.back()) make_sphere(&attribs, chain.back(), r, color8);
	}
}

void make_stitches_tristrip(
	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > *attribs_,
	std::vector< std::vector< glm::vec3 > > const &chains,
	std::vector< std::vector< ak::Stitch > > const &stitches,
	float r) {
	assert(attribs_);
	auto &attribs = *attribs_;

	std::vector< std::vector< glm::vec3 > > stitch_locations = interpolate_stitch_locations(chains, stitches);

	for (uint32_t ci = 0; ci < chains.size(); ++ci) {
		for (auto const &s : stitches[ci]) {
			glm::vec3 at = stitch_locations[ci][&s - &stitches[ci][0]];

			glm::u8vec4 col = glm::u8vec4(0xff, 0x00, 0xff, 0xff);
			float sr = 1.5f * r;
			if (s.flag == ak::Stitch::FlagDiscard) {
				col = glm::u8vec4(0xff, 0x00, 0x00, 0xff);
			} else if (s.flag == ak::Stitch::FlagLinkOne) {
				col = glm::u8vec4(0x88, 0x88, 0x88, 0xff);
			} else if (s.flag == ak::Stitch::FlagLinkAny) {
				col = glm::u8vec4(0xdd, 0xdd, 0xdd, 0xff);
			}
			make_sphere(&attribs, at, sr, col);
		}
	}
}

std::vector< std::vector< glm::vec3 > > interpolate_locations( ak::Model const &model, std::vector< std::vector< ak::EmbeddedVertex > > const &chains) {
	std::vector< std::vector< glm::vec3 > > locations;
	locations.reserve(chains.size());
	for (auto const &chain : chains) {
		locations.emplace_back();
		locations.back().reserve(chain.size());
		for (auto const &ev : chain) {
			locations.back().emplace_back(ev.interpolate(model.vertices));
		}
	}
	return locations;
}

std::vector< std::vector< glm::vec3 > > copy_locations( ak::Model const &model, std::vector< std::vector< uint32_t > > const &chains) {
	std::vector< std::vector< glm::vec3 > > locations;
	locations.reserve(chains.size());
	for (auto const &chain : chains) {
		locations.emplace_back();
		locations.back().reserve(chain.size());
		for (auto v : chain) {
			locations.back().emplace_back(model.vertices[v]);
		}
	}
	return locations;
}

void Interface::update_rowcol_graph_tristrip() {
	rowcol_graph_tristrip_dirty = false;

	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > attribs;

	std::vector< glm::vec3 > locations;
	locations.reserve(rowcol_graph.vertices.size());
	for (auto const &v : rowcol_graph.vertices) {
		locations.emplace_back(v.at.interpolate(constrained_model.vertices));
	}
	//vertices:
	float row_r = 0.075f * parameters.stitch_width_mm / parameters.model_units_mm;
	float col_r = row_r;
	float stitch_r = 1.5f * row_r;
	#define COL(H) glm::u8vec4(uint32_t(H) >> 24, uint32_t(H) >> 16, uint32_t(H) >> 8, uint32_t(H))
	for (uint32_t vi = 0; vi < rowcol_graph.vertices.size(); ++vi) {
		make_sphere(&attribs,
			locations[vi],
			stitch_r, glm::u8vec4(0x80, 0x80, 0x80, 0xff)
		);
		auto const &v = rowcol_graph.vertices[vi];
		if (v.row_in != -1U) {
			make_tube(&attribs,
				locations[vi],
				0.5f * (locations[v.row_in] + locations[vi]),
				row_r, COL(0xddd6bbff)
			);
		}
		if (v.row_out != -1U) {
			make_tube(&attribs,
				locations[vi],
				0.5f * (locations[v.row_out] + locations[vi]),
				row_r, COL(0xd0cab1ff)
			);
		}
		if (v.col_in[0] != -1U) {
			make_tube(&attribs,
				locations[vi],
				0.5f * (locations[v.col_in[0]] + locations[vi]),
				col_r, (v.col_in[1] == -1U ? COL(0xe1cfe6ff) : COL(0xffdc6aff))
			);
		}
		if (v.col_in[1] != -1U) {
			make_tube(&attribs,
				locations[vi],
				0.5f * (locations[v.col_in[1]] + locations[vi]),
				col_r, COL(0xffdc6aff)
			);
		}

		if (v.col_out[0] != -1U) {
			make_tube(&attribs,
				locations[vi],
				0.5f * (locations[v.col_out[0]] + locations[vi]),
				col_r, (v.col_out[1] == -1U ? COL(0xd2b7daff) : COL(0xfec200ff) )
			);
		}
		if (v.col_out[1] != -1U) {
			make_tube(&attribs,
				locations[vi],
				0.5f * (locations[v.col_out[1]] + locations[vi]),
				col_r, COL(0xfec200ff)
			);
		}

	}
	#undef COL

	rowcol_graph_tristrip.set(attribs, GL_STATIC_DRAW);
}



void Interface::update_active_chains_tristrip() {
	active_chains_tristrip_dirty = false;

	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > attribs;

	std::vector< std::vector< glm::vec3 > > active_locations = interpolate_locations(constrained_model, active_chains);
	make_chains_tristrip(&attribs,
		active_locations, active_stitches,
		glm::u8vec4(0x00, 0xff, 0x88, 0xff),
		0.01f);
	make_stitches_tristrip(&attribs,
		active_locations, active_stitches,
		0.01f);
	active_chains_tristrip.set(attribs, GL_STATIC_DRAW);
}

void Interface::update_slice_triangles() {
	slice_triangles_dirty = false;

	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > attribs;
	attribs.reserve(3 * constrained_model.triangles.size());

	for (auto const &t : slice.triangles) {
		glm::vec3 const &a = slice.vertices[t.x];
		glm::vec3 const &b = slice.vertices[t.y];
		glm::vec3 const &c = slice.vertices[t.z];

		glm::vec3 n = glm::normalize(glm::cross(b-a, c-a));

		glm::u8vec4 color = glm::u8vec4(0xff, 0xff, 0xdd, 0xff);

		glm::vec3 m = (a + b + c) / 3.0f;
		attribs.emplace_back(glm::mix(a, m, 0.1f), n, color);
		attribs.emplace_back(glm::mix(b, m, 0.1f), n, color);
		attribs.emplace_back(glm::mix(c, m, 0.1f), n, color);
	}

	slice_triangles.set(attribs, GL_STATIC_DRAW);
}

void Interface::update_slice_chains_tristrip() {
	slice_chains_tristrip_dirty = false;

	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > attribs;
	std::vector< std::vector< glm::vec3 > > active_locations = copy_locations(slice, slice_active_chains);
	make_chains_tristrip(&attribs,
		active_locations, active_stitches,
		glm::u8vec4(0x22, 0x22, 0xff, 0xff),
		0.01f);
	make_stitches_tristrip(&attribs,
		active_locations, active_stitches,
		0.01f);


	{ //handle chains separately for extra coloring:
		glm::u8vec4 regular = glm::u8vec4(0xff, 0x22, 0x22, 0xff);
		glm::u8vec4 boundary = glm::u8vec4(0xff, 0x55, 0x55, 0xff);
		float r = 0.01f;

		std::vector< std::vector< glm::vec3 > > chains = copy_locations(slice, slice_next_chains);
	
		for (auto const &chain : chains) {
			glm::u8vec4 color8 = (slice_next_used_boundary[&chain - &chains[0]] ? boundary : regular);

			//the chain itself:
			for (uint32_t pi = 0; pi + 1 < chain.size(); ++pi) {
				glm::vec3 a = chain[pi];
				glm::vec3 b = chain[pi+1];
				//TODO: some sort of direction indication?
	
				make_sphere(&attribs, a, r, color8);
				make_tube(&attribs, a, b, r, color8);
			}
			if (chain[0] != chain.back()) make_sphere(&attribs, chain.back(), r, color8);
		}
	}

	slice_chains_tristrip.set(attribs, GL_STATIC_DRAW);
}


void Interface::update_links_tristrip() {
	links_tristrip_dirty = false;

	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > attribs;

	std::vector< std::vector< glm::vec3 > > from = interpolate_stitch_locations(copy_locations(slice, slice_active_chains), active_stitches);
	std::vector< std::vector< glm::vec3 > > to = interpolate_stitch_locations(copy_locations(slice, slice_next_chains), next_stitches);

	for (auto const &link : links) {
		assert(link.from_chain < from.size());
		assert(link.from_stitch < from[link.from_chain].size());
		assert(link.to_chain < to.size());
		assert(link.to_stitch < to[link.to_chain].size());

		glm::vec3 const &a = from[link.from_chain][link.from_stitch];
		glm::vec3 const &b = to[link.to_chain][link.to_stitch];

		make_tube(&attribs, a, b, 0.02f, glm::u8vec4(0x00, 0xff, 0x00, 0xff));
	}

	make_stitches_tristrip(&attribs,
		copy_locations(slice, slice_next_chains), next_stitches,
		0.01f);

	links_tristrip.set(attribs, GL_STATIC_DRAW);
}

void Interface::update_next_active_chains_tristrip() {
	next_active_chains_tristrip_dirty = false;

	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > attribs;

	std::vector< std::vector< glm::vec3 > > next_active_locations = interpolate_locations(constrained_model, next_active_chains);
	make_chains_tristrip(&attribs,
		next_active_locations, next_active_stitches,
		glm::u8vec4(0x00, 0xff, 0x88, 0xff),
		0.01f);
	make_stitches_tristrip(&attribs,
		next_active_locations, next_active_stitches,
		0.01f);
	next_active_chains_tristrip.set(attribs, GL_STATIC_DRAW);
}

void Interface::update_traced_tristrip() {
	traced_tristrip_dirty = false;

	static std::vector< glm::u8vec4 > yarn_colors{
		glm::u8vec4(0xee, 0xbb, 0x55, 0xff),
		glm::u8vec4(0xbb, 0x55, 0xee, 0xff),
		glm::u8vec4(0x55, 0xee, 0xbb, 0xff),
		glm::u8vec4(0xbb, 0xee, 0x55, 0xff),
		glm::u8vec4(0x55, 0xbb, 0xee, 0xff),
		glm::u8vec4(0xee, 0x55, 0xbb, 0xff)
	};

	uint32_t yarn_color = 0;
	std::vector< GLAttribBuffer< glm::vec3, glm::vec3, glm::u8vec4 >::Vertex > attribs;
	for (uint32_t ti = 1; ti < traced.size(); ++ti) {
		if (traced[ti].yarn != traced[ti-1].yarn) {
			yarn_color = (yarn_color + 1) % yarn_colors.size();
			continue;
		}
		make_tube(&attribs,
			traced[ti-1].at, traced[ti].at,
			0.01f,
			yarn_colors[yarn_color]);
	}
	for (uint32_t ti = 0; ti < traced.size(); ++ti) {
		glm::vec3 const &at = traced[ti].at;
		if (traced[ti].ins[0] != -1U && traced[ti].ins[1] != -1U) {
			glm::vec3 const &in0 = traced[traced[ti].ins[0]].at;
			glm::vec3 const &in1 = traced[traced[ti].ins[1]].at;
			make_tube(&attribs, 0.5f * (at + in0), at, 0.005f, glm::u8vec4(0x55, 0x55, 0x55, 0xff));
			make_tube(&attribs, 0.5f * (at + in1), at, 0.005f, glm::u8vec4(0xcc, 0xcc, 0xcc, 0xff));
		}
		if (traced[ti].outs[0] != -1U && traced[ti].outs[1] != -1U) {
			glm::vec3 const &out0 = traced[traced[ti].outs[0]].at;
			glm::vec3 const &out1 = traced[traced[ti].outs[1]].at;
			make_tube(&attribs,at, 0.5f * (at + out0), 0.005f, glm::u8vec4(0x44, 0x44, 0x44, 0xff));
			make_tube(&attribs,at, 0.5f * (at + out1), 0.005f, glm::u8vec4(0xbb, 0xbb, 0xbb, 0xff));
		}
		if (traced[ti].ins[0] != -1U && traced[ti].ins[1] == -1U) {
			glm::vec3 const &in0 = traced[traced[ti].ins[0]].at;
			make_tube(&attribs, 0.5f * (at + in0), at, 0.005f, glm::u8vec4(0x88, 0x88, 0x88, 0xff));
		}
		if (traced[ti].outs[0] != -1U && traced[ti].outs[1] == -1U) {
			glm::vec3 const &out0 = traced[traced[ti].outs[0]].at;
			make_tube(&attribs, at, 0.5f * (at + out0), 0.005f, glm::u8vec4(0x77, 0x77, 0x77, 0xff));
		}
		assert(!(traced[ti].ins[0] == -1U && traced[ti].ins[1] != -1U));
		assert(!(traced[ti].outs[0] == -1U && traced[ti].outs[1] != -1U));
	}

	traced_tristrip.set(attribs, GL_STATIC_DRAW);
}




void Interface::DEBUG_test_linking(bool flip) {
	//save_constraints_file = ""; <-- should probably actually save constraints
	clear_constraints();

	glm::vec3 min = glm::vec3( std::numeric_limits< float >::infinity());
	glm::vec3 max = glm::vec3(-std::numeric_limits< float >::infinity());

	for (auto const &v : model.vertices) {
		min = glm::min(min, v);
		max = glm::max(max, v);
	}

	//build next-along-triangle map for directed edges:
	std::unordered_map< glm::uvec2, uint32_t > next;
	for (auto const &tri : model.triangles) {
		auto do_edge = [&](uint32_t a, uint32_t b, uint32_t c) {
			auto ret = next.insert(std::make_pair(glm::uvec2(a,b), c));
			assert(ret.second);
		};
		do_edge(tri.x, tri.y, tri.z);
		do_edge(tri.y, tri.z, tri.x);
		do_edge(tri.z, tri.x, tri.y);
	}

	std::vector< ak::Constraint > new_constraints;
	std::unordered_set< glm::uvec2 > visited;
	for (auto const &en : next) {
		//not a boundary:
		if (next.count(glm::uvec2(en.first.y, en.first.x))) continue;
		//already done:
		if (visited.count(en.first)) continue;

		std::vector< uint32_t > path;
		path.emplace_back(en.first.y);
		path.emplace_back(en.first.x);
		visited.insert(en.first);
		while (true) {
			glm::uvec2 e(path[path.size()-2], path[path.size()-1]);
			while (true) {
				auto f = next.find(glm::uvec2(e.y, e.x));
				if (f == next.end()) break;
				e = glm::uvec2(f->second, e.y);
				assert(e != glm::uvec2(path[path.size()-2], path[path.size()-1]));
			}
			assert(e != glm::uvec2(path[path.size()-2], path[path.size()-1]));
			assert(!next.count(glm::uvec2(e.y, e.x)));

			path.emplace_back(e.x);
			if (e.y == path[0] && e.x == path[1]) break;

			auto ret = visited.insert(e);
			assert(ret.second);
		}
		assert(path[0] == path[path.size()-2] && path[1] == path[path.size()-1]);
		path.pop_back();

		uint32_t start_votes = 0;
		uint32_t end_votes = 0;
		for (auto v : path) {
			if (model.vertices[v].z > 0.75f * (max.z + min.z)) {
				++end_votes;
			}
			if (model.vertices[v].z < 0.25f * (max.z + min.z)) {
				++start_votes;
			}
		}

		if (start_votes == end_votes) {
			std::cout << "Skipping constraint -- same number of start and end votes." << std::endl;
		}

		new_constraints.emplace_back();
		new_constraints.back().chain = path;
		if (flip) {
			new_constraints.back().value = (start_votes > end_votes ?  1.0f :-1.0f);
		} else {
			new_constraints.back().value = (start_votes > end_votes ? -1.0f : 1.0f);
		}
		new_constraints.back().radius = 0.0f;
	}

	std::cout << "Made " << new_constraints.size() << " constraints." << std::endl;

	set_constraints(new_constraints);

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
