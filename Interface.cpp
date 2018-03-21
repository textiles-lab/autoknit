#include "Interface.hpp"

#include <kit/GLProgram.hpp>
#include <kit/Load.hpp>
#include <kit/gl_errors.hpp>

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
		"in vec3 normal;\n"
		"in vec3 position;\n"
		"layout(location = 0) out vec4 fragColor;\n"
		"void main() {\n"
		"	vec3 n = normalize(normal);\n"
		"	vec3 l = vec3(0.0, 0.0, 1.0);\n"
		"	float nl = dot(n, l) * 0.5 + 0.5;\n"
		"	fragColor = vec4(nl * (0.5 * normal + 0.5), 1.0);\n"
		"}\n"
	);

	model_draw_p2c = ret->getUniformLocation("p2c", GLProgram::MissingIsError);
	model_draw_p2l = ret->getUniformLocation("p2l", GLProgram::MissingIsWarning);
	model_draw_n2l = ret->getUniformLocation("n2l", GLProgram::MissingIsWarning);

	return ret;
});

Interface::Interface() {
	model_triangles_for_model_draw = GLVertexArray::make_binding(model_draw->program, {
		{model_draw->getAttribLocation("Position", GLProgram::MissingIsError), model_triangles[0]},
		{model_draw->getAttribLocation("Normal", GLProgram::MissingIsWarning), model_triangles[1]}
	});
}

Interface::~Interface() {
}

void Interface::update(float elapsed) {
}

void Interface::draw() {
	glViewport(0, 0, kit::display.size.x, kit::display.size.y);

	glClearColor(0.9, 0.9, 0.9, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);

	//Position-to-clip matrix:
	glm::mat4 p2c =
		glm::infinitePerspective(camera.fovy, kit::display.aspect, 0.1f)
		* camera.mv()
		;
	//Position-to-light matrix:
	glm::mat4x3 p2l = glm::mat4(1.0f);
	//Normal-to-light matrix:
	glm::mat3 n2l = glm::inverse(glm::transpose(glm::mat3(p2l)));

	glUseProgram(model_draw->program);

	glUniformMatrix4fv(model_draw_p2c, 1, GL_FALSE, glm::value_ptr(p2c));
	glUniformMatrix4x3fv(model_draw_p2l, 1, GL_FALSE, glm::value_ptr(p2l));
	glUniformMatrix3fv(model_draw_n2l, 1, GL_FALSE, glm::value_ptr(n2l));

	glBindVertexArray(model_triangles_for_model_draw.array);

	glDrawArrays(GL_TRIANGLES, 0, model_triangles.count);

	glBindVertexArray(0);

	glUseProgram(0);

	GL_ERRORS();
}

void Interface::pointer_action(kit::PointerID pointer, kit::PointerAction action, kit::Pointer const &old_state, kit::Pointer const &new_state) {
	
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

void Interface::set_model(Model const &new_model) {
	model = new_model;

	std::vector< GLAttribBuffer< glm::vec3, glm::vec3 >::Vertex > attribs;
	attribs.reserve(3 * model.triangles.size());

	for (auto const &t : model.triangles) {
		glm::vec3 const &a = model.vertices[t.x];
		glm::vec3 const &b = model.vertices[t.y];
		glm::vec3 const &c = model.vertices[t.z];

		glm::vec3 n = glm::normalize(glm::cross(b-a, c-a));

		attribs.emplace_back(a, n);
		attribs.emplace_back(b, n);
		attribs.emplace_back(c, n);
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
