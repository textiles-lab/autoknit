#include "Interface.hpp"


Interface::Interface() {
}

Interface::~Interface() {
}

void Interface::update(float elapsed) {
}

void Interface::draw() {
	glm::mat4 mvp;
	{ //compute mvp from camera position:
		mvp = glm::infinitePerpective(camera.fovy, kit::display.aspect, 0.1f)
		    * glm::mat4_cast(glm::inverse(camera.rot));
		    * glm::mat4(
				glm::vec4(1.0f, 0.0f, 0.0f, 0.0f),
				glm::vec4(0.0f, 1.0f, 0.0f, 0.0f),
				glm::vec4(0.0f, 0.0f, 1.0f, 0.0f),
				glm::vec4(-camera.pos, 1.0f)
			);
	}
}

void Interface::pointer_action(PointerID pointer, PointerAction action, Pointer const &old_state, Pointer const &new_state) {
}


