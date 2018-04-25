
#include "Interface.hpp"

#include <kit/kit.hpp>
#include <kit/Load.hpp>

kit::Config kit_config() {
	kit::Config config;
	config.size = glm::uvec2(1000, 800);
	config.title = "autoknit";
	return config;
}

std::shared_ptr< kit::Mode > kit_mode() {
	kit::call_load_functions();

	std::string obj_file = "";
	std::string load_constraints_file = "";
	std::string save_constraints_file = "";
	for (auto const &arg : kit::args) {
		if (&arg == &kit::args[0]) continue;

		std::string tag = arg;
		std::string value = "";
		for (uint32_t i = 0; i < tag.size(); ++i) {
			if (tag[i] == ':') {
				value = tag.substr(i+1);
				tag = tag.substr(0,i);
				break;
			}
		}
		if (tag == "obj") {
			obj_file = value;
		} else if (tag == "load-constraints") {
			load_constraints_file = value;
		} else if (tag == "save-constraints") {
			save_constraints_file = value;
		} else {
			std::cerr << "Unrecognized tag:value argument pair '" << tag << ':' << value << "'." << std::endl;
			return nullptr;
		}
	}

	if (obj_file == "") {
		std::cerr << "Please pass an obj:file.obj parameter." << std::endl;
		return nullptr;
	}
	ak::Model model;
	ak::load_obj(obj_file, &model);

	if (model.triangles.empty()) {
		std::cerr << "ERROR: model is empty." << std::endl;
		return nullptr;
	}

	std::vector< ak::Constraint > constraints;
	if (load_constraints_file != "") {
		ak::load_constraints(model, load_constraints_file, &constraints);
	}

	std::shared_ptr< Interface > interface = std::make_shared< Interface >();

	interface->save_constraints_file = save_constraints_file;

	interface->set_model(model);
	interface->set_constraints(constraints);

	return interface;
}
