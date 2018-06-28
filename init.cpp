
#include "Interface.hpp"
#include "TaggedArguments.hpp"

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
	ak::Parameters parameters;
	{
		TaggedArguments args;
		args.emplace_back("obj", &obj_file, "input obj file (required)");
		args.emplace_back("obj-scale", &parameters.model_units_mm, "length of one unit in obj file (mm)");
		args.emplace_back("load-constraints", &load_constraints_file, "file to load time constraints from");
		args.emplace_back("save-constraints", &save_constraints_file, "file to save time constraints to");
		args.emplace_back("stitch-width", &parameters.stitch_width_mm, "stitch width (mm)");
		args.emplace_back("stitch-height", &parameters.stitch_height_mm, "stitch height (mm)");
		bool usage = !args.parse(kit::args);
		if (!usage && obj_file == "") {
			std::cerr << "ERROR: 'obj:' argument is required." << std::endl;
			usage = true;
		}
		if (usage) {
			std::cerr << "Usage:\n\t./interface [tag:value] [...]\n" << args.help_string() << std::endl;
			return nullptr;
		}
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

	interface->parameters = parameters;

	interface->save_constraints_file = save_constraints_file;

	interface->set_model(model);
	interface->set_constraints(constraints);

	return interface;
}
