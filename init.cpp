
#include "Interface.hpp"
#include "TaggedArguments.hpp"

#include "Stitch.hpp"

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
	std::string constraints_file = "";
	std::string save_traced_file = "";
	std::string load_traced_file = "";
	float constraint_tolerance = 0.01f;
	int32_t peel_test = 0;
	int32_t peel_step = 0;
	int32_t test_constraints = 0;
	ak::Parameters parameters;
	{
		TaggedArguments args;
		args.emplace_back("obj", &obj_file, "input obj file (required)");
		args.emplace_back("obj-scale", &parameters.model_units_mm, "length of one unit in obj file (mm)");
		args.emplace_back("test-constraints", &test_constraints, "if non-zero, generate linking-test-style constraints [+z boundaries to 1.0, -z boundaries to -1.0; flipped if negative");
		args.emplace_back("load-constraints", &load_constraints_file, "file to load time constraints from");
		args.emplace_back("save-constraints", &save_constraints_file, "file to save time constraints to");
		args.emplace_back("constraints", &constraints_file, "try to load constraints from the named file, and definitely save them to it (load_constraints_file or save_constraints_file will override)");
		args.emplace_back("constraint-tolerance", &constraint_tolerance, "maximum distance away from a model vertex that a loaded constraint vertex can be before it is considered missing");
		args.emplace_back("save-traced", &save_traced_file, "save traced stitches to this file");
		args.emplace_back("load-traced", &load_traced_file, "load and view traced stitches from this file (may break other interface actions)");
		args.emplace_back("stitch-width", &parameters.stitch_width_mm, "stitch width (mm)");
		args.emplace_back("stitch-height", &parameters.stitch_height_mm, "stitch height (mm)");
		args.emplace_back("peel-test", &peel_test, "run N rounds of peeling then quit (-1 to run until done)");
		args.emplace_back("peel-step", &peel_step, "run N rounds of peeling then show interface (-1 to run until done)");
		bool usage = !args.parse(kit::args);
		if (!usage && obj_file == "") {
			std::cerr << "ERROR: 'obj:' argument is required." << std::endl;
			usage = true;
		}
		if (!usage && (peel_step != 0 && peel_test != 0)) {
			std::cerr << "ERROR: Please specify only one of 'peel-test:' and 'peel-step:'" << std::endl;
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
		ak::load_constraints(model, load_constraints_file, constraint_tolerance, &constraints);
	} else if (constraints_file != "") {
		try {
			ak::load_constraints(model, constraints_file, constraint_tolerance, &constraints);
		} catch (std::exception const &e) {
			std::cerr << "WARNING: failed to load from '" << constraints_file << "' (" << e.what() << ")>" << std::endl;
		}
	}

	std::shared_ptr< Interface > interface = std::make_shared< Interface >();

	interface->parameters = parameters;

	if (save_constraints_file != "") {
		interface->save_constraints_file = save_constraints_file;
	} else if (constraints_file != "") {
		interface->save_constraints_file = constraints_file;
	}

	interface->set_model(model);
	interface->set_constraints(constraints);

	if (test_constraints != 0) {
		interface->DEBUG_test_linking(test_constraints < 0);
	}

	//hack to make it easier to view output ".st" files:
	if (load_traced_file != "") {
		std::vector< Stitch > stitches;
		if (!load_stitches(load_traced_file, &stitches)) {
			std::cerr << "ERROR: failed to load traced stitches from '" << load_traced_file << "'." << std::endl;
			return nullptr;
		}
		interface->traced.clear();
		interface->traced.reserve(stitches.size());
		for (auto const &stitch : stitches) {
			ak::TracedStitch ts;
			ts.yarn = stitch.yarn;
			ts.type = ak::TracedStitch::Type(stitch.type);
			ts.dir = ak::TracedStitch::Dir(stitch.direction);
			ts.ins[0] = stitch.in[0];
			ts.ins[1] = stitch.in[1];
			ts.outs[0] = stitch.out[0];
			ts.outs[1] = stitch.out[1];
			ts.at = stitch.at;
			interface->traced.emplace_back(ts);
		}

		std::cout << "Loaded traced stitches from '" << load_traced_file << "' -- these will probably be visible until you do anything with tracing." << std::endl;
		interface->show = Interface::ShowTraced;
		interface->traced_tristrip_dirty = true;
	}

	if (peel_test != 0 || peel_step != 0) {
		uint32_t target = (peel_test != 0 ? uint32_t(peel_test) : uint32_t(peel_step));
		interface->clear_peeling();
		while (interface->peel_step <= target) {
			if (!interface->step_peeling()) {
				std::cout << "--- NOTE: peeling finished ---" << std::endl;
				break;
			}
		}
		if (save_traced_file != "") {
			interface->save_traced_file = save_traced_file;
			interface->update_traced();
		}
		if (peel_test != 0) return nullptr;
	}

	interface->save_traced_file = save_traced_file;

	return interface;
}
