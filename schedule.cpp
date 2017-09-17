#include "Stitch.hpp"

#include "TaggedArguments.hpp"

int main(int argc, char **argv) {
	std::string in_st = "";
	std::string out_js = "";
	{ //parse arguments:
		TaggedArguments args;
		args.emplace_back("st", &in_st, "input stitches file (required)");
		args.emplace_back("js", &out_js, "output knitting file");
		bool usage = !args.parse(argc, argv);
		if (!usage && in_st == "") {
			std::cerr << "ERROR: 'st:' argument is required." << std::endl;
			usage = true;
		}
		if (usage) {
			std::cerr << "Usage:\n\t./schedule [tag:value] [...]\n" << args.help_string() << std::endl;
			return 1;
		}
	}

	std::vector< Stitch > stitches;
	if (!load_stitches(in_st, &stitches)) {
		std::cerr << "ERROR: failed to load stitches from '" << in_st << "'." << std::endl;
		return 1;
	}
	std::cout << "Read " << stitches.size() << " stitches from '" << in_st << "'." << std::endl;

	//New scheduling workflow:
	// (1) split into "steps" ==> tube-supported bits of knitting that will be done at once
	//    - each step eventually needs a shape + roll + offset for its output loops
	//      (implies a shape for input loops)
	// (2) pick a shape + roll for "interesting" steps
	//    - these are steps that take loops from more than one output or input
	// (2) figure out a layout (shape + roll) for each step

	return 0;
}
