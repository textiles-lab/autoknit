
#include "Interface.hpp"

#include <kit/kit.hpp>

kit::Config kit_config() {
	kit::Config config;
	config.size = glm::uvec2(1000, 800);
	config.title = "autoknit";
	return config;
}

std::shared_ptr< kit::Mode > kit_mode() {
	std::string model = "";
	for (auto const &arg : kit::args) {
		std::string tag = arg;
		std::string value = "";
		for (uint32_t i = 0; i < tag.size(); ++i) {
			if (tag[i] == ':') {
				value = tag.substr(i+1);
				tag = tag.substr(0,i);
				break;
			}
		}
		if (tag == "model") {
			model = value;
		} else {
			std::cerr << "Unrecognized tag:value argument pair '" << tag << ':' << value << "'." << std::endl;
			return nullptr;
		}
	}

	if (model == "") {
		std::cerr << "Please pass a model:file.obj parameter." << std::endl;
	}

	std::shared_ptr< Interface > interface = std::make_shared< Interface >();

	return std::make_shared< Interface >(model);

}
