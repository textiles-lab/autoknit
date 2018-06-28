#pragma once

#include <vector>
#include <string>
#include <functional>
#include <sstream>
#include <iostream>

struct TaggedArgument {
	std::string tag;
	std::function< bool(std::string) > parse;
	std::string help;

	TaggedArgument(std::string const &tag_, std::function< bool(std::string) > const &parse_, std::string const &help_) : tag(tag_), parse(parse_), help(help_) {
	}
	TaggedArgument(std::string const &tag_, float *value, std::string const &help_)
		: TaggedArgument(
			tag_,
			[value](std::string arg) -> bool {
				std::istringstream iss(arg);
				if (iss >> *value) return true;
				else return false;
			},
			tag_ + ":<float> (default:" + std::to_string(*value) + ") " + help_ ) {
	}
	TaggedArgument(std::string const &tag_, std::string *value, std::string const &help_)
		: TaggedArgument(
			tag_,
			[value](std::string arg) -> bool {
				*value = arg;
				return true;
			},
			tag_ + ":<string> " + (*value != "" ? "default: \"" + *value + "\") " : "") + help_
		) {
	}
};

struct TaggedArguments : std::vector< TaggedArgument > {
	bool parse(int argc, char **argv) {
		std::vector< std::string > args;
		for (int i = 0; i < argc; ++i) {
			args.emplace_back(argv[i]);
		}
		return parse(args);
	}
	bool parse(const std::vector< std::string > &args) {
		bool success = true;
		for (const std::string &arg : args) {
			if (&arg == &args[0]) continue;

			std::string tag = "";
			std::string value = arg;
			auto idx = value.find(':');
			if (idx != std::string::npos) {
				tag = value.substr(0, idx);
				value = value.substr(idx + 1);
			}
			TaggedArgument const *found = nullptr;
			for (auto const &ta : *this) {
				if (ta.tag == tag) {
					found = &ta;
					break;
				}
			}
			if (!found) {
				std::cerr << "ERROR: No tag match for argument '" << arg << "'." << std::endl;
				success = false;
			} else if (!found->parse(value)) {
				std::cerr << "ERROR: Failed to parse argument '" << arg << "'." << std::endl;
				success = false;
			}
		}
		return success;
	}
	std::string help_string() const {
		std::string ret = "Options:";
		for (auto const &ta : *this) {
			ret += "\n  " + ta.help;
		}
		return ret;
	}
};
