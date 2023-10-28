#pragma once

#include <vector>
#include <string>
#include <functional>

template< typename T >
void typeset_beds(
	std::vector< T > const &front,
	std::vector< T > const &back,
	std::function< std::string(T) > const &to_string,
	std::string const &gap,
	std::string *front_string_,
	std::string *back_string_ ) {

	std::string temp_a, temp_b;
	auto &front_string = (front_string_ ? *front_string_ : temp_a);
	auto &back_string = (back_string_ ? *back_string_ : temp_b);

	std::vector< std::vector< std::string > > columns;
	columns.resize(std::max(front.size(), back.size()));
	for (uint32_t i = 0; i < columns.size(); ++i) {
		if (i < back.size()) {
			columns[i].emplace_back(to_string(back[i]));
		} else {
			columns[i].emplace_back("");
		}
		if (i < front.size()) {
			columns[i].emplace_back(to_string(front[i]));
		} else {
			columns[i].emplace_back("");
		}
	}
	for (auto &column : columns) {
		uint32_t len = 0;
		for (auto const &elt : column) {
			len = std::max(len, uint32_t(elt.size()));
		}
		for (auto &elt : column) {
			while (elt.size() < len) {
				elt = " " + elt;
			}
		}
	}
	for (uint32_t i = 0; i < columns.size(); ++i) {
		if (i > 0) {
			back_string += gap;
			front_string += gap;
		}
		assert(columns[i].size() == 2);
		back_string += columns[i][0];
		front_string += columns[i][1];
	}
}
