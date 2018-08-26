#include "Shape.hpp"

#include <vector>
#include <limits>
#include <iostream>
#include <string>
#include <algorithm>

int main() {
	for (uint32_t count : {0,1,2,3,4,5,6,7,10,21,37,100,101}) {
		std::vector< Shape > shapes = Shape::make_shapes_for(count);

		assert(shapes.size() == Shape::count_shapes_for(count));

		std::vector< uint32_t > inds;
		inds.reserve(count);
		for (uint32_t i = 0; i < count; ++i) {
			inds.emplace_back(i);
		}

		for (auto const &shape : shapes) {
			//fits on the bed as expected:
			std::vector< uint32_t > front, back;
			shape.append_to_beds(inds, -1U, &front, &back);

			{ //dump for reference:
				auto str = [](uint32_t i) {
					std::string rep;
					if (i == -1U) rep = '.';
					else rep = std::to_string(i);
					while (rep.size() < 3) rep = ' ' + rep;
					return rep;
				};
				std::cout << count << "[" << &shape - &shapes[0] << "/" << shapes.size() << "]:\n";
				for (auto i : front) std::cout << str(i);
				std::cout << '\n';
				for (auto i : back) std::cout << str(i);
				std::cout << '\n';
				std::cout.flush();
			}

			//pack/unpack:
			assert(Shape::unpack(shape.pack()).pack() == shape.pack());

			//indexing:
			assert(shape.index_for(count) == &shape - &shapes[0]);


			{ //does everything appear?
				std::vector< bool > used(inds.size(), false);
				for (auto i : front) {
					if (i != -1U) {
						assert(i < count);
						assert(used[i] == false);
						used[i] = true;
					}
				}
				for (auto i : back) {
					if (i != -1U) {
						assert(i < count);
						assert(used[i] == false);
						used[i] = true;
					}
				}
				for (auto u : used) {
					assert(u);
				}
			}

			{ //is it "as-left-as-possible"? + "no-right-padding"?
				if (front.empty() && back.empty()) {
				} else if (front.empty() && !back.empty()) {
					assert(back[0] != -1U);
					assert(back.back() != -1U);
				} else if (!front.empty() && back.empty()) {
					assert(front[0] != -1U);
					assert(front.back() != -1U);
				} else { assert(!front.empty() && !back.empty());
					assert(front[0] != -1U || back[0] != -1U);
					assert(front.back() != -1U);
					assert(back.back() != -1U);
				}
			}

			{ //does it match the helper function results? (size_index_to_bed_needle)
				for (uint32_t i = 0; i < count; ++i) {
					char bed;
					int32_t needle;
					shape.size_index_to_bed_needle(count, i, &bed, &needle);
					if (bed == 'f') {
						assert(uint32_t(needle) < front.size());
						assert(front[needle] == i);
					} else { assert(bed == 'b');
						assert(uint32_t(needle) < back.size());
						assert(back[needle] == i);
					}
				}
			}

			{ //does it match the helper function results? (size_to_range)
				int32_t front_min, front_max, back_min, back_max;
				shape.size_to_range(count, &front_min, &front_max, &back_min, &back_max);

				int32_t real_front_min = std::numeric_limits< int32_t >::max();
				int32_t real_front_max = std::numeric_limits< int32_t >::min();
				int32_t real_back_min = std::numeric_limits< int32_t >::max();
				int32_t real_back_max = std::numeric_limits< int32_t >::min();
				for (uint32_t i = 0; i < front.size(); ++i) {
					if (front[i] != -1U) {
						real_front_min = std::min(real_front_min, int32_t(i));
						real_front_max = std::max(real_front_max, int32_t(i));
					}
				}
				for (uint32_t i = 0; i < back.size(); ++i) {
					if (back[i] != -1U) {
						real_back_min = std::min(real_back_min, int32_t(i));
						real_back_max = std::max(real_back_max, int32_t(i));
					}
				}

				assert(real_front_min == front_min);
				assert(real_front_max == front_max);
				assert(real_back_min == back_min);
				assert(real_back_max == back_max);
			}

		}
	}

	return 0;
}

