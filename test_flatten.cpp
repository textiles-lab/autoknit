#include <iostream>
#include <cassert>
#include <vector>

//'flatten' is defined in ak-link_chains.cpp
void flatten(std::vector< uint32_t > &closest, std::vector< float > const &lengths, bool is_loop);

int main(int argc, char **argv) {

	auto test = [](std::vector< uint32_t > closest, std::vector< float > lengths) {
		if (lengths.empty()) {
			lengths.assign(closest.size(), 1.0f);
		}
		assert(lengths.size() == closest.size());

		std::vector< uint32_t > before = closest;

		std::cout << "Flatten from:";
		for (auto c : closest) {
			std::cout << ' ' << c;
		};
		std::cout << std::endl;

		flatten(closest, lengths, true);

		std::cout << "          to:";
		for (auto c : closest) {
			std::cout << ' ' << c;
		};
		std::cout << std::endl;
		
	};

	test(std::vector< uint32_t >{5,5,5,5,5,5}, std::vector< float >{});
	test(std::vector< uint32_t >{5,5,5,2,2,5}, std::vector< float >{});
	test(std::vector< uint32_t >{1,1,2,1,2,2,2,1},
	        std::vector< float >{1,1,1,1,1,1,1,1});

	test(std::vector< uint32_t >{1,1,2,1,2,2,2,1},
	        std::vector< float >{1,1,1,5,1,1,1,1});
	
	test(std::vector< uint32_t >{5,1,1,1,2,2,5,5,3,4,4,4,4,1,1,3,3,5,5,2,2,5}, std::vector< float >{});

	return 0;
}
