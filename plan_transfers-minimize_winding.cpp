#include "plan_transfers-helpers.hpp"

#include <algorithm>

void minimize_winding(std::vector< int32_t > *winding_) {
	assert(winding_);
	auto &winding = *winding_;

	if (winding.empty()) return;
	//can be minimized by subtracting one of the multiples of two that bracket the median:
	std::vector< int32_t > temp = winding;

	//DEBUG: dump numbers
	//for (auto t : temp) std::cout << ' ' << t;
	//std::cout << '\n';

	std::sort(temp.begin(), temp.end());
	int32_t twice_median;
	if (temp.size() == 1) twice_median = 2 * temp[0];
	else twice_median = temp[temp.size()/2] + temp[(temp.size()+1)/2];

	int32_t close_multiple_a = (twice_median / 4) * 2;
	int32_t close_multiple_b = ((twice_median + (twice_median < 0 ? -4 : 4)) / 4) * 2;


	int32_t sum_a = 0;
	int32_t sum_b = 0;
	for (auto &w : winding) {
		sum_a += std::abs(w - close_multiple_a);
		sum_b += std::abs(w - close_multiple_b);
	}
	int32_t close_multiple = (sum_a <= sum_b ? close_multiple_a : close_multiple_b);

	//std::cout << "close multiples: " << close_multiple_a << " (sum:" << sum_a << "), " << close_multiple_b << " (sum: " << sum_b << ")" << std::endl;

	int32_t DEBUG_before = 0;
	int32_t DEBUG_after = 0;
	int32_t DEBUG_after_minus = 0;
	int32_t DEBUG_after_plus = 0;

	for (auto &w : winding) {
		DEBUG_before += std::abs(w);
		w -= close_multiple;
		DEBUG_after += std::abs(w);
		DEBUG_after_minus += std::abs(w - 2);
		DEBUG_after_plus += std::abs(w + 2);
	}
	assert(DEBUG_after <= DEBUG_before);
	assert(DEBUG_after <= DEBUG_after_minus);
	assert(DEBUG_after <= DEBUG_after_plus);
}
