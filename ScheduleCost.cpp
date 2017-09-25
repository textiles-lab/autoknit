#include "ScheduleCost.hpp"
#include "plan_transfers.hpp"

ScheduleCost ScheduleCost::shape_cost( Shape const &shape ) {
	ScheduleCost ret;
	ret.shape =
		  ((shape.nibbles & Shape::BackLeft) ? 1 : 0)
		| ((shape.nibbles & Shape::BackRight) ? 1 : 0)
		| ((shape.nibbles & Shape::FrontLeft) ? 1 : 0)
		| ((shape.nibbles & Shape::FrontRight) ? 1 : 0);
	return ret;
}

ScheduleCost ScheduleCost::transfer_cost(
		uint32_t from_count, Shape from_shape,
		uint32_t to_count, Shape to_shape,
		std::vector< uint32_t > const &map
	) {
	assert(map.size() == from_count);

	auto make_needles = [](uint32_t count, Shape shape) -> std::vector< BedNeedle > {
		std::vector< BedNeedle > ret(count);
		std::vector< uint32_t > inds;
		inds.reserve(count);
		for (uint32_t i = 0; i < count; ++i) {
			inds.emplace_back(i);
		}
		std::vector< uint32_t > front_inds, back_inds;
		shape.append_to_beds(inds, -1U, &front_inds, &back_inds);
		for (uint32_t i = 0; i < front_inds.size(); ++i) {
			if (front_inds[i] == -1U) continue;
			assert(front_inds[i] < ret.size());
			assert(ret[front_inds[i]].needle == 0);
			ret[front_inds[i]].bed = BedNeedle::Front;
			ret[front_inds[i]].needle = i + 1;
		}
		for (uint32_t i = 0; i < back_inds.size(); ++i) {
			if (back_inds[i] == -1U) continue;
			assert(back_inds[i] < ret.size());
			assert(ret[back_inds[i]].needle == 0);
			ret[back_inds[i]].bed = BedNeedle::Back;
			ret[back_inds[i]].needle = i + 1;
		}
		for (auto const &bn : ret) {
			assert(bn.needle != 0);
		}
		return ret;
	};
	std::vector< BedNeedle > from  = make_needles(from_count, from_shape);
	std::vector< BedNeedle > to_needles = make_needles(to_count, to_shape);
	std::vector< BedNeedle > to;
	to.reserve(from.size());
	for (uint32_t i = 0; i < from.size(); ++i) {
		assert(map[i] < to_count);
		to.emplace_back(to_needles[map[i]]);
	}
	assert(to.size() == from.size());

	//TODO: compute slack
	//TODO: actually call xfer planning?

	ScheduleCost cost;

	for (uint32_t i = 0; i < from.size(); ++i) {
		if (from[i].bed != to[i].bed) {
			cost.roll += 1;
		} else if (from[i].needle != to[i].needle) {
			cost.shift += 1;
		}
	}

	return cost;
}

