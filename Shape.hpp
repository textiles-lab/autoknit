#pragma once

//Shape represents the shape of a cycle of stitches on the bed.
// It is a map from stitch indices to needles, such that each stitch
//is on a needle adjacent to its index-wise neighbors.

//For a given count, shapes are controlled by a type and offset:
//Valid Shapes:
//   . . o o o o . . .
//   . . o o o o . . .

//   . . . o o . . . .
//   . . o o o o . . .

//   . . . o o o . . .
//   . . o o o . . . .

//   . . o o o . . . .
//   . . . o o o . . .

//   . . o o o o . . .
//   . . . o o . . . .

//   . . o o o o . . .
//   . . o o o . . . .

//   . . o o o . . . .
//   . . o o o o . . .

//   . . . o o o . . .
//   . . o o o o . . .

//   . . o o o o . . .
//   . . . o o o . . .

typedef uint32_t PackedShape;
struct Shape {
	uint32_t roll = 0;
	enum Nibble : uint8_t {
		BackLeft   = (1 << 0),
		BackRight  = (1 << 1),
		FrontLeft  = (1 << 2),
		FrontRight = (1 << 3),
	};
	uint8_t nibbles = 0;

	//------------------

	Shape(uint32_t roll_, uint8_t nibbles_) : roll(roll_), nibbles(nibbles_) {
	}

	//------------------

	static uint32_t count_shapes_for(uint32_t count) {
		if (count % 2 == 0) {
			return 5 * count;
		} else {
			return 4 * count;
		}
	}

	static std::vector< Shape > make_shapes_for(uint32_t count) {
		std::vector< Shape > ret;
		if (count % 2 == 0) {
			ret.reserve(5 * count);
			for (uint32_t i = 0; i < count; ++i) {
				ret.emplace_back(i, 0);
				ret.emplace_back(i, BackLeft | BackRight);
				ret.emplace_back(i, BackLeft | FrontRight);
				ret.emplace_back(i, FrontLeft | BackRight);
				ret.emplace_back(i, FrontLeft | FrontRight);
			}
			assert(ret.size() == count * 5);
		} else {
			ret.reserve(4 * count);
			for (uint32_t i = 0; i < count; ++i) {
				ret.emplace_back(i, BackLeft);
				ret.emplace_back(i, BackRight);
				ret.emplace_back(i, FrontLeft);
				ret.emplace_back(i, FrontRight);
			}
			assert(ret.size() == count * 4);
		}
		return ret;
	}

	//------------------

	PackedShape pack() const {
		return uint32_t(nibbles) | (roll << 4);
	}
	static Shape unpack(PackedShape ps) {
		return Shape(ps >> 4, ps & 0x0f);
	}
};
