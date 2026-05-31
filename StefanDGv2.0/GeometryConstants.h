#pragma once
#include <cstdint>

namespace constants
{
	const uint8_t DIMENSION_0D = 0;
	const uint8_t DIMENSION_1D = 1;
	const uint8_t DIMENSION_2D = 2;
	const uint8_t DIMENSION_3D = 1;

	namespace tetrahedron
	{
		const uint8_t N_FACES = 4;
		const uint8_t N_NODES = 4;
	}

	namespace triangle
	{
		const uint8_t N_NODES = 3;
		const uint8_t N_EDGES = 3;
	}

}
