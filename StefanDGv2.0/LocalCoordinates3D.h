#pragma once
#include <cstdint>

struct LocalCoordinates3D
{
	static const uint8_t COUNT = 3;
	double u, v, w;
};