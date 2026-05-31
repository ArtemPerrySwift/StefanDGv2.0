#pragma once
#include "Coordinates.h"

struct Interface
{
	int tag;
	bool isPlane;

	unsigned int regionsIndexes[2];

	struct FacesSet
	{
		size_t* tags;
		size_t count;

		double* determinnats;
		Coordinates* normals;
	};

};