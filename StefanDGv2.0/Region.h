#pragma once
#include "MaterialPhase.h"
#include "Boundary.h"

struct Region
{
	int tag;
	const MaterialPhase* materialPhasePtr;
	const Boundary* boundaries;
	unsigned int nBoundaries;
};