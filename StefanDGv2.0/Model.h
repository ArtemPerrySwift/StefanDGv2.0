#pragma once
#include "MaterialPhase.h"
#include "Region.h"
#include "Boundary.h"
#include "Interface.h"
#include "NonconformInterface.h"

#include <utility>

struct Model
{
public:
	enum class Error : uint8_t {NO_ERRORS, E_BOUNDARIES, E_CONDITIONS, E_MATERIAL_PHASES, E_REGIONS};
	Boundary::ValueCondition* valueConditions;
	Boundary::ConformCondition* conformConditions;
	Boundary* boundaries;
	NonconformInterface* nonconformInterfaces;

	MaterialPhase* materialPhases;
	const MaterialPhase** regionsMaterialPhases;

	int *regionsTags;
	int frontTag;

	unsigned int nRegions;
	unsigned int nBoundaries;
	unsigned int nNonconformInterfaces;
	
	static void movePoints(const int* pointTagIt,
						   const unsigned int nPoints,
						   const Coordinates* pointIt,
						   const Coordinates* movementIt,
						   Coordinates* movedPointIt);

	static void getFronts(int*& frontsPointsTags, Coordinates*& frontsPoints, unsigned int*& frontsPointsStartIndexes, unsigned int& nFronts);
	Error initilizeByCurrentGMSHModel();
	void clear();
};

