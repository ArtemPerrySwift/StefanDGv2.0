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
	Boundary::ValueCondition* valueConditions;
	Boundary::ConformCondition* conformConditions;
	Boundary* boundaries;
	NonconformInterface* nonconformInterfaces;

	MaterialPhase* materialPhases;
	const MaterialPhase** regionsMaterialPhases;

	unsigned int nRegions;
	unsigned int nBoundaries;
	
	bool initilizeByCurrentGMSHModel();
	void clear();
};

