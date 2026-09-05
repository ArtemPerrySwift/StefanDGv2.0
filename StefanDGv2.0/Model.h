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

	unsigned int nRegions;
	unsigned int nBoundaries;
	unsigned int nNonconformInterfaces;
	
	Error initilizeByCurrentGMSHModel();
	void clear();
};

