#pragma once
#include <stdio.h>
#include <cstdint>

struct MaterialPhase
{
	enum State : uint8_t { SOLID, LIQUID };

	State state;
	double thermalConductivity;
	double heatCapacity;
	double density;
};