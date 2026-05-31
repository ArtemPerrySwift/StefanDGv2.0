#pragma once
#include <cstdint>
#include "LocalCoordinates3D.h"
#include "LocalCoordinates2D.h"

class LinearLagrangeBasis
{
public:
	static const uint8_t N_FUNCTIONS = 4;
	static const uint8_t ORDER = 1;

	static double compute(const LocalCoordinates3D& localPoint, const double* dofs);
	static void compute(const LocalCoordinates3D& localPoint, const double* dofs, LocalCoordinates3D& localGradient);

	static void compute(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, double* valueIt);
	static void compute(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt);

	void computeOnFace0(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, double* valueIt);
	void computeOnFace1(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, double* valueIt);
	void computeOnFace2(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, double* valueIt);
	void computeOnFace3(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, double* valueIt);

	void computeOnFace0(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt);
	void computeOnFace1(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt);
	void computeOnFace2(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt);
	void computeOnFace3(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt);
	//static void compute(const LocalCoordinates* localPointIt, const size_t nPoints, double* gradientIt);
};

//const uint8_t LinearLagrangeBasis::N_FUNCTIONS = 4;
//const uint8_t LinearLagrangeBasis::ORDER = 1;

