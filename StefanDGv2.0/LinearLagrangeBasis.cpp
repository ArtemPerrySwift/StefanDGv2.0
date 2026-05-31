#include "LinearLagrangeBasis.h"

double LinearLagrangeBasis::compute(const LocalCoordinates3D& localPoint, const double* dofs)
{
	double startDOF = dofs[0];
	double sum = startDOF;
	sum += localPoint.u * (dofs[1] - startDOF);
	sum += localPoint.v * (dofs[2] - startDOF);
	sum += localPoint.w * (dofs[3] - startDOF);
	return sum;
}

void LinearLagrangeBasis::compute(const LocalCoordinates3D& localPoint, const double* dofs, LocalCoordinates3D& localGradient)
{
	localGradient.u = -dofs[0] + dofs[1];
	localGradient.v = -dofs[0] + dofs[2];
	localGradient.w = -dofs[0] + dofs[3];
}

void LinearLagrangeBasis::compute(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, double* valueIt)
{
	double* function1ValueIt = valueIt + nPoints;
	double* function2ValueIt = function1ValueIt + nPoints;
	double* function3ValueIt = function2ValueIt + nPoints;

	for (uint8_t i = 0; i < nPoints; ++i)
	{
		*valueIt = 1.0 - localPointIt->u - localPointIt->v - localPointIt->w;
		++valueIt;
		*function1ValueIt = localPointIt->u;
		++function1ValueIt;
		*function2ValueIt = localPointIt->v;
		++function2ValueIt;
		*function3ValueIt = localPointIt->w;
		++function3ValueIt;

		++localPointIt;
	}
}


void LinearLagrangeBasis::compute(const LocalCoordinates3D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt)
{
	LocalCoordinates3D* function1GradientIt = gradientIt + nPoints;
	LocalCoordinates3D* function2GradientIt = function1GradientIt + nPoints;
	LocalCoordinates3D* function3GradientIt = function2GradientIt + nPoints;

	for (uint8_t i = 0; i < nPoints; ++i)
	{
		gradientIt->u = -1.0;
		gradientIt->v = -1.0;
		gradientIt->w = -1.0;
		++gradientIt;

		function1GradientIt->u = 1.0;
		function1GradientIt->v = 0.0;
		function1GradientIt->w = 0.0;
		++function1GradientIt;

		function2GradientIt->u = 0.0;
		function2GradientIt->v = 1.0;
		function2GradientIt->w = 0.0;
		++function2GradientIt;


		function3GradientIt->u = 0.0;
		function3GradientIt->v = 0.0;
		function3GradientIt->w = 1.0;
		++function3GradientIt;

		++localPointIt;
	}
}

// u = v, v = u, w = 0
void LinearLagrangeBasis::computeOnFace0(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, double* valueIt)
{
	double* function1ValueIt = valueIt + nPoints;
	double* function2ValueIt = function1ValueIt + nPoints;
	double* function3ValueIt = function2ValueIt + nPoints;

	for (uint8_t i = 0; i < nPoints; ++i)
	{
		*valueIt = 1.0 - localPointIt->u - localPointIt->v;
		++valueIt;
		*function1ValueIt = localPointIt->v;
		++function1ValueIt;
		*function2ValueIt = localPointIt->u;
		++function2ValueIt;
		*function3ValueIt = 0.0;
		++function3ValueIt;

		++localPointIt;
	}
}

// u = u, v = 0, w = v
void LinearLagrangeBasis::computeOnFace1(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, double* valueIt)
{
	double* function1ValueIt = valueIt + nPoints;
	double* function2ValueIt = function1ValueIt + nPoints;
	double* function3ValueIt = function2ValueIt + nPoints;

	for (uint8_t i = 0; i < nPoints; ++i)
	{
		*valueIt = 1.0 - localPointIt->u - localPointIt->v;
		++valueIt;
		*function1ValueIt = localPointIt->u;
		++function1ValueIt;
		*function2ValueIt = 0.0;
		++function2ValueIt;
		*function3ValueIt = localPointIt->v;
		++function3ValueIt;

		++localPointIt;
	}
}

// u = 0, v = v, w = u
void LinearLagrangeBasis::computeOnFace2(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, double* valueIt)
{
	double* function1ValueIt = valueIt + nPoints;
	double* function2ValueIt = function1ValueIt + nPoints;
	double* function3ValueIt = function2ValueIt + nPoints;

	for (uint8_t i = 0; i < nPoints; ++i)
	{
		*valueIt = 1.0 - localPointIt->u - localPointIt->v;
		++valueIt;
		*function1ValueIt = 0.0;
		++function1ValueIt;
		*function2ValueIt = localPointIt->v;
		++function2ValueIt;
		*function3ValueIt = localPointIt->u;
		++function3ValueIt;

		++localPointIt;
	}
}

// u = u, v = v, w = 1.0 - u - v
void LinearLagrangeBasis::computeOnFace3(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, double* valueIt)
{
	double* function1ValueIt = valueIt + nPoints;
	double* function2ValueIt = function1ValueIt + nPoints;
	double* function3ValueIt = function2ValueIt + nPoints;

	for (uint8_t i = 0; i < nPoints; ++i)
	{
		*valueIt = 0.0;
		++valueIt;
		*function1ValueIt = localPointIt->u;
		++function1ValueIt;
		*function2ValueIt = localPointIt->v;
		++function2ValueIt;
		*function3ValueIt = 1.0 - localPointIt->u - localPointIt->v;
		++function3ValueIt;

		++localPointIt;
	}
}

void LinearLagrangeBasis::computeOnFace0(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt)
{
	LocalCoordinates3D* function1GradientIt = gradientIt + nPoints;
	LocalCoordinates3D* function2GradientIt = function1GradientIt + nPoints;
	LocalCoordinates3D* function3GradientIt = function2GradientIt + nPoints;

	for (uint8_t i = 0; i < nPoints; ++i)
	{
		gradientIt->u = -1.0;
		gradientIt->v = -1.0;
		gradientIt->w = -1.0;
		++gradientIt;

		function1GradientIt->u = 1.0;
		function1GradientIt->v = 0.0;
		function1GradientIt->w = 0.0;
		++function1GradientIt;

		function2GradientIt->u = 0.0;
		function2GradientIt->v = 1.0;
		function2GradientIt->w = 0.0;
		++function2GradientIt;


		function3GradientIt->u = 0.0;
		function3GradientIt->v = 0.0;
		function3GradientIt->w = 1.0;
		++function3GradientIt;

		++localPointIt;
	}
}

void LinearLagrangeBasis::computeOnFace1(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt)
{
	LocalCoordinates3D* function1GradientIt = gradientIt + nPoints;
	LocalCoordinates3D* function2GradientIt = function1GradientIt + nPoints;
	LocalCoordinates3D* function3GradientIt = function2GradientIt + nPoints;

	for (uint8_t i = 0; i < nPoints; ++i)
	{
		gradientIt->u = -1.0;
		gradientIt->v = -1.0;
		gradientIt->w = -1.0;
		++gradientIt;

		function1GradientIt->u = 1.0;
		function1GradientIt->v = 0.0;
		function1GradientIt->w = 0.0;
		++function1GradientIt;

		function2GradientIt->u = 0.0;
		function2GradientIt->v = 1.0;
		function2GradientIt->w = 0.0;
		++function2GradientIt;


		function3GradientIt->u = 0.0;
		function3GradientIt->v = 0.0;
		function3GradientIt->w = 1.0;
		++function3GradientIt;

		++localPointIt;
	}
}

void LinearLagrangeBasis::computeOnFace2(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt)
{
	LocalCoordinates3D* function1GradientIt = gradientIt + nPoints;
	LocalCoordinates3D* function2GradientIt = function1GradientIt + nPoints;
	LocalCoordinates3D* function3GradientIt = function2GradientIt + nPoints;

	for (uint8_t i = 0; i < nPoints; ++i)
	{
		gradientIt->u = -1.0;
		gradientIt->v = -1.0;
		gradientIt->w = -1.0;
		++gradientIt;

		function1GradientIt->u = 1.0;
		function1GradientIt->v = 0.0;
		function1GradientIt->w = 0.0;
		++function1GradientIt;

		function2GradientIt->u = 0.0;
		function2GradientIt->v = 1.0;
		function2GradientIt->w = 0.0;
		++function2GradientIt;


		function3GradientIt->u = 0.0;
		function3GradientIt->v = 0.0;
		function3GradientIt->w = 1.0;
		++function3GradientIt;

		++localPointIt;
	}
}

void LinearLagrangeBasis::computeOnFace3(const LocalCoordinates2D* localPointIt, const uint8_t nPoints, LocalCoordinates3D* gradientIt)
{
	LocalCoordinates3D* function1GradientIt = gradientIt + nPoints;
	LocalCoordinates3D* function2GradientIt = function1GradientIt + nPoints;
	LocalCoordinates3D* function3GradientIt = function2GradientIt + nPoints;

	for (uint8_t i = 0; i < nPoints; ++i)
	{
		gradientIt->u = -1.0;
		gradientIt->v = -1.0;
		gradientIt->w = -1.0;
		++gradientIt;

		function1GradientIt->u = 1.0;
		function1GradientIt->v = 0.0;
		function1GradientIt->w = 0.0;
		++function1GradientIt;

		function2GradientIt->u = 0.0;
		function2GradientIt->v = 1.0;
		function2GradientIt->w = 0.0;
		++function2GradientIt;


		function3GradientIt->u = 0.0;
		function3GradientIt->v = 0.0;
		function3GradientIt->w = 1.0;
		++function3GradientIt;

		++localPointIt;
	}
}


