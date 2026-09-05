#pragma once
#include "LocalCoordinates3D.h"
#include "Coordinates.h"
#include "NumericalIntegration.h"
#include "LLt.h"

namespace DG
{
	//Element local adjusmnets calculator
	template<class Basis>
	class ElementLAC
	{
		//using Basis = LinearLagrangeBasis;
		using NumericalIntegrationMethod = NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>;
	public:
		static const uint16_t N_BASIS_VALUES = Basis::N_FUNCTIONS * NumericalIntegrationMethod::nSteps;
		static const uint8_t N_LOCAL_MATRIX_ELEMENTS = Basis::N_FUNCTIONS * Basis::N_FUNCTIONS;

		static void init();
		static void finalize();
		static bool isInitialized();

		static const double* getMassMatrix();
		static const double* getMassVector();
		static const double* getIntegrationValues();
		static const LocalCoordinates3D* getIntegrationLocalGradients();

		static void computeStiffnessMatrix(const Coordinates gradients[N_BASIS_VALUES], double stiffnessMatrix[N_LOCAL_MATRIX_ELEMENTS]);
		static void computePowerVector(const double targetFunctionValues[NumericalIntegrationMethod::nSteps], double powerVector[Basis::N_FUNCTIONS]);

		static const LLt* getMassSLAESolverPtr();

	private:
		ElementLAC() = default;

		static bool _initializationState;

		static void computeMassMatrix();
		static void computeMassVector();

		static double* _memoryBuffer;
		static double* _massMatrix;
		static double* _massVector;

		static double* _values;
		static LocalCoordinates3D* _localGradients;

		static LLt* _massSLAESolverPtr;
	};
}


