#include "ElementLAC.h"
#include "Bases.h"

namespace DG
{
	template<class Basis>
	void ElementLAC<Basis>::computeMassMatrix()
	{
		double* di = _massMatrix;
		const double* iFunctionValues = _values;
		for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
		{
			*di = NumericalIntegrationMethod::integrateProduct(iFunctionValues, iFunctionValues);
			const double* jFunctionValues = iFunctionValues + NumericalIntegrationMethod::nSteps;
			double* ijElement = di + 1;
			double* jiElement = di + NumericalIntegrationMethod::nSteps;
			for (uint8_t j = i + 1; j < Basis::N_FUNCTIONS; ++j)
			{
				*ijElement = *jiElement = NumericalIntegrationMethod::integrateProduct(iFunctionValues, jFunctionValues);
				++ijElement;

				jiElement += Basis::N_FUNCTIONS;
				jFunctionValues += NumericalIntegrationMethod::nSteps;
			}

			iFunctionValues += NumericalIntegrationMethod::nSteps;
			di += Basis::N_FUNCTIONS + 1;
		}
	}

	template<class Basis>
	void ElementLAC<Basis>::computeMassVector()
	{
		double* massVectorElementIt = _massVector;
		const double* iBasisFunctionValues = _values;
		for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
		{
			*massVectorElementIt = NumericalIntegrationMethod::integrate(iBasisFunctionValues);

			iBasisFunctionValues += NumericalIntegrationMethod::nSteps;
			++massVectorElementIt;
		}
	}

	template<class Basis>
	void ElementLAC<Basis>::init()
	{
		_memoryBuffer = new double[N_LOCAL_MATRIX_ELEMENTS + Basis::N_FUNCTIONS + N_BASIS_VALUES * (1 + LocalCoordinates3D::COUNT)];
		_massMatrix = _memoryBuffer;
		_massVector = _massMatrix + N_LOCAL_MATRIX_ELEMENTS;
		_values = _massVector + Basis::N_FUNCTIONS;
		_localGradients = (LocalCoordinates3D*)(_values + N_BASIS_VALUES);

		Basis::compute(NumericalIntegrationMethod::localPoints, NumericalIntegrationMethod::nSteps, _values);
		Basis::compute(NumericalIntegrationMethod::localPoints, NumericalIntegrationMethod::nSteps, _localGradients);

		computeMassMatrix();
		computeMassVector();

		_massSLAESolverPtr = new LLt(_massMatrix, Basis::N_FUNCTIONS);

		_initializationState = true;
	}

	template<class Basis>
	void ElementLAC<Basis>::finalize()
	{
		delete[] _memoryBuffer;
		_memoryBuffer = nullptr;
		_massMatrix = nullptr;
		_massVector = nullptr;

		_values = nullptr;
		_localGradients = nullptr;

		delete _massSLAESolverPtr;
		_massSLAESolverPtr = nullptr;

		_initializationState = false;
	}

	template<class Basis>
	bool ElementLAC<Basis>::isInitialized()
	{
		return _initializationState;
	}

	template<class Basis>
	const double* ElementLAC<Basis>::getMassMatrix()
	{
		return _massMatrix;
	}

	template<class Basis>
	const double* ElementLAC<Basis>::getMassVector()
	{
		return _massVector;
	}

	template<class Basis>
	const double* ElementLAC<Basis>::getIntegrationValues()
	{
		return _values;
	}

	template<class Basis>
	const LocalCoordinates3D* ElementLAC<Basis>::getIntegrationLocalGradients()
	{
		return _localGradients;
	}

	template<class Basis>
	void ElementLAC<Basis>::computeStiffnessMatrix(const Coordinates gradients[N_BASIS_VALUES], double stiffnessMatrix[N_LOCAL_MATRIX_ELEMENTS])
	{
		double* di = stiffnessMatrix;
		const Coordinates* iFunctionGradients = gradients;
		for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
		{
			*di = NumericalIntegrationMethod::integrateProduct(iFunctionGradients, iFunctionGradients);
			const Coordinates* jFunctionGradients = iFunctionGradients + NumericalIntegrationMethod::nSteps;
			double* ijElement = di + 1;
			double* jiElement = di + NumericalIntegrationMethod::nSteps;
			for (uint8_t j = i + 1; j < Basis::N_FUNCTIONS; ++j)
			{
				*ijElement = *jiElement = NumericalIntegrationMethod::integrateProduct(iFunctionGradients, jFunctionGradients);

				++ijElement;
				jiElement += Basis::N_FUNCTIONS;
				jFunctionGradients += NumericalIntegrationMethod::nSteps;
			}

			iFunctionGradients += NumericalIntegrationMethod::nSteps;
			di += Basis::N_FUNCTIONS + 1;
		}
	}

	template<class Basis>
	void ElementLAC<Basis>::computePowerVector(const double targetFunctionValues[NumericalIntegrationMethod::nSteps], double powerVector[Basis::N_FUNCTIONS])
	{
		const double* iBasisFunctionValues = _values;
		for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
		{
			*powerVector = NumericalIntegrationMethod::integrateProduct(iBasisFunctionValues, targetFunctionValues);

			iBasisFunctionValues += NumericalIntegrationMethod::nSteps;
			++powerVector;
		}
	}

	template<class Basis>
	const LLt* ElementLAC<Basis>::getMassSLAESolverPtr()
	{
		return _massSLAESolverPtr;
	}

	template<class Basis>
	bool ElementLAC<Basis>::_initializationState = false;

	template<class Basis>
	double* ElementLAC<Basis>::_memoryBuffer = nullptr;

	template<class Basis>
	double* ElementLAC<Basis>::_massMatrix = nullptr;

	template<class Basis>
	double* ElementLAC<Basis>::_massVector = nullptr;

	template<class Basis>
	double* ElementLAC<Basis>::_values = nullptr;

	template<class Basis>
	LocalCoordinates3D* ElementLAC<Basis>::_localGradients = nullptr;

	template<class Basis>
	LLt* ElementLAC<Basis>::_massSLAESolverPtr = nullptr;

#define X(BasisName) template class ElementLAC<BasisName>;
	BASES
#undef X;

}
