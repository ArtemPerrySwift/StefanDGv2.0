#include "FaceDGCalculator.h"
#include "CoordinatesFunctions.h"
#include "Bases.h"

namespace DG
{
	template <class Basis>
	void FaceLAC<Basis>::computeFlowMatrix(const uint8_t faceIndex, const double normalDerivatives[N_BASIS_VALUES], double flowMatrix[N_LOCAL_MATRIX_ELEMENTS])
	{
		const double* iFunctionValues = _valuesByFace[faceIndex];
		double* ijElement = flowMatrix;
		for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
		{
			const double* jFunctionNormalDerivatives = normalDerivatives;
			for (uint8_t j = 0; j < Basis::N_FUNCTIONS; ++j)
			{
				*ijElement = NumericalIntegrationMethod::integrateProduct(iFunctionValues, jFunctionNormalDerivatives);
				++ijElement;

				jFunctionNormalDerivatives += NumericalIntegrationMethod::nSteps;
			}
			iFunctionValues += NumericalIntegrationMethod::nSteps;
		}
	}

	template <class Basis>
	void FaceLAC<Basis>::computeFlowVector(const double normalDerivatives[N_BASIS_VALUES], const double targetFunctionValues[NumericalIntegrationMethod::nSteps], double flowVector[Basis::N_FUNCTIONS])
	{
		for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
		{
			*flowVector = NumericalIntegrationMethod::integrateProduct(targetFunctionValues, normalDerivatives);

			normalDerivatives += NumericalIntegrationMethod::nSteps;
			++flowVector;
		}
	}

	template<class Basis>
	void FaceLAC<Basis>::computeFlowVector(const double normalDerivatives[N_BASIS_VALUES], double flowVector[Basis::N_FUNCTIONS])
	{
		for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
		{
			*flowVector = NumericalIntegrationMethod::integrate(normalDerivatives);

			normalDerivatives += NumericalIntegrationMethod::nSteps;
			++flowVector;
		}
	}

	template <class Basis>
	void FaceLAC<Basis>::computePowerVector(const uint8_t faceIndex, const double targetFunctionValues[NumericalIntegrationMethod::nSteps], double powerVector[Basis::N_FUNCTIONS])
	{
		const double* faceValuesIt = _valuesByFace[faceIndex];
		for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
		{
			*powerVector = NumericalIntegrationMethod::integrateProduct(faceValuesIt, targetFunctionValues);
			faceValuesIt += NumericalIntegrationMethod::nSteps;
			++powerVector;
		}
	}

	template<class Basis>
	void FaceLAC<Basis>::changeIntegrationValuesOrder(double values[N_BASIS_VALUES], const uint8_t nodesChangins[constants::triangle::N_NODES])
	{
		NumericalIntegrationMethod::changeIntegrationValuesOrder(values, nodesChangins);
		values += NumericalIntegrationMethod::nSteps;
		NumericalIntegrationMethod::changeIntegrationValuesOrder(values, nodesChangins);
		values += NumericalIntegrationMethod::nSteps;
		NumericalIntegrationMethod::changeIntegrationValuesOrder(values, nodesChangins);
		values += NumericalIntegrationMethod::nSteps;
		NumericalIntegrationMethod::changeIntegrationValuesOrder(values, nodesChangins);
	}

	template<class Basis>
	void FaceLAC<Basis>::changeIntegrationValuesOrder(const double values[N_BASIS_VALUES], const uint8_t nodesChangins[constants::triangle::N_NODES], double outValues[N_BASIS_VALUES])
	{
		NumericalIntegrationMethod::changeIntegrationValuesOrder(values, nodesChangins, outValues);
		values += NumericalIntegrationMethod::nSteps;
		NumericalIntegrationMethod::changeIntegrationValuesOrder(values, nodesChangins, outValues);
		values += NumericalIntegrationMethod::nSteps;
		NumericalIntegrationMethod::changeIntegrationValuesOrder(values, nodesChangins, outValues);
		values += NumericalIntegrationMethod::nSteps;
		NumericalIntegrationMethod::changeIntegrationValuesOrder(values, nodesChangins, outValues);
	}

	template <class Basis>
	const double(*FaceLAC<Basis>::getValuesByFace())[N_BASIS_VALUES]
	{
		return _valuesByFace;
	}

	template <class Basis>
	const LocalCoordinates3D(*FaceLAC<Basis>::getLocalGradientsByFace())[N_BASIS_VALUES]
	{
		return _localGradientsByFace;
	}

	template <class Basis>
	const double* const* FaceLAC<Basis>::getMassMatrixies()
	{
		return _massMatrixByFace;
	}

	template <class Basis>
	const double(*FaceLAC<Basis>::getCrossMassMatrixies())[constants::tetrahedron::N_FACES][N_LOCAL_MATRIX_ELEMENTS]
	{
		return _massCrossMatrixByFaces;
	}

		template<class Basis>
	const double(*FaceLAC<Basis>::getMassVectors())[Basis::N_FUNCTIONS]
	{
		return _massVectorByFace;
	}

	template <class Basis>
	void FaceLAC<Basis>::init()
	{
		_massCrossMatrixByFaces = new double[constants::tetrahedron::N_FACES][constants::tetrahedron::N_FACES][N_LOCAL_MATRIX_ELEMENTS];
		_massMatrixByFace = new double* [constants::tetrahedron::N_FACES];
		_valuesByFace = new double[constants::tetrahedron::N_FACES][N_BASIS_VALUES];
		_localGradientsByFace = new LocalCoordinates3D[constants::tetrahedron::N_FACES][N_BASIS_VALUES];
		_massVectorByFace = new double[constants::tetrahedron::N_FACES][Basis::N_FUNCTIONS];

		computeBasisDataByFace();
		computeMassMatrixies();
		computeMassVectors();

		_initializationState = true;
	}

	template <class Basis>
	void FaceLAC<Basis>::finalize()
	{
		delete[] _massCrossMatrixByFaces;
		delete[] _massMatrixByFace;
		delete[] _valuesByFace;
		delete[] _localGradientsByFace;
		delete[] _massVectorByFace;

		_massCrossMatrixByFaces = nullptr;
		_massMatrixByFace = nullptr;
		_valuesByFace = nullptr;
		_localGradientsByFace = nullptr;
		_massVectorByFace = nullptr;

		_initializationState;
	}

	template <class Basis>
	bool FaceLAC<Basis>::isInitialized()
	{
		return _initializationState;
	}

	template <class Basis>
	bool FaceLAC<Basis>::_initializationState = false;
	template <class Basis>
	double (*FaceLAC<Basis>::_massCrossMatrixByFaces)[constants::tetrahedron::N_FACES][N_LOCAL_MATRIX_ELEMENTS] = nullptr;
	template <class Basis>
	double** FaceLAC<Basis>::_massMatrixByFace = nullptr;
	template <class Basis>
	double (*FaceLAC<Basis>::_massVectorByFace)[Basis::N_FUNCTIONS] = nullptr;
	template <class Basis>
	double (*FaceLAC<Basis>::_valuesByFace)[N_BASIS_VALUES] = nullptr;
	template <class Basis>
	LocalCoordinates3D(*FaceLAC<Basis>::_localGradientsByFace)[N_BASIS_VALUES] = nullptr;

	template <class Basis>
	void FaceLAC<Basis>::computeBasisDataByFace()
	{
		LocalCoordinates3D localCoordinates3D[NumericalIntegrationMethod::nSteps];

		CoordinatesFunctions::translate0FacePoints(NumericalIntegrationMethod::localPoints, NumericalIntegrationMethod::nSteps, localCoordinates3D);
		Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _valuesByFace[0]);
		Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _localGradientsByFace[0]);

		CoordinatesFunctions::translate1FacePoints(NumericalIntegrationMethod::localPoints, NumericalIntegrationMethod::nSteps, localCoordinates3D);
		Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _valuesByFace[1]);
		Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _localGradientsByFace[1]);

		CoordinatesFunctions::translate2FacePoints(NumericalIntegrationMethod::localPoints, NumericalIntegrationMethod::nSteps, localCoordinates3D);
		Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _valuesByFace[2]);
		Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _localGradientsByFace[2]);

		CoordinatesFunctions::translate3FacePoints(NumericalIntegrationMethod::localPoints, NumericalIntegrationMethod::nSteps, localCoordinates3D);
		Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _valuesByFace[3]);
		Basis::compute(localCoordinates3D, NumericalIntegrationMethod::nSteps, _localGradientsByFace[3]);

	}

	template <class Basis>
	void FaceLAC<Basis>::computeMassMatrix(const double values[N_BASIS_VALUES], double massMatrix[N_LOCAL_MATRIX_ELEMENTS])
	{
		double* di = massMatrix;
		const double* iFunctionValues = values;
		for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
		{
			*di = NumericalIntegrationMethod::integrateProduct(iFunctionValues, iFunctionValues);
			const double* jFunctionValues = iFunctionValues;
			double* ijElement = di;
			double* jiElement = di;
			for (uint8_t j = i + 1; j < Basis::N_FUNCTIONS; ++j)
			{
				++ijElement;
				jiElement += Basis::N_FUNCTIONS;
				jFunctionValues += NumericalIntegrationMethod::nSteps;

				*ijElement = *jiElement = NumericalIntegrationMethod::integrateProduct(iFunctionValues, jFunctionValues);

			}

			iFunctionValues += NumericalIntegrationMethod::nSteps;
			di += Basis::N_FUNCTIONS + 1;
		}
	}

	template <class Basis>
	void FaceLAC<Basis>::computeCrossMassMatrix(const double values1[N_BASIS_VALUES], const double values2[N_BASIS_VALUES], double crossMassMatrix1[N_LOCAL_MATRIX_ELEMENTS], double crossMassMatrix2[N_LOCAL_MATRIX_ELEMENTS])
	{
		const double* iFunctionValues = values1;
		double* ijElementMatrix1 = crossMassMatrix1;
		double* matrix2InitRowElementIt = crossMassMatrix2;

		for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
		{
			const double* jFunctionValues = values2;
			double* jiElementMatrix2 = matrix2InitRowElementIt;
			for (uint8_t j = 0; j < Basis::N_FUNCTIONS; ++j)
			{
				*ijElementMatrix1 = *jiElementMatrix2 = NumericalIntegrationMethod::integrateProduct(iFunctionValues, jFunctionValues);
				++ijElementMatrix1;

				jiElementMatrix2 += Basis::N_FUNCTIONS;
				jFunctionValues += NumericalIntegrationMethod::nSteps;
			}
			iFunctionValues += NumericalIntegrationMethod::nSteps;
			++matrix2InitRowElementIt;
		}
	}

	template <class Basis>
	void FaceLAC<Basis>::computeMassMatrixies()
	{
		double** massMetrixIt = _massMatrixByFace;
		double (*iiMassMatrixPtr)[N_LOCAL_MATRIX_ELEMENTS] = *_massCrossMatrixByFaces;
		const double (*iFaceValuesIt)[N_BASIS_VALUES] = _valuesByFace;

		for (uint8_t i = 0; i < constants::tetrahedron::N_FACES; ++i)
		{
			*massMetrixIt = *iiMassMatrixPtr;
			computeMassMatrix(*iFaceValuesIt, *iiMassMatrixPtr);
			double (*ijMassMatrixIt)[N_LOCAL_MATRIX_ELEMENTS] = iiMassMatrixPtr;
			double (*jiMassMatrixIt)[N_LOCAL_MATRIX_ELEMENTS] = iiMassMatrixPtr;

			const double (*jFaceValuesIt)[N_BASIS_VALUES] = iFaceValuesIt;

			for (uint8_t j = i + 1; j < constants::tetrahedron::N_FACES; ++j)
			{
				++jFaceValuesIt;
				++ijMassMatrixIt;
				jiMassMatrixIt += constants::tetrahedron::N_FACES;

				computeCrossMassMatrix(*iFaceValuesIt, *jFaceValuesIt, *ijMassMatrixIt, *jiMassMatrixIt);
			}

			iiMassMatrixPtr += constants::tetrahedron::N_FACES + 1;
			++massMetrixIt;
		}
	}

	template <class Basis>
	void FaceLAC<Basis>::computeMassVector(const double* faceValueIt, double* massVectorElementIt)
	{
		for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
		{
			*massVectorElementIt = NumericalIntegrationMethod::integrate(faceValueIt);

			faceValueIt += NumericalIntegrationMethod::nSteps;
			++massVectorElementIt;
		}
	}

	template <class Basis>
	void FaceLAC<Basis>::computeMassVectors()
	{
		double (*faceValuesIt)[N_BASIS_VALUES] = _valuesByFace;
		double (*massVectorIt)[Basis::N_FUNCTIONS] = _massVectorByFace;

		for (uint8_t i = 0; i < constants::tetrahedron::N_FACES; ++i)
		{
			NumericalIntegrationMethod::integrateProduct(*faceValuesIt, *massVectorIt);

			++faceValuesIt;
			++massVectorIt;
		}
	}
}

#define X(BasisName) template class DG::FaceLAC<BasisName>;
	BASES
#undef X;
