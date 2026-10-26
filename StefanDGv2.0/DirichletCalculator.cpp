#include "NumericalIntegration.h"
#include "Bases.h"
#include "FaceDGCalculator.h"
#include "CoordinatesFunctions.h"
#include "Boundary.h"
#include "NodesFunctions.h"
#include <cmath>

//#define VALUE_CONDITION;

//#define ON_PLANE
//#define SHARED_BOUNDARY

namespace DG
{
    using Basis = LinearLagrangeBasis;

    template<size_t N_FUNCTIONS>
    static void addFaceLinearAdjusments(const double* powerVectorIt,
        const double* flowVectorElementIt,
        const double penaltyCoeff,
        const double commonCoeff,
        double* adjusmentIt)
    {
        for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
        {
            *adjusmentIt += commonCoeff * (penaltyCoeff * (*powerVectorIt) - *flowVectorElementIt);

            ++powerVectorIt;
            ++flowVectorElementIt;
            ++adjusmentIt;
        }
    }


    template<size_t N_FUNCTIONS>
    static void addFaceBilinearAdjusmets(const double* flowMatrixElementIt,
        const double* massMatrixElementIt,
        const double penaltyCoeff,
        const double detLambdaCoeff,
        double* adjusmentIt)
    {
        double* jiAdjusmentIt;
        double* ijAdjusmentIt;
        for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
        {
            massMatrixElementIt += i;
            *adjusmentIt += detLambdaCoeff * (-*flowMatrixElementIt - *flowMatrixElementIt + penaltyCoeff * (*massMatrixElementIt));
            jiAdjusmentIt = adjusmentIt;
            ijAdjusmentIt = adjusmentIt;

            const double* ijFlowMatrixElementPtr = flowMatrixElementIt;
            const double* jiFlowMatrixElementPtr = flowMatrixElementIt;
            for (uint8_t j = i + 1; j < N_FUNCTIONS; ++j)
            {
                ++massMatrixElementIt;

                ++ijAdjusmentIt;
                jiAdjusmentIt += N_FUNCTIONS;

                ++ijFlowMatrixElementPtr;
                jiFlowMatrixElementPtr += N_FUNCTIONS;

                double adjusment = detLambdaCoeff * (-*ijFlowMatrixElementPtr - *jiFlowMatrixElementPtr + penaltyCoeff * (*massMatrixElementIt));
                *ijAdjusmentIt += adjusment;
                *jiAdjusmentIt += adjusment;
            }

            flowMatrixElementIt += N_FUNCTIONS + 1;
            adjusmentIt += N_FUNCTIONS + 1;
            ++massMatrixElementIt;
        }

        adjusmentIt = nullptr;
    }


#ifdef VALUE_CONDITION
#ifdef NO_REAPEAT_VALUE_CONDITION
#define NO_REPEAT_FACE_ADJUSMENTS
#endif
#else
#ifdef NO_REPEAT_FUNCTION_CONDITION
#define NO_REPEAT_FACE_ADJUSMENTS
#endif
#endif

#ifndef NO_REPEAT_FACE_ADJUSMENTS
    static void addFaceAdjusments(const uint8_t localIndex,
                                  const double lambda,
                                  const double penalty,
                                  const double det,
                                  const Coordinates& unitNormal,
#ifdef VALUE_CONDITION
                                  const double conditionValue,
#else
                                  double conditionValues[NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps],
#endif
                                  const double transpLocalJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                  void* memoryBuffer,
                                  double bilinearAdjusments[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                  double linearAdjusments[Basis::N_FUNCTIONS])
    {
        double* flowMatrix = (double*)memoryBuffer;

        double* flowVector = flowMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
        double* normalDerivatives = flowVector + Basis::N_FUNCTIONS;

        Coordinates* gradients = (Coordinates*)(normalDerivatives + FaceLAC<Basis>::N_BASIS_VALUES);

#ifndef VALUE_CONDITION
        double* powerVector = (double*)(gradients + FaceLAC<Basis>::N_BASIS_VALUES);
#endif

        double commonMultiplier = det * lambda;

        CoordinatesFunctions::translate(transpLocalJacobianMatrix, FaceLAC<Basis>::getLocalGradientsByFace()[localIndex], FaceLAC<Basis>::N_BASIS_VALUES, gradients);

        CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, unitNormal, normalDerivatives);

        FaceLAC<Basis>::computeFlowMatrix(localIndex, normalDerivatives, flowMatrix);
        addFaceBilinearAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, FaceLAC<Basis>::getMassMatrixies()[localIndex], penalty, commonMultiplier, bilinearAdjusments);

        FaceLAC<Basis>::computeFlowVector(normalDerivatives, flowVector);

#ifdef VALUE_CONDITION
        addFaceLinearAdjusments<Basis::N_FUNCTIONS>(FaceLAC<Basis>::getMassVectors()[localIndex], flowVector, penalty, commonMultiplier * conditionValue, linearAdjusments);
#else
        FaceLAC<Basis>::computeFlowVector(normalDerivatives, conditionValues, flowVector);
        FaceLAC<Basis>::computePowerVector(localIndex, conditionValues, powerVector);

        addFaceLinearAdjusments<Basis::N_FUNCTIONS>(powerVector, flowVector, penalty, commonMultiplier, linearAdjusments);
#endif
    }
#ifdef VALUE_CONDITION
#define NO_REAPEAT_VALUE_CONDITION
#else
#define NO_REPEAT_FUNCTION_CONDITION
#endif
#endif

#ifdef NO_REPEAT_FACE_ADJUSMENTS
#undef NO_REPEAT_FACE_ADJUSMENTS
#endif


#ifdef ON_PLANE
    static void computeDirichletFacesAdjusmentsOnPlane(
#else
    static void computeDirichletFacesAdjusments(
#endif
        const size_t* baseFaceIndexesIt,
#ifdef SHARED_BOUNDARY
        const size_t* neigbourFaceIndexesIt,
#endif
        const size_t nFaces,
        const double lambda,
        const double penalty,
#ifdef VALUE_CONDITION
        const double conditionValue,
#endif
        const Coordinates nodes[],
        const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
        const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
        void* memoryBuffer,
        double baseBilinearTetrahedronsAdjusments[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
        double baseLinearTetrahedronsAdjusments[][Basis::N_FUNCTIONS]
#ifdef SHARED_BOUNDARY
        ,
        double neigbourBilinearTetrahedronsAdjusments[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
        double neigbourLinearTetrahedronsAdjusments[][Basis::N_FUNCTIONS]
#endif
    )
    {

        double* dirichletValues = (double*)memoryBuffer;

        memoryBuffer = (void*)(dirichletValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps);

        Coordinates normal, unitNormal, reverseUnitNormal;
        size_t elementIndex;
        uint8_t localIndex;

        size_t baseFaceNodesTags[constants::triangle::N_NODES];
#ifdef SHARED_BOUNDARY
        size_t neighbourFaceNodesTags[constants::triangle::N_NODES];
#endif
        double transpJacobian[LocalCoordinates2D::COUNT * Coordinates::COUNT];

#ifdef ON_PLANE
        elementIndex = *baseFaceIndexesIt >> 2;
        localIndex = *baseFaceIndexesIt & 3;

        size_t faceNodesTags[constants::triangle::N_NODES];
        NodesFunctions::extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], baseFaceNodesTags);
        CoordinatesFunctions::computeTranspJacobianToFace(nodes, baseFaceNodesTags, transpJacobian);
        CoordinatesFunctions::computeNormal(transpJacobian, normal);

        double det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
        unitNormal.x = normal.x / det;
        unitNormal.y = normal.y / det;
        unitNormal.z = normal.z / det;

        reverseUnitNormal.x = -unitNormal.x;
        reverseUnitNormal.y = -unitNormal.y;
        reverseUnitNormal.z = -unitNormal.z;
#endif
        for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
        {
            //double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

            elementIndex = *baseFaceIndexesIt >> 2;
            localIndex = *baseFaceIndexesIt & 3;

            NodesFunctions::extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], baseFaceNodesTags);

#ifdef VALUE_CONDITION
            CoordinatesFunctions::computeNormal(nodes, baseFaceNodesTags, normal);
#else
            CoordinatesFunctions::computeTranspJacobianToFace(nodes, baseFaceNodesTags, transpJacobian);
            Boundary::FunctionCondition::compute<Boundary::FunctionCondition::Type::DIRICHLET, NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps>(
                NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                transpJacobian,
                nodes[*baseFaceNodesTags],
                dirichletValues);
            CoordinatesFunctions::computeNormal(transpJacobian, normal);
#endif

            double det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

#ifndef ON_PLANE
            unitNormal.x = normal.x / det;
            unitNormal.y = normal.y / det;
            unitNormal.z = normal.z / det;
#endif

            addFaceAdjusments(localIndex,
                lambda,
                penalty,
                det,
                normal,
#ifdef VALUE_CONDITION
                conditionValue,
#else
                dirichletValues,
#endif
                transpLocalJacobianMatrix[elementIndex],
                memoryBuffer,
                baseBilinearTetrahedronsAdjusments[elementIndex],
                baseBilinearTetrahedronsAdjusments[elementIndex]);


#ifdef SHARED_BOUNDARY
            elementIndex = *neigbourFaceIndexesIt >> 2;
            localIndex = *neigbourFaceIndexesIt & 3;

#ifndef ON_PLANE
            reverseUnitNormal.x = -unitNormal.x;
            reverseUnitNormal.y = -unitNormal.y;
            reverseUnitNormal.z = -unitNormal.z;
#endif

            NodesFunctions::extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], neighbourFaceNodesTags);

            uint8_t changingIndexes = NodesFunctions::computeChangingIndexes(baseFaceNodesTags, neighbourFaceNodesTags);
            if (changingIndexes != 0)
            {
                NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::changeValuesOrder(dirichletValues, Basis::N_FUNCTIONS, changingIndexes);
            }

            addFaceAdjusments(localIndex,
                lambda,
                penalty,
                det,
                normal,
#ifdef VALUE_CONDITION
                conditionValue,
#else
                dirichletValues,
#endif
                transpLocalJacobianMatrix[elementIndex],
                memoryBuffer,
                neigbourBilinearTetrahedronsAdjusments[elementIndex],
                neigbourBilinearTetrahedronsAdjusments[elementIndex]);

            ++neigbourFaceIndexesIt;
#endif
            ++baseFaceIndexesIt;
        }
    }
}
