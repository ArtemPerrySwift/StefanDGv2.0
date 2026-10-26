#include "NumericalIntegration.h"
#include "Bases.h"
#include "FaceDGCalculator.h"
#include "CoordinatesFunctions.h"
#include "Boundary.h"
#include "NodesFunctions.h"
#include <cmath>

//#define SHARED_BOUNDARY
namespace DG
{
    using Basis = LinearLagrangeBasis;

    template<size_t N_FUNCTIONS>
    static void addFaceAdjusments(const double coeff, const double* newmanVectorElementIt, double* adjusmentIt)
    {
        for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
        {
            *adjusmentIt += coeff * (*newmanVectorElementIt);

            ++adjusmentIt;
            ++newmanVectorElementIt;
        }
    }

    static void computeNewmanFacesAdjusments(const size_t* baseFaceIndexesIt,
                                             const size_t* neigbourFaceIndexesIt,
                                             const size_t nFaces,
#ifdef VALUE_CONDITION
                                             const double newmanValue,
#endif
                                             const Coordinates nodes[],
                                             const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                             const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
#ifndef VALUE_CONDITION
                                             void* memoryBuffer,
#endif
                                             double baseLinearTetrahedronsAdjusments[][Basis::N_FUNCTIONS]
#ifdef SHARED_BOUNDARY
                                             ,
                                             double neighbourLinearTetrahedronsAdjusments[][Basis::N_FUNCTIONS]
#endif
    )
    {
        const double (*massVectors)[Basis::N_FUNCTIONS] = FaceLAC<Basis>::getMassVectors();
        size_t elementIndex;
        uint8_t localIndex;

#ifndef VALUE_CONDITION
        double* newmanValues = (double*)memoryBuffer;
        double* newmanVector = newmanValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;
        double transpJacobian[LocalCoordinates2D::COUNT * Coordinates::COUNT];
#endif

        size_t baseFaceNodesTags[constants::triangle::N_NODES];

#ifdef SHARED_BOUNDARY
        size_t neighbourFaceNodesTags[constants::triangle::N_NODES];
#endif

        Coordinates normal;

        for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
        {
            elementIndex = *baseFaceIndexesIt >> 2;
            localIndex = *baseFaceIndexesIt & 3;

            NodesFunctions::extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], baseFaceNodesTags);

#ifdef VALUE_CONDITION
            CoordinatesFunctions::computeNormal(nodes, baseFaceNodesTags, normal);
            double det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
            addFaceAdjusments<Basis::N_FUNCTIONS>(newmanValue * det, massVectors[localIndex], baseLinearTetrahedronsAdjusments[elementIndex]);
#else
            CoordinatesFunctions::computeTranspJacobianToFace(nodes, baseFaceNodesTags, transpJacobian);
            Boundary::FunctionCondition::compute<Boundary::FunctionCondition::Type::NEWMAN, NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps>(
                NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                transpJacobian,
                nodes[*(baseFaceNodesTags)],
                newmanValues);

            double det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

            FaceLAC<Basis>::computePowerVector(localIndex, newmanValues, newmanVector);
            addFaceAdjusments<Basis::N_FUNCTIONS>(det, newmanVector, baseLinearTetrahedronsAdjusments[elementIndex]);
#endif
#ifdef SHARED_BOUNDARY
            elementIndex = *neigbourFaceIndexesIt >> 2;
            localIndex = *neigbourFaceIndexesIt & 3;

#ifdef VALUE_CONDITION
            addFaceAdjusments<Basis::N_FUNCTIONS>(newmanValue * det, massVectors[localIndex], neighbourLinearTetrahedronsAdjusments[elementIndex]);
#else
            NodesFunctions::extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], baseFaceNodesTags);

            uint8_t changingIndexes = NodesFunctions::computeChangingIndexes(baseFaceNodesTags, neighbourFaceNodesTags);
            if (changingIndexes != 0)
            {
                NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::changeValuesOrder(newmanValues, Basis::N_FUNCTIONS, changingIndexes);
            }

#endif

            ++neigbourFaceIndexesIt;
#endif
            ++baseFaceIndexesIt;

        }
    }
}