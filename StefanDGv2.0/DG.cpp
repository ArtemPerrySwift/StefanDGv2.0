#include "DG.h"
#include "LinearLagrangeBasis.h"
#include "MaterialPhase.h"
#include "Boundary.h"
#include "NonconformInterface.h"
#include "ElementLAC.h"
#include "CoordinatesFunctions.h"
#include "FaceDGCalculator.h"
#include "ArrayFunctions.h"
#include <Eigen/IterativeLinearSolvers>
#include <gmsh.h>

namespace DG
{
	using Basis = LinearLagrangeBasis;
    namespace StefanTask
    {
        template<MaterialPhase::State>
        inline static double initialCondition(const Coordinates& point)
        {
            return 0.0;
        }

        template<>
        inline static double initialCondition<MaterialPhase::SOLID>(const Coordinates& point)
        {
            return 0.0;
        }

        template<>
        inline static double initialCondition<MaterialPhase::LIQUID>(const Coordinates& point)
        {
            return 0.0;
        }

        template<MaterialPhase::State state>
        static void computeInitialCondition(const LocalCoordinates3D *localPointIt,
                                            const double transpJacobian[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                            const Coordinates tetrahedronStartPoint,
                                            double* initlialConditionValueIt)
        {
            Coordinates point;
            for (uint8_t i = 0; i < NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps; ++i)
            {
                point = tetrahedronStartPoint;
                point.x += transpJacobian[0] * localPointIt->u;
                point.y += transpJacobian[1] * localPointIt->u;
                point.z += transpJacobian[2] * localPointIt->u;

                point.x += transpJacobian[3] * localPointIt->v;
                point.y += transpJacobian[4] * localPointIt->v;
                point.z += transpJacobian[5] * localPointIt->v;

                point.x += transpJacobian[6] * localPointIt->w;
                point.y += transpJacobian[7] * localPointIt->w;
                point.z += transpJacobian[8] * localPointIt->w;

                *initlialConditionValueIt = initialCondition<state>(point);
                ++initlialConditionValueIt;
                ++localPointIt;
            }
        }

        template<Boundary::FunctionCondition::Type conditionType>
        inline static double computeBoundaryCondition(const Coordinates& point)
        {
            return 0.0;
        }

        template<>
        inline static double computeBoundaryCondition<Boundary::FunctionCondition::Type::DIRICHLET>(const Coordinates& point)
        {
            return 0.0;
        }

        template<>
        inline static double computeBoundaryCondition<Boundary::FunctionCondition::Type::NEWMAN>(const Coordinates& point)
        {
            return 0.0;
        }

        template<Boundary::FunctionCondition::Type conditionType>
        static void computeBoundaryCodnition(const LocalCoordinates2D* localPointIt,
                                              const double transpJacobian[LocalCoordinates2D::COUNT * Coordinates::COUNT],
                                              const Coordinates triangleStartPoint,
                                              double* initlialConditionValueIt)
        {
            Coordinates point;
            for (uint8_t i = 0; i < NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps; ++i)
            {
                point = triangleStartPoint;
                point.x += transpJacobian[0] * localPointIt->u;
                point.y += transpJacobian[1] * localPointIt->u;
                point.z += transpJacobian[2] * localPointIt->u;

                point.x += transpJacobian[3] * localPointIt->v;
                point.y += transpJacobian[4] * localPointIt->v;
                point.z += transpJacobian[5] * localPointIt->v;

                *initlialConditionValueIt = computeBoundaryCondition<conditionType>(point);
                ++initlialConditionValueIt;
                ++localPointIt;
            }
        }

        template<size_t N_FUNCTIONS>
        static void computeTetrahedronAdjusments(const double* massMatrixElementIt,
                                                 const double gamma,
                                                 const double* stiffnessMatrixElementIt,
                                                 const double lambda,
                                                 const double* trackPowerVectorElementIt,
                                                 const double determinant,
                                                 double* bilinearAdjusmentIt,
                                                 double* linearAdjusmentIt)
        {
            double multiplier = determinant * gamma;
            for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
            {
                *linearAdjusmentIt = multiplier * (*trackPowerVectorElementIt);
                for (uint8_t j = 0; j < N_FUNCTIONS; ++j)
                {
                    *bilinearAdjusmentIt = determinant * (*massMatrixElementIt - *stiffnessMatrixElementIt);

                    ++massMatrixElementIt;
                    ++stiffnessMatrixElementIt;
                    ++bilinearAdjusmentIt;
                }

                ++trackPowerVectorElementIt;
                ++linearAdjusmentIt;
            }
        }


        template<MaterialPhase::State MATERIAL_STATE>
        void processTetrahedrons(const Coordinates nodes[],
                                 const size_t (*tetrahedronNodesIndexesIt)[constants::tetrahedron::N_NODES],
                                 const size_t nTetrahedrons,
                                 const double lambda,
                                 const double gamma,
                                 void* calculationMemoryBuffer,
                                 double(*localJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                 double (*initialXIt)[Basis::N_FUNCTIONS],
                                 double (*bilinearTetrahedronsAdjusmentsIt)[ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                 double (*linearTetrahedronsAdjusmentsIt)[Basis::N_FUNCTIONS])
        {
            const double* templateMassMatrix = ElementLAC<Basis>::getMassMatrix();
            double* stiffnessMatrix = (double*)calculationMemoryBuffer;
            double* localVector = stiffnessMatrix + ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;

            double* initlialConditionValues = localVector + Basis::N_FUNCTIONS;

            const LocalCoordinates3D* localGradients = ElementLAC<Basis>::getIntegrationLocalGradients();
            Coordinates* gradients = (Coordinates*)(initlialConditionValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps);

            const LLt* massSLAESolverPtr = ElementLAC<Basis>::getMassSLAESolverPtr();
            double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT];

            for (size_t i = 0; i < nTetrahedrons; ++i)
            {
                double det = CoordinatesFunctions::computeTranspJacobianMatrix(nodes, *tetrahedronNodesIndexesIt, transpJacobianMatrix);
                CoordinatesFunctions::computeLocalJacobianMatrix(transpJacobianMatrix, det, *localJacobianMatrixIt);

                CoordinatesFunctions::translate(*localJacobianMatrixIt, localGradients, Basis::N_FUNCTIONS * NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps, gradients);
                ElementLAC<Basis>::computeStiffnessMatrix(gradients, stiffnessMatrix);

                computeInitialCondition<MATERIAL_STATE>(NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::localPoints, transpJacobianMatrix, nodes[**tetrahedronNodesIndexesIt], initlialConditionValues);

                ElementLAC<Basis>::computePowerVector(initlialConditionValues, localVector);
                massSLAESolverPtr->solve(localVector, *initialXIt);

                computeTetrahedronAdjusments<Basis::N_FUNCTIONS>(templateMassMatrix,
                    gamma,
                    stiffnessMatrix,
                    lambda,
                    localVector,
                    det,
                    *bilinearTetrahedronsAdjusmentsIt,
                    *linearTetrahedronsAdjusmentsIt);

                ++initialXIt;
                ++localJacobianMatrixIt;

                ++bilinearTetrahedronsAdjusmentsIt;
                ++linearTetrahedronsAdjusmentsIt;

                ++tetrahedronNodesIndexesIt;

            }
        }

        template<size_t N_FUNCTIONS>
        static void addInteriorFaceAdjusmets(const double* flowMatrixElementIt, const double* massMatrixElementIt, const double penaltyCoeff, const double detLambdaCoeff, double* adjusmentIt)
        {
            for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
            {
                *adjusmentIt = detLambdaCoeff * (*flowMatrixElementIt - penaltyCoeff * (*massMatrixElementIt));
                const double* ijFlowMatrixElementPtr = flowMatrixElementIt;
                const double* jiFlowMatrixElementPtr = flowMatrixElementIt;
                for (uint8_t j = i + 1; j < N_FUNCTIONS; ++j)
                {
                    ++massMatrixElementIt;
                    ++ijFlowMatrixElementPtr;
                    jiFlowMatrixElementPtr += N_FUNCTIONS;

                    *adjusmentIt += detLambdaCoeff * ((*ijFlowMatrixElementPtr + *jiFlowMatrixElementPtr) / 2.0 - penaltyCoeff * (*massMatrixElementIt));
                }

                flowMatrixElementIt += N_FUNCTIONS + 1;
            }
        }

        template<size_t N_FUNCTIONS>
        void computeInteriorFaceCrossAdjusmets(const double* crossFlowMatrix1ElementIt,
            const double* crossFlowMatrix2ElementIt,
            const double* massMatrixElementIt,
            const double penaltyCoeff,
            const double detLambdaCoeff,
            double* adjusmentIt)
        {
            for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
            {
                *adjusmentIt = detLambdaCoeff * ((*crossFlowMatrix1ElementIt + *crossFlowMatrix2ElementIt) / 2.0 - penaltyCoeff * (*massMatrixElementIt));
                const double* jiFlowMatrix2ElementPtr = crossFlowMatrix2ElementIt;
                for (uint8_t j = 0; j < N_FUNCTIONS; ++j)
                {
                    ++massMatrixElementIt;
                    ++crossFlowMatrix1ElementIt;
                    jiFlowMatrix2ElementPtr += N_FUNCTIONS;

                    *adjusmentIt = detLambdaCoeff * ((*crossFlowMatrix1ElementIt + *jiFlowMatrix2ElementPtr) / 2.0 - penaltyCoeff * (*massMatrixElementIt));
                }

                ++crossFlowMatrix2ElementIt;
            }
        }

        template<size_t N_FUNCTIONS>
        void addInteriorFaceCrossAdjusmets(const double* crossFlowMatrix1ElementIt,
            const double* crossFlowMatrix2ElementIt,
            const double* massMatrixElementIt,
            const double penaltyCoeff,
            const double detLambdaCoeff,
            double* adjusmentIt)
        {
            for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
            {
                *adjusmentIt += detLambdaCoeff * ((*crossFlowMatrix1ElementIt + *crossFlowMatrix2ElementIt) / 2.0 - penaltyCoeff * (*massMatrixElementIt));
                const double* jiFlowMatrix2ElementPtr = crossFlowMatrix2ElementIt;
                for (uint8_t j = 0; j < N_FUNCTIONS; ++j)
                {
                    ++massMatrixElementIt;
                    ++crossFlowMatrix1ElementIt;
                    jiFlowMatrix2ElementPtr += N_FUNCTIONS;

                    *adjusmentIt += detLambdaCoeff * ((*crossFlowMatrix1ElementIt + *jiFlowMatrix2ElementPtr) / 2.0 - penaltyCoeff * (*massMatrixElementIt));
                }

                ++crossFlowMatrix2ElementIt;
            }
        }


        static void computeBilinearMatrixOnTriangle(const double rowsValues[FaceLAC<Basis>::N_BASIS_VALUES],
                                                    const double columnsValues[FaceLAC<Basis>::N_BASIS_VALUES],
                                                    double matrix[Basis::N_FUNCTIONS * Basis::N_FUNCTIONS])
        {
            const double* rawValues = rowsValues;
            for (uint8_t i = 0; i < Basis::N_FUNCTIONS; ++i)
            {
                const double* columnValues = columnsValues;
                for (uint8_t j = 0; j < Basis::N_FUNCTIONS; ++j)
                {
                    *matrix = NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::integrateProduct(rawValues, columnValues);
                    ++matrix;

                    columnValues += NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;
                }
                rawValues += NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;
            }
        }
        
        static void addInteriorFaceAdjusments(const uint8_t localIndex,
                                               const double commonMultiplier,
                                               const double penalty,
                                               const Coordinates& unitNormal,
                                               const double localJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                               Coordinates gradients[FaceLAC<Basis>::N_BASIS_VALUES],
                                               double flowMatrix[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                               double normalDerivatives[FaceLAC<Basis>::N_BASIS_VALUES],
                                               double bilinearAdjusments[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])
        {

            Coordinates* gradients = (Coordinates*)(normalDerivatives + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps);

            double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;

            CoordinatesFunctions::translate(localJacobianMatrix, FaceLAC<Basis>::getLocalGradientsByFace()[localIndex], FaceLAC<Basis>::N_BASIS_VALUES, gradients);

            CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, unitNormal, normalDerivatives);

            computeBilinearMatrixOnTriangle(FaceLAC<Basis>::getValuesByFace()[localIndex], normalDerivatives, flowMatrix);
            addInteriorFaceAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, FaceLAC<Basis>::getMassMatrixies()[localIndex], facePenalty, commonMultiplier, bilinearAdjusments);
        }

        static void computeOutwardNormal(const Coordinates nodes[], const size_t tetrahedronNodesTags[constants::tetrahedron::N_NODES], const uint8_t faceLocalIndex, Coordinates& outwardNormal)
        {
            Coordinates sideVector;
            double transpJacobian[LocalCoordinates2D::COUNT * Coordinates::COUNT];

            switch (faceLocalIndex)
            {
            case 0:
            {
                CoordinatesFunctions::computeTranspJacobianTo0Face(nodes, tetrahedronNodesTags, transpJacobian);
                CoordinatesFunctions::computeDiffrence(nodes[tetrahedronNodesTags[0]], nodes[tetrahedronNodesTags[1]], sideVector);
                break;
            }
            case 1:
            {
                CoordinatesFunctions::computeTranspJacobianTo0Face(nodes, tetrahedronNodesTags, transpJacobian);
                CoordinatesFunctions::computeDiffrence(nodes[tetrahedronNodesTags[1]], nodes[tetrahedronNodesTags[2]], sideVector);
                break;
            }
            case 2:
            {
                CoordinatesFunctions::computeTranspJacobianTo0Face(nodes, tetrahedronNodesTags, transpJacobian);
                CoordinatesFunctions::computeDiffrence(nodes[tetrahedronNodesTags[2]], nodes[tetrahedronNodesTags[3]], sideVector);
                break;
            }
            case 3:
            {
                CoordinatesFunctions::computeTranspJacobianTo0Face(nodes, tetrahedronNodesTags, transpJacobian);
                CoordinatesFunctions::computeDiffrence(nodes[tetrahedronNodesTags[3]], nodes[tetrahedronNodesTags[0]], sideVector);
                break;
            }
            }

            CoordinatesFunctions::computeNormal(transpJacobian, outwardNormal);
            if ((outwardNormal.x * sideVector.x + outwardNormal.y * sideVector.y + outwardNormal.z * sideVector.z) < 0.0)
            {
                outwardNormal.x = -outwardNormal.x;
                outwardNormal.y = -outwardNormal.y;
                outwardNormal.z = -outwardNormal.z;
            }
        }
        
        static void processInteriorFaces(const size_t* faceIndexesSide1It,
                                         const size_t* faceIndexesSide2It,
                                         const size_t nInteriorFaces,
                                         const double lambda,
                                         const double penalty,
                                         const Coordinates nodes[],
                                         const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                         void *calculationBuffer,
                                         double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                         double (*crossAdjusmentsIt)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                         double bilinearTetrahedronsAdjusmentsSide1[][ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                         double bilinearTetrahedronsAdjusmentsSide2[][ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])
        {
            const double* const* massMatrixByFace = FaceLAC<Basis>::getMassMatrixies();
            const double (*crossMassMatrixies)[constants::tetrahedron::N_FACES][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS] = FaceLAC<Basis>::getCrossMassMatrixies();

            const LocalCoordinates3D(*gradientsByFace)[FaceLAC<Basis>::N_BASIS_VALUES] = FaceLAC<Basis>::getLocalGradientsByFace();

            double* flowMatrix = (double*)calculationBuffer;
            double* crossFlowMatrix1 = flowMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* crossFlowMatrix2 = crossFlowMatrix1 + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* normalDerivatives = crossFlowMatrix2 + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;

            Coordinates* gradients = (Coordinates*)(normalDerivatives + FaceLAC<Basis>::N_BASIS_VALUES);
            Coordinates normal;

            for (size_t faceIndex = 0; faceIndex < nInteriorFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                size_t elementIndex = *faceIndexesSide1It >> 2;
                uint8_t localIndex0 = *faceIndexesSide1It && 3;
                ++faceIndexesSide1It;

                computeOutwardNormal(nodes, tetrahedronsNodesTags[elementIndex], localIndex0, normal);
                
                double det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

                normal.x /= det;
                normal.y /= det;
                normal.z /= det;

                double commonMultiplier = det * lambda;

                addInteriorFaceAdjusments(localIndex0,
                                          commonMultiplier,
                                          penalty,
                                          normal,
                                          localJacobianMatrix[elementIndex],
                                          gradients,
                                          flowMatrix,
                                          normalDerivatives,
                                          bilinearTetrahedronsAdjusmentsSide1[elementIndex]);

                elementIndex = *faceIndexesSide2It >> 2;
                uint8_t localIndex1 = *faceIndexesSide2It && 3;
                ++faceIndexesSide2It;

                FaceLAC<Basis>::computeFlowMatrix(localIndex1, normalDerivatives, crossFlowMatrix1);

                normal.x = -normal.x;
                normal.y = -normal.y;
                normal.z = -normal.z;

                addInteriorFaceAdjusments(localIndex1,
                                          commonMultiplier,
                                          penalty,
                                          normal,
                                          localJacobianMatrix[elementIndex],
                                          gradients,
                                          flowMatrix,
                                          normalDerivatives,
                                          bilinearTetrahedronsAdjusmentsSide2[elementIndex]);

                FaceLAC<Basis>::computeFlowMatrix(localIndex0, normalDerivatives, crossFlowMatrix2);

                computeInteriorFaceCrossAdjusmets<Basis::N_FUNCTIONS>(crossFlowMatrix2,
                                                                      crossFlowMatrix1,
                                                                      crossMassMatrixies[localIndex0][localIndex1],
                                                                      facePenalty,
                                                                      commonMultiplier,
                                                                      *crossAdjusmentsIt);

                ++crossAdjusmentsIt;
            }
        }


        static double computeFaceDeterminant(const uint8_t localIndex,
                                             const Coordinates nodes[],
                                             const size_t tetrahedronNodesTags[constants::tetrahedron::N_NODES])
        {
            Coordinates normal;
            double transpJacobian[LocalCoordinates2D::COUNT * Coordinates::COUNT];
            switch (localIndex)
            {
            case 0:
            {
                CoordinatesFunctions::computeTranspJacobianTo0Face(nodes, tetrahedronNodesTags, transpJacobian);
                break;
            }
            case 1:
            {
                CoordinatesFunctions::computeTranspJacobianTo1Face(nodes, tetrahedronNodesTags, transpJacobian);
                break;
            }
            case 2:
            {
                CoordinatesFunctions::computeTranspJacobianTo2Face(nodes, tetrahedronNodesTags, transpJacobian);
                break;
            }
            case 3:
            {
                CoordinatesFunctions::computeTranspJacobianTo3Face(nodes, tetrahedronNodesTags, transpJacobian);
                break;
            }
            }

            CoordinatesFunctions::computeNormal(transpJacobian, normal);
            return sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);;
        }

        static void processFacesOnPlaneInterface(const size_t* faceIndexesSide1It,
                                                 const size_t* faceIndexesSide2It,
                                                 const size_t nInteriorFaces,
                                                 const double lambda,
                                                 const double penalty,
                                                 const Coordinates nodes[],
                                                 const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                 void* calculationBuffer,
                                                 double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                 double (*crossAdjusmentsIt)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                 double bilinearTetrahedronsAdjusmentsSide1[][ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                 double bilinearTetrahedronsAdjusmentsSide2[][ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])
        {
            if (nInteriorFaces == 0)
            {
                return;
            }

            const double* const* massMatrixByFace = FaceLAC<Basis>::getMassMatrixies();
            const double (*crossMassMatrixies)[constants::tetrahedron::N_FACES][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS] = FaceLAC<Basis>::getCrossMassMatrixies();

            const LocalCoordinates3D(*gradientsByFace)[FaceLAC<Basis>::N_BASIS_VALUES] = FaceLAC<Basis>::getLocalGradientsByFace();

            double* flowMatrix = (double*)calculationBuffer;
            double* crossFlowMatrix1 = flowMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* crossFlowMatrix2 = crossFlowMatrix1 + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* normalDerivatives = crossFlowMatrix2 + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;

            Coordinates* gradients = (Coordinates*)(normalDerivatives + FaceLAC<Basis>::N_BASIS_VALUES);

            Coordinates unitNormal;
            Coordinates reverseUnitNormal;

            computeOutwardNormal(nodes, tetrahedronsNodesTags[*faceIndexesSide1It >> 2], *faceIndexesSide1It && 3, unitNormal);
            double det = sqrt(unitNormal.x * unitNormal.x + unitNormal.y * unitNormal.y + unitNormal.z * unitNormal.z);

            unitNormal.x /= det;
            unitNormal.y /= det;
            unitNormal.z /= det;

            reverseUnitNormal.x = -unitNormal.x;
            reverseUnitNormal.y = -unitNormal.y;
            reverseUnitNormal.z = -unitNormal.z;

            for (size_t faceIndex = 0; faceIndex < nInteriorFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                size_t elementIndex = *faceIndexesSide1It >> 2;
                uint8_t localIndex0 = *faceIndexesSide1It && 3;
                ++faceIndexesSide1It;

                det = computeFaceDeterminant(localIndex0, nodes, tetrahedronsNodesTags[elementIndex]);

                double commonMultiplier = det * lambda;

                addInteriorFaceAdjusments(localIndex0,
                    commonMultiplier,
                    penalty,
                    unitNormal,
                    localJacobianMatrix[elementIndex],
                    gradients,
                    flowMatrix,
                    normalDerivatives,
                    bilinearTetrahedronsAdjusmentsSide1[elementIndex]);

                elementIndex = *faceIndexesSide2It >> 2;
                uint8_t localIndex1 = *faceIndexesSide2It && 3;
                ++faceIndexesSide2It;

                FaceLAC<Basis>::computeFlowMatrix(localIndex1, normalDerivatives, crossFlowMatrix1);

                addInteriorFaceAdjusments(localIndex1,
                    commonMultiplier,
                    penalty,
                    reverseUnitNormal,
                    localJacobianMatrix[elementIndex],
                    gradients,
                    flowMatrix,
                    normalDerivatives,
                    bilinearTetrahedronsAdjusmentsSide2[elementIndex]);

                FaceLAC<Basis>::computeFlowMatrix(localIndex0, normalDerivatives, crossFlowMatrix2);

                computeInteriorFaceCrossAdjusmets<Basis::N_FUNCTIONS>(crossFlowMatrix2,
                    crossFlowMatrix1,
                    crossMassMatrixies[localIndex0][localIndex1],
                    facePenalty,
                    commonMultiplier,
                    *crossAdjusmentsIt);

                ++crossAdjusmentsIt;
            }
        }

        void extractFaceNodesTags(const uint8_t faceIndex, const size_t tetrahedronNodesTags[constants::tetrahedron::N_NODES], size_t faceNodes[constants::triangle::N_NODES])
        {
            switch (faceIndex)
            {
            case 0:
            {
                faceNodes[0] = tetrahedronNodesTags[1];
                faceNodes[1] = tetrahedronNodesTags[2];
                faceNodes[2] = tetrahedronNodesTags[3];
                break;
            }
            case 1:
            {
                faceNodes[0] = tetrahedronNodesTags[0];
                faceNodes[1] = tetrahedronNodesTags[2];
                faceNodes[2] = tetrahedronNodesTags[3];
                break;
            }
            case 2:
            {
                faceNodes[0] = tetrahedronNodesTags[0];
                faceNodes[1] = tetrahedronNodesTags[1];
                faceNodes[2] = tetrahedronNodesTags[3];
                break;
            }
            case 3:
            {
                faceNodes[0] = tetrahedronNodesTags[0];
                faceNodes[1] = tetrahedronNodesTags[1];
                faceNodes[2] = tetrahedronNodesTags[2];
                break;
            }
            }
        }

        static void addInteriorFaceFragmentAdjusments(const double commonMultiplier,
                                                      const double penalty,
                                                      const Coordinates& unitNormal,
                                                      const double localJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                      const double basisValues[],
                                                      const LocalCoordinates3D localGradients[],
                                                      Coordinates gradients[FaceLAC<Basis>::N_BASIS_VALUES],
                                                      double flowMatrix[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                      double massMatrix[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                      double normalDerivatives[FaceLAC<Basis>::N_BASIS_VALUES],
                                                      double bilinearAdjusments[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])
        {

            Coordinates* gradients = (Coordinates*)(normalDerivatives + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps);

            double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;

            CoordinatesFunctions::translate(localJacobianMatrix, localGradients, FaceLAC<Basis>::N_BASIS_VALUES, gradients);

            CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, unitNormal, normalDerivatives);

            computeBilinearMatrixOnTriangle(basisValues, basisValues, massMatrix);
            computeBilinearMatrixOnTriangle(basisValues, normalDerivatives, flowMatrix);

            addInteriorFaceAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, massMatrix, facePenalty, commonMultiplier, bilinearAdjusments);
        }


        static void processFacesOnNonconformInterface(const size_t faceIndexesSide1[],
                                                      const size_t faceIndexesSide2[],
                                                      const unsigned int* fragmentFaceSurfaceIndexSide1It,
                                                      const unsigned int* fragmentFaceSurfaceIndexSide2It,
                                                      const size_t nFragments,
                                                      const size_t* fragmentsTrianglesStartIndexIt,
                                                      const Coordinates fragmentsModelsNodes[],
                                                      const size_t (*triangleNodesTagsIt)[constants::triangle::N_NODES],
                                                      const Coordinates nodes[],
                                                      const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                      const double lambda,
                                                      const double penalty,
                                                      void* calculationBuffer,
                                                      double localJacobianMatrixes[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                      double (*crossAdjusmentsIt)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                      double bilinearTetrahedronsAdjusmentsSide1[][ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                      double bilinearTetrahedronsAdjusmentsSide2[][ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])
        {
            if (nFragments == 0)
            {
                return;
            }

            const LocalCoordinates3D(*gradientsByFace)[FaceLAC<Basis>::N_BASIS_VALUES] = FaceLAC<Basis>::getLocalGradientsByFace();

            double* flowMatrix = (double*)calculationBuffer;
            double* massMatrix = flowMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* crossMatrix = massMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;

            double* crossFlowMatrixes[2] = { crossMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS,  crossMatrix + 2 * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS };
            double* sidesBilinearAdjusments[2] = { crossFlowMatrixes[1] + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS, crossFlowMatrixes[1] + 2 * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS };

            LocalCoordinates2D* localPoints = (LocalCoordinates2D*)(sidesBilinearAdjusments[1] + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS);
            Coordinates* gradients = (Coordinates*)(localPoints + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps);
            LocalCoordinates3D* localGradients = (LocalCoordinates3D*)(gradients + FaceLAC<Basis>::N_BASIS_VALUES);

            double* basisValues[2];
            basisValues[0] = (double*)(localGradients + FaceLAC<Basis>::N_BASIS_VALUES);
            basisValues[1] = basisValues[0] + FaceLAC<Basis>::N_BASIS_VALUES;

            double* normalDerivatives[2];
            normalDerivatives[0] = basisValues[1] + FaceLAC<Basis>::N_BASIS_VALUES;
            normalDerivatives[1] = normalDerivatives[1] + FaceLAC<Basis>::N_BASIS_VALUES;

            Coordinates unitNormals[2];
            double translMatrix[2][4];
            size_t faceNodes[2][constants::triangle::N_NODES];
            LocalCoordinates2D fragmentPointsLocalCoordinates[2][constants::triangle::N_NODES];
            size_t elementsIndexes[2] = {SIZE_MAX, SIZE_MAX };
            uint8_t localIndexes[2];
            Coordinates facesBasePoints[2];

            double* sidesBilinearAdjusments[2];
            size_t faceIndex = faceIndexesSide1[*fragmentFaceSurfaceIndexSide1It];
            computeOutwardNormal(nodes, tetrahedronsNodesTags[faceIndex >> 2], faceIndex && 3, unitNormals[0]);
            unitNormals[1].x = -unitNormals[0].x;
            unitNormals[1].y = -unitNormals[0].y;
            unitNormals[1].z = -unitNormals[0].z;

            void (*computeTransitionMatrix2D)(const Coordinates nodes[], const size_t * triangleNodeIndexIt, double transitionMatrix[4]);
            void (*translatePointsCoordinatesToLocal2D)(const Coordinates & basePoint,
                const Coordinates points[],
                const size_t pointsIndexes[constants::triangle::N_NODES],
                const double transitionMatrix[4],
                LocalCoordinates2D localCoordinates[constants::triangle::N_NODES]);


            if (unitNormals[0].x > unitNormals[0].y && unitNormals[0].x > unitNormals[0].z)
            {
                computeTransitionMatrix2D = CoordinatesFunctions::computeTransitionMatrix2DExcludeX;
                translatePointsCoordinatesToLocal2D = CoordinatesFunctions::translatePointsCoordinatesToLocal2DExcludeX;
            }
            else if (unitNormals[0].y > unitNormals[0].z)
            {
                computeTransitionMatrix2D = CoordinatesFunctions::computeTransitionMatrix2DExcludeY;
                translatePointsCoordinatesToLocal2D = CoordinatesFunctions::translatePointsCoordinatesToLocal2DExcludeY;
            }
            else
            {
                computeTransitionMatrix2D = CoordinatesFunctions::computeTransitionMatrix2DExcludeZ;
                translatePointsCoordinatesToLocal2D = CoordinatesFunctions::translatePointsCoordinatesToLocal2DExcludeZ;
            }

            for (size_t iFragment = 0, iTriangle = *fragmentsTrianglesStartIndexIt; iFragment < nFragments; ++iFragment)
            {
                ++fragmentsTrianglesStartIndexIt;
                const size_t triangleEndIndex = *fragmentsTrianglesStartIndexIt;

                for (uint8_t i = 0; i < 2; ++i)
                {
                    const size_t* sideFaceIndexes = faceIndexesSide1;
                    const unsigned int* faceSurfaceIndex = fragmentFaceSurfaceIndexSide1It;
                    if (i == 0)
                    {
                        sideFaceIndexes = faceIndexesSide1;
                        faceSurfaceIndex = fragmentFaceSurfaceIndexSide1It;
                    }
                    else
                    {
                        sideFaceIndexes = faceIndexesSide2;
                        faceSurfaceIndex = fragmentFaceSurfaceIndexSide2It;
                    }
                    faceIndex = sideFaceIndexes[*faceSurfaceIndex];

                    if (elementsIndexes[i] != (faceIndex >> 2))
                    {
                        elementsIndexes[i] = faceIndex >> 2;
                        localIndexes[i] = faceIndex && 3;

                        extractFaceNodesTags(localIndexes[i], tetrahedronsNodesTags[elementsIndexes[i]], faceNodes[i]);
                        computeTransitionMatrix2D(nodes, faceNodes[i], translMatrix[i]);
                        facesBasePoints[i] = nodes[tetrahedronsNodesTags[elementsIndexes[i]][0]];
                    }

                    memset(sidesBilinearAdjusments[i], 0, ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS * sizeof(double));
                }

                memset(crossMatrix, 0, ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS * sizeof(double));


                for (; iTriangle < triangleEndIndex; ++iTriangle)
                {

                    double det = CoordinatesFunctions::computeTriangleDeterminant(fragmentsModelsNodes, *triangleNodesTagsIt);
                    double commonMultiplier = det * lambda;

                    for (uint8_t i = 0; i < 2; ++i)
                    {
                        translatePointsCoordinatesToLocal2D(facesBasePoints[0],
                                                            fragmentsModelsNodes,
                                                            *triangleNodesTagsIt,
                                                            translMatrix[0],
                                                            fragmentPointsLocalCoordinates[0]);
                        LocalCoordinates2D translCoefficients[2];
                        CoordinatesFunctions::computeTranslationCoefficients2D(fragmentPointsLocalCoordinates[i], translCoefficients);
                        CoordinatesFunctions::translateFragmentLocalPoints(fragmentPointsLocalCoordinates[i][0],
                            NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                            NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
                            translCoefficients,
                            localPoints);


                        switch (localIndexes[i])
                        {
                        case 0:
                        {
                            Basis::computeOnFace0(localPoints, NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps, basisValues[i]);
                            Basis::computeOnFace0(localPoints, NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps, localGradients);
                            break;
                        }
                        case 1:
                        {
                            Basis::computeOnFace1(localPoints, NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps, basisValues[i]);
                            Basis::computeOnFace1(localPoints, NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps, localGradients);
                            break;
                        }
                        case 2:
                        {
                            Basis::computeOnFace2(localPoints, NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps, basisValues[i]);
                            Basis::computeOnFace2(localPoints, NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps, localGradients);
                            break;
                        }
                        case 3:
                        {
                            Basis::computeOnFace3(localPoints, NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps, basisValues[i]);
                            Basis::computeOnFace3(localPoints, NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps, localGradients);
                            break;
                        }
                        }


                        CoordinatesFunctions::translate(localJacobianMatrixes[elementsIndexes[i]], localGradients, FaceLAC<Basis>::N_BASIS_VALUES, gradients);
                        CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, unitNormals[i], normalDerivatives[i]);

                        computeBilinearMatrixOnTriangle(basisValues[i], basisValues[i], massMatrix);
                        computeBilinearMatrixOnTriangle(basisValues[i], normalDerivatives[i], flowMatrix);

                        addInteriorFaceAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, massMatrix, penalty, commonMultiplier, sidesBilinearAdjusments[i]);
                    }

                    computeBilinearMatrixOnTriangle(basisValues[1], normalDerivatives[0], crossFlowMatrixes[0]);
                    computeBilinearMatrixOnTriangle(basisValues[0], normalDerivatives[1], crossFlowMatrixes[1]);
                    computeBilinearMatrixOnTriangle(basisValues[0], basisValues[1], massMatrix);

                    addInteriorFaceCrossAdjusmets<Basis::N_FUNCTIONS>(crossFlowMatrixes[1],
                                                                      crossFlowMatrixes[0],
                                                                      massMatrix,
                                                                      penalty,
                                                                      commonMultiplier,
                                                                      crossMatrix);

                    ++triangleNodesTagsIt;
                }

                double* side0BilinearAdjusments = sidesBilinearAdjusments[0];
                double* side1BilinearAdjusments = sidesBilinearAdjusments[1];

                double* side0GlobalBilinearAdjusments = bilinearTetrahedronsAdjusmentsSide1[elementsIndexes[0]];
                double* side1GlobalBilinearAdjusments = bilinearTetrahedronsAdjusmentsSide2[elementsIndexes[1]];
                double* crossGlobalAdjusments = *crossAdjusmentsIt;

                for (uint8_t i = 0; i < FaceLAC<Basis>::N_BASIS_VALUES; ++i)
                {
                    side0GlobalBilinearAdjusments[i] += side0BilinearAdjusments[i];
                    side1GlobalBilinearAdjusments[i] += side1BilinearAdjusments[i];
                    crossGlobalAdjusments[i] += crossMatrix[i];
                }

                ++crossAdjusmentsIt;
                ++fragmentFaceSurfaceIndexSide1It;
                ++fragmentFaceSurfaceIndexSide2It;
            }
        }

        static void processRegions(const unsigned int nRegions,
                                   const MaterialPhase* const* regionMaterialPhaseIt,
                                   const double dt,
                                   const double penalty,
                                   void* calculationBuffer,
                                   const size_t* regionStartTetrahedronIndexIt,
                                   const Coordinates nodes[],
                                   const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                   const size_t* regionInteriorFacesStartIndexIt,
                                   const size_t* baseFacesIndexes,
                                   const size_t* neighborFacesIndexes,
                                   double  localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                   double* bilinearAdjusmentsIt,
                                   double** regionsBilinearAdjasmentsIt,
                                   double* crossBilinearAdjusmentsIt,
                                   double* initialIt,
                                   double* fIt,
                                   double** regionsLinearAdjusments)
        {
            for (unsigned int i = 0; i < nRegions; ++i)
            {
                size_t nRegionTetrahedrons = *(regionStartTetrahedronIndexIt + 1) - *regionStartTetrahedronIndexIt;
                if ((*regionMaterialPhaseIt)->state == MaterialPhase::LIQUID)
                {
                    processTetrahedrons<MaterialPhase::LIQUID>(nodes,
                        tetrahedronsNodesTags + *regionStartTetrahedronIndexIt,
                        nRegionTetrahedrons,
                        (*regionMaterialPhaseIt)->thermalConductivity,
                        (*regionMaterialPhaseIt)->density * (*regionMaterialPhaseIt)->heatCapacity / dt,
                        calculationBuffer,
                        localJacobianMatrix + *regionStartTetrahedronIndexIt,
                        (double(*)[Basis::N_FUNCTIONS])initialIt,
                        (double(*)[ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])bilinearAdjusmentsIt,
                        (double(*)[Basis::N_FUNCTIONS])fIt);
                }
                else
                {
                    processTetrahedrons<MaterialPhase::LIQUID>(nodes,
                        tetrahedronsNodesTags + *regionStartTetrahedronIndexIt,
                        nRegionTetrahedrons,
                        (*regionMaterialPhaseIt)->thermalConductivity,
                        (*regionMaterialPhaseIt)->density * (*regionMaterialPhaseIt)->heatCapacity / dt,
                        calculationBuffer,
                        localJacobianMatrix + *regionStartTetrahedronIndexIt,
                        (double(*)[Basis::N_FUNCTIONS])initialIt,
                        (double(*)[ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])bilinearAdjusmentsIt,
                        (double(*)[Basis::N_FUNCTIONS])fIt);
                }

                *regionsBilinearAdjasmentsIt = bilinearAdjusmentsIt - (*regionStartTetrahedronIndexIt) * ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
                *regionsLinearAdjusments = fIt - (*regionStartTetrahedronIndexIt) * Basis::N_FUNCTIONS;

                size_t nInteriorFaces = *(regionInteriorFacesStartIndexIt + 1) - *regionInteriorFacesStartIndexIt;
                processInteriorFaces(baseFacesIndexes + *regionInteriorFacesStartIndexIt,
                                     neighborFacesIndexes + *regionInteriorFacesStartIndexIt,
                                     nInteriorFaces,
                                     (*regionMaterialPhaseIt)->thermalConductivity,
                                     penalty,
                                     nodes,
                                     tetrahedronsNodesTags,
                                     calculationBuffer,
                                     localJacobianMatrix,
                                     (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])crossBilinearAdjusmentsIt,
                                     (double(*)[ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])bilinearAdjusmentsIt - *regionStartTetrahedronIndexIt, 
                                     (double(*)[ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])bilinearAdjusmentsIt - *regionStartTetrahedronIndexIt);

                ++regionInteriorFacesStartIndexIt;

                initialIt += nRegionTetrahedrons * Basis::N_FUNCTIONS;
                fIt += nRegionTetrahedrons * Basis::N_FUNCTIONS;
                bilinearAdjusmentsIt += nRegionTetrahedrons * ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;

                crossBilinearAdjusmentsIt += nInteriorFaces * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;

                ++regionStartTetrahedronIndexIt;
                ++regionMaterialPhaseIt;
                ++regionsBilinearAdjasmentsIt;
            }
        }

        template<size_t N_FUNCTIONS>
        static void addDirichletFaceBilinearAdjusmets(const double* flowMatrixElementIt,
                                                      const double* massMatrixElementIt,
                                                      const double penaltyCoeff,
                                                      const double detLambdaCoeff,
                                                      double* adjusmentIt)
        {
            for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
            {
                *adjusmentIt = detLambdaCoeff * (*flowMatrixElementIt + *flowMatrixElementIt - penaltyCoeff * (*massMatrixElementIt));
                const double* ijFlowMatrixElementPtr = flowMatrixElementIt;
                const double* jiFlowMatrixElementPtr = flowMatrixElementIt;
                for (uint8_t j = i + 1; j < N_FUNCTIONS; ++j)
                {
                    ++massMatrixElementIt;
                    ++ijFlowMatrixElementPtr;
                    jiFlowMatrixElementPtr += N_FUNCTIONS;

                    *adjusmentIt += detLambdaCoeff * (*ijFlowMatrixElementPtr + *jiFlowMatrixElementPtr - penaltyCoeff * (*massMatrixElementIt));
                }

                flowMatrixElementIt += N_FUNCTIONS + 1;
            }
        }

        template<size_t N_FUNCTIONS>
        static void addDirichletFaceLinearAdjusments(const double* powerVectorIt,
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

        static double computeFaceConditionData(const uint8_t localIndex,
                                               const Coordinates nodes[],
                                               const size_t tetrahedronNodesTags[constants::tetrahedron::N_NODES],
                                               double dirichletValues[],
                                               Coordinates &unitNormal)
        {
            Coordinates sideVector;
            double transpJacobian[LocalCoordinates2D::COUNT * Coordinates::COUNT];
            switch (localIndex)
            {
            case 0:
            {
                CoordinatesFunctions::computeTranspJacobianTo0Face(nodes, tetrahedronNodesTags, transpJacobian);
                CoordinatesFunctions::computeDiffrence(nodes[tetrahedronNodesTags[0]], nodes[tetrahedronNodesTags[1]], sideVector);
                computeBoundaryCodnition<Boundary::FunctionCondition::Type::DIRICHLET>(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                    transpJacobian,
                    nodes[tetrahedronNodesTags[1]],
                    dirichletValues);
                break;
            }
            case 1:
            {
                CoordinatesFunctions::computeTranspJacobianTo1Face(nodes, tetrahedronNodesTags, transpJacobian);
                CoordinatesFunctions::computeDiffrence(nodes[tetrahedronNodesTags[1]], nodes[tetrahedronNodesTags[2]], sideVector);
                computeBoundaryCodnition<Boundary::FunctionCondition::Type::DIRICHLET>(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                    transpJacobian,
                    nodes[tetrahedronNodesTags[0]],
                    dirichletValues);
                break;
            }
            case 2:
            {
                CoordinatesFunctions::computeTranspJacobianTo2Face(nodes, tetrahedronNodesTags, transpJacobian);
                CoordinatesFunctions::computeDiffrence(nodes[tetrahedronNodesTags[2]], nodes[tetrahedronNodesTags[3]], sideVector);
                computeBoundaryCodnition<Boundary::FunctionCondition::Type::DIRICHLET>(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                    transpJacobian,
                    nodes[tetrahedronNodesTags[0]],
                    dirichletValues);
                break;
            }
            case 3:
            {
                CoordinatesFunctions::computeTranspJacobianTo3Face(nodes, tetrahedronNodesTags, transpJacobian);
                CoordinatesFunctions::computeDiffrence(nodes[tetrahedronNodesTags[3]], nodes[tetrahedronNodesTags[0]], sideVector);
                computeBoundaryCodnition<Boundary::FunctionCondition::Type::DIRICHLET>(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                    transpJacobian,
                    nodes[tetrahedronNodesTags[0]],
                    dirichletValues);
                break;
            }
            }

            CoordinatesFunctions::computeNormal(transpJacobian, unitNormal);
            if ((unitNormal.x * sideVector.x + unitNormal.y * sideVector.y + unitNormal.z * sideVector.z) < 0.0)
            {
                unitNormal.x = -unitNormal.x;
                unitNormal.y = -unitNormal.y;
                unitNormal.z = -unitNormal.z;
            }

            double det = sqrt(unitNormal.x * unitNormal.x + unitNormal.y * unitNormal.y + unitNormal.z * unitNormal.z);
            unitNormal.x /= det;
            unitNormal.y /= det;
            unitNormal.z /= det;

            return det;
        }

        template<Boundary::FunctionCondition::Type ConditionType>
        static double computeFaceConditionData(const uint8_t localIndex,
                                               const Coordinates nodes[],
                                               const size_t tetrahedronNodesTags[constants::tetrahedron::N_NODES],
                                               double conditionValues[])
        {
            Coordinates normal;
            double transpJacobian[LocalCoordinates2D::COUNT * Coordinates::COUNT];
            switch (localIndex)
            {
            case 0:
            {
                CoordinatesFunctions::computeTranspJacobianTo0Face(nodes, tetrahedronNodesTags, transpJacobian);
                computeBoundaryCodnition<ConditionType>(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                    transpJacobian,
                    nodes[tetrahedronNodesTags[1]],
                    conditionValues);
                break;
            }
            case 1:
            {
                CoordinatesFunctions::computeTranspJacobianTo1Face(nodes, tetrahedronNodesTags, transpJacobian);
                computeBoundaryCodnition<ConditionType>(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                    transpJacobian,
                    nodes[tetrahedronNodesTags[0]],
                    conditionValues);
                break;
            }
            case 2:
            {
                CoordinatesFunctions::computeTranspJacobianTo2Face(nodes, tetrahedronNodesTags, transpJacobian);
                computeBoundaryCodnition<ConditionType>(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                    transpJacobian,
                    nodes[tetrahedronNodesTags[0]],
                    conditionValues);
                break;
            }
            case 3:
            {
                CoordinatesFunctions::computeTranspJacobianTo3Face(nodes, tetrahedronNodesTags, transpJacobian);
                computeBoundaryCodnition<ConditionType>(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                    transpJacobian,
                    nodes[tetrahedronNodesTags[0]],
                    conditionValues);
                break;
            }
            }
            return sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);;
        }

        static double computeFaceUnitNormal(const uint8_t localIndex, const Coordinates nodes[], const size_t tetrahedronNodesTags[constants::tetrahedron::N_NODES], Coordinates& unitNormal)
        {
            double transpJacobian[LocalCoordinates2D::COUNT * Coordinates::COUNT];
            Coordinates sideVector, normal;
            switch (localIndex)
            {
            case 0:
            {
                CoordinatesFunctions::computeTranspJacobianTo0Face(nodes, tetrahedronNodesTags, transpJacobian);
                CoordinatesFunctions::computeDiffrence(nodes[tetrahedronNodesTags[0]], nodes[tetrahedronNodesTags[1]], sideVector);
                break;
            }
            case 1:
            {
                CoordinatesFunctions::computeTranspJacobianTo1Face(nodes, tetrahedronNodesTags, transpJacobian);
                CoordinatesFunctions::computeDiffrence(nodes[tetrahedronNodesTags[1]], nodes[tetrahedronNodesTags[2]], sideVector);
                break;
            }
            case 2:
            {
                CoordinatesFunctions::computeTranspJacobianTo2Face(nodes, tetrahedronNodesTags, transpJacobian);
                CoordinatesFunctions::computeDiffrence(nodes[tetrahedronNodesTags[2]], nodes[tetrahedronNodesTags[3]], sideVector);
                break;
            }
            case 3:
            {
                CoordinatesFunctions::computeTranspJacobianTo3Face(nodes, tetrahedronNodesTags, transpJacobian);
                CoordinatesFunctions::computeDiffrence(nodes[tetrahedronNodesTags[3]], nodes[tetrahedronNodesTags[0]], sideVector);
                break;
            }
            }

            CoordinatesFunctions::computeNormal(transpJacobian, normal);
            if ((normal.x * sideVector.x + normal.y * sideVector.y + normal.z * sideVector.z) < 0.0)
            {
                normal.x = -normal.x;
                normal.y = -normal.y;
                normal.z = -normal.z;
            }

            double det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

            unitNormal.x = normal.x / det;
            unitNormal.y = normal.y / det;
            unitNormal.z = normal.z / det;

            return det;
        }


        static void addDirichletFaceAdjusments(const uint8_t localIndex,
                                               const double lambda,
                                               const double penalty,
                                               const double det,
                                               const Coordinates& unitNormal,
                                               const double dirichletValues[NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps],
                                               const double localJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                               void* memoryBuffer,
                                               double bilinearAdjusments[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                               double linearAdjusments[Basis::N_FUNCTIONS])
        {
            double* flowMatrix = (double*)memoryBuffer;

            double* flowVector = flowMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* powerVector = flowVector + Basis::N_FUNCTIONS;
            double* normalDerivatives = powerVector + Basis::N_FUNCTIONS;

            Coordinates* gradients = (Coordinates*)(normalDerivatives + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps);

            double commonMultiplier = det * lambda;
            double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;

            CoordinatesFunctions::translate(localJacobianMatrix, FaceLAC<Basis>::getLocalGradientsByFace()[localIndex], NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps, gradients);

            CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, unitNormal, normalDerivatives);

            FaceLAC<Basis>::computeFlowMatrix(localIndex, normalDerivatives, flowMatrix);
            addDirichletFaceBilinearAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, FaceLAC<Basis>::getMassMatrixies()[localIndex], facePenalty, commonMultiplier, bilinearAdjusments);

            FaceLAC<Basis>::computeFlowVector(normalDerivatives, dirichletValues, flowVector);
            FaceLAC<Basis>::computePowerVector(localIndex, dirichletValues, powerVector);

            addDirichletFaceLinearAdjusments<Basis::N_FUNCTIONS>(powerVector, flowVector, facePenalty, commonMultiplier, linearAdjusments);
        }

        static void addDirichletFaceAdjusments(const uint8_t localIndex,
                                               const double lambda,
                                               const double penalty,
                                               const double det,
                                               const Coordinates& unitNormal,
                                               const double conditionValue,
                                               const double localJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                               void* memoryBuffer,
                                               double bilinearAdjusments[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                               double linearAdjusments[Basis::N_FUNCTIONS])
        {
            double* flowMatrix = (double*)memoryBuffer;

            double* flowVector = flowMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* normalDerivatives = flowVector + Basis::N_FUNCTIONS;

            Coordinates* gradients = (Coordinates*)(normalDerivatives + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps);

            double commonMultiplier = det * lambda;
            double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;

            CoordinatesFunctions::translate(localJacobianMatrix, FaceLAC<Basis>::getLocalGradientsByFace()[localIndex], NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps, gradients);

            CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, unitNormal, normalDerivatives);

            FaceLAC<Basis>::computeFlowMatrix(localIndex, normalDerivatives, flowMatrix);
            addDirichletFaceBilinearAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, FaceLAC<Basis>::getMassMatrixies()[localIndex], facePenalty, commonMultiplier, bilinearAdjusments);

            FaceLAC<Basis>::computeFlowVector(normalDerivatives, flowVector);

            addDirichletFaceLinearAdjusments<Basis::N_FUNCTIONS>(FaceLAC<Basis>::getMassVectors()[localIndex], flowVector, facePenalty, commonMultiplier * conditionValue, linearAdjusments);
        }

        static void computeDirichletFacesAdjusments(const size_t* faceIndexIt,
                                                    const size_t nFaces,
                                                    const double lambda,
                                                    const double penalty,
                                                    const Coordinates nodes[],
                                                    const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                    const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                    void* memoryBuffer,
                                                    double bilinearTetrahedronsAdjusments[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                    double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
        {

            double* dirichletValues = (double*)memoryBuffer;

            memoryBuffer = (void*)(dirichletValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps);

            Coordinates unitNormal;

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                size_t elementIndex = *faceIndexIt >> 2;
                uint8_t localIndex = *faceIndexIt && 3;

                double det = computeFaceConditionData(localIndex, nodes, tetrahedronsNodesTags[elementIndex], dirichletValues, unitNormal);

                addDirichletFaceAdjusments(localIndex,
                                           lambda,
                                           penalty,
                                           det,
                                           unitNormal,
                                           dirichletValues,
                                           localJacobianMatrix[elementIndex],
                                           memoryBuffer,
                                           bilinearTetrahedronsAdjusments[elementIndex],
                                           linearTetrahedronsAdjusments[elementIndex]);

                ++faceIndexIt;
            }
        }

        static void computeDirichletFacesAdjusments(const size_t* faceIndexIt,
                                                    const size_t nFaces,
                                                    const double lambda,
                                                    const double penalty,
                                                    const double dirichletValue,
                                                    const Coordinates nodes[],
                                                    const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                    const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                    void* memoryBuffer,
                                                    double bilinearTetrahedronsAdjusments[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                    double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
        {
            Coordinates unitNormal;

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                size_t elementIndex = *faceIndexIt >> 2;
                uint8_t localIndex = *faceIndexIt && 3;

                double det = computeFaceUnitNormal(localIndex, nodes, tetrahedronsNodesTags[elementIndex],unitNormal);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValue,
                    localJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments[elementIndex],
                    linearTetrahedronsAdjusments[elementIndex]);

                ++faceIndexIt;
            }
        }

        static void computeDirichletFacesAdjusments(const size_t* baseFaceIndexesIt,
                                                    const size_t* neigbourFaceIndexesIt,
                                                    const size_t nFaces,
                                                    const double lambda,
                                                    const double penalty,
                                                    const Coordinates nodes[],
                                                    const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                    const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                    void* memoryBuffer,
                                                    double bilinearTetrahedronsAdjusments1Area[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                    double bilinearTetrahedronsAdjusments2Area[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                    double linearTetrahedronsAdjusments1Area[][Basis::N_FUNCTIONS],
                                                    double linearTetrahedronsAdjusments2Area[][Basis::N_FUNCTIONS])
        {

            double* dirichletValues = (double*)memoryBuffer;

            memoryBuffer = (void*)(dirichletValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps);

            Coordinates unitNormal;
            size_t elementIndex;
            uint8_t localIndex;

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                //double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                elementIndex = *baseFaceIndexesIt >> 2;
                localIndex = *baseFaceIndexesIt && 3;

                double det = computeFaceConditionData(localIndex, nodes, tetrahedronsNodesTags[elementIndex], dirichletValues, unitNormal);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValues,
                    localJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments1Area[elementIndex],
                    linearTetrahedronsAdjusments1Area[elementIndex]);

                elementIndex = *neigbourFaceIndexesIt >> 2;
                localIndex = *neigbourFaceIndexesIt && 3;

                unitNormal.x = -unitNormal.x;
                unitNormal.y = -unitNormal.y;
                unitNormal.z = -unitNormal.z;

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValues,
                    localJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments1Area[elementIndex],
                    linearTetrahedronsAdjusments1Area[elementIndex]);

                ++baseFaceIndexesIt;
                ++neigbourFaceIndexesIt;
            }
        }

        static void computeDirichletFacesAdjusments(const size_t* baseFaceIndexesIt,
                                                    const size_t* neigbourFaceIndexesIt,
                                                    const size_t nFaces,
                                                    const double lambda,
                                                    const double penalty,
                                                    const double dirichletValue,
                                                    const Coordinates nodes[],
                                                    const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                    const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                    void* memoryBuffer,
                                                    double bilinearTetrahedronsAdjusments1Area[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                    double bilinearTetrahedronsAdjusments2Area[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                    double linearTetrahedronsAdjusments1Area[][Basis::N_FUNCTIONS],
                                                    double linearTetrahedronsAdjusments2Area[][Basis::N_FUNCTIONS])
        {
            Coordinates unitNormal;

            size_t elementIndex;
            uint8_t localIndex;


            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                //double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                elementIndex = *baseFaceIndexesIt >> 2;
                localIndex = *baseFaceIndexesIt && 3;

                double det = computeFaceUnitNormal(localIndex, nodes, tetrahedronsNodesTags[elementIndex], unitNormal);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValue,
                    localJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments1Area[elementIndex],
                    linearTetrahedronsAdjusments1Area[elementIndex]);

                elementIndex = *neigbourFaceIndexesIt >> 2;
                localIndex = *neigbourFaceIndexesIt && 3;

                unitNormal.x = -unitNormal.x;
                unitNormal.y = -unitNormal.y;
                unitNormal.z = -unitNormal.z;

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValue,
                    localJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments1Area[elementIndex],
                    linearTetrahedronsAdjusments1Area[elementIndex]);

                ++baseFaceIndexesIt;
                ++neigbourFaceIndexesIt;
            }
        }

        static void computeDirichletFacesAdjusmentsOnPlane(const size_t *faceIndexIt,
                                                           const size_t nFaces,
                                                           const double lambda,
                                                           const double penalty,
                                                           const Coordinates nodes[],
                                                           const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                           const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                           void* memoryBuffer,
                                                           double bilinearTetrahedronsAdjusments[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                           double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
        {
            if (nFaces == 0)
            {
                return;
            }

            double* dirichletValues = (double*)memoryBuffer;

            memoryBuffer = (void*)(dirichletValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps);

            size_t elementIndex = *faceIndexIt >> 2;
            uint8_t localIndex = *faceIndexIt && 3;

            Coordinates unitNormal;
            computeFaceUnitNormal(localIndex, nodes, tetrahedronsNodesTags[elementIndex], unitNormal);

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                elementIndex = *faceIndexIt >> 2;
                localIndex = *faceIndexIt && 3;

                double det = computeFaceConditionData<Boundary::FunctionCondition::Type::DIRICHLET>(localIndex, nodes, tetrahedronsNodesTags[elementIndex], dirichletValues);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValues,
                    localJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments[elementIndex],
                    linearTetrahedronsAdjusments[elementIndex]);

                ++faceIndexIt;
            }
        }

        static void computeDirichletFacesAdjusmentsOnPlane(const size_t* faceIndexIt,
                                                           const size_t nFaces,
                                                           const double lambda,
                                                           const double penalty,
                                                           const double dirichletValue,
                                                           const Coordinates nodes[],
                                                           const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                           const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                           void* memoryBuffer,
                                                           double bilinearTetrahedronsAdjusments[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                           double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
        {
            if (nFaces == 0)
            {
                return;
            }

            size_t elementIndex = *faceIndexIt >> 2;
            uint8_t localIndex = *faceIndexIt && 3;

            Coordinates unitNormal;
            computeFaceUnitNormal(localIndex, nodes, tetrahedronsNodesTags[elementIndex], unitNormal);

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                elementIndex = *faceIndexIt >> 2;
                localIndex = *faceIndexIt && 3;

                double det = computeFaceDeterminant(localIndex, nodes, tetrahedronsNodesTags[elementIndex]);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValue,
                    localJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments[elementIndex],
                    linearTetrahedronsAdjusments[elementIndex]);

                ++faceIndexIt;
            }
        }

        static void computeDirichletFacesAdjusmentsOnPlane(const size_t* baseFaceIndexesIt,
                                                           const size_t* neigbourFaceIndexesIt,
                                                           const size_t nFaces,
                                                           const double lambda,
                                                           const double penalty,
                                                           const Coordinates nodes[],
                                                           const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                           const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                           void* memoryBuffer,
                                                           double bilinearTetrahedronsAdjusments1Area[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                           double bilinearTetrahedronsAdjusments2Area[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                           double linearTetrahedronsAdjusments1Area[][Basis::N_FUNCTIONS],
                                                           double linearTetrahedronsAdjusments2Area[][Basis::N_FUNCTIONS])
        {
            if (nFaces == 0)
            {
                return;
            }

            double* dirichletValues = (double*)memoryBuffer;

            memoryBuffer = (void*)(dirichletValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps);

            size_t elementIndex = *baseFaceIndexesIt >> 2;
            uint8_t localIndex = *baseFaceIndexesIt && 3;

            Coordinates unitNormal;
            computeFaceUnitNormal(localIndex, nodes, tetrahedronsNodesTags[elementIndex], unitNormal);

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                elementIndex = *baseFaceIndexesIt >> 2;
                localIndex = *baseFaceIndexesIt && 3;

                double det = computeFaceConditionData<Boundary::FunctionCondition::Type::DIRICHLET>(localIndex, nodes, tetrahedronsNodesTags[elementIndex], dirichletValues);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValues,
                    localJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments1Area[elementIndex],
                    linearTetrahedronsAdjusments1Area[elementIndex]);

                elementIndex = *neigbourFaceIndexesIt >> 2;
                localIndex = *neigbourFaceIndexesIt && 3;

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValues,
                    localJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments2Area[elementIndex],
                    linearTetrahedronsAdjusments2Area[elementIndex]);

                ++baseFaceIndexesIt;
                ++neigbourFaceIndexesIt;
            }
        }

        static void computeDirichletFacesAdjusmentsOnPlane(const size_t* baseFaceIndexesIt,
                                                           const size_t* neigbourFaceIndexesIt,
                                                           const size_t nFaces,
                                                           const double lambda,
                                                           const double penalty,
                                                           const double dirichletValue,
                                                           const Coordinates nodes[],
                                                           const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                           const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                           void* memoryBuffer,
                                                           double bilinearTetrahedronsAdjusments1Area[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                           double bilinearTetrahedronsAdjusments2Area[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                           double linearTetrahedronsAdjusments1Area[][Basis::N_FUNCTIONS],
                                                           double linearTetrahedronsAdjusments2Area[][Basis::N_FUNCTIONS])
        {
            if (nFaces == 0)
            {
                return;
            }

            size_t elementIndex = *baseFaceIndexesIt >> 2;
            uint8_t localIndex = *baseFaceIndexesIt && 3;

            Coordinates unitNormal;
            computeFaceUnitNormal(localIndex, nodes, tetrahedronsNodesTags[elementIndex], unitNormal);

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                elementIndex = *baseFaceIndexesIt >> 2;
                localIndex = *baseFaceIndexesIt && 3;

                double det = computeFaceDeterminant(localIndex, nodes, tetrahedronsNodesTags[elementIndex]);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValue,
                    localJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments1Area[elementIndex],
                    linearTetrahedronsAdjusments1Area[elementIndex]);

                elementIndex = *neigbourFaceIndexesIt >> 2;
                localIndex = *neigbourFaceIndexesIt && 3;

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValue,
                    localJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments2Area[elementIndex],
                    linearTetrahedronsAdjusments2Area[elementIndex]);

                ++baseFaceIndexesIt;
                ++neigbourFaceIndexesIt;
            }
        }

        template<size_t N_FUNCTIONS>
        static void addNemanFaceAdjusments(const double determinant, const double* newmanVectorElementIt, double* adjusmentIt)
        {
            for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
            {
                *adjusmentIt += determinant * (*newmanVectorElementIt);

                ++adjusmentIt;
                ++newmanVectorElementIt;
            }
        }

        static void computeNewmanFacesAdjusments(const size_t* baseFaceIndexesIt,
                                                 const size_t* neigbourFaceIndexesIt,
                                                 const size_t nFaces,
                                                 const Coordinates nodes[],
                                                 const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                 const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                 void* memoryBuffer,
                                                 double linearTetrahedronsAdjusments1Area[][Basis::N_FUNCTIONS],
                                                 double linearTetrahedronsAdjusments2Area[][Basis::N_FUNCTIONS])
        {
            double* newmanValues = (double*)memoryBuffer;
            double* newmanVector = newmanValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;

            size_t elementIndex;
            uint8_t localIndex;

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                elementIndex = *baseFaceIndexesIt >> 2;
                localIndex = *baseFaceIndexesIt && 3;
                double det = computeFaceConditionData<Boundary::FunctionCondition::Type::DIRICHLET>(localIndex, nodes, tetrahedronsNodesTags[elementIndex], newmanValues);

                FaceLAC<Basis>::computePowerVector(localIndex, newmanValues, newmanVector);
                addNemanFaceAdjusments<Basis::N_FUNCTIONS>(det, newmanVector, linearTetrahedronsAdjusments1Area[elementIndex]);

                elementIndex = *neigbourFaceIndexesIt >> 2;
                localIndex = *neigbourFaceIndexesIt && 3;

                FaceLAC<Basis>::computePowerVector(localIndex, newmanValues, newmanVector);
                addNemanFaceAdjusments<Basis::N_FUNCTIONS>(det, newmanVector, linearTetrahedronsAdjusments2Area[elementIndex]);

                ++baseFaceIndexesIt;
                ++neigbourFaceIndexesIt;
            }
        }

        static void computeNewmanFacesAdjusments(const size_t* faceIndexIt,
                                                 const size_t nFaces,
                                                 const Coordinates nodes[],
                                                 const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                 const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                 void* memoryBuffer,
                                                 double linearTetrahedronsAdjusmentsArea[][Basis::N_FUNCTIONS])
        {
            double* newmanValues = (double*)memoryBuffer;
            double* newmanVector = newmanValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;

            size_t elementIndex;
            uint8_t localIndex;

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                elementIndex = *faceIndexIt >> 2;
                localIndex = *faceIndexIt && 3;
                double det = computeFaceConditionData<Boundary::FunctionCondition::Type::DIRICHLET>(localIndex, nodes, tetrahedronsNodesTags[elementIndex], newmanValues);

                FaceLAC<Basis>::computePowerVector(localIndex, newmanValues, newmanVector);
                addNemanFaceAdjusments<Basis::N_FUNCTIONS>(det, newmanVector, linearTetrahedronsAdjusmentsArea[elementIndex]);

                ++faceIndexIt;
            }
        }

        template<size_t N_FUNCTIONS>
        static void addNemanFaceAdjusments(const double newmanValue, const double determinant, const double* massVectorElementIt, double* adjusmentIt)
        {
            for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
            {
                *adjusmentIt += newmanValue * determinant * (*massVectorElementIt);

                ++adjusmentIt;
                ++massVectorElementIt;
            }
        }

        static void computeNewmanFacesAdjusments(const size_t* baseFaceIndexesIt,
                                                 const size_t* neigbourFaceIndexesIt,
                                                 const size_t nFaces,
                                                 const double nerwmanValue,
                                                 const Coordinates nodes[],
                                                 const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                 const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                 double linearTetrahedronsAdjusments1Area[][Basis::N_FUNCTIONS],
                                                 double linearTetrahedronsAdjusments2Area[][Basis::N_FUNCTIONS])
        {
            const double (*massVectors)[Basis::N_FUNCTIONS] = FaceLAC<Basis>::getMassVectors();
            size_t elementIndex;
            uint8_t localIndex;

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                elementIndex = *baseFaceIndexesIt >> 2;
                localIndex = *baseFaceIndexesIt && 3;
                double det = computeFaceDeterminant(localIndex, nodes, tetrahedronsNodesTags[elementIndex]);;

                addNemanFaceAdjusments<Basis::N_FUNCTIONS>(nerwmanValue, det, massVectors[localIndex], linearTetrahedronsAdjusments1Area[elementIndex]);

                elementIndex = *neigbourFaceIndexesIt >> 2;
                localIndex = *neigbourFaceIndexesIt && 3;

                addNemanFaceAdjusments<Basis::N_FUNCTIONS>(nerwmanValue, det, massVectors[localIndex], linearTetrahedronsAdjusments2Area[elementIndex]);

                ++baseFaceIndexesIt;
                ++neigbourFaceIndexesIt;
            }
        }

        static void computeNewmanFacesAdjusments(const size_t* faceIndexIt,
                                                 const size_t nFaces,
                                                 const double nerwmanValue,
                                                 const Coordinates nodes[],
                                                 const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                 const double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                 double linearTetrahedronsAdjusmentsArea[][Basis::N_FUNCTIONS])
        {
            const double (*massVectors)[Basis::N_FUNCTIONS] = FaceLAC<Basis>::getMassVectors();
            size_t elementIndex;
            uint8_t localIndex;

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                elementIndex = *faceIndexIt >> 2;
                localIndex = *faceIndexIt && 3;
                double det = computeFaceDeterminant(localIndex, nodes, tetrahedronsNodesTags[elementIndex]);;

                addNemanFaceAdjusments<Basis::N_FUNCTIONS>(nerwmanValue, det, massVectors[localIndex], linearTetrahedronsAdjusmentsArea[elementIndex]);

                ++faceIndexIt;
            }
        }

        static void processInteriorFacesOnInterface(const size_t* baseFaceIndexesIt,
                                                    const size_t* neigbourFaceIndexesIt,
                                                    const size_t nInteriorFaces,
                                                    const double lambda,
                                                    const double penalty,
                                                    const Coordinates nodes[],
                                                    const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                    void* calculationBuffer,
                                                    double localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                    double (*crossAdjusmentsIt)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                    double bilinearTetrahedronsAdjusments1[][ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                    double bilinearTetrahedronsAdjusments2[][ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])
        {
            const double* const* massMatrixByFace = FaceLAC<Basis>::getMassMatrixies();
            const double (*crossMassMatrixies)[constants::tetrahedron::N_FACES][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS] = FaceLAC<Basis>::getCrossMassMatrixies();

            const LocalCoordinates3D(*gradientsByFace)[FaceLAC<Basis>::N_BASIS_VALUES] = FaceLAC<Basis>::getLocalGradientsByFace();

            double* flowMatrix = (double*)calculationBuffer;
            double* crossFlowMatrix1 = flowMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* crossFlowMatrix2 = crossFlowMatrix1 + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* normalDerivatives = crossFlowMatrix2 + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* basisValuesBuffer = normalDerivatives + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps;

            Coordinates* gradients = (Coordinates*)(basisValuesBuffer + NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps);

            double transpJacobian[LocalCoordinates2D::COUNT * Coordinates::COUNT];
            Coordinates normal, sideVector;

            for (size_t faceIndex = 0; faceIndex < nInteriorFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                size_t elementIndex0 = *baseFaceIndexesIt >> 2;
                uint8_t localIndex0 = *baseFaceIndexesIt && 3;
                ++baseFaceIndexesIt;

                switch (localIndex0)
                {
                case 0:
                {
                    CoordinatesFunctions::computeTranspJacobianTo0Face(nodes, tetrahedronsNodesTags[elementIndex0], transpJacobian);
                    CoordinatesFunctions::computeDiffrence(nodes[tetrahedronsNodesTags[elementIndex0][0]], nodes[tetrahedronsNodesTags[elementIndex0][1]], sideVector);
                    break;
                }
                case 1:
                {
                    CoordinatesFunctions::computeTranspJacobianTo1Face(nodes, tetrahedronsNodesTags[elementIndex0], transpJacobian);
                    CoordinatesFunctions::computeDiffrence(nodes[tetrahedronsNodesTags[elementIndex0][1]], nodes[tetrahedronsNodesTags[elementIndex0][2]], sideVector);
                    break;
                }
                case 2:
                {
                    CoordinatesFunctions::computeTranspJacobianTo2Face(nodes, tetrahedronsNodesTags[elementIndex0], transpJacobian);
                    CoordinatesFunctions::computeDiffrence(nodes[tetrahedronsNodesTags[elementIndex0][2]], nodes[tetrahedronsNodesTags[elementIndex0][3]], sideVector);
                    break;
                }
                case 3:
                {
                    CoordinatesFunctions::computeTranspJacobianTo3Face(nodes, tetrahedronsNodesTags[elementIndex0], transpJacobian);
                    CoordinatesFunctions::computeDiffrence(nodes[tetrahedronsNodesTags[elementIndex0][3]], nodes[tetrahedronsNodesTags[elementIndex0][0]], sideVector);
                    break;
                }
                }

                CoordinatesFunctions::computeNormal(transpJacobian, normal);
                if ((normal.x * sideVector.x + normal.y * sideVector.y + normal.z * sideVector.z) < 0.0)
                {
                    normal.x = -normal.x;
                    normal.y = -normal.y;
                    normal.z = -normal.z;
                }

                double commonMultiplier = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z) * lambda;

                CoordinatesFunctions::translate(localJacobianMatrix[elementIndex0],
                    gradientsByFace[localIndex0],
                    NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
                    gradients);

                CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, normal, normalDerivatives);

                FaceLAC<Basis>::computeFlowMatrix(localIndex0, normalDerivatives, flowMatrix);
                addInteriorFaceAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, massMatrixByFace[localIndex0], penalty, commonMultiplier, bilinearTetrahedronsAdjusments1[elementIndex0]);


                size_t elementIndex1 = *neigbourFaceIndexesIt >> 2;
                uint8_t localIndex1 = *neigbourFaceIndexesIt && 3;
                ++neigbourFaceIndexesIt;

                FaceLAC<Basis>::computeFlowMatrix(localIndex1, normalDerivatives, crossFlowMatrix1);

                normal.x = -normal.x;
                normal.y = -normal.y;
                normal.z = -normal.z;

                CoordinatesFunctions::translate(localJacobianMatrix[elementIndex1],
                    gradientsByFace[localIndex1],
                    NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
                    gradients);

                CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, normal, normalDerivatives);

                FaceLAC<Basis>::computeFlowMatrix(localIndex1, normalDerivatives, flowMatrix);
                addInteriorFaceAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, massMatrixByFace[localIndex1], penalty, commonMultiplier, bilinearTetrahedronsAdjusments2[elementIndex1]);

                FaceLAC<Basis>::computeFlowMatrix(localIndex0, normalDerivatives, crossFlowMatrix2);

                computeInteriorFaceCrossAdjusmets<Basis::N_FUNCTIONS>(crossFlowMatrix2,
                    crossFlowMatrix1,
                    crossMassMatrixies[localIndex0][localIndex1],
                    facePenalty,
                    commonMultiplier,
                    *crossAdjusmentsIt);

                ++crossAdjusmentsIt;
            }
        }


        static void proceedBoundaries(const Boundary* boundaryIt,
                                      const unsigned int nBoundaries,
                                      const size_t baseFaceIndexes[],
                                      const size_t neigbourFaceIndexes[],
                                      size_t *surfaceFacesStartIndexeIt,
                                      const MaterialPhase* const regionMaterialPhases[],
                                      const double penalty,
                                      void* calculationBuffer,
                                      const size_t* regionStartTetrahedronIndexIt,
                                      const Coordinates nodes[],
                                      const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                      double  localJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                      double* const regionsBilinearAdjusments[],
                                      double* const regionsLinearAdjusmenst[])
        {

            for (unsigned int i = 0; i < nBoundaries; ++i)
            {
                const size_t nFaces = *(surfaceFacesStartIndexeIt + 1) - *surfaceFacesStartIndexeIt;
                const size_t* boundaryBasesFacesIndexes = baseFaceIndexes + *surfaceFacesStartIndexeIt;
                if (boundaryIt->regionsIndexes[1] == UINT_MAX)
                {

                    switch (boundaryIt->condition->macroType)
                    {
                    case  Boundary::Condition::MacroType::STEFAN_V:
                    case  Boundary::Condition::MacroType::DIRICHLET_V:
                    {
                        if (boundaryIt->type == Boundary::PLANE_SURFACE)
                        {
                            computeDirichletFacesAdjusmentsOnPlane(boundaryBasesFacesIndexes,
                                nFaces,
                                regionMaterialPhases[boundaryIt->regionsIndexes[0]]->thermalConductivity,
                                penalty,
                                ((Boundary::ValueCondition*)(boundaryIt->condition))->value,
                                nodes,
                                tetrahedronsNodesTags,
                                localJacobianMatrix,
                                calculationBuffer,
                                (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[boundaryIt->regionsIndexes[0]],
                                (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]]);
                        }
                        else
                        {
                            computeDirichletFacesAdjusments(boundaryBasesFacesIndexes,
                                nFaces,
                                regionMaterialPhases[boundaryIt->regionsIndexes[0]]->thermalConductivity,
                                penalty,
                                ((Boundary::ValueCondition*)(boundaryIt->condition))->value,
                                nodes,
                                tetrahedronsNodesTags,
                                localJacobianMatrix,
                                calculationBuffer,
                                (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[boundaryIt->regionsIndexes[0]],
                                (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]]);
                        }
                        break;
                    }
                    case  Boundary::Condition::MacroType::DIRICHLET_F:
                    {
                        if (boundaryIt->type == Boundary::PLANE_SURFACE)
                        {
                            computeDirichletFacesAdjusmentsOnPlane(boundaryBasesFacesIndexes,
                                nFaces,
                                regionMaterialPhases[boundaryIt->regionsIndexes[0]]->thermalConductivity,
                                penalty,
                                nodes,
                                tetrahedronsNodesTags,
                                localJacobianMatrix,
                                calculationBuffer,
                                (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[boundaryIt->regionsIndexes[0]],
                                (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]]);
                        }
                        else
                        {
                            computeDirichletFacesAdjusments(boundaryBasesFacesIndexes,
                                nFaces,
                                regionMaterialPhases[boundaryIt->regionsIndexes[0]]->thermalConductivity,
                                penalty,
                                nodes,
                                tetrahedronsNodesTags,
                                localJacobianMatrix,
                                calculationBuffer,
                                (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[boundaryIt->regionsIndexes[0]],
                                (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]]);
                        }
                        break;
                    }
                    case  Boundary::Condition::MacroType::NEWMAN_V:
                    {
                        computeNewmanFacesAdjusments(boundaryBasesFacesIndexes,
                                                     nFaces,
                                                     ((Boundary::ValueCondition*)(boundaryIt->condition))->value,
                                                     nodes,
                                                     tetrahedronsNodesTags,
                                                     localJacobianMatrix,
                                                     (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]]);

                        break;
                    }
                    case  Boundary::Condition::MacroType::NEWMAN_F:
                    {
                        computeNewmanFacesAdjusments(boundaryBasesFacesIndexes,
                            nFaces,
                            nodes,
                            tetrahedronsNodesTags,
                            localJacobianMatrix,
                            calculationBuffer,
                            (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]]);

                        break;
                    }
                    default:
                    {
                        break;
                    }
                    }
                }
                else
                {
                    const size_t* boundaryNeighbourFacesIndexes = neigbourFaceIndexes + *surfaceFacesStartIndexeIt;
                    switch (boundaryIt->condition->macroType)
                    {
                    case  Boundary::Condition::MacroType::STEFAN_V:
                    case  Boundary::Condition::MacroType::DIRICHLET_V:
                    {
                        if (boundaryIt->type == Boundary::PLANE_SURFACE)
                        {
                            computeDirichletFacesAdjusmentsOnPlane(boundaryBasesFacesIndexes,
                                boundaryNeighbourFacesIndexes,
                                nFaces,
                                regionMaterialPhases[boundaryIt->regionsIndexes[0]]->thermalConductivity,
                                penalty,
                                ((Boundary::ValueCondition*)(boundaryIt->condition))->value,
                                nodes,
                                tetrahedronsNodesTags,
                                localJacobianMatrix,
                                calculationBuffer,
                                (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[boundaryIt->regionsIndexes[0]],
                                (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[boundaryIt->regionsIndexes[1]],
                                (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]],
                                (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[1]]);
                        }
                        else
                        {
                            computeDirichletFacesAdjusments(boundaryBasesFacesIndexes,
                                boundaryNeighbourFacesIndexes,
                                nFaces,
                                regionMaterialPhases[boundaryIt->regionsIndexes[0]]->thermalConductivity,
                                penalty,
                                ((Boundary::ValueCondition*)(boundaryIt->condition))->value,
                                nodes,
                                tetrahedronsNodesTags,
                                localJacobianMatrix,
                                calculationBuffer,
                                (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[boundaryIt->regionsIndexes[0]],
                                (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[boundaryIt->regionsIndexes[1]],
                                (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]],
                                (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[1]]);
                        }
                        break;
                    }
                    case  Boundary::Condition::MacroType::DIRICHLET_F:
                    {
                        if (boundaryIt->type == Boundary::PLANE_SURFACE)
                        {
                            computeDirichletFacesAdjusmentsOnPlane(boundaryBasesFacesIndexes,
                                boundaryNeighbourFacesIndexes,
                                nFaces,
                                regionMaterialPhases[boundaryIt->regionsIndexes[0]]->thermalConductivity,
                                penalty,
                                nodes,
                                tetrahedronsNodesTags,
                                localJacobianMatrix,
                                calculationBuffer,
                                (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[boundaryIt->regionsIndexes[0]],
                                (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[boundaryIt->regionsIndexes[1]],
                                (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]],
                                (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[1]]);
                        }
                        else
                        {
                            computeDirichletFacesAdjusments(boundaryBasesFacesIndexes,
                                boundaryNeighbourFacesIndexes,
                                nFaces,
                                regionMaterialPhases[boundaryIt->regionsIndexes[0]]->thermalConductivity,
                                penalty,
                                nodes,
                                tetrahedronsNodesTags,
                                localJacobianMatrix,
                                calculationBuffer,
                                (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[boundaryIt->regionsIndexes[0]],
                                (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[boundaryIt->regionsIndexes[1]],
                                (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]],
                                (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[1]]);
                        }
                        break;
                    }
                    case  Boundary::Condition::MacroType::NEWMAN_V:
                    {
                        computeNewmanFacesAdjusments(boundaryBasesFacesIndexes,
                            boundaryNeighbourFacesIndexes,
                            nFaces,
                            ((Boundary::ValueCondition*)(boundaryIt->condition))->value,
                            nodes,
                            tetrahedronsNodesTags,
                            localJacobianMatrix,
                            (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]],
                            (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[1]]);
                        break;
                    }
                    case  Boundary::Condition::MacroType::NEWMAN_F:
                    {
                        computeNewmanFacesAdjusments(boundaryBasesFacesIndexes,
                            boundaryNeighbourFacesIndexes,
                            nFaces,
                            nodes,
                            tetrahedronsNodesTags,
                            localJacobianMatrix,
                            calculationBuffer,
                            (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]],
                            (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[1]]);
                        break;
                    }

                    case  Boundary::Condition::MacroType::CONFORM_INTERFACE:
                    {
                        printf("Conform interafce condition is not implemented yet");
                        break;
                    }
                    default:
                    {
                        break;
                    }
                    }
                }
            }
        }

        static size_t proceedNoncofromInterfaces(const NonconformInterface nonconformInterfaces[],
                                               const unsigned int nNonconformInterfaces,
                                               const size_t baseFaceIndexes[],
                                               const size_t neigbourFaceIndexes[],
                                               size_t surfaceFacesStartIndexes[],
                                               const double penalty,
                                               void* calculationBuffer,
                                               const size_t* regionStartTetrahedronIndexIt,
                                               const Coordinates nodes[],
                                               const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                               double  localJacobianMatrixes[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                               double* const regionsBilinearAdjusments[],
                                               unsigned int** interfaceFragmentsTrianglesIndexesIt,
                                               unsigned int* interfacFragmentsCountIt,
                                               double** crossAdjusmentsIt)
        {
            size_t nAllFragments = 0;
            for (unsigned int i = 0; i < nNonconformInterfaces; ++i)
            {
                const size_t* boundaryBasesFacesIndexes = baseFaceIndexes + surfaceFacesStartIndexes[nonconformInterfaces->sideIndexes[0]];
                const size_t* boundaryNeigbourFacesIndexes = baseFaceIndexes + surfaceFacesStartIndexes[nonconformInterfaces->sideIndexes[1]];

                unsigned int nFragments = gmsh::model::surfaces::createFragmentModel(nonconformInterfaces->sidesTags, "FragmentModel", interfaceFragmentsTrianglesIndexesIt);
                nAllFragments += nFragments;
                *interfacFragmentsCountIt = nFragments;

                size_t nNodes = gmsh::model::mesh::getNodesCount();
                Coordinates* fragmentsModelsNodes = (Coordinates*)malloc((nNodes + 1) * sizeof(Coordinates));
                size_t* trianglesStartIndexes = (size_t*)malloc(nFragments * sizeof(size_t));
                gmsh::model::mesh::getSurfacesTrianglesStartIndexes(trianglesStartIndexes);
                size_t* triangleNodesTagsIt = (size_t*)(trianglesStartIndexes[nFragments] * constants::triangle::N_NODES * sizeof(size_t));
                gmsh::model::mesh::getTriangles(triangleNodesTagsIt);

                *crossAdjusmentsIt = (double*)malloc(nFragments * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS * sizeof(double));
                processFacesOnNonconformInterface(boundaryBasesFacesIndexes,
                    boundaryNeigbourFacesIndexes,
                    *interfaceFragmentsTrianglesIndexesIt,
                    *interfaceFragmentsTrianglesIndexesIt + nFragments,
                    nFragments,
                    trianglesStartIndexes,
                    fragmentsModelsNodes,
                    (const size_t(*)[constants::triangle::N_NODES])triangleNodesTagsIt,
                    nodes,
                    tetrahedronsNodesTags,
                    nonconformInterfaces->thermalConductivity,
                    penalty,
                    calculationBuffer,
                    localJacobianMatrixes,
                    (double(*)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS]) * crossAdjusmentsIt,
                    (double(*)[ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[nonconformInterfaces->regionsIndexes[0]],
                    (double(*)[ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])regionsBilinearAdjusments[nonconformInterfaces->regionsIndexes[1]]);

                gmsh::model::remove();
                free(fragmentsModelsNodes);
                free(trianglesStartIndexes);
                free(triangleNodesTagsIt);

                ++interfaceFragmentsTrianglesIndexesIt;
                ++crossAdjusmentsIt;
                ++interfacFragmentsCountIt;
            }
            
            return nAllFragments;
        }

        Eigen::Triplet<double, Eigen::Index>* computeAdjusmentsTriplets(const double* tetrahedronsAdjusmentsIt,
            const size_t nTetrahedrons,
            const uint8_t nLocalDOFs,
            const Eigen::Index startIndex,
            Eigen::Triplet<double, Eigen::Index>* tripletsIt)
        {
            Eigen::Index iIndex = startIndex;
            for (size_t k = 0; k < nTetrahedrons; ++k)
            {
                for (uint8_t i = 0; i < nLocalDOFs; ++i)
                {
                    tetrahedronsAdjusmentsIt += i;
                    *tripletsIt = { iIndex, iIndex, *tetrahedronsAdjusmentsIt };
                    ++tripletsIt;
                    ++tetrahedronsAdjusmentsIt;

                    Eigen::Index jIndex = iIndex + 1;
                    for (uint8_t j = i + 1; j < nLocalDOFs; ++j)
                    {
                        *tripletsIt = { iIndex, jIndex, *tetrahedronsAdjusmentsIt };
                        ++tripletsIt;
                        *tripletsIt = { jIndex, iIndex, *tetrahedronsAdjusmentsIt };
                        ++tripletsIt;
                        ++tetrahedronsAdjusmentsIt;

                        ++jIndex;
                    }
                    ++iIndex;
                }
            }

            return tripletsIt;
        }

        Eigen::Triplet<double, Eigen::Index>* computeAdjusmentsTriplets(const double* crossAdjusmentsIt,
            const size_t nCrossElements,
            const size_t* baseFaceIndexIt,
            const size_t* neigbourFaceIndexIt,
            const uint8_t nLocalDOFs,
            const Eigen::Index startIndex,
            Eigen::Triplet<double, Eigen::Index>* tripletsIt)
        {
            for (size_t k = 0; k < nCrossElements; ++k)
            {
                Eigen::Index jInitIndex = startIndex + (*neigbourFaceIndexIt >> 2) * nLocalDOFs;

                Eigen::Index iIndex = startIndex + (*baseFaceIndexIt >> 2) * nLocalDOFs;
                for (uint8_t i = 0; i < nLocalDOFs; ++i)
                {
                    Eigen::Index jIndex = jInitIndex;
                    for (uint8_t j = 0; j < nLocalDOFs; ++j)
                    {
                        *tripletsIt = { iIndex, jIndex, *crossAdjusmentsIt };
                        ++tripletsIt;
                        *tripletsIt = { jIndex, iIndex, *crossAdjusmentsIt };
                        ++tripletsIt;
                        ++jIndex;
                        ++crossAdjusmentsIt;
                    }
                    ++iIndex;
                }

                ++baseFaceIndexIt;
                ++neigbourFaceIndexIt;
            }

            return tripletsIt;
        }

        Eigen::Triplet<double, Eigen::Index>* computeAdjusmentsTriplets(const double* crossAdjusmentsIt,
                                                                        const unsigned int nCrossElements,
                                                                        const size_t* baseFaceIndexes,
                                                                        const size_t* neigbourFaceIndexes,
                                                                        const unsigned int* interfaceFragmentTriangleIndexes,
                                                                        const uint8_t nLocalDOFs,
                                                                        const Eigen::Index startIndexes[2],
                                                                        Eigen::Triplet<double, Eigen::Index>* tripletsIt)
        {
            const unsigned int* interfaceFragmentSide1TriangleIndexIt = interfaceFragmentTriangleIndexes;
            const unsigned int* interfaceFragmentSide2TriangleIndexIt = interfaceFragmentTriangleIndexes;

            for (unsigned int k = 0; k < nCrossElements; ++k)
            {
                Eigen::Index jInitIndex = startIndexes[1] + (neigbourFaceIndexes[*interfaceFragmentSide2TriangleIndexIt] >> 2) * nLocalDOFs;

                Eigen::Index iIndex = startIndexes[0] + (baseFaceIndexes[*interfaceFragmentSide1TriangleIndexIt] >> 2) * nLocalDOFs;
                for (uint8_t i = 0; i < nLocalDOFs; ++i)
                {
                    Eigen::Index jIndex = jInitIndex;
                    for (uint8_t j = 0; j < nLocalDOFs; ++j)
                    {
                        *tripletsIt = { iIndex, jIndex, *crossAdjusmentsIt };
                        ++tripletsIt;
                        *tripletsIt = { jIndex, iIndex, *crossAdjusmentsIt };
                        ++tripletsIt;
                        ++jIndex;
                        ++crossAdjusmentsIt;
                    }
                    ++iIndex;
                }

                ++interfaceFragmentSide2TriangleIndexIt;
                ++interfaceFragmentSide1TriangleIndexIt;
            }

            return tripletsIt;
        }

        static void buildSLAE(const unsigned int nRegions,
                              double* bilinearAdjusments,
                              double* interiorCrossAdjusments,
                              const size_t regionsStartTetrahedronsIndexes[],
                              const size_t regionsInteriorFacesStartIndexes[],
                              const size_t surfacesFacesStartIndexes[],
                              const size_t tetBaseFacesIndexes[],
                              const size_t tetNeigbourFacesIndexes[],
                              const NonconformInterface nonconformInterfaces[],
                              const unsigned int nNonconformInterfaces,
                              unsigned int** interfacesFragmentsTrianglesIndexes,
                              unsigned int* interfacFragmentsCount,
                              double** nonconformCrossAdjusments,
                              Eigen::Index regionsDOFsStartIndexes[],
                              Eigen::Triplet<double, Eigen::Index>* tripletIt)
        {
            Eigen::Index volumeDOFsStartIndex = 0;
            for (unsigned int i = 0; i < nRegions; ++i)
            {   
                regionsDOFsStartIndexes[i] = volumeDOFsStartIndex - regionsStartTetrahedronsIndexes[i] * LinearLagrangeBasis::N_FUNCTIONS;
                size_t nTetrahedrons = regionsStartTetrahedronsIndexes[i + 1] - regionsStartTetrahedronsIndexes[i];
                tripletIt = computeAdjusmentsTriplets(bilinearAdjusments,
                    nTetrahedrons,
                    LinearLagrangeBasis::N_FUNCTIONS,
                    volumeDOFsStartIndex,
                    tripletIt);

                size_t nVolumeDOFs = nTetrahedrons * LinearLagrangeBasis::N_FUNCTIONS;

                size_t nInteriorFaces = regionsInteriorFacesStartIndexes[i + 1] - regionsInteriorFacesStartIndexes[i];
                tripletIt = computeAdjusmentsTriplets(interiorCrossAdjusments,
                                                      nInteriorFaces,
                                                      tetBaseFacesIndexes + regionsInteriorFacesStartIndexes[i],
                                                      tetNeigbourFacesIndexes + regionsInteriorFacesStartIndexes[i],
                                                      LinearLagrangeBasis::N_FUNCTIONS,
                                                      regionsDOFsStartIndexes[i],
                                                      tripletIt);
                
                bilinearAdjusments += nTetrahedrons * ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
                interiorCrossAdjusments += nInteriorFaces * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
                volumeDOFsStartIndex += nVolumeDOFs;
            }

            Eigen::Index startIndexes[2];
            for (unsigned int i = 0; i < nNonconformInterfaces; ++i)
            {
                startIndexes[0] = regionsDOFsStartIndexes[nonconformInterfaces[i].regionsIndexes[0]];
                startIndexes[1] = regionsDOFsStartIndexes[nonconformInterfaces[i].regionsIndexes[1]];
                tripletIt = computeAdjusmentsTriplets(nonconformCrossAdjusments[i],
                                                      interfacFragmentsCount[i],
                                                      tetBaseFacesIndexes + surfacesFacesStartIndexes[nonconformInterfaces[i].sideIndexes[0]],
                                                      tetNeigbourFacesIndexes + surfacesFacesStartIndexes[nonconformInterfaces[i].sideIndexes[1]],
                                                      interfacesFragmentsTrianglesIndexes[i],
                                                      LinearLagrangeBasis::N_FUNCTIONS,
                                                      startIndexes,
                                                      tripletIt);
            }

        }

        const int MAX_ITERATIONS = 1'000;
        void solveInitialIteration(const unsigned int nRegions,
                                          const MaterialPhase* const regionsMaterialPhases[],
                                          const Boundary boundaries[],
                                          const unsigned int nBoundaries,
                                          const NonconformInterface nonconformInterfaces[],
                                          const unsigned int nNonconformInterfaces,
                                          const double dt,
                                          const double penalty,
                                          void* calculationBuffer,
                                          void* additionalBuffer,
                                          Solution* solutionIt)
        {
            size_t* regionsStartTetrahedronsIndexes = (size_t*)additionalBuffer;
            size_t* surfacesFacesStartIndexes = regionsStartTetrahedronsIndexes + nRegions;
            size_t* regionsStartTetrahedronsTags = surfacesFacesStartIndexes + nBoundaries;
            size_t* regionsInteriorFacesStartIndexes = regionsStartTetrahedronsTags + nRegions;

            double** regionsBilinearAdjusments = (double**)(regionsInteriorFacesStartIndexes + nRegions);
            double** regionsLinearAdjusments = (double**)(regionsBilinearAdjusments + nRegions);

            unsigned int** interfacesFragmentsTrianglesIndexes = (unsigned int**)(regionsBilinearAdjusments + nRegions);
            unsigned int* interfacFragmentsCount = (unsigned int*)(interfacesFragmentsTrianglesIndexes + nNonconformInterfaces);
            double** nonconformCrossAdjusments = (double**)(interfacFragmentsCount + nNonconformInterfaces);


            size_t nNodes = gmsh::model::mesh::getNodesCount();
            Coordinates* nodes = (Coordinates*)malloc((nNodes + 1) * sizeof(Coordinates));
            gmsh::model::mesh::getNodes((double(*)[3])nodes);

            gmsh::model::mesh::getRegionsTetrahedronsStartIndexes(regionsStartTetrahedronsIndexes);
            size_t nTetrahedrons = regionsStartTetrahedronsIndexes[nRegions];
            double* bilinearAdjusments = (double*)malloc(nTetrahedrons * ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS * sizeof(double));
            size_t nTetFaces = nTetrahedrons << 2;
            size_t nTetNodes = nTetrahedrons << 2;
            size_t* tetFacesIndexes = (size_t*)malloc(nTetFaces * sizeof(size_t));
            size_t* tetBaseFacesIndexes = tetFacesIndexes;
            size_t* tetNeigbourFacesIndexes = tetFacesIndexes;
            size_t* tetNodesTag = (size_t*)malloc(nTetNodes * sizeof(size_t));

            gmsh::model::mesh::getTetrahedrons(surfacesFacesStartIndexes, regionsStartTetrahedronsTags, tetNodesTag, regionsInteriorFacesStartIndexes, tetBaseFacesIndexes, tetNeigbourFacesIndexes);

            size_t nUniqueFaces = regionsInteriorFacesStartIndexes[nRegions];
            size_t (*facesTetrahedronsIndexes)[2] = (size_t(*)[2])malloc((nUniqueFaces << 1) * sizeof(size_t));
            memset(facesTetrahedronsIndexes, ~0, (nUniqueFaces << 1) * sizeof(size_t));

            size_t nInteriorFaces = regionsInteriorFacesStartIndexes[nRegions] - regionsInteriorFacesStartIndexes[0];

            size_t nBilinearAdjusments = nTetrahedrons * ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS + nInteriorFaces * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;

            size_t nDOFs = nTetrahedrons * Basis::N_FUNCTIONS;

            Eigen::VectorXd initial(nDOFs);
            Eigen::VectorXd f(nDOFs);

            double* initialPtr = initial.data();
            double* fPtr = f.data();

            double (*localJacobianMatrixies)[LocalCoordinates3D::COUNT * Coordinates::COUNT] = 
                (double(*)[LocalCoordinates3D::COUNT * Coordinates::COUNT])malloc(nTetrahedrons * LocalCoordinates3D::COUNT * Coordinates::COUNT * sizeof(double));

            double* bilinearAdjusmentsIt = bilinearAdjusments;
            double* initialIt = initial.data();
            double* fIt = f.data();


            const MaterialPhase* const* materialPhasePtrIt = regionsMaterialPhases;
            const size_t *regionStartTetrahedronIndexIt = regionsStartTetrahedronsIndexes;
            double (*localJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT] = localJacobianMatrixies;

            double* crossBilinearAdjusments = (double*)malloc(nInteriorFaces * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS * sizeof(double));
            processRegions(nRegions,
                           regionsMaterialPhases,
                           dt,
                           penalty,
                           calculationBuffer,
                           regionStartTetrahedronIndexIt,
                           nodes,
                           (size_t(*)[constants::tetrahedron::N_NODES])tetNodesTag,
                           regionsInteriorFacesStartIndexes,
                           tetBaseFacesIndexes,
                           tetNeigbourFacesIndexes,
                           localJacobianMatrixIt,
                           bilinearAdjusmentsIt,
                           regionsBilinearAdjusments,
                           crossBilinearAdjusments,
                           initialIt,
                           fIt,
                           regionsLinearAdjusments);

            proceedBoundaries(boundaries,
                             nBoundaries,
                             tetBaseFacesIndexes,
                             tetNeigbourFacesIndexes,
                             surfacesFacesStartIndexes,
                             regionsMaterialPhases,
                             penalty,
                             calculationBuffer,
                             regionStartTetrahedronIndexIt,
                             nodes,
                             (size_t(*)[constants::tetrahedron::N_NODES])tetNodesTag,
                             localJacobianMatrixIt,
                             regionsBilinearAdjusments,
                             regionsLinearAdjusments);

            size_t nFragments = proceedNoncofromInterfaces(nonconformInterfaces,
                                       nNonconformInterfaces,
                                       tetBaseFacesIndexes,
                                       tetNeigbourFacesIndexes,
                                       surfacesFacesStartIndexes,
                                       penalty,
                                       calculationBuffer,
                                       regionStartTetrahedronIndexIt,
                                       nodes,
                                       (size_t(*)[constants::tetrahedron::N_NODES])tetNodesTag,
                                       localJacobianMatrixIt,
                                       regionsBilinearAdjusments,
                                       interfacesFragmentsTrianglesIndexes,
                                       interfacFragmentsCount,
                                       nonconformCrossAdjusments);

            nBilinearAdjusments += nFragments * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;

            Eigen::Index* regionsDOFsStartIndexes = (Eigen::Index*)malloc(nRegions * sizeof(Eigen::Index));
            Eigen::Triplet<double, Eigen::Index>* triplets = (Eigen::Triplet<double, Eigen::Index>*)malloc(nBilinearAdjusments * sizeof(Eigen::Triplet<double, Eigen::Index>));

            buildSLAE(nRegions,
                     bilinearAdjusments,
                     crossBilinearAdjusments,
                     regionStartTetrahedronIndexIt,
                     regionsInteriorFacesStartIndexes,
                     surfacesFacesStartIndexes,
                     tetBaseFacesIndexes,
                     tetNeigbourFacesIndexes,
                     nonconformInterfaces,
                     nNonconformInterfaces,
                     interfacesFragmentsTrianglesIndexes,
                     interfacFragmentsCount,
                     nonconformCrossAdjusments,
                     regionsDOFsStartIndexes, 
                     triplets);

            Eigen::SparseMatrix<double, 0, Eigen::Index> A(nDOFs, nDOFs);
            A.setFromTriplets(triplets, triplets + nBilinearAdjusments);
            delete[] triplets;
            printf("Building matrix finished\n");

            printf("\nSolving slae...\n");

#if DEBUG_SHOW_SLAE
            std::cout << A << std::endl;
            std::cout << b << std::endl;
#endif
            Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double, 0, Eigen::Index>> solver;
            //Eigen::BiCGSTAB<Eigen::SparseMatrix<double, 0, Eigen::Index>> solver;
            solver.setMaxIterations(MAX_ITERATIONS);
            //Eigen::Index maxIter = solver.maxIterations();
            //double tol = solver.tolerance();
            solver.compute(A);
            solutionIt->_DOFs = solver.solveWithGuess(f, initial);
            solutionIt->_DOFsPtr = solutionIt->_DOFs.data() - gmsh::model::mesh::getFirstTetrahedronTag();

            printf("Solving slae finished\n");
            printf("Estimated error %e; Iterations count %zu\n", solver.error(), solver.iterations());

#if DEBUG_SHOW_SLAE
            std::cout << "Solution" << std::endl;
            std::cout << x << std::endl;
#endif

            for (unsigned int i = 0; i < nNonconformInterfaces; ++i)
            {
                free(interfacesFragmentsTrianglesIndexes[i]);
                free(nonconformCrossAdjusments);
            }

            free(nodes);
            free(tetFacesIndexes);
            free(tetNodesTag);
            free(facesTetrahedronsIndexes);
            free(bilinearAdjusments);
            free(localJacobianMatrixies);
            free(crossBilinearAdjusments);
            free(triplets);
        }


        static void relocateFrontNodes(const Solution& solution,
                                       int frontTag,
                                       int solidRegionTag,
                                       int liquidRegionTag,
                                       MaterialPhase solidMaterialsPhase,
                                       MaterialPhase liquidMaterialsPhase,
                                       const double dt,
                                       const double latentHeat,
                                       Coordinates* frontNodesIt,
                                       Coordinates* normals,
                                       const unsigned int nFrontNodes)
        {
            LocalCoordinates3D localGradient;
            Coordinates gradient;

            double commonMultiplier = dt / (latentHeat * (liquidMaterialsPhase.density + solidMaterialsPhase.density) / 2.0);
            LocalCoordinates3D tetrahedronMainLocalPoints[4] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0} };
            uint8_t iLocalNode;

            uint8_t nElementsTagsMax = 10;
            size_t* elementsTags = (size_t*)malloc(nElementsTagsMax * sizeof(size_t));

            size_t tetrahedronNodeTags[constants::tetrahedron::N_NODES];
            Coordinates tetrahedronNodes[constants::tetrahedron::N_NODES];
            size_t nodeTag;

            uint8_t nElements;
            int entityTag;

            LocalCoordinates3D localGradient;
            double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT];
            double localJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT];

            Coordinates gradientSolidSideSum, gradientLiquidSideSum;
            uint8_t nElementSolidSide, nElementLiquidSide;

            gmsh::model::getNormals(frontTag, (double*)frontNodesIt, nFrontNodes, (double*)normals);
            
            for (size_t i = 0; i < nFrontNodes; ++i)
            {
                gmsh::model::mesh::getNodeByCoordinates(frontNodesIt->x, frontNodesIt->y, frontNodesIt->z, &nodeTag);

                if (gmsh::model::mesh::get3DElementsByCoordinates(frontNodesIt->x, frontNodesIt->y, frontNodesIt->z, nElementsTagsMax, elementsTags, &nElements) == 1)
                {
                    free(elementsTags);
                    nElementsTagsMax = nElements;
                    elementsTags = (size_t*)malloc(nElementsTagsMax * sizeof(size_t));
                    gmsh::model::mesh::get3DElementsByCoordinates(frontNodesIt->x, frontNodesIt->y, frontNodesIt->z, nElementsTagsMax, elementsTags, &nElements);
                }

                gradientSolidSideSum.x = 0;
                gradientSolidSideSum.y = 0;
                gradientSolidSideSum.z = 0;

                gradientLiquidSideSum.x = 0;
                gradientLiquidSideSum.y = 0;
                gradientLiquidSideSum.z = 0;
                

                for (uint8_t j = 0; j < nElements; ++j)
                {
                    gmsh::model::mesh::get3DElement(elementsTags[j], tetrahedronNodeTags, (double*)tetrahedronNodes, &entityTag);
                    for (iLocalNode = 0; tetrahedronNodeTags[iLocalNode] != 0; ++iLocalNode);

                    solution.compute(elementsTags[j], tetrahedronMainLocalPoints[iLocalNode], localGradient);

                    double det = CoordinatesFunctions::computeTranspJacobianMatrix(tetrahedronNodes, transpJacobianMatrix);
                    CoordinatesFunctions::computeLocalJacobianMatrix(transpJacobianMatrix, det, localJacobianMatrix);

                    CoordinatesFunctions::translate(localJacobianMatrix, localGradient, gradient);

                    if (entityTag == solidRegionTag)
                    {
                        gradientSolidSideSum.x += gradient.x;
                        gradientSolidSideSum.y += gradient.y;
                        gradientSolidSideSum.z += gradient.x;

                        ++nElementSolidSide;
                    }
                    else
                    {
                        gradientLiquidSideSum.x += gradient.x;
                        gradientLiquidSideSum.y += gradient.y;
                        gradientLiquidSideSum.z += gradient.x;

                        ++nElementLiquidSide;
                    }
                }

                double flowSolidSide = gradientSolidSideSum.x * normals->x + gradientSolidSideSum.y * normals->y + gradientSolidSideSum.z * normals->z;
                double flowLiquidSide = gradientLiquidSideSum.x * normals->x + gradientLiquidSideSum.y * normals->y + gradientLiquidSideSum.z * normals->z;
                double moveCoeff = (flowSolidSide / nElementSolidSide + flowLiquidSide / nElementLiquidSide) / 2;

                frontNodesIt->x += moveCoeff * normals->x;
                frontNodesIt->y += moveCoeff * normals->y;
                frontNodesIt->z += moveCoeff * normals->z;
            }

            free(elementsTags);
        }
        }

        void solve(const unsigned int nRegions,
                   const MaterialPhase* const regionsMaterialPhases[],
                   const Boundary boundaries[],
                   const unsigned int nBoundaries,
                   const NonconformInterface nonconformInterfaces[],
                   const unsigned int nNonconformInterfaces,
                   const double tMin,
                   const double tMax,
                   const size_t nTSteps,
                   Solution* solutionIt)
        {
            const size_t nTIntervals = nTSteps - 1;

            if (nTIntervals == 0)
            {
                return;
            }

            double dt = (tMax - tMin) / nTIntervals;

            size_t* modelMemoryBuffer = (size_t*)malloc((3 * (nRegions + 1) + nBoundaries + 1) * sizeof(size_t));
            size_t* regionsStartTetrahedronsIndexes = modelMemoryBuffer;
            size_t* surfacesFacesStartIndexes = regionsStartTetrahedronsIndexes + nRegions + 1;
            size_t* regionsStartTetrahedronsTags = surfacesFacesStartIndexes + nBoundaries + 1;
            size_t* regionsInteriorFacesStartIndexes = regionsStartTetrahedronsTags + nRegions;

            double* calculationBaffer = (double*)malloc(1024 * sizeof(double));
            ElementLAC<Basis>::init();
            solveInitialIteration(nRegions,
                                  regionsMaterialPhases,
                                  boundaries,
                                  nBoundaries,
                                  nonconformInterfaces,
                                  nNonconformInterfaces,
                                  dt,
                                  calculationBaffer,
                                  regionsStartTetrahedronsIndexes,
                                  surfacesFacesStartIndexes,
                                  regionsStartTetrahedronsTags,
                                  regionsInteriorFacesStartIndexes,
                                  solutionIt);

            ++solutionIt;

            const uint8_t DEBUG_DOFS_SIZE = 4;

            for (size_t i = 0; i < nTIntervals; ++i)
            {


                ++solutionIt;
            }

            free(modelMemoryBuffer);
            free(calculationBaffer);
        }

    }

	double DG::Solution::compute(const size_t elementTag, const LocalCoordinates3D LocalPoint3D) const
	{
		return LinearLagrangeBasis::compute(LocalPoint3D, _DOFsPtr + LinearLagrangeBasis::N_FUNCTIONS * elementTag);
	}

    void DG::Solution::compute(const size_t elementTag, const LocalCoordinates3D LocalPoint3D, LocalCoordinates3D& gradient) const
    {
        LinearLagrangeBasis::compute(LocalPoint3D, _DOFsPtr + LinearLagrangeBasis::N_FUNCTIONS * elementTag, gradient);
    }
    
}

