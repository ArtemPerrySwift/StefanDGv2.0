#include "DG.h"
#include "Model.h"
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
#include <set>
#include <iostream>

#define DEBUG_SHOW_SLAE false;
#define POWER_FUNCTION false;

template<uint8_t MATRIX_SIZE>
static void printMatrix(const double* matrixElementsIt)
{
    for (uint8_t i = 0; i < MATRIX_SIZE; ++i)
    {
        printf("%lf", *matrixElementsIt);
        ++matrixElementsIt;
        for (uint8_t j = 1; j < MATRIX_SIZE; ++j)
        {
            printf(" %lf", *matrixElementsIt);
            ++matrixElementsIt;
        }
        printf("\n");
    }
}
double T = 0.0;
namespace DG
{
	using Basis = LinearLagrangeBasis;
    namespace StefanTask
    {
        // Фаза 1: Лёд(твёрдая фаза)
        const double K_S     = 2.23;       // Вт / (м * К) - теплопроводность льда
        const double RHO_S   = 916.4;    // кг / м3 - плотность льда
        const double C_S     = 2'046.0;      // Дж / (кг * К) - удельная теплоемкость льда
        const double ALPHA_S = sqrt(K_S / (RHO_S * C_S)); // м / с ^ 0.5 - температуропроводность льда
        const double DBL_ALPHA_S = 2 * ALPHA_S;

        // Фаза 2: Вода(жидкая фаза)
        const double K_L     = 0.584;      // Вт / (м * К) - теплопроводность воды
        const double RHO_L   = 999.5;   // кг / м3 - плотность воды(для простоты Stefan eq берет плотность фазы фронта)
        const double C_L     = 4'200.0;     // Дж / (кг * К) - удельная теплоемкость воды
        const double ALPHA_L = sqrt(K_L / (RHO_L * C_L)); // м / с ^ 0.5
        const double DBL_ALPHA_L = 2 * ALPHA_L;

        // Скрытая теплота плавления и температуры
        const double L       = 334'000.0;     // Дж / кг
        const double T_0     = -2.0;       // град.C(температура на стенке)
        const double T_M     = 0.0;        // град.C(температура плавления)
        const double T_INF   = 20.0;      // град.C(температура воды вдали)
        const double RHO     = RHO_S;      // плотность твердой фазы(выделяющей / поглощающей теплоту)
        const double F_SPEED = 5.7750988819576507e-05;
        const double TIME_ADD = 1 / (4 * F_SPEED * F_SPEED);
        
        const double DIV_S = F_SPEED / ALPHA_S;
        const double DIV_L = F_SPEED / ALPHA_L;

        const double A_S =   (T_M - T_0) / erf(DIV_S);
        const double A_L = (T_M - T_INF) / erfc(DIV_L);

        const double DS = A_S / (ALPHA_S * sqrt(M_PI)) * exp(-(DIV_S * DIV_S)) * 2.0 * F_SPEED;
        const double DL = -A_L / (ALPHA_L * sqrt(M_PI)) * exp(-(DIV_L * DIV_L)) * 2.0 * F_SPEED;

        template<MaterialPhase::State state>
        inline static double analyticalSolution(const Coordinates& point, double t)
        {
            return 0;
        }


        template<>
        inline static double analyticalSolution<MaterialPhase::SOLID>(const Coordinates& point, double t)
        {
            return T_0 + A_S * erf(point.x / (DBL_ALPHA_S * sqrt(TIME_ADD + t)));
        }

        template<>
        inline static double analyticalSolution<MaterialPhase::LIQUID>(const Coordinates& point, double t)
        {
            return T_INF + A_L * erfc(point.x / (DBL_ALPHA_L * sqrt(TIME_ADD + t)));
        }

        template<MaterialPhase::State>
        inline static double initialCondition(const Coordinates& point)
        {
            return 0.0;
        }

        template<>
        inline static double initialCondition<MaterialPhase::SOLID>(const Coordinates& point)
        {
            return T_0 + A_S * erf(DIV_S * point.x);
            //return 271.15 + point.x * 2;
            //return point.x;
            //return 0.0;
        }

        template<>
        inline static double initialCondition<MaterialPhase::LIQUID>(const Coordinates& point)
        {
            return T_INF + A_L * erfc(DIV_L * point.x);
            //return 293.15 - (2 - point.x)*20;
            //return 273.15 + (point.x - 1) * 20;
            //return point.x;
            //return 0.0;
        }

        template<MaterialPhase::State state>
        static void computeInitialCondition(const LocalCoordinates3D* localPointIt,
                                            const double jacobian[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                            const Coordinates tetrahedronStartPoint,
                                            double* initlialConditionValueIt)
        {
            Coordinates point;
            for (uint8_t i = 0; i < NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps; ++i)
            {
                CoordinatesFunctions::translate(*localPointIt, tetrahedronStartPoint, jacobian, point);
                *initlialConditionValueIt = initialCondition<state>(point);
                ++initlialConditionValueIt;
                ++localPointIt;
            }
        }

        template<MaterialPhase::State>
        inline static void initialConditionGradient(const Coordinates& point, Coordinates &gradient)
        {
            return  gradient = { 0, 0, 0 };;
        }

        template<>
        inline static void initialConditionGradient<MaterialPhase::SOLID>(const Coordinates& point, Coordinates& gradient)
        {
            //gradient = { -2, 0, 0 };
            gradient = { DS, 0, 0 };
        }

        template<>
        inline static void initialConditionGradient<MaterialPhase::LIQUID>(const Coordinates& point, Coordinates& gradient)
        {
            //gradient = { -20, 0, 0 };
            gradient = { DL, 0, 0 };
        }

#if POWER_FUNCTION
        inline static double powerFunction(const Coordinates& point)
        {
            //return 293.15 - (point.x + 1)*20;
            return point.x;
        }

        static void computePowerValues(const LocalCoordinates3D* localPointIt,
                                       const double jacobian[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                       const Coordinates tetrahedronStartPoint,
                                       double* powerValueIt)
        {
            Coordinates point;
            for (uint8_t i = 0; i < NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps; ++i)
            {
                CoordinatesFunctions::translate(*localPointIt, tetrahedronStartPoint, jacobian, point);
                *powerValueIt = powerFunction(point);
                ++powerValueIt;
                ++localPointIt;
            }
        }
#endif

        /*
        inline static double analyticalSolution(const Coordinates& point, double t)
        {
            return (point.x);
        }
        */

        template<MaterialPhase::State state>
        static double computeErrors(const LocalCoordinates3D* localPointIt,
                                    const double jacobian[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                    const Coordinates tetrahedronStartPoint,
                                    const double* tetrahedronDOFs,
                                    const double t,
                                    double* errorIt)
        {
            Coordinates point;
            double err, analSolution;
            for (uint8_t i = 0; i < NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps; ++i)
            {
                CoordinatesFunctions::translate(*localPointIt, tetrahedronStartPoint, jacobian, point);
                analSolution = analyticalSolution<state>(point, t);
                err = (analSolution - LinearLagrangeBasis::compute(*localPointIt, tetrahedronDOFs)) / analSolution;
                *errorIt = err * err;
                ++errorIt;
                ++localPointIt;
            }
            return 0.0;
        }

        static double computeEnergyError(const double (*DOFs)[Basis::N_FUNCTIONS],
                                         const unsigned int nRegions,
                                         const MaterialPhase* const* regionMaterialPhaseIt,
                                         const size_t* regionStartTetrahedronIndexIt,
                                         Coordinates nodes[],
                                         const size_t(*tetrahedronNodesIndexesIt)[constants::tetrahedron::N_NODES],
                                         double t)
        {
            double jacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT];
            double errors[NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps];
            double detSum = 0.0;
            double error = 0.0;

            size_t j = *regionStartTetrahedronIndexIt;
            for (unsigned int i = 0; i < nRegions; ++i)
            {
                ++regionStartTetrahedronIndexIt;
                size_t iRegionTetrahedronsEnd = *regionStartTetrahedronIndexIt;

                if ((*regionMaterialPhaseIt)->state == MaterialPhase::SOLID)
                {
                    for (; j < iRegionTetrahedronsEnd; ++j)
                    {
                        double det = CoordinatesFunctions::computeJacobianMatrix(nodes, *tetrahedronNodesIndexesIt, jacobianMatrix);
                        computeErrors<MaterialPhase::SOLID>(NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::localPoints, jacobianMatrix, nodes[**tetrahedronNodesIndexesIt], *DOFs, t, errors);
                        error += det * NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::integrate(errors);
                        detSum += det;

                        ++tetrahedronNodesIndexesIt;
                        ++DOFs;
                    }
                }
                else
                {
                    for (; j < iRegionTetrahedronsEnd; ++j)
                    {
                        double det = CoordinatesFunctions::computeJacobianMatrix(nodes, *tetrahedronNodesIndexesIt, jacobianMatrix);
                        computeErrors<MaterialPhase::LIQUID>(NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::localPoints, jacobianMatrix, nodes[**tetrahedronNodesIndexesIt], *DOFs, t, errors);
                        error += det * NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::integrate(errors);
                        detSum += det;

                        ++tetrahedronNodesIndexesIt;
                        ++DOFs;
                    }
                }
                ++regionMaterialPhaseIt;
            }


            return sqrt(error) / detSum * 6.0;
        }

        template<Boundary::FunctionCondition::Type conditionType>
        inline static double computeBoundaryCondition(const Coordinates& point)
        {
            return 0.0;
        }

        template<>
        inline static double computeBoundaryCondition<Boundary::FunctionCondition::Type::DIRICHLET>(const Coordinates& point)
        {
            //return 0.0;
            return point.x;
            //return (point.x) * T;
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

        double* computeInitialDOFs(const unsigned int nRegions, const MaterialPhase* const regionsMaterialPhases[], size_t* regionsStartTetrahedronsIndexes)
        {
            gmsh::model::mesh::getRegionsTetrahedronsStartIndexes(regionsStartTetrahedronsIndexes);
            size_t nNodes = gmsh::model::mesh::getNodesCount();

            size_t nTetNodes = 4 * regionsStartTetrahedronsIndexes[nRegions];
            double* DOFs = (double*)malloc(nTetNodes * sizeof(double));

            void* memoryBufffer = malloc(nTetNodes * sizeof(size_t) + (nNodes + 1) * sizeof(Coordinates));
            size_t* tetNodesTags = (size_t*)memoryBufffer;
            gmsh::model::mesh::getTetrahedrons(tetNodesTags);

            Coordinates* nodes = (Coordinates*)(tetNodesTags + nTetNodes);
            gmsh::model::mesh::getNodes((double(*)[3])nodes);

            size_t* tetNodesTagIt = tetNodesTags;
            double* DOFIt = DOFs;
            size_t iTet = *regionsStartTetrahedronsIndexes * constants::tetrahedron::N_NODES;
            size_t regionEndTetNodeIndex;
            for (unsigned int i = 0; i < nRegions; ++i)
            {
                ++regionsStartTetrahedronsIndexes;
                regionEndTetNodeIndex = *regionsStartTetrahedronsIndexes * constants::tetrahedron::N_NODES;
                if ((*regionsMaterialPhases)->state == MaterialPhase::SOLID)
                {
                    for (; iTet < regionEndTetNodeIndex; ++iTet)
                    {
                        *DOFIt = initialCondition<MaterialPhase::SOLID>(nodes[*tetNodesTagIt]);
                        ++DOFIt;
                        ++tetNodesTagIt;
                    }

                }
                else
                {
                    for (; iTet < regionEndTetNodeIndex; ++iTet)
                    {
                        *DOFIt = initialCondition<MaterialPhase::LIQUID>(nodes[*tetNodesTagIt]);
                        ++DOFIt;
                        ++tetNodesTagIt;
                    }
                }
                ++regionsMaterialPhases;
            }

            free(memoryBufffer);
            return DOFs;
        }

        const uint8_t ELEMENTS_TAGS_BUFFER_MAX_SIZE = 10;
        void computePreviousSolutionValues(const Solution& previousSolution,
            const LocalCoordinates3D* localPointIt,
            const double jacobian[LocalCoordinates3D::COUNT * Coordinates::COUNT],
            const Coordinates tetrahedronStartPoint,
            double bufferTranspJacobian[LocalCoordinates3D::COUNT * Coordinates::COUNT],
            size_t elementsTagsBuffer[ELEMENTS_TAGS_BUFFER_MAX_SIZE],
            Coordinates tetrahedronNodes[constants::tetrahedron::N_NODES],
            double* initlialConditionValueIt)
        {
            size_t* elementsTags = elementsTagsBuffer;
            uint8_t nElementsMax = ELEMENTS_TAGS_BUFFER_MAX_SIZE;
            uint8_t nElements;

            size_t previousSolutionElementTag = *elementsTags;
            Coordinates point;
            LocalCoordinates3D localCoordinates;

            int entityTag;

            double jacobianMatrixBuffer[LocalCoordinates3D::COUNT * Coordinates::COUNT];

            for (uint8_t i = 0; i < NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps; ++i)
            {
                CoordinatesFunctions::translate(*localPointIt, tetrahedronStartPoint, jacobian, point);

                if (gmsh::model::mesh::get3DElementsByCoordinates(point.x, point.y, point.z, nElementsMax, elementsTags, &nElements))
                {
                    if (elementsTags != elementsTagsBuffer)
                    {
                        free(elementsTags);
                    }
                    elementsTags = (size_t*)malloc(nElements * sizeof(size_t));
                    gmsh::model::mesh::get3DElementsByCoordinates(point.x, point.y, point.z, nElementsMax, elementsTags, &nElements);
                }

                double value = 0;

                if (*elementsTags != previousSolutionElementTag)
                {
                    gmsh::model::mesh::get3DElement(*elementsTags, (double*)tetrahedronNodes, &entityTag);
                    double det = CoordinatesFunctions::computeJacobianMatrix(tetrahedronNodes, jacobianMatrixBuffer);
                    CoordinatesFunctions::computeTranspLocalJacobianMatrix(jacobianMatrixBuffer, det, bufferTranspJacobian);
                    previousSolutionElementTag = *elementsTags;
                }

                Coordinates* basePoint = tetrahedronNodes;
                Coordinates pointsDiff = { point.x - basePoint->x, point.y - basePoint->y, point.z - basePoint->z };
                localCoordinates.u = bufferTranspJacobian[0] * pointsDiff.x + bufferTranspJacobian[3] * pointsDiff.y + bufferTranspJacobian[6] * pointsDiff.z;
                localCoordinates.v = bufferTranspJacobian[1] * pointsDiff.x + bufferTranspJacobian[4] * pointsDiff.y + bufferTranspJacobian[7] * pointsDiff.z;
                localCoordinates.w = bufferTranspJacobian[2] * pointsDiff.x + bufferTranspJacobian[5] * pointsDiff.y + bufferTranspJacobian[8] * pointsDiff.z;

                value = previousSolution.compute(*elementsTags, localCoordinates);

                for (uint8_t j = 1; j < nElements; ++j)
                {
                    gmsh::model::mesh::get3DElement(elementsTags[j], (double*)tetrahedronNodes, &entityTag);
                    double det = CoordinatesFunctions::computeJacobianMatrix(tetrahedronNodes, jacobianMatrixBuffer);
                    CoordinatesFunctions::computeTranspLocalJacobianMatrix(jacobianMatrixBuffer, det, bufferTranspJacobian);

                    Coordinates* basePoint = tetrahedronNodes;
                    Coordinates pointsDiff = { point.x - basePoint->x, point.y - basePoint->y, point.z - basePoint->z };
                    localCoordinates.u = bufferTranspJacobian[0] * pointsDiff.x + bufferTranspJacobian[1] * pointsDiff.y + bufferTranspJacobian[2] * pointsDiff.z;
                    localCoordinates.v = bufferTranspJacobian[3] * pointsDiff.x + bufferTranspJacobian[4] * pointsDiff.y + bufferTranspJacobian[5] * pointsDiff.z;
                    localCoordinates.w = bufferTranspJacobian[6] * pointsDiff.x + bufferTranspJacobian[7] * pointsDiff.y + bufferTranspJacobian[8] * pointsDiff.z;

                    value += previousSolution.compute(elementsTags[j], localCoordinates);
                }

                if (nElements > 1)
                {
                    value /= nElements;
                    previousSolutionElementTag = elementsTags[nElements - 1];
                }

                *initlialConditionValueIt = value;
                ++initlialConditionValueIt;
                ++localPointIt;
            }

            if (elementsTags != elementsTagsBuffer)
            {
                free(elementsTags);
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
                    *bilinearAdjusmentIt = determinant * (gamma*(*massMatrixElementIt) + lambda*(*stiffnessMatrixElementIt));

                    ++massMatrixElementIt;
                    ++stiffnessMatrixElementIt;
                    ++bilinearAdjusmentIt;
                }

                ++trackPowerVectorElementIt;
                ++linearAdjusmentIt;
            }
        }

enum IterationType {INITIAL, GENERAL};

        template<const MaterialPhase::State MATERIAL_STATE>
        void processTetrahedrons(const Coordinates nodes[],
                                 const size_t (*tetrahedronNodesIndexesIt)[constants::tetrahedron::N_NODES],
                                 const size_t nTetrahedrons,
                                 const double lambda,
                                 const double gamma,
                                 void* calculationMemoryBuffer,
                                 double(*transpLocalJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
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
            double jacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT];
#if POWER_FUNCTION
            double powerVector[Basis::N_FUNCTIONS];
#endif
            for (size_t i = 0; i < nTetrahedrons; ++i)
            {
                double det = CoordinatesFunctions::computeJacobianMatrix(nodes, *tetrahedronNodesIndexesIt, jacobianMatrix);
                CoordinatesFunctions::computeTranspLocalJacobianMatrix(jacobianMatrix, det, *transpLocalJacobianMatrixIt);

                CoordinatesFunctions::translate(*transpLocalJacobianMatrixIt, localGradients, Basis::N_FUNCTIONS * NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps, gradients);
                ElementLAC<Basis>::computeStiffnessMatrix(gradients, stiffnessMatrix);

                computeInitialCondition<MATERIAL_STATE>(NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::localPoints, jacobianMatrix, nodes[**tetrahedronNodesIndexesIt], initlialConditionValues);

                ElementLAC<Basis>::computePowerVector(initlialConditionValues, localVector);
                massSLAESolverPtr->solve(localVector, *initialXIt);

#if POWER_FUNCTION
                double* powerValues = initlialConditionValues;
                computePowerValues(NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::localPoints, jacobianMatrix, nodes[**tetrahedronNodesIndexesIt], powerValues);
                ElementLAC<Basis>::computePowerVector(powerValues, powerVector);

                for (uint8_t j = 0; j < Basis::N_FUNCTIONS; ++j)
                {
                    localVector[j] += powerVector[j];
                }
#endif

                computeTetrahedronAdjusments<Basis::N_FUNCTIONS>(templateMassMatrix,
                    gamma,
                    stiffnessMatrix,
                    lambda,
                    localVector,
                    det,
                    *bilinearTetrahedronsAdjusmentsIt,
                    *linearTetrahedronsAdjusmentsIt);

                ++initialXIt;
                ++transpLocalJacobianMatrixIt;

                ++bilinearTetrahedronsAdjusmentsIt;
                ++linearTetrahedronsAdjusmentsIt;

                ++tetrahedronNodesIndexesIt;
            }
        }
        
        void processTetrahedrons(const Coordinates nodes[],
                                 const size_t(*tetrahedronNodesIndexesIt)[constants::tetrahedron::N_NODES],
                                 const size_t nTetrahedrons,
                                 const double lambda,
                                 const double gamma,
                                 const Solution &solution,
                                 void* calculationMemoryBuffer,
                                 double(*transptranspLocalJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
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
            double jacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT];

            double bufferTranspJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT];
            size_t elementTagsBuffer[ELEMENTS_TAGS_BUFFER_MAX_SIZE];
            *elementTagsBuffer = ~size_t(0);
            Coordinates tetrahedronNodesBuffer[constants::tetrahedron::N_NODES];

#if POWER_FUNCTION
            double powerVector[Basis::N_FUNCTIONS];
#endif

            for (size_t i = 0; i < nTetrahedrons; ++i)
            {
                double det = CoordinatesFunctions::computeJacobianMatrix(nodes, *tetrahedronNodesIndexesIt, jacobianMatrix);
                CoordinatesFunctions::computeTranspLocalJacobianMatrix(jacobianMatrix, det, *transptranspLocalJacobianMatrixIt);

                CoordinatesFunctions::translate(*transptranspLocalJacobianMatrixIt, localGradients, Basis::N_FUNCTIONS * NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps, gradients);
                ElementLAC<Basis>::computeStiffnessMatrix(gradients, stiffnessMatrix);

                computePreviousSolutionValues(solution,
                                              NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::localPoints,
                                              jacobianMatrix,
                                              nodes[**tetrahedronNodesIndexesIt],
                                              bufferTranspJacobianMatrix,
                                              elementTagsBuffer,
                                              tetrahedronNodesBuffer,
                                              initlialConditionValues);
               
                ElementLAC<Basis>::computePowerVector(initlialConditionValues, localVector);
                massSLAESolverPtr->solve(localVector, *initialXIt);

#if POWER_FUNCTION
                double* powerValues = initlialConditionValues;
                computePowerValues(NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::localPoints, jacobianMatrix, nodes[**tetrahedronNodesIndexesIt], powerValues);
                ElementLAC<Basis>::computePowerVector(powerValues, powerVector);

                for (uint8_t j = 0; j < Basis::N_FUNCTIONS; ++j)
                {
                    localVector[j] += powerVector[j];
                }

#endif
                computeTetrahedronAdjusments<Basis::N_FUNCTIONS>(templateMassMatrix,
                                                                 gamma,
                                                                 stiffnessMatrix,
                                                                 lambda,
                                                                 localVector,
                                                                 det,
                                                                 *bilinearTetrahedronsAdjusmentsIt,
                                                                 *linearTetrahedronsAdjusmentsIt);

                ++initialXIt;
                ++transptranspLocalJacobianMatrixIt;

                ++bilinearTetrahedronsAdjusmentsIt;
                ++linearTetrahedronsAdjusmentsIt;

                ++tetrahedronNodesIndexesIt;

            }
        }
        

        template<size_t N_FUNCTIONS>
        static void addInteriorFaceAdjusments(const double* flowMatrixElementIt, const double* massMatrixElementIt, const double penaltyCoeff, const double detLambdaCoeff, double* adjusmentIt)
        {
            double* ijAdjusmentPtr;
            double* jiAdjusmentPtr;
            for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
            {
                massMatrixElementIt += i;
                *adjusmentIt += detLambdaCoeff * (-(*flowMatrixElementIt) + penaltyCoeff * (*massMatrixElementIt));

                ijAdjusmentPtr = adjusmentIt;
                jiAdjusmentPtr = adjusmentIt;

                const double* ijFlowMatrixElementPtr = flowMatrixElementIt;
                const double* jiFlowMatrixElementPtr = flowMatrixElementIt;
                for (uint8_t j = i + 1; j < N_FUNCTIONS; ++j)
                {
                    ++massMatrixElementIt;

                    ++ijAdjusmentPtr;
                    jiAdjusmentPtr += N_FUNCTIONS;

                    ++ijFlowMatrixElementPtr;
                    jiFlowMatrixElementPtr += N_FUNCTIONS;

                    double adjuement = detLambdaCoeff * (-(*ijFlowMatrixElementPtr + *jiFlowMatrixElementPtr) / 2.0 + penaltyCoeff * (*massMatrixElementIt));
                    *ijAdjusmentPtr += adjuement;
                    *jiAdjusmentPtr += adjuement;
                }

                ++massMatrixElementIt;
                flowMatrixElementIt += N_FUNCTIONS + 1;
                adjusmentIt += N_FUNCTIONS + 1;
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
            double* adj = adjusmentIt;
            for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
            {
                const double* jiFlowMatrix2ElementPtr = crossFlowMatrix2ElementIt;
                for (uint8_t j = 0; j < N_FUNCTIONS; ++j)
                {
                    *adjusmentIt = detLambdaCoeff * ((*crossFlowMatrix1ElementIt + *jiFlowMatrix2ElementPtr) / 2.0 - penaltyCoeff * (*massMatrixElementIt));

                    ++adjusmentIt;
                    ++massMatrixElementIt;
                    ++crossFlowMatrix1ElementIt;
                    jiFlowMatrix2ElementPtr += N_FUNCTIONS;
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
                const double* jiFlowMatrix2ElementPtr = crossFlowMatrix2ElementIt;
                for (uint8_t j = 0; j < N_FUNCTIONS; ++j)
                {
                    *adjusmentIt += detLambdaCoeff * ((*crossFlowMatrix1ElementIt + *jiFlowMatrix2ElementPtr) / 2.0 - penaltyCoeff * (*massMatrixElementIt));

                    ++massMatrixElementIt;
                    ++crossFlowMatrix1ElementIt;
                    jiFlowMatrix2ElementPtr += N_FUNCTIONS;
                    ++adjusmentIt;
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
                                               const double transpLocalJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                               Coordinates gradients[FaceLAC<Basis>::N_BASIS_VALUES],
                                               double flowMatrix[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                               double normalDerivatives[FaceLAC<Basis>::N_BASIS_VALUES],
                                               double bilinearAdjusments[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])
        {
            double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;

            CoordinatesFunctions::translate(transpLocalJacobianMatrix, FaceLAC<Basis>::getLocalGradientsByFace()[localIndex], FaceLAC<Basis>::N_BASIS_VALUES, gradients);

            CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, unitNormal, normalDerivatives);

            FaceLAC<Basis>::computeFlowMatrix(localIndex, normalDerivatives, flowMatrix);
            //computeBilinearMatrixOnTriangle(FaceLAC<Basis>::getValuesByFace()[localIndex], normalDerivatives, flowMatrix);
            addInteriorFaceAdjusments<Basis::N_FUNCTIONS>(flowMatrix, FaceLAC<Basis>::getMassMatrixies()[localIndex], facePenalty, commonMultiplier, bilinearAdjusments);
        }

        static void extractFaceNodesTags(const uint8_t localIndex,
                                         const size_t tetrahedronNodes[constants::tetrahedron::N_NODES],
                                         size_t faceNodes[constants::triangle::N_NODES])
        {
            switch (localIndex)
            {
            case 0:
            {
                faceNodes[0] = tetrahedronNodes[0];
                faceNodes[1] = tetrahedronNodes[2];
                faceNodes[2] = tetrahedronNodes[1];
                break;
            }
            case 1:
            {
                faceNodes[0] = tetrahedronNodes[0];
                faceNodes[1] = tetrahedronNodes[1];
                faceNodes[2] = tetrahedronNodes[3];
                break;
            }
            case 2:
            {
                faceNodes[0] = tetrahedronNodes[0];
                faceNodes[1] = tetrahedronNodes[3];
                faceNodes[2] = tetrahedronNodes[2];
                break;
            }
            case 3:
            {
                faceNodes[0] = tetrahedronNodes[3];
                faceNodes[1] = tetrahedronNodes[1];
                faceNodes[2] = tetrahedronNodes[2];
                break;
            }
            }
        }

        uint8_t computeChangingIndexes(const size_t original[constants::triangle::N_NODES], const size_t comparison[constants::triangle::N_NODES])
        {
            uint8_t changingIndexes = 0;
            if (*original == *comparison)
            {
                ++original;
                ++comparison;
                if (*original == *comparison)
                {
                    changingIndexes |= 1 << 2;
                    changingIndexes |= 2 << 4;

                    return changingIndexes;
                }

                changingIndexes |= 2 << 2;
                changingIndexes |= 1 << 4;
                return changingIndexes;
            }

            ++comparison;
            if (*original == *comparison)
            {
                changingIndexes |= 1;
                ++original;
                ++comparison;

                if (*original == *comparison)
                {
                    changingIndexes |= 2 << 2;
                    return changingIndexes;
                }

                changingIndexes |= 2 << 4;
                return changingIndexes;

            }

            changingIndexes |= 2;
            ++original;
            if (*original == *comparison)
            {
                changingIndexes |= 1 << 2;
                return changingIndexes;
            }

            changingIndexes |= 1 << 4;
            return changingIndexes;
        }
        
        static void processInteriorFaces(const size_t* faceIndexesSide1It,
                                         const size_t* faceIndexesSide2It,
                                         const size_t nInteriorFaces,
                                         const double lambda,
                                         const double penalty,
                                         const Coordinates nodes[],
                                         const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                         void *calculationBuffer,
                                         double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                         double (*crossAdjusmentsIt)[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                         double bilinearTetrahedronsAdjusmentsSide1[][ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                         double bilinearTetrahedronsAdjusmentsSide2[][ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])
        {
            const double* const* massMatrixByFace = FaceLAC<Basis>::getMassMatrixies();
            const double (*crossMassMatrixies)[constants::tetrahedron::N_FACES][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS] = FaceLAC<Basis>::getCrossMassMatrixies();

            const LocalCoordinates3D(*gradientsByFace)[FaceLAC<Basis>::N_BASIS_VALUES] = FaceLAC<Basis>::getLocalGradientsByFace();

            double* flowMatrix = (double*)calculationBuffer;
            double* crossMassMatrix = flowMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* crossFlowMatrix1 = crossMassMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* crossFlowMatrix2 = crossFlowMatrix1 + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* normalDerivatives = crossFlowMatrix2 + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* basisValues = normalDerivatives + FaceLAC<Basis>::N_BASIS_VALUES;

            Coordinates* gradients = (Coordinates*)(basisValues + FaceLAC<Basis>::N_BASIS_VALUES);
            Coordinates normal;

            size_t facesNodes[2][constants::triangle::N_NODES];

            for (size_t faceIndex = 0; faceIndex < nInteriorFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                size_t elementIndex = *faceIndexesSide1It >> 2;
                uint8_t localIndex0 = *faceIndexesSide1It & 3;
                ++faceIndexesSide1It;

                extractFaceNodesTags(localIndex0, tetrahedronsNodesTags[elementIndex], facesNodes[0]);
                CoordinatesFunctions::computeNormal(nodes, facesNodes[0], normal);
                
                double det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

                normal.x /= det;
                normal.y /= det;
                normal.z /= det;

                double commonMultiplier = det * lambda;

                addInteriorFaceAdjusments(localIndex0,
                                          commonMultiplier,
                                          penalty,
                                          normal,
                                          transpLocalJacobianMatrix[elementIndex],
                                          gradients,
                                          flowMatrix,
                                          normalDerivatives,
                                          bilinearTetrahedronsAdjusmentsSide1[elementIndex]);

                elementIndex = *faceIndexesSide2It >> 2;
                uint8_t localIndex1 = *faceIndexesSide2It & 3;
                ++faceIndexesSide2It;

                extractFaceNodesTags(localIndex1, tetrahedronsNodesTags[elementIndex], facesNodes[1]);
                uint8_t changingIndexes = computeChangingIndexes(facesNodes[0], facesNodes[1]);
                if (changingIndexes != 0)
                {
                    NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::changeValuesOrder(normalDerivatives, Basis::N_FUNCTIONS, changingIndexes);
                }

                FaceLAC<Basis>::computeFlowMatrix(localIndex1, normalDerivatives, crossFlowMatrix1);

                normal.x = -normal.x;
                normal.y = -normal.y;
                normal.z = -normal.z;

                addInteriorFaceAdjusments(localIndex1,
                                          commonMultiplier,
                                          penalty,
                                          normal,
                                          transpLocalJacobianMatrix[elementIndex],
                                          gradients,
                                          flowMatrix,
                                          normalDerivatives,
                                          bilinearTetrahedronsAdjusmentsSide2[elementIndex]);

                changingIndexes = computeChangingIndexes(facesNodes[1], facesNodes[0]);
                if (changingIndexes != 0)
                {
                    NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::changeValuesOrder(normalDerivatives, Basis::N_FUNCTIONS, changingIndexes);
                    NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::changeValuesOrder(FaceLAC<Basis>::getValuesByFace()[localIndex1], Basis::N_FUNCTIONS, changingIndexes, basisValues);
                    computeBilinearMatrixOnTriangle(FaceLAC<Basis>::getValuesByFace()[localIndex0], basisValues, crossMassMatrix);
                }

                FaceLAC<Basis>::computeFlowMatrix(localIndex0, normalDerivatives, crossFlowMatrix2);

                computeInteriorFaceCrossAdjusmets<Basis::N_FUNCTIONS>(crossFlowMatrix2,
                                                                      crossFlowMatrix1,
                                                                      crossMassMatrix,
                                                                      facePenalty,
                                                                      commonMultiplier,
                                                                      *crossAdjusmentsIt);

                ++crossAdjusmentsIt;
            }
        }


        static void processFacesOnPlaneInterface(const size_t* faceIndexesSide1It,
                                                 const size_t* faceIndexesSide2It,
                                                 const size_t nInteriorFaces,
                                                 const double lambda,
                                                 const double penalty,
                                                 const Coordinates nodes[],
                                                 const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                 void* calculationBuffer,
                                                 double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
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

            size_t facesNodes[2][constants::triangle::N_NODES];
            Coordinates normal, unitNormal, reverseUnitNormal;

            extractFaceNodesTags(*faceIndexesSide1It & 3, tetrahedronsNodesTags[*faceIndexesSide1It >> 2], facesNodes[0]);
            CoordinatesFunctions::computeNormal(nodes, facesNodes[0], normal);

            double det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

            unitNormal.x = normal.x / det;
            unitNormal.y = normal.x / det;
            unitNormal.z = normal.x / det;

            reverseUnitNormal.x = -unitNormal.x;
            reverseUnitNormal.y = -unitNormal.y;
            reverseUnitNormal.z = -unitNormal.z;

            for (size_t faceIndex = 0; faceIndex < nInteriorFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                size_t elementIndex = *faceIndexesSide1It >> 2;
                uint8_t localIndex0 = *faceIndexesSide1It & 3;
                ++faceIndexesSide1It;

                extractFaceNodesTags(localIndex0, tetrahedronsNodesTags[elementIndex], facesNodes[0]);
                CoordinatesFunctions::computeNormal(nodes, facesNodes[0], normal);
                det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

                double commonMultiplier = det * lambda;

                addInteriorFaceAdjusments(localIndex0,
                    commonMultiplier,
                    penalty,
                    unitNormal,
                    transpLocalJacobianMatrix[elementIndex],
                    gradients,
                    flowMatrix,
                    normalDerivatives,
                    bilinearTetrahedronsAdjusmentsSide1[elementIndex]);

                elementIndex = *faceIndexesSide2It >> 2;
                uint8_t localIndex1 = *faceIndexesSide2It & 3;
                ++faceIndexesSide2It;

                extractFaceNodesTags(localIndex1, tetrahedronsNodesTags[elementIndex], facesNodes[1]);
                uint8_t changingIndexes = computeChangingIndexes(facesNodes[1], facesNodes[0]);
                if (changingIndexes != 0)
                {
                    NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::changeValuesOrder(normalDerivatives, Basis::N_FUNCTIONS, changingIndexes);
                }


                FaceLAC<Basis>::computeFlowMatrix(localIndex1, normalDerivatives, crossFlowMatrix1);

                addInteriorFaceAdjusments(localIndex1,
                    commonMultiplier,
                    penalty,
                    reverseUnitNormal,
                    transpLocalJacobianMatrix[elementIndex],
                    gradients,
                    flowMatrix,
                    normalDerivatives,
                    bilinearTetrahedronsAdjusmentsSide2[elementIndex]);

                changingIndexes = computeChangingIndexes(facesNodes[0], facesNodes[1]);
                if (changingIndexes != 0)
                {
                    NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::changeValuesOrder(normalDerivatives, Basis::N_FUNCTIONS, changingIndexes);
                }

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

        static void addInteriorFaceFragmentAdjusments(const double commonMultiplier,
                                                      const double penalty,
                                                      const Coordinates& unitNormal,
                                                      const double transpLocalJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                      const double basisValues[],
                                                      const LocalCoordinates3D localGradients[],
                                                      Coordinates gradients[FaceLAC<Basis>::N_BASIS_VALUES],
                                                      double flowMatrix[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                      double massMatrix[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                      double normalDerivatives[FaceLAC<Basis>::N_BASIS_VALUES],
                                                      double bilinearAdjusments[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])
        {
            double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;

            CoordinatesFunctions::translate(transpLocalJacobianMatrix, localGradients, FaceLAC<Basis>::N_BASIS_VALUES, gradients);

            CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, unitNormal, normalDerivatives);

            computeBilinearMatrixOnTriangle(basisValues, basisValues, massMatrix);
            computeBilinearMatrixOnTriangle(basisValues, normalDerivatives, flowMatrix);

            addInteriorFaceAdjusments<Basis::N_FUNCTIONS>(flowMatrix, massMatrix, facePenalty, commonMultiplier, bilinearAdjusments);
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
                                                      double transpLocalJacobianMatrixes[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
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
            normalDerivatives[1] = normalDerivatives[0] + FaceLAC<Basis>::N_BASIS_VALUES;

            Coordinates unitNormals[2];
            double translMatrix[2][4];
            size_t faceNodes[2][constants::triangle::N_NODES];
            LocalCoordinates2D fragmentPointsLocalCoordinates[2][constants::triangle::N_NODES];
            size_t elementsIndexes[2] = {SIZE_MAX, SIZE_MAX};
            uint8_t localIndexes[2];
            Coordinates facesBasePoints[2];

            size_t faceIndex = faceIndexesSide1[*fragmentFaceSurfaceIndexSide1It];
            extractFaceNodesTags(faceIndex & 3, tetrahedronsNodesTags[faceIndex >> 2], faceNodes[0]);
            CoordinatesFunctions::computeNormal(nodes, faceNodes[0], unitNormals[0]);

            double det = (unitNormals->x * unitNormals->x + unitNormals->y * unitNormals->y + unitNormals->z * unitNormals->z);

            unitNormals->x /= sqrt(det);
            unitNormals->y /= sqrt(det);
            unitNormals->z /= sqrt(det);

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

                faceIndex = faceIndexesSide1[*fragmentFaceSurfaceIndexSide1It];
                for (uint8_t i = 0; i < 2; ++i)
                {
                    if (elementsIndexes[i] != (faceIndex >> 2))
                    {
                        elementsIndexes[i] = faceIndex >> 2;
                        localIndexes[i] = faceIndex & 3;

                        extractFaceNodesTags(localIndexes[i], tetrahedronsNodesTags[elementsIndexes[i]], faceNodes[i]);

                        computeTransitionMatrix2D(nodes, faceNodes[i], translMatrix[i]);
                        facesBasePoints[i] = nodes[faceNodes[i][0]];
                    }

                    memset(sidesBilinearAdjusments[i], 0, ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS * sizeof(double));
                    faceIndex = faceIndexesSide2[*fragmentFaceSurfaceIndexSide2It];
                }

                memset(crossMatrix, 0, ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS * sizeof(double));

                for (; iTriangle < triangleEndIndex; ++iTriangle)
                {
                    double det = CoordinatesFunctions::computeTriangleDeterminant(fragmentsModelsNodes, *triangleNodesTagsIt);
                    double commonMultiplier = det * lambda;

                    for (uint8_t i = 0; i < 2; ++i)
                    {
                        translatePointsCoordinatesToLocal2D(facesBasePoints[i],
                                                            fragmentsModelsNodes,
                                                            *triangleNodesTagsIt,
                                                            translMatrix[i],
                                                            fragmentPointsLocalCoordinates[i]);
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


                        CoordinatesFunctions::translate(transpLocalJacobianMatrixes[elementsIndexes[i]], localGradients, FaceLAC<Basis>::N_BASIS_VALUES, gradients);
                        CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, unitNormals[i], normalDerivatives[i]);

                        computeBilinearMatrixOnTriangle(basisValues[i], basisValues[i], massMatrix);
                        computeBilinearMatrixOnTriangle(basisValues[i], normalDerivatives[i], flowMatrix);

                        addInteriorFaceAdjusments<Basis::N_FUNCTIONS>(flowMatrix, massMatrix, penalty, commonMultiplier, sidesBilinearAdjusments[i]);
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

                for (uint8_t i = 0; i < FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS; ++i)
                {
                    side0GlobalBilinearAdjusments[i] += side0BilinearAdjusments[i];
                    side1GlobalBilinearAdjusments[i] += side1BilinearAdjusments[i];
                    crossGlobalAdjusments[i] = crossMatrix[i];
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
                                   const char* currentModelName,
                                   const char* previousModelName,
                                   const Solution &previousSolution,
                                   void* calculationBuffer,
                                   const size_t* regionStartTetrahedronIndexIt,
                                   const Coordinates nodes[],
                                   const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                   const size_t* regionInteriorFacesStartIndexIt,
                                   const size_t* baseFacesIndexes,
                                   const size_t* neighborFacesIndexes,
                                   double  transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                   double* bilinearAdjusmentsIt,
                                   double** regionsBilinearAdjasmentsIt,
                                   double* crossBilinearAdjusmentsIt,
                                   double* initialIt,
                                   double* fIt,
                                   double** regionsLinearAdjusments)
        {
            gmsh::model::setCurrent(previousModelName);
            for (unsigned int i = 0; i < nRegions; ++i)
            {
                size_t nRegionTetrahedrons = *(regionStartTetrahedronIndexIt + 1) - *regionStartTetrahedronIndexIt;
                processTetrahedrons(nodes,
                    tetrahedronsNodesTags + *regionStartTetrahedronIndexIt,
                    nRegionTetrahedrons,
                    (*regionMaterialPhaseIt)->thermalConductivity,
                    (*regionMaterialPhaseIt)->density * (*regionMaterialPhaseIt)->heatCapacity / dt,
                    previousSolution,
                    calculationBuffer,
                    transpLocalJacobianMatrix + *regionStartTetrahedronIndexIt,
                    (double(*)[Basis::N_FUNCTIONS])initialIt,
                    (double(*)[ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])bilinearAdjusmentsIt,
                    (double(*)[Basis::N_FUNCTIONS])fIt);

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
                    transpLocalJacobianMatrix,
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
                ++regionsLinearAdjusments;
            }
            gmsh::model::setCurrent(currentModelName);
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
                                   double  transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
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
                        transpLocalJacobianMatrix + *regionStartTetrahedronIndexIt,
                        (double(*)[Basis::N_FUNCTIONS])initialIt,
                        (double(*)[ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS])bilinearAdjusmentsIt,
                        (double(*)[Basis::N_FUNCTIONS])fIt);
                }
                else
                {
                    processTetrahedrons<MaterialPhase::SOLID>(nodes,
                        tetrahedronsNodesTags + *regionStartTetrahedronIndexIt,
                        nRegionTetrahedrons,
                        (*regionMaterialPhaseIt)->thermalConductivity,
                        (*regionMaterialPhaseIt)->density * (*regionMaterialPhaseIt)->heatCapacity / dt,
                        calculationBuffer,
                        transpLocalJacobianMatrix + *regionStartTetrahedronIndexIt,
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
                                     transpLocalJacobianMatrix,
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
                ++regionsLinearAdjusments;
            }
        }

        template<size_t N_FUNCTIONS>
        static void addDirichletFaceBilinearAdjusmets(const double* flowMatrixElementIt,
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
                *adjusmentIt += detLambdaCoeff * (- *flowMatrixElementIt - *flowMatrixElementIt + penaltyCoeff * (*massMatrixElementIt));
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

                    double adjusment = detLambdaCoeff * (- *ijFlowMatrixElementPtr - *jiFlowMatrixElementPtr + penaltyCoeff * (*massMatrixElementIt));
                    *ijAdjusmentIt += adjusment;
                    *jiAdjusmentIt += adjusment;
                }

                flowMatrixElementIt += N_FUNCTIONS + 1;
                adjusmentIt += N_FUNCTIONS + 1;
                ++massMatrixElementIt;
            }

            adjusmentIt = nullptr;
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


        static double computeFaceConditionData(const Coordinates nodes[],
                                               const size_t triangleNodesTags[constants::triangle::N_NODES],
                                               double dirichletValues[],
                                               Coordinates &unitNormal)
        {
            double transpJacobian[LocalCoordinates2D::COUNT * Coordinates::COUNT];
            CoordinatesFunctions::computeTranspJacobianToFace(nodes, triangleNodesTags, transpJacobian);
            computeBoundaryCodnition<Boundary::FunctionCondition::Type::DIRICHLET>(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                transpJacobian,
                nodes[triangleNodesTags[0]],
                dirichletValues);

            CoordinatesFunctions::computeNormal(transpJacobian, unitNormal);

            double det = sqrt(unitNormal.x * unitNormal.x + unitNormal.y * unitNormal.y + unitNormal.z * unitNormal.z);
            unitNormal.x /= det;
            unitNormal.y /= det;
            unitNormal.z /= det;

            return det;
        }

        template<Boundary::FunctionCondition::Type ConditionType>
        static double computeFaceConditionData(const Coordinates nodes[],
                                               const size_t triangleNodesTags[constants::tetrahedron::N_NODES],
                                               double conditionValues[])
        {
            Coordinates normal;
            double transpJacobian[LocalCoordinates2D::COUNT * Coordinates::COUNT];
            CoordinatesFunctions::computeTranspJacobianToFace(nodes, triangleNodesTags, transpJacobian);
            computeBoundaryCodnition<ConditionType>(NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::localPoints,
                transpJacobian,
                nodes[*triangleNodesTags],
                conditionValues);

            CoordinatesFunctions::computeNormal(transpJacobian, normal);

            return sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);;
        }

        static double computeFaceUnitNormal(const Coordinates nodes[], const size_t triangleNodesTags[constants::tetrahedron::N_NODES], Coordinates& unitNormal)
        {
            Coordinates normal;
            CoordinatesFunctions::computeNormal(nodes, triangleNodesTags, normal);

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
                                               const double transpLocalJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                               void* memoryBuffer,
                                               double bilinearAdjusments[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                               double linearAdjusments[Basis::N_FUNCTIONS])
        {
            double* flowMatrix = (double*)memoryBuffer;

            double* flowVector = flowMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* powerVector = flowVector + Basis::N_FUNCTIONS;
            double* normalDerivatives = powerVector + Basis::N_FUNCTIONS;

            Coordinates* gradients = (Coordinates*)(normalDerivatives + FaceLAC<Basis>::N_BASIS_VALUES);

            double commonMultiplier = det * lambda;
            double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;

            CoordinatesFunctions::translate(transpLocalJacobianMatrix, FaceLAC<Basis>::getLocalGradientsByFace()[localIndex], FaceLAC<Basis>::N_BASIS_VALUES, gradients);

            CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, unitNormal, normalDerivatives);

            FaceLAC<Basis>::computeFlowMatrix(localIndex, normalDerivatives, flowMatrix);
            addDirichletFaceBilinearAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, FaceLAC<Basis>::getMassMatrixies()[localIndex], facePenalty, commonMultiplier, bilinearAdjusments);

            FaceLAC<Basis>::computeFlowVector(normalDerivatives, dirichletValues, flowVector);
            FaceLAC<Basis>::computePowerVector(localIndex, dirichletValues, powerVector);

            addDirichletFaceLinearAdjusments<Basis::N_FUNCTIONS>(powerVector, flowVector, facePenalty, commonMultiplier, linearAdjusments);
        }

        /*
        template<size_t N_FUNCTIONS>
        static void addDirichletFaceLinearAdjusments(const double* massVectorIt,
                                                     const double* flowVectorElementIt,
                                                     const double penaltyCoeff,
                                                     const double commonCoeff,
                                                     double* adjusmentIt)
        {
            for (uint8_t i = 0; i < N_FUNCTIONS; ++i)
            {
                *adjusmentIt += commonCoeff * (penaltyCoeff * (*massVectorIt) - *flowVectorElementIt);

                ++massVectorIt;
                ++flowVectorElementIt;
                ++adjusmentIt;
            }
        }
        */

        static void addDirichletFaceAdjusments(const uint8_t localIndex,
                                               const double lambda,
                                               const double penalty,
                                               const double det,
                                               const Coordinates& unitNormal,
                                               const double conditionValue,
                                               const double transpLocalJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                               void* memoryBuffer,
                                               double bilinearAdjusments[FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                               double linearAdjusments[Basis::N_FUNCTIONS])
        {
            double* flowMatrix = (double*)memoryBuffer;

            double* flowVector = flowMatrix + FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;
            double* normalDerivatives = flowVector + Basis::N_FUNCTIONS;

            Coordinates* gradients = (Coordinates*)(normalDerivatives + FaceLAC<Basis>::N_BASIS_VALUES);

            double commonMultiplier = det * lambda;
            double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(transpMatrixies_[iTriangle], determinant)*/;

            CoordinatesFunctions::translate(transpLocalJacobianMatrix, FaceLAC<Basis>::getLocalGradientsByFace()[localIndex], FaceLAC<Basis>::N_BASIS_VALUES, gradients);

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
                                                    const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                    void* memoryBuffer,
                                                    double bilinearTetrahedronsAdjusments[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                    double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
        {

            double* dirichletValues = (double*)memoryBuffer;

            memoryBuffer = (void*)(dirichletValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps);

            Coordinates unitNormal;
            size_t faceNodesTags[constants::triangle::N_NODES];

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                size_t elementIndex = *faceIndexIt >> 2;
                uint8_t localIndex = *faceIndexIt & 3;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], faceNodesTags);
                double det = computeFaceConditionData(nodes, faceNodesTags, dirichletValues, unitNormal);

                addDirichletFaceAdjusments(localIndex,
                                           lambda,
                                           penalty,
                                           det,
                                           unitNormal,
                                           dirichletValues,
                                           transpLocalJacobianMatrix[elementIndex],
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
                                                    const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                    void* memoryBuffer,
                                                    double bilinearTetrahedronsAdjusments[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                    double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
        {
            Coordinates unitNormal;
            size_t faceNodesTags[constants::triangle::N_NODES];

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                size_t elementIndex = *faceIndexIt >> 2;
                uint8_t localIndex = *faceIndexIt & 3;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], faceNodesTags);
                double det = computeFaceUnitNormal(nodes, faceNodesTags, unitNormal);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValue,
                    transpLocalJacobianMatrix[elementIndex],
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
                                                    const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
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

            size_t facesNodesTags[2][constants::triangle::N_NODES];

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                //double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                elementIndex = *baseFaceIndexesIt >> 2;
                localIndex = *baseFaceIndexesIt & 3;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], facesNodesTags[0]);
                double det = computeFaceConditionData(nodes, facesNodesTags[0], dirichletValues, unitNormal);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValues,
                    transpLocalJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments1Area[elementIndex],
                    linearTetrahedronsAdjusments1Area[elementIndex]);

                elementIndex = *neigbourFaceIndexesIt >> 2;
                localIndex = *neigbourFaceIndexesIt & 3;

                unitNormal.x = -unitNormal.x;
                unitNormal.y = -unitNormal.y;
                unitNormal.z = -unitNormal.z;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], facesNodesTags[1]);

                uint8_t changingIndexes = computeChangingIndexes(facesNodesTags[0], facesNodesTags[1]);
                if (changingIndexes != 0)
                {
                    NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::changeValuesOrder(dirichletValues, Basis::N_FUNCTIONS, changingIndexes);
                }


                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValues,
                    transpLocalJacobianMatrix[elementIndex],
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
                                                    const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                    void* memoryBuffer,
                                                    double bilinearTetrahedronsAdjusments1Area[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                    double bilinearTetrahedronsAdjusments2Area[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                    double linearTetrahedronsAdjusments1Area[][Basis::N_FUNCTIONS],
                                                    double linearTetrahedronsAdjusments2Area[][Basis::N_FUNCTIONS])
        {
            Coordinates unitNormal;

            size_t elementIndex;
            uint8_t localIndex;

            size_t faceNodesTags[constants::triangle::N_NODES];

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                //double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                elementIndex = *baseFaceIndexesIt >> 2;
                localIndex = *baseFaceIndexesIt & 3;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], faceNodesTags);
                double det = computeFaceUnitNormal(nodes, faceNodesTags, unitNormal);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValue,
                    transpLocalJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments1Area[elementIndex],
                    linearTetrahedronsAdjusments1Area[elementIndex]);

                elementIndex = *neigbourFaceIndexesIt >> 2;
                localIndex = *neigbourFaceIndexesIt & 3;

                unitNormal.x = -unitNormal.x;
                unitNormal.y = -unitNormal.y;
                unitNormal.z = -unitNormal.z;

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValue,
                    transpLocalJacobianMatrix[elementIndex],
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
                                                           const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
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
            uint8_t localIndex = *faceIndexIt & 3;

            Coordinates unitNormal;
            size_t faceNodesTags[constants::triangle::N_NODES];
            extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], faceNodesTags);
            computeFaceUnitNormal(nodes, faceNodesTags, unitNormal);

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                elementIndex = *faceIndexIt >> 2;
                localIndex = *faceIndexIt & 3;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], faceNodesTags);
                double det = computeFaceConditionData<Boundary::FunctionCondition::Type::DIRICHLET>(nodes, faceNodesTags, dirichletValues);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValues,
                    transpLocalJacobianMatrix[elementIndex],
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
                                                           const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                           void* memoryBuffer,
                                                           double bilinearTetrahedronsAdjusments[][FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS],
                                                           double linearTetrahedronsAdjusments[][Basis::N_FUNCTIONS])
        {
            if (nFaces == 0)
            {
                return;
            }

            size_t elementIndex = *faceIndexIt >> 2;
            uint8_t localIndex = *faceIndexIt & 3;

            Coordinates unitNormal, normal;
            size_t faceNodesTags[constants::triangle::N_NODES];
            extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], faceNodesTags);
            computeFaceUnitNormal(nodes, faceNodesTags, unitNormal);

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                elementIndex = *faceIndexIt >> 2;
                localIndex = *faceIndexIt & 3;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], faceNodesTags);
                CoordinatesFunctions::computeNormal(nodes, faceNodesTags, normal);
                double det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValue,
                    transpLocalJacobianMatrix[elementIndex],
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
                                                           const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
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
            uint8_t localIndex = *baseFaceIndexesIt & 3;

            Coordinates unitNormal, reverseUnitNormal;
            size_t facesNodesTags[2][constants::triangle::N_NODES];
            extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], facesNodesTags[0]);
            computeFaceUnitNormal(nodes, facesNodesTags[0], unitNormal);
            
            reverseUnitNormal.x = -unitNormal.x;
            reverseUnitNormal.y = -unitNormal.y;
            reverseUnitNormal.z = -unitNormal.z;

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                elementIndex = *baseFaceIndexesIt >> 2;
                localIndex = *baseFaceIndexesIt & 3;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], facesNodesTags[0]);
                double det = computeFaceConditionData<Boundary::FunctionCondition::Type::DIRICHLET>(nodes, facesNodesTags[0], dirichletValues);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValues,
                    transpLocalJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments1Area[elementIndex],
                    linearTetrahedronsAdjusments1Area[elementIndex]);

                elementIndex = *neigbourFaceIndexesIt >> 2;
                localIndex = *neigbourFaceIndexesIt & 3;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], facesNodesTags[1]);

                uint8_t changingIndexes = computeChangingIndexes(facesNodesTags[0], facesNodesTags[1]);
                if (changingIndexes != 0)
                {
                    NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::changeValuesOrder(dirichletValues, 1, changingIndexes);
                }


                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    reverseUnitNormal,
                    dirichletValues,
                    transpLocalJacobianMatrix[elementIndex],
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
                                                           const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
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
            uint8_t localIndex = *baseFaceIndexesIt & 3;

            Coordinates unitNormal, reverseUnitNormal;
            Coordinates normal;
            size_t facesNodesTags[2][constants::triangle::N_NODES];
            extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], facesNodesTags[0]);
            computeFaceUnitNormal(nodes, facesNodesTags[0], unitNormal);

            reverseUnitNormal.x = -unitNormal.x;
            reverseUnitNormal.y = -unitNormal.y;
            reverseUnitNormal.z = -unitNormal.z;

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                double facePenalty = penalty /*/ geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt*/;

                elementIndex = *baseFaceIndexesIt >> 2;
                localIndex = *baseFaceIndexesIt & 3;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], facesNodesTags[0]);
                CoordinatesFunctions::computeNormal(nodes, facesNodesTags[0], normal);
                double det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    unitNormal,
                    dirichletValue,
                    transpLocalJacobianMatrix[elementIndex],
                    memoryBuffer,
                    bilinearTetrahedronsAdjusments1Area[elementIndex],
                    linearTetrahedronsAdjusments1Area[elementIndex]);

                elementIndex = *neigbourFaceIndexesIt >> 2;
                localIndex = *neigbourFaceIndexesIt & 3;

                addDirichletFaceAdjusments(localIndex,
                    lambda,
                    penalty,
                    det,
                    reverseUnitNormal,
                    dirichletValue,
                    transpLocalJacobianMatrix[elementIndex],
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
                                                 const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                 void* memoryBuffer,
                                                 double linearTetrahedronsAdjusments1Area[][Basis::N_FUNCTIONS],
                                                 double linearTetrahedronsAdjusments2Area[][Basis::N_FUNCTIONS])
        {
            double* newmanValues = (double*)memoryBuffer;
            double* newmanVector = newmanValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;

            size_t elementIndex;
            uint8_t localIndex;

            size_t facesNodesTags[2][constants::triangle::N_NODES];

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                elementIndex = *baseFaceIndexesIt >> 2;
                localIndex = *baseFaceIndexesIt & 3;
                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], facesNodesTags[0]);
                double det = computeFaceConditionData<Boundary::FunctionCondition::Type::NEWMAN>(nodes, facesNodesTags[0], newmanValues);

                FaceLAC<Basis>::computePowerVector(localIndex, newmanValues, newmanVector);
                addNemanFaceAdjusments<Basis::N_FUNCTIONS>(det, newmanVector, linearTetrahedronsAdjusments1Area[elementIndex]);

                elementIndex = *neigbourFaceIndexesIt >> 2;
                localIndex = *neigbourFaceIndexesIt & 3;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], facesNodesTags[1]);

                uint8_t changingIndexes = computeChangingIndexes(facesNodesTags[0], facesNodesTags[1]);
                if (changingIndexes != 0)
                {
                    NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::changeValuesOrder(newmanValues, Basis::N_FUNCTIONS, changingIndexes);
                }

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
                                                 const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                 void* memoryBuffer,
                                                 double linearTetrahedronsAdjusmentsArea[][Basis::N_FUNCTIONS])
        {
            double* newmanValues = (double*)memoryBuffer;
            double* newmanVector = newmanValues + NumericalIntegration::Tetrahedron::Gauss<Basis::ORDER + 1>::nSteps;

            size_t elementIndex;
            uint8_t localIndex;

            size_t faceNodesTags[constants::triangle::N_NODES];

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                elementIndex = *faceIndexIt >> 2;
                localIndex = *faceIndexIt & 3;
                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], faceNodesTags);
                double det = computeFaceConditionData<Boundary::FunctionCondition::Type::DIRICHLET>(nodes, faceNodesTags, newmanValues);

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
                                                 const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                 double linearTetrahedronsAdjusments1Area[][Basis::N_FUNCTIONS],
                                                 double linearTetrahedronsAdjusments2Area[][Basis::N_FUNCTIONS])
        {
            const double (*massVectors)[Basis::N_FUNCTIONS] = FaceLAC<Basis>::getMassVectors();
            size_t elementIndex;
            uint8_t localIndex;

            size_t faceNodesTags[constants::triangle::N_NODES];

            Coordinates normal;

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                elementIndex = *baseFaceIndexesIt >> 2;
                localIndex = *baseFaceIndexesIt & 3;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], faceNodesTags);
                CoordinatesFunctions::computeNormal(nodes, faceNodesTags, normal);
                double det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

                addNemanFaceAdjusments<Basis::N_FUNCTIONS>(nerwmanValue, det, massVectors[localIndex], linearTetrahedronsAdjusments1Area[elementIndex]);

                elementIndex = *neigbourFaceIndexesIt >> 2;
                localIndex = *neigbourFaceIndexesIt & 3;

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
                                                 const double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                                 double linearTetrahedronsAdjusmentsArea[][Basis::N_FUNCTIONS])
        {
            const double (*massVectors)[Basis::N_FUNCTIONS] = FaceLAC<Basis>::getMassVectors();
            size_t elementIndex;
            uint8_t localIndex;

            size_t faceNodesTags[constants::triangle::N_NODES];

            Coordinates normal;

            for (size_t faceIndex = 0; faceIndex < nFaces; ++faceIndex)
            {
                elementIndex = *faceIndexIt >> 2;
                localIndex = *faceIndexIt & 3;

                extractFaceNodesTags(localIndex, tetrahedronsNodesTags[elementIndex], faceNodesTags);
                CoordinatesFunctions::computeNormal(nodes, faceNodesTags, normal);
                double det = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

                addNemanFaceAdjusments<Basis::N_FUNCTIONS>(nerwmanValue, det, massVectors[localIndex], linearTetrahedronsAdjusmentsArea[elementIndex]);

                ++faceIndexIt;
            }
        }

        /*
        static void processInteriorFacesOnInterface(const size_t* baseFaceIndexesIt,
                                                    const size_t* neigbourFaceIndexesIt,
                                                    const size_t nInteriorFaces,
                                                    const double lambda,
                                                    const double penalty,
                                                    const Coordinates nodes[],
                                                    const size_t tetrahedronsNodesTags[][constants::tetrahedron::N_NODES],
                                                    void* calculationBuffer,
                                                    double transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
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
                double facePenalty = penalty // geom_funct::computeTriangleDiametr(*transpMatrixiesIt, determinant); ++transpMatrixiesIt/;

                size_t elementIndex0 = *baseFaceIndexesIt >> 2;
                uint8_t localIndex0 = *baseFaceIndexesIt & 3;
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

                CoordinatesFunctions::translate(transpLocalJacobianMatrix[elementIndex0],
                    gradientsByFace[localIndex0],
                    NumericalIntegration::Triangle::Gauss<Basis::ORDER + 1>::nSteps,
                    gradients);

                CoordinatesFunctions::coomputeDirectionalDerivative(gradients, FaceLAC<Basis>::N_BASIS_VALUES, normal, normalDerivatives);

                FaceLAC<Basis>::computeFlowMatrix(localIndex0, normalDerivatives, flowMatrix);
                addInteriorFaceAdjusmets<Basis::N_FUNCTIONS>(flowMatrix, massMatrixByFace[localIndex0], penalty, commonMultiplier, bilinearTetrahedronsAdjusments1[elementIndex0]);


                size_t elementIndex1 = *neigbourFaceIndexesIt >> 2;
                uint8_t localIndex1 = *neigbourFaceIndexesIt & 3;
                ++neigbourFaceIndexesIt;

                FaceLAC<Basis>::computeFlowMatrix(localIndex1, normalDerivatives, crossFlowMatrix1);

                normal.x = -normal.x;
                normal.y = -normal.y;
                normal.z = -normal.z;

                CoordinatesFunctions::translate(transpLocalJacobianMatrix[elementIndex1],
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
        */


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
                                      double  transpLocalJacobianMatrix[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                      double* const regionsBilinearAdjusments[],
                                      double* const regionsLinearAdjusmenst[])
        {

            for (unsigned int i = 0; i < nBoundaries; ++i)
            {
                if (boundaryIt->condition == nullptr)
                {
                    continue;
                }
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
                                transpLocalJacobianMatrix,
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
                                transpLocalJacobianMatrix,
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
                                transpLocalJacobianMatrix,
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
                                transpLocalJacobianMatrix,
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
                                                     transpLocalJacobianMatrix,
                                                     (double(*)[Basis::N_FUNCTIONS])regionsLinearAdjusmenst[boundaryIt->regionsIndexes[0]]);

                        break;
                    }
                    case  Boundary::Condition::MacroType::NEWMAN_F:
                    {
                        computeNewmanFacesAdjusments(boundaryBasesFacesIndexes,
                            nFaces,
                            nodes,
                            tetrahedronsNodesTags,
                            transpLocalJacobianMatrix,
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
                                transpLocalJacobianMatrix,
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
                                transpLocalJacobianMatrix,
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
                                transpLocalJacobianMatrix,
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
                                transpLocalJacobianMatrix,
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
                            transpLocalJacobianMatrix,
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
                            transpLocalJacobianMatrix,
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

                ++boundaryIt;
                ++surfaceFacesStartIndexeIt;
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
                                                 double  transpLocalJacobianMatrixes[][LocalCoordinates3D::COUNT * Coordinates::COUNT],
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

                unsigned int nFragments;
                gmsh::model::surfaces::createFragmentModel(nonconformInterfaces->sidesTags, "FragmentModel", interfaceFragmentsTrianglesIndexesIt, &nFragments);
                nAllFragments += nFragments;
                *interfacFragmentsCountIt = nFragments;

                gmsh::model::mesh::generate(constants::DIMENSION_2D);
                size_t nNodes = gmsh::model::mesh::getNodesCount();
                Coordinates* fragmentsModelsNodes = (Coordinates*)malloc((nNodes + 1) * sizeof(Coordinates));
                gmsh::model::mesh::getNodes((double(*)[3])fragmentsModelsNodes);
                size_t* trianglesStartIndexes = (size_t*)malloc((nFragments + 1) * sizeof(size_t));
                gmsh::model::mesh::getSurfacesTrianglesStartIndexes(trianglesStartIndexes);
                size_t* triangleNodesTagsIt = (size_t*)malloc(trianglesStartIndexes[nFragments] * constants::triangle::N_NODES * sizeof(size_t));
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
                    transpLocalJacobianMatrixes,
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
            const unsigned int* interfaceFragmentSide2TriangleIndexIt = interfaceFragmentTriangleIndexes + nCrossElements;

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
                              unsigned int* interfacesFragmentsCount,
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
                                                      interfacesFragmentsCount[i],
                                                      tetBaseFacesIndexes + surfacesFacesStartIndexes[nonconformInterfaces[i].sideIndexes[0]],
                                                      tetBaseFacesIndexes + surfacesFacesStartIndexes[nonconformInterfaces[i].sideIndexes[1]],
                                                      interfacesFragmentsTrianglesIndexes[i],
                                                      LinearLagrangeBasis::N_FUNCTIONS,
                                                      startIndexes,
                                                      tripletIt);
            }

        }

        static size_t findFirstStartTetrahedronTag(size_t* regionStartTetrahedronsTagIt)
        {
            while (*regionStartTetrahedronsTagIt == 0);
            return *regionStartTetrahedronsTagIt;
        }

        const int MAX_ITERATIONS = 10'000;

        void initSolver()
        {
            ElementLAC<Basis>::init();
            FaceLAC<Basis>::init();
        }

        void finalizeSolver()
        {
            ElementLAC<Basis>::finalize();
            FaceLAC<Basis>::finalize();
        }

        double solveInitialIteration(const unsigned int nRegions,
                                     const MaterialPhase* const regionsMaterialPhases[],
                                     const Boundary boundaries[],
                                     const unsigned int nBoundaries,
                                     const NonconformInterface nonconformInterfaces[],
                                     const unsigned int nNonconformInterfaces,
                                     const double t,
                                     const double dt,
                                     const double penalty,
                                     void* calculationBuffer,
                                     void* additionalBuffer,
                                     Solution* solutionIt)
        {
            T = t;
            size_t* regionsStartTetrahedronsIndexes = (size_t*)additionalBuffer;
            size_t* surfacesFacesStartIndexes = regionsStartTetrahedronsIndexes + nRegions + 1;
            size_t* regionsStartTetrahedronsTags = surfacesFacesStartIndexes + nBoundaries + 1;
            size_t* regionsInteriorFacesStartIndexes = regionsStartTetrahedronsTags + nRegions + 1;

            double** regionsBilinearAdjusments = (double**)(regionsInteriorFacesStartIndexes + nRegions + 1);
            double** regionsLinearAdjusments = regionsBilinearAdjusments + nRegions;

            unsigned int** interfacesFragmentsTrianglesIndexes = (unsigned int**)(regionsLinearAdjusments + nRegions);
            unsigned int* interfacFragmentsCount = (unsigned int*)(interfacesFragmentsTrianglesIndexes + nNonconformInterfaces);
            double** nonconformCrossAdjusments = (double**)(interfacFragmentsCount + nNonconformInterfaces);

            size_t nNodes = gmsh::model::mesh::getNodesCount();
            Coordinates* nodes = (Coordinates*)malloc((nNodes + 1) * sizeof(Coordinates));
            gmsh::model::mesh::getNodes((double(*)[3])nodes);

            gmsh::model::mesh::getSurfacesTrianglesStartIndexes(surfacesFacesStartIndexes);
            gmsh::model::mesh::getRegionsTetrahedronsStartIndexes(regionsStartTetrahedronsIndexes);
            size_t nTetrahedrons = regionsStartTetrahedronsIndexes[nRegions];
            double* bilinearAdjusments = (double*)malloc(nTetrahedrons * ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS * sizeof(double));
            size_t nTetNodes = nTetrahedrons * constants::tetrahedron::N_NODES;
            size_t tetFecesBufferSize = nTetrahedrons * constants::tetrahedron::N_FACES + surfacesFacesStartIndexes[nBoundaries];
            size_t* tetFacesIndexes = (size_t*)malloc(tetFecesBufferSize * sizeof(size_t));
            size_t* tetBaseFacesIndexes = tetFacesIndexes;
            size_t* tetNeigbourFacesIndexes = tetFacesIndexes + tetFecesBufferSize / 2;
            size_t* tetNodesTag = (size_t*)malloc(nTetNodes * sizeof(size_t));

            gmsh::model::mesh::getTetrahedrons(regionsStartTetrahedronsTags, tetNodesTag, regionsInteriorFacesStartIndexes, tetBaseFacesIndexes, tetNeigbourFacesIndexes);

            size_t nInteriorFaces = regionsInteriorFacesStartIndexes[nRegions] - regionsInteriorFacesStartIndexes[0];

            size_t nBilinearAdjusments = nTetrahedrons * ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS + nInteriorFaces * 2 * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;

            size_t nDOFs = nTetrahedrons * Basis::N_FUNCTIONS;

            Eigen::VectorXd initial(nDOFs);
            Eigen::VectorXd f(nDOFs);

            double* initialPtr = initial.data();
            double* fPtr = f.data();

            double (*transpLocalJacobianMatrixies)[LocalCoordinates3D::COUNT * Coordinates::COUNT] = 
                (double(*)[LocalCoordinates3D::COUNT * Coordinates::COUNT])malloc(nTetrahedrons * LocalCoordinates3D::COUNT * Coordinates::COUNT * sizeof(double));

            double* bilinearAdjusmentsIt = bilinearAdjusments;
            double* initialIt = initial.data();
            double* fIt = f.data();


            const MaterialPhase* const* materialPhasePtrIt = regionsMaterialPhases;
            const size_t *regionStartTetrahedronIndexIt = regionsStartTetrahedronsIndexes;
            double (*transpLocalJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT] = transpLocalJacobianMatrixies;

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
                           transpLocalJacobianMatrixIt,
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
                             transpLocalJacobianMatrixIt,
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
                                       transpLocalJacobianMatrixIt,
                                       regionsBilinearAdjusments,
                                       interfacesFragmentsTrianglesIndexes,
                                       interfacFragmentsCount,
                                       nonconformCrossAdjusments);

            nBilinearAdjusments += nFragments * 2 * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;

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
            free(triplets);
            printf("Building matrix finished\n");

            printf("\nSolving slae...\n");

#if DEBUG_SHOW_SLAE
            std::cout << A << std::endl;
            std::cout << "Intitial" << std::endl;
            std::cout << initial << std::endl;
            std::cout << "Power vector" << std::endl;
            std::cout << f << std::endl;
#endif
            Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double, 0, Eigen::Index>> solver;
            //Eigen::BiCGSTAB<Eigen::SparseMatrix<double, 0, Eigen::Index>> solver;
            solver.setMaxIterations(MAX_ITERATIONS);
            //Eigen::Index maxIter = solver.maxIterations();
            //double tol = solver.tolerance();
            solver.compute(A);
            solutionIt->_DOFs = solver.solveWithGuess(f, initial);
            solutionIt->_DOFsPtr = solutionIt->_DOFs.data() - findFirstStartTetrahedronTag(regionsStartTetrahedronsTags) * Basis::N_FUNCTIONS;

            printf("Solving slae finished\n");
            printf("Estimated error %e; Iterations count %zu\n", solver.error(), solver.iterations());

#if DEBUG_SHOW_SLAE
            std::cout << "Solution" << std::endl;
            std::cout << solutionIt->_DOFs << std::endl;
#endif

            //std::cout << "Solution" << std::endl;
            //std::cout << solutionIt->_DOFs << std::endl;
            //std::cout << "Solution" << std::endl;
            //std::cout << solutionIt->_DOFs << std::endl;

            double err = computeEnergyError((double(*)[Basis::N_FUNCTIONS])(solutionIt->_DOFs.data()),
                                            nRegions,
                                            regionsMaterialPhases,
                                            regionsStartTetrahedronsIndexes,
                                            nodes,
                                            (size_t(*)[constants::tetrahedron::N_NODES])tetNodesTag,
                                            t);

            printf("Energy error: %lf\n", err);

            for (unsigned int i = 0; i < nNonconformInterfaces; ++i)
            {
                free(interfacesFragmentsTrianglesIndexes[i]);
                free(nonconformCrossAdjusments[i]);
            }

            free(nodes);
            free(tetFacesIndexes);
            free(tetNodesTag);
            free(bilinearAdjusments);
            free(transpLocalJacobianMatrixies);
            free(crossBilinearAdjusments);

            return err;
        }

        double solveIteration(const unsigned int nRegions,
                              const MaterialPhase* const regionsMaterialPhases[],
                              const Boundary boundaries[],
                              const unsigned int nBoundaries,
                              const NonconformInterface nonconformInterfaces[],
                              const unsigned int nNonconformInterfaces,
                              const double t,
                              const double dt,
                              const double penalty,
                              const char* currentSolutionModelName,
                              const Solution& previousSolution,
                              const char* previousSolutionModelName,
                              void* calculationBuffer,
                              void* additionalBuffer,
                              Solution* solutionIt)
        {
            T = t;
            size_t* regionsStartTetrahedronsIndexes = (size_t*)additionalBuffer;
            size_t* surfacesFacesStartIndexes = regionsStartTetrahedronsIndexes + nRegions + 1;
            size_t* regionsStartTetrahedronsTags = surfacesFacesStartIndexes + nBoundaries + 1;
            size_t* regionsInteriorFacesStartIndexes = regionsStartTetrahedronsTags + nRegions + 1;

            double** regionsBilinearAdjusments = (double**)(regionsInteriorFacesStartIndexes + nRegions + 1);
            double** regionsLinearAdjusments = regionsBilinearAdjusments + nRegions;

            unsigned int** interfacesFragmentsTrianglesIndexes = (unsigned int**)(regionsLinearAdjusments + nRegions);
            unsigned int* interfacFragmentsCount = (unsigned int*)(interfacesFragmentsTrianglesIndexes + nNonconformInterfaces);
            double** nonconformCrossAdjusments = (double**)(interfacFragmentsCount + nNonconformInterfaces);

            size_t nNodes = gmsh::model::mesh::getNodesCount();
            Coordinates* nodes = (Coordinates*)malloc((nNodes + 1) * sizeof(Coordinates));
            gmsh::model::mesh::getNodes((double(*)[3])nodes);

            gmsh::model::mesh::getSurfacesTrianglesStartIndexes(surfacesFacesStartIndexes);
            gmsh::model::mesh::getRegionsTetrahedronsStartIndexes(regionsStartTetrahedronsIndexes);
            size_t nTetrahedrons = regionsStartTetrahedronsIndexes[nRegions];
            double* bilinearAdjusments = (double*)malloc(nTetrahedrons * ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS * sizeof(double));
            size_t nTetNodes = nTetrahedrons * constants::tetrahedron::N_NODES;
            size_t tetFecesBufferSize = nTetrahedrons * constants::tetrahedron::N_FACES + surfacesFacesStartIndexes[nBoundaries];
            size_t* tetFacesIndexes = (size_t*)malloc(tetFecesBufferSize * sizeof(size_t));
            size_t* tetBaseFacesIndexes = tetFacesIndexes;
            size_t* tetNeigbourFacesIndexes = tetFacesIndexes + tetFecesBufferSize / 2;
            size_t* tetNodesTag = (size_t*)malloc(nTetNodes * sizeof(size_t));

            gmsh::model::mesh::getTetrahedrons(regionsStartTetrahedronsTags, tetNodesTag, regionsInteriorFacesStartIndexes, tetBaseFacesIndexes, tetNeigbourFacesIndexes);

            size_t nInteriorFaces = regionsInteriorFacesStartIndexes[nRegions] - regionsInteriorFacesStartIndexes[0];

            size_t nBilinearAdjusments = nTetrahedrons * ElementLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS + nInteriorFaces * 2 * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;

            size_t nDOFs = nTetrahedrons * Basis::N_FUNCTIONS;

            Eigen::VectorXd initial(nDOFs);
            Eigen::VectorXd f(nDOFs);

            double* initialPtr = initial.data();
            double* fPtr = f.data();

            double (*transpLocalJacobianMatrixies)[LocalCoordinates3D::COUNT * Coordinates::COUNT] =
                (double(*)[LocalCoordinates3D::COUNT * Coordinates::COUNT])malloc(nTetrahedrons * LocalCoordinates3D::COUNT * Coordinates::COUNT * sizeof(double));

            double* bilinearAdjusmentsIt = bilinearAdjusments;
            double* initialIt = initial.data();
            double* fIt = f.data();


            const MaterialPhase* const* materialPhasePtrIt = regionsMaterialPhases;
            const size_t* regionStartTetrahedronIndexIt = regionsStartTetrahedronsIndexes;
            double (*transpLocalJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT] = transpLocalJacobianMatrixies;

            double* crossBilinearAdjusments = (double*)malloc(nInteriorFaces * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS * sizeof(double));

            //gmsh::model::setCurrent(previousSolutionModelName);

            processRegions(nRegions,
                           regionsMaterialPhases,
                           dt,
                           penalty,
                           currentSolutionModelName,
                           previousSolutionModelName,
                           previousSolution,
                           calculationBuffer,
                           regionStartTetrahedronIndexIt,
                           nodes,
                           (size_t(*)[constants::tetrahedron::N_NODES])tetNodesTag,
                           regionsInteriorFacesStartIndexes,
                           tetBaseFacesIndexes,
                           tetNeigbourFacesIndexes,
                           transpLocalJacobianMatrixIt,
                           bilinearAdjusmentsIt,
                           regionsBilinearAdjusments,
                           crossBilinearAdjusments,
                           initialIt,
                           fIt,
                           regionsLinearAdjusments);

            //gmsh::model::setCurrent(currentSolutionModelName);

            //std::cout << "Power vector" << std::endl;
            //std::cout << f << std::endl;

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
                transpLocalJacobianMatrixIt,
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
                transpLocalJacobianMatrixIt,
                regionsBilinearAdjusments,
                interfacesFragmentsTrianglesIndexes,
                interfacFragmentsCount,
                nonconformCrossAdjusments);

            nBilinearAdjusments += nFragments * 2 * FaceLAC<Basis>::N_LOCAL_MATRIX_ELEMENTS;

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
            free(triplets);
            printf("Building matrix finished\n");

            printf("\nSolving slae...\n");

#if DEBUG_SHOW_SLAE
            std::cout << A << std::endl;
            std::cout << "Intitial" << std::endl;
            std::cout << initial << std::endl;
            std::cout << "Power vector" << std::endl;
            std::cout << f << std::endl;
#endif
            //std::cout << "Intitial" << std::endl;
            //std::cout << initial << std::endl;

            Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double, 0, Eigen::Index>> solver;
            //Eigen::BiCGSTAB<Eigen::SparseMatrix<double, 0, Eigen::Index>> solver;
            solver.setMaxIterations(MAX_ITERATIONS);
            //Eigen::Index maxIter = solver.maxIterations();
            //double tol = solver.tolerance();
            solver.compute(A);
            solutionIt->_DOFs = solver.solveWithGuess(f, initial);
            solutionIt->_DOFsPtr = solutionIt->_DOFs.data() - findFirstStartTetrahedronTag(regionsStartTetrahedronsTags) * Basis::N_FUNCTIONS;

            printf("Solving slae finished\n");
            printf("Estimated error %e; Iterations count %zu\n", solver.error(), solver.iterations());

#if DEBUG_SHOW_SLAE
            std::cout << "Solution" << std::endl;
            std::cout << solutionIt->_DOFs << std::endl;
#endif

            //std::cout << "Solution" << std::endl;
            //std::cout << solutionIt->_DOFs << std::endl;

            double err = computeEnergyError((double(*)[Basis::N_FUNCTIONS])(solutionIt->_DOFs.data()),
                                            nRegions,
                                            regionsMaterialPhases,
                                            regionsStartTetrahedronsIndexes,
                                            nodes,
                                            (size_t(*)[constants::tetrahedron::N_NODES])tetNodesTag,
                                            t);

            printf("Energy error: %lf\n", err);

            for (unsigned int i = 0; i < nNonconformInterfaces; ++i)
            {
                free(interfacesFragmentsTrianglesIndexes[i]);
                free(nonconformCrossAdjusments[i]);
            }

            free(nodes);
            free(tetFacesIndexes);
            free(tetNodesTag);
            free(bilinearAdjusments);
            free(transpLocalJacobianMatrixies);
            free(crossBilinearAdjusments);

            return err;
        }

        void getFrontNodes(const int frontTag, int** frontNodesTags, unsigned int* nFrontNodes)
        {
            unsigned int nEdges, nEmbeddedEdges, nEmbeddedNodes;
            gmsh::model::getSurfaceChildrenCounts(frontTag, &nEdges, &nEmbeddedEdges, &nEmbeddedNodes);
            int* nodesBuffer = (int*)malloc(((nEdges + nEmbeddedEdges) * 2 + nEmbeddedNodes) * sizeof(int));
            gmsh::model::getSurfaceChildrenNodes(frontTag, nodesBuffer + nEmbeddedNodes + 2 * nEmbeddedEdges, nodesBuffer + nEmbeddedNodes, nodesBuffer);

            std::set<int> nodesTagSet(nodesBuffer, nodesBuffer + (nEdges + nEmbeddedEdges) * 2 + nEmbeddedNodes);
            unsigned int _nFrontNodes = unsigned int(nodesTagSet.size());
            nodesBuffer = (int*)realloc(nodesBuffer, _nFrontNodes * sizeof(int));

            std::set<int>::const_iterator nodeTagSetIt = nodesTagSet.begin();
            int* nodeTagOutIt = nodesBuffer;

            for (unsigned int i = 0; i < _nFrontNodes; ++i)
            {
                *nodeTagOutIt = *nodeTagSetIt;
                ++nodeTagOutIt;
                ++nodeTagSetIt;
            }

            *nFrontNodes = _nFrontNodes;
            *frontNodesTags = nodesBuffer;
        }

        void computeFrontNodesInitialMovements(MaterialPhase solidMaterialsPhase,
                                               MaterialPhase liquidMaterialsPhase,
                                               const double dt,
                                               const double latentHeat,
                                               const Coordinates* normalIt,
                                               const unsigned int nFrontNodes,
                                               const Coordinates* frontNodesIt,
                                               Coordinates* movementIt)
        {
            double commonMultiplier = dt / latentHeat;

            Coordinates gradientSolidSide, gradientLiquidSide;
            for (size_t i = 0; i < nFrontNodes; ++i)
            {
                initialConditionGradient<MaterialPhase::SOLID>(*frontNodesIt, gradientSolidSide);
                initialConditionGradient<MaterialPhase::LIQUID>(*frontNodesIt, gradientLiquidSide);
                double flowSolidSide = (gradientSolidSide.x * normalIt->x + gradientSolidSide.y * normalIt->y + gradientSolidSide.z * normalIt->z);
                double flowLiquidSide = -(gradientLiquidSide.x * normalIt->x + gradientLiquidSide.y * normalIt->y + gradientLiquidSide.z * normalIt->z);
                double moveCoeff = commonMultiplier * (solidMaterialsPhase.thermalConductivity * flowSolidSide + liquidMaterialsPhase.thermalConductivity * flowLiquidSide);
                movementIt->x = moveCoeff * normalIt->x;
                movementIt->y = moveCoeff * normalIt->y;
                movementIt->z = moveCoeff * normalIt->z;

                ++movementIt;
                ++frontNodesIt;
                ++normalIt;
            }
        }

        uint8_t findClosestTetrahedronPointIndex(const Coordinates& point, const Coordinates* tetrahedronPointIt)
        {
            double pointsMetric, pointsMinMetric;
            pointsMinMetric = abs(point.x - tetrahedronPointIt->x) + abs(point.y - tetrahedronPointIt->y) + abs(point.z - tetrahedronPointIt->z);
            uint8_t iClosestPoint = 0;

            for (uint8_t i = 1; i < constants::tetrahedron::N_NODES; ++i)
            {
                ++tetrahedronPointIt;

                pointsMetric = abs(point.x - tetrahedronPointIt->x) + abs(point.y - tetrahedronPointIt->y) + abs(point.z - tetrahedronPointIt->z);
                if (pointsMetric < pointsMinMetric)
                {
                    pointsMinMetric = pointsMetric;
                    iClosestPoint = i;
                }
            }

            return iClosestPoint;
        }
        void computeFrontNodesMovements(const Solution& solution,
                                        int solidRegionTag,
                                        int liquidRegionTag,
                                        MaterialPhase solidMaterialsPhase,
                                        MaterialPhase liquidMaterialsPhase,
                                        const double dt,
                                        const double latentHeat,
                                        const Coordinates* normalIt,
                                        const unsigned int nFrontNodes,
                                        const Coordinates* frontNodesIt,
                                        Coordinates* movements)
        {
            LocalCoordinates3D localGradient;
            Coordinates gradient;

            double commonMultiplier = dt / latentHeat;
            LocalCoordinates3D tetrahedronMainLocalPoints[4] = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0} };
            uint8_t iLocalNode;

            size_t elementsTagsStaticBuffer[ELEMENTS_TAGS_BUFFER_MAX_SIZE];
            size_t* elementsTags = elementsTagsStaticBuffer;
            uint8_t elementsTagsBufferSize = ELEMENTS_TAGS_BUFFER_MAX_SIZE;
            bool usingDynamicBuffer = false;

            size_t tetrahedronNodeTags[constants::tetrahedron::N_NODES];
            Coordinates tetrahedronNodes[constants::tetrahedron::N_NODES];
            size_t nodeTag;

            uint8_t nElements;
            int entityTag;

            double jacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT];
            double transpLocalJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT];

            Coordinates gradientSolidSideSum, gradientLiquidSideSum;
            uint8_t nElementSolidSide, nElementLiquidSide;

            for (size_t i = 0; i < nFrontNodes; ++i)
            {
                gmsh::model::mesh::getNodeByCoordinates(frontNodesIt->x, frontNodesIt->y, frontNodesIt->z, &nodeTag);

                if (gmsh::model::mesh::get3DElementsByCoordinates(frontNodesIt->x, frontNodesIt->y, frontNodesIt->z, elementsTagsBufferSize, elementsTags, &nElements) == 1)
                {
                    if (usingDynamicBuffer)
                    {
                        free(elementsTags);
                    }
                    else
                    {
                        usingDynamicBuffer = true;
                    }
                    elementsTagsBufferSize = nElements;
                    elementsTags = (size_t*)malloc(ELEMENTS_TAGS_BUFFER_MAX_SIZE * sizeof(size_t));
                    gmsh::model::mesh::get3DElementsByCoordinates(frontNodesIt->x, frontNodesIt->y, frontNodesIt->z, elementsTagsBufferSize, elementsTags, &nElements);
                }

                gradientSolidSideSum.x = 0;
                gradientSolidSideSum.y = 0;
                gradientSolidSideSum.z = 0;

                gradientLiquidSideSum.x = 0;
                gradientLiquidSideSum.y = 0;
                gradientLiquidSideSum.z = 0;
                
                nElementSolidSide = 0;
                nElementLiquidSide = 0;

                for (uint8_t j = 0; j < nElements; ++j)
                {
                    gmsh::model::mesh::get3DElement(elementsTags[j], tetrahedronNodeTags, (double*)tetrahedronNodes, &entityTag);
                    iLocalNode = findClosestTetrahedronPointIndex(*frontNodesIt, tetrahedronNodes);

                    solution.compute(elementsTags[j], tetrahedronMainLocalPoints[iLocalNode], localGradient);

                    double det = CoordinatesFunctions::computeJacobianMatrix(tetrahedronNodes, jacobianMatrix);
                    //CoordinatesFunctions::computeTranspLocalJacobianMatrix(jacobianMatrix, det, transpLocalJacobianMatrix);
                    CoordinatesFunctions::computeTranspLocalJacobianMatrixTemplate(jacobianMatrix, transpLocalJacobianMatrix);

                    CoordinatesFunctions::translate(transpLocalJacobianMatrix, localGradient, gradient);
                    gradient.x /= det;
                    gradient.y /= det;
                    gradient.z /= det;

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

                double flowSolidSide = (gradientSolidSideSum.x * normalIt->x + gradientSolidSideSum.y * normalIt->y + gradientSolidSideSum.z * normalIt->z) / nElementSolidSide;
                double flowLiquidSide = -(gradientLiquidSideSum.x * normalIt->x + gradientLiquidSideSum.y * normalIt->y + gradientLiquidSideSum.z * normalIt->z) / nElementLiquidSide;
                
                double moveCoeff = commonMultiplier * (solidMaterialsPhase.thermalConductivity * flowSolidSide + liquidMaterialsPhase.thermalConductivity * flowLiquidSide);

                //moveCoeff = 0.1;
                movements->x += moveCoeff * normalIt->x;
                movements->y += moveCoeff * normalIt->y;
                movements->z += moveCoeff * normalIt->z;

                ++movements;
                ++frontNodesIt;
                ++normalIt;
            }
            if (usingDynamicBuffer)
            {
                free(elementsTags);
            }
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

    void DG::Solution::clear()
    {
        _DOFs.resize(0);
        _DOFsPtr = nullptr;
    }

    const double* DG::Solution::getDOFs() const
    {
        return _DOFs.data();
    }
    
}

