#define _CRT_SECURE_NO_WARNINGS
#define GMSH_DLL
#include "gmsh.h"
#include <iostream>
#include <cstdio>
#include "DG.h"
#include "Boundary.h"
#include "Model.h"

const char timeMeshFileName[] = "TimeMesh.txt";
const char gmshInitialModelName[] = "Stefan";
//const char gmshInitialModelFileName[] = "./models/Stefan.geo";
//const char gmshInitialModelFileName[] = "./models/Stefan_conform.xao";

const char gmshBaseModelFileName[] = "./models/FrontPoints.geo";
const char gmshMergeModelFileName[] = "./models/FrontPointsMerge.geo";

enum Error {
	OPENING_TIME_FILE = 1,
	READING_TIME_PARAMETERS,
	OPENING_PARAMETER_FILE,
	READING_PARAMETER,
	READING_MODEL
};

static uint8_t readTimeParameters(double* tMin, double* tMax, size_t* nTIntervals)
{
	FILE* timeFile = fopen(timeMeshFileName, "r");

	if (timeFile == NULL)
	{
		printf("Error: unable to open file %s \n", timeMeshFileName);
		return Error::OPENING_TIME_FILE;
	}

	if (fscanf(timeFile, "%lf %lf %zu", tMin, tMax, nTIntervals) != 3)
	{
		printf("Error: unable to read time parameters in file %s\n", timeMeshFileName);
		return Error::READING_TIME_PARAMETERS;
	}

	fclose(timeFile);
	return 0;
}

static uint8_t readParameter(const char parameterName[], const char parameterFileName[], double* parameter)
{
	FILE* penaltyFile = fopen(parameterFileName, "r");

	if (penaltyFile == NULL)
	{
		printf("Error: unable to open file %s \n", parameterFileName);
		return Error::OPENING_PARAMETER_FILE;
	}

	if (fscanf(penaltyFile, "%lf", parameter) != 1)
	{
		printf("Error: unable to read %s in file %s\n", parameterName, parameterFileName);
		return Error::READING_PARAMETER;
	}

	fclose(penaltyFile);
	return 0;
}

static uint8_t test()
{
	double penalty;
	uint8_t err = readParameter("Penalty", "Penalty.txt", &penalty);
	if (err)
	{
		return err;
	}

	gmsh::initialize();
	//const char gmshFilePath[] = "./models/two_elements_nonconformal.geo";
	//const char gmshFilePath[] = "./models/Stefan_conform.geo";
	//const char gmshFilePath[] = "./models/BoxGeo.geo";
	const char gmshFilePath[] = "./models/StefanConfromOCC.geo";
	//const char gmshFilePath[] = "./models/box_occ.geo";
	//const char gmshFilePath[] = "./models/one_element_occ.geo";
	const char gmshModelName[] = "Test";
	//gmsh::open(gmshFilePath);
	gmsh::model::add(gmshModelName);

	gmsh::merge(gmshFilePath);
	//gmsh::open("./Models/one_element_occ.geo");
	//gmsh::open("./Models/two_elements_occ.geo");
	//gmsh::open("./Models/box_occ.geo");
	//gmsh::open("./models/two_elements_nonconformal.geo");
	

	Model model;
	if (model.initilizeByCurrentGMSHModel() != Model::Error::NO_ERRORS)
	{
		gmsh::finalize();
		return 1;
	}
	//gmsh::model::occ::synchronize();
	//gmsh::model::geo::synchronize();
	gmsh::vectorpair vectp = { {0, 7}, {0, 2}, {0, 6}, {0, 3} };
	//gmsh::vectorpair vectp = { {2, 1}};
	/*vectp.resize(nFrontNodes);
	for (size_t j = 0; j < nFrontNodes; ++j)
	{
		vectp[j] = { 0, frontNodeTags[j] };
	}*/
	//int frontNodeTags[] = {7, 2, 6, 3};
	//Coordinates frontNodes[] = { {1.1, 0, 1}, {1.1, 0, 0}, {1.1, 1, 1}, {1.1, 1, 0} };
	//gmsh::model::setNodes(frontNodeTags, 4, (double(*)[3])frontNodes);
	//gmsh::model::occ::translate(vectp, 0.1, 0.0, 0.0);
	//gmsh::model::occ::synchronize();
	//gmsh::model::mesh::generate();
	
	size_t* modelMemoryBuffer = (size_t*)malloc((3 * (model.nRegions + 1) + model.nBoundaries + 1) * sizeof(size_t) +
		2 * model.nRegions * sizeof(double*) +
		model.nNonconformInterfaces * sizeof(unsigned int*) +
		model.nNonconformInterfaces * sizeof(unsigned int) +
		model.nNonconformInterfaces * sizeof(double*));

	double* calculationBuffer = (double*)malloc(1024 * sizeof(double));

	DG::StefanTask::initSolver();
	gmsh::model::mesh::generate();

	DG::Solution solution;
	DG::StefanTask::solveInitialIteration(model.nRegions,
		model.regionsMaterialPhases,
		model.boundaries,
		model.nBoundaries,
		model.nonconformInterfaces,
		model.nNonconformInterfaces,
		1.0,
		1.0,
		penalty,
		calculationBuffer,
		modelMemoryBuffer,
		&solution);

	
	const char gmshModel1Name[] = "Test1";
	gmsh::model::add(gmshModel1Name);
	gmsh::merge(gmshFilePath);
	gmsh::model::mesh::generate();
	
	DG::Solution solution2;
	DG::StefanTask::solveIteration(model.nRegions,
								   model.regionsMaterialPhases,
								   model.boundaries,
								   model.nBoundaries,
								   model.nonconformInterfaces,
								   model.nNonconformInterfaces,
								   2.0,
								   1.0,
								   penalty,
								   gmshModelName,
								   solution,
								   gmshModel1Name,
								   calculationBuffer,
								   modelMemoryBuffer,
								   &solution2);
	
	solution2.clear();
	//*/
	//solution.clear();


	//gmsh::fltk::run();
	gmsh::finalize();
}

int main()
{
	//test();
	//return 0;
	double tMin, tMax;
	size_t nTIntervals;

	uint8_t err = readTimeParameters(&tMin, &tMax, &nTIntervals);
	if (err)
	{
		return err;
	}

	double penalty;
	err = readParameter("Penalty", "Penalty.txt", &penalty);
	if (err)
	{
		return err;
	}

	double latentHeat;
	err = readParameter("LatentHeat", "LatentHeat.txt", &latentHeat);
	if (err)
	{
		return err;
	}

	gmsh::initialize();
	//test();
	//return 0;
	gmsh::model::add(gmshInitialModelName);
	gmsh::merge(gmshBaseModelFileName);

	int* frontsPointsTags;
	Coordinates* baseFrontsPoints;
	unsigned int* frontsPointsStartIndexes;
	unsigned int nFronts;
	Model::getFronts(frontsPointsTags, baseFrontsPoints, frontsPointsStartIndexes, nFronts);

	unsigned int nFrontsNodes = frontsPointsStartIndexes[nFronts];
	
	Coordinates* frontsBuffer = (Coordinates*)malloc(3 * nFrontsNodes * sizeof(Coordinates));
	Coordinates* frontPoints = frontsBuffer;
	Coordinates* frontPointsMovements = frontPoints + nFrontsNodes;
	Coordinates* frontNormals = frontPointsMovements + nFrontsNodes;

	gmsh::merge(gmshMergeModelFileName);

	Model model;
	if (model.initilizeByCurrentGMSHModel() != Model::Error::NO_ERRORS)
	{
		gmsh::finalize();
		return 1;
	}

	const size_t nTSteps = nTIntervals + 1;

	if (nTIntervals == 0)
	{
		return 0;
	}

	double dt = (tMax - tMin) / nTIntervals;

	const Boundary* frontBoundary = model.boundaries;
	while (frontBoundary->condition->macroType != Boundary::Condition::MacroType::STEFAN_V)
	{
		++frontBoundary;
	}

	int solidRegionTag, liquidRegionTag;
	const MaterialPhase* solidMaterialsPhase, * liquidMaterialsPhase;

	if (model.regionsMaterialPhases[frontBoundary->regionsIndexes[0]]->state == MaterialPhase::SOLID)
	{
		solidRegionTag = model.regionsTags[frontBoundary->regionsIndexes[0]];
		liquidRegionTag = model.regionsTags[frontBoundary->regionsIndexes[1]];
		solidMaterialsPhase = model.regionsMaterialPhases[frontBoundary->regionsIndexes[0]];
		liquidMaterialsPhase = model.regionsMaterialPhases[frontBoundary->regionsIndexes[1]];
	}
	else
	{
		solidRegionTag = model.regionsTags[frontBoundary->regionsIndexes[1]];
		liquidRegionTag = model.regionsTags[frontBoundary->regionsIndexes[0]];
		solidMaterialsPhase = model.regionsMaterialPhases[frontBoundary->regionsIndexes[1]];
		liquidMaterialsPhase = model.regionsMaterialPhases[frontBoundary->regionsIndexes[0]];
	}

	size_t memoryBufferSize = (3 * (model.nRegions + 1) + model.nBoundaries + 1) * sizeof(size_t) +
							  2 * model.nRegions * sizeof(double*) +
							  model.nNonconformInterfaces * sizeof(unsigned int*) +
							  model.nNonconformInterfaces * sizeof(unsigned int) +
							  model.nNonconformInterfaces * sizeof(double*);
	size_t* modelMemoryBuffer = (size_t*)malloc((3 * (model.nRegions + 1) + model.nBoundaries + 1) * sizeof(size_t) + 
												 2* model.nRegions * sizeof(double*) + 
												 model.nNonconformInterfaces * sizeof(unsigned int*) +
												 model.nNonconformInterfaces * sizeof(unsigned int) +
												 model.nNonconformInterfaces * sizeof(double*));

	double* calculationBuffer = (double*)malloc(1024 * sizeof(double));
	gmsh::model::mesh::generate();
	double* initialDOFs = DG::StefanTask::computeInitialDOFs(model.nRegions, model.regionsMaterialPhases, modelMemoryBuffer);
	int viewTag = gmsh::view::add(gmshInitialModelName);
	gmsh::view::addTetrahedronsNodesData(viewTag, 0, gmshInitialModelName, initialDOFs, tMin);
	free(initialDOFs);

	//int* frontNodeTags;
	//unsigned int nFrontNodes;
	//DG::StefanTask::getFrontNodes(model.frontTag, &frontNodeTags, &nFrontNodes);
	//Coordinates* frontNodes = (Coordinates*)malloc(nFrontNodes * sizeof(Coordinates));
	//gmsh::model::getNodes(frontNodeTags, nFrontNodes, (double(*)[3])frontNodes);
	//Coordinates* frontNormals = (Coordinates*)malloc(nFrontNodes * sizeof(Coordinates));

	char* gmshPreviousModelName = (char*)malloc(sizeof(gmshInitialModelName) + (int)floor(log10(nTIntervals)) + 1);
	char* gmshModelName = (char*)malloc(sizeof(gmshInitialModelName) + (int)floor(log10(nTIntervals)) + 1);
	//char* gmshModelFileName = (char*)malloc(sizeof(gmshInitialModelFileName) + (int)floor(log10(nTIntervals)) + sizeof("_.geo_unrolled.xao") + 1);
	//char* gmshModelFileName = (char*)malloc(sizeof(gmshInitialModelFileName) + (int)floor(log10(nTIntervals)) + sizeof("_.xao") + 1);

	strcpy(gmshPreviousModelName, gmshInitialModelName);
	strcpy(gmshModelName, gmshInitialModelName);
	//strcpy(gmshModelFileName, gmshInitialModelFileName);

	//strcpy(gmshModelFileName + sizeof(gmshInitialModelFileName) - 5, "_1.xao");
	strcpy(gmshModelName + sizeof(gmshInitialModelName) - 1, "_1");

	gmsh::model::getNormals(model.frontTag, (double*)baseFrontsPoints, nFrontsNodes, (double*)frontNormals);
	DG::StefanTask::computeFrontNodesInitialMovements(*solidMaterialsPhase, *liquidMaterialsPhase, dt, latentHeat, frontNormals, nFrontsNodes, baseFrontsPoints, frontPointsMovements);

	gmsh::model::add(gmshModelName);
	gmsh::merge(gmshBaseModelFileName);
	Model::movePoints(frontsPointsTags, nFrontsNodes, baseFrontsPoints, frontPointsMovements, frontPoints);
	gmsh::merge(gmshMergeModelFileName);
	//gmsh::model::setNodes(frontNodeTags, nFrontNodes, (double(*)[3])frontNodes);
	//gmsh::write(gmshModelFileName);

	gmsh::model::mesh::generate();

	//gmsh::fltk::run();
	//gmsh::finalize();
	//return 0;

	DG::StefanTask::initSolver();

	//initialDOFs = DG::StefanTask::computeInitialDOFs(model.nRegions, model.regionsMaterialPhases, modelMemoryBuffer);
	//free(initialDOFs);

	DG::Solution solutions[2];

	DG::Solution* previousSolution = solutions;
	DG::Solution* currentSolution = solutions + 1;

	DG::StefanTask::solveInitialIteration(model.nRegions,
										  model.regionsMaterialPhases,
										  model.boundaries,
										  model.nBoundaries,
										  model.nonconformInterfaces,
										  model.nNonconformInterfaces,
										  tMin + dt,
										  dt,
										  penalty,
										  calculationBuffer,
										  modelMemoryBuffer,
										  currentSolution);

	//viewTag = gmsh::view::add(gmshModelName);
	gmsh::view::addTetrahedronsNodesData(viewTag, 1, gmshModelName, currentSolution->getDOFs(), tMin + dt);

	strcpy(gmshPreviousModelName + sizeof(gmshInitialModelName) - 1, "_0");
	for (size_t i = 2; i < nTSteps; ++i)
	{
		std::swap(previousSolution, currentSolution);
		currentSolution->clear();

		std::swap(gmshModelName, gmshPreviousModelName);
		double t = tMin + i * dt;
		gmsh::model::getNormals(model.frontTag, (double*)baseFrontsPoints, nFrontsNodes, (double*)frontNormals);
		//gmsh::model::getNormals(model.frontTag, (double*)frontNodes, nFrontNodes, (double*)frontNormals);

		DG::StefanTask::computeFrontNodesMovements(*previousSolution,
										   solidRegionTag,
										   liquidRegionTag,
										   *solidMaterialsPhase,
										   *liquidMaterialsPhase,
										   dt,
										   latentHeat,
										   frontNormals,
										   nFrontsNodes,
										   frontPoints,
										   frontPointsMovements);

		sprintf(gmshModelName + sizeof(gmshInitialModelName), "%zu", i);
		//gmsh::open(gmshModelFileName);
		gmsh::model::add(gmshModelName);

		gmsh::merge(gmshBaseModelFileName);
		//gmsh::merge(gmshInitialModelFileName);
		Model::movePoints(frontsPointsTags, nFrontsNodes, baseFrontsPoints, frontPointsMovements, frontPoints);
		gmsh::merge(gmshMergeModelFileName);

		//gmsh::model::setNodes(frontNodeTags, nFrontNodes, (double(*)[3])frontNodes);
		//gmsh::model::occ::synchronize();
		//{
		//	//gmsh::vectorpair vectpdim;
		//	//gmsh::model::occ::getEntities(vectpdim, 0);
		//	gmsh::vectorpair vectp = { {3, 1} };
		//	/*vectp.resize(nFrontNodes);
		//	for (size_t j = 0; j < nFrontNodes; ++j)
		//	{
		//		vectp[j] = { 0, frontNodeTags[j] };
		//	}*/

		//	gmsh::model::occ::translate(vectp, 0.1, 0.0, 0.0);
		//	gmsh::model::occ::synchronize();
		//}

		//sprintf(gmshModelFileName + sizeof(gmshInitialModelFileName) - 4, "%zu.xao", i);
		//gmsh::write(gmshModelFileName);

		//gmsh::model::remove();
		//gmsh::model::add(gmshModelName);

		//gmsh::merge(gmshModelFileName);
		//gmsh::open(gmshModelFileName);

		gmsh::model::mesh::generate();

		DG::StefanTask::solveIteration(model.nRegions,
									   model.regionsMaterialPhases,
									   model.boundaries,
									   model.nBoundaries,
									   model.nonconformInterfaces,
									   model.nNonconformInterfaces,
									   t,
									   dt,
									   penalty,
									   gmshModelName,
			                           *previousSolution,
									   gmshPreviousModelName,
									   calculationBuffer,
									   modelMemoryBuffer,
									   currentSolution);

		/*DG::StefanTask::solveInitialIteration(model.nRegions,
			model.regionsMaterialPhases,
			model.boundaries,
			model.nBoundaries,
			model.nonconformInterfaces,
			model.nNonconformInterfaces,
			t,
			dt,
			penalty,
			calculationBuffer,
			modelMemoryBuffer,
			currentSolution);*/

		//viewTag = gmsh::view::add(gmshModelName);
		gmsh::view::addTetrahedronsNodesData(viewTag, i, gmshModelName, currentSolution->getDOFs(), t);

		previousSolution->clear();
		std::swap(previousSolution, currentSolution);
	}

	free(frontsPointsTags);
	free(baseFrontsPoints);
	free(frontsPointsStartIndexes);
	free(frontsBuffer);
	//free(gmshModelName);
	//free(gmshModelFileName);
	free(modelMemoryBuffer);
	free(calculationBuffer);

	model.clear();
	DG::StefanTask::finalizeSolver();

	gmsh::fltk::run();
	gmsh::finalize();

	return 0;
}

const uint8_t N_CHECK_POINTS = 4;
const uint8_t N_CHECK_POINTS_DECR = N_CHECK_POINTS - 1;
LocalCoordinates3D checkLocalPoints[N_CHECK_POINTS] = { {0.0, 0.0, 0.0 },  {1.0, 0.0, 0.0 }, {0.0, 0.1, 0.0 }, {0.0, 0.0, 0.1 } };

	//for (size_t i = 0; i < nTSteps; ++i)
	//{
		//printf("Results of DG soltion #%zu\n", i);
		//printf("%lf %lf %lf %lf\n", solutionIt->compute(0, checkLocalPoints[0]), solutionIt->compute(0, checkLocalPoints[1]), solutionIt->compute(0, checkLocalPoints[2]), solutionIt->compute(0, checkLocalPoints[3]));
	//}

	//free(solutions);
	/*
	gmsh::open("./models/Stefan.geo");
	gmsh::model::mesh::generate();


	/*gmsh::vectorpair entityesDimTags;
	gmsh::model::getEntities(entityesDimTags);

	for (auto dimTag : entityesDimTags)
	{
		std::cout << "Entity dim: " << dimTag.first << "; tag: " << dimTag.second << "; is orphan: " << (gmsh::model::isEntityOrphan(dimTag.first, dimTag.second) ? true : false) << std::endl;
	}*/

	//gmsh::model::mesh::renumberElements();

#ifdef ELEMENT_DEBUG
	std::vector<int> elementTypes;
	std::vector<std::vector<size_t>> elementTags, nodesTags;
	gmsh::model::mesh::getElements(elementTypes, elementTags, nodesTags);

	int pointType = gmsh::model::mesh::getElementType("Point", 1);
	int lineType = gmsh::model::mesh::getElementType("Line", 1);
	int triangleType = gmsh::model::mesh::getElementType("Triangle", 1);
	int tetrahedronType = gmsh::model::mesh::getElementType("Tetrahedron", 1);
	size_t elementMaxTag;
	gmsh::model::mesh::getMaxElementTag(elementMaxTag);

	size_t elementsEnd = 4281;
	for (size_t i = 1; i < elementsEnd; ++i)
	{
		int elementType, dim, entityTag;
		std::vector<size_t> elementNodesTags;
		gmsh::model::mesh::getElement(i, elementType, elementNodesTags, dim, entityTag);

		std::cout << "Element tag: " << i << "; Dimension: " << dim << "; Element type:" << elementType << "; Entity tag " << entityTag << ";" << std::endl;
	}

#endif
	/*
	gmsh::finalize();

	return 0;
}
*/