#define _CRT_SECURE_NO_WARNINGS
#define GMSH_DLL
#include "gmsh.h"
#include <iostream>
#include <cstdio>
#include "DG.h"
#include "Boundary.h"
#include "Model.h"

const char timeMeshFileName[] = "TimeMesh.txt";
const char penaltyFileName[] = "Penalty.txt";
const char latentHeatFileName[] = "LatentHeat.txt";
const char gmshInitialModelName[] = "Stefan";
const char gmshInitialModelFileName[] = "./Models/Stefan.geo";

enum Error {
	OPENING_TIME_FILE = 1,
	READING_TIME_PARAMETERS,
	OPENING_PENALTY_FILE,
	READING_PENALTY_PARAMETER,
	OPENING_LATENT_HEAT_FILE,
	READING_LATENT_HEAT_PARAMETER,
	READING_MODEL
};

static uint8_t readTimeParameters(double* tMin, double* tMax, size_t* nTSteps)
{
	FILE* timeFile = fopen(timeMeshFileName, "r");

	if (timeFile == NULL)
	{
		printf("Error: unable to open file %s \n", timeMeshFileName);
		return Error::OPENING_TIME_FILE;
	}

	if (fscanf(timeFile, "%lf %lf %zu", tMin, tMax, nTSteps) != 3)
	{
		printf("Error: unable to read time parameters in file %s\n", timeMeshFileName);
		return Error::READING_TIME_PARAMETERS;
	}

	fclose(timeFile);
	return 0;
}

static uint8_t readPenalty(double* penalty)
{
	FILE* penaltyFile = fopen(penaltyFileName, "r");

	if (penaltyFile == NULL)
	{
		printf("Error: unable to open file %s \n", penaltyFileName);
		return Error::OPENING_PENALTY_FILE;
	}

	if (fscanf(penaltyFile, "%lf", penalty) != 1)
	{
		printf("Error: unable to read penalty parameter in file %s\n", penaltyFileName);
		return Error::READING_PENALTY_PARAMETER;
	}

	fclose(penaltyFile);
	return 0;
}

static uint8_t readLatentHeat(double* latentHeat)
{
	FILE* latentHeatFile = fopen(penaltyFileName, "r");

	if (latentHeatFile == NULL)
	{
		printf("Error: unable to open file %s \n", penaltyFileName);
		return Error::OPENING_LATENT_HEAT_FILE;
	}

	if (fscanf(latentHeatFile, "%lf", latentHeat) != 1)
	{
		printf("Error: unable to latent heat parameter in file %s\n", penaltyFileName);
		return Error::READING_LATENT_HEAT_PARAMETER;
	}

	fclose(latentHeatFile);
	return 0;
}



int main()
{
	gmsh::initialize();

	double tMin, tMax;
	size_t nTSteps;

	uint8_t err = readTimeParameters(&tMin, &tMax, &nTSteps);
	if (err)
	{
		return err;
	}

	double penalty;
	err = readPenalty(&penalty);
	if (err)
	{
		return err;
	}

	double latentHeat;
	err = readLatentHeat(&latentHeat);
	if (err)
	{
		return err;
	}

	gmsh::initialize();
	gmsh::open(gmshInitialModelFileName);

	gmsh::model::mesh::generate();

	DG::Solution solutions[2];

	DG::Solution* previousSolution = solutions;
	DG::Solution* currentSolution = solutions + 1;

	gmsh::model::mesh::generate();

	Model model;
	if (model.initilizeByCurrentGMSHModel() != Model::Error::NO_ERRORS)
	{
		gmsh::finalize();
		return 1;
	}

	const size_t nTIntervals = nTSteps - 1;

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

	int* frontNodeTags;
	unsigned int nFrontNodes;
	DG::StefanTask::getFrontNodes(model.frontTag, &frontNodeTags, &nFrontNodes);
	Coordinates* frontNodes = (Coordinates*)malloc(nFrontNodes * sizeof(Coordinates));
	gmsh::model::getNodes(frontNodeTags, nFrontNodes, (double(*)[3])frontNodes);
	Coordinates* frontNormals = (Coordinates*)malloc(nFrontNodes * sizeof(Coordinates));
	gmsh::model::getNormals(model.frontTag, (double*)frontNodes, nFrontNodes, (double*)frontNormals);

	DG::StefanTask::relocateInitialFrontNodes(*solidMaterialsPhase, *liquidMaterialsPhase,dt, latentHeat, frontNormals, nFrontNodes, frontNodes);

	char* gmshModelName = (char*)malloc(sizeof(gmshInitialModelName) + (int)floor(log10(nTIntervals)));
	char* gmshModelFileName = (char*)malloc(sizeof(gmshInitialModelFileName) + (int)floor(log10(nTIntervals)));

	strcpy(gmshModelName, gmshInitialModelName);
	strcpy(gmshModelFileName, gmshInitialModelFileName);

	strcpy(gmshModelFileName + sizeof(gmshInitialModelFileName) - 1, "_1.geo");
	strcpy(gmshModelName + sizeof(gmshInitialModelName) - 1, "_1");

	gmsh::model::add(gmshModelName);
	gmsh::merge(gmshInitialModelFileName);
	gmsh::model::setNodes(frontNodeTags, nFrontNodes, (double(*)[3])frontNodes);
	gmsh::write(gmshModelFileName);

	gmsh::model::mesh::generate();

	size_t* modelMemoryBuffer = (size_t*)malloc((3 * (model.nRegions + 1) + model.nBoundaries + 1) * sizeof(size_t));
	double* calculationBuffer = (double*)malloc(1024 * sizeof(double));

	DG::StefanTask::initSolver();

	DG::StefanTask::solveInitialIteration(model.nRegions,
										  model.regionsMaterialPhases,
										  model.boundaries,
										  model.nBoundaries,
										  model.nonconformInterfaces,
										  model.nNonconformInterfaces,
										  dt,
										  penalty,
										  calculationBuffer,
										  modelMemoryBuffer,
										  currentSolution);



	const uint8_t DEBUG_DOFS_SIZE = 4;

	for (size_t i = 2; i < nTSteps; ++i)
	{
		gmsh::model::getNormals(model.frontTag, (double*)frontNodes, nFrontNodes, (double*)frontNormals);

		DG::StefanTask::relocateFrontNodes(*previousSolution,
										   solidRegionTag,
										   liquidRegionTag,
										   *solidMaterialsPhase,
										   *liquidMaterialsPhase,
										   dt,
										   latentHeat,
										   frontNormals,
										   nFrontNodes,
										   frontNodes);

		sprintf(gmshModelName + sizeof(gmshInitialModelName), "%zu.geo", i);
		gmsh::model::add(gmshModelName);

		gmsh::merge(gmshModelFileName);

		gmsh::model::setNodes(frontNodeTags, nFrontNodes, (double(*)[3])frontNodes);
		sprintf(gmshModelFileName + sizeof(gmshInitialModelFileName), "%zu.geo", i);
		gmsh::write(gmshModelFileName);

		gmsh::model::mesh::generate();

		DG::StefanTask::solveIteration(model.nRegions,
									   model.regionsMaterialPhases,
									   model.boundaries,
									   model.nBoundaries,
									   model.nonconformInterfaces,
									   model.nNonconformInterfaces,
									   dt,
									   penalty,
			                           *previousSolution,
									   calculationBuffer,
									   modelMemoryBuffer,
									   currentSolution);

		int viewTag = gmsh::view::add(gmshModelName);
		gmsh::view::addTetrahedronsNodesData(viewTag, i, gmshModelName, currentSolution->getDOFs(), tMin + i * dt);

		previousSolution->clear();
		std::swap(previousSolution, currentSolution);
	}

	free(gmshModelName);
	free(gmshModelFileName);
	free(modelMemoryBuffer);
	free(calculationBuffer);

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