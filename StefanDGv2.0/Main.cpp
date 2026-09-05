#define _CRT_SECURE_NO_WARNINGS
#define GMSH_DLL
#include "gmsh.h"
#include <iostream>
#include <cstdio>
#include "DG.h"
#include "Boundary.h"
#include "Model.h"

const char timeMeshFileName[] = "TimeMesh.txt";
const int ERROR_OF_OPENING_TIME_FILE = 1;
const int ERROR_OF_READING_TIME_PARAMETERS = 1;

int main()
{
	gmsh::initialize();
	gmsh::open("./Models/Stefan.geo");
	gmsh::model::mesh::generate();

	/*
	unsigned int nRegions = gmsh::model::regions::getCount();
	unsigned int nSurfaces = gmsh::model::surfaces::getCount();

	size_t nRegionsStartIndexes = nRegions + 1;
	size_t* memoryBuffer = (size_t*)malloc((3 * (nRegionsStartIndexes) + nSurfaces + 1) * sizeof(size_t));
	size_t* regionsStartTetrahedronsIndexes = memoryBuffer;
	gmsh::model::mesh::getRegionsTetrahedronsStartIndexes(regionsStartTetrahedronsIndexes);

	size_t* surfacesFacesStartIndexes = regionsStartTetrahedronsIndexes + nRegionsStartIndexes;
	size_t* regionsStartTetrahedronsTags = surfacesFacesStartIndexes + nSurfaces + 1;
	size_t* regionsInteriorFacesStartIndexes = regionsStartTetrahedronsTags + nRegions;

	size_t nTetFaces = regionsStartTetrahedronsIndexes[nRegions] << 2;
	size_t nTetNodes = regionsStartTetrahedronsIndexes[nRegions] << 2;
	size_t* tetFacesIndexes = (size_t*)malloc(nTetFaces * sizeof(size_t));
	size_t* tetNodesTag = (size_t*)malloc(nTetNodes * sizeof(size_t));

	gmsh::model::mesh::getTetrahedrons(surfacesFacesStartIndexes, regionsStartTetrahedronsTags, tetNodesTag, regionsInteriorFacesStartIndexes, tetFacesIndexes);

	size_t nFaces = regionsInteriorFacesStartIndexes[nRegions];

	free(memoryBuffer);
	free(tetFacesIndexes);
	free(tetNodesTag);
	gmsh::finalize();
	return 0;
	*/

	double tMin, tMax;
	size_t nTSteps;

	FILE* timeFile = fopen(timeMeshFileName, "r");

	if (timeFile == NULL)
	{
		printf("Error: unable to open file %s \n", timeMeshFileName);
		return ERROR_OF_OPENING_TIME_FILE;
	}

	if (fscanf(timeFile, "%lf %lf %zu", &tMin, &tMax, &nTSteps) != 3)
	{
		printf("Error: unable to read time parameters in file %s\n", timeMeshFileName);
		return ERROR_OF_READING_TIME_PARAMETERS;
	}

	fclose(timeFile);
	gmsh::initialize();
	gmsh::open("./Models/Stefan.geo");

	Model model;
	if (model.initilizeByCurrentGMSHModel() != Model::Error::NO_ERRORS)
	{
		gmsh::finalize();
		return 1;
	}

	gmsh::model::mesh::generate();

	DG::Solution* solutions = (DG::Solution*)malloc((nTSteps - 1)* sizeof(DG::Solution));
	
	DG::StefanTask::solve(model.nRegions,
						  model.regionsMaterialPhases,
						  model.boundaries,
						  model.nBoundaries,
						  model.nonconformInterfaces,
						  model.nNonconformInterfaces,
						  tMin,
						  tMax,
						  nTSteps,
						  solutions);

	const DG::Solution* solutionIt = solutions;
	const uint8_t N_CHECK_POINTS = 4;
	const uint8_t N_CHECK_POINTS_DECR = N_CHECK_POINTS - 1;
	LocalCoordinates3D checkLocalPoints[N_CHECK_POINTS] = { {0.0, 0.0, 0.0 },  {1.0, 0.0, 0.0 }, {0.0, 0.1, 0.0 }, {0.0, 0.0, 0.1 } };

	//for (size_t i = 0; i < nTSteps; ++i)
	//{
		//printf("Results of DG soltion #%zu\n", i);
		//printf("%lf %lf %lf %lf\n", solutionIt->compute(0, checkLocalPoints[0]), solutionIt->compute(0, checkLocalPoints[1]), solutionIt->compute(0, checkLocalPoints[2]), solutionIt->compute(0, checkLocalPoints[3]));
	//}

	free(solutions);
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
	gmsh::finalize();

	return 0;
}