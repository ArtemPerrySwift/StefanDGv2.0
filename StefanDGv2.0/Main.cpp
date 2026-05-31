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

	/*
	int* surfacesTags, * surfacesPGsTags;
	unsigned char* surfaceTypes;
	unsigned int nSurfaces, * surfacesNPGs;
	gmsh::model::getSurfaces(surfacesTags,
		nSurfaces,
		surfaceTypes,
		surfacesPGsTags,
		surfacesNPGs);
		*/

	/*
	int* regionsTags, * regionsBoundariesTags, * regionsPGsTags;
	unsigned int nRegions, * regionsNBoundariesTags, * regionsNPGs;
	gmsh::model::getRegions(regionsTags, nRegions, regionsBoundariesTags,
		regionsNBoundariesTags, regionsPGsTags, regionsNPGs);


	int* surfacesTags, * surfacesPGsTags;
	unsigned int nSurfaces, *surfacesNPGs;
	gmsh::model::getSurfaces(surfacesTags, nSurfaces, surfacesPGsTags, surfacesNPGs);

	int* PGsTags;
	unsigned int dimensionsNTags[4];
	char* names;
	unsigned int* namesSizes;
	gmsh::model::getPhysicalGroups(PGsTags, dimensionsNTags, names, namesSizes);
	*/
	Model model;
	model.initilizeByCurrentGMSHModel();

	DG::Solution* solutions = DG::solveStefanTask(tMin, tMax, nTSteps);

	const DG::Solution* solutionIt = solutions;
	const uint8_t N_CHECK_POINTS = 4;
	const uint8_t N_CHECK_POINTS_DECR = N_CHECK_POINTS - 1;
	LocalCoordinates3D checkLocalPoints[N_CHECK_POINTS] = { {0.0, 0.0, 0.0 },  {1.0, 0.0, 0.0 }, {0.0, 0.1, 0.0 }, {0.0, 0.0, 0.1 } };

	for (size_t i = 0; i < nTSteps; ++i)
	{
		printf("Results of DG soltion #%zu\n", i);
		printf("%lf %lf %lf %lf\n", solutionIt->compute(0, checkLocalPoints[0]), solutionIt->compute(0, checkLocalPoints[1]), solutionIt->compute(0, checkLocalPoints[2]), solutionIt->compute(0, checkLocalPoints[3]));
	}
	delete[] solutions;
	
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