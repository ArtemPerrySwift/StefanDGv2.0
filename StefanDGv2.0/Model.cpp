#define GMSH_DLL
#include "Model.h"
#include "gmsh.h"
#include <climits>
#include <utility>
#include <unordered_map>
#include <stdexcept>

typedef gmsh::model::regions::RegionApiResult GMSHRegionApiResult;
typedef gmsh::model::surfaces::SurfaceApiResult GMSHSurfaceApiResult;
typedef gmsh::model::physical_groups::PGApiResult GMSHPGApiResult;

void addToMap(std::pair<int, const Boundary::Condition*>* const conditionByPGTagBegin,
              const std::pair<int, const Boundary::Condition*>* const conditionByPGTagEnd,
              const unsigned int nPGs,
              const int tag,
              const Boundary::Condition* condition)
{
    std::pair<int, const Boundary::Condition*>* conditionByPGTagIt = conditionByPGTagBegin + tag % nPGs;
    while (conditionByPGTagIt->first != 0 && conditionByPGTagIt != conditionByPGTagEnd)
    {
        ++conditionByPGTagIt;
    }

    if (conditionByPGTagIt->first != 0)
    {
        conditionByPGTagIt = conditionByPGTagBegin;
        while (conditionByPGTagIt->first != 0)
        {
            ++conditionByPGTagIt;
        }
    }

    conditionByPGTagIt->first = tag;
    conditionByPGTagIt->second = condition;
}

unsigned int getConditions(const int* PGsTags,
				   const unsigned int nPGs,
				   const char* names,
				   const unsigned int* namesStartIndexes,
				   Boundary::ValueCondition* valueConditionMemoryIt,
				   std::pair<int, const Boundary::Condition*>* conditionByPGTag,
                   unsigned int &nNoncofnormInterfaces)
{

    const std::pair<int, const Boundary::Condition*>* const conditionByPGTagEndIt = conditionByPGTag + nPGs;

    unsigned int iNonconformInterface = 0;
    unsigned int nValuesConditions = 0;

	for (unsigned int i = 0; i < nPGs; ++i)
	{
		const char* name = names + *namesStartIndexes;
		switch (*name)
		{
        case 'D':
        {
            ++name;
            if (*name == ':')
            {
                ++name;

                char* endptr;
                double value = strtod(name, &endptr);
                if (endptr == name)
                {
                    throw std::runtime_error("Dirichlet value is not set in physical name");
                }
                new (valueConditionMemoryIt) Boundary::ValueCondition(Boundary::ValueCondition::Type::DIRICHLET, value);
                addToMap(conditionByPGTag, conditionByPGTagEndIt, nPGs, *PGsTags, valueConditionMemoryIt);

                ++valueConditionMemoryIt;
                ++nValuesConditions;
            }
            else
            {
                addToMap(conditionByPGTag, conditionByPGTagEndIt, nPGs, *PGsTags, &Boundary::DirichletFunctionCondition);
            }
            break;
        }
        case 'N':
        {
            ++name;
            if (*name == ':')
            {
                ++name;
                char* endptr;
                double value = strtod(name, &endptr);
                if (endptr == name)
                {
                    throw std::runtime_error("Newman value is not set in physical name");
                }
                new (valueConditionMemoryIt) Boundary::ValueCondition(Boundary::ValueCondition::Type::NEWMAN, value);
                addToMap(conditionByPGTag, conditionByPGTagEndIt, nPGs, *PGsTags, valueConditionMemoryIt);

                ++valueConditionMemoryIt;
                ++nValuesConditions;
            }
            else
            {
                addToMap(conditionByPGTag, conditionByPGTagEndIt, nPGs, *PGsTags, &Boundary::NewmanFunctionCondition);
            }
            break;
        }
        case 'C':
        {
            ++iNonconformInterface;
            addToMap(conditionByPGTag, conditionByPGTagEndIt, nPGs, *PGsTags, &Boundary::NonconformInterfaceCondition);
            break;
        }
        case 'S':
        {
            ++name;
            if (*name == ':')
            {
                ++name;
                char* endptr;
                double value = strtod(name, &endptr);
                if (endptr == name)
                {
                    throw std::runtime_error("Stefan value is not set in physical name");
                }
                new (valueConditionMemoryIt) Boundary::ValueCondition(Boundary::ValueCondition::Type::STEFAN, value);
                addToMap(conditionByPGTag, conditionByPGTagEndIt, nPGs, *PGsTags, valueConditionMemoryIt);

                ++valueConditionMemoryIt;
                ++nValuesConditions;
            }
            else
            {
                throw std::runtime_error("Stefan name should contain \':\' after \'S\'");
            }
            break;
        }

        default:
        {
            throw std::runtime_error("Unknown type of condition");
        }
		}

        ++namesStartIndexes;
        ++PGsTags;
	}

    nNoncofnormInterfaces = iNonconformInterface;

    return nValuesConditions;
}

enum READ_CONDITIONS_ERRORS : uint8_t { VALUE_ERR = 1, AFTER_VALUE_ERR };

uint8_t readValueCondition(Boundary::ValueCondition::Type conditionType, Boundary::ValueCondition* valueConditionMemoryIt, const char* strBegin, const char* strEnd)
{
    char* endptr;
    double value = strtod(strBegin, &endptr);
    if (endptr == strBegin)
    {
        return VALUE_ERR;
    }

    while (*endptr == ' ' || *endptr == ';')
    {
        ++endptr;
    }

    if (endptr != strEnd)
    {
        return AFTER_VALUE_ERR;
    }
    new (valueConditionMemoryIt) Boundary::ValueCondition(conditionType, value);
    return 0;
}

bool getConditions(const unsigned int *nameStartIndexIt,
                   const char* namesStart,
                   const unsigned int nPGs,
                   const Boundary::Condition** conditionPtrIt,
                   Boundary::ValueCondition* valueConditionMemoryIt,
                   unsigned int& nValueConditions,
                   unsigned int& nNonconformInterfaces)
{
    unsigned int iValueCondition = 0, iNonconformInterafce = 0;
    bool isNameError = false, isNamesErrors = false;
    const char* name, *nameSymbolIt, *nameEnd;

    nameEnd = namesStart + *nameStartIndexIt;

    for (unsigned int i = 0; i < nPGs; ++i)
    {
        name = nameEnd;
        nameSymbolIt = name;

        ++nameStartIndexIt;
        nameEnd = namesStart + *nameStartIndexIt;

        switch (*nameSymbolIt)
        {
        case 'D':
        {
            ++nameSymbolIt;
            if (*nameSymbolIt == ':')
            {
                ++nameSymbolIt;
                uint8_t err = readValueCondition(Boundary::ValueCondition::Type::DIRICHLET, valueConditionMemoryIt, nameSymbolIt, nameEnd);

                if (err != 0)
                {
                    isNameError = true;
                    break;
                }

                *conditionPtrIt = valueConditionMemoryIt;

                ++valueConditionMemoryIt;
                ++iValueCondition;
            }
            else
            {
                while (*nameSymbolIt == ' ' || *nameSymbolIt == ';')
                {
                    ++nameSymbolIt;
                }

                if (nameSymbolIt != nameEnd)
                {
                    isNameError = true;
                    break;
                }

                *conditionPtrIt = &Boundary::DirichletFunctionCondition;

            }
            break;
        }
        case 'N':
        {
            ++nameSymbolIt;
            if (*nameSymbolIt == ':')
            {
                ++nameSymbolIt;
                uint8_t err = readValueCondition(Boundary::ValueCondition::Type::NEWMAN, valueConditionMemoryIt, nameSymbolIt, nameEnd);

                if (err != 0)
                {
                    isNameError = true;
                    break;
                }
                *conditionPtrIt = valueConditionMemoryIt;

                ++valueConditionMemoryIt;
                ++iValueCondition;
            }
            else
            {
                while (*nameSymbolIt == ' ' || *nameSymbolIt == ';')
                {
                    ++nameSymbolIt;
                }

                if (nameSymbolIt != nameEnd)
                {
                    isNameError = true;
                    break;
                }

                *conditionPtrIt = &Boundary::NewmanFunctionCondition;
            }
            break;
        }
        case 'C':
        {
            ++nameSymbolIt;
            if (*nameSymbolIt != '_')
            {
                isNameError = true;
                break;
            }

            ++iNonconformInterafce;
            *conditionPtrIt = &Boundary::NonconformInterfaceCondition;
            break;
        }
        case 'S':
        {
            ++nameSymbolIt;
            if (*nameSymbolIt == ':')
            {
                ++nameSymbolIt;
                uint8_t err = readValueCondition(Boundary::ValueCondition::Type::STEFAN, valueConditionMemoryIt, nameSymbolIt, nameEnd);

                if (err != 0)
                {
                    isNameError = true;
                    break;
                }

                *conditionPtrIt = valueConditionMemoryIt;

                ++valueConditionMemoryIt;
                ++iValueCondition;
                break;
            }
            else
            {
                isNameError = true;
            }
            break;
        }

        default:
        {
            isNameError = true;
        }
        }

        if (isNameError)
        {
            conditionPtrIt = nullptr;
            printf("Error: Invalid condition name ");
            nameSymbolIt = name;
            while (nameSymbolIt != nameEnd)
            {
                printf("%c", *nameSymbolIt);
                ++nameSymbolIt;
            }
            printf("\n");
            isNamesErrors = true;

            isNameError = false;
        }
        ++conditionPtrIt;
    }

    nValueConditions = iValueCondition;
    nNonconformInterfaces = iNonconformInterafce;

    return isNamesErrors;
}

bool getMaterialPhases(const unsigned int* nameStartIndexIt,
                       const char* namesStart,
                       const unsigned int nPGs,
                       MaterialPhase* materialPhaseIt)
{
    bool errFl = false;
    const char* nameEnd = namesStart + *nameStartIndexIt;
    const char* name, *nameSymbolIt;
    for (unsigned int i = 0; i < nPGs; ++i)
    {
        name = nameEnd;

        ++nameStartIndexIt;
        nameEnd = namesStart + *nameStartIndexIt;

        int nReadSymbols;
        char state;
#pragma warning(suppress : 4996)
        uint8_t nReadValues = sscanf(name, "M:%c;%lf;%lf;%lf;%n", &state, &materialPhaseIt->thermalConductivity, &materialPhaseIt->heatCapacity, &materialPhaseIt->density, &nReadSymbols);

        nameSymbolIt = name + nReadSymbols;
        if (nReadValues != 4 || nameSymbolIt != nameEnd)
        {
            printf("Error: Invalid material name %s", name);

            nameSymbolIt = name;
            while (nameSymbolIt != nameEnd)
            {
                printf("%c", *nameSymbolIt);
                ++nameSymbolIt;
            }
            printf("\n");
            errFl = true;
        }

        if (state == 'L')
        {
            materialPhaseIt->state = MaterialPhase::LIQUID;
        }
        else if (state == 'S')
        {
            materialPhaseIt->state = MaterialPhase::SOLID;
        }
        else
        {
            printf("Error: Invalid material name %s", name);
            errFl = true;
        }


        ++materialPhaseIt;
    }

    return errFl;
}

const unsigned char PLANE_SURFACE = 16;

unsigned int bsearch(const int key, const int arr[], const unsigned int n)
{
    unsigned int l = 0, r = n, i;
    while (arr[l] != key)
    {
        i = (l + r) >> 1;
        if (arr[i] > key)
        {
            r = i;
        }
        else
        {
            l = i;
        }
    }

    return l;
}

bool initilizeBoundaries(const int orderedSurfacesPGsTags[],
                          const unsigned int nSurfacesPGs,
                          const Boundary::Condition* conditions[],
                          unsigned int* surfacesPGsStartIndexIt,
                          int* surfacesPGTagIt,
                          const unsigned int nSurfaces,
                          Boundary* boundaryIt)
{
    
    for (unsigned int i = 0; i < nSurfaces; ++i)
    {
        unsigned int nSurfacePGs = *(surfacesPGsStartIndexIt + 1) - *surfacesPGsStartIndexIt;
        if (nSurfacePGs > 1)
        {
            printf("Error: boundary have several conditions");
            return true;
        }

        if (nSurfacePGs == 0)
        {
            boundaryIt->condition = nullptr;
        }
        else
        {
            unsigned int PGIndex = bsearch(*surfacesPGTagIt, orderedSurfacesPGsTags, nSurfacesPGs);
            boundaryIt->condition = conditions[PGIndex];

            ++surfacesPGTagIt;
        }

        boundaryIt->regionsIndexes[0] = UINT_MAX;
        boundaryIt->regionsIndexes[1] = UINT_MAX;

        ++surfacesPGsStartIndexIt;
        ++boundaryIt;
    }

    return false;
}

bool processRegions(const unsigned int nRegions,
                    const unsigned int* regionsBoundaryStartIndexIt,
                    const int* regionsBoundaryTagIt,
                    const unsigned int* regionsPGStartIndexIt,
                    const int* regionsPGTagIt,
                    const int orderedBoundariesTags[],
                    const unsigned int nBoundaries,
                    const int orderedUniqueRegionsPGsTags[],
                    const MaterialPhase materialPhases[],
                    const unsigned int nUniqueRegionsPGs,
                    Boundary boundaries[],
                    const MaterialPhase** regionsMaterialPhases)
{
    unsigned int iBoundary = *regionsBoundaryStartIndexIt;

    for (unsigned int i = 0; i < nRegions; ++i)
    {
        unsigned int nRegionPGs = *(regionsPGStartIndexIt + 1) - *regionsPGStartIndexIt;
        if (nRegionPGs != 1)
        {
            printf("Error: region should have only 1 material phase");
            return true;
        }

        unsigned int tagIndex = bsearch(*regionsPGTagIt, orderedUniqueRegionsPGsTags, nUniqueRegionsPGs);
        *regionsMaterialPhases = materialPhases + tagIndex;

        ++regionsBoundaryStartIndexIt;
        unsigned int regionBoundariesEndIndex = *regionsBoundaryStartIndexIt;
        while (iBoundary < regionBoundariesEndIndex)
        {
            tagIndex = bsearch(*regionsBoundaryTagIt, orderedBoundariesTags, nBoundaries);
            Boundary* boundary = boundaries + tagIndex;
            if (*boundary->regionsIndexes == UINT_MAX)
            {
                *boundary->regionsIndexes = i;
            }
            else
            {
                *(boundary->regionsIndexes + 1) = i;
            }

            ++regionsBoundaryTagIt;
            ++iBoundary;
        }

        ++regionsMaterialPhases;
        ++regionsPGStartIndexIt;
        ++regionsPGTagIt;
    }

    return false;
}

void extractNonconformPGsTags(const int *surfacePGTagIt,
    const unsigned int nSurfacesPGs,
    const Boundary::Condition** conditionIt,
    int* nonconformConditionPGTagIt,
    bool* nonconformConditionPGTagIsMetIt)
{
    for (unsigned int i = 0; i < nSurfacesPGs; ++i)
    {
        if ((*conditionIt)->macroType == Boundary::Condition::MacroType::NONCONFORM_INTERFACE)
        {
            *nonconformConditionPGTagIt = *surfacePGTagIt;
            *nonconformConditionPGTagIsMetIt = false;

            ++nonconformConditionPGTagIt;
            ++nonconformConditionPGTagIsMetIt;
        }

        ++conditionIt;
        ++surfacePGTagIt;
    }
}

unsigned int processBoundariesFinal(const int surfacesTags[],
                                    const unsigned int nSurfaces,
                                    const int surfacesUniquePGsTags[],   
                                    const unsigned int nSurfacesPGs,
                                    const Boundary::Condition** conditionIt, 
                                    const int surfacesPGsTags[],
                                    const unsigned int surfacesPGsTagsStartTags[],
                                    const unsigned int nNonconformInterfaces,
                                    const MaterialPhase* const regionsMaterialPhases[],
                                    Boundary* boundaryIt,
                                    Boundary::ConformCondition* conformConditionMemoryIt,
                                    NonconformInterface nonconfromInterfaces[])
{
    void* memoryBuffer = malloc(nNonconformInterfaces * (sizeof(int) + sizeof(bool)));
    int* nonconformConditionPGTags = (int*)memoryBuffer;
    bool* nonconformConditionPGTagIsMet = (bool*)(nonconformConditionPGTags + nNonconformInterfaces);

    unsigned int iConformCondition = 0;

    extractNonconformPGsTags(surfacesUniquePGsTags, nSurfacesPGs, conditionIt, nonconformConditionPGTags, nonconformConditionPGTagIsMet);

    for (unsigned int i = 0; i < nSurfaces; ++i)
    {
        if (boundaryIt->condition == NULL)
        {
            if (boundaryIt->regionsIndexes[1])
            {
                boundaryIt->condition = &Boundary::HomogeneousNemanCondition;
            }
            else
            {
                if (regionsMaterialPhases[boundaryIt->regionsIndexes[0]] == regionsMaterialPhases[boundaryIt->regionsIndexes[1]])
                {
                    new (conformConditionMemoryIt) Boundary::ConformCondition(regionsMaterialPhases[boundaryIt->regionsIndexes[0]]->thermalConductivity);
                }
                else
                {
                    new (conformConditionMemoryIt) Boundary::ConformCondition((regionsMaterialPhases[boundaryIt->regionsIndexes[0]]->thermalConductivity +
                        regionsMaterialPhases[boundaryIt->regionsIndexes[1]]->thermalConductivity) / 2.0);
                }

                boundaryIt->condition = conformConditionMemoryIt;
                ++conformConditionMemoryIt;
                ++iConformCondition;
            }
        }
        else if (boundaryIt->condition->macroType == Boundary::Condition::MacroType::NONCONFORM_INTERFACE)
        {
            unsigned int tagIndex = bsearch(surfacesPGsTags[surfacesPGsTagsStartTags[i]], nonconformConditionPGTags, nNonconformInterfaces);
            NonconformInterface* nonconformInterface = nonconfromInterfaces + tagIndex;

            if (nonconformConditionPGTagIsMet[tagIndex])
            {
                nonconformInterface->sidesTags[1] = surfacesTags[i];
                nonconformInterface->sideIndexes[1] = i;
                nonconformInterface->regionsIndexes[1] = boundaryIt->regionsIndexes[0];

                if (regionsMaterialPhases[nonconformInterface->regionsIndexes[0]] == regionsMaterialPhases[nonconformInterface->regionsIndexes[1]])
                {
                    nonconformInterface->thermalConductivity = regionsMaterialPhases[nonconformInterface->regionsIndexes[0]]->thermalConductivity;
                }
                else
                {
                    nonconformInterface->thermalConductivity = (regionsMaterialPhases[nonconformInterface->regionsIndexes[0]]->thermalConductivity +
                        regionsMaterialPhases[nonconformInterface->regionsIndexes[1]]->thermalConductivity) / 2.0;
                }
            }
            else
            {
               *(nonconformInterface->sidesTags) = surfacesTags[i];
               *(nonconformInterface->sideIndexes) = i;
               *(nonconformInterface->regionsIndexes) = boundaryIt->regionsIndexes[0];

               nonconformConditionPGTagIsMet[tagIndex] = true;
            }
        }

        ++boundaryIt;
    }

    free(memoryBuffer);

    return iConformCondition;
}

bool Model::initilizeByCurrentGMSHModel()
{
    unsigned int nRegions = gmsh::model::regions::getCount();
    unsigned int nSurfaces = gmsh::model::surfaces::getCount();

    this->nRegions = nRegions;
    this->nBoundaries = nSurfaces;

    Boundary* boundaries = (Boundary*)malloc(nSurfaces * sizeof(Boundary));
    const MaterialPhase** regionsMaterialPhases = (const MaterialPhase**)malloc(nRegions * sizeof(MaterialPhase*));

    this->boundaries = boundaries;
    this->regionsMaterialPhases = regionsMaterialPhases;

    void* regionsMemoryPull = malloc(((nRegions + 1) << 1 + nRegions) * sizeof(int));
    int* regionsTags = (int*)regionsMemoryPull;
    unsigned int* regionsBoundariesStartIndexes = (unsigned int*)(regionsTags + nRegions);
    unsigned int* regionsPGsStartIndexes = regionsBoundariesStartIndexes + nRegions + 1;

    gmsh::model::regions::getStartIndexes(regionsBoundariesStartIndexes, regionsPGsStartIndexes);

    void* regionsAdditionMemoryPull = malloc((regionsBoundariesStartIndexes[nRegions] + regionsPGsStartIndexes[nRegions]) * sizeof(int));
    int* regionsBoundariesTags = (int *)regionsAdditionMemoryPull;
    int* regionsPGsTags = regionsBoundariesTags + regionsBoundariesStartIndexes[nRegions];

    gmsh::model::regions::getRegions(regionsTags, regionsBoundariesTags, regionsPGsTags);


    void* surfacesMemoryPull = malloc(nSurfaces * (sizeof(int) + sizeof(char)) + (nSurfaces + 1) * sizeof(unsigned int));
    int* surfacesTags = (int*)surfacesMemoryPull;
    unsigned char* surfacesTypes = (unsigned char*)(surfacesTags + nSurfaces);
    unsigned int* surfacesPGsStartIndexes = (unsigned int*)(surfacesTypes + nSurfaces);

    gmsh::model::surfaces::getStartIndexes(surfacesPGsStartIndexes);
    int* surfacesPGs = (int*)malloc(surfacesPGsStartIndexes[nSurfaces] * sizeof(int));

    gmsh::model::surfaces::getSurfaces(surfacesTags, surfacesTypes, surfacesPGs);

    int nPGs = gmsh::model::physical_groups::getCount();
    void* PGsMemoryPull = malloc((nPGs << 2 + 1) * sizeof(int));

    int* PGsTags = (int*)PGsMemoryPull;
    unsigned int* PGsNamesStartIndexes = (unsigned int*)(PGsTags + nPGs);
    unsigned int dimensionsPGsStartIndexes[5];

    gmsh::model::physical_groups::getStartIndexes(dimensionsPGsStartIndexes, PGsNamesStartIndexes);
    
    char* PGsNames = (char*)malloc(PGsNamesStartIndexes[nPGs] * sizeof(char));
    gmsh::model::physical_groups::getPhysicalGroups(PGsTags, PGsNames);

    unsigned int nSurfacesPGs = dimensionsPGsStartIndexes[3] - dimensionsPGsStartIndexes[2];
    const Boundary::Condition** conditions = (const Boundary::Condition**)malloc(nSurfacesPGs * sizeof(Boundary::Condition*));
    Boundary::ValueCondition* valueConditions = (Boundary::ValueCondition*)malloc(nSurfacesPGs * sizeof(Boundary::ValueCondition));

    unsigned int nValueConditions, nNonconformInterfaces;
    bool err = getConditions(PGsNamesStartIndexes + dimensionsPGsStartIndexes[2], PGsNames, nSurfacesPGs, conditions, valueConditions, nValueConditions, nNonconformInterfaces);

    if (err)
    {
        free(boundaries);
        free(regionsMaterialPhases);
        free(regionsMemoryPull);
        free(regionsAdditionMemoryPull);
        free(surfacesMemoryPull);
        free(surfacesPGs);
        free(PGsMemoryPull);
        free(PGsNames);
        free(valueConditions);
        free(conditions);

        return true;
    }

    realloc(valueConditions, nValueConditions * sizeof(Boundary::ValueCondition));
    this->valueConditions = valueConditions;

    err = initilizeBoundaries(PGsTags + dimensionsPGsStartIndexes[2], nSurfacesPGs, conditions, surfacesPGsStartIndexes, surfacesPGs, nSurfaces, boundaries);

    if (err)
    {
        free(boundaries);
        free(regionsMaterialPhases);
        free(regionsMemoryPull);
        free(regionsAdditionMemoryPull);
        free(surfacesMemoryPull);
        free(surfacesPGs);
        free(PGsMemoryPull);
        free(PGsNames);
        free(valueConditions);
        free(conditions);

        return true;
    }

    unsigned int nRegionsPGs = dimensionsPGsStartIndexes[4] - dimensionsPGsStartIndexes[3];
    MaterialPhase* materialPhases = (MaterialPhase*)malloc(nRegionsPGs * sizeof(MaterialPhase));

    err = getMaterialPhases(PGsNamesStartIndexes + dimensionsPGsStartIndexes[3], PGsNames, nRegionsPGs, materialPhases);

    if (err)
    {
        free(boundaries);
        free(regionsMaterialPhases);
        free(regionsMemoryPull);
        free(regionsAdditionMemoryPull);
        free(surfacesMemoryPull);
        free(surfacesPGs);
        free(PGsMemoryPull);
        free(PGsNames);
        free(valueConditions);
        free(conditions);
        free(materialPhases);

        return true;
    }

    this->materialPhases = materialPhases;
    free(PGsNames);

    err = processRegions(nRegions,
                         regionsBoundariesStartIndexes,
                         regionsBoundariesTags,
                         regionsPGsStartIndexes,
                         regionsPGsTags,
                         surfacesTags,
                         nSurfaces,
                         PGsTags + dimensionsPGsStartIndexes[3],
                         materialPhases,
                         nRegionsPGs,
                         boundaries,
                         regionsMaterialPhases);


    if (err)
    {
        free(boundaries);
        free(regionsMaterialPhases);
        free(regionsMemoryPull);
        free(regionsAdditionMemoryPull);
        free(surfacesMemoryPull);
        free(surfacesPGs);
        free(PGsMemoryPull);
        free(valueConditions);
        free(conditions);
        free(materialPhases);

        return true;
    }

    this->regionsMaterialPhases = regionsMaterialPhases;

    NonconformInterface* noncofnormInterfaces = (NonconformInterface*)malloc(nNonconformInterfaces * sizeof(NonconformInterface));
    Boundary::ConformCondition* conformConditions = (Boundary::ConformCondition*)malloc((nSurfaces - nSurfacesPGs) * sizeof(Boundary::ConformCondition));

    unsigned int nConformConditions = processBoundariesFinal(surfacesTags,
                                                             nSurfaces,
                                                             PGsTags + dimensionsPGsStartIndexes[2],
                                                             nSurfacesPGs,
                                                             conditions,
                                                             surfacesPGs,
                                                             surfacesPGsStartIndexes,
                                                             nNonconformInterfaces,
                                                             regionsMaterialPhases,
                                                             boundaries,
                                                             conformConditions,
                                                             noncofnormInterfaces);

    realloc(conformConditions, nConformConditions * sizeof(Boundary::ConformCondition));

    this->nonconformInterfaces = noncofnormInterfaces;
    this->conformConditions = conformConditions;

    free(regionsMemoryPull);
    free(regionsAdditionMemoryPull);
    free(surfacesMemoryPull);
    free(surfacesPGs);
    free(PGsMemoryPull);
    free(conditions);

    return false;
}
