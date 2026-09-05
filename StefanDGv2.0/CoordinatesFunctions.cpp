#include "CoordinatesFunctions.h"
#include <cmath>

namespace CoordinatesFunctions
{
    void computeDiffrence(const Coordinates& point1, const Coordinates& point2, Coordinates& vector)
    {
        vector.x = point1.x - point2.x;
        vector.y = point1.y - point2.y;
        vector.z = point1.z - point2.z;
    }

    void computeTransitionMatrix2DExcludeX(const Coordinates nodes[], const size_t* triangleNodeIndexIt, double transitionMatrix[4])
    {
        ++triangleNodeIndexIt;
        const Coordinates* tetrahedroBeginNodeIt = nodes + *triangleNodeIndexIt;
        ++triangleNodeIndexIt;
        const Coordinates* directionNodesPtr = nodes + *triangleNodeIndexIt;

        double transitionMatrixBuf[4];

        transitionMatrixBuf[0] = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        transitionMatrixBuf[2] = directionNodesPtr->z - tetrahedroBeginNodeIt->z;

        ++triangleNodeIndexIt;
        directionNodesPtr = nodes + *triangleNodeIndexIt;
        transitionMatrixBuf[1] = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        transitionMatrixBuf[3] = directionNodesPtr->z - tetrahedroBeginNodeIt->z;

        double det = transitionMatrixBuf[0] * transitionMatrixBuf[3] - transitionMatrixBuf[1] * transitionMatrixBuf[2];
        
        transitionMatrix[0] =  transitionMatrixBuf[3] / det;
        transitionMatrix[1] = -transitionMatrixBuf[2] / det;
        transitionMatrix[2] = -transitionMatrixBuf[1] / det;
        transitionMatrix[3] =  transitionMatrixBuf[0] / det;

    }


    void computeTransitionMatrix2DExcludeY(const Coordinates nodes[], const size_t* triangleNodeIndexIt, double transitionMatrix[4])
    {
        ++triangleNodeIndexIt;
        const Coordinates* tetrahedroBeginNodeIt = nodes + *triangleNodeIndexIt;
        ++triangleNodeIndexIt;
        const Coordinates* directionNodesPtr = nodes + *triangleNodeIndexIt;

        double transitionMatrixBuf[4];

        transitionMatrixBuf[0] = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        transitionMatrixBuf[2] = directionNodesPtr->z - tetrahedroBeginNodeIt->z;

        ++triangleNodeIndexIt;
        directionNodesPtr = nodes + *triangleNodeIndexIt;
        transitionMatrixBuf[1] = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        transitionMatrixBuf[3] = directionNodesPtr->z - tetrahedroBeginNodeIt->z;

        double det = transitionMatrixBuf[0] * transitionMatrixBuf[3] - transitionMatrixBuf[1] * transitionMatrixBuf[2];

        transitionMatrix[0] = transitionMatrixBuf[3] / det;
        transitionMatrix[1] = -transitionMatrixBuf[2] / det;
        transitionMatrix[2] = -transitionMatrixBuf[1] / det;
        transitionMatrix[3] = transitionMatrixBuf[0] / det;

    }

    void computeTransitionMatrix2DExcludeZ(const Coordinates nodes[], const size_t* triangleNodeIndexIt, double transitionMatrix[4])
    {
        ++triangleNodeIndexIt;
        const Coordinates* tetrahedroBeginNodeIt = nodes + *triangleNodeIndexIt;
        ++triangleNodeIndexIt;
        const Coordinates* directionNodesPtr = nodes + *triangleNodeIndexIt;

        double transitionMatrixBuf[4];

        transitionMatrixBuf[0] = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        transitionMatrixBuf[2] = directionNodesPtr->y - tetrahedroBeginNodeIt->y;

        ++triangleNodeIndexIt;
        directionNodesPtr = nodes + *triangleNodeIndexIt;
        transitionMatrixBuf[1] = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        transitionMatrixBuf[3] = directionNodesPtr->y - tetrahedroBeginNodeIt->y;

        double det = transitionMatrixBuf[0] * transitionMatrixBuf[3] - transitionMatrixBuf[1] * transitionMatrixBuf[2];

        transitionMatrix[0] = transitionMatrixBuf[3] / det;
        transitionMatrix[1] = -transitionMatrixBuf[2] / det;
        transitionMatrix[2] = -transitionMatrixBuf[1] / det;
        transitionMatrix[3] = transitionMatrixBuf[0] / det;

    }

    void translatePointsCoordinatesToLocal2DExcludeX(const Coordinates& basePoint,
                                                     const Coordinates points[], 
                                                     const size_t pointsIndexes[constants::triangle::N_NODES],
                                                     const double transitionMatrix[4], 
                                                     LocalCoordinates2D localCoordinates[constants::triangle::N_NODES])
    {
        double y, z;
        const Coordinates* point;
        for (uint8_t i = 0; i < constants::triangle::N_NODES; ++i)
        {
            point = points + pointsIndexes[i];
            y = point->y - basePoint.y;
            z = point->z - basePoint.z;

            localCoordinates->u = transitionMatrix[0] * y + transitionMatrix[1] * z;
            localCoordinates->v = transitionMatrix[2] * y + transitionMatrix[3] * z;

            ++points;
            ++localCoordinates;
        }
    }

    void translatePointsCoordinatesToLocal2DExcludeY(const Coordinates& basePoint,
                                                     const Coordinates points[],
                                                     const size_t pointsIndexes[constants::triangle::N_NODES],
                                                     const double transitionMatrix[4],
                                                     LocalCoordinates2D localCoordinates[constants::triangle::N_NODES])
    {
        double x, z;
        const Coordinates* point;
        for (uint8_t i = 0; i < constants::triangle::N_NODES; ++i)
        {
            point = points + pointsIndexes[i];
            x = point->x - basePoint.x;
            z = point->z - basePoint.z;

            localCoordinates->u = transitionMatrix[0] * x + transitionMatrix[1] * z;
            localCoordinates->v = transitionMatrix[2] * x + transitionMatrix[3] * z;

            ++points;
            ++localCoordinates;
        }
    }

    void translatePointsCoordinatesToLocal2DExcludeZ(const Coordinates& basePoint,
                                                     const Coordinates points[],
                                                     const size_t pointsIndexes[constants::triangle::N_NODES],
                                                     const double transitionMatrix[4],
                                                     LocalCoordinates2D localCoordinates[constants::triangle::N_NODES])
    {
        double x, y;
        const Coordinates* point;
        for (uint8_t i = 0; i < constants::triangle::N_NODES; ++i)
        {
            point = points + pointsIndexes[i];
            x = point->x - basePoint.x;
            y = point->y - basePoint.y;

            localCoordinates->u = transitionMatrix[0] * x + transitionMatrix[1] * y;
            localCoordinates->v = transitionMatrix[2] * x + transitionMatrix[3] * y;

            ++points;
            ++localCoordinates;
        }
    }

    void computeTranslationCoefficients2D(LocalCoordinates2D triangleLocalNodes[constants::triangle::N_NODES], LocalCoordinates2D translCoeff[2])
    {
        translCoeff->u = triangleLocalNodes[1].u - triangleLocalNodes[0].u;
        translCoeff->v = triangleLocalNodes[1].v - triangleLocalNodes[0].v;

        ++translCoeff;
        translCoeff->u = triangleLocalNodes[1].u - triangleLocalNodes[0].u;
        translCoeff->v = triangleLocalNodes[1].v - triangleLocalNodes[0].v;
    }

    void translateFragmentLocalPoints(const LocalCoordinates2D fragmentBasePoint,
                                     const LocalCoordinates2D* fragmentPointIt,
                                     const uint8_t nPoints,
                                     const LocalCoordinates2D translCoeff[2],
                                     LocalCoordinates2D* trianglePointIt)
    {
        for (uint8_t i = 0; i < nPoints; ++i)
        {
            trianglePointIt->u = fragmentBasePoint.u;
            trianglePointIt->v = fragmentBasePoint.v;

            trianglePointIt->u += fragmentPointIt->u * translCoeff[0].u + fragmentPointIt->v * translCoeff[1].u;
            trianglePointIt->v += fragmentPointIt->u * translCoeff[0].v + fragmentPointIt->v * translCoeff[1].v;

            ++trianglePointIt;
            ++fragmentPointIt;
        }
    }


    void computeTranspJacobianTo0Face(const Coordinates nodes[], const size_t* tetrahedronNodeIndexIt, double transpJacobianMatrix[LocalCoordinates2D::COUNT * Coordinates::COUNT])
    {
        ++tetrahedronNodeIndexIt;
        const Coordinates* tetrahedroBeginNodeIt = nodes + *tetrahedronNodeIndexIt;
        ++tetrahedronNodeIndexIt;
        const Coordinates* directionNodesPtr = nodes + *tetrahedronNodeIndexIt;
        double* transpJacobianMatrixElementIt = transpJacobianMatrix;

        *transpJacobianMatrixElementIt = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedroBeginNodeIt->z;

        ++tetrahedronNodeIndexIt;
        directionNodesPtr = nodes + *tetrahedronNodeIndexIt;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedroBeginNodeIt->z;
    }

    void computeTranspJacobianTo1Face(const Coordinates nodes[], const size_t* tetrahedronNodeIndexIt, double transpJacobianMatrix[LocalCoordinates2D::COUNT * Coordinates::COUNT])
    {
        const Coordinates* tetrahedroBeginNodeIt = nodes + *tetrahedronNodeIndexIt;
        ++tetrahedronNodeIndexIt;
        ++tetrahedronNodeIndexIt;
        const Coordinates* directionNodesPtr = nodes + *tetrahedronNodeIndexIt;
        double* transpJacobianMatrixElementIt = transpJacobianMatrix;

        *transpJacobianMatrixElementIt = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedroBeginNodeIt->z;

        ++tetrahedronNodeIndexIt;
        directionNodesPtr = nodes + *tetrahedronNodeIndexIt;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedroBeginNodeIt->z;
    }

    void computeTranspJacobianTo2Face(const Coordinates nodes[], const size_t* tetrahedronNodeIndexIt, double transpJacobianMatrix[LocalCoordinates2D::COUNT * Coordinates::COUNT])
    {
        const Coordinates* tetrahedroBeginNodeIt = nodes + *tetrahedronNodeIndexIt;
        ++tetrahedronNodeIndexIt;
        const Coordinates* directionNodesPtr = nodes + *tetrahedronNodeIndexIt;
        double* transpJacobianMatrixElementIt = transpJacobianMatrix;

        *transpJacobianMatrixElementIt = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedroBeginNodeIt->z;

        ++tetrahedronNodeIndexIt;
        ++tetrahedronNodeIndexIt;
        directionNodesPtr = nodes + *tetrahedronNodeIndexIt;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedroBeginNodeIt->z;
    }

    void computeTranspJacobianTo3Face(const Coordinates nodes[], const size_t* tetrahedronNodeIndexIt, double transpJacobianMatrix[LocalCoordinates2D::COUNT * Coordinates::COUNT])
    {
        const Coordinates* tetrahedroBeginNodeIt = nodes + *tetrahedronNodeIndexIt;
        ++tetrahedronNodeIndexIt;
        const Coordinates* directionNodesPtr = nodes + *tetrahedronNodeIndexIt;
        double* transpJacobianMatrixElementIt = transpJacobianMatrix;

        *transpJacobianMatrixElementIt = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedroBeginNodeIt->z;

        ++tetrahedronNodeIndexIt;
        directionNodesPtr = nodes + *tetrahedronNodeIndexIt;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedroBeginNodeIt->z;
    }

    void computeNormal(const double transpJacobianMatrix[LocalCoordinates2D::COUNT * Coordinates::COUNT], Coordinates& normal)
    {
        normal.x = transpJacobianMatrix[1] * transpJacobianMatrix[5] - transpJacobianMatrix[4] * transpJacobianMatrix[2];
        normal.y = transpJacobianMatrix[0] * transpJacobianMatrix[5] - transpJacobianMatrix[3] * transpJacobianMatrix[2];
        normal.z = transpJacobianMatrix[0] * transpJacobianMatrix[4] - transpJacobianMatrix[3] * transpJacobianMatrix[1];
    }

    double computeTranspJacobianMatrix(const Coordinates* tetrahedronNodeIt, double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT])
    {
        const Coordinates* tetrahedroBeginNodeIt = tetrahedronNodeIt;
        ++tetrahedronNodeIt;
        double* transpJacobianMatrixElementIt = transpJacobianMatrix;

        *transpJacobianMatrixElementIt = tetrahedronNodeIt->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = tetrahedronNodeIt->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = tetrahedronNodeIt->z - tetrahedroBeginNodeIt->z;

        ++tetrahedronNodeIt;
        *(++transpJacobianMatrixElementIt) = tetrahedronNodeIt->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = tetrahedronNodeIt->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = tetrahedronNodeIt->z - tetrahedroBeginNodeIt->z;

        ++tetrahedronNodeIt;
        *(++transpJacobianMatrixElementIt) = tetrahedronNodeIt->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = tetrahedronNodeIt->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = tetrahedronNodeIt->z - tetrahedroBeginNodeIt->z;

        double det = transpJacobianMatrix[0] * transpJacobianMatrix[4] * transpJacobianMatrix[8];
        det += transpJacobianMatrix[2] * transpJacobianMatrix[3] * transpJacobianMatrix[7];
        det += transpJacobianMatrix[1] * transpJacobianMatrix[5] * transpJacobianMatrix[6];
        det -= transpJacobianMatrix[2] * transpJacobianMatrix[4] * transpJacobianMatrix[6];
        det -= transpJacobianMatrix[0] * transpJacobianMatrix[5] * transpJacobianMatrix[7];
        det -= transpJacobianMatrix[1] * transpJacobianMatrix[3] * transpJacobianMatrix[8];

        return det;
    }

    double computeTranspJacobianMatrix(const Coordinates nodes[], const size_t* tetrahedronNodeIndexIt, double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT])
    {
        const Coordinates* tetrahedroBeginNodeIt = nodes + *tetrahedronNodeIndexIt;
        ++tetrahedronNodeIndexIt;
        const Coordinates* directionNodesPtr = nodes + *tetrahedronNodeIndexIt;
        double* transpJacobianMatrixElementIt = transpJacobianMatrix;

        *transpJacobianMatrixElementIt = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedroBeginNodeIt->z;

        ++tetrahedronNodeIndexIt;
        directionNodesPtr = nodes + *tetrahedronNodeIndexIt;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedroBeginNodeIt->z;

        ++tetrahedronNodeIndexIt;
        directionNodesPtr = nodes + *tetrahedronNodeIndexIt;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->x - tetrahedroBeginNodeIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedroBeginNodeIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedroBeginNodeIt->z;

        double det = transpJacobianMatrix[0] * transpJacobianMatrix[4] * transpJacobianMatrix[8];
        det += transpJacobianMatrix[2] * transpJacobianMatrix[3] * transpJacobianMatrix[7];
        det += transpJacobianMatrix[1] * transpJacobianMatrix[5] * transpJacobianMatrix[6];
        det -= transpJacobianMatrix[2] * transpJacobianMatrix[4] * transpJacobianMatrix[6];
        det -= transpJacobianMatrix[0] * transpJacobianMatrix[5] * transpJacobianMatrix[7];
        det -= transpJacobianMatrix[1] * transpJacobianMatrix[3] * transpJacobianMatrix[8];

        return det;
    }


    double computeTranspJacobianMatrix(const Coordinates* tetrahedronNodeBeginIt, double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT])
    {
        const Coordinates* directionNodesPtr = tetrahedronNodeBeginIt + 1;
        double* transpJacobianMatrixElementIt = transpJacobianMatrix;
        *transpJacobianMatrixElementIt = directionNodesPtr->x - tetrahedronNodeBeginIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedronNodeBeginIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedronNodeBeginIt->z;

        ++directionNodesPtr;
        *transpJacobianMatrixElementIt = directionNodesPtr->x - tetrahedronNodeBeginIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedronNodeBeginIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedronNodeBeginIt->z;

        ++directionNodesPtr;
        *transpJacobianMatrixElementIt = directionNodesPtr->x - tetrahedronNodeBeginIt->x;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->y - tetrahedronNodeBeginIt->y;
        *(++transpJacobianMatrixElementIt) = directionNodesPtr->z - tetrahedronNodeBeginIt->z;

        double det = transpJacobianMatrix[0] * transpJacobianMatrix[4] * transpJacobianMatrix[8];
        det += transpJacobianMatrix[2] * transpJacobianMatrix[3] * transpJacobianMatrix[7];
        det += transpJacobianMatrix[1] * transpJacobianMatrix[5] * transpJacobianMatrix[6];
        det -= transpJacobianMatrix[2] * transpJacobianMatrix[4] * transpJacobianMatrix[6];
        det -= transpJacobianMatrix[0] * transpJacobianMatrix[5] * transpJacobianMatrix[7];
        det -= transpJacobianMatrix[1] * transpJacobianMatrix[3] * transpJacobianMatrix[8];

        return det;
    }
    void computeLocalJacobianMatrix(const double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT], const double determinant, double* localJacobianMatrixElementIt)
    {
        *localJacobianMatrixElementIt = (transpJacobianMatrix[4] * transpJacobianMatrix[8] - transpJacobianMatrix[5] * transpJacobianMatrix[7]) / determinant;
        ++localJacobianMatrixElementIt;
        *localJacobianMatrixElementIt = (transpJacobianMatrix[5] * transpJacobianMatrix[6] - transpJacobianMatrix[3] * transpJacobianMatrix[8]) / determinant;
        ++localJacobianMatrixElementIt;
        *localJacobianMatrixElementIt = (transpJacobianMatrix[3] * transpJacobianMatrix[7] - transpJacobianMatrix[4] * transpJacobianMatrix[6]) / determinant;
        ++localJacobianMatrixElementIt;
        
        *localJacobianMatrixElementIt = (transpJacobianMatrix[2] * transpJacobianMatrix[7] - transpJacobianMatrix[1] * transpJacobianMatrix[8]) / determinant;
        ++localJacobianMatrixElementIt;
        *localJacobianMatrixElementIt = (transpJacobianMatrix[0] * transpJacobianMatrix[8] - transpJacobianMatrix[2] * transpJacobianMatrix[6]) / determinant;
        ++localJacobianMatrixElementIt;
        *localJacobianMatrixElementIt = (transpJacobianMatrix[1] * transpJacobianMatrix[6] - transpJacobianMatrix[0] * transpJacobianMatrix[7]) / determinant;
        ++localJacobianMatrixElementIt;
        
        *localJacobianMatrixElementIt = (transpJacobianMatrix[1] * transpJacobianMatrix[5] - transpJacobianMatrix[2] * transpJacobianMatrix[4]) / determinant;
        ++localJacobianMatrixElementIt;
        *localJacobianMatrixElementIt = (transpJacobianMatrix[2] * transpJacobianMatrix[3] - transpJacobianMatrix[0] * transpJacobianMatrix[5]) / determinant;
        ++localJacobianMatrixElementIt;
        *localJacobianMatrixElementIt = (transpJacobianMatrix[0] * transpJacobianMatrix[4] - transpJacobianMatrix[1] * transpJacobianMatrix[3]) / determinant;
        ++localJacobianMatrixElementIt;
    }

    void copyTetrahedronsFacePoints(const uint8_t faceLocalIndex,
                                    const Coordinates tetrahedronPoints[constants::tetrahedron::N_NODES],
                                    Coordinates facePoints[constants::tetrahedron::N_NODES])
    {
        switch (faceLocalIndex)
        {
        case 0:
        {
            *facePoints = tetrahedronPoints[0];
            *(++facePoints) = tetrahedronPoints[2];
            *(++facePoints) = tetrahedronPoints[1];
            break;
        }
        case 1:
        {
            *facePoints = tetrahedronPoints[0];
            *(++facePoints) = tetrahedronPoints[1];
            *(++facePoints) = tetrahedronPoints[3];
            break;
        }
        case 2:
        {
            *facePoints = tetrahedronPoints[0];
            *(++facePoints) = tetrahedronPoints[3];
            *(++facePoints) = tetrahedronPoints[2];
            break;
        }
        case 3:
        {
            *facePoints = tetrahedronPoints[3];
            *(++facePoints) = tetrahedronPoints[1];
            *(++facePoints) = tetrahedronPoints[2];
            break;
        }
        }
    }

    void computeTetrahedronFaceDirections(const uint8_t faceLocalIndex, const Coordinates tetrahedronPoints[constants::tetrahedron::N_NODES], Coordinates directions[2])
    {
        switch (faceLocalIndex)
        {
        case 0:
        {
            computeDiffrence(tetrahedronPoints[2], tetrahedronPoints[0], directions[0]);
            computeDiffrence(tetrahedronPoints[1], tetrahedronPoints[0], directions[1]);
            break;
        }
        case 1:
        {
            computeDiffrence(tetrahedronPoints[1], tetrahedronPoints[0], directions[0]);
            computeDiffrence(tetrahedronPoints[3], tetrahedronPoints[0], directions[1]);
            break;
        }
        case 2:
        {
            computeDiffrence(tetrahedronPoints[3], tetrahedronPoints[0], directions[0]);
            computeDiffrence(tetrahedronPoints[2], tetrahedronPoints[0], directions[1]);
            break;
        }
        case 3:
        {
            computeDiffrence(tetrahedronPoints[1], tetrahedronPoints[3], directions[0]);
            computeDiffrence(tetrahedronPoints[2], tetrahedronPoints[3], directions[1]);
            break;
        }
        }
    }

    double computeTriagnleJacobianDet(const Coordinates triangleDirections[2])
    {
        const Coordinates* triangleDirection2 = triangleDirections + 1;
        double e = triangleDirections->x * triangleDirections->x + triangleDirections->y * triangleDirections->y + triangleDirections->z * triangleDirections->z;
        double g = triangleDirection2->x * triangleDirection2->x + triangleDirection2->y * triangleDirection2->y + triangleDirection2->z * triangleDirection2->z;
        double f = triangleDirections->x * triangleDirection2->x + triangleDirections->y * triangleDirection2->y + triangleDirections->z * triangleDirection2->z;

        return sqrt(e * g - f * f);
    }

    void computeNormalViaTriangleDirections(const Coordinates triangleDirections[2], Coordinates& normal)
    {
        const Coordinates* triangleDirection2 = triangleDirections + 1;
        normal.x = triangleDirections->y * triangleDirection2->z - triangleDirections->z * triangleDirection2->y;
        normal.y = triangleDirections->z * triangleDirection2->x - triangleDirections->x * triangleDirection2->z;
        normal.z = triangleDirections->x * triangleDirection2->y - triangleDirections->y * triangleDirection2->x;

        double norm = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

        normal.x /= norm;
        normal.y /= norm;
        normal.z /= norm;
    }

    void computeLocalJacobianMatrixies(const double(*transpJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
                                       const size_t nJacobianMatrixies,
                                       const double* determinnatIt,
                                       double(*localJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT])
    { 
        for (size_t i = 0; i < nJacobianMatrixies; ++i)
        {
            computeLocalJacobianMatrix(*transpJacobianMatrixIt, *determinnatIt, *localJacobianMatrixIt);
            
            ++transpJacobianMatrixIt;
            ++determinnatIt;
            ++localJacobianMatrixIt;
        }
    }

    void translate(const double localJacobian[LocalCoordinates3D::COUNT * Coordinates::COUNT], const LocalCoordinates3D &localGradient, Coordinates &gradientIt)
    {
        gradientIt.x = localJacobian[0] * localGradient.u + localJacobian[1] * localGradient.v + localJacobian[2] * localGradient.w;
        gradientIt.y = localJacobian[3] * localGradient.u + localJacobian[4] * localGradient.v + localJacobian[5] * localGradient.w;
        gradientIt.z = localJacobian[6] * localGradient.u + localJacobian[7] * localGradient.v + localJacobian[8] * localGradient.w;
    }


    void translate(const double localJacobian[LocalCoordinates3D::COUNT * Coordinates::COUNT], const LocalCoordinates3D* localGradientIt, const uint16_t nPoints, Coordinates* gradientIt)
	{
		for (uint16_t i = 0; i < nPoints; ++i)
		{
			gradientIt->x = localJacobian[0] * localGradientIt->u + localJacobian[1] * localGradientIt->v + localJacobian[2] * localGradientIt->w;
			gradientIt->y = localJacobian[3] * localGradientIt->u + localJacobian[4] * localGradientIt->v + localJacobian[5] * localGradientIt->w;
			gradientIt->z = localJacobian[6] * localGradientIt->u + localJacobian[7] * localGradientIt->v + localJacobian[8] * localGradientIt->w;

			++gradientIt;
			++localGradientIt;
		}
	}

    void translate0FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints)
    {
        for (uint8_t i = 0; i < nPoints; ++i)
        {
            LocalCoordinates3D& localPoint3D = tetrahedraLocalPoints[i];
            const LocalCoordinates2D& templateLocalPoint2D = templateFaceLocalPoints[i];
            localPoint3D.u = templateLocalPoint2D.v;
            localPoint3D.v = templateLocalPoint2D.u;
            localPoint3D.w = 0.0;
        }
    }

    void translate1FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints)
    {
        for (uint8_t i = 0; i < nPoints; ++i)
        {
            LocalCoordinates3D& localPoint3D = tetrahedraLocalPoints[i];
            const LocalCoordinates2D& templateLocalPoint2D = templateFaceLocalPoints[i];
            localPoint3D.u = templateLocalPoint2D.u;
            localPoint3D.v = 0.0;
            localPoint3D.w = templateLocalPoint2D.v;
        }
    }

    void translate2FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints)
    {
        for (uint8_t i = 0; i < nPoints; ++i)
        {
            LocalCoordinates3D& localPoint3D = tetrahedraLocalPoints[i];
            const LocalCoordinates2D& templateLocalPoint2D = templateFaceLocalPoints[i];
            localPoint3D.u = 0.0;
            localPoint3D.v = templateLocalPoint2D.v;
            localPoint3D.w = templateLocalPoint2D.u;
        }
    }

    void translate3FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints)
    {
        for (uint8_t i = 0; i < nPoints; ++i)
        {
            LocalCoordinates3D& localPoint3D = tetrahedraLocalPoints[i];
            const LocalCoordinates2D& templateLocalPoint2D = templateFaceLocalPoints[i];
            localPoint3D.u = templateLocalPoint2D.u;
            localPoint3D.v = templateLocalPoint2D.v;
            localPoint3D.w = 1.0 - templateLocalPoint2D.u - templateLocalPoint2D.v;
        }
    }

    void translate(const LocalCoordinates3D localPoints[],
                   const uint8_t pointsCount,
                   const Coordinates& initPoint,
                   const double transpJacobianMatrix[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                   Coordinates points[])
    {
        for (uint8_t i = 0; i < pointsCount; ++i)
        {
            Coordinates& point = points[i];
            point = initPoint;
            const LocalCoordinates3D& localPoint = localPoints[i];

            point.x += transpJacobianMatrix[0] * localPoint.u;
            point.y += transpJacobianMatrix[1] * localPoint.u;
            point.z += transpJacobianMatrix[2] * localPoint.u;

            point.x += transpJacobianMatrix[3] * localPoint.v;
            point.y += transpJacobianMatrix[4] * localPoint.v;
            point.z += transpJacobianMatrix[5] * localPoint.v;

            point.x += transpJacobianMatrix[6] * localPoint.w;
            point.y += transpJacobianMatrix[7] * localPoint.w;
            point.z += transpJacobianMatrix[8] * localPoint.w;
        }
    }

    void translate(const LocalCoordinates2D localPoints[],
                   const uint8_t pointsCount,
                   const Coordinates& initPoint,
                   const double transpJacobianMatrix[Coordinates::COUNT * LocalCoordinates3D::COUNT],
                   Coordinates points[])
    {
        for (uint8_t i = 0; i < pointsCount; ++i)
        {
            Coordinates& point = points[i];
            point = initPoint;
            const LocalCoordinates2D& localPoint = localPoints[i];

            point.x += transpJacobianMatrix[0] * localPoint.u;
            point.y += transpJacobianMatrix[1] * localPoint.u;
            point.z += transpJacobianMatrix[2] * localPoint.u;

            point.x += transpJacobianMatrix[3] * localPoint.v;
            point.y += transpJacobianMatrix[4] * localPoint.v;
            point.z += transpJacobianMatrix[5] * localPoint.v;
        }
    }

    void coomputeDirectionalDerivative(const Coordinates* gradientIt, const uint16_t nGradients, const Coordinates direction, double* directionalDerivativeIt)
    {
        for (uint8_t i = 0; i < nGradients; ++i)
        {
            *directionalDerivativeIt = gradientIt->x * direction.x + gradientIt->y * direction.y + gradientIt->z * direction.z;
            
            ++gradientIt;
            ++directionalDerivativeIt;
        }
    }

    double computeTriangleDeterminant(const Coordinates points[], const size_t trianglePointsIndexes[constants::triangle::N_NODES])
    {
        Coordinates vector1, vector2, normal;

        computeDiffrence(points[trianglePointsIndexes[1]], points[trianglePointsIndexes[0]], vector1);
        computeDiffrence(points[trianglePointsIndexes[2]], points[trianglePointsIndexes[0]], vector2);

        normal.x = vector1.y * vector2.z - vector1.z * vector2.y;
        normal.y = vector1.z * vector2.x - vector1.x * vector2.z;
        normal.z = vector1.x * vector2.y - vector1.y * vector2.x;

        return sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    }

    void computeNormal(const Coordinates trianglePoints[constants::triangle::N_NODES], Coordinates &normal)
    {
        Coordinates vector1, vector2;

        computeDiffrence(trianglePoints[1], trianglePoints[0], vector1);
        computeDiffrence(trianglePoints[2], trianglePoints[0], vector2);

        normal.x = vector1.y * vector2.z - vector1.z * vector2.y;
        normal.y = vector1.z * vector2.x - vector1.x * vector2.z;
        normal.z = vector1.x * vector2.y - vector1.y * vector2.x;
    }

}

