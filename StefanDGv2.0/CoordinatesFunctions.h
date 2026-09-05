#pragma once
#include "LocalCoordinates2D.h"
#include "LocalCoordinates3D.h"
#include "Coordinates.h"
#include "GeometryConstants.h"

namespace CoordinatesFunctions
{
	void computeDiffrence(const Coordinates& point1, const Coordinates& point2, Coordinates& vector);

	void computeTransitionMatrix2DExcludeX(const Coordinates nodes[], const size_t* triangleNodeIndexIt, double transitionMatrix[4]);
	void computeTransitionMatrix2DExcludeY(const Coordinates nodes[], const size_t* triangleNodeIndexIt, double transitionMatrix[4]);
	void computeTransitionMatrix2DExcludeZ(const Coordinates nodes[], const size_t* triangleNodeIndexIt, double transitionMatrix[4]);

	void translatePointsCoordinatesToLocal2DExcludeX(const Coordinates& basePoint,
												     const Coordinates points[],
												     const size_t pointsIndexes[constants::triangle::N_NODES],
												     const double transitionMatrix[4],
												     LocalCoordinates2D localCoordinates[constants::triangle::N_NODES]);
	void translatePointsCoordinatesToLocal2DExcludeY(const Coordinates& basePoint,
													 const Coordinates points[],
													 const size_t pointsIndexes[constants::triangle::N_NODES],
													 const double transitionMatrix[4],
													 LocalCoordinates2D localCoordinates[constants::triangle::N_NODES]);
	void translatePointsCoordinatesToLocal2DExcludeZ(const Coordinates& basePoint,
													 const Coordinates points[],
													 const size_t pointsIndexes[constants::triangle::N_NODES],
													 const double transitionMatrix[4],
													 LocalCoordinates2D localCoordinates[constants::triangle::N_NODES]);

	void computeTranslationCoefficients2D(LocalCoordinates2D triangleLocalNodes[constants::triangle::N_NODES], LocalCoordinates2D translCoeff[2]);
	
	void translateFragmentLocalPoints(const LocalCoordinates2D fragmentBasePoint,
									  const LocalCoordinates2D* fragmentPointIt,
									  const uint8_t nPoints,
									  const LocalCoordinates2D translCoeff[2],
									  LocalCoordinates2D* trianglePointIt);


	void computeTranspJacobianTo0Face(const Coordinates nodes[], const size_t* tetrahedronNodeIndexIt, double transpJacobianMatrix[LocalCoordinates2D::COUNT * Coordinates::COUNT]);
	void computeTranspJacobianTo1Face(const Coordinates nodes[], const size_t* tetrahedronNodeIndexIt, double transpJacobianMatrix[LocalCoordinates2D::COUNT * Coordinates::COUNT]);
	void computeTranspJacobianTo2Face(const Coordinates nodes[], const size_t* tetrahedronNodeIndexIt, double transpJacobianMatrix[LocalCoordinates2D::COUNT * Coordinates::COUNT]);
	void computeTranspJacobianTo3Face(const Coordinates nodes[], const size_t* tetrahedronNodeIndexIt, double transpJacobianMatrix[LocalCoordinates2D::COUNT * Coordinates::COUNT]);

	void computeNormal(const double transpJacobianMatrix[LocalCoordinates2D::COUNT * Coordinates::COUNT], Coordinates& normal);

	double computeTranspJacobianMatrix(const Coordinates* tetrahedronNodeIt, double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT]);
	double computeTranspJacobianMatrix(const Coordinates nodes[], const size_t* tetrahedronNodeIndexIt, double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT]);
	double computeTranspJacobianMatrix(const Coordinates* tetrahedronNodeBeginIt, double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT]);

	void computeLocalJacobianMatrix(const double transpJacobianMatrix[LocalCoordinates3D::COUNT * Coordinates::COUNT],
		const double determinant,
		double* localJacobianMatrixElementIt);

	void computeLocalJacobianMatrixies(const double(*transpJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT],
		const size_t nJacobianMatrixies,
		const double* determinnatIt,
		double(*localJacobianMatrixIt)[LocalCoordinates3D::COUNT * Coordinates::COUNT]);

	void copyTetrahedronsFacePoints(const uint8_t faceLocalIndex,
		const Coordinates tetrahedronPoints[constants::tetrahedron::N_NODES],
		Coordinates facePoints[constants::tetrahedron::N_NODES]);

	void computeTetrahedronFaceDirections(const uint8_t faceLocalIndex, const Coordinates tetrahedronPoints[constants::tetrahedron::N_NODES], Coordinates directions[2]);
	double computeTriagnleJacobianDet(const Coordinates triangleDirections[2]);
	void computeNormalViaTriangleDirections(const Coordinates triangleDirections[2], Coordinates& normal);

	void translate(const double localJacobian[LocalCoordinates3D::COUNT * Coordinates::COUNT], const LocalCoordinates3D& localGradient, Coordinates& gradientIt);
	void translate(const double localJacobian[LocalCoordinates3D::COUNT * Coordinates::COUNT], const LocalCoordinates3D* localGradientIt, const uint16_t nPoints, Coordinates* gradientIt);

	void translate0FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints);
	void translate1FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints);
	void translate2FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints);
	void translate3FacePoints(const LocalCoordinates2D* templateFaceLocalPoints, const uint8_t nPoints, LocalCoordinates3D* tetrahedraLocalPoints);

	void translate(const LocalCoordinates3D localPoints[],
				   const uint8_t pointsCount,
				   const Coordinates& initPoint,
				   const double transpJacobianMatrix[Coordinates::COUNT * LocalCoordinates3D::COUNT],
				   Coordinates points[]);

	void translate(const LocalCoordinates2D localPoints[],
				   const uint8_t pointsCount,
				   const Coordinates& initPoint,
				   const double transpJacobianMatrix[Coordinates::COUNT * LocalCoordinates3D::COUNT],
				   Coordinates points[]);

	void coomputeDirectionalDerivative(const Coordinates* gradientIt, const uint16_t nGradients, const Coordinates direction, double* directionalDerivativeIt);

	double computeTriangleDeterminant(const Coordinates points[], const size_t trianglePointsIndexes[constants::triangle::N_NODES]);
	void computeNormal(const Coordinates trianglePoints[constants::triangle::N_NODES], Coordinates& normal);
};

