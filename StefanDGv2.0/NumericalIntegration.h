#pragma once
#include "LocalCoordinates3D.h"
#include "LocalCoordinates2D.h"
#include "Coordinates.h"
#include "GeometryConstants.h"
#include <cstdint>


namespace NumericalIntegration
{
	namespace Tetrahedron
	{
		template <size_t ORDER> requires (ORDER == 2)
		class Gauss
		{
		public:
			Gauss() = default;
			static const LocalCoordinates3D localPoints[];
			static const double weights[];
			static const uint8_t nSteps;
		};

		template<>
		class Gauss<2>
		{
		public:
			static const uint8_t nSteps = 4;
			static const LocalCoordinates3D localPoints[nSteps];
			static const double weights[nSteps];

			static double integrate(const double values[nSteps]);
			static double integrateProduct(const double multipliers1[nSteps], const double multipliers2[nSteps]);
			static double integrateProduct(const Coordinates vectors1[nSteps], const Coordinates vectors2[nSteps]);

		private:
			Gauss<2>() = default;
			static const double coord1;
			static const double coord2;
			static const double weight;
		};

	}

	namespace Triangle
	{
		template <size_t ORDER> requires (ORDER == 2)
			class Gauss
		{
		public:
			Gauss() = default;
			static const LocalCoordinates2D localPoints[];
			static const double weights[];
			static const uint8_t nSteps;
		};

		template<>
		class Gauss<2>
		{
		public:
			static const uint8_t nSteps = 3;
			static const LocalCoordinates2D localPoints[nSteps];
			static const double weights[nSteps];

			static double integrate(const double values[nSteps]);
			static double integrateProduct(const double multipliers1[nSteps], const double multipliers2[nSteps]);

			static void changeValuesOrder(double values[], const uint8_t nFunctions, const uint8_t nodesChangins[constants::triangle::N_NODES]);
			static void changeValuesOrder(double values[nSteps], const uint8_t nodesChangins[constants::triangle::N_NODES]);
			static void changeValuesOrder(const double values[nSteps], const uint8_t nodesChangins[constants::triangle::N_NODES], double outValues[nSteps]);

			static void changeValuesOrder(const double values[nSteps], const uint8_t nodesChangins, double outValues[nSteps]);
			static void changeValuesOrder(double values[], const uint8_t nFunctions, const uint8_t nodesChangins);

		private:
			Gauss<2>() = default;
			static const double coord1;
			static const double coord2;
			static const double weight;
		};

	}


}
