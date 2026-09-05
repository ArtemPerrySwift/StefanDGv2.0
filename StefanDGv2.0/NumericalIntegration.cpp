#include "NumericalIntegration.h"
#include <cmath>

namespace NumericalIntegration
{
	namespace Tetrahedron
	{
		const double Gauss<2>::coord1 = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
		const double Gauss<2>::coord2 = (5.0 - sqrt(5.0)) / 20.0;

		const double Gauss<2>::weight = 1.0 / 24.0;

		const LocalCoordinates3D Gauss<2>::localPoints[] = { {coord1, coord2, coord2},
													 {coord2, coord1, coord2},
													 {coord2, coord2, coord1},
													 {coord2, coord2, coord2} };

		const double Gauss<2>::weights[] = { weight, weight, weight, weight };

		double Gauss<2>::integrate(const double values[nSteps])
		{
			double sum = *values;
			sum += *(++values);
			sum += *(++values);
			sum += *(++values);
			return sum * weight;
		}

		double Gauss<2>::integrateProduct(const double multipliers1[nSteps], const double multipliers2[nSteps])
		{
			double sum = (*multipliers1) * (*multipliers2);
			sum += (*(++multipliers1)) * (*(++multipliers2));
			sum += (*(++multipliers1)) * (*(++multipliers2));
			sum += (*(++multipliers1)) * (*(++multipliers2));

			return sum * weight;
		}

		double Gauss<2>::integrateProduct(const Coordinates vectors1[nSteps], const Coordinates vectors2[nSteps])
		{
			double sum = vectors1->x * vectors2->x + vectors1->y * vectors2->y + vectors1->z * vectors2->z;
			++vectors1;
			++vectors2;

			sum += vectors1->x * vectors2->x + vectors1->y * vectors2->y + vectors1->z * vectors2->z;
			++vectors1;
			++vectors2;

			sum += vectors1->x * vectors2->x + vectors1->y * vectors2->y + vectors1->z * vectors2->z;
			return sum * weight;
		}
	}

	namespace Triangle
	{
		const double Gauss<2>::coord1 = 1.0 / 6.0;
		const double Gauss<2>::coord2 = 2.0 / 3.0;
		const double Gauss<2>::weight = 1.0 / 6.0;

		const LocalCoordinates2D Gauss<2>::localPoints[] = { {coord1, coord1},
													 {coord2, coord1},
													 {coord1, coord2}};

		const double Gauss<2>::weights[] = { weight, weight, weight};

		double Gauss<2>::integrate(const double values[nSteps])
		{
			double sum = *values;
			sum += *(++values);
			sum += *(++values);

			return sum * weight;
		}

		double Gauss<2>::integrateProduct(const double multipliers1[nSteps], const double multipliers2[nSteps])
		{
			double sum = (*multipliers1) * (*multipliers2);
			sum += (*(++multipliers1)) * (*(++multipliers2));
			sum += (*(++multipliers1)) * (*(++multipliers2));

			return sum * weight;
		}

		void Gauss<2>::changeValuesOrder(double values[], const uint8_t nFunctions, const uint8_t nodesChangins[constants::triangle::N_NODES])
		{
			double outValues[constants::triangle::N_NODES];
			for (uint8_t i = 0; i < nFunctions; ++i)
			{
				changeValuesOrder(values, nodesChangins, outValues);

				*values = outValues[0];
				*(++values) = outValues[1];
				*(++values) = outValues[2];
				++values;
			}

		}

		void Gauss<2>::changeValuesOrder(double values[nSteps], const uint8_t nodesChangins[constants::triangle::N_NODES])
		{
			double outValues[constants::triangle::N_NODES];

			changeValuesOrder(values, nodesChangins, outValues);

			values[0] = outValues[0];
			values[1] = outValues[1];
			values[2] = outValues[2];
		}

		void Gauss<2>::changeValuesOrder(const double values[nSteps], const uint8_t nodesChangins[constants::triangle::N_NODES], double outValues[nSteps])
		{
			outValues[*nodesChangins] = *values;
			outValues[*(++nodesChangins)] = *(++values);
			outValues[*(++nodesChangins)] = *(++values);
		}

		void Gauss<2>::changeValuesOrder(const double values[nSteps], const uint8_t nodesChangins, double outValues[nSteps])
		{
			outValues[nodesChangins & 0b11] = *values;
			outValues[(nodesChangins >> 2) & 0b11] = *(++values);
			outValues[nodesChangins >> 4] = *(++values);
		}

		void Gauss<2>::changeValuesOrder(double values[], const uint8_t nFunctions, const uint8_t nodesChangins)
		{
			double outValues[constants::triangle::N_NODES];
			for (uint8_t i = 0; i < nFunctions; ++i)
			{
				changeValuesOrder(values, nodesChangins, outValues);

				*values = outValues[0];
				*(++values) = outValues[1];
				*(++values) = outValues[2];
				++values;
			}

		}
	}
}