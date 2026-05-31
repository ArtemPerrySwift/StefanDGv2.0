#include "DG.h"
#include "LinearLagrangeBasis.h"

namespace DG
{
	Solution* solveStefanTask(const double tMin, const double tMax, const size_t nTSteps)
	{
		Solution* solutions = new Solution[nTSteps];

		const size_t nTIntervals = nTSteps - 1;
		double dt = (tMax - tMin) / nTIntervals;

		Solution* solutionIt = solutions;
		const uint8_t DEBUG_DOFS_SIZE = 4;

		for (size_t i = 0; i < nTIntervals; ++i)
		{
			solutionIt->_elementStartTag = 0;
			solutionIt->_DOFs.resize(DEBUG_DOFS_SIZE);
			double* DOFIt = solutions->_DOFs.data();
			*DOFIt = 0.0;
			++DOFIt;
			*DOFIt = 1.0;
			++DOFIt;
			*DOFIt = 1.0;
			++DOFIt;
			*DOFIt = 1.0;

			++solutionIt;
		}

		solutionIt->_elementStartTag = 0;
		solutionIt->_DOFs.resize(DEBUG_DOFS_SIZE);
		double* DOFIt = solutions->_DOFs.data();
		*DOFIt = 0.0;
		++DOFIt;
		*DOFIt = 1.0;
		++DOFIt;
		*DOFIt = 1.0;
		++DOFIt;
		*DOFIt = 1.0;

		++solutionIt;

		return solutions;
	}

	double Solution::compute(const size_t elementTag, const LocalCoordinates3D LocalPoint3D) const
	{
		return LinearLagrangeBasis::compute(LocalPoint3D, _DOFs.data() + LinearLagrangeBasis::N_FUNCTIONS * (elementTag - _elementStartTag));
	}
}

