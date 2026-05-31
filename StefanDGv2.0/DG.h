#include <Eigen/Dense>
#include "LocalCoordinates3D.h"

namespace DG
{
	class Solution;

	Solution* solveStefanTask(const double tMin, const double tMax, const size_t nTSteps);

	class Solution
	{
	public:
		friend Solution* solveStefanTask(const double tMin, const double tMax, const size_t nTSteps);
		double compute(const size_t elementTag, const LocalCoordinates3D LocalPoint3D) const;

	private:
		Solution() = default;
		Eigen::VectorXd _DOFs;
		size_t _elementStartTag;
	};
}


