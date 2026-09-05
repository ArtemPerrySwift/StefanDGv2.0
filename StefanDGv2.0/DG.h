#include <Eigen/Dense>
#include "LocalCoordinates3D.h"
#include "MaterialPhase.h"
#include "Boundary.h"
#include "NonconformInterface.h"

namespace DG
{
	class Solution;

	namespace StefanTask
	{
		void solve(const unsigned int nRegions,
			const MaterialPhase* const regionsMaterialPhases[],
			const Boundary boundaries[],
			const unsigned int nBoundaries,
			const NonconformInterface nonconformInterfaces[],
			const unsigned int nNonconformInterfaces,
			const double tMin,
			const double tMax,
			const size_t nTSteps,
			Solution* solutionIt);

		void solveInitialIteration(const unsigned int nRegions,
			const MaterialPhase* const regionsMaterialPhases[],
			const Boundary boundaries[],
			const unsigned int nBoundaries,
			const NonconformInterface nonconformInterfaces[],
			const unsigned int nNonconformInterfaces,
			const double dt,
			const double penalty,
			void* calculationBuffer,
			void* additionalBuffer,
			Solution* solutionIt);

	}


	class Solution
	{
	public:
		friend 	void solveStefanTask(const unsigned int nRegions,
			const MaterialPhase* const regionsMaterialPhases[],
			const Boundary boundaries[],
			const unsigned int nBoundaries,
			const NonconformInterface nonconformInterfaces[],
			const unsigned int nNonconformInterfaces,
			const double tMin,
			const double tMax,
			const size_t nTSteps,
			Solution* solutionIt);

		friend void StefanTask::solveInitialIteration(const unsigned int nRegions,
			const MaterialPhase* const regionsMaterialPhases[],
			const Boundary boundaries[],
			const unsigned int nBoundaries,
			const NonconformInterface nonconformInterfaces[],
			const unsigned int nNonconformInterfaces,
			const double dt,
			const double penalty,
			void* calculationBuffer,
			void* additionalBuffer,
			Solution* solutionIt);

		double compute(const size_t elementTag, const LocalCoordinates3D LocalPoint3D) const;
		void compute(const size_t elementTag, const LocalCoordinates3D LocalPoint3D, LocalCoordinates3D& gradient) const;

	private:
		Solution() = default;
		Eigen::VectorXd _DOFs;
		double* _DOFsPtr;
	};
}


