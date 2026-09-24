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
		void initSolver();
		void finalizeSolver();

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

		double* computeInitialDOFs(const unsigned int nRegions, const MaterialPhase* const regionsMaterialPhases[], size_t* regionsStartTetrahedronsIndexes);

		double solveInitialIteration(const unsigned int nRegions,
								   const MaterialPhase* const regionsMaterialPhases[],
								   const Boundary boundaries[],
								   const unsigned int nBoundaries,
								   const NonconformInterface nonconformInterfaces[],
								   const unsigned int nNonconformInterfaces,
								   const double t,
								   const double dt,
								   const double penalty,
								   void* calculationBuffer,
								   void* additionalBuffer,
								   Solution* solutionIt);

		double solveIteration(const unsigned int nRegions,
						    const MaterialPhase* const regionsMaterialPhases[],
						    const Boundary boundaries[],
						    const unsigned int nBoundaries,
						    const NonconformInterface nonconformInterfaces[],
						    const unsigned int nNonconformInterfaces,
							const double t,
						    const double dt,
						    const double penalty,
						    const Solution& solution,
						    void* calculationBuffer,
						    void* additionalBuffer,
						    Solution* solutionIt);

		void relocateInitialFrontNodes(MaterialPhase solidMaterialsPhase,
									   MaterialPhase liquidMaterialsPhase,
									   const double dt,
									   const double latentHeat,
									   const Coordinates* normalIt,
									   const unsigned int nFrontNodes,
									   Coordinates* frontNodesIt);

		void relocateFrontNodes(const Solution& solution,
								int solidRegionTag,
								int liquidRegionTag,
								MaterialPhase solidMaterialsPhase,
								MaterialPhase liquidMaterialsPhase,
								const double dt,
								const double latentHeat,
								const Coordinates* normalIt,
								const unsigned int nFrontNodes,
								Coordinates* frontNodesIt);

		void getFrontNodes(const int frontTag, int** frontNodesTags,unsigned int* nFrontNodes);
	}


	class Solution
	{
	public:
		friend double StefanTask::solveInitialIteration(const unsigned int nRegions,
													  const MaterialPhase* const regionsMaterialPhases[],
													  const Boundary boundaries[],
													  const unsigned int nBoundaries,
													  const NonconformInterface nonconformInterfaces[],
													  const unsigned int nNonconformInterfaces,
												      const double t,
													  const double dt,
													  const double penalty,
													  void* calculationBuffer,
													  void* additionalBuffer,
													  Solution* solutionIt);

		friend double StefanTask::solveIteration(const unsigned int nRegions,
											   const MaterialPhase* const regionsMaterialPhases[],
											   const Boundary boundaries[],
											   const unsigned int nBoundaries,
											   const NonconformInterface nonconformInterfaces[],
											   const unsigned int nNonconformInterfaces,
											   const double t,
											   const double dt,
											   const double penalty,
											   const Solution& solution,
											   void* calculationBuffer,
											   void* additionalBuffer,
											   Solution* solutionIt);

		double compute(const size_t elementTag, const LocalCoordinates3D LocalPoint3D) const;
		void compute(const size_t elementTag, const LocalCoordinates3D LocalPoint3D, LocalCoordinates3D& gradient) const;
		void clear();
		const double* getDOFs() const;

		Solution() = default;
	private:
		Eigen::VectorXd _DOFs;
		double* _DOFsPtr;
	};
}


