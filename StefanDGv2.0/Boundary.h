#pragma once
#include "Coordinates.h"
#include "LocalCoordinates2D.h"
#include "LocalCoordinates3D.h"
#include "GeometryConstants.h"

struct Boundary
{
	static const unsigned char PLANE_SURFACE = 16;

	struct ImplicitCondition;
	struct FunctionCondition;
	struct ValueCondition;
	struct ConformCondition;

	struct Condition
	{
		friend struct ImplicitCondition;
		friend struct FunctionCondition;
		friend struct ValueCondition;
		friend struct ConformCondition;

		enum class MacroType: uint8_t { DIRICHLET_F, DIRICHLET_V, NEWMAN_F, NEWMAN_V, STEFAN_V, CONFORM_INTERFACE, NONCONFORM_INTERFACE, HOMOGENEOUS_NEWMAN};
		const MacroType macroType;

	private:
		Condition(MacroType macroType);
	};

	struct ImplicitCondition : public Condition
	{
		enum class Type : uint8_t { NONCONFORM_INTERFACE = MacroType::NONCONFORM_INTERFACE, 
										   HOMOGENEOUS_NEWMAN = MacroType::HOMOGENEOUS_NEWMAN};

		ImplicitCondition(Type type);
	};

	struct ConformCondition : public Condition
	{
		double thermalConductivity;
		ConformCondition(double thermalConductivity);
	};

	struct FunctionCondition : public Condition
	{
		enum class Type : uint8_t { DIRICHLET = MacroType::DIRICHLET_F, NEWMAN = MacroType::NEWMAN_F};
		FunctionCondition(Type type);

		template<Type conditionType>
		inline static double compute(const Coordinates& point)
		{
			return 0.0;
		}

		template<>
		inline static double compute<Type::DIRICHLET>(const Coordinates& point)
		{
			//return 0.0;
			return point.x;
			//return (point.x) * T;
		}

		template<>
		inline static double compute<Type::NEWMAN>(const Coordinates& point)
		{
			return 0.0;
		}

		template<Type conditionType, const uint8_t N_POINTS>
		static void compute(const LocalCoordinates2D* localPointIt,
						    const double transpJacobian[LocalCoordinates2D::COUNT * Coordinates::COUNT],
						    const Coordinates triangleStartPoint,
						    double* conditionValueIt)
		{
			Coordinates point;
			for (uint8_t i = 0; i < N_POINTS; ++i)
			{
				point = triangleStartPoint;
				point.x += transpJacobian[0] * localPointIt->u;
				point.y += transpJacobian[1] * localPointIt->u;
				point.z += transpJacobian[2] * localPointIt->u;

				point.x += transpJacobian[3] * localPointIt->v;
				point.y += transpJacobian[4] * localPointIt->v;
				point.z += transpJacobian[5] * localPointIt->v;

				*conditionValueIt = compute<conditionType>(point);
				++conditionValueIt;
				++localPointIt;
			}
		}
	};

	struct ValueCondition : public Condition
	{
		enum class Type : uint8_t { DIRICHLET = MacroType::DIRICHLET_V, NEWMAN = MacroType::NEWMAN_V, STEFAN = MacroType::STEFAN_V };
		double value;

		ValueCondition(Type type, double value);
	};

	static const FunctionCondition NewmanFunctionCondition;
	static const FunctionCondition DirichletFunctionCondition;
	static const ImplicitCondition HomogeneousNemanCondition;
	static const ImplicitCondition NonconformInterfaceCondition;

	//int tag;
	//bool isPlane;
	//bool isInterface;
	//unsigned int interfaceIndex;
	const Condition* condition;
	unsigned int regionsIndexes[2];
	unsigned char type;
};

