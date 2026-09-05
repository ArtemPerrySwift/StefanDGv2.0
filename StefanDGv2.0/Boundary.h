#pragma once
#include "Coordinates.h"
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

		static enum class MacroType: uint8_t { DIRICHLET_F, DIRICHLET_V, NEWMAN_F, NEWMAN_V, STEFAN_V, CONFORM_INTERFACE, NONCONFORM_INTERFACE, HOMOGENEOUS_NEWMAN};
		const MacroType macroType;

	private:
		Condition(MacroType macroType);
	};

	struct ImplicitCondition : public Condition
	{
		static enum class Type : uint8_t { NONCONFORM_INTERFACE = MacroType::NONCONFORM_INTERFACE, 
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
		static enum class Type : uint8_t { DIRICHLET = MacroType::DIRICHLET_F, NEWMAN = MacroType::NEWMAN_F};
		FunctionCondition(Type type);

		static void computeDirichletValues(const Coordinates points[], uint8_t nPoints, double values[])
		{
			for (uint8_t i = 0; i < nPoints; ++i)
			{
				*values = 0.0;

				++points;
				++values;
			}
		}

		static void computeNewmanValues(const Coordinates points[], uint8_t nPoints, double values[])
		{
			for (uint8_t i = 0; i < nPoints; ++i)
			{
				*values = 0.0;

				++points;
				++values;
			}
		}
	};

	struct ValueCondition : public Condition
	{
		static enum class Type : uint8_t { DIRICHLET = MacroType::DIRICHLET_V, NEWMAN = MacroType::NEWMAN_V, STEFAN = MacroType::STEFAN_V };
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

