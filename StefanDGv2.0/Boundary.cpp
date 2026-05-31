#include "Boundary.h"

const Boundary::FunctionCondition Boundary::NewmanFunctionCondition = { FunctionCondition::Type::NEWMAN };
const Boundary::FunctionCondition Boundary::DirichletFunctionCondition = { FunctionCondition::Type::DIRICHLET };
const Boundary::ImplicitCondition Boundary::HomogeneousNemanCondition = { ImplicitCondition::Type::HOMOGENEOUS_NEWMAN };
const Boundary::ImplicitCondition Boundary::NonconformInterfaceCondition = { ImplicitCondition::Type::NONCONFORM_INTERFACE };

Boundary::Condition::Condition(MacroType macroType) : macroType{ macroType } 
{
}

Boundary::ImplicitCondition::ImplicitCondition(Type type) : Condition((MacroType)type)
{
}

Boundary::FunctionCondition::FunctionCondition(Type type) : Condition((MacroType)type)
{
}

Boundary::ValueCondition::ValueCondition(Type type, double value) : Condition((MacroType)type), value{value}
{
}

Boundary::ConformCondition::ConformCondition(double thermalConductivity) : Condition(MacroType::CONFORM_INTERFACE), thermalConductivity{ thermalConductivity }
{
}
