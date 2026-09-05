#pragma once
#include <cinttypes>

namespace array_functions
{
	template<uint8_t ARRAY_SIZE>
	inline void multiply(const double multiplier, double* elementIt)
	{
		for (uint8_t i = 0; i < ARRAY_SIZE; ++i)
		{
			*elementIt *= multiplier;
			++elementIt;
		}
	}

	template<uint8_t ARRAY_SIZE>
	inline void add(const double* addendumIt, double* resultIt)
	{
		for (uint8_t i = 0; i < ARRAY_SIZE; ++i)
		{
			*resultIt += *addendumIt;

			++addendumIt;
			++resultIt;
		}
	}
}