#include "LLt.h"
#include <cmath>

LLt::LLt(double matrix[], const size_t matrixSize) : _matrixSize{ matrixSize }
{
	_elements = new double[matrixSize * (matrixSize - 1) / 2 + matrixSize];
	_di = _elements;
	_l = _di + matrixSize;
	_il = new double*[matrixSize];

	double** _ilIt = _il;
	double* lElementIt = _l;
	double* diIt = _di;

	size_t nElementsUntillRowEnd = matrixSize;
	for (size_t i = 0; i < matrixSize; ++i)
	{
		diIt = _di;
		*_ilIt = lElementIt;

		double* iRow = lElementIt;
		double** jColumnIt = _il;

		double diSum = 0;
		for (size_t j = 0; j < i; ++j)
		{
			double sum = 0;
			double* ikIt = iRow;
			double* jkIt = *jColumnIt;
			for (size_t k = 0; k < j; ++k)
			{
				sum += (*ikIt) * (*jkIt);
				++ikIt;
				++jkIt;
			}

			double element = (*matrix - sum) / *diIt;
			*lElementIt = element;
			diSum += element * element;

			++matrix;
			++lElementIt;
			++diIt;
			++jColumnIt;
		}

		*diIt = sqrt(*matrix - diSum);
		matrix += nElementsUntillRowEnd;
		--nElementsUntillRowEnd;
		++_ilIt;
	}
}

LLt::~LLt()
{
	delete[] _elements;
	delete[] _il;
}

void LLt::solve(const double b[], double x[]) const
{
	directStep(b, x);
	reverseStep(x);
}

void LLt::directStep(const double b[], double x[]) const
{
	double* diIt = _di;
	double* iyIt = x;
	const double* lIt = _l;

	for (size_t i = 0; i < _matrixSize; ++i)
	{
		double sum = 0;
		const double* jyIt = x;
		for (size_t j = 0; j < i; ++j)
		{
			sum += (*lIt) * (*jyIt);
			++lIt;
			++jyIt;
		}

		*iyIt = (*b - sum) / *diIt;

		++diIt;
		++b;
		++iyIt;
	}
}

void LLt::reverseStep(double x[]) const
{
	const double* const* ilRevIt = _il + _matrixSize;
	const double* diRevIt = _di + _matrixSize;
	double* yRevIt = x + _matrixSize;
	double* xRevIt = yRevIt;
	size_t diIndex = _matrixSize;
	for (size_t i = 0; i < _matrixSize; ++i)
	{
		--ilRevIt;
		--diRevIt;
		--diIndex;
		--xRevIt;
		--yRevIt;

		double xElement = (*yRevIt) / (*diRevIt);
		*xRevIt = xElement;
		const double* columnIt = *ilRevIt;

		double* yIt = x;
		for (size_t j = 0; j < diIndex; ++j)
		{
			*yIt -= (*columnIt) * xElement;

			++columnIt;
			++yIt;
		}

	}
}
