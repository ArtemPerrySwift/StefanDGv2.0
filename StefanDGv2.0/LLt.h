#pragma once

class LLt
{
public:
	LLt(double matrix[], const size_t matrixSize);
	~LLt();

	void solve(const double b[], double x[]) const;
private:

	void directStep(const double b[], double x[]) const;
	void reverseStep(double x[]) const;


	double* _elements;
	double* _di;
	double* _l;
	double** _il;
	size_t _matrixSize;
};

