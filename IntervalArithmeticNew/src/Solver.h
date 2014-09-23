/*
 * Solver.h
 *
 *  Created on: 25-01-2014
 *      Author: thof
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "Utils.h"
#include "IntervalArithmetic.h"
#include "BoundaryConditions.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <time.h>
#include <boost/lexical_cast.hpp>

using namespace std;

namespace interval_arithmetic
{

template<typename T>
class Solver
{
protected:
	bool _initparams;
	bool _estimateMN;
	Parameters<T> params;
	Interval<T>* X;
	BoundaryConditions<T>* bc;
	long double** u;
	long double maxM;
	long double maxN;
	vector<long double> vecConstM;
	vector<long double> vecConstN;
public:
	Solver();
	virtual ~Solver();
	void SetParameters(Parameters<T>& p);
	void WriteFPResultsToFile();
	void WriteIntervalResultsToFile();
	void WriteConstMResults();
	void WriteResults();
	void Execute();
	int ConstMExperiment();
	int SolveInterval();
	bool SetEstimateMN(bool b);
	bool GetEstimateMN();
	long double GetMaxM();
	long double GetMaxN();
	virtual int SolveFP();
	virtual int SolvePIA();
	virtual int SolveDIA();
	virtual int SetExample(int eid);
};

} /* namespace intervalarth */
#endif /* SOLVER_H_ */
