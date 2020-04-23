/*
 * Solver.h
 *
 *  Created on: 25-01-2014
 *      Author: thof
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "Utils.h"
#include "BoundaryConditions.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <time.h>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>


using namespace std;

namespace fs = boost::filesystem ;

namespace interval_arithmetic
{

template<typename T>
class Solver
{
protected:
	bool _initparams;
	bool _estimateMN;
	Parameters<long double> params;
	BoundaryConditions<T>* bc;
	long double** u;
	long double maxP;
	long double maxQ;
	long double maxM;
	long double maxN;
	vector<long double> vecConstP;
	vector<long double> vecConstQ;
public:
	Interval<T>* X;
	Solver();
	virtual ~Solver();
	void SetParameters(Parameters<long double>& p);
	void WriteFPResultsToFile();
	void WriteIntervalResultsToFile();
	void WriteFPResultsToCsv();
	void WriteIntervalResultsToCsv();
	void WriteConstMResults();
	void WriteResults();
	void InitializeX(int m, int n);
	void Execute();
	int ConstMExperiment();
	int SolveInterval();
	bool SetEstimateMN(bool b);
	bool GetEstimateMN();
	long double GetMaxP();
	long double GetMaxQ();
	virtual int SolveFP();
	virtual int SolvePIA();
	virtual int SolveDIA();
	virtual int SetExample(int eid);
};

} /* namespace intervalarth */



#endif /* SOLVER_H_ */
