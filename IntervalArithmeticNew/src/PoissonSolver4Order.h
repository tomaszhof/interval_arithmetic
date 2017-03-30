/*
 * PoissonSolver.h
 *
 *  Created on: 25-01-2014
 *      Author: thof
 */

#ifndef POISSON4SOLVER_H_
#define POISSON4SOLVER_H_

#include "Solver.h"
#include "THashMap.h"
//#include "BoundaryConditions.h"
//#include "Example01.h"
#include "Example02.h"
#include "Example03.h"
#include "Example04.h"
#include "Example05.h"
#include "Example06.h"
#include "Example07.h"
#include "Example08.h"
#include "Example09.h"
#include <iostream>
#include <fstream>

using namespace interval_arithmetic;

namespace interval_arithmetic
{

template<typename T>
class PoissonSolver4Order: public Solver<T>
{
using Solver<T>::bc;
using Solver<T>::params;
using Solver<T>::u;
//using Solver<T>::X;
using Solver<T>::maxM;
using Solver<T>::maxN;
using Solver<T>::_initparams;

public:
	PoissonSolver4Order();
	virtual ~PoissonSolver4Order();
	int SolveFP();
	int SolvePIA();
	int SolveDIA();
	int SetExample(int eid);
	void WriteFPResultsToCsv();
	void WriteIntervalResultsToCsv();
};

} /* namespace intervalarth */

template<typename T>
inline void interval_arithmetic::PoissonSolver4Order<T>::WriteFPResultsToCsv() {
}

template<typename T>
inline void interval_arithmetic::PoissonSolver4Order<T>::WriteIntervalResultsToCsv() {
}

#endif /* POISSON4SOLVER_H_ */
