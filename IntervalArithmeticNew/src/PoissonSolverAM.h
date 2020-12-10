/*
 * PoissonSolverAM.h
 *
 *  Created on: 13 lis 2017
 *      Author: tomhof
 */

#ifndef POISSONSOLVERAM_H_
#define POISSONSOLVERAM_H_

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
#include "Example10.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class PoissonSolverAM: public Solver<T> {
	using Solver<T>::bc;
	using Solver<T>::params;
	using Solver<T>::u;
	//using Solver<T>::X;
	using Solver<T>::maxM;
	using Solver<T>::maxN;
	using Solver<T>::_initparams;

public:
	PoissonSolverAM();
	virtual ~PoissonSolverAM();

	int SolveFP();
	int SolvePIA();
	int SolveDIA();
	int SetExample(int eid);
	void WriteFPResultsToCsv();
	void WriteIntervalResultsToCsv();
};

}

#endif /* POISSONSOLVERAM_H_ */
