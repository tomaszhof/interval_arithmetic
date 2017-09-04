/*
 * PoissonSolver4OrderAM.h
 *
 *  Created on: 4 wrz 2017
 *      Author: tomhof
 */

#ifndef POISSONSOLVER4ORDERAM_H_
#define POISSONSOLVER4ORDERAM_H_

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
class PoissonSolver4OrderAM : public Solver<T>{
using Solver<T>::bc;
using Solver<T>::params;
using Solver<T>::u;
//using Solver<T>::X;
using Solver<T>::maxM;
using Solver<T>::maxN;
using Solver<T>::_initparams;

public:
	PoissonSolver4OrderAM();
	virtual ~PoissonSolver4OrderAM();

	int SolveFP();
	int SolvePIA();
	int SolveDIA();
	int SetExample(int eid);
	void WriteFPResultsToCsv();
	void WriteIntervalResultsToCsv();
};
template<typename T>
inline void interval_arithmetic::PoissonSolver4OrderAM<T>::WriteFPResultsToCsv() {
}

template<typename T>
inline void interval_arithmetic::PoissonSolver4OrderAM<T>::WriteIntervalResultsToCsv() {
}
}

#endif /* POISSONSOLVER4ORDERAM_H_ */
