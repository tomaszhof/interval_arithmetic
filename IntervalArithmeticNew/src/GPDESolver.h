/*
 * GPDESolver.h
 *
 *  Created on: 25-01-2014
 *      Author: thof
 */

#ifndef GPDESOLVER_H_
#define GPDESOLVER_H_

#include "Solver.h"
#include "THashMap.h"
#include "BoundaryConditions.h"
//#include "Example01.h"
//#include "Example02.h"
#include "Example03.h"
//#include "Example04.h"

using namespace interval_arithmetic;

namespace interval_arithmetic
{

template<typename T>
class GPDESolver: public Solver<T>
{
using Solver<T>::bc;
using Solver<T>::params;
using Solver<T>::u;
//using Solver<T>::X;
using Solver<T>::maxM;
using Solver<T>::maxN;
using Solver<T>::_initparams;

public:
	GPDESolver();
	virtual ~GPDESolver();
	int SolveFP();
	int SolvePIA();
	int SolveDIA();
	int SetExample(int eid);
};

} /* namespace intervalarth */
#endif /* GPDESOLVER_H_ */
