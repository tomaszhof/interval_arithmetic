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
#include "ExampleGPE01.h"
#include "ExampleGPE02.h"

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

T alphax(T xi, T yj);
T alphay(T xi, T yj);
T betax(T xi, T yj);
T betay(T xi, T yj);

public:
	GPDESolver();
	virtual ~GPDESolver();
	int SolveFP();
	int SolvePIA();
	int SolveDIA();
	int SetExample(int eid);
	void WriteFPResultsToCsv();
	void WriteIntervalResultsToCsv();
};

} /* namespace intervalarth */

template<typename T>
inline void interval_arithmetic::GPDESolver<T>::WriteFPResultsToCsv() {
}

template<typename T>
inline void interval_arithmetic::GPDESolver<T>::WriteIntervalResultsToCsv() {
}

#endif /* GPDESOLVER_H_ */
