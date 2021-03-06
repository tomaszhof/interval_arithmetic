/*
 * GPDESolver.h
 *
 *  Created on: 25-01-2014
 *      Author: thof
 */

#ifndef GPE5CSOLVER_H_
#define GPE5CSOLVER_H_

#include "Solver.h"
#include "THashMap.h"
#include "ExampleGPE01.h"
#include "ExampleGPE02.h"
#include "ExampleGPE03.h"
#include "ExampleGPE04.h"
#include "ExampleGPE05.h"
#include "ExampleGPE06.h"
#include "ExampleGPE08.h"
#include "ExampleGPE11.h"

using namespace interval_arithmetic;

namespace interval_arithmetic
{

template<typename T>
class GPE5CSolver: public Solver<T>
{
using Solver<T>::bc;
using Solver<T>::params;
using Solver<T>::u;
//using Solver<T>::X;
using Solver<T>::maxP;
using Solver<T>::maxQ;
using Solver<T>::maxR;
using Solver<T>::maxS;
using Solver<T>::maxT;
using Solver<T>::_initparams;

using Solver<T>::ih;
using Solver<T>::ik;

const Interval<T> im1 = {-1.0L, -1.0L};
const Interval<T> i1 = {1.0L, 1.0L};
const Interval<T> i2 = {2.0L, 2.0L};
const Interval<T> i6 = {6.0L, 6.0L};
const Interval<T> i12 = {12.0L, 12.0L};

T alphax(T xi, T yj);
T alphay(T xi, T yj);
T betax(T xi, T yj);
T betay(T xi, T yj);

Interval<T> IAlphaX(Interval<T> xi, Interval<T> yj);
Interval<T> IAlphaY(Interval<T> xi, Interval<T> yj);
Interval<T> IBetaX(Interval<T> xi, Interval<T> yj);
Interval<T> IBetaY(Interval<T> xi, Interval<T> yj);

public:
	GPE5CSolver();
	virtual ~GPE5CSolver();
	int SolveFP();
	int SolvePIA();
	int SolveDIA();
	int SetExample(int eid);
	void WriteFPResultsToCsv();
	void WriteIntervalResultsToCsv();
};

} /* namespace intervalarth */

template<typename T>
inline void interval_arithmetic::GPE5CSolver<T>::WriteFPResultsToCsv() {
}

template<typename T>
inline void interval_arithmetic::GPE5CSolver<T>::WriteIntervalResultsToCsv() {
}

#endif /* GPE5CSOLVER_H_ */
