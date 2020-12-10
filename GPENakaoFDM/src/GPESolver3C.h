/*
 * GPESolver3C.h
 *
 *  Created on: 01.05.2019
 *      Author: numeric
 */

#ifndef GPESOLVER3C_H_
#define GPESOLVER3C_H_

#include "Solver.h"
#include "THashMap.h"
#include "ExampleGPE01.h"
#include "ExampleGPE02.h"
#include "ExampleGPE03.h"
#include "ExampleGPE04.h"
#include "ExampleGPE05.h"
#include "ExampleGPE06.h"
#include "ExampleGPE07.h"
#include "ExampleGPE08.h"
#include "ExampleGPE11.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class GPESolver3C: public Solver<T> {
	using Solver<T>::bc;
	using Solver<T>::params;
	using Solver<T>::u;
	//using Solver<T>::X;
	using Solver<T>::maxP;
	using Solver<T>::maxQ;
	using Solver<T>::maxR;
	using Solver<T>::maxS;
	using Solver<T>::_initparams;

	using Solver<T>::ih;
	using Solver<T>::ik;

	const Interval<T> im2 = { -2.0L, -2.0L };
	const Interval<T> im1 = { -1.0L, -1.0L };
	const Interval<T> i1 = { 1.0L, 1.0L };
	const Interval<T> i2 = { 2.0L, 2.0L };
	const Interval<T> i6 = { 6.0L, 6.0L };
	const Interval<T> i12 = { 12.0L, 12.0L };

	T alphax(T xi, T yj);
	T alphay(T xi, T yj);
	T betax(T xi, T yj);
	T betay(T xi, T yj);
	T gammaxy(T xi, T yj);

	Interval<T> IAlphaX(Interval<T> xi, Interval<T> yj);
	Interval<T> IAlphaY(Interval<T> xi, Interval<T> yj);
	Interval<T> IBetaX(Interval<T> xi, Interval<T> yj);
	Interval<T> IBetaY(Interval<T> xi, Interval<T> yj);
	Interval<T> IGammaXY(Interval<T> xi, Interval<T> yj);

	Interval<T> DIAlphaX(Interval<T> xi, Interval<T> yj);
	Interval<T> DIAlphaY(Interval<T> xi, Interval<T> yj);
	Interval<T> DIBetaX(Interval<T> xi, Interval<T> yj);
	Interval<T> DIBetaY(Interval<T> xi, Interval<T> yj);
	Interval<T> DIGammaXY(Interval<T> xi, Interval<T> yj);

public:
	GPESolver3C();
	virtual ~GPESolver3C();
	int SolveFP();
	int SolvePIA();
	int SolveDIA();
	int SetExample(int eid);
	void WriteFPResultsToCsv();
	void WriteIntervalResultsToCsv();
};

}

#endif /* GPESOLVER3C_H_ */
