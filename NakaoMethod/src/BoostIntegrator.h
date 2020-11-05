/*
 * BoostIntegrator.h
 *
 *  Created on: 2 lis 2020
 *      Author: tomhof
 */

#ifndef SRC_BOOSTINTEGRATOR_H_
#define SRC_BOOSTINTEGRATOR_H_

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <iostream>

using namespace std;

class BoostIntegrator {
public:
	BoostIntegrator();
	virtual ~BoostIntegrator();
	long double integrate(long double (*f)(long double x, long double y), long double a1, long double b1, long double a2, long double b2, long double* int_err);
};

#endif /* SRC_BOOSTINTEGRATOR_H_ */
