/*
 * BoostIntegrator.cpp
 *
 *  Created on: 2 lis 2020
 *      Author: tomhof
 */

#include "BoostIntegrator.h"

BoostIntegrator::BoostIntegrator() {
	// TODO Auto-generated constructor stub

}

BoostIntegrator::~BoostIntegrator() {
	// TODO Auto-generated destructor stub
}

long double BoostIntegrator::integrate(long double (*f)(long double x, long double y), long double a1, long double b1, long double a2, long double b2, long double* int_err) {
	using namespace boost::math::quadrature;

	    auto f1 = [&](long double x) {
	        auto g = [&](long double y) {
	            return f(x, y);
	        };
	        return gauss_kronrod<long double, 61>::integrate(g, a2, b2, 3);
	    };

	    long double Q = gauss_kronrod<long double, 61>::integrate(f1, a1, b1, 3, 1e-15, int_err);
	    //std::cout << Q << ", error estimated at " << *int_err <<std::endl;
	    return Q;
}
