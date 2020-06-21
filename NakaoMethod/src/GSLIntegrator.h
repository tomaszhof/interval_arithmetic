/*
 * GSLIntegrator.h
 *
 *  Created on: 6 cze 2020
 *      Author: tomhof
 */

#ifndef SRC_GSLINTEGRATOR_H_
#define SRC_GSLINTEGRATOR_H_

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

class GSLIntegrator {
public:
	GSLIntegrator();
	virtual ~GSLIntegrator();
    double integrate(double (*f)(double * x_array, size_t dim, void * params), double a1, double b1, double a2, double b2, double* int_err);
    void display_results(char *title, double result, double error);
};

#endif /* SRC_GSLINTEGRATOR_H_ */
