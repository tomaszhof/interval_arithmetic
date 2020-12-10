/*
 * GSLIntegrator.cpp
 *
 *  Created on: 6 cze 2020
 *      Author: tomhof
 */

#include "GSLIntegrator.h"

GSLIntegrator::GSLIntegrator() {
	// TODO Auto-generated constructor stub

}

GSLIntegrator::~GSLIntegrator() {
	// TODO Auto-generated destructor stub
}

void GSLIntegrator::display_results(char *title, double result, double error){

}
double GSLIntegrator::integrate(double (*f)(double * x_array, size_t dim, void * params), double a1, double b1, double a2, double b2, double* int_err){
	double res, err;

		double xl[2] = { a1, a2 };
		double xu[2] = { b1, b2 };

		const gsl_rng_type *T;
		gsl_rng *r;

		gsl_monte_function G = { f, 2, 0 };

		gsl_rng_env_setup();

		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
//
		{
			gsl_monte_plain_state *s = gsl_monte_plain_alloc(2);
			gsl_monte_plain_integrate(&G, xl, xu, 2, calls, r, s, &res, &err);
			gsl_monte_plain_free(s);

//			display_results("plain", res, err);
		}

//		{
//			gsl_monte_miser_state *s = gsl_monte_miser_alloc(2);
//			gsl_monte_miser_integrate(&G, xl, xu, 2, calls, r, s, &res, &err);
//			gsl_monte_miser_free(s);
//			display_results("miser", res, err);
//			*int_err = err;
//		}


//		  {
//		    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);
//
//		    gsl_monte_vegas_integrate (&G, xl, xu, 2, 10000, r, s,
//		                               &res, &err);
//
////		    display_results ("vegas warm-up", res, err);
////		    printf ("converging...\n");
//
//		    double chisq = 0.0;
//		    do
//		      {
//		        gsl_monte_vegas_integrate (&G, xl, xu, 2, 100000, r, s,
//		                                   &res, &err);
//		        chisq = gsl_monte_vegas_chisq (s);
////		        printf ("result = % .6f sigma = % .6f "
////		                "chisq/dof = %.1f\n", res, err, chisq);
//		      }
//		    while ((fabs(chisq) > 0.0) && (fabs(chisq - 1.0) > 0.5));
//
////		    display_results ("vegas final", res, err);
//
//		    gsl_monte_vegas_free (s);
//		  }

		  gsl_rng_free (r);

		return res;
}
