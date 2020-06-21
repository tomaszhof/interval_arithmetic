/*
 * NakaoExperiment2DApprox.h
 *
 *  Created on: Jul 21, 2018
 *      Author: numeric
 */

#ifndef NAKAOEXPERIMENT2D_APPROX_H_
#define NAKAOEXPERIMENT2D_APPROX_H_

#include "Interval.h"
#include "GSLIntegrator.h"
#include <string>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <boost/lexical_cast.hpp>
#include <stdlib.h>


using namespace interval_arithmetic;
using namespace std;

const long double PI = 3.141592653589793238L;

class NakaoExperiment2DApprox {

private:
	GSLIntegrator* integrator;
	static int i, j;
	static long double h;
	int error, it, j1, jh, k, kh, l, l1, l2,
	lh, n, n1, n2, p, q, rh;
	long double a, abs_alpha, abs_alphaij, b, b_dash,
	beta, c, d, delta, epsilon, exact,
	max, norm_u, norm_uij, s, sph;
	Interval<long double> ia, ia1, ialpha, ib, ib1, ibeta, ic, ic1,
	id, id1, ih,ih2, imax, interval_s, ipi, ipi12, iz, iierr;
	int *r;
	long double *a1, *b1, *x;
	Interval<long double> *interval_a1, *interval_b1, *interval_x;
	long double **alpha_k, **alpha_km1, **u;
	Interval<long double> **iu_k, **iu_km1;
	char answer, output;
	bool alpha_OK, finish, u_OK;
	string file_name, left, right, st, tmpstr;
	fstream results;

	const Interval<long double> im1 = {-1.0, -1.0};
	const Interval<long double> i0 = {0.0, 0.0};
	const Interval<long double> i1 = {1.0, 1.0};
	const Interval<long double> i2 = {2.0, 2.0};
	const Interval<long double> i4 = {4.0, 4.0};
	const Interval<long double> i6 = {6.0, 6.0};
	const Interval<long double> i12 = {12.0, 12.0};
	const Interval<long double> im11 = {-1.0, 1.0};


public:
	NakaoExperiment2DApprox();
	virtual ~NakaoExperiment2DApprox();
//	static long double phi(int i, int j, long double x, long double y);
//	static long double f(long double x, long double y);
//	static double g_f_phi(double *k, size_t dim, void *params);
//	double g_int_c_ij1(double *k, size_t dim, void *params);
//	double g_int_c_ij2(double *k, size_t dim, void *params);
//	double g_int_c_ij3(double *k, size_t dim, void *params);
//	double g_int_c_ij4(double *k, size_t dim, void *params);
//	double g_int_c_ij5(double *k, size_t dim, void *params);
//	double g_int_c_ij6(double *k, size_t dim, void *params);
//	double g_int_c_ij7(double *k, size_t dim, void *params);
	void execute();

	static long double f(long double x, long double y){
		//(1.0-2.0*M_PI)*
		return sin(M_PI*x)*sin(M_PI*y);
	}

	static long double fe(long double x, long double y){
			//
			return (1.0-2.0*M_PI)*sin(M_PI*x)*sin(M_PI*y);
	}

	static long double ce(long double x, long double y){
			return M_PI;
	}

	static long double phi(int i, int j, long double x, long double y){
			if ((i*h <= x) && (x <= (i+1.0)*h) && (j*h <= y) && (y <= (i+j+1.0)*h-x)){
				return 1.0 + i + j - (1.0/h)*(x+y);
			}

			if (((i-1.0)*h <= x) && (x <= i*h) && ((i+j)*h-x <= y) && (y <= (j+1.0)*h)){
				return 1.0 + j - y/h;
			}

			if (((i-1.0)*h <= x) && (x <= i*h) && (j*h <= y ) && (y <=(i+j)*h -x)){
				return 1.0 - i + x/h;
			}

			if (((i-1.0)*h <= x) && (x <= i*h) && ((i + j - 1.0)*h - x <= y) && (y <= j*h)){
				return 1.0 - i - j + (1.0/h)*(x+y);
			}

			if ((i*h <= x) && (x <= (i+1.0)*h) && ((j-1.0)*h <= y) && (y <= (i+j)*h-x)){
				return 1.0 - j + y/h;
			}

			if ((i*h <= x) && (x <= (i+1.0)*h) && ((i+j)*h-x <= y) && (y <= j*h)){
				return 1.0 + i - x/h;
			}
			return 0.0;
	}

	static double g_f_phi(double *k, size_t dim, void *params){
		double x = k[0];
		double y = k[1];

		return fe(x,y) * phi(i, j, x, y);

	}

	static double g_int_c_ij1(double *k, size_t dim, void *params){
		double x = k[0];
		double y = k[1];

		return ce(x,y)*phi(i, j, x, y)*phi(i, j, x, y);

	}

	static double g_int_c_ij2(double *k, size_t dim, void *params){
		double x = k[0];
		double y = k[1];

		return ce(x,y)*phi(i, j, x, y)*phi(i, j-1, x, y);

	}

	static double g_int_c_ij3(double *k, size_t dim, void *params){
		double x = k[0];
		double y = k[1];

		return ce(x,y)*phi(i, j, x, y)*phi(i, j+1, x, y);

	}

	static double g_int_c_ij4(double *k, size_t dim, void *params){
		double x = k[0];
		double y = k[1];

		return ce(x,y)*phi(i, j, x, y)*phi(i-1, j, x, y);

	}

	static double g_int_c_ij5(double *k, size_t dim, void *params){
		double x = k[0];
		double y = k[1];

		return ce(x,y)*phi(i, j, x, y)*phi(i+1, j, x, y);

	}

	static double g_int_c_ij6(double *k, size_t dim, void *params){
		double x = k[0];
		double y = k[1];

		return ce(x,y)*phi(i, j, x, y)*phi(i-1, j+1, x, y);

	}


	static double g_int_c_ij7(double *k, size_t dim, void *params){
		double x = k[0];
		double y = k[1];

		return ce(x,y)*phi(i, j, x, y)*phi(i+1, j-1, x, y);

	}

};


#endif /* NAKAOEXPERIMENT2D_APPROX_H_ */
