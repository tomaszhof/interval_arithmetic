/*
 * NakaoExperiment2DApprox.h
 *
 *  Created on: Jun 10, 2020
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
#include <chrono>


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
	long double abs_alpha, abs_alphaij,
	beta, c, delta, epsilon, exact, width,
	max, norm_u, norm_uij, s, d, sph;
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

	long double aij, aij_l, aij_u, bij, bij_l, b_dash_ij, b_dash_ij_u, hpm2, mhpm2;
	long int duration = 0;

	const Interval<long double> im1 = {-1.0, -1.0};
	const Interval<long double> i0 = {0.0, 0.0};
	const Interval<long double> i1 = {1.0, 1.0};
	const Interval<long double> i2 = {2.0, 2.0};
	const Interval<long double> i4 = {4.0, 4.0};
	const Interval<long double> i6 = {6.0, 6.0};
	const Interval<long double> i12 = {12.0, 12.0};
	const Interval<long double> im11 = {-1.0, 1.0};

	size_t integrator_calls = 500000;


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


//-------------------------------------------------------------
//-------------------Nakao example 1---------------------------
//-------------------------------------------------------------
//	static long double u_exact(long double x, long double y){
//		return (1.0 / M_PI) * sin(M_PI * x) * sin(M_PI * y);
//	}
//
//	static long double fe(long double x, long double y){
//			//
//			return (1.0-2.0*M_PI)*sin(M_PI*x)*sin(M_PI*y);
//	}
//
//	static long double ce(long double x, long double y){
//			return  M_PI; //5.0/4.0*M_PI*M_PI;
//	}

////-------------------------------------------------------------
////-------------------Nakao example 2---------------------------
////-------------------------------------------------------------
//	static long double u_exact(long double x, long double y){
//		return 0.0;
//	}
//
//	static long double fe(long double x, long double y){
//			return (1.0-2.0*M_PI)*sin(M_PI*x)*sin(M_PI*y);
//	}
//
//	static long double ce(long double x, long double y){
//			return  20.0*x*y;
//	}

////-------------------------------------------------------------
////--------------------TH example-------------------------------
////-------------------------------------------------------------
	static long double u_exact(long double x, long double y){
		return x*cos(M_PI/2.0 * x)*sin(M_PI*y);
	}


	static long double fe(long double x, long double y){
			//
			return (-1.0*M_PI)*sin(M_PI/2.0*x)*sin(M_PI*y);
	}

	static long double ce(long double x, long double y){
			return  5.0/4.0*M_PI*M_PI;
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
