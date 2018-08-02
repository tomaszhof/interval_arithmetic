/*
 * NakaoExperiment2D.h
 *
 *  Created on: Jul 21, 2018
 *      Author: numeric
 */

#ifndef NAKAOEXPERIMENT2D_H_
#define NAKAOEXPERIMENT2D_H_

#include "Interval.h"
#include <string>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <boost/lexical_cast.hpp>
#include <stdlib.h>


using namespace interval_arithmetic;
using namespace std;

const long double PI = 3.141592653589793238L;

class NakaoExperiment2D {

private:
	int error, i, it, j, j1, jh, k, kh, l, l1, l2,
	lh, n, n1, n2, p, q, rh;
	long double a, abs_alpha, abs_alphaij, b, b_dash,
	beta, c, d, delta, epsilon, exact, h,
	max, norm_u, norm_uij, s, sph;
	Interval<long double> ia, ia1, ialpha, ib, ib1, ibeta, ic, ic1,
	id, id1, ih, imax, interval_s, ipi, iz;
	int *r;
	long double *a1, *b1, *x;
	Interval<long double> *interval_a1, *interval_b1, *interval_x;
	long double **alpha_k, **alpha_km1, **u;
	Interval<long double> **iu_k, **iu_km1;
	char answer, output;
	bool alpha_OK, finish, u_OK;
	string file_name, left, right, st, tmpstr;
	fstream results;

	const Interval<long double> i0 = {0.0, 0.0};
	const Interval<long double> i1 = {1.0, 1.0};
	const Interval<long double> i6 = {6.0, 6.0};
	const Interval<long double> im11 = {-1.0, 1.0};


public:
	NakaoExperiment2D();
	virtual ~NakaoExperiment2D();
	void execute();
};

#endif /* NAKAOEXPERIMENT2D_H_ */
