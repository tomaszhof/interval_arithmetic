/*
 * NakaoExperiment2DCxy.h
 *
 *  Created on: Jul 21, 2018
 *      Author: numeric
 */

#ifndef NAKAOEXPERIMENT2DCXY_H_
#define NAKAOEXPERIMENT2DCXY_H_

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

class NakaoExperiment2DCxy {

private:
	int error, i, it, j, j1, jh, k, kh, l, l1, l2,
	lh, n, n1, n2, p, q, rh;
	long double a, abs_alpha, abs_alphaij, b,
	beta, c, d, delta, epsilon, exact, h,
	max, norm_u, norm_uij, s, sph, ij_over_12, i_over_24, j_over_24, twenty_h2;
	Interval<long double> ia, ia1, ialpha, ib, ib1, ibeta, ibeta1, ibeta2, ic, ic1,
	id, id1, ih,ih2, imax, interval_s, isph, ipi, ipi12, iz, i20h4,
	i_plus_j_ih, i_minus_j_ih, interval_i, interval_j, interval_ij, ipih, i2pih,iz1,
	interval_ij_over_12, interval_i_over_24, interval_j_over_24;
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
	const Interval<long double> i2 = {2.0, 2.0};
	const Interval<long double> i3 = {3.0, 3.0};
	const Interval<long double> i4 = {4.0, 4.0};
	const Interval<long double> i5 = {5.0, 5.0};
	const Interval<long double> i6 = {6.0, 6.0};
	const Interval<long double> i7 = {7.0, 7.0};
	const Interval<long double> i12 = {12.0, 12.0};
	const Interval<long double> i15 = {15.0, 15.0};
	const Interval<long double> i20 = {20.0, 20.0};
	const Interval<long double> i24 = {24.0, 24.0};
	const Interval<long double> i36 = {36.0, 36.0};
	const Interval<long double> i45 = {45.0, 45.0};
	const Interval<long double> i120 = {120.0, 120.0};
	const Interval<long double> i180 = {180.0, 180.0};
	const Interval<long double> i360 = {360.0, 360.0};
	const Interval<long double> im1 = {-1.0, -1.0};
	const Interval<long double> im11 = {-1.0, 1.0};


public:
	NakaoExperiment2DCxy();
	virtual ~NakaoExperiment2DCxy();
	void execute();
};

#endif /* NAKAOEXPERIMENT2DCXY_H_ */
