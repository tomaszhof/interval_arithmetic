/*
0 * NakaoExperiment.h
 *
 *  Created on: Feb 25, 2018
 *      Author: numeric
 */

#ifndef NAKAOEXPERIMENT_H_
#define NAKAOEXPERIMENT_H_
#include "Interval.h"
#include <string>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <boost/lexical_cast.hpp>

using namespace interval_arithmetic;
using namespace std;

class NakaoExperiment {
private:
	int error, i, k, n;
	long double a, abs_alpha, abs_alphai, alpha, b,
	beta, c, d, delta, epsilon, exact,
	h, norm_u, norm_ui;
	Interval<long double> ia, ia1, ib, ib1, ibeta, ic, ic1,
	id, id1, ialpha, ih, iPi, iSin;
	long double* alpha_k, *alpha_km1, *u, *y, *z;
	Interval<long double>* iu_k, *iu_km1, *iy, *iz;
	char answer, output;
	bool alpha_OK, finish, u_OK;
	string file_name, left, right, s;
	fstream results;

	Interval<long double> i4 = {4.0, 4.0};
	Interval<long double> i6 = {6.0, 6.0};
	Interval<long double> im11 = {-1.0, 1.0};


public:
	NakaoExperiment();
	virtual ~NakaoExperiment();
	void execute();
};

#endif /* NAKAOEXPERIMENT_H_ */
