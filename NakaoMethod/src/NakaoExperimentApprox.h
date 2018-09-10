/*
0 * NakaoExperimentApprox.h
 *
 *  Created on: Feb 25, 2018
 *      Author: numeric
 */

#ifndef NAKAOEXPERIMENTAPPROX_H_
#define NAKAOEXPERIMENTAPPROX_H_
#include "Interval.h"
#include <string>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <boost/lexical_cast.hpp>

using namespace interval_arithmetic;
using namespace std;



class NakaoExperimentApprox {
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


	std::vector<long double> GL_X4 = {-0.906180, -0.538469, 0., 0.538469, 0.906180};
	std::vector<long double> GL_A4 = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927};

	Interval<long double> i4 = {4.0, 4.0};
	Interval<long double> i6 = {6.0, 6.0};
	Interval<long double> im11 = {-1.0, 1.0};

	friend long double phi(int i, long double x);
	friend long double phi_p(int i, long double x);

	friend long double f(long double x);

	void initialize_vx(long double a, long double b, int n);

public:
	NakaoExperimentApprox();
	virtual ~NakaoExperimentApprox();
	void initialize();
	void execute();
};

#endif /* NAKAOEXPERIMENTAPPROX_H_ */
