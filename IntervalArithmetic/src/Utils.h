/*
 * Utils.h
 *
 *  Created on: 26 kwi 2014
 *      Author: thof
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <string>

using namespace std;

namespace intervalarth
{
enum IAMode {DINT_MODE, PINT_MODE};
enum ExperimentMode {CONST_M_EXP, CLASSICAL_EXP, INTERVAL_EXP};
enum Solvers {GPDE_SOLVER};

struct Parameters
{
public:
	Solvers selected_solver;
	IAMode ia_mode;
	ExperimentMode exp_mode;
	int m;
	int n;
	int example_id;
	long double alpha;
	long double beta;
	long double gamma;
	long double delta;
	long double eps;
	string file_name;
};
}



#endif /* UTILS_H_ */
