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

namespace interval_arithmetic
{
enum IAMode {DINT_MODE, PINT_MODE};
enum ExperimentMode {CONST_M_EXP, CLASSICAL_EXP, INTERVAL_EXP};
enum Solvers {GPDE_SOLVER};

template<typename T>
struct Parameters
{
public:
	Solvers selected_solver;
	IAMode ia_mode;
	ExperimentMode exp_mode;
	int m;
	int n;
	int example_id;
	T alpha;
	T beta;
	T gamma;
	T delta;
	T eps;
	string file_name;
	bool print_csv;
};
}



#endif /* UTILS_H_ */
