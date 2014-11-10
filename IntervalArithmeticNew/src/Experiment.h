/*
 * Experiment.h
 *
 *  Created on: Jan 5, 2013
 *      Author: tomaszhof
 */

#ifndef EXPERIMENT_H_
#define EXPERIMENT_H_

#include "Interval.h"
//#include "DiffPoisson.h"
#include "BoundaryConditions.h"
//#include "Example01.h"
//#include "Example02.h"
#include "Example03.h"
#include "Example04.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <float.h>
#include <boost/lexical_cast.hpp>
#include "ConfigReader.h"
#include "Solver.h"
#include "GPDESolver.h"
#include "Utils.h"

using namespace std;
using namespace interval_arithmetic;

namespace interval_arithmetic
{

template<typename T>
class Experiment
{
private:
	BoundaryConditions<T>* _example;
	Parameters<T> parameters;
	Solver<T>* solver;
	bool _param_initialized;
	bool _solver_initialized;
public:
	Experiment();
	Experiment(int ac, char *av[]);
	virtual ~Experiment();
	void SetExample(int eid, int arth_mode);
	void SetSolver(Parameters<T> p);
	void SetExampleForSolver(int eid);
	void Initialize();
	void SetParameters(Parameters<T> p);
	void Execute();
	static void SetMode(IAMode m);
};


}

#endif /* EXPERIMENT_H_ */
