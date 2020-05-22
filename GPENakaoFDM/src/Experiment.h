/*
 * Experiment.h
 *
 *  Created on: Jan 5, 2013
 *      Author: tomaszhof
 */

#ifndef EXPERIMENT_H_
#define EXPERIMENT_H_

#include "Utils.h"
#include "BoundaryConditions.h"
#include "ExampleGPE01.h"
#include "ExampleGPE02.h"
#include "ExampleGPE03.h"
#include "ExampleGPE04.h"
#include "ExampleGPE05.h"
#include "ExampleGPE06.h"
//#include "Example09.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <float.h>
#include <boost/lexical_cast.hpp>
#include "ConfigReader.h"
#include "GPDESolver.h"
#include "GPESolver3C.h"


using namespace std;
using namespace interval_arithmetic;

namespace interval_arithmetic
{

template<typename T>
class Experiment
{
private:
	BoundaryConditions<T>* _example;
	Parameters<long double> parameters;
	Solver<T>* solver;
	bool _param_initialized;
	bool _solver_initialized;
public:
	Experiment();
	Experiment(int ac, char *av[]);
	virtual ~Experiment();
	void SetExample(int eid, int arth_mode);
	void SetSolver(Parameters<long double> p);
	void SetExampleForSolver(int eid);
	void Initialize();
	void SetParameters(Parameters<long double> p);
	void Execute();
	static void SetMode(IAMode m);
};


}

#endif /* EXPERIMENT_H_ */
