/*
 * Experiment.cpp
 *
 *  Created on: Jan 5, 2013
 *      Author: tomaszhof
 */

#include "Experiment.h"

namespace interval_arithmetic {

template<typename T>
Experiment<T>::Experiment() {
	_param_initialized = false;
	_solver_initialized = false;
	_example = NULL;
	solver = NULL;
}

template<typename T>
Experiment<T>::~Experiment() {
}

template<typename T>
void Experiment<T>::SetExample(int eid, int arth_mode) {
	switch (eid) {
	case 1:
		//_example = new Example01();
		break;
	case 2:
		//_example = new Example02();
		break;
	case 3:
		//_example = new Example03<T>();
		break;
	case 4:
		_example = new Example04<T>();
		break;
	case 9:
		_example = new Example09<T>();
	default:
		_example = NULL;
		break;
	}
}

template<typename T>
Experiment<T>::Experiment(int ac, char* av[]) {
	_example = NULL;
	solver = NULL;
	_param_initialized = false;
	_solver_initialized = false;
	ConfigReader<long double>* cf = new ConfigReader<long double>(ac, av);
	this->SetParameters(cf->GetParameters());
	delete cf;
}

template<typename T>
void Experiment<T>::SetSolver(Parameters<long double> p) {
	switch (p.selected_solver) {
	case GPDE_SOLVER:
		solver = new GPDESolver<T>();
		break;
	case POISSON:
		solver = new PoissonSolver<T>();
		break;
	case POISSONAM:
		solver = new PoissonSolverAM<T>();
		break;

	case POISSON4:
		solver = new PoissonSolver4Order<T>();
		break;
	case POISSON4AM:
		solver = new PoissonSolver4OrderAM<T>();
		break;
	default:
		solver = NULL;
		break;
	}
	_solver_initialized = (solver != NULL);

	if (_solver_initialized) {
		solver->SetParameters(parameters);
	}
}

template<typename T>
void Experiment<T>::SetParameters(Parameters<long double> p) {
	this->parameters = p;
	_param_initialized = true;
}

template<typename T>
void Experiment<T>::Initialize() {
	Interval<T>::Initialize();
	if (!_param_initialized)
		return;
	this->SetSolver(parameters);
}

template<typename T>
void Experiment<T>::SetExampleForSolver(int eid) {
	if (_solver_initialized)
		solver->SetExample(eid);
}

template<typename T>
void Experiment<T>::Execute() {
	if (_solver_initialized) {
		solver->Execute();
		solver->WriteResults();
	}
}

//The explicit instantiation part
template class Experiment<long double> ;
template class Experiment<mpreal> ;
}

