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
		_example = new ExampleGPE01<T>();
		break;
	case 2:
		_example = new ExampleGPE02<T>();
		break;
	case 3:
		_example = new ExampleGPE03<T>();
		break;
	case 4:
		_example = new ExampleGPE04<T>();
		break;
	case 5:
		_example = new ExampleGPE05<T>();
		break;
	case 6:
		_example = new ExampleGPE06<T>();
		break;
	case 7:
		_example = new ExampleGPE07<T>();
		break;
	case 8:
		_example = new ExampleGPE08<T>();
		break;
	case 9:
		//_example = new Example09<T>();
		//break;
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
	case GPE_SOLVER5C:
		solver = new GPE5CSolver<T>();
		break;
	case GPE_SOLVER3C:
		solver = new GPESolver3C<T>();
		break;
//	case POISSON:
//		solver = new PoissonSolver<T>();
//		break;
//	case POISSONAM:
//		solver = new PoissonSolverAM<T>();
//		break;
//
//	case POISSON4:
//		solver = new PoissonSolver4Order<T>();
//		break;
//	case POISSON4AM:
//		solver = new PoissonSolver4OrderAM<T>();
//		break;
//	case POISSON6:
//		solver = new PoissonSolver6Order<T>();
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
	Interval<T>::SetMode(p.ia_mode);
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
//template class Experiment<mpreal> ;
}

