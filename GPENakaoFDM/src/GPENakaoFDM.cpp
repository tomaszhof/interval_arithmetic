//============================================================================
// Name        : GPENakaoFDM.cpp
// Author      : TH
// Version     :
// Copyright   : 
// Description : GPENakaoFDM  - Elliptic PDE FDM Solvers
//============================================================================

#include <iostream>
#include <chrono>

#include "Interval.h"
#include "Experiment.h"

using namespace std;

int main(int ac, char *av[]) {
	cout << "GPENakaoExperiment  - Elliptic PDE FDM Solvers" << endl;

	//experiment for the Generalized Poisson Equation
	Experiment<long double>* exper = new Experiment<long double>(ac, av);
	exper->Initialize();

	auto t1 = std::chrono::high_resolution_clock::now();
	exper->Execute();
	auto t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
	cout << duration;

	delete exper;
	return 0;
}
