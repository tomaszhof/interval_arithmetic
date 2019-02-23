//============================================================================
// Name        : GPENakaoFDM.cpp
// Author      : TH
// Version     :
// Copyright   : 
// Description : GPENakaoFDM  - Elliptic PDE FDM Solvers
//============================================================================

#include <iostream>

#include "Interval.h"
#include "Experiment.h"

using namespace std;

int main(int ac, char *av[]) {
	cout << "GPENakaoExperiment  - Elliptic PDE FDM Solvers" << endl;

	//experiment for the Generalized Poisson Equation
	Experiment<long double>* exper = new Experiment<long double>(ac, av);
	exper->Initialize();
	exper->Execute();
	delete exper;
	return 0;
}
