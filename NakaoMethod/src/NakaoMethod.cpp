//============================================================================
// Name        : NakaoMethod.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
//#include "Interval.h"
//#include "NakaoExperiment.h"
//#include "NakaoExperiment2D.h"
//#include "Integration.h"

#include "NakaoExperimentApprox.h"

using namespace std;

int main() {

//	int st;
//	long double iex = GaussLegendreTH(e, 4, GL_X4, GL_A4, -3.0, 3.0, st);
//
//	std::cout << std::setprecision(10);
//	std::cout << "Integrating ExpTH(X) over [-3, 3]: " << iex << '\n';
//    std::cout << "Actual value:                    " << e(3) -e(-3) << '\n';
//
//	    GaussLegendreQuadrature<5> gl5;
//
//
//
//	    gl5.print_roots_and_weights(std::cout);
//	    std::cout << "Integrating Exp(X) over [-3, 3]: " << gl5.integrate(-3., 3., e) << '\n';
//	    std::cout << "Actual value:                    " << e(3) -e(-3) << '\n';

//	NakaoExperiment* exper1D = new NakaoExperiment();
//	exper1D->execute();
//	delete exper1D;

//	NakaoExperiment2D* exper = new NakaoExperiment2D();
//	exper->execute();
//	delete exper;

	    	NakaoExperimentApprox* exper1D = new NakaoExperimentApprox();
	    	exper1D->execute();
	    	delete exper1D;

	return 0;
}
