//============================================================================
// Name        : IntervalArithmeticTest.cpp
// Author      : Tomasz Hoffmann
// Version     :
// Copyright   : TH
// Description : IntervalAritmetic module test project
//============================================================================

#include <iostream>
#include "Tester.h"
#include "Experiment.h"
#include "Interval.h"

using namespace std;
using namespace interval_arithmetic;

int main(int ac, char *av[]) {
	cout << "\n \n ------------------------------------------------------- \n \n";
	cout << "\n \n ------TEST IMPLEMENTATION OF INTERVAL ARITHMETIC------- \n \n";
	cout << "\n \n ------------------------------------------------------- \n \n";

	//simple tester for basic interval operations
//	Tester* tester = new Tester();
//	tester->ArithmeticTestNew();
//	delete tester;

	Interval<long double> ione = {1, 1};
	Interval<long double> a = {-2, -1};
	Interval<long double> b = {-3, -1};
	Interval<long double> c, d;

	c = ione /a;
	d = ione /b;

	//experiment for the Generalized Poisson Equation
	Experiment<long double>* exper = new Experiment<long double>(ac, av);
	exper->Initialize();
	exper->Execute();
	delete exper;

	return 0;
}
