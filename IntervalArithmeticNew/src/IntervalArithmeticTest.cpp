//============================================================================
// Name        : IntervalArithmeticTest.cpp
// Author      : Tomasz Hoffmann
// Version     :
// Copyright   : TH
// Description : IntervalAritmetic module test project
//============================================================================

#include <iostream>
//#include "Tester.h"
//#include "Experiment.h"
//#include "Interval.h"
#include "IntervalMpreal.h"

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

//	Interval<mpreal> ione = {1, 1};
//	Interval<mpreal> a = {-2, -1};
//	Interval<mpreal> b = {-3, -1};
//	Interval<mpreal> c, d;

	Interval<mpreal> ione = {1, 1};
	Interval<mpreal> a = {-2, -1};
	Interval<mpreal> b = {-3, -1};
	Interval<mpreal> c, d;

	c = ione /a;
	d = ione /b;

	//experiment for the Generalized Poisson Equation
//	Experiment<long double>* exper = new Experiment<long double>(ac, av);
//	exper->Initialize();
//	exper->Execute();
//	delete exper;

	return 0;
}
