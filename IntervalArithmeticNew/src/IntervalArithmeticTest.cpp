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

using namespace std;
using namespace interval_arithmetic;

int main(int ac, char *av[]) {
	cout << "IntervalArtithmetic CPP test project" << endl;

	cout << "\n \n ------------------------------------------------------- \n \n";
	cout << "\n \n -----------------NEW IMPLEMENTATION-------------------- \n \n";
	cout << "\n \n ------------------------------------------------------- \n \n";

	Tester* tester = new Tester();
	tester->ArithmeticTestNew();
	delete tester;

	Experiment<long double>* exper = new Experiment<long double>(ac, av);
	exper->Initialize();
	exper->Execute();
	delete exper;

	return 0;
}
