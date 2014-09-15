//============================================================================
// Name        : IntervalArithmeticTest.cpp
// Author      : Tomasz Hoffmann
// Version     :
// Copyright   : TH
// Description : IntervalAritmetic module test project
//============================================================================

#include <iostream>
#include "Tester.h"

using namespace std;
using namespace intervalarth;

int main() {
	cout << "IntervalArtithmetic CPP test project" << endl;

	Tester* tester = new Tester();
	tester->ArithmeticTest();
	cout << "\n \n ------------------------------------------------------- \n \n";
	cout << "\n \n -----------------NEW IMPLEMENTATION-------------------- \n \n";
	cout << "\n \n ------------------------------------------------------- \n \n";
	tester->ArithmeticTestNew();
	delete tester;

	return 0;
}
