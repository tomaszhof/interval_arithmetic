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
	delete tester;

	return 0;
}
