//============================================================================
// Name        : NakaoMethod.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Interval.h"
#include "NakaoExperiment.h"
#include "NakaoExperiment2D.h"

using namespace std;

int main() {
	NakaoExperiment2D* exper = new NakaoExperiment2D();
	exper->execute();
	delete exper;

	return 0;
}
