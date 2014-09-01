/*
 * Tester.h
 *
 *  Created on: 12-04-2014
 *      Author: thof
 */

#ifndef TESTER_H_
#define TESTER_H_

#include <string>
#include <ieee754.h>
#include <bitset>

#include "IntervalArithmetic.h"

using namespace std;

namespace intervalarth
{

class Tester
{
public:
	Tester();
	virtual ~Tester();
	void PrintBinary(long double x);
	void IEndsToStrings(interval i, string& left, string& right);
	void ArithmeticTest();
	void StuffTest();
};

} /* namespace intervalarth */
#endif /* TESTER_H_ */
