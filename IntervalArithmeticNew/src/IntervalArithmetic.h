/*
 * IntervalArithmetic.h
 *
 *  Created on: 20-11-2012
 *      Author: Tomasz Hoffmann and Andrzej Marciniak
 *
 *  Before you start use module, please install libraries:
 *  	- .boost ( http://www.boost.org/ )
 *  	- GNU MPFR ( http://www.mpfr.org/ )
 */

#ifndef INTERVALARITHMETIC_H_
#define INTERVALARITHMETIC_H_
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <fenv.h>
#include <stdlib.h>
#include <stdint.h>
#include <cmath>
#include <mpfr.h>
#include <boost/lexical_cast.hpp>

#include <iomanip>
#include <fstream>
#include <float.h>


using namespace std;

namespace intervalarth
{

struct interval
{
	long double a, b;
};

//set precision
const int curr_precision = 80;

// store the original rounding mode
const int originalRounding = fegetround();
const char DecimalSeparator = '.';

const int DMODE = 1;
const int PMODE = 2;

class IntervalArithmetic
{
private:
	int mode;
public:
	IntervalArithmetic();
	virtual ~IntervalArithmetic();
	long double IntWidth(const interval& x);
	interval IAdd(const interval& x, const interval& y);
	interval ISub(const interval& x, const interval& y);
	interval IMul(const interval& x, const interval& y);
	interval IDiv(const interval& x, const interval& y);
	long double DIntWidth(const interval& x);
	interval Projection(const interval& x);
	interval Opposite(const interval& x);
	interval Inverse(const interval& x);
	interval DIAdd(const interval& x, const interval& y);
	interval DISub(const interval& x, const interval& y);
	interval DIMul(const interval& x, const interval& y);
	interval DIDiv(const interval& x, const interval& y);
	interval IntRead(const string& sa);
	long double LeftRead(const string& sa);
	long double RightRead(const string& sa);
	interval ISin(const interval& x, int& st);
	interval ICos(const interval& x, int& st);
	interval IExp(const interval& x, int& st);
	interval ISqr(const interval& x, int& st);
	interval DISin(const interval& x, int& st);
	interval DICos(const interval& x, int& st);
	interval DIExp(const interval& x, int& st);
	interval DISqr(const interval& x, int& st);
	interval ISqr2();
	interval ISqr3();
	interval IPi();
	void IEndsToStrings(const interval& i, string& left, string& right);
	void SetArthMode(int m);
	int GetArthMode();
};

} /* namespace interval */
#endif /* INTERVALARITHMETIC_H_ */
