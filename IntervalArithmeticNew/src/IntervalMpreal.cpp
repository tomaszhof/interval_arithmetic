/*
 * IntervalMpreal.cpp
 *
 *  Created on: 13 wrz 2014
 *      Author: tomhof
 */
#include "Interval.h"
#include <mpfr.h>
#include <mpreal.h>

using namespace std;
using namespace mpfr;

namespace interval_arithmetic {
template<>
class Interval<mpreal> {
private:
	mpreal a;
	mpreal b;
public:
	Interval();
	virtual ~Interval();
	void Abc(const Interval<mpreal>& i);
	//void operator+(const Interval<long double>& i);
};

Interval<mpreal>::Interval()
{
	this->a = 0.0;
	this->b = 0.0;
}
Interval<mpreal>::~Interval()
{
}

void Interval<mpreal>::Abc(const Interval<mpreal>& i){
	Interval<mpreal> r;
	//fesetround(FE_DOWNWARD);
	r.a = this->a + i.a;
	//fesetround(FE_UPWARD);
	r.b = this->b + i.b;
	//fesetround(FE_TONEAREST);
//	/return r;
}

}
