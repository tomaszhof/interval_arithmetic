/*
 * Interval.h
 *
 *  Created on: 7 wrz 2014
 *      Author: tomhof
 */

#ifndef INTERVAL_H_
#define INTERVAL_H_

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
#include <string.h>
#include <iomanip>
#include <fstream>
#include <float.h>
#include <typeinfo>

using namespace std;


namespace interval_arithmetic {

enum IAMode {
	DINT_MODE, PINT_MODE
};

enum IAPrecision {
	LONGDOUBLE_PREC=80, FLOAT_PREC=60
};

enum IAOutDigits
{
	LONGDOUBLE_DIGITS=17, FLOAT_DIGITS=17
};

template<typename T> class Interval {
private:
	static IAMode mode;
	static IAPrecision precision;
	static IAOutDigits outdigits;
public:
	T a;
	T b;
	Interval();
	Interval(T a, T b);
	virtual ~Interval();
	Interval<T> operator=(const Interval<T>& i);
	Interval<T> operator+(const Interval<T>& i);
	Interval<T> operator-(const Interval<T>& i);
	Interval<T> operator*(const Interval<T>& i);
	Interval<T> operator/(const Interval<T>& i);
	static void SetMode(IAMode m);
	static IAMode GetMode(IAMode m);
	static void SetPrecision(IAPrecision p);
	static IAPrecision GetPrecision(IAPrecision p);
	static void SetOutDigits(IAOutDigits o);
	static IAOutDigits GetOutDigits(IAOutDigits o);
	static Interval<T> IntRead(const string & sa);
	void IEndsToStrings(string & left, string & right);
	T LeftRead(const string& sa);
	T RightRead(const string& sa);
};

template<typename T>
inline Interval<T>::~Interval() {
}

template<typename T>
Interval<T>::Interval() {
	this->a = 0;
	this->b = 0;
}

template<typename T>
inline Interval<T>::Interval(T a, T b) {
	this->a = a;
	this->b = b;
}

template<typename T>
inline void Interval<T>::SetMode(IAMode m) {
	Interval<T>::mode = m;
}

template<typename T>
inline Interval<T> Interval<T>::operator /(const Interval<T>& y) {
	Interval<T> r;
	Interval<T> x(this->a, this->b);
	T x1y1, x1y2, x2y1, t;

	if ((y.a <= 0) && (y.b >= 0)) {
		throw runtime_error("Division by an interval containing 0.");
	} else {
		fesetround(FE_DOWNWARD);
		x1y1 = x.a / y.a;
		x1y2 = x.a / y.b;
		x2y1 = x.b / y.a;
		r.a = x.b / y.b;
		t = r.a;
		if (x2y1 < t)
			r.a = x2y1;
		if (x1y2 < t)
			r.a = x1y2;
		if (x1y1 < t)
			r.a = x1y1;

		fesetround(FE_UPWARD);
		x1y1 = x.a / y.a;
		x1y2 = x.a / y.b;
		x2y1 = x.b / y.a;

		r.b = x.b / y.b;
		t = r.b;
		if (x2y1 > t)
			r.b = x2y1;
		if (x1y2 > t)
			r.b = x1y2;
		if (x1y1 > t)
			r.b = x1y1;

	}
	fesetround(FE_TONEAREST);
	return r;
}

template<typename T> IAMode Interval<T>::mode = PINT_MODE;
template<typename T> IAPrecision Interval<T>::precision = LONGDOUBLE_PREC;
template<typename T> IAOutDigits Interval<T>::outdigits = LONGDOUBLE_DIGITS;

template<typename T>
inline IAMode Interval<T>::GetMode(IAMode m) {
	return Interval<T>::mode;
}

template<typename T>
inline Interval<T> Interval<T>::operator =(const Interval<T>& i) {
	this->a = i.a;
	this->b = i.b;
}

template<typename T>
inline Interval<T> Interval<T>::operator +(const Interval<T>& i) {
	Interval<T> r;
	fesetround(FE_DOWNWARD);
	r.a = this->a + i.a;
	fesetround(FE_UPWARD);
	r.b = this->b + i.b;
	fesetround(FE_TONEAREST);
	return r;
}

template<typename T>
inline Interval<T> Interval<T>::operator -(const Interval<T>& i) {
	Interval<T> r;
	fesetround(FE_DOWNWARD);
	r.a = this->a - i.a;
	fesetround(FE_UPWARD);
	r.b = this->b - i.b;
	fesetround(FE_TONEAREST);
	return r;
}

template<typename T>
inline Interval<T> Interval<T>::operator *(const Interval<T>& y) {
	Interval<T> r(0, 0);
	T x1y1, x1y2, x2y1;

	Interval<T> x(this->a, this->b);
	fesetround(FE_DOWNWARD);
	x1y1 = x.a * y.a;
	x1y2 = x.a * y.b;
	x2y1 = x.b * y.a;
	r.a = x.b * y.b;
	if (x2y1 < r.a)
		r.a = x2y1;
	if (x1y2 < r.a)
		r.a = x1y2;
	if (x1y1 < r.a)
		r.a = x1y1;

	fesetround(FE_UPWARD);
	x1y1 = x.a * y.a;
	x1y2 = x.a * y.b;
	x2y1 = x.b * y.a;

	r.b = x.b * y.b;
	if (x2y1 > r.b)
		r.b = x2y1;
	if (x1y2 > r.b)
		r.b = x1y2;
	if (x1y1 > r.b)
		r.b = x1y1;
	fesetround(FE_TONEAREST);
	return r;
}

template<typename T>
inline void Interval<T>::SetPrecision(IAPrecision p) {
	Interval<T>::precision = p;
}

template<typename T>
inline IAPrecision Interval<T>::GetPrecision(IAPrecision p) {
	return Interval<T>::precision;
}

template<typename T>
inline void Interval<T>::SetOutDigits(IAOutDigits o) {
	Interval<T>::outdigits = 0;
}

template<typename T>
inline IAOutDigits Interval<T>::GetOutDigits(IAOutDigits o) {
	return Interval<T>::outdigits;
}

template<typename T>
inline Interval<T> Interval<T>::IntRead(const string& sa) {
	Interval<T> r;
		mpfr_t rop;
		mpfr_init2(rop, precision);
		mpfr_set_str(rop, sa.c_str(), 10, MPFR_RNDD);
		T le = 0.0;
		if (strcmp(typeid(T).name(), typeid(long double).name()) == 0)
		{
			le = mpfr_get_ld(rop, MPFR_RNDD);
		}
		if (strcmp(typeid(T).name(), typeid(double).name()) == 0)
		{
			le = mpfr_get_d(rop, MPFR_RNDD);
		}
		if (strcmp(typeid(T).name(), typeid(float).name()) == 0)
		{
			le = mpfr_get_flt(rop, MPFR_RNDD);
		}

		mpfr_set_str(rop, sa.c_str(), 10, MPFR_RNDU);
		T re = 0.0;
		if (strcmp(typeid(T).name(), typeid(long double).name()) == 0)
		{
			re = mpfr_get_ld(rop, MPFR_RNDU);
		}
		if (strcmp(typeid(T).name(), typeid(double).name()) == 0)
		{
			re = mpfr_get_d(rop, MPFR_RNDU);
		}
		if (strcmp(typeid(T).name(), typeid(float).name()) == 0)
		{
			re = mpfr_get_flt(rop, MPFR_RNDU);
		}
		fesetround(FE_TONEAREST);

		r.a = le;
		r.b = re;
		return r;
}

template<typename T>
inline void Interval<T>::IEndsToStrings(string& left, string& right) {
	mpfr_t rop;
	mpfr_exp_t exponent;
	mpfr_init2(rop, precision);
	char* str = NULL;
	char *buffer = new char(precision + 3);
	mpfr_set_ld(rop, this->a, MPFR_RNDD);

	mpfr_get_str(buffer, &exponent, 10, outdigits, rop, MPFR_RNDD);
	str = buffer;

	stringstream ss;
	int prec = std::numeric_limits<T>::digits10;
	ss.setf(std::ios_base::scientific);
	bool minus = (str[0] == '-');
	int splitpoint = minus ? 1 : 0;
	string sign = minus ? "-" : "";

	ss << std::setprecision(prec) << sign << str[splitpoint] << "."
			<< &str[splitpoint + 1] << "E" << exponent - 1;
	left = ss.str();
	ss.str(std::string());

	mpfr_set_ld(rop, this->b, MPFR_RNDU);
	mpfr_get_str(buffer, &exponent, 10, outdigits, rop, MPFR_RNDU);
	str = buffer;
	splitpoint = (str[0] == '-') ? 1 : 0;
	ss << std::setprecision(prec) << sign << str[splitpoint] << "."
			<< &str[splitpoint + 1] << "E" << exponent - 1;
	right = ss.str();
	ss.clear();
}

template<typename T>
inline T Interval<T>::LeftRead(const string& sa) {
	Interval<T> int_number;
	int_number = IntRead(sa);
	return int_number.a;
}



template<typename T>
inline T Interval<T>::RightRead(const string& sa) {
	Interval<T> int_number;
	int_number = IntRead(sa);
	return int_number.b;
}

} /* namespace interval_arithmetic */

#endif /* INTERVAL_H_ */
