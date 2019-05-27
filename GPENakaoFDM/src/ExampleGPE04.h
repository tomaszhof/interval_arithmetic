/*
 * Example09.h
 *
 *  Created on: 11-11-2014
 *      Author: thof
 */

#ifndef EXAMPLEGPE04_H_
#define EXAMPLEGPE04_H_

#include "Utils.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class ExampleGPE04: public BoundaryConditions<T> {
	const long double PI = 3.141592653589793238L;
public:
	ExampleGPE04();
	virtual ~ExampleGPE04();
	Interval<T> F(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> PHI1(const Interval<T>& iy, int& st);
	Interval<T> PHI2(const Interval<T>& ix, int& st);
	Interval<T> PHI3(const Interval<T>& iy, int& st);
	Interval<T> PHI4(const Interval<T>& ix, int& st);
	Interval<T> PSI(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> OMEGA(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> PSI1(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> PSI2(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> PSI3(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> PSI4(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> A1(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> A2(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> DA1DX(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> DA1DY(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> DA2DX(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> DA2DY(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> D2A1DX2(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> D2A1DY2(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> D2A2DX2(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> D2A2DY2(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> C(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> DCDX(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> DCDY(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> D2CDX2(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> D2CDY2(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> DFDX(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> DFDY(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> D2FDX2(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> D2FDY2(const Interval<T>& ix, const Interval<T>& iy, int& st);

	//classical functions
	long double f(const long double& x, const long double& y);
	long double phi1(const long double& y);
	long double phi2(const long double& x);
	long double phi3(const long double& y);
	long double phi4(const long double& x);
	long double a1(const long double& x, const long double& y);
	long double a2(const long double& x, const long double& y);
	long double da1dx(const long double& x, const long double& y);
	long double da1dy(const long double& x, const long double& y);
	long double da2dx(const long double& x, const long double& y);
	long double da2dy(const long double& x, const long double& y);
	long double d2a1dx2(const long double& x, const long double& y);
	long double d2a1dy2(const long double& x, const long double& y);
	long double d2a2dx2(const long double& x, const long double& y);
	long double d2a2dy2(const long double& x, const long double& y);
	long double c(const long double& x, const long double& y);
	long double dcdx(const long double& x, const long double& y);
	long double dcdy(const long double& x, const long double& y);
	long double d2cdx2(const long double& x, const long double& y);
	long double d2cdy2(const long double& x, const long double& y);

	int boundconds_classic(const long double& b1, const long double& b2,
			const long double eps);

	long double ExactSol(long double x, long double y);
	long double GetConstM();
	long double GetConstN();
	long double GetConstP();
	long double GetConstQ();
	long double GetConstR();
	long double GetConstS();
	void SetArithmeticMode(IAMode mode);

	const Interval<long double> i0 = { 0.0, 0.0 };
	const Interval<long double> i1 = { 1.0, 1.0 };
	const Interval<long double> i2 = { 2.0, 2.0 };
	const Interval<long double> i3 = { 3.0, 3.0 };
	const Interval<long double> i4 = { 4.0, 4.0 };
	const Interval<long double> i5 = { 5.0, 5.0 };
	const Interval<long double> i6 = { 6.0, 6.0 };
	const Interval<long double> i7 = { 7.0, 7.0 };
	const Interval<long double> i10 = { 10.0, 10.0 };
	const Interval<long double> i12 = { 12.0, 12.0 };
	const Interval<long double> i15 = { 15.0, 15.0 };
	const Interval<long double> i20 = { 20.0, 20.0 };
	const Interval<long double> i24 = { 24.0, 24.0 };
	const Interval<long double> i36 = { 36.0, 36.0 };
	const Interval<long double> i45 = { 45.0, 45.0 };
	const Interval<long double> i120 = { 120.0, 120.0 };
	const Interval<long double> i180 = { 180.0, 180.0 };
	const Interval<long double> i360 = { 360.0, 360.0 };
	const Interval<long double> im1 = { -1.0, -1.0 };
	const Interval<long double> im11 = { -1.0, 1.0 };
	const Interval<long double> ipi = Interval<long double>::IPi();
};

template<typename T>
ExampleGPE04<T>::ExampleGPE04() {
	// TODO Auto-generated constructor stub

}

template<typename T>
ExampleGPE04<T>::~ExampleGPE04() {
	// TODO Auto-generated destructor stub
}

template<typename T>
long double ExampleGPE04<T>::f(const long double& x, const long double& y) {
//	long double r = x*exp(x*y);
////	cout << "r=" << r << endl;
//	long double rest = (x*x*x*y*sin(PI*y) + x*x*y*y*sin(PI*x) + 2*x*x*sin(PI*y) - y);
////	cout << "rest=" << rest << endl;
//	r = r * rest;
////	cout << "r*rest=" << r << endl;
	long double r = x * y * exp(x * y)
			* (x * y * sin(PI * y) + x * y * sin(PI * x) - 1);
	return r;
}

template<typename T>
long double ExampleGPE04<T>::phi1(const long double& y) {
	return 1.0L;
}

template<typename T>
long double ExampleGPE04<T>::phi2(const long double& x) {
	return 1.0L;
}

template<typename T>
long double ExampleGPE04<T>::phi3(const long double& y) {
	return exp(y);
}

template<typename T>
long double ExampleGPE04<T>::phi4(const long double& x) {
	return exp(x);
}

template<typename T>
long double ExampleGPE04<T>::c(const long double& x, const long double& y) {
	return -x * y;
}

template<typename T>
long double ExampleGPE04<T>::dcdx(const long double& x, const long double& y) {
	return -y;
}

template<typename T>
long double ExampleGPE04<T>::dcdy(const long double& x, const long double& y) {
	return -x;
}

template<typename T>
long double ExampleGPE04<T>::d2cdx2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE04<T>::d2cdy2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE04<T>::a1(const long double& x, const long double& y) {
	return (x * x) * sin(PI * y);
}

template<typename T>
long double ExampleGPE04<T>::da1dx(const long double& x, const long double& y) {
	return (2.0 * x) * sin(PI * y);
}

template<typename T>
long double ExampleGPE04<T>::da1dy(const long double& x, const long double& y) {
	return (PI * x * x) * cos(PI * y);
}

template<typename T>
long double ExampleGPE04<T>::d2a1dx2(const long double& x,
		const long double& y) {
	return 2.0L * sin(PI * y);
}

template<typename T>
long double ExampleGPE04<T>::d2a1dy2(const long double& x,
		const long double& y) {
	return -PI * PI * (x * x) * sin(PI * y);
}

template<typename T>
long double ExampleGPE04<T>::a2(const long double& x, const long double& y) {
	return (y * y) * sin(PI * x);
}

template<typename T>
long double ExampleGPE04<T>::da2dx(const long double& x, const long double& y) {
	return (PI * y * y) * cos(PI * x);
}

template<typename T>
long double ExampleGPE04<T>::da2dy(const long double& x, const long double& y) {
	return (2.0 * y) * sin(PI * x);
}

template<typename T>
long double ExampleGPE04<T>::d2a2dx2(const long double& x,
		const long double& y) {
	return -PI * PI * (y * y) * sin(PI * x);
}

template<typename T>
long double ExampleGPE04<T>::d2a2dy2(const long double& x,
		const long double& y) {
	return 2.0L * sin(PI * x);
}

template<typename T>
Interval<T> ExampleGPE04<T>::F(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//f(x,y)= x*exp(x*y)*(x*x*x*y*sin(PI*y) + x*x*y*y*sin(PI*x) + 2*x*x*sin(PI*y) - y);
	//x*exp(x*y)*(x*x*x*y*sin(PI*y) + x*x*y*y*sin(PI*x) + 2*x*x*sin(PI*y) - y);
	Interval<T> r = { 0.0L, 0.0L };
//	r = (ix * IExp(ix * iy));
//	cout << "[" << ix.a << ", " << ix.b << "] " << "ix*IExp(ix*iy): r= [" << r.a << " ; " << r.b << "]" << endl;
//	Interval<T> rest = (ix*ix*ix*iy*ISin(ipi * iy) + ix*ix*iy*iy*ISin(ipi*ix) + i2*ix*ix*ISin(ipi*iy) - iy);
//	cout << "rest= [" << rest.a << " ; " << rest.b << "]" << endl;
//	r = r * rest;
//	cout << "r*rest= [" << r.a << " ; " << r.b << "]" << endl;
	r = ix * iy * IExp(ix * iy);
	r = r * (ix * iy * ISin(ipi * iy) + ix * iy * ISin(ipi * ix) - i1);
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE04<T>::PHI1(const Interval<T>& iy, int& st) {
	//phi1(y) = 0.0;
	Interval<T> r = { 1.0, 1.0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE04<T>::PHI2(const Interval<T>& ix, int& st) {
	//phi2(x) = 0.0;
	Interval<T> r = { 1.0, 1.0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE04<T>::PHI3(const Interval<T>& iy, int& st) {
	//phi3(y) = 0.0;
	Interval<T> r = { 0, 0 };

	r = IExp(iy);
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE04<T>::PHI4(const Interval<T>& ix, int& st) {
	//phi3(y) = 0.0;
	Interval<T> r = { 0, 0 };
	st = 0;
	r = IExp(ix);
	return r;
}

template<typename T>
Interval<T> ExampleGPE04<T>::PSI(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//psi(x,y) = 0
	Interval<T> r = { 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE04<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//omega(x,y) = 0
	Interval<T> r = { 0, 0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE04<T>::PSI1(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE04<T>::PSI2(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE04<T>::PSI3(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE04<T>::PSI4(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE04<T>::C(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//c(x,y) = x * y;

	return im1 * ix * iy;
}

template<typename T>
long double ExampleGPE04<T>::ExactSol(long double x, long double y) {
	return exp(x * y);
}

template<typename T>
long double ExampleGPE04<T>::GetConstM() {
	long double constM = 0.0;

	return constM;
}

template<typename T>
long double ExampleGPE04<T>::GetConstN() {
	long double constN = 0.0;

	return constN;
}

template<typename T>
long double ExampleGPE04<T>::GetConstP() {
	long double constP = 8.16L;

	return constP;
}

template<typename T>
long double ExampleGPE04<T>::GetConstQ() {
	long double constQ = 8.16L;

	return constQ;
}

template<typename T>
long double ExampleGPE04<T>::GetConstR() {
	long double constR = 8.16L;

	return constR;
}

template<typename T>
long double ExampleGPE04<T>::GetConstS() {
	long double constS = 19.03L;

	return constS;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::A1(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {

	//a1(x,y) = (x*x)*sin(PI*y);
	Interval<T> r = { 0.0L, 0.0L };
	st = 0;

	r = ix * ix * ISin(PI * iy);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::A2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	//a2(x,y) = (y*y)*sin(PI*x);
	Interval<T> r = { 0.0L, 0.0L };
	st = 0;
	r = iy * iy * ISin(PI * ix);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::D2A1DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	//2.0L*sin(PI*y);
	Interval<T> r = { 0, 0 };
	st = 0;
	r = i2 * ISin(ipi * iy);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::D2A1DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	//-PI*PI*(x*x)*sin(PI*y);
	Interval<T> r = { 0, 0 };
	st = 0;
	r = im1 * ipi * ipi * ix * ix * ISin(ipi * iy);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::D2A2DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	//-PI*PI*(y*y)*sin(PI*x);
	Interval<T> r = { 0, 0 };
	st = 0;

	r = im1 * ipi * ipi * iy * iy * ISin(ipi * ix);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::D2A2DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	//2.0L*sin(PI*x);
	Interval<T> r = { 0, 0 };
	st = 0;

	r = i2 * ISin(ipi * ix);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::DCDX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	r = im1 * iy;
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::DCDY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	r = im1 * ix;
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::D2CDX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::D2CDY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::DFDX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
//	Interval<T> iexpXY = IExp(ix * iy);
//	Interval<T> isinPIX = ISin(ipi * ix);
//	Interval<T> isinPIY = IExp(ipi * iy);
//	Interval<T> icosPIX = IExp(ipi * ix);
//
//	r = iy * iexpXY
//			* (ix * ix * iy * iy * isinPIX + ix * ix * iy * iy * isinPIX
//					+ ix * ix * iy * iy * isinPIY + ipi * ix * ix * iy * icosPIX
//					- ix * iy + i2 * ix * iy * isinPIX + i2 * ix * iy * isinPIY
//					- i1);

	Interval<T> s = ISin(ipi*(ix+iy)/i2);
	Interval<T> c = ICos(ipi*(ix+iy)/i2);
	Interval<T> xy = ix*iy;

	r = i2 * xy * (xy + i2) * s * c;

	c = ICos(ipi*ix);
	r = r + xy * (ipi * ix * c - i1) - i1;
	r = iy * IExp(xy)*r;
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::DFDY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
//	Interval<T> iexpXY = IExp(ix * iy);
//	Interval<T> isinPIX = ISin(ipi * ix);
//	Interval<T> isinPIY = IExp(ipi * iy);
//	Interval<T> icosPIY = IExp(ipi * iy);
//
//	st = 0;
//
//	r = ix * iexpXY
//			* (ix * ix * iy * iy * isinPIX + ix * ix * iy * iy * isinPIY
//					+ ipi * ix * iy * iy * icosPIY - ix * iy
//					+ i2 * ix * iy * isinPIX + i2 * ix * iy * isinPIY - i1);

	Interval<T> s = ISin(ipi*(ix+iy)/i2);
	Interval<T> c = ICos(ipi*(ix+iy)/i2);
	Interval<T> xy = ix*iy;

	r = i2 * xy * (xy + i2) * s * c;

	c = ICos(ipi*iy);
	r = r + xy * (ipi * iy * c - i1) - i1;
	r = ix * IExp(xy)*r;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::D2FDX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
//	Interval<T> iexpXY = IExp(ix * iy);
//	Interval<T> isinPIX = ISin(ipi * ix);
//	Interval<T> isinPIY = IExp(ipi * iy);
//	Interval<T> icosPIX = IExp(ipi * ix);
//
//	r = iy * (ix * iy * iy * iexpXY)
//			* (ix * iy * isinPIX + ix * iy * isinPIY - i1);
//	r = r
//			+ ix * iy * iexpXY
//					* (i2 * ipi * iy * icosPIX - ipi * ipi * ix * iy * isinPIX);
//	r = r
//			+ i2 * iy * (ix * iy * iexpXY + iexpXY)
//					* (iy * isinPIX + ipi * ix * iy * icosPIX + iy * isinPIY);


	Interval<T> s = ISin(ipi*(ix + iy)/i2);
	Interval<T> c = ICos(ipi*(ix - iy)/i2);
	Interval<T> xy = ix * iy;
	r = i2 * (xy * (xy + i4) + i2) * s * c;

	c = ICos(ipi*ix);
	s = ISin(ipi*ix);

	r = r + ipi*ix*(i2*(xy + i2)*c - ipi*ix*s) - xy - i2;
	r = iy *iy * IExp(xy)*r;

	st = 0;


	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::D2FDY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
//	Interval<T> iexpXY = IExp(ix * iy);
//	Interval<T> isinPIX = ISin(ipi * ix);
//	Interval<T> isinPIY = IExp(ipi * iy);
//	Interval<T> icosPIX = IExp(ipi * ix);
//	Interval<T> icosPIY = IExp(ipi * iy);
//
//	r = ix * ix * iexpXY
//			* (ix * ix * isinPIX + ix * ix * iy * iy * isinPIY
//					+ i2 * ipi * ix * iy * iy * icosPIY - ix * iy
//					+ i4 * ix * iy * isinPIX + i4 * ix * iy * isinPIY
//					+ i2 * isinPIX - ipi * ipi * iy * iy * isinPIY
//					+ i2 * isinPIY + i4 * ipi * iy * icosPIY - i2);

	Interval<T> s = ISin(ipi*(ix + iy)/i2);
	Interval<T> c = ICos(ipi*(ix - iy)/i2);
	Interval<T> xy = ix * iy;
	r = i2 * (xy * (xy + i4) + i2) * s * c;

	c = ICos(ipi*iy);
	s = ISin(ipi*iy);

	r = r + ipi*iy*(i2*(xy + i2)*c - ipi*iy*s) - xy - i2;
	r = ix *ix * IExp(xy)*r;

	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::DA1DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	r = i2 * ix * ISin(ipi * iy);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::DA1DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	r = ipi * ix * ix * ICos(ipi * iy);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::DA2DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	r = ipi * iy * iy * ICos(ipi * ix);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE04<T>::DA2DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	r = i2 * iy * ISin(PI * ix);
	return r;
}
template<typename T>
void ExampleGPE04<T>::SetArithmeticMode(IAMode mode) {
	Interval<T>::SetMode(mode);
}

template<typename T>
int ExampleGPE04<T>::boundconds_classic(const long double& b1,
		const long double& b2, const long double eps) {
	if ((b1 != 0) && (b2 != 0)) {
		if (abs(b1 - b2) / abs(b1) >= eps)
			return 4;
		else
			return 0;
	} else if (b1 == 0) {
		if (abs(b2) >= eps)
			return 4;
		else
			return 0;
	} else if (b2 == 0) {
		if (abs(b1) >= eps)
			return 4;
		else
			return 0;
	} else
		return 0;

}

} /* namespace interval_arithmetic */
#endif /* EXAMPLEGPE03_H_ */
