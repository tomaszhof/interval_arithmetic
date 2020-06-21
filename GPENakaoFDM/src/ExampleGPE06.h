/*
 * Example09.h
 *
 *  Created on: 11-11-2014
 *      Author: thof
 */

#ifndef EXAMPLEGPE06_H_
#define EXAMPLEGPE06_H_

#include "Utils.h"
//#include "Interval.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class ExampleGPE06: public BoundaryConditions<T> {
	const long double PI = 3.141592653589793238L;
public:
	ExampleGPE06();
	virtual ~ExampleGPE06();
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
	long double d2a1dx2(const long double& x, const long double& y);
	long double d2a1dy2(const long double& x, const long double& y);
	long double d2a2dx2(const long double& x, const long double& y);
	long double d2a2dy2(const long double& x, const long double& y);
	long double c(const long double& x, const long double& y);
	long double dcdx(const long double& x, const long double& y);
	long double dcdy(const long double& x, const long double& y);

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
	const Interval<long double> i12 = { 12.0, 12.0 };
	const Interval<long double> i15 = { 15.0, 15.0 };
	const Interval<long double> i20 = { 20.0, 20.0 };
	const Interval<long double> i24 = { 24.0, 24.0 };
	const Interval<long double> i36 = { 36.0, 36.0 };
	const Interval<long double> i40 = { 40.0, 40.0 };
	const Interval<long double> i45 = { 45.0, 45.0 };
	const Interval<long double> i120 = { 120.0, 120.0 };
	const Interval<long double> i180 = { 180.0, 180.0 };
	const Interval<long double> i360 = { 360.0, 360.0 };
	const Interval<long double> im1 = { -1.0, -1.0 };
	const Interval<long double> im11 = { -1.0, 1.0 };
	const Interval<long double> ipi = Interval<long double>::IPi();
};

template<typename T>
ExampleGPE06<T>::ExampleGPE06() {
	// TODO Auto-generated constructor stub

}

template<typename T>
ExampleGPE06<T>::~ExampleGPE06() {
	// TODO Auto-generated destructor stub
}

template<typename T>
long double ExampleGPE06<T>::f(const long double& x, const long double& y) {

	//20*(1-x)*(1-y)*(1-%e^(x*y))*sin(%pi*x*y)-(1-x)*(1-y)*y^2*%e^(x*y)+2*(1-y)*y*%e^(x*y)-(1-x)*x^2*(1-y)*%e^(x*y)+2*(1-x)*x*%e^(x*y)

	return 20*(1-x)*(1-y)*(1-exp(x*y))*sin(PI*x*y)-(1-x)*(1-y)*pow(y,2.0)*exp(x*y)+2*(1-y)*y*exp(x*y)-(1-x)*pow(x,2.)*(1-y)*exp(x*y)+2*(1-x)*x*exp(x*y);
}

template<typename T>
long double ExampleGPE06<T>::phi1(const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE06<T>::phi2(const long double& x) {
	return 0.0;
}

template<typename T>
long double ExampleGPE06<T>::phi3(const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE06<T>::phi4(const long double& x) {
	return 0.0;
}

template<typename T>
long double ExampleGPE06<T>::c(const long double& x, const long double& y) {
	return 20*sin(PI*x*y);
}

template<typename T>
long double ExampleGPE06<T>::dcdx(const long double& x, const long double& y) {
	//20*%pi*y*cos(%pi*x*y)
	return 20*PI*y*cos(PI*x*y);
}

template<typename T>
long double ExampleGPE06<T>::dcdy(const long double& x, const long double& y) {
	//20*%pi*y*cos(%pi*x*y)
	return 20*PI*x*cos(PI*x*y);
}

template<typename T>
long double ExampleGPE06<T>::a1(const long double& x, const long double& y) {
	return 1.0;
}

template<typename T>
long double ExampleGPE06<T>::d2a1dx2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE06<T>::d2a1dy2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE06<T>::a2(const long double& x, const long double& y) {
	return 1.0;
}

template<typename T>
long double ExampleGPE06<T>::d2a2dx2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE06<T>::d2a2dy2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
Interval<T> ExampleGPE06<T>::F(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//f(x,y)= 20*(1-x)*(1-y)*(1-e^(x*y))*sin(pi*x*y)-(1-x)*(1-y)*y^2*e^(x*y)+2*(1-y)*y*e^(x*y)-(1-x)*x^2*(1-y)*e^(x*y)+2*(1-x)*x*e^(x*y)
	Interval<T> r = { 0.0L, 0.0L };
	Interval<T> iexpXY = IExp(ix * iy);
	Interval<T> isinPIXY = ISin(ipi*ix*iy);
	r = i20 * (i1 - ix)*(i1 - iy)*(i1-iexpXY)*isinPIXY - (i1 -ix)*(i1-iy)*iy*iy*iexpXY + i2*(i1-iy)*iy*iexpXY-(i1-ix)*ix*ix*(i1-iy)*iexpXY+i2*(i1-ix)*ix*iexpXY;
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE06<T>::PHI1(const Interval<T>& iy, int& st) {
	//phi1(y) = 0.0;
	Interval<T> r = { 0, 0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE06<T>::PHI2(const Interval<T>& ix, int& st) {
	//phi2(x) = 0.0;
	Interval<T> r = { 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE06<T>::PHI3(const Interval<T>& iy, int& st) {
	//phi3(y) = 0.0;
	Interval<T> r = { 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE06<T>::PHI4(const Interval<T>& ix, int& st) {
	//phi3(y) = 0.0;
	Interval<T> r = { 0, 0 };
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE06<T>::PSI(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//psi(x,y) = 0
	Interval<T> r = { 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE06<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//omega(x,y) = 0
	Interval<T> r = { 0, 0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE06<T>::PSI1(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE06<T>::PSI2(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE06<T>::PSI3(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE06<T>::PSI4(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE06<T>::C(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//c(x,y) = 20*sin(PI*x*y);

	return i20 * ISin(ipi*ix*iy);
}

template<typename T>
long double ExampleGPE06<T>::ExactSol(long double x, long double y) {

	//u=(1-x)*(1-y)*(1-%e^(x*y))
	return (1-x)*(1-y)*(1-exp(x*y));
}

template<typename T>
long double ExampleGPE06<T>::GetConstM() {
	long double constM = 1627;

	return constM;
}

template<typename T>
long double ExampleGPE06<T>::GetConstN() {
	long double constN = 1627;

	return constN;
}


//m=n=50
//P = 3.94667
//Q = 3.94667
//R = 0.973341
//S = 17.8609

template<typename T>
long double ExampleGPE06<T>::GetConstP() {
	long double constP = 4.0L;

	return constP;
}

template<typename T>
long double ExampleGPE06<T>::GetConstQ() {
	long double constQ = 4.0L;

	return constQ;
}

template<typename T>
long double ExampleGPE06<T>::GetConstR() {
	long double constP = 1.0L;

	return constP;
}

template<typename T>
long double ExampleGPE06<T>::GetConstS() {
	long double constQ = 18L;

	return constQ;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::A1(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = i1;
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::A2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = i1;
	st = 0;

	return r;
}



template<typename T>
inline Interval<T> ExampleGPE06<T>::D2A1DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::D2A1DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::D2A2DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::D2A2DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}



template<typename T>
inline Interval<T> ExampleGPE06<T>::DCDX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	//20*PI*y*cos(PI*x*y)
	Interval<T> r = { 0, 0 };
	st = 0;
	r = i20 * ipi * iy * ICos(ipi*ix*iy);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::DCDY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	r = i20 * ipi * ix * ICos(ipi*ix*iy);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::DFDX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	//e^(x y) (-x^3 (y - 1) y + x^2 (y^2 - 6 y + 3) + x (-y^4 + y^3 + 4 y - 6) + y^4 - 4 y^3 + 3 y^2 + 2)
	//- 20 (y - 1) (e^(x y) ((x - 1) y + 1) - 1) sin(π x y) - 20 π (x - 1) (y - 1) y (e^(x y) - 1) cos(π x y)

	Interval<T> iypow4 = iy*iy*iy*iy;
	Interval<T> iypow3 = iy*iy*iy;
	Interval<T> ixpow3 = ix*ix*ix;
	Interval<T> iexpXY = IExp(ix*iy);

	r = iexpXY*(ixpow3*(iy-i1)*iy + ix*ix*(iy*iy-i6*iy+i3) + ix*(im1*iypow4 + iypow3 + i4*iy - i6) + iypow4 -i4*iypow3 + i3*iy*iy + i2);
	r = r - i20*(iy-i1)*(iexpXY *((ix-i1)*iy + i1)-i1)*ISin(ipi*ix*iy) - i20*ipi*(ix-i1)*(iy-i1)*iy*(iexpXY - i1)*ICos(ipi*ix*iy);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::DFDY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	Interval<T> ixpow4 = ix*ix*ix*ix;
	Interval<T> ixpow3 = ix*ix*ix;
	Interval<T> iypow3 = iy*iy*iy;
	Interval<T> iexpXY = IExp(iy*ix);

	r = iexpXY*(iypow3*(ix-i1)*ix + iy*iy*(ix*ix-i6*ix+i3) + iy*(im1*ixpow4 + ixpow3 + i4*ix - i6) + ixpow4 -i4*ixpow3 + i3*ix*ix + i2);
	r = r - i20*(ix-i1)*(iexpXY *((iy-i1)*ix + i1)-i1)*ISin(ipi*iy*ix) - i20*ipi*(iy-i1)*(ix-i1)*ix*(iexpXY - i1)*ICos(ipi*iy*ix);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::D2FDX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	//-x^2 (y - 1) y e^(x y) ((x - 1) y + 2) - 2 (y - 1) y^3 e^(x y) - (y - 1) y^3 e^(x y) ((x - 1) y + 2) - 2 (x - 1) x y^2 e^(x y)
	//+ 20 π^2 (x - 1) (y - 1) y^2 (e^(x y) - 1) sin(π x y) + 4 (1 - 2 x) y e^(x y) - 4 e^(x y) - 2 (x - 1) (y - 1) e^(x y)
	//- 4 x (y - 1) e^(x y) ((x - 1) y + 1) - 20 (y - 1) y e^(x y) ((x - 1) y + 2) sin(π x y) - 40 π (y - 1) y (e^(x y) ((x - 1) y + 1) - 1) cos(π x y)
	Interval<T> r = { 0, 0 };
	Interval<T> iexpXY = IExp(ix * iy);
	Interval<T> isinPIXY = ISin(ipi*ix*iy);
	st = 0;
	r = im1*(ix*ix) *(iy - i1) * iy * iexpXY *((ix - i1)*iy +i2) -i2*(iy-i1)*(iy*iy*iy)*iexpXY - (iy - i1)*(iy *iy *iy) * iexpXY *((ix -i1)*iy + i2) - i2 * (ix - i1)*ix*iy*iy*iexpXY;
	r = r + i20 * ipi *ipi * (ix -i1)*(iy-i1)*iy*iy*(iexpXY - i1) * isinPIXY + i4 *(i1 - i2 *ix)*iy*iexpXY - i4*iexpXY -i2*(ix-i1)*iexpXY;
	r = r  - i4*ix*(iy-i1)*iexpXY*((ix-i1)*iy + i1) - i20*(iy-i1)*iy*iexpXY*((ix-i1)*iy + i2)*isinPIXY - i40*ipi*(iy-i1)*iy*(iexpXY*((ix-i1)+i1) - i1) * ICos(ipi*ix*iy);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::D2FDY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {

    //(x - 1) x^4 (y - 1) (-e^(x y)) - 4 (x - 1) x^3 e^(x y) - (x - 1) x^2 (y - 1) y^2 e^(x y) - 2 x^2 (y - 1) y e^(x y)
	//- 20 (x - 1) x^2 (y - 1) e^(x y) sin(π x y) + 4 x (1 - 2 y) e^(x y) - 2 (x - 1) x y (3 y - 2) e^(x y) - 4 e^(x y)
	//- 2 (x - 1) (3 y - 1) e^(x y) - 40 (x - 1) x e^(x y) (sin(π x y) + π x (y - 1) cos(π x y))
	//+ 20 π (x - 1) x (e^(x y) - 1) (π x (y - 1) sin(π x y) - 2 cos(π x y))

	Interval<T> r = { 0, 0 };
	Interval<T> iexpXY = IExp(ix * iy);
	Interval<T> isinPIXY = ISin(ipi*ix*iy);
	Interval<T> icosPIXY = ICos(ipi*ix*iy);
	st = 0;
	r = im1*(iy*iy) *(ix - i1) * ix * iexpXY *((iy - i1)*ix +i2) -i2*(ix-i1)*(ix*ix*ix)*iexpXY - (ix - i1)*(ix *ix *ix) * iexpXY *((iy -i1)*ix + i2) - i2 * (iy - i1)*iy*ix*ix*iexpXY;
	r = r + i20 * ipi *ipi * (iy -i1)*(ix-i1)*ix*ix*(iexpXY - i1) * isinPIXY + i4 *(i1 - i2 *iy)*ix*iexpXY - i4*iexpXY -i2*(iy-i1)*iexpXY;
	r = r  - i4*iy*(ix-i1)*iexpXY*((iy-i1)*ix + i1) - i20*(ix-i1)*ix*iexpXY*((iy-i1)*ix + i2)*isinPIXY - i40*ipi*(ix-i1)*ix*(iexpXY*((iy-i1)+i1) - i1) * ICos(ipi*iy*ix);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::DA1DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::DA1DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::DA2DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE06<T>::DA2DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}
template<typename T>
void ExampleGPE06<T>::SetArithmeticMode(IAMode mode) {
	Interval<T>::SetMode(mode);
}

template<typename T>
int ExampleGPE06<T>::boundconds_classic(const long double& b1,
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
#endif /* EXAMPLEGPE06_H_ */
