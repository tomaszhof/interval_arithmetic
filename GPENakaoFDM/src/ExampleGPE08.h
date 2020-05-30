/*
 * Example09.h
 *
 *  Created on: 30-05-2020
 *      Author: thof
 */

#ifndef EXAMPLEGPE08_H_
#define EXAMPLEGPE08_H_

#include "Utils.h"
//#include "Interval.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class ExampleGPE08: public BoundaryConditions<T> {
	const long double PI = 3.141592653589793238L;
public:
	ExampleGPE08();
	virtual ~ExampleGPE08();
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
ExampleGPE08<T>::ExampleGPE08() {
	// TODO Auto-generated constructor stub

}

template<typename T>
ExampleGPE08<T>::~ExampleGPE08() {
	// TODO Auto-generated destructor stub
}

template<typename T>
long double ExampleGPE08<T>::f(const long double& x, const long double& y) {

	//-(Pi Sin[(Pi x)/2] Sin[Pi y])
	return -PI*sin((PI*x)/2.0)*sin(PI*y);
}

template<typename T>
long double ExampleGPE08<T>::phi1(const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE08<T>::phi2(const long double& x) {
	return 0.0;
}

template<typename T>
long double ExampleGPE08<T>::phi3(const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE08<T>::phi4(const long double& x) {
	return 0.0;
}

template<typename T>
long double ExampleGPE08<T>::c(const long double& x, const long double& y) {
	return (5.0/4.0)*PI*PI;
}

template<typename T>
long double ExampleGPE08<T>::dcdx(const long double& x, const long double& y) {
	//0.0;
	return 0.0;
}

template<typename T>
long double ExampleGPE08<T>::dcdy(const long double& x, const long double& y) {
	//return 0.0
	return 0.0;
}

template<typename T>
long double ExampleGPE08<T>::d2cdx2(const long double& x, const long double& y) {
	//0.0
	return 0.0;
}

template<typename T>
long double ExampleGPE08<T>::d2cdy2(const long double& x, const long double& y) {
	//return 0.0;
	return 0.0;
}

template<typename T>
long double ExampleGPE08<T>::a1(const long double& x, const long double& y) {
	return 1.0;
}

template<typename T>
long double ExampleGPE08<T>::d2a1dx2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE08<T>::d2a1dy2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE08<T>::a2(const long double& x, const long double& y) {
	return 1.0;
}

template<typename T>
long double ExampleGPE08<T>::d2a2dx2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE08<T>::d2a2dy2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
Interval<T> ExampleGPE08<T>::F(const Interval<T>& ix, const Interval<T>& iy,

		int& st) {
	//-PI*sin((PI*x)/2.0)*sin(PI*y);
	Interval<T> r = { 0.0L, 0.0L };
	st = 0;
	r = im1*ipi*ISin((ipi*ix)/i2)*ISin(ipi*iy);
	return r;
}

template<typename T>
Interval<T> ExampleGPE08<T>::PHI1(const Interval<T>& iy, int& st) {
	//phi1(y) = 0.0;
	Interval<T> r = { 0, 0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE08<T>::PHI2(const Interval<T>& ix, int& st) {
	//phi2(x) = 0.0;
	Interval<T> r = { 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE08<T>::PHI3(const Interval<T>& iy, int& st) {
	//phi3(y) = 0.0;
	Interval<T> r = { 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE08<T>::PHI4(const Interval<T>& ix, int& st) {
	//phi3(y) = 0.0;
	Interval<T> r = { 0, 0 };
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE08<T>::PSI(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//psi(x,y) = 0
	Interval<T> r = { 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE08<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//omega(x,y) = 0
	Interval<T> r = { 0, 0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE08<T>::PSI1(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE08<T>::PSI2(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE08<T>::PSI3(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE08<T>::PSI4(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE08<T>::C(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	return (i5/i4) * ipi *ipi;
}

template<typename T>
long double ExampleGPE08<T>::GetConstM() {
	long double constM = 1627;

	return constM;
}

template<typename T>
long double ExampleGPE08<T>::GetConstN() {
	long double constN = 1627;

	return constN;
}


//m=n=50
//P = 3.94667
//Q = 3.94667
//R = 0.973341
//S = 17.8609

template<typename T>
long double ExampleGPE08<T>::GetConstP() {
	long double constP = 15.0L;

	return constP;
}

template<typename T>
long double ExampleGPE08<T>::GetConstQ() {
	long double constQ = 15.0L;

	return constQ;
}

template<typename T>
long double ExampleGPE08<T>::GetConstR() {
	long double constP = 7.5L;

	return constP;
}

template<typename T>
long double ExampleGPE08<T>::GetConstS() {
	long double constQ = 35L;

	return constQ;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::A1(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = i1;
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::A2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = i1;
	st = 0;

	return r;
}



template<typename T>
inline Interval<T> ExampleGPE08<T>::D2A1DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::D2A1DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::D2A2DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::D2A2DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}



template<typename T>
inline Interval<T> ExampleGPE08<T>::DCDX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::DCDY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::DFDX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	//d/dx(sin((π x)/2) sin(π y)) = 1/2 π cos((π x)/2) sin(π y)
	r = (i1/i2)*ipi*ICos((ipi*ix)/i2)*ISin(ipi*iy);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::DFDY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	//d/dy(sin((π x)/2) sin(π y)) = π sin((π x)/2) cos(π y)
	r = ipi*ISin((ipi*ix)/i2)*ICos(ipi*iy);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::D2FDX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	// d^2/dx^2(sin((π x)/2) sin(π y)) = -1/4 π^2 sin((π x)/2) sin(π y)
	r = (im1/i4)*ipi*ipi*ISin((ipi*ix)/i2)*ISin(ipi*iy);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::D2FDY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
    // d^2/dy^2(sin((π x)/2) sin(π y)) = -π^2 sin((π x)/2) sin(π y)
	r = im1*ipi*ipi*ISin((ipi*ix)/i2)*ISin(ipi*iy);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::DA1DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::DA1DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::DA2DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE08<T>::DA2DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}
template<typename T>
void ExampleGPE08<T>::SetArithmeticMode(IAMode mode) {
	Interval<T>::SetMode(mode);
}

template<typename T>
int ExampleGPE08<T>::boundconds_classic(const long double& b1,
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

template<typename T>
long double ExampleGPE08<T>::ExactSol(long double x, long double y) {

	//u=x*cos(pi/2.0 * x)*sin(pi*y)
	return x*cos(PI/2.0 * x)*sin(PI*y);
}

} /* namespace interval_arithmetic */
#endif /* ExampleGPE08_H_ */
