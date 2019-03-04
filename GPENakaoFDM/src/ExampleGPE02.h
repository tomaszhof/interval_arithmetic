/*
 * Example09.h
 *
 *  Created on: 11-11-2014
 *      Author: thof
 */

#ifndef EXAMPLEGPE02_H_
#define EXAMPLEGPE02_H_

#include "Utils.h"
//#include "Interval.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class ExampleGPE02: public BoundaryConditions<T> {
	const long double PI = 3.141592653589793238L;
public:
	ExampleGPE02();
	virtual ~ExampleGPE02();
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
	const Interval<long double> i45 = { 45.0, 45.0 };
	const Interval<long double> i120 = { 120.0, 120.0 };
	const Interval<long double> i180 = { 180.0, 180.0 };
	const Interval<long double> i360 = { 360.0, 360.0 };
	const Interval<long double> im1 = { -1.0, -1.0 };
	const Interval<long double> im11 = { -1.0, 1.0 };
	const Interval<long double> ipi = Interval<long double>::IPi();
};

template<typename T>
ExampleGPE02<T>::ExampleGPE02() {
	// TODO Auto-generated constructor stub

}

template<typename T>
ExampleGPE02<T>::~ExampleGPE02() {
	// TODO Auto-generated destructor stub
}

template<typename T>
long double ExampleGPE02<T>::f(const long double& x, const long double& y) {
	return (2.0 * PI - 1.0) * sin(PI * x) * sin(PI * y);
}

template<typename T>
long double ExampleGPE02<T>::phi1(const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::phi2(const long double& x) {
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::phi3(const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::phi4(const long double& x) {
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::c(const long double& x, const long double& y) {
	return -PI;
}

template<typename T>
long double ExampleGPE02<T>::dcdx(const long double& x, const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::dcdy(const long double& x, const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::a1(const long double& x, const long double& y) {
	return -1.0;
}

template<typename T>
long double ExampleGPE02<T>::d2a1dx2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::d2a1dy2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::a2(const long double& x, const long double& y) {
	return -1.0;
}

template<typename T>
long double ExampleGPE02<T>::d2a2dx2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::d2a2dy2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
Interval<T> ExampleGPE02<T>::F(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//f(x,y)= (2.0 * PI - 1.0) *sin(PI * x)*sin(PI * y);
	Interval<T> r = { 0.0L, 0.0L };
	r = (i2 * ipi - i1) * ISin(ipi * ix) * ISin(ipi * iy);
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PHI1(const Interval<T>& iy, int& st) {
	//phi1(y) = 0.0;
	Interval<T> r = { 0, 0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PHI2(const Interval<T>& ix, int& st) {
	//phi2(x) = 0.0;
	Interval<T> r = { 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PHI3(const Interval<T>& iy, int& st) {
	//phi3(y) = 0.0;
	Interval<T> r = { 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PHI4(const Interval<T>& ix, int& st) {
	//phi3(y) = 0.0;
	Interval<T> r = { 0, 0 };
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PSI(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//psi(x,y) = 0
	Interval<T> r = { 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//omega(x,y) = 0
	Interval<T> r = { 0, 0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PSI1(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PSI2(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PSI3(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PSI4(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::C(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//c(x,y) = -PI;

	return im1 * ipi;
}

template<typename T>
long double ExampleGPE02<T>::ExactSol(long double x, long double y) {
	return (1.0 / PI) * sin(PI * x) * sin(PI * y);
}

template<typename T>
long double ExampleGPE02<T>::GetConstM() {
	long double constM = 1627;

	return constM;
}

template<typename T>
long double ExampleGPE02<T>::GetConstN() {
	long double constN = 1627;

	return constN;
}

template<typename T>
long double ExampleGPE02<T>::GetConstP() {
	long double constP = 10.0L;

	return constP;
}

template<typename T>
long double ExampleGPE02<T>::GetConstQ() {
	long double constQ = 10.0L;

	return constQ;
}

template<typename T>
long double ExampleGPE02<T>::GetConstR() {
	long double constP = 31.0L;

	return constP;
}

template<typename T>
long double ExampleGPE02<T>::GetConstS() {
	long double constQ = 10L;

	return constQ;
}

template<typename T>
inline Interval<T> ExampleGPE02<T>::A1(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = im1;
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE02<T>::A2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = im1;
	st = 0;

	return r;
}



template<typename T>
inline Interval<T> ExampleGPE02<T>::D2A1DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE02<T>::D2A1DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE02<T>::D2A2DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE02<T>::D2A2DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}



template<typename T>
inline Interval<T> ExampleGPE02<T>::DCDX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE02<T>::DCDY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE02<T>::D2FDX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	r = im1*(ipi*ipi) *(i2 * ipi - i1) * ISin(ipi * ix) * ISin(ipi * iy);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE02<T>::D2FDY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	r = im1*(ipi*ipi) *(i2 * ipi - i1) * ISin(ipi * ix) * ISin(ipi * iy);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE02<T>::DA1DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE02<T>::DA1DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE02<T>::DA2DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE02<T>::DA2DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}
template<typename T>
void ExampleGPE02<T>::SetArithmeticMode(IAMode mode) {
	Interval<T>::SetMode(mode);
}

template<typename T>
int ExampleGPE02<T>::boundconds_classic(const long double& b1,
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
#endif /* EXAMPLEGPE01_H_ */
