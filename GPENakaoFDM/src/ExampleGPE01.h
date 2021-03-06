/*
 * Example09.h
 *
 *  Created on: 11-11-2014
 *      Author: thof
 */

#ifndef EXAMPLEGPE01_H_
#define EXAMPLEGPE01_H_

#include "Utils.h"
//#include "Interval.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class ExampleGPE01: public BoundaryConditions<T> {
	const Interval<T> i0 = { 0.0L, 0.0L };
	const Interval<T> i1 = { 1.0L, 1.0L };

public:
	ExampleGPE01();
	virtual ~ExampleGPE01();
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
	virtual long double a1(const long double& x, const long double& y);
	virtual long double a2(const long double& x, const long double& y);
	virtual long double d2a1dx2(const long double& x, const long double& y);
	virtual long double d2a1dy2(const long double& x, const long double& y);
	virtual long double d2a2dx2(const long double& x, const long double& y);
	virtual long double d2a2dy2(const long double& x, const long double& y);
	virtual long double c(const long double& x, const long double& y);
	virtual long double dcdx(const long double& x, const long double& y);
	virtual long double dcdy(const long double& x, const long double& y);

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
};

template<typename T>
ExampleGPE01<T>::ExampleGPE01() {
	// TODO Auto-generated constructor stub

}

template<typename T>
ExampleGPE01<T>::~ExampleGPE01() {
	// TODO Auto-generated destructor stub
}

template<typename T>
long double ExampleGPE01<T>::f(const long double& x, const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE01<T>::phi1(const long double& y) {
	return cos(3.0 * y);
}

template<typename T>
long double ExampleGPE01<T>::phi2(const long double& x) {
	return exp(3.0 * x);
}

template<typename T>
long double ExampleGPE01<T>::phi3(const long double& y) {
	return exp(3.0) * cos(3.0 * y);
}

template<typename T>
long double ExampleGPE01<T>::phi4(const long double& x) {
	return exp(3.0 * x) * cos(3.0);
}

template<typename T>
long double ExampleGPE01<T>::c(const long double& x, const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE01<T>::dcdx(const long double& x, const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE01<T>::dcdy(const long double& x, const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE01<T>::a1(const long double& x, const long double& y) {
	return 1.0;
}

template<typename T>
long double ExampleGPE01<T>::d2a1dx2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE01<T>::d2a1dy2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE01<T>::a2(const long double& x, const long double& y) {
	return 1.0;
}

template<typename T>
long double ExampleGPE01<T>::d2a2dx2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE01<T>::d2a2dy2(const long double& x,
		const long double& y) {
	return 0.0;
}

template<typename T>
Interval<T> ExampleGPE01<T>::F(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//f(x,y)= 0;
	Interval<T> r;
	r.a = 0;
	r.b = 0;
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE01<T>::PHI1(const Interval<T>& iy, int& st) {
	//phi1(y) = cos(3*y);
	Interval<T> r = { 0, 0 };
	Interval<T> ithree = { 3, 3 };

	st = 0;

	r = ICos(ithree * iy);
	return r;
}

template<typename T>
Interval<T> ExampleGPE01<T>::PHI2(const Interval<T>& ix, int& st) {
	//phi2(x) = exp(3*x);
	Interval<T> r = { 0, 0 };
	Interval<T> ithree = { 3, 3 };

	st = 0;

	r = IExp(ithree * ix);
	return r;
	return r;
}

template<typename T>
Interval<T> ExampleGPE01<T>::PHI3(const Interval<T>& iy, int& st) {
	//phi3(y) = exp(3.0) * cos(3.0*y);
	Interval<T> r = { 0, 0 };
	Interval<T> ithree = { 3, 3 };

	st = 0;

	r = IExp(ithree) * ICos(ithree * iy);
	return r;
}

template<typename T>
Interval<T> ExampleGPE01<T>::PHI4(const Interval<T>& ix, int& st) {
	//phi3(y) = exp(3.0*x) * cos(3.0);
	Interval<T> r = { 0, 0 };
	Interval<T> ithree = { 3, 3 };
	st = 0;
	r = IExp(ithree * ix) * ICos(ithree);
	return r;
}

template<typename T>
Interval<T> ExampleGPE01<T>::PSI(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//psi(x,y) = 0

	Interval<T> r = { 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE01<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//omega(x,y) = 0
	Interval<T> r = { 0, 0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE01<T>::PSI1(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE01<T>::PSI2(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE01<T>::PSI3(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE01<T>::PSI4(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE01<T>::C(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//c(x,y) = 1.0;
	Interval<T> ione = { 1, 1 };

	return ione;
}

template<typename T>
long double ExampleGPE01<T>::ExactSol(long double x, long double y) {
	return exp(3.0 * x) * cos(3.0 * y);
}

template<typename T>
long double ExampleGPE01<T>::GetConstM() {
	long double constM = 1627;

	return constM;
}

template<typename T>
long double ExampleGPE01<T>::GetConstN() {
	long double constN = 1627;

	return constN;
}

template<typename T>
long double ExampleGPE01<T>::GetConstP() {
	long double constP = 14643;

	return constP;
}

template<typename T>
long double ExampleGPE01<T>::GetConstQ() {
	long double constQ = 14643;

	return constQ;
}

template<typename T>
long double ExampleGPE01<T>::GetConstR() {
	long double constR = 0;

	return constR;
}

template<typename T>
long double ExampleGPE01<T>::GetConstS() {
	long double constS = 0;

	return constS;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::A1(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::A2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::D2A1DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::DCDX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::DCDY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::D2A1DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::D2A2DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::D2A2DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::D2FDY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::D2FDX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::DA1DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::DA1DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::DA2DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE01<T>::DA2DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
void ExampleGPE01<T>::SetArithmeticMode(IAMode mode) {
	Interval<T>::SetMode(mode);
}

template<typename T>
int ExampleGPE01<T>::boundconds_classic(const long double& b1,
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
