/*
 * ExampleGPE05.h
 *
 *  Created on: 31.03.2019
 *      Author: numeric
 */

#ifndef EXAMPLEGPE05_H_
#define EXAMPLEGPE05_H_
#include "Utils.h"
#include "BoundaryConditions.h"
using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class ExampleGPE05 : public BoundaryConditions<T>{
	const long double PI = 3.141592653589793238L;
public:
	ExampleGPE05();
	virtual ~ExampleGPE05();
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
ExampleGPE05<T>::ExampleGPE05() {
	// TODO Auto-generated constructor stub

}

template<typename T>
ExampleGPE05<T>::~ExampleGPE05() {
	// TODO Auto-generated destructor stub
}

template<typename T>
long double ExampleGPE05<T>::f(const long double& x, const long double& y) {
    //sin(π x) sin(π y) (-π^2 cos(π x^2 y) - π^2 cos(π x y^2) + x y)
	return sin(PI*x) * sin(PI*y) * (-PI*PI * cos(PI*PI*y) - PI*PI*cos(PI*x*y*y) + 1.0L); //+ x*y);
}

template<typename T>
long double ExampleGPE05<T>::phi1(const long double& y) {
	return 0.0L;
}

template<typename T>
long double ExampleGPE05<T>::phi2(const long double& x) {
	return 0.0L;
}

template<typename T>
long double ExampleGPE05<T>::phi3(const long double& y) {
	return sin(PI*y);
}

template<typename T>
long double ExampleGPE05<T>::phi4(const long double& x) {
	return sin(PI*x);
}

template<typename T>
long double ExampleGPE05<T>::c(const long double& x, const long double& y) {
	return 1.0L;
}

template<typename T>
long double ExampleGPE05<T>::dcdx(const long double& x, const long double& y) {
	return 0.0L;
}

template<typename T>
long double ExampleGPE05<T>::dcdy(const long double& x, const long double& y) {
	return 0.0L;
}

template<typename T>
long double ExampleGPE05<T>::a1(const long double& x, const long double& y) {
	return cos(PI*x*x*y);
}

template<typename T>
long double ExampleGPE05<T>::d2a1dx2(const long double& x,
		const long double& y) {
	//-2 π y (2 π x^2 y cos(π x^2 y) + sin(π x^2 y))
	return -2*PI*y*(2*PI*x*x*y*cos(PI*x*x*y) + sin(PI*x*x*y));
}

template<typename T>
long double ExampleGPE05<T>::d2a1dy2(const long double& x,
		const long double& y) {
	//-π^2 x^4 cos(π x^2 y)
	return -PI*PI*(x*x*x*x)*sin(PI*x*x*y);
}

template<typename T>
long double ExampleGPE05<T>::a2(const long double& x, const long double& y) {
	return cos(PI*x*y*y);
}

template<typename T>
long double ExampleGPE05<T>::d2a2dx2(const long double& x,
		const long double& y) {
	//-π^2 y^4 cos(π x y^2)
	return -PI*PI*(y*y*y*y)*sin(PI*x*y*y);
}

template<typename T>
long double ExampleGPE05<T>::d2a2dy2(const long double& x,
		const long double& y) {
	//-2 π x (2 π y^2 x cos(π y^2 x) + sin(π y^2 x))
	return -2*PI*x*(2*PI*y*y*x*cos(PI*y*y*x) + sin(PI*y*y*x));
}

template<typename T>
Interval<T> ExampleGPE05<T>::F(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//f(x,y)= x*exp(x*y)*(x*x*x*y*sin(PI*y) + x*x*y*y*sin(PI*x) + 2*x*x*sin(PI*y) - y);
	Interval<T> r = { 0.0L, 0.0L };
	r = (ix * IExp(ix * iy)) * (ix*ix*ix*iy*ISin(ipi * iy) + ix*ix*iy*iy*ISin(ipi*ix) + i2*ix*ix*ISin(ipi*iy) - iy);
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE05<T>::PHI1(const Interval<T>& iy, int& st) {
	//phi1(y) = 0.0;
	Interval<T> r = { 1.0, 1.0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE05<T>::PHI2(const Interval<T>& ix, int& st) {
	//phi2(x) = 0.0;
	Interval<T> r = { 1.0, 1.0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE05<T>::PHI3(const Interval<T>& iy, int& st) {
	//phi3(y) = 0.0;
	Interval<T> r = { 0, 0 };

	r = IExp(iy);
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE05<T>::PHI4(const Interval<T>& ix, int& st) {
	//phi3(y) = 0.0;
	Interval<T> r = { 0, 0 };
	st = 0;
	r = IExp(ix);
	return r;
}

template<typename T>
Interval<T> ExampleGPE05<T>::PSI(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//psi(x,y) = 0
	Interval<T> r = { 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE05<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//omega(x,y) = 0
	Interval<T> r = { 0, 0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE05<T>::PSI1(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE05<T>::PSI2(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE05<T>::PSI3(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE05<T>::PSI4(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE05<T>::C(const Interval<T>& ix, const Interval<T>& iy,
		int& st) {
	//c(x,y) = x * y;

	return ix * iy;
}

template<typename T>
long double ExampleGPE05<T>::ExactSol(long double x, long double y) {
	return sin(PI*x)*sin(PI*y);
}

template<typename T>
long double ExampleGPE05<T>::GetConstM() {
	long double constM = 0.0;

	return constM;
}

template<typename T>
long double ExampleGPE05<T>::GetConstN() {
	long double constN = 0.0;

	return constN;
}

template<typename T>
long double ExampleGPE05<T>::GetConstP() {
	long double constP = 80.0L;

	return constP;
}

template<typename T>
long double ExampleGPE05<T>::GetConstQ() {
	long double constQ = 400.0L;

	return constQ;
}

template<typename T>
long double ExampleGPE05<T>::GetConstR() {
	long double constP = 4500.0L;

	return constP;
}

template<typename T>
long double ExampleGPE05<T>::GetConstS() {
	long double constQ = 200L;

	return constQ;
}

template<typename T>
inline Interval<T> ExampleGPE05<T>::A1(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {

	//a1(x,y) = (x*x)*sin(PI*y);
	Interval<T> r = {0.0L, 0.0L};
	st = 0;

	r = ix * ix * ISin(PI * iy);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE05<T>::A2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	//a2(x,y) = (y*y)*sin(PI*x);
	Interval<T> r = {0.0L, 0.0L};
	st = 0;
	r = iy * iy * ISin(PI*ix);

	return r;
}



template<typename T>
inline Interval<T> ExampleGPE05<T>::D2A1DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	//2.0L*sin(PI*y);
	Interval<T> r = { 0, 0 };
	st = 0;
	r = i2 * ISin(ipi * iy);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE05<T>::D2A1DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	//-PI*PI*(x*x)*sin(PI*y);
	Interval<T> r = { 0, 0 };
	st = 0;
	r = im1 * ipi * ipi * ix * ix * ISin(ipi * iy);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE05<T>::D2A2DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	//-PI*PI*(y*y)*sin(PI*x);
	Interval<T> r = { 0, 0 };
	st = 0;

	r = im1 * ipi * ipi * iy * iy * ISin(ipi * ix);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE05<T>::D2A2DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	//2.0L*sin(PI*x);
	Interval<T> r = { 0, 0 };
	st = 0;

	r = i2 * ISin(ipi * ix);
	return r;
}



template<typename T>
inline Interval<T> ExampleGPE05<T>::DCDX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	r = im1 * iy;
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE05<T>::DCDY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	r = im1 * ix;
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE05<T>::D2FDX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	r = ix*ix*ix*ix*iy*iy*iy*ISin(ipi * iy);
	r = r + ix*ix*ix*iy*iy*iy*iy*ISin(ipi * ix);
	r = r + i2*ipi*ix*ix*ix*iy*iy*iy*ICos(ipi * ix);
	r = r - ipi*ipi*ix*ix*ix*iy*iy*ISin(ipi*ix);
	r = r + i10*ix*ix*ix*iy*iy*ISin(ipi*iy) + i6*ix*ix*iy*iy*iy*ISin(ipi * ix);
	r = r + i6*ipi*ix*ix*iy*iy*ICos(ipi*ix);
	r = r + i24*ix*ix*iy*ISin(ipi*iy);
	r = r - ix*iy*iy*iy + i6*ix*iy*iy*ISin(ipi*ix);
	r = r + i12*ix*ISin(ipi*iy);
	r = r - i2*iy*iy;
	r = r * IExp(ix*iy);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE05<T>::D2FDY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	r = ix*ix*ix*ix*iy*ISin(ipi*iy)+ix*ix*ix*iy*iy*ISin(ipi*ix)+i4*ix*ix*ix*ISin(ipi*iy)+i2*ipi*ix*ix*ix*iy*ICos(ipi*iy);
	r = r + i4*ix*ix*iy*ISin(ipi*ix) - ipi*ipi*ix*ix*iy*ISin(ipi*iy) + i6*ipi*ix*ix*ICos(ipi*iy) - ix*iy;
	r = r - i2*ipi*ipi*ix*ISin(ipi*iy)+i2*ix*ISin(ipi*ix)-i2;
	r = r * ix * ix * IExp(ix*iy);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE05<T>::DA1DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	r = i2 * ix * ISin(PI*iy);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE05<T>::DA1DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	r = ipi * ix * ix * ICos(ipi*ix);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE05<T>::DA2DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	r = ipi * iy * iy * ICos(ipi*ix);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE05<T>::DA2DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	r = i2 * iy * ISin(PI*ix);
	return r;
}
template<typename T>
void ExampleGPE05<T>::SetArithmeticMode(IAMode mode) {
	Interval<T>::SetMode(mode);
}

template<typename T>
int ExampleGPE05<T>::boundconds_classic(const long double& b1,
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
#endif /* EXAMPLEGPE05_H_ */
