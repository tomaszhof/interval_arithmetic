/*
 * Example09.h
 *
 *  Created on: 11-11-2014
 *      Author: thof
 */

#ifndef EXAMPLEGPE11_H_
#define EXAMPLEGPE11_H_

#include "Utils.h"
//#include "Interval.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class ExampleGPE11: public BoundaryConditions<T> {
	const Interval<T> i0 = { 0.0L, 0.0L };
	const Interval<T> i1 = { 1.0L, 1.0L };

public:
	ExampleGPE11();
	virtual ~ExampleGPE11();
	Interval<T> F(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> PHI1(const Interval<T> &iy, int &st);
	Interval<T> PHI2(const Interval<T> &ix, int &st);
	Interval<T> PHI3(const Interval<T> &iy, int &st);
	Interval<T> PHI4(const Interval<T> &ix, int &st);
	Interval<T> PSI(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> OMEGA(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> PSI1(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> PSI2(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> PSI3(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> PSI4(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> A1(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> A2(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> DA1DX(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> DA1DY(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> DA2DX(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> DA2DY(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> D2A1DX2(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> D2A1DY2(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> D2A2DX2(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> D2A2DY2(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> C(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> DCDX(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> DCDY(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> D2FDX2(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> D2FDY2(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> DFDX(const Interval<T> &ix, const Interval<T> &iy, int &st);
	Interval<T> DFDY(const Interval<T> &ix, const Interval<T> &iy, int &st);

	//classical functions
	long double f(const long double &x, const long double &y);
	long double phi1(const long double &y);
	long double phi2(const long double &x);
	long double phi3(const long double &y);
	long double phi4(const long double &x);
	virtual long double a1(const long double &x, const long double &y);
	virtual long double a2(const long double &x, const long double &y);
	virtual long double d2a1dx2(const long double &x, const long double &y);
	virtual long double d2a1dy2(const long double &x, const long double &y);
	virtual long double d2a2dx2(const long double &x, const long double &y);
	virtual long double d2a2dy2(const long double &x, const long double &y);
	virtual long double c(const long double &x, const long double &y);
	virtual long double dcdx(const long double &x, const long double &y);
	virtual long double dcdy(const long double &x, const long double &y);

	int boundconds_classic(const long double &b1, const long double &b2,
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
ExampleGPE11<T>::ExampleGPE11() {
	// TODO Auto-generated constructor stub

}

template<typename T>
ExampleGPE11<T>::~ExampleGPE11() {
	// TODO Auto-generated destructor stub
}

template<typename T>
long double ExampleGPE11<T>::f(const long double &x, const long double &y) {
	long double sinXsinY = std::sin(M_PI * x) * std::sin(M_PI * y);
	return -1.0 * M_PI * M_PI * x * y
			* (std::exp(x) * sinXsinY + std::exp(y) * sinXsinY);
}

template<typename T>
long double ExampleGPE11<T>::phi1(const long double &y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE11<T>::phi2(const long double &x) {
	return 0.0;
}

template<typename T>
long double ExampleGPE11<T>::phi3(const long double &y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE11<T>::phi4(const long double &x) {
	return 0.0;
}

template<typename T>
long double ExampleGPE11<T>::c(const long double &x, const long double &y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE11<T>::dcdx(const long double &x, const long double &y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE11<T>::dcdy(const long double &x, const long double &y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE11<T>::a1(const long double &x, const long double &y) {
	return x * y * std::exp(y);
}

template<typename T>
long double ExampleGPE11<T>::d2a1dx2(const long double &x,
		const long double &y) {
	return 0.0;
}

template<typename T>
long double ExampleGPE11<T>::d2a1dy2(const long double &x,
		const long double &y) {
	long double expy = std::exp(y);
	return x * (expy * y + 2.0 * expy);
}

template<typename T>
long double ExampleGPE11<T>::a2(const long double &x, const long double &y) {
	return x * y * std::exp(x);
}

template<typename T>
long double ExampleGPE11<T>::d2a2dx2(const long double &x,
		const long double &y) {
	long double expx = std::exp(x);
	return y * (expx * x + 2.0 * expx);
}

template<typename T>
long double ExampleGPE11<T>::d2a2dy2(const long double &x,
		const long double &y) {
	return 0.0;
}

template<typename T>
Interval<T> ExampleGPE11<T>::F(const Interval<T> &ix, const Interval<T> &iy,
		int &st) {
	//f(x,y)=-1.0*M_PI*M_PI*x*y*(std::exp(x)*sinXsinY + std::exp(y)*sinXsinY);

	Interval<T> r, sinx, siny, iexpx, iexpy;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> imone(-1, -1);
	Interval<T> ipi = Interval<T>::IPi();

	if (Interval<T>::GetMode() == PINT_MODE) {
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		iexpx = IExp(ix);
		iexpy = IExp(iy);
	} else {
		sinx = DISin(ipi * ix);
		siny = DISin(ipi * iy);
		iexpx = DIExp(ix);
		iexpy = DIExp(iy);
	}

	r = imone * ipi * ipi * ix * iy
			* (iexpx * sinx * siny + iexpy * sinx * siny);
	return r;
}

template<typename T>
Interval<T> ExampleGPE11<T>::PHI1(const Interval<T> &iy, int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE11<T>::PHI2(const Interval<T> &ix, int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE11<T>::PHI3(const Interval<T> &iy, int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE11<T>::PHI4(const Interval<T> &ix, int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE11<T>::PSI(const Interval<T> &ix, const Interval<T> &iy,
		int &st) {
	//	d2f/dx2 = -\pi^2*y*(2e^x*sin(pi*x)*sin(pi*y)+2*pi*e^x*cos(pi*x)*sin(pi*y)
	//						 +2*pi*e^y*cos(pi*x)*sin(pi*y)+e^x*x*sin(pi*x)*sin(pi*y)
	//						 +2*pi*e^x*x*cos(pi*x)*sin(pi*y)
	//						 -pi^2*e^x*x*sin(pi*x)*sin(pi*y)
	//						 -pi^2*e^y*x*sin(pi x)*sin(pi*y))

	Interval<T> r, sinx, siny, cosx, cosy, iex, iey, impi2y, ip1, ip2, ip3, ip4,
			ip5;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> itwo(2, 2);
	Interval<T> im1(-1, -1);
	Interval<T> ipi = Interval<T>::IPi();
	Interval<T> ipi2 = ipi * ipi;
	Interval<T> ipi3 = ipi * ipi * ipi;

	if (Interval<T>::GetMode() == PINT_MODE) {
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		cosx = ICos(ipi * ix);
		cosy = ICos(ipi * iy);
		iex = IExp(ix);
		iey = IExp(iy);

	} else {
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		cosx = ICos(ipi * ix);
		cosy = ICos(ipi * iy);
		iex = IExp(ix);
		iey = IExp(iy);
	}

	impi2y = im1 * ipi * ipi * iy;

	ip1 = itwo * iex * sinx * siny + itwo * ipi * iex * cosx * siny;
	ip2 = itwo * ipi * iey * cosx * siny + iex * ix * sinx * siny;
	ip3 = itwo * ipi * iex * ix * cosx * siny;
	ip4 = im1 * ipi2 * iex * ix * sinx * siny;
	ip5 = im1 * ipi2 * iey * ix * sinx * siny;

	r = impi2y * (ip1 + ip2 + ip3 + ip4 + ip5);

	return r;
}

template<typename T>
Interval<T> ExampleGPE11<T>::OMEGA(const Interval<T> &ix, const Interval<T> &iy,
		int &st) {
	Interval<T> r, sinx, siny, cosx, cosy, iex, iey, impi2x, ip1, ip2, ip3, ip4,
			ip5;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> itwo(2, 2);
	Interval<T> im1(-1, -1);
	Interval<T> ipi = Interval<T>::IPi();
	Interval<T> ipi2 = ipi * ipi;
	Interval<T> ipi3 = ipi * ipi * ipi;

	if (Interval<T>::GetMode() == PINT_MODE) {
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		cosx = ICos(ipi * ix);
		cosy = ICos(ipi * iy);
		iex = IExp(ix);
		iey = IExp(iy);

	} else {
		sinx = DISin(ipi * ix);
		siny = DISin(ipi * iy);
		cosx = DICos(ipi * ix);
		cosy = DICos(ipi * iy);
		iex = DIExp(ix);
		iey = DIExp(iy);
	}

	impi2x = im1 * ipi * ipi * ix;

	ip1 = itwo * iey * sinx * siny + itwo * ipi * iex * sinx * cosy;	//
	ip2 = itwo * ipi * iey * iy * sinx * cosy + iey * iy * sinx * siny;	//
	ip3 = itwo * ipi * iey * sinx * cosy;	//
	ip4 = im1 * ipi2 * iex * iy * sinx * siny;	//
	ip5 = im1 * ipi2 * iey * iy * sinx * siny;	//

	r = impi2x * (ip1 + ip2 + ip3 + ip4 + ip5);
	return r;
}

template<typename T>
Interval<T> ExampleGPE11<T>::PSI1(const Interval<T> &ix, const Interval<T> &iy,
		int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE11<T>::PSI2(const Interval<T> &ix, const Interval<T> &iy,
		int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE11<T>::PSI3(const Interval<T> &ix, const Interval<T> &iy,
		int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE11<T>::PSI4(const Interval<T> &ix, const Interval<T> &iy,
		int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE11<T>::C(const Interval<T> &ix, const Interval<T> &iy,
		int &st) {
	//c(x,y) = 0;
	Interval<T> izero = { 0, 0 };

	return izero;
}

template<typename T>
long double ExampleGPE11<T>::ExactSol(long double x, long double y) {
	long double exact = std::sin(M_PI * x) * std::sin(M_PI * y);
	return exact;
}

template<typename T>
long double ExampleGPE11<T>::GetConstM() {
	long double constM = 162;

	return constM;
}

template<typename T>
long double ExampleGPE11<T>::GetConstN() {
	long double constN = 162;

	return constN;
}

template<typename T>
long double ExampleGPE11<T>::GetConstP() {
	long double constP = 32;

	return constP;
}

template<typename T>
long double ExampleGPE11<T>::GetConstQ() {
	long double constQ = 32;

	return constQ;
}

template<typename T>
long double ExampleGPE11<T>::GetConstR() {
	long double constR = 32;

	return constR;
}

template<typename T>
long double ExampleGPE11<T>::GetConstS() {
	long double constS = 97.5;

	return constS;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::A1(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> r, iey;
	r.a = 0;
	r.b = 0;
	st = 0;

	if (Interval<T>::GetMode() == PINT_MODE) {
		iey = IExp(iy);

	} else {
		iey = DIExp(iy);
	}
	st = 0;
	r = ix * iy * iey;
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::A2(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> r, iex;
	r.a = 0;
	r.b = 0;
	st = 0;

	if (Interval<T>::GetMode() == PINT_MODE) {
		iex = IExp(ix);

	} else {
		iex = DIExp(ix);
	}
	st = 0;
	r = ix * iy * iex;
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::DCDX(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::DCDY(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::D2A1DX2(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::D2A1DY2(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> r, iey;
	Interval<T> itwo = { 2.0, 2.0 };

	r.a = 0;
	r.b = 0;
	st = 0;

	if (Interval<T>::GetMode() == PINT_MODE) {
		iey = IExp(iy);

	} else {
		iey = DIExp(iy);
	}
	st = 0;
	r = ix * (iey * iy + itwo * iey);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::D2A2DX2(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> r, iex;
	Interval<T> itwo = { 2.0, 2.0 };
	r.a = 0;
	r.b = 0;
	st = 0;

	if (Interval<T>::GetMode() == PINT_MODE) {
		iex = IExp(ix);

	} else {
		iex = DIExp(ix);
	}
	st = 0;
	r = iy * (iex * ix + itwo * iex);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::D2A2DY2(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::DFDY(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> sinx, siny, icosy, sinxsiny, isinxcosy, iexpx, iexpy;
	Interval<T> r = { 0, 0 };
	st = 0;
	Interval<T> imone(-1, -1);
	Interval<T> ipi = Interval<T>::IPi();
	Interval<T> ipi2 = ipi * ipi;

	if (Interval<T>::GetMode() == PINT_MODE) {
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		icosy = ICos(ipi * iy);
		iexpx = IExp(ix);
		iexpy = IExp(iy);
	} else {
		sinx = DISin(ipi * ix);
		siny = DISin(ipi * iy);
		icosy = DICos(ipi * iy);
		iexpx = DIExp(ix);
		iexpy = DIExp(iy);
	}
	sinxsiny = sinx * siny;
	isinxcosy = sinx * icosy;

	r =
			imone * ipi2 * ix * (iexpx * sinxsiny + iexpy * sinxsiny)
					- ipi2 * ix * iy
							* (iexpy * sinxsiny + iexpx * isinxcosy
									+ iexpy * isinxcosy);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::DFDX(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> sinx, siny, icosx, sinxsiny, icosxsiny, iexpx, iexpy;
	Interval<T> r = { 0, 0 };
	st = 0;
	Interval<T> imone(-1, -1);
	Interval<T> ipi = Interval<T>::IPi();
	Interval<T> ipi2 = ipi * ipi;

	if (Interval<T>::GetMode() == PINT_MODE) {
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		icosx = ICos(ipi * ix);
		iexpx = IExp(ix);
		iexpy = IExp(iy);
	} else {
		sinx = DISin(ipi * ix);
		siny = DISin(ipi * iy);
		icosx = DICos(ipi * ix);
		iexpx = DIExp(ix);
		iexpy = DIExp(iy);
	}
	sinxsiny = sinx * siny;
	icosxsiny = icosx * siny;

	r =
			imone * ipi2 * iy * (iexpx * sinxsiny + iexpy * sinxsiny)
					- ipi2 * ix * iy
							* (iexpx * sinxsiny + iexpx * icosxsiny
									+ iexpy * icosxsiny);
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::D2FDY2(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> sinx, siny, icosy, sinxsiny, isinxcosy, iexpx, iexpy;
	Interval<T> r = { 0, 0 };
	st = 0;
	Interval<T> imone(-1, -1);
	Interval<T> itwo(2, 2);
	Interval<T> ipi = Interval<T>::IPi();
	Interval<T> ipi2 = ipi * ipi;

	if (Interval<T>::GetMode() == PINT_MODE) {
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		icosy = ICos(ipi * iy);
		iexpx = IExp(ix);
		iexpy = IExp(iy);
	} else {
		sinx = DISin(ipi * ix);
		siny = DISin(ipi * iy);
		icosy = DICos(ipi * iy);
		iexpx = DIExp(ix);
		iexpy = DIExp(iy);
	}
	sinxsiny = sinx * siny;
	isinxcosy = sinx * icosy;

	r = imone * itwo * ipi2 * ix
			* (itwo
					* (iexpx * ipi * isinxcosy + iexpy * ipi * isinxcosy
							+ iexpy * sinxsiny)
					+ iy
							* (itwo * iexpy * ipi * isinxcosy + iexpy * sinxsiny
									- iexpx * ipi2 * sinxsiny
									- iexpy * ipi2 * sinxsiny));

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::D2FDX2(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> sinx, siny, icosx, sinxsiny, icosxsiny, iexpx, iexpy;
	Interval<T> r = { 0, 0 };
	st = 0;
	Interval<T> imone(-1, -1);
	Interval<T> itwo(2, 2);
	Interval<T> ipi = Interval<T>::IPi();
	Interval<T> ipi2 = ipi * ipi;

	if (Interval<T>::GetMode() == PINT_MODE) {
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		icosx = ICos(ipi * ix);
		iexpx = IExp(ix);
		iexpy = IExp(iy);
	} else {
		sinx = DISin(ipi * ix);
		siny = DISin(ipi * iy);
		icosx = DICos(ipi * ix);
		iexpx = DIExp(ix);
		iexpy = DIExp(iy);
	}
	sinxsiny = sinx * siny;
	icosxsiny = icosx * siny;

	r = imone * itwo * ipi2 * iy
			* (iexpx * ipi * icosxsiny + iexpy * ipi * icosxsiny
					+ iexpx * sinxsiny
					- ipi2 * ix * iy
							* (imone * (iexpy * ipi2 * sinxsiny)
									+ (itwo * iexpx * ipi * icosx + iexpx * sinx
											- iexpx * ipi2 * sinx) * siny));
	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::DA1DX(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> r = { 0, 0 };
	st = 0;
	Interval<T> iey;

	if (Interval<T>::GetMode() == PINT_MODE) {
		iey = IExp(iy);

	} else {
		iey = DIExp(iy);
	}
	st = 0;
	r = iey * iy;

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::DA1DY(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> iey;
	Interval<T> r = { 0, 0 };
	st = 0;
	Interval<T> itwo = { 2.0, 2.0 };

	r.a = 0;
	r.b = 0;
	st = 0;

	if (Interval<T>::GetMode() == PINT_MODE) {
		iey = IExp(iy);

	} else {
		iey = DIExp(iy);
	}
	st = 0;
	r = ix * (iey * iy + iey);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::DA2DX(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {

	Interval<T> r, iex;
	r.a = 0;
	r.b = 0;
	st = 0;

	if (Interval<T>::GetMode() == PINT_MODE) {
		iex = IExp(ix);

	} else {
		iex = DIExp(ix);
	}
	st = 0;
	r = iy * (iex * ix + iex);

	return r;
}

template<typename T>
inline Interval<T> ExampleGPE11<T>::DA2DY(const Interval<T> &ix,
		const Interval<T> &iy, int &st) {
	Interval<T> r, iex;
	r.a = 0;
	r.b = 0;
	st = 0;

	if (Interval<T>::GetMode() == PINT_MODE) {
		iex = IExp(ix);

	} else {
		iex = DIExp(ix);
	}
	st = 0;
	r = iex * ix;

	return r;
}

template<typename T>
void ExampleGPE11<T>::SetArithmeticMode(IAMode mode) {
	Interval<T>::SetMode(mode);
}

template<typename T>
int ExampleGPE11<T>::boundconds_classic(const long double &b1,
		const long double &b2, const long double eps) {
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
