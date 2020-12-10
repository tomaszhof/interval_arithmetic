/*
 * Example11.h
 *
 *  Created on: 04-05-2013
 *      Author: thof
 */

#ifndef EXAMPLE11_H_
#define EXAMPLE11_H_

#include "Utils.h"
//#include "Interval.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class Example11 : public BoundaryConditions<T> {
public:
	Example11();
	virtual ~Example11();
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
	Interval<T> A(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> C(const Interval<T>& ix, const Interval<T>& iy, int& st);

	//classical functions
	long double f(const long double& x, const long double& y);
	long double phi1(const long double& y);
	long double phi2(const long double& x);
	long double phi3(const long double& y);
	long double phi4(const long double& x);
	long double a(const long double& x, const long double& y);
	long double c(const long double& x, const long double& y);
	int boundconds_classic(const long double& b1, const long double& b2,
	    		const long double eps);

	long double ExactSol(long double x, long double y);
	long double GetConstM();
	long double GetConstN();
	void SetArithmeticMode(IAMode mode);
};

template<typename T>
Example11<T>::Example11()
{
	// TODO Auto-generated constructor stub

}

template<typename T>
Example11<T>::~Example11()
{
	// TODO Auto-generated destructor stub
}

template<typename T>
long double Example11<T>::f(const long double& x, const long double& y)
{
	long double sinXsinY = std::sin(M_PI*x) * std::sin(M_PI*y);
	return -1.0*M_PI*M_PI*x*y*(std::exp(x)*sinXsinY + std::exp(y)*sinXsinY);
}

template<typename T>
long double Example11<T>::phi1(const long double& y)
{
	return 0.0;
}

template<typename T>
long double Example11<T>::phi2(const long double& x)
{
	return 0.0;
}

template<typename T>
long double Example11<T>::phi3(const long double& y)
{
	return 0.0;
}

template<typename T>
long double Example11<T>::phi4(const long double& x)
{
	return 0.0;
}

template<typename T>
long double Example11<T>::a(const long double& x, const long double& y)
{
	return x*y*std::exp(y);
}

template<typename T>
long double Example11<T>::c(const long double& x, const long double& y)
{
	return x*y*std::exp(x);
}

template<typename T>
Interval<T> Example11<T>::F(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//f(x,y)=-1.0*M_PI*M_PI*x*y*(std::exp(x)*sinXsinY + std::exp(y)*sinXsinY);

	Interval<T> r, sinx, siny, iexpx, iexpy;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> imone (-1,-1);
	Interval<T> ipi = Interval<T>::IPi();

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		iexpx = IExp(ix);
		iexpy = IExp(iy);
	}
	else
	{
		sinx = DISin(ipi * ix);
		siny = DISin(ipi * iy);
		iexpx = DIExp(ix);
		iexpy = DIExp(iy);
	}

	r = imone*ipi*ipi*ix*iy *(iexpx*sinx*siny + iexpy*sinx*siny);
	return r;
}

template<typename T>
Interval<T> Example11<T>::PHI1(const Interval<T>& iy, int& st)
{
	Interval<T> r =
		{ 0, 0 };
		return r;
}

template<typename T>
Interval<T> Example11<T>::PHI2(const Interval<T>& ix, int& st)
{
	Interval<T> r =
		{ 0, 0 };

		return r;
}

template<typename T>
Interval<T> Example11<T>::PHI3(const Interval<T>& iy, int& st)
{
	Interval<T> r =
		{ 0, 0 };
		return r;
}

template<typename T>
Interval<T> Example11<T>::PHI4(const Interval<T>& ix, int& st)
{
	Interval<T> r =
		{ 0, 0 };
		return r;
}

template<typename T>
Interval<T> Example11<T>::PSI(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
//	d2f/dx2 = -\pi^2*y*(2e^x*sin(pi*x)*sin(pi*y)+2*pi*e^x*cos(pi*x)*sin(pi*y)
//						 +2*pi*e^y*cos(pi*x)*sin(pi*y)+e^x*x*sin(pi*x)*sin(pi*y)
//						 +2*pi*e^x*x*cos(pi*x)*sin(pi*y)
//						 -pi^2*e^x*x*sin(pi*x)*sin(pi*y)
//						 -pi^2*e^y*x*sin(pi x)*sin(pi*y))

	Interval<T> r, sinx, siny,cosx, cosy, iex, iey, impi2y, ip1, ip2,ip3, ip4, ip5;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> itwo (2,2);
	Interval<T> im1 (-1,-1);
	Interval<T> ipi = Interval<T>::IPi();
	Interval<T> ipi2 = ipi * ipi;
	Interval<T> ipi3 = ipi * ipi * ipi;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		cosx = ICos(ipi * ix);
		cosy = ICos(ipi * iy);
		iex = IExp(ix);
		iey=IExp(iy);

	}
	else
	{
		sinx = ISin(ipi * ix);
				siny = ISin(ipi * iy);
				cosx = ICos(ipi * ix);
				cosy = ICos(ipi * iy);
				iex = IExp(ix);
				iey=IExp(iy);
	}

	impi2y=im1*ipi*ipi*iy;

	ip1 = itwo*iex*sinx*siny + itwo*ipi*iex*cosx*siny;
	ip2 = itwo*ipi*iey*cosx*siny + iex*ix*sinx*siny;
	ip3 = itwo*ipi*iex*ix*cosx*siny;
	ip4 = im1*ipi2*iex*ix*sinx*siny;
	ip5 = im1*ipi2*iey*ix*sinx*siny;

	r = impi2y*(ip1+ip2+ip3+ip4+ip5);

	return r;
}

template<typename T>
Interval<T> Example11<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r, sinx, siny,cosx, cosy, iex, iey, impi2x, ip1, ip2,ip3, ip4, ip5;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> itwo (2,2);
	Interval<T> im1 (-1,-1);
	Interval<T> ipi = Interval<T>::IPi();
	Interval<T> ipi2 = ipi * ipi;
	Interval<T> ipi3 = ipi * ipi * ipi;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		cosx = ICos(ipi * ix);
		cosy = ICos(ipi * iy);
		iex = IExp(ix);
		iey=IExp(iy);

	}
	else
	{
		sinx = DISin(ipi * ix);
		siny = DISin(ipi * iy);
		cosx = DICos(ipi * ix);
		cosy = DICos(ipi * iy);
		iex = DIExp(ix);
		iey= DIExp(iy);
	}

	impi2x=im1*ipi*ipi*ix;

	ip1 = itwo*iey*sinx*siny + itwo*ipi*iex*sinx*cosy;//
	ip2 = itwo*ipi*iey*iy*sinx*cosy + iey*iy*sinx*siny;//
	ip3 = itwo*ipi*iey*sinx*cosy;//
	ip4 = im1*ipi2*iex*iy*sinx*siny;//
	ip5 = im1*ipi2*iey*iy*sinx*siny;//

	r = impi2x*(ip1+ip2+ip3+ip4+ip5);
	return r;
}

template<typename T>
Interval<T> Example11<T>::PSI1(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r, sinx, siny;
	r.a = 0;
	r.b = 0;
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example11<T>::PSI2(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r, sinx, siny;
	r.a = 0;
	r.b = 0;
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example11<T>::PSI3(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r, sinx, siny;
	r.a = 0;
	r.b = 0;
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example11<T>::PSI4(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example11<T>::A(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r, iey;
	r.a = 0;
	r.b = 0;
	st = 0;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		iey=IExp(iy);

	}
	else
	{
		iey= DIExp(iy);
	}
	st = 0;
	r = ix*iy*iey;
	return r;
}

template<typename T>
Interval<T> Example11<T>::C(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r, iex;
	r.a = 0;
	r.b = 0;
	st = 0;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		iex=IExp(ix);

	}
	else
	{
		iex= DIExp(ix);
	}
	st = 0;
	r = ix*iy*iex;
	return r;
}

template<typename T>
long double Example11<T>::ExactSol(long double x, long double y)
{
	long double exact = std::sin(M_PI * x) * std::sin(M_PI * y);
	return exact;
}

template<typename T>
long double Example11<T>::GetConstM()
{
	long double constM = 97.5;

	return constM;
}

template<typename T>
long double Example11<T>::GetConstN()
{
	long double constN = 97.5;

	return constN;
}

template<typename T>
void Example11<T>::SetArithmeticMode(IAMode mode)
{
	Interval<T>::SetMode(mode);
}

template<typename T>
int Example11<T>::boundconds_classic(const long double& b1, const long double& b2,
		const long double eps)
{
	if ((b1 != 0) && (b2 != 0))
	{
		if (abs(b1 - b2) / abs(b1) >= eps)
			return 4;
		else
			return 0;
	}
	else if (b1 == 0)
	{
		if (abs(b2) >= eps)
			return 4;
		else
			return 0;
	}
	else if (b2 == 0)
	{
		if (abs(b1) >= eps)
			return 4;
		else
			return 0;
	}
	else
		return 0;

}

} /* namespace interval_arithmetic */
#endif /* EXAMPLE10_H_ */
