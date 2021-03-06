/*
 * Example10.h
 *
 *  Created on: 04-05-2013
 *      Author: thof
 */

#ifndef EXAMPLE10_H_
#define EXAMPLE10_H_

#include "Utils.h"
//#include "Interval.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class Example10 : public BoundaryConditions<T> {
public:
	Example10();
	virtual ~Example10();
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
Example10<T>::Example10()
{
	// TODO Auto-generated constructor stub

}

template<typename T>
Example10<T>::~Example10()
{
	// TODO Auto-generated destructor stub
}

template<typename T>
long double Example10<T>::f(const long double& x, const long double& y)
{
	return -2.0 * M_PI*(std::sin(M_PI*x) * std::sin(M_PI*y));
}

template<typename T>
long double Example10<T>::phi1(const long double& y)
{
	return 0.0;
}

template<typename T>
long double Example10<T>::phi2(const long double& x)
{
	return 0.0;
}

template<typename T>
long double Example10<T>::phi3(const long double& y)
{
	return 0.0;
}

template<typename T>
long double Example10<T>::phi4(const long double& x)
{
	return 0.0;
}

template<typename T>
long double Example10<T>::a(const long double& x, const long double& y)
{
	return 0.0;
}

template<typename T>
long double Example10<T>::c(const long double& x, const long double& y)
{
	return 0.0;
}

template<typename T>
Interval<T> Example10<T>::F(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//f(x,y)=-2*PI*sin(PI*x)*sin(PI*y)

	Interval<T> r, sinx, siny;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> imtwo (-2,-2);
	Interval<T> ipi = Interval<T>::IPi();

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		r = imtwo * ipi *sinx * siny;
	}
	else
	{
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		r = imtwo * ipi *sinx * siny;
	}
	return r;
}

template<typename T>
Interval<T> Example10<T>::PHI1(const Interval<T>& iy, int& st)
{
	Interval<T> r =
		{ 0, 0 };
		return r;
}

template<typename T>
Interval<T> Example10<T>::PHI2(const Interval<T>& ix, int& st)
{
	Interval<T> r =
		{ 0, 0 };

		return r;
}

template<typename T>
Interval<T> Example10<T>::PHI3(const Interval<T>& iy, int& st)
{
	Interval<T> r =
		{ 0, 0 };
		return r;
}

template<typename T>
Interval<T> Example10<T>::PHI4(const Interval<T>& ix, int& st)
{
	Interval<T> r =
		{ 0, 0 };
		return r;
}

template<typename T>
Interval<T> Example10<T>::PSI(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//f(x,y)=-2*PI*sin(PI*x)*sin(PI*y)
	Interval<T> r, sinx, siny;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> itwo (2,2);
	Interval<T> ipi = Interval<T>::IPi();
	Interval<T> ipi3 = ipi * ipi * ipi;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		r = itwo * ipi3 * sinx * siny;
	}
	else
	{
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		r = itwo * ipi3 * sinx * siny;
	}
	return r;
}

template<typename T>
Interval<T> Example10<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r, sinx, siny;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> itwo (2,2);
	Interval<T> ipi = Interval<T>::IPi();
	Interval<T> ipi3 = ipi * ipi * ipi;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		r = itwo * ipi3 * sinx * siny;
	}
	else
	{
		sinx = ISin(ipi * ix);
		siny = ISin(ipi * iy);
		r = itwo * ipi3 * sinx * siny;
	}
	return r;
}

template<typename T>
Interval<T> Example10<T>::PSI1(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r, sinx, siny;
		r.a = 0;
		r.b = 0;
		st = 0;

		return r;
}

template<typename T>
Interval<T> Example10<T>::PSI2(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r, sinx, siny;
	r.a = 0;
	r.b = 0;
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example10<T>::PSI3(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r, sinx, siny;
	r.a = 0;
	r.b = 0;
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example10<T>::PSI4(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example10<T>::A(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
		{ 1, 1 };
		st = 0;

		return r;
}

template<typename T>
Interval<T> Example10<T>::C(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
		{ 1, 1 };
		st = 0;

		return r;
}

template<typename T>
long double Example10<T>::ExactSol(long double x, long double y)
{
	long double exact = std::sin(M_PI * x) * std::sin(M_PI * y);
	return exact;
}

template<typename T>
long double Example10<T>::GetConstM()
{
	long double constM = 97.5;

	return constM;
}

template<typename T>
long double Example10<T>::GetConstN()
{
	long double constN = 97.5;

	return constN;
}

template<typename T>
void Example10<T>::SetArithmeticMode(IAMode mode)
{
	Interval<T>::SetMode(mode);
}

template<typename T>
int Example10<T>::boundconds_classic(const long double& b1, const long double& b2,
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
