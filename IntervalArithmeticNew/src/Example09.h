/*
 * Example09.h
 *
 *  Created on: 11-11-2014
 *      Author: thof
 */

#ifndef EXAMPLE09_H_
#define EXAMPLE09_H_

#include "Utils.h"
//#include "Interval.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class Example09 : public BoundaryConditions<T> {
public:
	Example09();
	virtual ~Example09();
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
Example09<T>::Example09()
{
	// TODO Auto-generated constructor stub

}

template<typename T>
Example09<T>::~Example09()
{
	// TODO Auto-generated destructor stub
}

template<typename T>
long double Example09<T>::f(const long double& x, const long double& y)
{
	return 0.0;
}

template<typename T>
long double Example09<T>::phi1(const long double& y)
{
	return cos(3.0*y);
}

template<typename T>
long double Example09<T>::phi2(const long double& x)
{
	return exp(3.0*x);
}

template<typename T>
long double Example09<T>::phi3(const long double& y)
{
	return exp(3.0) * cos(3.0*y);
}

template<typename T>
long double Example09<T>::phi4(const long double& x)
{
	return exp(3.0*x) * cos(3.0);
}

template<typename T>
long double Example09<T>::a(const long double& x, const long double& y)
{
	return 1.0;
}

template<typename T>
long double Example09<T>::c(const long double& x, const long double& y)
{
	return 1.0;
}

template<typename T>
Interval<T> Example09<T>::F(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//f(x,y)= 0;
	Interval<T> r;
	r.a = 0;
	r.b = 0;
	st = 0;
	return r;
}

template<typename T>
Interval<T> Example09<T>::PHI1(const Interval<T>& iy, int& st)
{
	//phi1(y) = cos(3*y);
	Interval<T> r =
	{ 0, 0 };
	Interval<T> ithree =
	{ 3, 3};

	st=0;

	r = ICos(ithree * iy);
	return r;
}

template<typename T>
Interval<T> Example09<T>::PHI2(const Interval<T>& ix, int& st)
{
	//phi2(x) = exp(3*x);
	Interval<T> r =
	{ 0, 0 };
	Interval<T> ithree =
	{ 3, 3 };

	st=0;

	r = IExp(ithree * ix);
	return r;
	return r;
}

template<typename T>
Interval<T> Example09<T>::PHI3(const Interval<T>& iy, int& st)
{
	//phi3(y) = exp(3.0) * cos(3.0*y);
	Interval<T> r =
	{ 0, 0 };
	Interval<T> ithree =
	{ 3, 3};

	st=0;

	r = IExp(ithree) * ICos(ithree * iy);
	return r;
}

template<typename T>
Interval<T> Example09<T>::PHI4(const Interval<T>& ix, int& st)
{
	//phi3(y) = exp(3.0*x) * cos(3.0);
	Interval<T> r =
	{ 0, 0 };
	Interval<T> ithree =
	{ 3, 3};
	st=0;
	r = IExp(ithree*ix) * ICos(ithree);
	return r;
}

template<typename T>
Interval<T> Example09<T>::PSI(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//psi(x,y) = 14*y*x^2 - 6*y - 3*y^2

	Interval<T> powx2, powy2, tmp;
	Interval<T> ifourteen =
	{ 14, 14 };
	Interval<T> isix =
	{ 6, 6 };
	Interval<T> ithree =
	{ 3, 3 };
	Interval<T> r =
	{ 0, 0 };

	st = 0;
	int status = 0;
	int mode = Interval<T>::GetMode();

	if (mode == PINT_MODE)
	{
		powx2 = (ix * ix);
		powy2 = (iy * iy);
		r = ifourteen *(iy * powx2);
		tmp = (isix * iy);
		r = (r - tmp);
		tmp = (ithree * powy2);
		r = (r - tmp);
	}
	else
	{
		powx2 = (ix * ix);
		powy2 = (iy * iy);
		r = (ifourteen* (iy * powx2));
		tmp = (isix * iy);
		r = (r - tmp);
		tmp = (ithree * powy2);
		r = (r - tmp);
	}
	return r;
}

template<typename T>
Interval<T> Example09<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//omega(x,y) = 14*x*y^2 - 6*x - 3*x^2
	Interval<T> powx2, powy2, tmp;
	Interval<T> ifourteen =
	{ 14, 14 };
	Interval<T> isix =
	{ 6, 6 };
	Interval<T> ithree =
	{ 3, 3 };
	Interval<T> r =
	{ 0, 0 };

	st = 0;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		powx2 = (ix * ix);
		powy2 = (iy * iy);
		r = (ifourteen *(ix * powy2));
		tmp = (isix * ix);
		r = (r - tmp);
		tmp = (ithree * powx2);
		r = (r - tmp);
	}
	else
	{
		powx2 = (ix * ix);
		powy2 = (iy * iy);
		r = (ifourteen * (ix * powy2));
		tmp = (isix * ix);
		r = (r - tmp);
		tmp = (ithree * powx2);
		r = (r - tmp);
	}
	return r;
}

template<typename T>
Interval<T> Example09<T>::PSI1(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example09<T>::PSI2(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example09<T>::PSI3(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example09<T>::PSI4(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example09<T>::A(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	// a(x,y) = 1.0;
	Interval<T> ione =
	{ 1, 1 };
	return ione;
}

template<typename T>
Interval<T> Example09<T>::C(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//c(x,y) = 1.0;
	Interval<T> ione =
	{ 1, 1 };

	return ione;
}

template<typename T>
long double Example09<T>::ExactSol(long double x, long double y)
{
	return exp(3.0*x)*cos(3.0*y);
}

template<typename T>
long double Example09<T>::GetConstM()
{
	long double constM = 1627;

	return constM;
}

template<typename T>
long double Example09<T>::GetConstN()
{
	long double constN = 1627;

	return constN;
}

template<typename T>
void Example09<T>::SetArithmeticMode(IAMode mode)
{
	Interval<T>::SetMode(mode);
}

template<typename T>
int Example09<T>::boundconds_classic(const long double& b1, const long double& b2,
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
#endif /* EXAMPLE09_H_ */
