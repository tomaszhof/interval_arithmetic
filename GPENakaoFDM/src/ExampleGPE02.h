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
class ExampleGPE02 : public BoundaryConditions<T> {
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
	Interval<T> A(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> C(const Interval<T>& ix, const Interval<T>& iy, int& st);

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
	void SetArithmeticMode(IAMode mode);
};

template<typename T>
ExampleGPE02<T>::ExampleGPE02()
{
	// TODO Auto-generated constructor stub

}

template<typename T>
ExampleGPE02<T>::~ExampleGPE02()
{
	// TODO Auto-generated destructor stub
}

template<typename T>
long double ExampleGPE02<T>::f(const long double& x, const long double& y)
{
	return (2.0 * PI - 1.0) *sin(PI * x)*sin(PI * y);
}

template<typename T>
long double ExampleGPE02<T>::phi1(const long double& y)
{
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::phi2(const long double& x)
{
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::phi3(const long double& y)
{
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::phi4(const long double& x)
{
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::c(const long double& x, const long double& y)
{
	return -PI;
}

template<typename T>
long double ExampleGPE02<T>::dcdx(const long double& x, const long double& y)
{
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::dcdy(const long double& x, const long double& y)
{
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::a1(const long double& x, const long double& y)
{
	return -1.0;
}

template<typename T>
long double ExampleGPE02<T>::d2a1dx2(const long double& x, const long double& y)
{
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::d2a1dy2(const long double& x, const long double& y)
{
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::a2(const long double& x, const long double& y)
{
	return -1.0;
}

template<typename T>
long double ExampleGPE02<T>::d2a2dx2(const long double& x, const long double& y)
{
	return 0.0;
}

template<typename T>
long double ExampleGPE02<T>::d2a2dy2(const long double& x, const long double& y)
{
	return 0.0;
}

template<typename T>
Interval<T> ExampleGPE02<T>::F(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//f(x,y)= 0;
	Interval<T> r;
	r.a = 0;
	r.b = 0;
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PHI1(const Interval<T>& iy, int& st)
{
	//phi1(y) = 0.0;
	Interval<T> r =
	{ 0, 0 };

	st=0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PHI2(const Interval<T>& ix, int& st)
{
	//phi2(x) = 0.0;
	Interval<T> r =
	{ 0, 0 };

	st=0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PHI3(const Interval<T>& iy, int& st)
{
	//phi3(y) = 0.0;
	Interval<T> r =
		{ 0, 0 };

	st=0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PHI4(const Interval<T>& ix, int& st)
{
	//phi3(y) = 0.0;
	Interval<T> r =
	{ 0, 0 };
	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PSI(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//psi(x,y) = 0
	Interval<T> r =
	{ 0, 0 };

	st = 0;
	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//omega(x,y) = 0
	Interval<T> r =
	{ 0, 0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PSI1(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PSI2(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PSI3(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::PSI4(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> ExampleGPE02<T>::A(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	// a(x,y) = 1.0;
	Interval<T> ione =
	{ 1, 1 };
	return ione;
}

template<typename T>
Interval<T> ExampleGPE02<T>::C(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//c(x,y) = 1.0;
	Interval<T> ione =
	{ 1, 1 };

	return ione;
}

template<typename T>
long double ExampleGPE02<T>::ExactSol(long double x, long double y)
{
	return (1.0 / PI) * sin(PI*x)*sin(PI*y);
}

template<typename T>
long double ExampleGPE02<T>::GetConstM()
{
	long double constM = 1627;

	return constM;
}

template<typename T>
long double ExampleGPE02<T>::GetConstN()
{
	long double constN = 1627;

	return constN;
}

template<typename T>
long double ExampleGPE02<T>::GetConstP()
{
	long double constP = 14643;

	return constP;
}

template<typename T>
long double ExampleGPE02<T>::GetConstQ()
{
	long double constQ = 14643;

	return constQ;
}

template<typename T>
void ExampleGPE02<T>::SetArithmeticMode(IAMode mode)
{
	Interval<T>::SetMode(mode);
}

template<typename T>
int ExampleGPE02<T>::boundconds_classic(const long double& b1, const long double& b2,
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
#endif /* EXAMPLEGPE01_H_ */
