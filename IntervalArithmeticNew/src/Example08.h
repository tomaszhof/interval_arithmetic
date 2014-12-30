/*
 * Example08.h
 *
 *  Created on: 11-11-2014
 *      Author: thof
 */

#ifndef EXAMPLE08_H_
#define EXAMPLE08_H_

#include "Utils.h"
#include "Interval.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class Example08 : public BoundaryConditions<T> {
public:
	Example08();
	virtual ~Example08();
	Interval<T> F(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> PHI1(const Interval<T>& iy, int& st);
	Interval<T> PHI2(const Interval<T>& ix, int& st);
	Interval<T> PHI3(const Interval<T>& iy, int& st);
	Interval<T> PHI4(const Interval<T>& ix, int& st);
	Interval<T> PSI(const Interval<T>& ix, const Interval<T>& iy, int& st);
	Interval<T> OMEGA(const Interval<T>& ix, const Interval<T>& iy, int& st);
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
Example08<T>::Example08()
{
	// TODO Auto-generated constructor stub

}

template<typename T>
Example08<T>::~Example08()
{
	// TODO Auto-generated destructor stub
}

template<typename T>
long double Example08<T>::f(const long double& x, const long double& y)
{
	return (1/(x*y))*cos(y*M_PI/2)*cos(x*M_PI/2);//x*y*exp(cos(x*M_PI/2))*exp(cos(y*M_PI/2));//pow(x,2)*pow(y,3)-pow(x,3)*pow(y,2);//x*pow(y,2)/(pow(x,2) + pow(y,4));
}

template<typename T>
long double Example08<T>::phi1(const long double& y)
{
	//long double pi = 4*atan(1.0);
	return sin((y-1)*M_PI); //y-1; //exp(y-1);
}

template<typename T>
long double Example08<T>::phi2(const long double& x)
{
	//long double pi = 4*atan(1.0);
	//return exp(sin((x-1)*(pi/2.0)));
	return sin((x-1)*M_PI); //x-1; //exp(x-1);
}

template<typename T>
long double Example08<T>::phi3(const long double& y)
{
	return sin((2-y)*M_PI); //2-y;//exp(2-y);
}

template<typename T>
long double Example08<T>::phi4(const long double& x)
{
	return sin((2-x)*M_PI); //2-x; //exp(2-x);
}

template<typename T>
long double Example08<T>::a(const long double& x, const long double& y)
{
	return exp(cos(x*M_PI)-cos(y*M_PI));//x*(y-1)*exp(x*x);
}

template<typename T>
long double Example08<T>::c(const long double& x, const long double& y)
{
	return exp(x - y);//(x-1)*y*exp(y*y);//sin((2*x-1)*M_PI)*cos((2*y-3)*M_PI);
}

template<typename T>
Interval<T> Example08<T>::F(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//f(x,y)=(1/(x*y))*cos(y*M_PI/2)*cos(x*M_PI/2);
	Interval<T> r, xy, xpy, xym3, i3y, xpow2, ypow2;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> itwo (2,2);
	Interval<T> ione (1,1);
	Interval<T> ionexy = ione/(ix * iy);
	Interval<T> xcos = ICos(ix*Interval<T>::IPi()/itwo);
	Interval<T> ycos = ICos(iy*Interval<T>::IPi()/itwo);

	r = ionexy*ycos*xcos;
	return r;
}

template<typename T>
Interval<T> Example08<T>::PHI1(const Interval<T>& iy, int& st)
{
	//phi1(y) = sin((y-1)*M_PI)
	Interval<T> r =
	{ 0, 0 };
	Interval<T> ione =
	{ 1, 1 };
	r = ISin((iy-ione)*Interval<T>::IPi());

	return r;
}

template<typename T>
Interval<T> Example08<T>::PHI2(const Interval<T>& ix, int& st)
{
	//phi2(x) = sin((x-1)*M_PI);
	Interval<T> ione =
		{ 1, 1 };
	Interval<T> r =
	{ 0, 0 };

	r = ISin((ix-ione)*Interval<T>::IPi());
	return r;
}

template<typename T>
Interval<T> Example08<T>::PHI3(const Interval<T>& iy, int& st)
{
	//phi3(y) = sin((2-y)*M_PI);
	Interval<T> r =
	{ 0, 0 };
	Interval<T> itwo =
		{ 2, 2 };

	r = ISin((itwo - iy)*Interval<T>::IPi());
	return r;
}

template<typename T>
Interval<T> Example08<T>::PHI4(const Interval<T>& ix, int& st)
{
	//phi4(y) = sin((2-x)*M_PI);
	Interval<T> r =
	{ 0, 0 };
	Interval<T> itwo =
	{ 2, 2 };

	r = ISin((itwo - ix)*Interval<T>::IPi());
	return r;
}

template<typename T>
Interval<T> Example08<T>::PSI(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//psi(x,y) = 0
	Interval<T> r =
	{ 0, 0 };

	return r;
}

template<typename T>
Interval<T> Example08<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//omega(x,y) = 0

	Interval<T> r =
	{ 0, 0 };
	return r;
}

template<typename T>
Interval<T> Example08<T>::A(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	// a(x,y) =  exp(cos(x*M_PI)-cos(y*M_PI));
	Interval<T> r =
	{ 0, 0 };
	Interval<T> xcos = ICos(ix*Interval<T>::IPi());
	Interval<T> ycos = ICos(iy*Interval<T>::IPi());

	r = IExp(xcos-ycos);

	return r;
}

template<typename T>
Interval<T> Example08<T>::C(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//c(x,y) = exp(x - y);
	Interval<T> r =
		{ 0, 0 };
	r = IExp(ix-iy);

	return r;
}

template<typename T>
long double Example08<T>::ExactSol(long double x, long double y)
{
	return 0.0;
}

template<typename T>
long double Example08<T>::GetConstM()
{
	long double constM = 640;

	return constM;
}

template<typename T>
long double Example08<T>::GetConstN()
{
	long double constN = 650;

	return constN;
}

template<typename T>
void Example08<T>::SetArithmeticMode(IAMode mode)
{
	Interval<T>::SetMode(mode);
}

template<typename T>
int Example08<T>::boundconds_classic(const long double& b1, const long double& b2,
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
#endif /* EXAMPLE08_H_ */
