/*
 * Example06.h
 *
 *  Created on: 11-11-2014
 *      Author: thof
 */

#ifndef EXAMPLE06_H_
#define EXAMPLE06_H_

#include "Utils.h"
//#include "Interval.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class Example06 : public BoundaryConditions<T> {
public:
	Example06();
	virtual ~Example06();
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
Example06<T>::Example06()
{
	// TODO Auto-generated constructor stub

}

template<typename T>
Example06<T>::~Example06()
{
	// TODO Auto-generated destructor stub
}

template<typename T>
long double Example06<T>::f(const long double& x, const long double& y)
{
	return x * y * (y*std::pow(x, 2) - 3 * x + x*std::pow(y, 2) - 3 * y);
}

template<typename T>
long double Example06<T>::phi1(const long double& y)
{
	return (1-y*std::cos(y*(M_PI / 2)))*std::pow((M_PI),y);
}

template<typename T>
long double Example06<T>::phi2(const long double& x)
{
	return std::pow(M_PI,x);
}

template<typename T>
long double Example06<T>::phi3(const long double& y)
{
	return std::pow(M_PI,(1-y))*exp(y);
}

template<typename T>
long double Example06<T>::phi4(const long double& x)
{
	return std::pow(M_PI,(1-x))*exp(sin(x*M_PI/2));
}

template<typename T>
long double Example06<T>::a(const long double& x, const long double& y)
{
	return (std::sin(y*M_PI/2) + std::cos(x*M_PI/2)) * std::exp(std::pow(x, 2)/ 2);
}

template<typename T>
long double Example06<T>::c(const long double& x, const long double& y)
{
	return std::exp((std::pow(x, 2) + std::pow(y, 2)) / 2);
}

template<typename T>
Interval<T> Example06<T>::F(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//f(x,y)=x * y * (y*pow(x, 2) - 3 * x + x*pow(y, 2) - 3 * y);

	Interval<T> r, xy, xpy, xym3, i3y, xpow2, ypow2;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> ithree (3,3);
	xpow2 = ix * ix;
	ypow2 = iy * iy;

	st=0;
	r = ix * iy * (iy*xpow2 - ithree*ix + ix*ypow2 - ithree*iy);
	return r;
}

template<typename T>
Interval<T> Example06<T>::PHI1(const Interval<T>& iy, int& st)
{
	//phi1(y) = (1-y*cos(y*(M_PI / 2)))*pow((M_PI),y);

	Interval<T> powy2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> izero =
	{ 0, 0 };
	Interval<T> ione =
	{ 1, 1 };
	Interval<T> itwo =
	{ 2, 2 };

	st=0;
	//r = (ione - iy*ICos(iy*Interval<T>::IPi()/2))*IPow(Interval<T>::IPi(), iy);

	return r;
}

template<typename T>
Interval<T> Example06<T>::PHI2(const Interval<T>& ix, int& st)
{
	//phi2(x) = pow(M_PI,x);

	Interval<T> powx2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };

	st=0;
	//r = IPow(Interval<T>::IPi(),ix);

	return r;
}

template<typename T>
Interval<T> Example06<T>::PHI3(const Interval<T>& iy, int& st)
{
	//phi3(y) = pow(M_PI,(1-y))*exp(y);
	Interval<T> powy2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> ione =
	{ 1, 1 };

	st=0;
	//r = Interval<T>::IPow(Interval<T>::IPi(), (ione - iy)) * Interval<T>::IExp(iy);
	return r;
}

template<typename T>
Interval<T> Example06<T>::PHI4(const Interval<T>& ix, int& st)
{
	//phi4(x) = pow(M_PI,(1-x))*exp(sin(x*M_PI/2));

	Interval<T> powx2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> izero =
	{ 0, 0 };
	Interval<T> ione =
	{ 1, 1 };
	Interval<T> itwo =
	{ 2, 2 };

	st=0;
	//r = IPow(Interval<T>::IPi(), (ione - ix))*IExp(ISin(ix*Interval<T>::IPi()/itwo));
	return r;
}

template<typename T>
Interval<T> Example06<T>::PSI(const Interval<T>& ix, const Interval<T>& iy, int& st)
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
Interval<T> Example06<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy, int& st)
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
Interval<T> Example06<T>::PSI1(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example06<T>::PSI2(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example06<T>::PSI3(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example06<T>::PSI4(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example06<T>::A(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	// a(x,y) = (sin(y*M_PI/2) + cos(x*M_PI/2)) * exp(pow(x, 2)/ 2);

	Interval<T> powx2, powy2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> itwo =
	{ 2, 2 };

	st = 0;
	Interval<T> xpow2 = ix * ix;

	r = (ISin(iy*(Interval<T>::IPi()/itwo)) + ICos(ix*(Interval<T>::IPi()/itwo)))*IExp(xpow2/itwo);

	return r;
}

template<typename T>
Interval<T> Example06<T>::C(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//c(x,y) = exp((pow(x,2)+pow(y,2))/2);
	Interval<T> powx2, powy2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> itwo =
	{ 2, 2 };
	int status = 0;

	//b:=iadd(imul(ix,ix),imul(iy,iy));
	//b:=iexp(idiv(b,itwo),st);
	//Result:=imul(ix,b)


	if (Interval<T>::GetMode() == PINT_MODE)
	{
		powy2 = (iy * iy);
		powx2 = (ix * ix);
		tmp = (powx2 + powy2);
		tmp = (tmp / itwo);
		iexp = IExp(tmp);
		if (status == 0)
		{
			st = 0;
			//r = (ix * iexp);
			return r;
		}
		else
		{
			st = 3;
			return r;
		}
	}
	else
	{
		//b:=iadd(imul(ix,ix),imul(iy,iy));
	    //b:=iexp(idiv(b,itwo),st);
		//Result:=imul(ix,b)
		r = ((ix * ix) + (iy *iy));
		r = IExp((r / itwo));
		//r = (ix * r);

//		powy2 = ia.DIMul(iy, iy);
//		powx2 = ia.DIMul(ix, ix);
//		tmp = ia.DIAdd(powx2, powy2);
//		tmp = ia.DIDiv(tmp, itwo);
//		iexp = ia.DIExp(tmp, status);
//		if (status == 0)
//		{
//			st = 0;
//			r = ia.DIMul(ix, iexp);
//			return r;
//		}
//		else
//		{
//			st = 3;
//			return r;
//		}
	}

	return r;
}

template<typename T>
long double Example06<T>::ExactSol(long double x, long double y)
{
	//long double exact = x * y * exp((-1) * (pow(x, 2) + pow(y, 2)) / 2);
	return 0.0; //exact;
}

template<typename T>
long double Example06<T>::GetConstM()
{
	long double constM = 2.2073;

	return constM;
}

template<typename T>
long double Example06<T>::GetConstN()
{
	long double constN = 2.2073;

	return constN;
}

template<typename T>
void Example06<T>::SetArithmeticMode(IAMode mode)
{
	Interval<T>::SetMode(mode);
}

template<typename T>
int Example06<T>::boundconds_classic(const long double& b1, const long double& b2,
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
#endif /* EXAMPLE06_H_ */
