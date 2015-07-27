/*
 * Example03.h
 *
 *  Created on: 04-05-2013
 *      Author: thof
 */

#ifndef EXAMPLE03_H_
#define EXAMPLE03_H_

#include "Utils.h"
#include "Interval.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class Example03 : public BoundaryConditions<T> {
public:
	Example03();
	virtual ~Example03();
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
Example03<T>::Example03()
{
	// TODO Auto-generated constructor stub

}

template<typename T>
Example03<T>::~Example03()
{
	// TODO Auto-generated destructor stub
}

template<typename T>
long double Example03<T>::f(const long double& x, const long double& y)
{
	return x * y * (y*std::pow(x, 2) - 3 * x + x*std::pow(y, 2) - 3 * y);
}

template<typename T>
long double Example03<T>::phi1(const long double& y)
{
	return y * std::exp((-1)*(1 + std::pow(y, 2)) / 2);
}

template<typename T>
long double Example03<T>::phi2(const long double& x)
{
	return x * std::exp((-1)*(1 + std::pow(x, 2)) / 2);
}

template<typename T>
long double Example03<T>::phi3(const long double& y)
{
	return 2 * y * std::exp((-1)*(4 + std::pow(y, 2)) / 2);
}

template<typename T>
long double Example03<T>::phi4(const long double& x)
{
	return 2 * x * std::exp((-1)*(4 + std::pow(x, 2)) / 2);
}

template<typename T>
long double Example03<T>::a(const long double& x, const long double& y)
{
	return y * std::exp((std::pow(x, 2) + std::pow(y, 2)) / 2);
}

template<typename T>
long double Example03<T>::c(const long double& x, const long double& y)
{
	return x * std::exp((std::pow(x, 2) + std::pow(y, 2)) / 2);
}

template<typename T>
Interval<T> Example03<T>::F(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//f(x,y)=x*y*(x+y)*(x*y-3)

	Interval<T> r, xy, xpy, xym3, i3y;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> ithree (3,3);

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		xy = (ix * iy);
		xpy = (ix + iy);
		xym3 = (xy - ithree);
		i3y = (ithree * iy);
		r = (xy * xpy);
		r = (r * xym3);
	}
	else
	{
//		f:=isub(imul(ix,iy),ithree);
//		f:=imul(iadd(ix,iy),f);
//		st:=0;
//		Result:=imul(imul(ix,iy),f);
		r = ((ix * iy) - ithree);
		r = ((ix + iy) * r);
		r =((ix * iy) * r);


//		xy = ia.DIMul(ix, iy);
//		xpy = ia.DIAdd(ix,iy);
//		xym3 = ia.DISub(xy,ithree);
//		i3y = ia.DIMul(ithree, iy);
//		r = ia.DIMul(xy, xpy);
//		r = ia.DIMul(r, xym3);
	}
	return r;
}

template<typename T>
Interval<T> Example03<T>::PHI1(const Interval<T>& iy, int& st)
{
	//phi1(y) = y*exp(-(1+pow(y,2))/2);

	Interval<T> powy2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> izero =
	{ 0, 0 };
	Interval<T> ione =
	{ 1, 1 };
	Interval<T> itwo =
	{ 2, 2 };
	int status = 0;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		powy2 = iy * iy;
		tmp = ione + powy2;
		tmp = tmp / itwo;
		tmp = izero - tmp;
		iexp = IExp(tmp);
		if (status == 0)
		{
			st = 0;
			r = iy * iexp;
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
//		phi1:=iadd(imul(iy,iy),ione);
//		phi1:=iexp(idiv(phi1,itwo),st);
//		Result:=idiv(iy,phi1);

		r = (iy * iy) + ione;
		r = IExp((r / itwo));
		r = iy / r;

//		powy2 = ia.DIMul(iy, iy);
//		tmp = ia.DIAdd(ione, powy2);
//		tmp = ia.DIDiv(tmp, itwo);
//		tmp = ia.DISub(izero, tmp);
//		iexp = ia.DIExp(tmp, status);
//		if (status == 0)
//		{
//			st = 0;
//			r = ia.DIMul(iy, iexp);
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
Interval<T> Example03<T>::PHI2(const Interval<T>& ix, int& st)
{
	//phi2(x) = x*exp(-(1+pow(x,2))/2);

	Interval<T> powx2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> izero =
	{ 0, 0 };
	Interval<T> ione =
	{ 1, 1 };
	Interval<T> itwo =
	{ 2, 2 };
	int status = 0;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		powx2 = ix * ix;
		tmp = ione + powx2;
		tmp = tmp / itwo;
		tmp = izero - tmp;
		iexp = IExp(tmp);
		if (status == 0)
		{
			st = 0;
			r = ix * iexp;
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
//		phi2:=iadd(imul(ix,ix),ione);
//		phi2:=iexp(idiv(phi2,itwo),st);
//		Result:=idiv(ix,phi2);
		r = (ix * ix) + ione;
		r = IExp((r / itwo));
		r = ix / r;

//		powx2 = ia.DIMul(ix, ix);
//		tmp = ia.DIAdd(ione, powx2);
//		tmp = ia.DIDiv(tmp, itwo);
//		tmp = ia.DISub(izero, tmp);
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
Interval<T> Example03<T>::PHI3(const Interval<T>& iy, int& st)
{
	//phi3(y) = 2 * y * exp(-(4 + pow(y, 2)) / 2);
	Interval<T> powy2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> izero =
	{ 0, 0 };
	Interval<T> itwo =
	{ 2, 2 };
	Interval<T> ifour =
	{ 4, 4 };
	int status = 0;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		powy2 = (iy * iy);
		tmp = ifour + powy2;
		tmp = tmp / itwo;
		tmp = izero - tmp;
		iexp = IExp(tmp);
		if (status == 0)
		{
			st = 0;
			tmp = itwo * iy;
			r = tmp * iexp;
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
//		phi3:=iadd(dimul(iy,iy),ifour);
//		phi3:=iexp(idiv(phi3,itwo),st);
//		Result:=imul(itwo,idiv(iy,phi3));
		r = (iy * iy) + ifour;
		r = IExp((r / itwo));
		r = itwo * (iy / r);

//		powy2 = ia.DIMul(iy, iy);
//		tmp = ia.DIAdd(ifour, powy2);
//		tmp = ia.DIDiv(tmp, itwo);
//		tmp = ia.DISub(izero, tmp);
//		iexp = ia.DIExp(tmp, status);
//		if (status == 0)
//		{
//			st = 0;
//			tmp = ia.DIMul(itwo, iy);
//			r = ia.DIMul(tmp, iexp);
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
Interval<T> Example03<T>::PHI4(const Interval<T>& ix, int& st)
{
	//phi4(x) = 2 * x * exp(-(4 + pow(x, 2)) / 2);

	Interval<T> powx2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> izero =
	{ 0, 0 };
	Interval<T> itwo =
	{ 2, 2 };
	Interval<T> ifour =
	{ 4, 4 };
	int status = 0;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		powx2 = (ix * ix);
		tmp = (ifour + powx2);
		tmp = (tmp / itwo);
		tmp = (izero - tmp);
		iexp = IExp(tmp);
		if (status == 0)
		{
			st = 0;
			tmp = (itwo * ix);
			r = (tmp * iexp);
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

//		phi4:=iadd(imul(ix,ix),ifour);
//		phi4:=iexp(idiv(phi4,itwo),st);
//		Result:=imul(itwo,idiv(ix,phi4));
		r = ((ix * ix) + ifour);
		r = IExp((r / itwo));
		r = itwo * (ix / r);

//		powx2 = ia.DIMul(ix, ix);
//		tmp = ia.DIAdd(ifour, powx2);
//		tmp = ia.DIDiv(tmp, itwo);
//		tmp = ia.DISub(izero, tmp);
//		iexp = ia.DIExp(tmp, status);
//		if (status == 0)
//		{
//			st = 0;
//			tmp = ia.DIMul(itwo, ix);
//			r = ia.DIMul(tmp, iexp);
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
Interval<T> Example03<T>::PSI(const Interval<T>& ix, const Interval<T>& iy, int& st)
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
Interval<T> Example03<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy, int& st)
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
Interval<T> Example03<T>::A(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	// a(x,y) = y * exp((pow(x, 2) + pow(y, 2)) / 2);
//
//	a:=iadd(imul(ix,ix),imul(iy,iy));
//	a:=iexp(idiv(a,itwo),st);
//	Result:=imul(iy,a)

	Interval<T> powx2, powy2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> itwo =
	{ 2, 2 };
	int status = 0;
	int mode = Interval<T>::GetMode();

	if (mode == PINT_MODE)
	{
		powy2 = (iy * iy);
		powx2 = (ix * ix);
		tmp = (powx2 + powy2);
		tmp = (tmp / itwo);
		iexp = IExp(tmp);
		if (status == 0)
		{
			st = 0;
			r = (iy * iexp);
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
		//	a:=iadd(imul(ix,ix),imul(iy,iy));
		//	a:=iexp(idiv(a,itwo),st);
		//	Result:=imul(iy,a)
		r = ((ix * ix)+ (iy* iy));
		r = IExp((r / itwo));
		r = (iy * r);

//		powy2 = ia.DIMul(iy, iy);
//		powx2 = ia.DIMul(ix, ix);
//		tmp = ia.DIAdd(powx2, powy2);
//		tmp = ia.DIDiv(tmp, itwo);
//		iexp = ia.DIExp(tmp, status);
//		if (status == 0)
//		{
//			st = 0;
//			r = ia.DIMul(iy, iexp);
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
Interval<T> Example03<T>::C(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//c(x,y) = x*exp((pow(x,2)+pow(y,2))/2);
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
			r = (ix * iexp);
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
		r = (ix * r);

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
long double Example03<T>::ExactSol(long double x, long double y)
{
	long double exact = x * y * std::exp((-1) * (std::pow(x, 2) + std::pow(y, 2)) / 2);
	return exact;
}

template<typename T>
long double Example03<T>::GetConstM()
{
	long double constM = 2.2073;

	return constM;
}

template<typename T>
long double Example03<T>::GetConstN()
{
	long double constN = 2.2073;

	return constN;
}

template<typename T>
void Example03<T>::SetArithmeticMode(IAMode mode)
{
	Interval<T>::SetMode(mode);
}

template<typename T>
int Example03<T>::boundconds_classic(const long double& b1, const long double& b2,
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
#endif /* EXAMPle03_H_ */
