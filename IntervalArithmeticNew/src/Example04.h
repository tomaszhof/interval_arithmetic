/*
 * Example04.h
 *
 *  Created on: 04-05-2013
 *      Author: thof
 */

#ifndef EXAMPLE04_H_
#define EXAMPLE04_H_

#include "Utils.h"
//#include "Interval.h"
#include "BoundaryConditions.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {



template <typename T>
class Example04 : public BoundaryConditions<T> {
public:
	Example04();
	virtual ~Example04();
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

template <typename T>
Example04<T>::Example04()
{
	// TODO Auto-generated constructor stub

}

template <typename T>
Example04<T>::~Example04()
{
	// TODO Auto-generated destructor stub
}

template <typename T>
long double Example04<T>::f(const long double& x, const long double& y)
{
	return std::pow(x, 2) * std::pow(y, 2)
			* (3 * std::pow(y, 2) + 2 * std::pow(x, 2) * std::pow(y, 2) - 3 * std::pow(x, 2));
}

template <typename T>
long double Example04<T>::phi1(const long double& y)
{
	return y * std::exp((1 - std::pow(y, 2)) / 2);
}

template <typename T>
long double Example04<T>::phi2(const long double& x)
{
	return x * std::exp((std::pow(x, 2) - 1) / 2);
}

template <typename T>
long double Example04<T>::phi3(const long double& y)
{
	return 2 * y * std::exp((4 - std::pow(y, 2)) / 2);
}

template <typename T>
long double Example04<T>::phi4(const long double& x)
{
	return 2 * x * std::exp((std::pow(x, 2) - 4) / 2);
}

template <typename T>
long double Example04<T>::a(const long double& x, const long double& y)
{
	return x * std::pow(y, 3) * std::exp((-1) * (std::pow(x, 2) - std::pow(y, 2)) / 2);
}

template <typename T>
long double Example04<T>::c(const long double& x, const long double& y)
{
	return y * std::pow(x, 3) * std::exp((-1) * (std::pow(x, 2) - std::pow(y, 2)) / 2);
}

template <typename T>
Interval<T> Example04<T>::F(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//f(x,y)= pow(x, 2) * pow(y, 2) * (3 * pow(y, 2) + 2 * pow(x, 2) * pow(y, 2) - 3 * pow(x, 2))

	Interval<T> r, powx2, powy2, powx2powy2;
	r.a = 0;
	r.b = 0;
	st = 0;

	Interval<T> itwo =
	{ 2, 2 };
	Interval<T> ithree =
	{ 3, 3 };

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		powx2 = ix * ix;
		powy2 = iy * iy;
		powx2powy2 = powx2 * powy2;

		r = ithree * powy2;
		r = r + (itwo * powx2powy2);
		r = r - (ithree * powx2);
		r = r * powx2powy2;
	}
	else
	{
		powx2 =(ix * ix);
		powy2 = (iy * iy);
		powx2powy2 = (powx2 * powy2);

		r = (ithree * powy2);
		r = r + (itwo * powx2powy2);
		r = r - (ithree * powx2);
		r = powx2powy2 * r;
	}
	return r;
}

template <typename T>
Interval<T> Example04<T>::PHI1(const Interval<T>& iy, int& st)
{
	//phi1(y) = y * exp((1 - pow(y, 2)) / 2);

	Interval<T> powy2, iexp;
	Interval<T> ione =
	{ 1, 1 };
	Interval<T> itwo =
	{ 2, 2 };
	Interval<T> r =
	{ 0, 0 };

	int status = 0;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		powy2 = (iy * iy);
		iexp = IExp((ione - powy2) / itwo);
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
		powy2 = (iy * iy);
		iexp = DIExp((ione - powy2) / itwo);
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

	return r;
}

template <typename T>
Interval<T> Example04<T>::PHI2(const Interval<T>& ix, int& st)
{
	//phi2(x) = x * exp((pow(x, 2) - 1) / 2);

	Interval<T> powx2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> ione =
	{ 1, 1 };
	Interval<T> itwo =
	{ 2, 2 };

	int status = 0;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		powx2 = ix * ix;
		iexp = IExp(((powx2 - ione) / itwo));
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
		powx2 = (ix * ix);
		iexp = DIExp((powx2 - ione) / itwo);
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

	return r;
}

template <typename T>
Interval<T> Example04<T>::PHI3(const Interval<T>& iy, int& st)
{
	//phi3(y) = 2 * y * exp((4 - pow(y, 2)) /2 );
	Interval<T> powy2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> itwo =
	{ 2, 2 };
	Interval<T> ifour =
	{ 4, 4 };
	int status = 0;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		powy2 = (iy * iy);
		iexp = IExp((ifour - powy2) / itwo);
		if (status == 0)
		{
			st = 0;
			r = ((itwo * iy) * iexp);
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
		powy2 = (iy * iy);
		iexp = DIExp((ifour - powy2) / itwo);
		if (status == 0)
		{
			st = 0;
			r = ((itwo * iy) * iexp);
			return r;
		}
		else
		{
			st = 3;
			return r;
		}
	}

	return r;
}

template <typename T>
Interval<T> Example04<T>::PHI4(const Interval<T>& ix, int& st)
{
	//phi4(x) = 2 * x * exp((pow(x, 2) - 4) / 2);

	Interval<T> powx2, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> itwo =
	{ 2, 2 };
	Interval<T> ifour =
	{ 4, 4 };
	int status = 0;

	if (Interval<T>::GetMode() == PINT_MODE)
	{
		powx2 = (ix * ix);
		iexp = IExp((powx2 - ifour) / itwo);
		if (status == 0)
		{
			st = 0;
			r = ((itwo * ix) * iexp);
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
		powx2 = (ix * ix);
		iexp = DIExp((powx2 - ifour) / itwo);
		if (status == 0)
		{
			st = 0;
			r = ((itwo * ix) * iexp);
			return r;
		}
		else
		{
			st = 3;
			return r;
		}
	}

	return r;
}

template <typename T>
Interval<T> Example04<T>::PSI(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };

	st = 0;

	return r;
}

template <typename T>
Interval<T> Example04<T>::OMEGA(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };

	st = 0;

	return r;
}

template<typename T>
Interval<T> Example04<T>::PSI1(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example04<T>::PSI2(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example04<T>::PSI3(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> Example04<T>::PSI4(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	Interval<T> r =
	{ 0, 0 };
	st = 0;

	return r;
}

template <typename T>
Interval<T> Example04<T>::A(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	// a(x,y) = x * pow(y, 3) * exp((-1) * (pow(x, 2) - pow(y, 2)) / 2);

	Interval<T> powx2, powy2, powy3, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };

	Interval<T> iminusone =
	{ -1, -1 };
	Interval<T> itwo =
	{ 2, 2 };
	int status = 0;
	int mode = Interval<T>::GetMode();

	if (mode == PINT_MODE)
	{
		powy2 = (iy * iy);
		powx2 = (ix * ix);
		powy3 = (iy * powy2);
		tmp = (powx2 - powy2);
		tmp = (iminusone * (tmp / itwo));
		iexp = IExp(tmp);
		if (status == 0)
		{
			st = 0;
			r = (powy3 * iexp);
			r = (ix * r);
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
		powy2 = (iy * iy);
		powx2 = (ix * ix);
		powy3 = (iy * powy2);
		tmp = (powx2 - powy2);
		tmp = (iminusone * (tmp / itwo));
		iexp = DIExp(tmp);
		if (status == 0)
		{
			st = 0;
			r = (powy3 * iexp);
			r = (ix * r);
			return r;
		}
		else
		{
			st = 3;
			return r;
		}
	}

	return r;
}

template <typename T>
Interval<T> Example04<T>::C(const Interval<T>& ix, const Interval<T>& iy, int& st)
{
	//c(x,y) = y * pow(x, 3) * exp((-1) * (pow(x, 2) - pow(y, 2)) / 2);
	Interval<T> powx2, powy2, powx3, iexp, tmp;
	Interval<T> r =
	{ 0, 0 };
	Interval<T> iminusone =
	{ -1, -1 };
	Interval<T> itwo =
	{ 2, 2 };
	int status = 0;
	int mode = Interval<T>::GetMode();

	if (mode == PINT_MODE)
	{
		powy2 = (iy * iy);
		powx2 = (ix * ix);
		powx3 = (ix * powx2);
		tmp = (powx2 - powy2);
		tmp = (iminusone * (tmp / itwo));
		iexp = IExp(tmp);
		if (status == 0)
		{
			st = 0;
			r = (powx3 * iexp);
			r = (iy * r);
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
		powy2 = (iy * iy);
		powx2 = (ix * ix);
		powx3 = (ix * powx2);
		tmp = (powx2 - powy2);
		tmp = (iminusone * (tmp / itwo));
		iexp = DIExp(tmp);
		if (status == 0)
		{
			st = 0;
			r = (powx3 * iexp);
			r = (iy * r);
			return r;
		}
		else
		{
			st = 3;
			return r;
		}
	}

	return r;
}

template <typename T>
long double Example04<T>::ExactSol(long double x, long double y)
{
	long double exact = x * y * std::exp((std::pow(x, 2) - std::pow(y, 2)) / 2);
	return exact;
}

template <typename T>
long double Example04<T>::GetConstM()
{
	long double constM = 636.4;

	return constM;
}

template <typename T>
long double Example04<T>::GetConstN()
{
	long double constM = 53.79;

	return constM;
}

template <typename T>
void Example04<T>::SetArithmeticMode(IAMode mode)
{
	Interval<T>::SetMode(mode);
}

template <typename T>
int Example04<T>::boundconds_classic(const long double& b1, const long double& b2,
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
#endif /* EXAMPLE04_H_ */
