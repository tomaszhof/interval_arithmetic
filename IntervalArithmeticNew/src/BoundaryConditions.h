/*
 * BoundaryConditions.h
 *
 *  Created on: 04-05-2013
 *      Author: thof
 */

#ifndef BOUNDARYCONDITIONS_H_
#define BOUNDARYCONDITIONS_H_
#include "Interval.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
class BoundaryConditions {
public:
	BoundaryConditions();
	virtual ~BoundaryConditions();

	//Interval<T> arithmetic functions
	virtual Interval<T> F(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> PHI1(const Interval<T>& iy, int& st);
	virtual Interval<T> PHI2(const Interval<T>& ix, int& st);
	virtual Interval<T> PHI3(const Interval<T>& iy, int& st);
	virtual Interval<T> PHI4(const Interval<T>& ix, int& st);
	virtual Interval<T> PSI(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> OMEGA(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> THETA(const Interval<T>& ix, const Interval<T>& iy,
				int& st);
	virtual Interval<T> KSI(const Interval<T>& ix, const Interval<T>& iy,
				int& st);
	virtual Interval<T> PSI1(const Interval<T>& ix, const Interval<T>& iy,
				int& st);
	virtual Interval<T> PSI2(const Interval<T>& ix, const Interval<T>& iy,
				int& st);
	virtual Interval<T> PSI3(const Interval<T>& ix, const Interval<T>& iy,
				int& st);
	virtual Interval<T> PSI4(const Interval<T>& ix, const Interval<T>& iy,
				int& st);
	virtual Interval<T> A(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> C(const Interval<T>& ix, const Interval<T>& iy,
			int& st);

	//floating-point arithmetic functions
	virtual long double f(const T& x, const T& y);
	virtual long double phi1(const T& y);
	virtual long double phi2(const T& x);
	virtual long double phi3(const T& y);
	virtual long double phi4(const T& x);
	virtual long double psi(const T& x, const T& y);
	virtual long double omega(const long double& x, const long double& y);
	virtual long double a(const long double& x, const long double& y);
	virtual long double c(const long double& x, const long double& y);
	virtual int boundconds_classic(const long double& b1, const long double& b2,
			const long double eps);
	virtual int boundconds(const Interval<T>& B1, const Interval<T>& B2,
			const long double eps);

	virtual long double ExactSol(long double x, long double y);
	virtual long double GetConstM();
	virtual long double GetConstN();
	virtual long double GetConstP();
	virtual long double GetConstQ();
	virtual void SetArithmeticMode(int mode);
};

template<typename T>
BoundaryConditions<T>::BoundaryConditions() {
	// TODO Auto-generated constructor stub

}

template<typename T>
BoundaryConditions<T>::~BoundaryConditions() {
	// TODO Auto-generated destructor stub
}

template<typename T>
long double BoundaryConditions<T>::f(const T& x,
		const T& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::phi1(const T& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::phi2(const T& x) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::phi3(const T& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::phi4(const T& x) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::psi(const T& x,
		const T& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::omega(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::a(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::c(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
Interval<T> BoundaryConditions<T>::F(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> BoundaryConditions<T>::PHI1(const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> BoundaryConditions<T>::PHI2(const Interval<T>& ix, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> BoundaryConditions<T>::PHI3(const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> BoundaryConditions<T>::PHI4(const Interval<T>& ix, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> BoundaryConditions<T>::PSI(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> BoundaryConditions<T>::OMEGA(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> BoundaryConditions<T>::THETA(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> BoundaryConditions<T>::KSI(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}


template<typename T>
Interval<T> BoundaryConditions<T>::PSI1(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> BoundaryConditions<T>::PSI2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> BoundaryConditions<T>::PSI3(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> BoundaryConditions<T>::PSI4(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}
template<typename T>
Interval<T> BoundaryConditions<T>::A(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
Interval<T> BoundaryConditions<T>::C(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
long double BoundaryConditions<T>::ExactSol(long double x, long double y) {
	long double exact = 0;
	return exact;
}

template<typename T>
long double BoundaryConditions<T>::GetConstM() {
	long double constM = 0;
	return constM;
}

template<typename T>
long double BoundaryConditions<T>::GetConstN() {
	long double constN = 0;
	return constN;
}

template<typename T>
long double BoundaryConditions<T>::GetConstP() {
	long double constP = 0;
	return constP;
}

template<typename T>
long double BoundaryConditions<T>::GetConstQ() {
	long double constQ = 0;
	return constQ;
}

template<typename T>
int BoundaryConditions<T>::boundconds_classic(const long double& b1,
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

template<typename T>
int BoundaryConditions<T>::boundconds(const Interval<T>& B1,
		const Interval<T>& B2, const long double eps) {
	int r = 0;
	if ((B1.a != 0) && (B2.a != 0)) {
		if (abs(B1.a - B2.a) / abs(B1.a) >= eps)
			r = 4;
	} else {
		if (B1.a == 0) {
			if (abs(B2.a) >= eps)
				r = 4;
		} else {
			if (B2.a == 0) {
				if (abs(B1.a) >= eps)
					r = 4;
			}
		}
	}
	if ((B1.b != 0) && (B2.b != 0)) {
		if (abs(B1.b - B2.b) / abs(B1.b) >= eps)
			r = 4;
	} else {
		if (B1.b == 0) {
			if (abs(B2.b) >= eps)
				r = 4;
		} else {
			if (B2.b == 0) {
				if (abs(B1.b) >= eps)
					r = 4;
			}
		}
	}

	return r;
}

template<typename T>
void BoundaryConditions<T>::SetArithmeticMode(int mode) {

}

} /* namespace interval_arithmetic */
#endif /* BOUNDARYCONDITIONS_H_ */
