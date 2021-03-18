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
	virtual Interval<T> A1(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> A2(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> DA1DX(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> DA1DY(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> DA2DX(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> DA2DY(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> D2A1DX2(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> D2A1DY2(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> D2A2DX2(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> D2A2DY2(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> C(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> DCDX(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> DCDY(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> D2CDX2(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> D2CDY2(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> DFDX(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> DFDY(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> D2FDX2(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	virtual Interval<T> D2FDY2(const Interval<T>& ix, const Interval<T>& iy,
			int& st);
	//floating-point arithmetic functions
	virtual long double f(const T& x, const T& y);
	virtual long double phi1(const T& y);
	virtual long double phi2(const T& x);
	virtual long double phi3(const T& y);
	virtual long double phi4(const T& x);
	virtual long double psi(const T& x, const T& y);
	virtual long double omega(const long double& x, const long double& y);
	virtual long double a1(const long double& x, const long double& y);
	virtual long double a2(const long double& x, const long double& y);
	virtual long double da1dx(const long double& x, const long double& y);
	virtual long double da2dx(const long double& x, const long double& y);
	virtual long double da1dy(const long double& x, const long double& y);
	virtual long double da2dy(const long double& x, const long double& y);
	virtual long double d2a1dx2(const long double& x, const long double& y);
	virtual long double d2a1dy2(const long double& x, const long double& y);
	virtual long double d2a2dx2(const long double& x, const long double& y);
	virtual long double d2a2dy2(const long double& x, const long double& y);
	virtual long double c(const long double& x, const long double& y);
	virtual long double dcdx(const long double& x, const long double& y);
	virtual long double dcdy(const long double& x, const long double& y);
	virtual long double d2cdx2(const long double& x, const long double& y);
	virtual long double d2cdy2(const long double& x, const long double& y);

	virtual int boundconds_classic(const long double& b1, const long double& b2,
			const long double eps);
	virtual int boundconds(const Interval<T>& B1, const Interval<T>& B2,
			const long double eps);

	virtual long double ExactSol(long double x, long double y);
	virtual long double GetConstM();
	virtual long double GetConstN();
	virtual long double GetConstP();
	virtual long double GetConstQ();
	virtual long double GetConstR();
	virtual long double GetConstS();
	virtual long double GetConstT();
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
long double BoundaryConditions<T>::f(const T& x, const T& y) {
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
long double BoundaryConditions<T>::psi(const T& x, const T& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::omega(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::a1(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::da1dx(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::da1dy(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::d2a1dx2(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::d2a1dy2(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::a2(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::da2dx(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::da2dy(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::d2a2dx2(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::d2a2dy2(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::c(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::dcdx(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::dcdy(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::d2cdx2(const long double& x,
		const long double& y) {
	return 0;
}

template<typename T>
long double BoundaryConditions<T>::d2cdy2(const long double& x,
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
Interval<T> BoundaryConditions<T>::A1(const Interval<T>& ix,
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
inline long double BoundaryConditions<T>::GetConstR() {
	return 0.0L;
}

template<typename T>
inline long double BoundaryConditions<T>::GetConstS() {
	return 0.0L;
}

template<typename T>
inline long double BoundaryConditions<T>::GetConstT() {
	return 0.0L;
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
inline Interval<T> BoundaryConditions<T>::A2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::D2A1DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::D2A1DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::D2A2DX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::D2A2DY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::DCDX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::DCDY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::D2CDX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::D2CDY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::DFDX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::DFDY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::D2FDX2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::D2FDY2(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::DA1DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::DA1DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::DA2DX(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
inline Interval<T> BoundaryConditions<T>::DA2DY(const Interval<T>& ix,
		const Interval<T>& iy, int& st) {
	Interval<T> r = { 0, 0 };
	st = 0;

	return r;
}

template<typename T>
void BoundaryConditions<T>::SetArithmeticMode(int mode) {

}

} /* namespace interval_arithmetic */
#endif /* BOUNDARYCONDITIONS_H_ */
