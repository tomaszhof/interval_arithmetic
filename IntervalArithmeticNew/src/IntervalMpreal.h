///*
// * IntervalMpreal.h
// *
// *  Created on: 28 lip 2015
// *      Author: tomhof
// */
//#include "Interval.h"
//#include <mpfr.h>
//#include <mpreal.h>
//
//using namespace std;
//using namespace mpfr;
//
//namespace interval_arithmetic {
//
//template<> IAMode Interval<mpreal>::mode = PINT_MODE;
//template<> IAPrecision Interval<mpreal>::precision = LONGDOUBLE_PREC;
//template<> IAOutDigits Interval<mpreal>::outdigits = LONGDOUBLE_DIGITS;
//
//
//template<>
//Interval<mpreal>::~Interval()
//{
//}
//
//template<>
//inline Interval<mpreal>::Interval(mpreal a, mpreal b) {
//	this->a = a;
//	this->b = b;
//}
//
//template<>
//inline IAMode Interval<mpreal>::GetMode() {
//	return Interval<mpreal>::mode;
//}
//
//template<>
//mpreal DIntWidth(const Interval<mpreal>& x) {
//	mpreal w1, w2;
//
//	mpreal::set_default_rnd(MPFR_RNDU);
//	w1 = x.b - x.a;
//	if (w1 < 0)
//		w1 = -w1;
//	mpreal::set_default_rnd(MPFR_RNDD);
//	w2 = x.b - x.a;
//	if (w2 < 0)
//		w2 = -w2;
//	mpreal::set_default_rnd(MPFR_RNDN);
//	if (w1 > w2)
//		return w1;
//	else
//		return w2;
//}
//
//template<>
//mpreal IntWidth(const Interval<mpreal>& x) {
//	mpreal::set_default_rnd(MPFR_RNDU);
//	mpreal w = x.b - x.a;
//	mpreal::set_default_rnd(MPFR_RNDN);
//	return w;
//}
//
//template<>
//inline mpreal Interval<mpreal>::GetWidth() {
//	Interval<mpreal> x(this->a, this->b);
//	switch (mode) {
//	case PINT_MODE:
//		return IntWidth(x);
//	case DINT_MODE:
//		return DIntWidth(x);
//	default:
//		return IntWidth(x);
//	}
//}
//
//template<>
//inline Interval<mpreal> Interval<mpreal>::operator =(const Interval<mpreal>& i) {
//	this->a = i.a;
//	this->b = i.b;
//
//	return *this;
//}
//
//template<>
//Interval<mpreal> IAdd(const Interval<mpreal>& x, const Interval<mpreal>& y) {
//	Interval<mpreal> r;
//	mpreal::set_default_rnd(MPFR_RNDD);
//	r.a = x.a + y.a;
//	mpreal::set_default_rnd(MPFR_RNDU);
//	r.b = x.b + y.b;
//	mpreal::set_default_rnd(MPFR_RNDN);
//	return r;
//}
//
//
//template<>
//Interval<mpreal> IDiv(const Interval<mpreal>& x, const Interval<mpreal>& y) {
//	Interval<mpreal> r;
//	mpreal x1y1, x1y2, x2y1, t;
//
//	if ((y.a <= 0) && (y.b >= 0)) {
//		throw runtime_error("Division by an interval containing 0.");
//	} else {
//		mpreal::set_default_rnd(MPFR_RNDD);
//		x1y1 = x.a / y.a;
//		x1y2 = x.a / y.b;
//		x2y1 = x.b / y.a;
//		r.a = x.b / y.b;
//		t = r.a;
//		if (x2y1 < t)
//			r.a = x2y1;
//		if (x1y2 < t)
//			r.a = x1y2;
//		if (x1y1 < t)
//			r.a = x1y1;
//
//		mpreal::set_default_rnd(MPFR_RNDU);
//		x1y1 = x.a / y.a;
//		x1y2 = x.a / y.b;
//		x2y1 = x.b / y.a;
//
//		r.b = x.b / y.b;
//		t = r.b;
//		if (x2y1 > t)
//			r.b = x2y1;
//		if (x1y2 > t)
//			r.b = x1y2;
//		if (x1y1 > t)
//			r.b = x1y1;
//
//	}
//	mpreal::set_default_rnd(MPFR_RNDN);
//	return r;
//}
//
//template<>
//Interval<mpreal> ISub(const Interval<mpreal>& x, const Interval<mpreal>& y) {
//	Interval<mpreal> r;
//	mpreal::set_default_rnd(MPFR_RNDD);
//	r.a = x.a - y.b;
//	mpreal::set_default_rnd(MPFR_RNDU);
//	r.b = x.b - y.a;
//	mpreal::set_default_rnd(MPFR_RNDN);
//	return r;
//}
//
//template<>
//Interval<mpreal> IMul(const Interval<mpreal>& x, const Interval<mpreal>& y) {
//	Interval<mpreal> r(0, 0);
//	mpreal x1y1, x1y2, x2y1;
//
//	mpreal::set_default_rnd(MPFR_RNDD);
//	x1y1 = x.a * y.a;
//	x1y2 = x.a * y.b;
//	x2y1 = x.b * y.a;
//	r.a = x.b * y.b;
//	if (x2y1 < r.a)
//		r.a = x2y1;
//	if (x1y2 < r.a)
//		r.a = x1y2;
//	if (x1y1 < r.a)
//		r.a = x1y1;
//
//	mpreal::set_default_rnd(MPFR_RNDU);
//	x1y1 = x.a * y.a;
//	x1y2 = x.a * y.b;
//	x2y1 = x.b * y.a;
//
//	r.b = x.b * y.b;
//	if (x2y1 > r.b)
//		r.b = x2y1;
//	if (x1y2 > r.b)
//		r.b = x1y2;
//	if (x1y1 > r.b)
//		r.b = x1y1;
//	mpreal::set_default_rnd(MPFR_RNDN);
//	return r;
//}
//
//template<>
//Interval<mpreal> ISin(const Interval<mpreal>& x) {
//	bool is_even, finished;
//	int k;
//	int st = 0;
//	Interval<mpreal> d, s, w, w1, x2;
//	mpreal eps = nextabove(0);
//
//	if (x.a > x.b)
//		st = 1;
//	else {
//		s = x;
//		w = x;
//		x2 = IMul(x, x);
//		k = 1;
//		is_even = true;
//		finished = false;
//		st = 0;
//
//		do {
//			d.a = (k + 1) * (k + 2);
//			d.b = d.a;
//			s = IMul(s, IDiv(x2, d));
//			if (is_even)
//				w1 = ISub(w, s);
//			else
//				w1 = IAdd(w, s);
//			if ((w.a != 0) && (w.b != 0)) {
//				if ((abs(w.a - w1.a) / abs(w.a) < eps)
//						&& (abs(w.b - w1.b) / abs(w.b) < eps))
//					finished = true;
//				else
//					;
//			} else if ((w.a == 0) && (w.b != 0)) {
//				if ((abs(w.a - w1.a) < eps)
//						&& (abs(w.b - w1.b) / abs(w.b) < eps))
//					finished = true;
//				else
//					;
//			} else if (w.a != 0) {
//				if ((abs(w.a - w1.a) / abs(w.a) < eps)
//						& (abs(w.b - w1.b) < eps))
//					finished = true;
//				else if ((abs(w.a - w1.a) < eps) & (abs(w.b - w1.b) < eps))
//					finished = true;
//			}
//
//			if (finished) {
//				if (w1.b > 1) {
//					w1.b = 1;
//					if (w1.a > 1)
//						w1.a = 1;
//				}
//				if (w1.a < -1) {
//					w1.a = -1;
//					if (w1.b < -1)
//						w1.b = -1;
//				}
//				return w1;
//			} else {
//				w = w1;
//				k = k + 2;
//				is_even = !is_even;
//				if ((w.a <= 0.0)&&(w.b >=0.0))
//				{
//					finished = true;
//					w = {0,0};
//					return w;
//				}
//			}
//		} while (!(finished || (k > INT_MAX / 2)));
//	}
//	if (!finished)
//		st = 2;
//
//	Interval<mpreal> r;
//	r.a = 0;
//	r.b = 0;
//	return r;
//}
//
//template<>
//Interval<mpreal> ICos(const Interval<mpreal>& x) {
//	bool is_even, finished;
//	int k;
//	int st = 0;
//	Interval<mpreal> d, c, w, w1, x2;
//	mpreal eps = nextabove(0);
//
//	if (x.a > x.b)
//		st = 1;
//	else {
//		c.a = 1;
//		c.b = 1;
//		w = c;
//		x2 = IMul(x, x);
//		k = 1;
//		is_even = true;
//		finished = false;
//		st = 0;
//
//		do {
//			d.a = k * (k + 1);
//			d.b = d.a;
//			c = IMul(c, IDiv(x2, d));
//			if (is_even)
//				w1 = ISub(w, c);
//			else
//				w1 = IAdd(w, c);
//
//			if ((w.a != 0) && (w.b != 0)) {
//				if ((abs(w.a - w1.a) / abs(w.a) < eps)
//						&& (abs(w.b - w1.b) / abs(w.b) < eps))
//					finished = true;
//				else
//					;
//			} else if ((w.a == 0) && (w.b != 0)) {
//				if ((abs(w.a - w1.a) < eps)
//						&& (abs(w.b - w1.b) / abs(w.b) < eps))
//					finished = true;
//				else
//					;
//			}
//
//			else if (w.a != 0) {
//				if ((abs(w.a - w1.a) / abs(w.a) < eps)
//						& (abs(w.b - w1.b) < eps))
//					finished = true;
//				else if ((abs(w.a - w1.a) < eps) & (abs(w.b - w1.b) < eps))
//					finished = true;
//			}
//
//			if (finished) {
//				if (w1.b > 1) {
//					w1.b = 1;
//					if (w1.a > 1)
//						w1.a = 1;
//				}
//				if (w1.a < -1) {
//					w1.a = -1;
//					if (w1.b < -1)
//						w1.b = -1;
//				}
//				return w1;
//			} else {
//				w = w1;
//				k = k + 2;
//				is_even = !is_even;
//			}
//		} while (!(finished || (k > INT_MAX / 2)));
//	}
//	if (!finished)
//		st = 2;
//
//	Interval<mpreal> r;
//	r.a = 0;
//	r.b = 0;
//	return r;
//}
//
//template<>
//Interval<mpreal> IExp(const Interval<mpreal>& x) {
//	bool finished;
//	int k;
//	int st = 0;
//	mpreal eps = nextabove(0);
//
//	Interval<mpreal> d, e, w, w1;
//	if (x.a > x.b)
//		st = 1;
//	else {
//		e.a = 1;
//		e.b = 1;
//		w = e;
//		k = 1;
//		finished = false;
//		st = 0;
//		do {
//			d.a = k;
//			d.b = k;
//			e = IMul(e, IDiv(x, d));
//			w1 = IAdd(w, e);
//			if ((abs(w.a - w1.a) / abs(w.a) < eps)
//					&& (abs(w.b - w1.b) / abs(w.b) < eps)) {
//				finished = true;
//				return w1;
//			} else {
//				w = w1;
//				k = k + 1;
//			}
//		} while (!(finished || (k > INT_MAX / 2)));
//		if (!finished)
//			st = 2;
//	}
//	Interval<mpreal> r;
//	r.a = 0;
//	r.b = 0;
//	return r;
//}
//
//template<>
//Interval<mpreal> ISqr(const Interval<mpreal>& x, int & st) {
//	mpreal minx, maxx;
//	Interval<mpreal> r;
//	r.a = 0;
//	r.b = 0;
//	if (x.a > x.b)
//		st = 1;
//	else {
//		st = 0;
//		if ((x.a <= 0) && (x.b >= 0))
//			minx = 0;
//		else if (x.a > 0)
//			minx = x.a;
//		else
//			minx = x.b;
//		if (abs(x.a) > mpfr::abs(x.b))
//			maxx = mpfr::abs(x.a);
//		else
//			maxx = mpfr::abs(x.b);
//		mpreal::set_default_rnd(MPFR_RNDD);
//		r.a = minx * minx;
//		mpreal::set_default_rnd(MPFR_RNDU);
//		r.b = maxx * maxx;
//		mpreal::set_default_rnd(MPFR_RNDN);
//	}
//	return r;
//}
//
//Interval<mpreal> DIAdd(const Interval<mpreal>& x, const Interval<mpreal>& y) {
//	Interval<mpreal> z1, z2;
//	if ((x.a <= x.b) && (y.a <= y.b)) {
//		return IAdd(x, y);
//	} else {
//		mpreal::set_default_rnd(MPFR_RNDD);
//		z1.a = x.a + y.a;
//		z2.b = x.b + y.b;
//		mpreal::set_default_rnd(MPFR_RNDU);
//		z1.b = x.b + y.b;
//		z2.a = x.a + y.a;
//		mpreal::set_default_rnd(MPFR_RNDN);
//		if (z1.GetWidth() >= z2.GetWidth())
//			return z1;
//		else
//			return z2;
//	}
//}
//
//template<>
//Interval<mpreal> DISub(const Interval<mpreal>& x, const Interval<mpreal>& y) {
//	Interval<mpreal> z1, z2;
//	if ((x.a <= x.b) && (y.a <= y.b)) {
//		return ISub(x, y);
//	} else {
//		mpreal::set_default_rnd(MPFR_RNDD);
//		z1.a = x.a - y.b;
//		z2.b = x.b - y.a;
//		mpreal::set_default_rnd(MPFR_RNDU);
//		z1.b = x.b - y.a;
//		z2.a = x.a - y.b;
//		mpreal::set_default_rnd(MPFR_RNDN);
//		if (z1.GetWidth() >= z2.GetWidth())
//			return z1;
//		else
//			return z2;
//	}
//}
//
//template<>
//Interval<mpreal> DIMul(const Interval<mpreal>& x, const Interval<mpreal>& y) {
//	Interval<mpreal> z1, z2, r;
//	mpreal z;
//	bool xn, xp, yn, yp, zero;
//
//	if ((x.a <= x.b) && (y.a <= y.b))
//		r = IMul(x, y);
//	else {
//		xn = (x.a < 0) and (x.b < 0);
//		xp = (x.a > 0) and (x.b > 0);
//		yn = (y.a < 0) and (y.b < 0);
//		yp = (y.a > 0) and (y.b > 0);
//		zero = false;
//		// A, B in H-T
//		if ((xn || xp) && (yn || yp))
//			if (xp && yp) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.a * y.a;
//				z2.b = x.b * y.b;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.b * y.b;
//				z2.a = x.a * y.a;
//			} else if (xp && yn) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.b * y.a;
//				z2.b = x.a * y.b;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.a * y.b;
//				z2.a = x.b * y.a;
//			} else if (xn && yp) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.a * y.b;
//				z2.b = x.b * y.a;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.b * y.a;
//				z2.a = x.a * y.b;
//			} else {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.b * y.b;
//				z2.b = x.a * y.a;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.a * y.a;
//				z2.a = x.b * y.b;
//			}
//		// A in H-T, B in T
//		else if ((xn || xp)
//				&& (((y.a <= 0) && (y.b >= 0)) || ((y.a >= 0) && (y.b <= 0))))
//			if (xp && (y.a <= y.b)) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.b * y.a;
//				z2.b = x.b * y.b;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.b * y.b;
//				z2.a = x.b * y.a;
//			} else if (xp && (y.a > y.b)) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.a * y.a;
//				z2.b = x.a * y.b;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.a * y.b;
//				z2.a = x.a * y.a;
//			} else if (xn && (y.a <= y.b)) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.a * y.b;
//				z2.b = x.a * y.a;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.a * y.a;
//				z2.a = x.a * y.b;
//			} else {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.b * y.b;
//				z2.b = x.b * y.a;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.b * y.a;
//				z2.a = x.b * y.b;
//			}
//		// A in T, B in H-T
//		else if ((((x.a <= 0) && (x.b >= 0)) || ((x.a >= 0) && (x.b <= 0)))
//				&& (yn || yp))
//			if ((x.a <= x.b) && yp) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.a * y.b;
//				z2.b = x.b * y.b;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.b * y.b;
//				z2.a = x.a * y.b;
//			} else if ((x.a <= 0) && yn) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.b * y.a;
//				z2.b = x.a * y.a;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.a * y.a;
//				z2.a = x.b * y.a;
//			} else if ((x.a > x.b) && yp) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.a * y.a;
//				z2.b = x.b * y.a;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.b * y.a;
//				z2.a = x.a * y.a;
//			} else {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.b * y.b;
//				z2.b = x.a * y.b;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.a * y.b;
//				z2.a = x.b * y.b;
//			}
//		// A, B in Z-
//		else if ((x.a >= 0) && (x.b <= 0) && (y.a >= 0) && (y.b <= 0)) {
//			mpreal::set_default_rnd(MPFR_RNDD);
//			z1.a = x.a * y.a;
//			z = x.b * y.b;
//			if (z1.a < z)
//				z1.a = z;
//			z2.b = x.a * y.b;
//			z = x.b * y.a;
//			if (z < z2.b)
//				z2.b = z;
//			mpreal::set_default_rnd(MPFR_RNDU);
//			z1.b = x.a * y.b;
//			z = x.b * y.a;
//			if (z < z1.b)
//				z1.b = z;
//			z2.a = x.a * y.a;
//			z = x.b * y.b;
//			if (z2.a < z)
//				z2.a = z;
//		}
//		// A in Z and B in Z- or A in Z- and B in Z
//		else
//			zero = true;
//		if (zero) {
//			r.a = 0;
//			r.b = 0;
//		} else if (z1.GetWidth() >= z2.GetWidth())
//			r = z1;
//		else
//			r = z2;
//	}
//
//	mpreal::set_default_rnd(MPFR_RNDN);
//	return r;
//}
//
//template<>
//Interval<mpreal> DIDiv(const Interval<mpreal>& x, const Interval<mpreal>& y) {
//	Interval<mpreal> z1, z2, r;
//	bool xn, xp, yn, yp, zero;
//
//	if ((x.a <= x.b) && (y.a <= y.b))
//		r = IDiv(x, y);
//	else {
//		xn = (x.a < 0) && (x.b < 0);
//		xp = (x.a > 0) && (x.b > 0);
//		yn = (y.a < 0) && (y.b < 0);
//		yp = (y.a > 0) && (y.b > 0);
//		zero = false;
//		// A, B in H-T
//		if ((xn || xp) && (yn || yp))
//			if (xp && yp) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.a / y.b;
//				z2.b = x.b / y.a;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.b / y.a;
//				z2.a = x.a / y.b;
//			} else if (xp && yn) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.b / y.b;
//				z2.b = x.a / y.a;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.a / y.a;
//				z2.a = x.b / y.b;
//			} else if (xn && yp) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.a / y.a;
//				z2.b = x.b / y.b;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.b / y.b;
//				z2.a = x.a / y.a;
//			} else {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.b / y.a;
//				z2.b = x.a / y.b;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.a / y.b;
//				z2.a = x.b / y.a;
//			}
//		// A in T, B in H-T
//		else if (((x.a <= 0) && (x.b >= 0))
//				|| (((x.a >= 0) && (x.b <= 0)) && (yn || yp)))
//			if ((x.a <= x.b) && yp) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.a / y.a;
//				z2.b = x.b / y.a;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.b / y.a;
//				z2.a = x.a / y.a;
//			} else if ((x.a <= x.b) && yn) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.b / y.b;
//				z2.b = x.a / y.b;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.a / y.b;
//				z2.a = x.b / y.b;
//			} else if ((x.a > x.b) && yp) {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.a / y.b;
//				z2.b = x.b / y.b;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.b / y.b;
//				z2.a = x.a / y.b;
//			} else {
//				mpreal::set_default_rnd(MPFR_RNDD);
//				z1.a = x.b / y.a;
//				z2.b = x.a / y.a;
//				mpreal::set_default_rnd(MPFR_RNDU);
//				z1.b = x.a / y.a;
//				z2.a = x.b / y.a;
//			}
//		else
//			zero = true;
//		if (zero)
//			throw runtime_error("Division by an interval containing 0.");
//		else if (z1.GetWidth() >= z2.GetWidth())
//			r = z1;
//		else
//			r = z2;
//		mpreal::set_default_rnd(MPFR_RNDN);
//	}
//	return r;
//}
//
//template<>
//Interval<mpreal> DISin(const Interval<mpreal>& x) {
//	bool is_even, finished;
//	int k;
//	int st = 0;
//	Interval<mpreal> d, s, w, w1, x2;
//	mpreal eps = nextabove(0);
//	if (x.a > x.b)
//		st = 1;
//	else {
//		s = x;
//		w = x;
//		x2 = DIMul(x, x);
//		k = 1;
//		is_even = true;
//		finished = false;
//		st = 0;
//
//		do {
//			d.a = (k + 1) * (k + 2);
//			d.b = d.a;
//			s = DIMul(s, DIDiv(x2, d));
//			if (is_even)
//				w1 = DISub(w, s);
//			else
//				w1 = DIAdd(w, s);
//			if ((w.a != 0) && (w.b != 0)) {
//				if ((abs(w.a - w1.a) / abs(w.a) < eps)
//						&& (abs(w.b - w1.b) / abs(w.b) < eps))
//					finished = true;
//				else
//					;
//			} else if ((w.a == 0) && (w.b != 0)) {
//				if ((abs(w.a - w1.a) < eps)
//						&& (abs(w.b - w1.b) / abs(w.b) < eps))
//					finished = true;
//				else
//					;
//			}
//
//			else if (w.a != 0) {
//				if ((abs(w.a - w1.a) / abs(w.a) < eps)
//						& (abs(w.b - w1.b) < eps))
//					finished = true;
//				else if ((abs(w.a - w1.a) < eps) & (abs(w.b - w1.b) < eps))
//					finished = true;
//			}
//
//			if (finished) {
//				if (w1.b > 1) {
//					w1.b = 1;
//					if (w1.a > 1)
//						w1.a = 1;
//				}
//				if (w1.a < -1) {
//					w1.a = -1;
//					if (w1.b < -1)
//						w1.b = -1;
//				}
//				return w1;
//			} else {
//				w = w1;
//				k = k + 2;
//				is_even = !is_even;
//			}
//		} while (!(finished || (k > INT_MAX / 2)));
//	}
//	if (!finished)
//		st = 2;
//
//	Interval<mpreal> r;
//	r.a = 0;
//	r.b = 0;
//	return r;
//}
//
//template<>
//Interval<mpreal> DICos(const Interval<mpreal>& x, int & st) {
//	bool is_even, finished;
//	int k;
//	Interval<mpreal> d, c, w, w1, x2;
//	mpreal eps = nextabove(0);
//
//	if (x.a > x.b)
//		st = 1;
//	else {
//		c.a = 1;
//		c.b = 1;
//		w = c;
//		x2 = DIMul(x, x);
//		k = 1;
//		is_even = true;
//		finished = false;
//		st = 0;
//
//		do {
//			d.a = k * (k + 1);
//			d.b = d.a;
//			c = DIMul(c, DIDiv(x2, d));
//			if (is_even)
//				w1 = DISub(w, c);
//			else
//				w1 = DIAdd(w, c);
//
//			if ((w.a != 0) && (w.b != 0)) {
//				if ((abs(w.a - w1.a) / abs(w.a) < eps)
//						&& (abs(w.b - w1.b) / abs(w.b) < eps))
//					finished = true;
//				else
//					;
//			} else if ((w.a == 0) && (w.b != 0)) {
//				if ((abs(w.a - w1.a) < eps)
//						&& (abs(w.b - w1.b) / abs(w.b) < eps))
//					finished = true;
//				else
//					;
//			}
//
//			else if (w.a != 0) {
//				if ((abs(w.a - w1.a) / abs(w.a) < eps)
//						& (abs(w.b - w1.b) < eps))
//					finished = true;
//				else if ((abs(w.a - w1.a) < eps) & (abs(w.b - w1.b) < eps))
//					finished = true;
//			}
//
//			if (finished) {
//				if (w1.b > 1) {
//					w1.b = 1;
//					if (w1.a > 1)
//						w1.a = 1;
//				}
//				if (w1.a < -1) {
//					w1.a = -1;
//					if (w1.b < -1)
//						w1.b = -1;
//				}
//				return w1;
//			} else {
//				w = w1;
//				k = k + 2;
//				is_even = !is_even;
//			}
//		} while (!(finished || (k > INT_MAX / 2)));
//	}
//	if (!finished)
//		st = 2;
//
//	Interval<mpreal> r;
//	r.a = 0;
//	r.b = 0;
//	return r;
//}
//
//template<>
//Interval<mpreal> DIExp(const Interval<mpreal>& x) {
//	bool finished;
//	int k;
//	int st = 0;
//	Interval<mpreal> d, e, w, w1;
//	mpreal eps = nextabove(0);
//
//	if (x.a > x.b)
//		st = 1;
//	else {
//		e.a = 1;
//		e.b = 1;
//		w = e;
//		k = 1;
//		finished = false;
//		st = 0;
//		do {
//			d.a = k;
//			d.b = k;
//			e = IMul(e, DIDiv(x, d));
//			w1 = DIAdd(w, e);
//			if ((abs(w.a - w1.a) / abs(w.a) < eps)
//					&& (abs(w.b - w1.b) / abs(w.b) < eps)) {
//				finished = true;
//				return w1;
//			} else {
//				w = w1;
//				k = k + 1;
//			}
//		} while (!(finished || (k > INT_MAX / 2)));
//		if (!finished)
//			st = 2;
//	}
//	Interval<mpreal> r;
//	r.a = 0;
//	r.b = 0;
//	return r;
//}
//
//template<>
//Interval<mpreal> DISqr(const Interval<mpreal>& x) {
//	mpreal minx, maxx;
//	int st = 0;
//	Interval<mpreal> r;
//
//	r.a = 0;
//	r.b = 0;
//	if (x.a > x.b)
//		st = 1;
//	else {
//		st = 0;
//		if ((x.a <= 0) && (x.b >= 0))
//			minx = 0;
//		else if (x.a > 0)
//			minx = x.a;
//		else
//			minx = x.b;
//		if (abs(x.a) > abs(x.b))
//			maxx = abs(x.a);
//		else
//			maxx = abs(x.b);
//		mpreal::set_default_rnd(MPFR_RNDD);
//		r.a = minx * minx;
//		mpreal::set_default_rnd(MPFR_RNDU);
//		r.b = maxx * maxx;
//		mpreal::set_default_rnd(MPFR_RNDN);
//	}
//	return r;
//}
//
//template<>
//inline Interval<mpreal> Interval<mpreal>::operator +(const Interval<mpreal>& y) {
//	Interval<mpreal> x(this->a, this->b);
//	Interval<mpreal> r = {0, 0};
//	switch (mode) {
//	case PINT_MODE:
//		r = IAdd(x, y);
//		break;
//	case DINT_MODE:
//		r = DIAdd(x, y);
//		break;
//	default:
//		r = IAdd(x, y);
//		break;
//	}
//
//	return r;
//}
//
//template<>
//inline Interval<mpreal> operator +(Interval<mpreal> x, const Interval<mpreal>& y) {
//	switch (Interval<mpreal>::mode) {
//	case PINT_MODE:
//		return  IAdd(x, y);
//	case DINT_MODE:
//		return  DIAdd(x, y);
//	default:
//		return IAdd(x, y);
//	}
//}
//
//template<>
//inline Interval<mpreal> Interval<mpreal>::operator -(const Interval<mpreal>& y) {
//	Interval<mpreal> x(this->a, this->b);
//	Interval<mpreal> r = {0, 0};
//	switch (mode) {
//	case PINT_MODE:
//		r = ISub(x, y);
//		break;
//	case DINT_MODE:
//		r = DISub(x, y);
//		break;
//	default:
//		r = ISub(x, y);
//		break;
//	}
//
//	return r;
//}
//
//template<>
//inline Interval<mpreal> operator -(Interval<mpreal> x, const Interval<mpreal>& y) {
//	switch (Interval<mpreal>::mode) {
//	case PINT_MODE:
//		return  ISub(x, y);
//	case DINT_MODE:
//		return  DISub(x, y);
//	default:
//		return ISub(x, y);
//	}
//}
//
//template<>
//inline Interval<mpreal> Interval<mpreal>::operator *(const Interval<mpreal>& y) {
//	Interval<mpreal> x(this->a, this->b);
//	Interval<mpreal> r = {0, 0};
//	switch (mode) {
//	case PINT_MODE:
//		r = IMul(x, y);
//		break;
//	case DINT_MODE:
//		r = DIMul(x, y);
//		break;
//	default:
//		r = IMul(x, y);
//		break;
//	}
//
//	return r;
//}
//
//template<>
//inline Interval<mpreal> operator *(Interval<mpreal> x, const Interval<mpreal>& y) {
//	switch (Interval<mpreal>::mode) {
//	case PINT_MODE:
//		return IMul(x, y);
//	case DINT_MODE:
//		return DIMul(x, y);
//	default:
//		return IMul(x, y);
//	}
//}
//
//template<>
//inline Interval<mpreal> Interval<mpreal>::operator /(const Interval<mpreal>& y) {
//	Interval<mpreal> x(this->a, this->b);
//	Interval<mpreal> r = {0, 0};
//	switch (mode) {
//	case PINT_MODE:
//		r = IDiv(x, y);
//		break;
//	case DINT_MODE:
//		r = DIDiv(x, y);
//		break;
//	default:
//		r = IDiv(x, y);
//		break;
//	}
//
//	return r;
//}
//
//template<>
//inline Interval<mpreal> operator /(Interval<mpreal> x, const Interval<mpreal>& y) {
//	switch (Interval<mpreal>::mode) {
//	case PINT_MODE:
//		return IDiv(x, y);
//	case DINT_MODE:
//		return DIDiv(x, y);
//	default:
//		return IDiv(x, y);
//	}
//}
//template<>
//inline void Interval<mpreal>::SetPrecision(IAPrecision p) {
//	Interval<mpreal>::precision = p;
//}
//
//template<>
//inline IAPrecision Interval<mpreal>::GetPrecision(IAPrecision p) {
//	return Interval<mpreal>::precision;
//}
//
//template<>
//inline void Interval<mpreal>::SetOutDigits(IAOutDigits o) {
//	Interval<mpreal>::outdigits = LONGDOUBLE_DIGITS;
//}
//
//template<>
//inline IAOutDigits Interval<mpreal>::GetOutDigits(IAOutDigits o) {
//	return Interval<mpreal>::outdigits;
//}
//
//template<>
//inline Interval<mpreal> Interval<mpreal>::IntRead(const string& sa) {
//	Interval<mpreal> r;
//	mpfr_t ropl;
//	mpfr_init2(ropl, precision);
//	mpfr_set_str(ropl, sa.c_str(), 10, MPFR_RNDD);
//
//	mpfr_t ropr;
//	mpfr_init2(ropr, precision);
//	mpfr_set_str(ropr, sa.c_str(), 10, MPFR_RNDU);
//
//	r.a = ropl;
//	r.b = ropr;
//	return r;
//}
//
//// TODO: perform changes form mpreal type template specialization
//
//template<>
//inline void Interval<mpreal>::IEndsToStrings(string& left, string& right) {
//	mpfr_t rop;
//	mpfr_exp_t exponent;
//	mpfr_init2(rop, precision);
//	char* str = NULL;
//	char *buffer = new char(precision + 3);
//	mpfr_set(rop, this->a.mpfr_ptr(), MPFR_RNDD);
//	mpfr_get_str(buffer, &exponent, 10, outdigits, rop, MPFR_RNDD);
//	str = buffer;
//
//	stringstream ss;
//	int prec = std::numeric_limits<mpreal>::digits10();
//	ss.setf(std::ios_base::scientific);
//	bool minus = (str[0] == '-');
//	int splitpoint = minus ? 1 : 0;
//	string sign = minus ? "-" : "";
//
//	ss << std::setprecision(prec) << sign << str[splitpoint] << "."
//			<< &str[splitpoint + 1] << "E" << exponent - 1;
//	left = ss.str();
//	ss.str(std::string());
//
//	mpfr_set(rop, this->b.mpfr_ptr(), MPFR_RNDU);
//	mpfr_get_str(buffer, &exponent, 10, outdigits, rop, MPFR_RNDU);
//	str = buffer;
//	splitpoint = (str[0] == '-') ? 1 : 0;
//	ss << std::setprecision(prec) << sign << str[splitpoint] << "."
//			<< &str[splitpoint + 1] << "E" << exponent - 1;
//	right = ss.str();
//	ss.clear();
//}
//
//template<>
//inline Interval<mpreal> Interval<mpreal>::Projection() {
//	Interval<mpreal> x(this->a, this->b);
//	Interval<mpreal> r;
//	r = x;
//	if (x.a > x.b) {
//		r.a = x.b;
//		r.b = x.a;
//	}
//	return r;
//}
//
//template<>
//inline Interval<mpreal> Interval<mpreal>::Opposite() {
//	Interval<mpreal> x(this->a, this->b);
//	Interval<mpreal> r;
//	r.a = -x.a;
//	r.b = -x.b;
//	return r;
//}
//
//template<>
//inline Interval<mpreal> Interval<mpreal>::Inverse() {
//	Interval<mpreal> x(this->a, this->b);
//	Interval<mpreal> z1, z2;
//
//	mpreal::set_default_rnd(MPFR_RNDD);
//	z1.a = 1 / x.a;
//	z2.b = 1 / x.b;
//	mpreal::set_default_rnd(MPFR_RNDU);
//	z1.b = 1 / x.b;
//	z2.a = 1 / x.a;
//	mpreal::set_default_rnd(MPFR_RNDN);
//	if (DIntWidth(z1) >= DIntWidth(z2))
//		return z1;
//	else
//		return z2;
//}
//
//template<>
//inline mpreal Interval<mpreal>::LeftRead(const string& sa) {
//	Interval<mpreal> int_number;
//	int_number = IntRead(sa);
//	return int_number.a;
//}
//
//template<>
//inline mpreal Interval<mpreal>::RightRead(const string& sa) {
//	Interval<mpreal> int_number;
//	int_number = IntRead(sa);
//	return int_number.b;
//}
//
//template<>
//inline Interval<mpreal> Interval<mpreal>::ISqr2() {
//	string i2;
//	Interval<mpreal> r;
//	i2 = "1.414213562373095048";
//	r.a = LeftRead(i2);
//	i2 = "1.414213562373095049";
//	r.b = RightRead(i2);
//	return r;
//}
//
//template<>
//inline Interval<mpreal> Interval<mpreal>::ISqr3() {
//	string i2;
//	Interval<mpreal> r;
//	i2 = "1.732050807568877293";
//	r.a = LeftRead(i2);
//	i2 = "1.732050807568877294";
//	r.b = RightRead(i2);
//	return r;
//}
//
//template<>
//inline Interval<mpreal> Interval<mpreal>::IPi() {
//	string i2;
//	Interval<mpreal> r;
//	i2 = "3.141592653589793238";
//	r.a = LeftRead(i2);
//	i2 = "3.141592653589793239";
//	r.b = RightRead(i2);
//	return r;
//}
//
//template<>
//inline void Interval<mpreal>::Initialize() {
//	Interval<mpreal>::precision = LONGDOUBLE_PREC;
//	Interval<mpreal>::outdigits = LONGDOUBLE_DIGITS;
//}
//
//}
