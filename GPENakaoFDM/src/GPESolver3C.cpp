/*
 * GPESolver3C.cpp
 *
 *  Created on: 01.05.2019
 *      Author: numeric
 */

#include "GPESolver3C.h"

namespace interval_arithmetic {

template<typename T>
GPESolver3C<T>::GPESolver3C() {

}

template<typename T>
GPESolver3C<T>::~GPESolver3C() {
	// TODO Auto-generated destructor stub
}

template<typename T>
int GPESolver3C<T>::SetExample(int eid) {
	switch (eid) {
	case 1:
		bc = new ExampleGPE01<T>();
		break;
	case 2:
		bc = new ExampleGPE02<T>();
		break;
	case 3:
		bc = new ExampleGPE03<T>();
		break;
	case 4:
		bc = new ExampleGPE04<T>();
		break;
	case 5:
		bc = new ExampleGPE05<T>();
		break;
	case 6:
		bc = new ExampleGPE06<T>();
		break;
	case 7:
		bc = new ExampleGPE07<T>();
		break;
	case 8:
		bc = new ExampleGPE08<T>();
		break;
	case 11:
		bc = new ExampleGPE11<T>();
		break;
	default:
		bc = NULL;
		break;
	}

//		if (_example != NULL)
//			_example->SetArithmeticMode(arth_mode);
	return 0;
}

template<typename T>
T GPESolver3C<T>::alphax(T xi, T yj) {
	T result = 0.0;
	T h2d12 = this->h * this->h / 12.0;
	T k2d12 = this->k * this->k / 12.0;
	T tmp2A1M1 = 2.0 * 1.0 / bc->a1(xi, yj);
	T tmp2A2M1 = 2.0 * 1.0 / bc->a2(xi, yj);

	result = bc->a1(xi, yj)
			- h2d12
					* (tmp2A1M1 * bc->da1dx(xi, yj) * bc->da1dx(xi, yj)
							- bc->d2a1dx2(xi, yj) - bc->c(xi, yj));
	result = result
			- k2d12
					* (tmp2A2M1 * bc->da2dy(xi, yj) * bc->da1dy(xi, yj)
							- bc->d2a1dy2(xi, yj));
	return result;
}

template<typename T>
T GPESolver3C<T>::alphay(T xi, T yj) {
	T result = 0.0;
	T h2d12 = this->h * this->h / 12.0;
	T k2d12 = this->k * this->k / 12.0;
	T tmp2A1M1 = 2.0 * 1.0 / bc->a1(xi, yj);
	T tmp2A2M1 = 2.0 * 1.0 / bc->a2(xi, yj);

	result = bc->a2(xi, yj)
			- h2d12
					* (tmp2A1M1 * bc->da1dx(xi, yj) * bc->da2dx(xi, yj)
							- bc->d2a2dx2(xi, yj));
	result = result
			- k2d12
					* (tmp2A2M1 * bc->da2dy(xi, yj) * bc->da2dy(xi, yj)
							- bc->d2a2dy2(xi, yj) - bc->c(xi, yj));
	return result;
}

template<typename T>
T GPESolver3C<T>::betax(T xi, T yj) {
	T result = 0.0;
	T h2d6 = (-1.0) * this->h * this->h / 6.0;
	T tmpA1M1 = 1.0 / bc->a1(xi, yj);
	result = h2d6
			* (tmpA1M1 * bc->da1dx(xi, yj) * bc->c(xi, yj) - bc->dcdx(xi, yj));
	return result;
}

template<typename T>
T GPESolver3C<T>::betay(T xi, T yj) {
	T result = 0.0;
	T k2d6 = (-1.0) * this->k * this->k / 6.0;
	T tmpA2M1 = 1.0 / bc->a2(xi, yj);
	result = k2d6
			* (tmpA2M1 * bc->da2dy(xi, yj) * bc->c(xi, yj) - bc->dcdy(xi, yj));

	return result;
}

template<typename T>
T GPESolver3C<T>::gammaxy(T xi, T yj) {
	T result = 0.0;
	T h2d12 = this->h * this->h / 12.0;
	T k2d12 = this->k * this->k / 12.0;
	T tmpA1M1 = 2.0 / bc->a1(xi, yj);
	T tmpA2M1 = 2.0 / bc->a2(xi, yj);

	result = bc->c(xi, yj)
			- h2d12
					* (tmpA1M1 * bc->da1dx(xi, yj) * bc->dcdx(xi, yj)
							- bc->d2cdx2(xi, yj));
	result = result
			- k2d12
					* (tmpA2M1 * bc->da2dy(xi, yj) * bc->dcdy(xi, yj)
							- bc->d2cdy2(xi, yj));
	return result;
}

template<typename T>
Interval<T> GPESolver3C<T>::IAlphaX(Interval<T> xi, Interval<T> yj) {
	Interval<T> result = { 0.0L, 0.0L };
	int st = 0;
	Interval<T> h2d12 = this->ih2 / this->i12;
	Interval<T> k2d12 = this->ik2 / this->i12;
	Interval<T> tmp2A1M1 = i2 / bc->A1(xi, yj, st);
	Interval<T> tmp2A2M1 = i2 / bc->A2(xi, yj, st);

	result = bc->A1(xi, yj, st)
			- h2d12
					* (tmp2A1M1 * bc->DA1DX(xi, yj, st) * bc->DA1DX(xi, yj, st)
							- bc->D2A1DX2(xi, yj, st) - bc->C(xi, yj, st));
	result = result
			- k2d12
					* (tmp2A2M1 * bc->DA2DY(xi, yj, st) * bc->DA1DY(xi, yj, st)
							- bc->D2A1DY2(xi, yj, st));

	return result;
}

template<typename T>
Interval<T> GPESolver3C<T>::IAlphaY(Interval<T> xi, Interval<T> yj) {
	Interval<T> result = { 0.0L, 0.0L };
	int st = 0;
	Interval<T> h2d12 = this->ih2 / this->i12;
	Interval<T> k2d12 = this->ik2 / this->i12;
	Interval<T> tmp2A1M1 = i2 / bc->A1(xi, yj, st);
	Interval<T> tmp2A2M1 = i2 / bc->A2(xi, yj, st);

	Interval<T> DY2COMP = tmp2A2M1 * bc->DA2DY(xi, yj, st)
			* bc->DA2DY(xi, yj, st) - bc->D2A2DY2(xi, yj, st)
			- bc->C(xi, yj, st);
	Interval<T> DA1A2DX = tmp2A1M1 * bc->DA1DX(xi, yj, st)
			* bc->DA2DX(xi, yj, st) - bc->D2A2DX2(xi, yj, st);
	Interval<T> T_DA1DX = bc->DA1DX(xi, yj, st);
	Interval<T> T_DA2DX = bc->DA2DX(xi, yj, st);
	Interval<T> T_D2A2DX2 = bc->D2A2DX2(xi, yj, st);

	result = bc->A2(xi, yj, st)
			- h2d12
					* (tmp2A1M1 * bc->DA1DX(xi, yj, st) * bc->DA2DX(xi, yj, st)
							- bc->D2A2DX2(xi, yj, st));
	result = result
			- k2d12
					* (tmp2A2M1 * bc->DA2DY(xi, yj, st) * bc->DA2DY(xi, yj, st)
							- bc->D2A2DY2(xi, yj, st) - bc->C(xi, yj, st));
	return result;
}

template<typename T>
Interval<T> GPESolver3C<T>::IBetaX(Interval<T> xi, Interval<T> yj) {

	Interval<T> result = { 0.0L, 0.0L };
	int st = 0;
	Interval<T> h2d6 = this->ih * this->ih / i6;
	Interval<T> tmpA1M1 = i1 / bc->A1(xi, yj, st);
	result = h2d6
			* (bc->DCDX(xi, yj, st)
					- bc->DA1DX(xi, yj, st) * bc->C(xi, yj, st) * tmpA1M1);
	return result;
}

template<typename T>
Interval<T> GPESolver3C<T>::IBetaY(Interval<T> xi, Interval<T> yj) {
	Interval<T> result = { 0.0L, 0.0L };
	int st;
	Interval<T> k2d6 = this->ik * this->ik / this->i6;
	Interval<T> tmpA2M1 = i1 / bc->A2(xi, yj, st);
	result = k2d6
			* (bc->DCDY(xi, yj, st)
					- bc->DA2DY(xi, yj, st) * bc->C(xi, yj, st) * tmpA2M1);

	return result;
}

template<typename T>
Interval<T> GPESolver3C<T>::IGammaXY(Interval<T> xi, Interval<T> yj) {
	Interval<T> result = { 0.0L, 0.0L };
	int st;
	Interval<T> h2d12 = im1 * this->ih * this->ih / i12;
	Interval<T> k2d12 = im1 * this->ik * this->ik / i12;
	Interval<T> tmpA1M1 = i2 / bc->A1(xi, yj, st);
	Interval<T> tmpA2M1 = i2 / bc->A2(xi, yj, st);

	Interval<T> tmpW5 = im1 * h2d12
			* (tmpA1M1 * bc->DA1DX(xi, yj, st) * bc->DCDX(xi, yj, st)
					- bc->D2CDX2(xi, yj, st));
	Interval<T> tmpW6 = im1 * k2d12
			* (tmpA2M1 * bc->DA2DY(xi, yj, st) * bc->DCDY(xi, yj, st)
					- bc->D2CDY2(xi, yj, st));
//	cout << "W5_ij= [" << tmpW5.a << " ; " << tmpW5.b << "]" << endl;
//	cout << "W6_ij= [" << tmpW6.a << " ; " << tmpW6.b << "]" << endl;

	result = bc->C(xi, yj, st)
			+ h2d12
					* (tmpA1M1 * bc->DA1DX(xi, yj, st) * bc->DCDX(xi, yj, st)
							- bc->D2CDX2(xi, yj, st));
	result = result
			+ k2d12
					* (tmpA2M1 * bc->DA2DY(xi, yj, st) * bc->DCDY(xi, yj, st)
							- bc->D2CDY2(xi, yj, st));
	return result;
}

template<typename T>
Interval<T> GPESolver3C<T>::DIAlphaX(Interval<T> xi, Interval<T> yj) {
	Interval<T> result = { 0.0L, 0.0L };
	int st = 0;
	Interval<T> h2d12 = this->ih2 / this->i12;
	Interval<T> k2d12 = this->ik2 / this->i12;
	Interval<T> tmp2A1M1 = i2 / bc->A1(xi, yj, st);
	Interval<T> tmp2A2M1 = i2 / bc->A2(xi, yj, st);

	result = bc->A1(xi, yj, st)
			- h2d12
					* (tmp2A1M1 * bc->DA1DX(xi, yj, st) * bc->DA1DX(xi, yj, st)
							- bc->D2A1DX2(xi, yj, st) - bc->C(xi, yj, st));
	result = result
			- k2d12
					* (tmp2A2M1 * bc->DA2DY(xi, yj, st) * bc->DA1DY(xi, yj, st)
							- bc->D2A1DY2(xi, yj, st));

	return result;
}

template<typename T>
Interval<T> GPESolver3C<T>::DIAlphaY(Interval<T> xi, Interval<T> yj) {
	Interval<T> result = { 0.0L, 0.0L };
	int st = 0;
	Interval<T> h2d12 = this->ih2 / this->i12;
	Interval<T> k2d12 = this->ik2 / this->i12;
	Interval<T> tmp2A1M1 = i2 / bc->A1(xi, yj, st);
	Interval<T> tmp2A2M1 = i2 / bc->A2(xi, yj, st);

	Interval<T> DY2COMP = tmp2A2M1 * bc->DA2DY(xi, yj, st)
			* bc->DA2DY(xi, yj, st) - bc->D2A2DY2(xi, yj, st)
			- bc->C(xi, yj, st);
	Interval<T> DA1A2DX = tmp2A1M1 * bc->DA1DX(xi, yj, st)
			* bc->DA2DX(xi, yj, st) - bc->D2A2DX2(xi, yj, st);
	Interval<T> T_DA1DX = bc->DA1DX(xi, yj, st);
	Interval<T> T_DA2DX = bc->DA2DX(xi, yj, st);
	Interval<T> T_D2A2DX2 = bc->D2A2DX2(xi, yj, st);

	result = bc->A2(xi, yj, st)
			+ (h2d12
					* (tmp2A1M1 * bc->DA1DX(xi, yj, st) * bc->DA2DX(xi, yj, st)
							- bc->D2A2DX2(xi, yj, st))).Opposite();
	result =
			result
					+ (k2d12
							* (tmp2A2M1 * bc->DA2DY(xi, yj, st)
									* bc->DA2DY(xi, yj, st)
									- bc->D2A2DY2(xi, yj, st)
									- bc->C(xi, yj, st))).Opposite();
	return result;
}

template<typename T>
Interval<T> GPESolver3C<T>::DIBetaX(Interval<T> xi, Interval<T> yj) {

	Interval<T> result = { 0.0L, 0.0L };
	int st = 0;
	Interval<T> h2d6 = im1 * this->ih * this->ih / i6;
	Interval<T> tmpA1M1 = i1 / bc->A1(xi, yj, st);
	result = h2d6
			* (bc->DCDX(xi, yj, st)
					- bc->DA1DX(xi, yj, st) * bc->C(xi, yj, st) * tmpA1M1);
	return result.Opposite();
}

template<typename T>
Interval<T> GPESolver3C<T>::DIBetaY(Interval<T> xi, Interval<T> yj) {
	Interval<T> result = { 0.0L, 0.0L };
	int st;
	Interval<T> k2d6 = im1 * this->ik * this->ik / this->i6;
	Interval<T> tmpA2M1 = i1 / bc->A2(xi, yj, st);
	result = k2d6
			* (bc->DCDY(xi, yj, st)
					- bc->DA2DY(xi, yj, st) * bc->C(xi, yj, st) * tmpA2M1);

	return result.Opposite();
}

template<typename T>
Interval<T> GPESolver3C<T>::DIGammaXY(Interval<T> xi, Interval<T> yj) {
	Interval<T> result = { 0.0L, 0.0L };
	int st;
	Interval<T> h2d12 = this->ih * this->ih / i12;
	Interval<T> k2d12 = this->ik * this->ik / i12;
	Interval<T> tmpA1M1 = i2 / bc->A1(xi, yj, st);
	Interval<T> tmpA2M1 = i2 / bc->A2(xi, yj, st);

	Interval<T> tmpW5 = im1 * h2d12
			* (tmpA1M1 * bc->DA1DX(xi, yj, st) * bc->DCDX(xi, yj, st)
					- bc->D2CDX2(xi, yj, st));
	Interval<T> tmpW6 = im1 * k2d12
			* (tmpA2M1 * bc->DA2DY(xi, yj, st) * bc->DCDY(xi, yj, st)
					- bc->D2CDY2(xi, yj, st));
//	cout << "W5_ij= [" << tmpW5.a << " ; " << tmpW5.b << "]" << endl;
//	cout << "W6_ij= [" << tmpW6.a << " ; " << tmpW6.b << "]" << endl;

	result = bc->C(xi, yj, st)
			+ (h2d12
					* (tmpA1M1 * bc->DA1DX(xi, yj, st) * bc->DCDX(xi, yj, st)
							- bc->D2CDX2(xi, yj, st))).Opposite();
	result = result
			+ (k2d12
					* (tmpA2M1 * bc->DA2DY(xi, yj, st) * bc->DCDY(xi, yj, st)
							- bc->D2CDY2(xi, yj, st))).Opposite();
	return result;
}

template<typename T>
int GPESolver3C<T>::SolveFP() {
	int i, j, jh, j1, k, kh, l, lh, l1, l2, n1, n2, p, q, rh, st;
	long double h1, k1, h2, k2, hh1, kk1, max, s, tmpP, tmpQ, tmpR, tmpS;

	if (!Solver<T>::_initparams)
		throw runtime_error("Parameters not initialized!");

	int m = params.m;
	int n = params.n;
	long double alpha = params.alpha;
	long double beta = params.beta;
	long double delta = params.delta;
	long double gamma = params.gamma;
	long double eps = params.eps;

	st = 0;
	if ((n < 2) || (m < 2))
		st = 1;
//	else if ((alpha < 0) || (beta < 0))
//		st = 2;
	if (st == 0) {
		st = this->bc->boundconds_classic(bc->phi1(beta), bc->phi2(alpha), eps);
		if (st == 0) {
			st = bc->boundconds_classic(bc->phi2(gamma), bc->phi3(beta), eps);
			if (st == 0) {
				st = bc->boundconds_classic(bc->phi3(delta), bc->phi4(gamma),
						eps);
				if (st == 0)
					st = bc->boundconds_classic(bc->phi4(alpha),
							bc->phi1(delta), eps);
			}
		}
	}
	if (st == 0) {
		h1 = (gamma - alpha) / n;
		k1 = (delta - beta) / m;
		h2 = h1 * h1;
		k2 = k1 * k1;
		n1 = (n - 1) * (m - 1);
		n2 = n1 + 1;
		p = n2;

		long double* a1 = new long double[n1 + 1];
		long double* b1 = new long double[n1 + 1];
		int* r = new int[n1 + 1];
		long double* x = new long double[(n * m - n - m + 4)
				* (n * m - n - m + 4) / 4];

		for (int i = 1; i <= n2; i++)
			r[i - 1] = 0;
		k = 0;
		j = 0;

		do {
			k = k + 1;
			for (int i = 1; i <= n1; i++)
				a1[i - 1] = 0;
			j = j + 1;
			i = ((k - 1) / (m - 1)) + 1;

			l1 = (i - 2) * (m - 1) + j;
			l2 = l1 + m - 1;
			hh1 = alpha + i * h1;
			kk1 = beta + j * k1;

			//parameters functions in order to generalize to elliptic PDE
			if (i > 1) {
				//u_i-1_j
				a1[l1 - 1] = alphax(hh1, kk1) / h2
						- betax(hh1, kk1) / (2.0 * h1);
			}

			//u_i_j
			a1[l2 - 1] = gammaxy(hh1, kk1)
					- 2.0 * (alphax(hh1, kk1) / h2 + alphay(hh1, kk1) / k2);

			if (j > 1) {
				//u_i_j-1
				a1[l2 - 2] = alphay(hh1, kk1) / k2
						- betay(hh1, kk1) / (2.0 * k1);
			}

			if (j < m - 1) {
				//u_i_j+1
				a1[l2] = alphay(hh1, kk1) / k2 + betay(hh1, kk1) / (2.0 * k1);
			}

			l1 = l2 + m - 1; //update index to next row position i -> i+1

			if (i < n - 1) {
				//u_i+1_j
				a1[l1 - 1] = alphax(hh1, kk1) / h2
						+ betax(hh1, kk1) / (2.0 * h1);
			}

			s = bc->f(hh1, kk1);
			//cout << k << "S = F = [" << s << "]" << endl;

			if (i == 1) {
				s = s
						- (alphax(hh1, kk1) / h2 - betax(hh1, kk1) / (2.0 * h1))
								* bc->phi1(kk1);
				if (j == 1)
					s = s
							- (alphay(hh1, kk1) / k2
									- betay(hh1, kk1) / (2.0 * k1))
									* bc->phi2(hh1);
				if (j == m - 1)
					s = s
							- (alphay(hh1, kk1) / k2
									+ betay(hh1, kk1) / (2.0 * k1))
									* bc->phi4(hh1);
			} else if (i == n - 1) {
				s = s
						- (alphax(hh1, kk1) / h2 + betax(hh1, kk1) / (2.0 * h1))
								* bc->phi3(kk1);
				if (j == 1)
					s = s
							- (alphay(hh1, kk1) / k2
									- betay(hh1, kk1) / (2.0 * k1))
									* bc->phi2(hh1);
				if (j == m - 1)
					s = s
							- (alphay(hh1, kk1) / k2
									+ betay(hh1, kk1) / (2.0 * k1))
									* bc->phi4(hh1);
			} else {
				if (j == 1)
					s = s
							- (alphay(hh1, kk1) / k2
									- betay(hh1, kk1) / (2.0 * k1))
									* bc->phi2(hh1);
				if (j == m - 1)
					s = s
							- (alphay(hh1, kk1) / k2
									+ betay(hh1, kk1) / (2.0 * k1))
									* bc->phi4(hh1);
			}
			//cout << k << ">> S(BC=XY)= [" << s << "]" << endl;
			//cout << "--------------------------------------------" << endl
			//		<< endl;

//			cout << k << " S= [" << s << "]" << endl;
//			if (k == 81){
//				cout << k << " S= [" << s << "]" << endl;
//			}

			a1[n2 - 1] = s;
			for (int i = 1; i <= n1; i++) {
				rh = r[i - 1];
				if (rh != 0)
					b1[rh - 1] = a1[i - 1];
			}
			kh = k - 1;
			l = 0;
			max = 0;

			for (int j1 = 1; j1 <= n2; j1++) {
				if (r[j1 - 1] == 0) {
					s = a1[j1 - 1];
					l = l + 1;
					q = l;
					for (int i = 1; i <= kh; i++) {
						s = s - b1[i - 1] * x[q - 1];
						q = q + p;
					}
					a1[l - 1] = s;
					s = abs(s);
					if ((j1 < n2) && (s > max)) {
						max = s;
						jh = j1;
						lh = l;
					}
				}
			}

			if (max == 0)
				st = 5;
			else {
				max = 1 / a1[lh - 1];
				r[jh - 1] = k;
				for (int i = 1; i <= p; i++)
					a1[i - 1] = max * a1[i - 1];
				jh = 0;
				q = 0;
				for (j1 = 1; j1 <= kh; j1++) {
					s = x[q + lh - 1];
					for (int i = 1; i <= p; i++) {
						if (i != lh) {
							jh = jh + 1;
							x[jh - 1] = x[q + i - 1] - s * a1[i - 1];
						}
					}
					q = q + p;
				}
				for (int i = 1; i <= p; i++) {
					if (i != lh) {
						jh = jh + 1;
						x[jh - 1] = a1[i - 1];
					}
				}
				p = p - 1;
			}
			if (j == m - 1)
				j = 0;
		} while ((k != n1) && (st != 5));

		delete[] a1;
		delete[] b1;
		if (st == 0) {
			u = new long double*[n + 1];
			for (int i = 0; i <= n; i++)
				u[i] = new long double[m + 1];
			for (k = 1; k <= n1; k++) {
				rh = r[k - 1];
				if ((rh != k) && (rh != 0)) {
					s = x[k - 1];
					x[k - 1] = x[rh - 1];
					i = r[rh - 1];
					while (i != k) {
						x[rh - 1] = x[i - 1];
						r[rh - 1] = rh;
						rh = i;
						i = r[rh - 1];
					}
					x[rh - 1] = s;
					r[rh - 1] = rh;
				}
			}
			for (int i = 1; i <= n - 1; i++)
				for (int j = 1; j <= m - 1; j++)
					u[i][j] = x[(i - 1) * (m - 1) + j - 1];
			for (int i = 1; i <= n - 1; i++) {
				hh1 = alpha + i * h1;
				u[i][0] = bc->phi2(hh1);
				u[i][m] = bc->phi4(hh1);
			}
			for (int j = 0; j <= m; j++) {
				kk1 = beta + j * k1;
				u[0][j] = bc->phi1(kk1);
				u[n][j] = bc->phi3(kk1);
			}

			maxP = 0;
			maxQ = 0;
			maxR = 0;
			maxS = 0;
			if (this->_estimateMN) {

				for (int i = 3; i <= n - 3; i++)
					for (int j = 3; j <= m - 3; j++) {
//						tmpP = u[i + 2][j + 1] - 2 * u[i][j + 1]
//										+ u[i - 2][j + 1] - u[i + 2][j - 1]
//										+ 2 * u[i][j - 1] - u[i - 2][j - 1];
//						tmpP = (1.0 / (8 * h1 * h1 * k1)) * tmpP;

						tmpP = u[i + 1][j + 1] - u[i + 1][j - 1]
								- 2.0 * u[i][j + 1] + 2.0 * u[i][j - 1]
								+ u[i - 1][j + 1] - u[i - 1][j - 1];
						tmpP = (1.0 / (2.0 * h1 * h1 * k1)) * tmpP;
						if (abs(tmpP) > maxP)
							maxP = abs(tmpP);

//						tmpQ = u[i + 1][j + 2] - 2 * u[i + 1][j]
//								+ u[i + 1][j - 2] - u[i - 1][j + 2]
//								+ 2 * u[i - 1][j] - u[i - 1][j - 2];
//						tmpQ = (1.0 / (8 * h1 * k1 * k1)) * tmpQ;
						tmpQ = u[i + 1][j + 1] + u[i + 1][j - 1]
								- 2.0 * u[i + 1][j] + 2.0 * u[i - 1][j]
								- u[i - 1][j + 1] - u[i - 1][j - 1];
						tmpQ = (1.0 / (2.0 * h1 * k1 * k1)) * tmpQ;
						if (abs(tmpQ) > maxQ)
							maxQ = abs(tmpQ);

//						tmpR = u[i + 2][j + 2] - 2 * u[i][j + 2]
//								+ u[i - 2][j + 2] - 2 * u[i + 2][j]
//								+ 4 * u[i][j] - 2 * u[i - 2][j]
//								+ u[i + 2][j - 2] - 2 * u[i][j - 2]
//								+ u[i - 2][j - 2];
//						tmpR = (1.0 / (16 * h1 * h1 * k1 * k1)) * tmpR;
						tmpR = u[i + 2][j] - 2.0 * u[i + 1][j]
								+ 2.0 * u[i - 1][j] - u[i - 2][j];
//						cout << "tmpR = " << tmpR << endl;
						tmpS = (1.0 / (2.0 * h1 * h1 * h1));
//						cout << "1/(2*h^3) = " << tmpS << endl;
						tmpR = tmpS * tmpR;
//						cout << "[1/(2*h^3)] * tmpR = " << tmpR << endl;
						if (abs(tmpR) > maxR)
							maxR = abs(tmpR);

//						tmpS = -u[i - 2][j] + 2 * u[i - 1][j] - 2 * u[i + 1][j]
//								+ u[i + 2][j];
//						tmpS = (1.0 / (2.0 * h1 * h1 * h1)) * tmpS;
						tmpS = u[i + 1][j + 1] - 2 * u[i][j + 1]
								+ u[i - 1][j + 1] - 2 * u[i + 1][j]
								+ 4 * u[i][j] - 2 * u[i - 1][j]
								+ u[i + 1][j - 1] - 2 * u[i][j - 1]
								+ u[i - 1][j - 1];
//						cout << "tmpS = " << tmpS << endl;
						tmpR = (1.0 / (h1 * h1 * k1 * k1));
//						cout << "1/(h^2*k^2) = " << tmpR << endl;
						tmpS = tmpR * tmpS;
//						cout << "[1/(h^2*k^2)] * tmpS = " << tmpS << endl;
						if (abs(tmpS) > maxS)
							maxS = abs(tmpS);
					}

//				for (int i = 0; i <= n; i++)
//					delete[] u[i];
				delete[] u;
			}

		}
		delete[] x;
		delete[] r;
	}

	return 0;
}

template<typename T>
int GPESolver3C<T>::SolvePIA() {
	if (!_initparams)
		throw runtime_error("Parameters not initialized!");

	int st = 0;
	int n = params.n;
	int m = params.m;
	Interval<T> ALPHA = { params.alpha, params.alpha };
	Interval<T> BETA = { params.beta, params.beta };
	Interval<T> GAMMA = { params.gamma, params.gamma };
	Interval<T> DELTA = { params.delta, params.delta };
	long double eps = params.eps;
	long double sigma = 1e-3;

	const Interval<T> izero = { 0, 0 };
	const Interval<T> ione = { 1, 1 };
	const Interval<T> itwo = { 2, 2 };
	const Interval<T> ithree = { 3, 3 };
	const Interval<T> itwelve = { 12, 12 };
	const Interval<T> isigma = { -sigma, sigma };

	Interval<T> tmpi = { 0, 0 };
	int i, j, jh, j1, k, kh, l, lh, l1, l2, n1, n2, n3, p, q, rh;
	Interval<T> AF, BB0, BB1, CF, H1, HH, HH1, MHH, II, JJ, K1, KK, KK1, MKK,
			MAX, MM, AA, CC, Pconst, Qconst, Rconst, Sconst, NN, S, S1, S2, S3,
			S4, S5, ERR, H1POW2, K1POW2, H1POW2K1POW2, A1F, A2F, H2D12, K2D12;
	Interval<T> aij;
	int* r;
	T z;
	THashMap<T> bm;

	st = 0;
	NN.a = n;
	NN.b = n;
	MM.a = m;
	MM.b = m;
	Pconst.a = -bc->GetConstP();
	Pconst.b = bc->GetConstP();
	Qconst.a = -bc->GetConstQ();
	Qconst.b = bc->GetConstQ();
	Rconst.a = -bc->GetConstR();
	Rconst.b = bc->GetConstR();
	Sconst.a = -bc->GetConstT();
	Sconst.b = bc->GetConstT();

	if ((n < 2) || (m < 2))
		st = 1;
	else if ((ALPHA.a < 0) || (BETA.a < 0))
		st = 2;

	if (st == 0)

	{
		S = bc->PHI1(ALPHA, st);
		if (st == 0)

		{
			S1 = bc->PHI2(BETA, st);
			if (st == 0)

			{
				st = bc->boundconds(S, S1, eps);
				if (st == 0)

				{
					S = bc->PHI2(GAMMA, st);
					if (st == 0)

					{
						S1 = bc->PHI3(BETA, st);
						if (st == 0)

						{
							st = bc->boundconds(S, S1, eps);
							if (st == 0)

							{
								S = bc->PHI3(DELTA, st);
								if (st == 0)

								{
									S1 = bc->PHI4(GAMMA, st);
									if (st == 0)

										//<--- 20

										st = bc->boundconds(S, S1, eps);
									if (st == 0)

									{
										S = bc->PHI4(ALPHA, st);
										if (st == 0)

										{
											S1 = bc->PHI1(DELTA, st);
											if (st == 0)
												st = bc->boundconds(S, S1, eps);
										}
									}
								}
								//---> 20
							}
						}
					}
				}
			}
		}
	}

	if (st == 0) {
		H1 = (GAMMA - ALPHA) / NN;
		K1 = (DELTA - BETA) / MM;
		HH.a = -H1.b;
		HH.b = H1.b;
		H1POW2 = H1 * H1;
		K1POW2 = K1 * K1;
		H1POW2K1POW2 = H1POW2 * K1POW2;
		KK.a = -K1.b;
		KK.b = K1.b;
		n1 = (n - 1) * (m - 1);
		n2 = n1 + 1;
		p = n2;

		r = new int[n1 + 1];
		for (i = 1; i <= n2; i++)
			r[i - 1] = 0;
		k = 0;
		j = 0;
		n3 = (n * m - n - m + 4) * (n * m - n - m + 4) / 4;

		this->X = new Interval<T> [n3];
		for (i = 1; i <= n3; i++)
			this->X[i - 1] = izero;

		do {
			k = k + 1;
			j = j + 1;
			i = trunc((k - 1) / (m - 1)) + 1;
			l1 = (i - 2) * (m - 1) + j;
			l2 = l1 + m - 1;
			II.a = i;
			II.b = i;
			JJ.a = j;
			JJ.b = j;
			HH1 = ALPHA + (II * H1);
			KK1 = BETA + (JJ * K1);

			//parameters functions in order to generalize to elliptic PDE
			A1F = bc->A1(HH1, KK1, st);
			A2F = bc->A2(HH1, KK1, st);

			if (i > 1) {
				//u_i-1_j
				S1 = IAlphaX(HH1, KK1) / this->ih2
						- IBetaX(HH1, KK1) / (i2 * H1);
				bm.ToMap(l1 - 1, S1);
//				cout << k << "S1(l1-1)= [" << S1.a << " ; " << S1.b << "]"
//						<< endl;

			}

			//u_i_j
			S1 = IGammaXY(HH1, KK1)
					- i2
							* (this->IAlphaX(HH1, KK1) / this->ih2
									+ this->IAlphaY(HH1, KK1) / this->ik2);
			bm.ToMap(l2 - 1, S1);
//			cout << k << "S1(l2-1)= [" << S1.a << " ; " << S1.b << "]" << endl;

			if (j > 1) {
				//u_i_j-1
				S1 = this->IAlphaY(HH1, KK1) / this->ik2
						- this->IBetaY(HH1, KK1) / (i2 * K1);
				bm.ToMap(l2 - 2, S1);
//				cout << k << "S1(l2-2)= [" << S1.a << " ; " << S1.b << "]"
//						<< endl;
			}

			if (j < m - 1) {
				//u_i_j+1
				S1 = IAlphaY(HH1, KK1) / this->ik2
						+ this->IBetaY(HH1, KK1) / (i2 * K1);
				bm.ToMap(l2, S1);
//				cout << k << "S1(l2)= [" << S1.a << " ; " << S1.b << "]"
//						<< endl;
			}

			l1 = l2 + m - 1; //update index to next row position i -> i+1

			if (i < n - 1) {
				//u_i+1_j
				S1 = IAlphaX(HH1, KK1) / this->ih2
						+ IBetaX(HH1, KK1) / (i2 * H1);
				bm.ToMap(l1 - 1, S1);
//				cout << k << "S1(l1-1)= [" << S1.a << " ; " << S1.b << "]"
//						<< endl;
			}

			S = bc->F(HH1, KK1, st);
//			cout << k << "S = F = [" << S.a << " ; " << S.b << "]" << endl;
			ERR = (this->ih2 / i12)
					* (bc->D2FDX2(HH1 + HH, KK1, st)
							+ im2 * bc->DA1DX(HH1, KK1, st)
									* bc->DFDX(HH1 + HH, KK1, st)
									/ bc->A1(HH1, KK1, st)
							+ im2
									* (bc->DA2DX(HH1, KK1, st)
											- bc->A2(HH1, KK1, st)
													* bc->DA1DX(HH1, KK1, st)
													/ bc->A1(HH1, KK1, st))
									* Pconst - bc->A2(HH1, KK1, st) * Sconst);

//			cout << "(i,j) = (" << i << ", " << j << ")" << endl;
//			cout << k << "ADD DX2 B: S= [" << S.a << " ; " << S.b << "]" << endl;
//			cout << k << "ADD1 DX2: ERR= [" << ERR.a << " ; " << ERR.b << "]"
//					<< endl;

			S = S + ERR;

			ERR = (this->ik2 / i12)
					* (bc->D2FDY2(HH1, KK1 + KK, st)
							+ im2 * bc->DA2DY(HH1, KK1, st)
									* bc->DFDY(HH1, KK1 + KK, st)
									/ bc->A2(HH1, KK1, st)
							+ im2
									* (bc->DA1DY(HH1, KK1, st)
											- bc->A1(HH1, KK1, st)
													* bc->DA2DY(HH1, KK1, st)
													/ bc->A2(HH1, KK1, st))
									* Qconst - bc->A1(HH1, KK1, st) * Sconst);

//			cout << "(i,j) = (" << i << ", " << j << ")" << endl;
//			cout << k << "ADD2 DY2: S= [" << S.a << " ; " << S.b << "]" << endl;
//			cout << k << "ADD2 DY2: ERR= [" << ERR.a << " ; " << ERR.b << "]"
//					<< endl;

			S = S + ERR;
//			cout << k << "S [+ERR2] = [" << S.a << " ; " << S.b << "]" << endl;
			S = S + isigma;
//			cout << k << "S(not BC)= [" << S.a << " ; " << S.b << "]" << endl;

			if (i == 1) {
				S = S
						- (IAlphaX(HH1, KK1) / this->ih2
								- IBetaX(HH1, KK1) / (i2 * H1))
								* bc->PHI1(KK1, st);

			} else if (i == n - 1) {
				S = S
						- (IAlphaX(HH1, KK1) / this->ih2
								+ IBetaX(HH1, KK1) / (i2 * H1))
								* bc->PHI3(KK1, st);
			}
			if (j == 1)
				S = S
						- (IAlphaY(HH1, KK1) / this->ik2
								- IBetaY(HH1, KK1) / (i2 * K1))
								* bc->PHI2(HH1, st);
			if (j == m - 1)
				S = S
						- (IAlphaY(HH1, KK1) / this->ik2
								+ IBetaY(HH1, KK1) / (i2 * K1))
								* bc->PHI4(HH1, st);

//			cout << k << ">> S(BC=XY)= [" << S.a << " ; " << S.b << "]" << endl;
//			cout << "--------------------------------------------" << endl
//					<< endl;
//			filestr << k << " E: S= [" << S.a << " ; " << S.b << "]" << endl;
			if (S.a > S.b)
				runtime_error("Wrong F value!");
			if (st == 0) {
				bm.ToMap(n2 - 1, S);

				for (int i = 1; i <= n1; i++) {
					rh = r[i - 1];
					if ((k > 1) && (rh == k - 1)) {
						BB1 = bm.FromMap(i - 1);
					}
					if ((k > m - 1) && (rh == k - m + 1)) {
						BB0 = bm.FromMap(i - 1);
					}
				}
				kh = k - 1;
				l = 0;
				MAX.a = 0;
				MAX.b = 0;
				for (j1 = 1; j1 <= n2; j1++) {
					if (r[j1 - 1] == 0) {
						S = bm.FromMap(j1 - 1);
						l = l + 1;
						q = l;
						for (int i = 1; i <= kh; i++) {
							if ((k > 1) && (i == k - 1))
								S = S - (BB1 * this->X[q - 1]);
							if ((k > m - 1) && (i - 1 == k - m))
								S = S - (BB0 * this->X[q - 1]);
							q = q + p;
						}
						if (!((S.a == 0) && (S.b == 0)))
							bm.ToMap(l - 1, S);
						else
							bm.Erase(l - 1);

						if (S.a < 0)
							S.a = abs(S.a);
						if (S.b < 0)
							S.b = abs(S.b);

						if (S.b < S.a) {
							z = S.a;
							S.a = S.b;
							S.b = z;
						}
						if ((j1 < n2) && (S.b > MAX.a)) {
							MAX = S;
							jh = j1;
							lh = l;
						}
					}
				}
				//cout << "MAX= [" << MAX.a << " ; " << MAX.b << "]" << endl;
				if ((MAX.a == 0) && (MAX.b == 0))
					st = 5;
				else {
					tmpi = bm.FromMap(lh - 1);
					MAX = (ione / tmpi);
					r[jh - 1] = k;
					for (int i = 1; i <= p; i++) {
						tmpi = bm.FromMap(i - 1);
						S = (MAX * tmpi);
						if (!((S.a == 0) && (S.b == 0)))
							bm.ToMap(i - 1, S);
						else
							bm.Erase(i - 1);
					}
					jh = 0;
					q = 0;
					for (j1 = 1; j1 <= kh; j1++) {
						S = this->X[q + lh - 1];
						for (int i = 1; i <= p; i++)
							if (i != lh) {
								jh = jh + 1;
								S1 = bm.FromMap(i - 1);
								this->X[jh - 1] = this->X[q + i - 1] - (S * S1);
							}
						q = q + p;
					}
					for (int i = 1; i <= p; i++) {
						if (i != lh) {
							jh = jh + 1;
							tmpi = bm.FromMap(i - 1);
							this->X[jh - 1] = tmpi;
							//filestr << jh-1 << ": X= [" << this->X[jh - 1].a << " ; " << this->X[jh - 1].b << "]" << endl;
						}
					}
					p = p - 1;
				}

				if (j == m - 1) {
					j = 0;
				}
			}

			bm.Clear();

		} while (!((k == n1) || (st == 5)));

		if (st == 0) {
			for (int k = 1; k <= n1; k++) {
				rh = r[k - 1];
				if (rh != k) {
					S = this->X[k - 1];
					this->X[k - 1] = this->X[rh - 1];
					i = r[rh - 1];
					while (i != k) {
						this->X[rh - 1] = this->X[i - 1];
						r[rh - 1] = rh;
						rh = i;
						i = r[rh - 1];
					}
					this->X[rh - 1] = S;
					r[rh - 1] = rh;
				}
			}
		}
	}
//	filestr.close();
	return 0;
}

template<typename T>
int GPESolver3C<T>::SolveDIA() {
	if (!_initparams)
			throw runtime_error("Parameters not initialized!");

		int st = 0;
		int n = params.n;
		int m = params.m;
		Interval<T> ALPHA = { params.alpha, params.alpha };
		Interval<T> BETA = { params.beta, params.beta };
		Interval<T> GAMMA = { params.gamma, params.gamma };
		Interval<T> DELTA = { params.delta, params.delta };
		long double eps = params.eps;
		long double sigma = 0.0; //1e-3;

		const Interval<T> izero = { 0, 0 };
		const Interval<T> ione = { 1, 1 };
		const Interval<T> itwo = { 2, 2 };
		const Interval<T> ithree = { 3, 3 };
		const Interval<T> itwelve = { 12, 12 };
		const Interval<T> isigma = { -sigma, sigma };

		Interval<T> tmpi = { 0, 0 };
		int i, j, jh, j1, k, kh, l, lh, l1, l2, n1, n2, n3, p, q, rh;
		Interval<T> AF, BB0, BB1, CF, H1, HH, HH1, MHH, II, JJ, K1, KK, KK1, MKK,
				MAX, MM, AA, CC, Pconst, Qconst, Rconst, Sconst, NN, S, S1, S2, S3,
				S4, S5, ERR, H1POW2, K1POW2, H1POW2K1POW2, A1F, A2F, H2D12, K2D12;
		Interval<T> aij;
		int* r;
		T z;
		THashMap<T> bm;

		st = 0;
		NN.a = n;
		NN.b = n;
		MM.a = m;
		MM.b = m;
		Pconst.a = bc->GetConstP();
		Pconst.b = -bc->GetConstP();
		Qconst.a = bc->GetConstQ();
		Qconst.b = -bc->GetConstQ();
		Rconst.a = bc->GetConstR();
		Rconst.b = -bc->GetConstR();
		Sconst.a = bc->GetConstT();
		Sconst.b = -bc->GetConstT();

		if ((n < 2) || (m < 2))
			st = 1;
		else if ((ALPHA.a < 0) || (BETA.a < 0))
			st = 2;

		if (st == 0)

		{
			S = bc->PHI1(ALPHA, st);
			if (st == 0)

			{
				S1 = bc->PHI2(BETA, st);
				if (st == 0)

				{
					st = bc->boundconds(S, S1, eps);
					if (st == 0)

					{
						S = bc->PHI2(GAMMA, st);
						if (st == 0)

						{
							S1 = bc->PHI3(BETA, st);
							if (st == 0)

							{
								st = bc->boundconds(S, S1, eps);
								if (st == 0)

								{
									S = bc->PHI3(DELTA, st);
									if (st == 0)

									{
										S1 = bc->PHI4(GAMMA, st);
										if (st == 0)

											//<--- 20

											st = bc->boundconds(S, S1, eps);
										if (st == 0)

										{
											S = bc->PHI4(ALPHA, st);
											if (st == 0)

											{
												S1 = bc->PHI1(DELTA, st);
												if (st == 0)
													st = bc->boundconds(S, S1, eps);
											}
										}
									}
									//---> 20
								}
							}
						}
					}
				}
			}
		}

		if (st == 0) {
			H1 = (GAMMA - ALPHA) / NN;
			K1 = (DELTA - BETA) / MM;
			HH.a = -H1.b;
			HH.b = H1.b;
			H1POW2 = H1 * H1;
			K1POW2 = K1 * K1;
			H1POW2K1POW2 = H1POW2 * K1POW2;
			KK.a = -K1.b;
			KK.b = K1.b;
			n1 = (n - 1) * (m - 1);
			n2 = n1 + 1;
			p = n2;

			r = new int[n1 + 1];
			for (i = 1; i <= n2; i++)
				r[i - 1] = 0;
			k = 0;
			j = 0;
			n3 = (n * m - n - m + 4) * (n * m - n - m + 4) / 4;

			this->X = new Interval<T> [n3];
			for (i = 1; i <= n3; i++)
				this->X[i - 1] = izero;

			do {
				k = k + 1;
				j = j + 1;
				i = trunc((k - 1) / (m - 1)) + 1;
				l1 = (i - 2) * (m - 1) + j;
				l2 = l1 + m - 1;
				II.a = i;
				II.b = i;
				JJ.a = j;
				JJ.b = j;
				HH1 = ALPHA + (II * H1);
				KK1 = BETA + (JJ * K1);

				//parameters functions in order to generalize to elliptic PDE
				A1F = bc->A1(HH1, KK1, st);
				A2F = bc->A2(HH1, KK1, st);

				if (i > 1) {
					//u_i-1_j
					S1 = IAlphaX(HH1, KK1) / this->ih2
							- IBetaX(HH1, KK1) / (i2 * H1);
					bm.ToMap(l1 - 1, S1);
//					cout << k << "S1(l1-1)= [" << S1.a << " ; " << S1.b << "]"
//							<< endl;

				}

				//u_i_j
				S1 = IGammaXY(HH1, KK1)
						- i2
								* (this->IAlphaX(HH1, KK1) / this->ih2
										+ this->IAlphaY(HH1, KK1) / this->ik2);
				bm.ToMap(l2 - 1, S1);
//				cout << k << "S1(l2-1)= [" << S1.a << " ; " << S1.b << "]" << endl;

				if (j > 1) {
					//u_i_j-1
					S1 = this->IAlphaY(HH1, KK1) / this->ik2
							- this->IBetaY(HH1, KK1) / (i2 * K1);
					bm.ToMap(l2 - 2, S1);
//					cout << k << "S1(l2-2)= [" << S1.a << " ; " << S1.b << "]"
//							<< endl;
				}

				if (j < m - 1) {
					//u_i_j+1
					S1 = IAlphaY(HH1, KK1) / this->ik2
							+ this->IBetaY(HH1, KK1) / (i2 * K1);
					bm.ToMap(l2, S1);
//					cout << k << "S1(l2)= [" << S1.a << " ; " << S1.b << "]"
//							<< endl;
				}

				l1 = l2 + m - 1; //update index to next row position i -> i+1

				if (i < n - 1) {
					//u_i+1_j
					S1 = IAlphaX(HH1, KK1) / this->ih2
							+ IBetaX(HH1, KK1) / (i2 * H1);
					bm.ToMap(l1 - 1, S1);
//					cout << k << "S1(l1-1)= [" << S1.a << " ; " << S1.b << "]"
//							<< endl;
				}

				S = bc->F(HH1, KK1, st);
//				cout << k << "S = F = [" << S.a << " ; " << S.b << "]" << endl;
				ERR = (this->ih2 / i12)
						* (bc->D2FDX2(HH1 + HH, KK1, st)
								+ im2 * bc->DA1DX(HH1, KK1, st)
										* bc->DFDX(HH1 + HH, KK1, st)
										/ bc->A1(HH1, KK1, st)
								+ i2
										* (bc->DA2DX(HH1, KK1, st)
												- bc->A2(HH1, KK1, st)
														* bc->DA1DX(HH1, KK1, st)
														/ bc->A1(HH1, KK1, st))
										* Pconst - bc->A2(HH1, KK1, st) * Sconst);

				S = S + ERR;

				ERR = (this->ik2 / i12)
						* (bc->D2FDY2(HH1, KK1 + KK, st)
								+ im2 * bc->DA2DY(HH1, KK1, st)
										* bc->DFDY(HH1, KK1 + KK, st)
										/ bc->A2(HH1, KK1, st)
								+ i2
										* (bc->DA1DY(HH1, KK1, st)
												- bc->A1(HH1, KK1, st)
														* bc->DA2DY(HH1, KK1, st)
														/ bc->A2(HH1, KK1, st))
										* Qconst - bc->A1(HH1, KK1, st) * Sconst);
				S = S + ERR;
//				cout << k << "S(not BC)= [" << S.a << " ; " << S.b << "]" << endl;

				if (i == 1) {
					S = S
							- (IAlphaX(HH1, KK1) / this->ih2
									- IBetaX(HH1, KK1) / (i2 * H1))
									* bc->PHI1(KK1, st);

				} else if (i == n - 1) {
					S = S
							- (IAlphaX(HH1, KK1) / this->ih2
									+ IBetaX(HH1, KK1) / (i2 * H1))
									* bc->PHI3(KK1, st);
				}
				if (j == 1)
					S = S
							- (IAlphaY(HH1, KK1) / this->ik2
									- IBetaY(HH1, KK1) / (i2 * K1))
									* bc->PHI2(HH1, st);
				if (j == m - 1)
					S = S
							- (IAlphaY(HH1, KK1) / this->ik2
									+ IBetaY(HH1, KK1) / (i2 * K1))
									* bc->PHI4(HH1, st);

//				cout << k << ">> S(BC=XY)= [" << S.a << " ; " << S.b << "]" << endl;
//				cout << "--------------------------------------------" << endl
//						<< endl;
	//			filestr << k << " E: S= [" << S.a << " ; " << S.b << "]" << endl;
				if (st == 0) {
					bm.ToMap(n2 - 1, S);

					for (int i = 1; i <= n1; i++) {
						rh = r[i - 1];
						if ((k > 1) && (rh == k - 1)) {
							BB1 = bm.FromMap(i - 1);
						}
						if ((k > m - 1) && (rh == k - m + 1)) {
							BB0 = bm.FromMap(i - 1);
						}
					}
					kh = k - 1;
					l = 0;
					MAX.a = 0;
					MAX.b = 0;
					for (j1 = 1; j1 <= n2; j1++) {
						if (r[j1 - 1] == 0) {
							S = bm.FromMap(j1 - 1);
							l = l + 1;
							q = l;
							for (int i = 1; i <= kh; i++) {
								if ((k > 1) && (i == k - 1))
									S = S  - (BB1 * this->X[q - 1]);
								if ((k > m - 1) && (i - 1 == k - m))
									S = S - (BB0 * this->X[q - 1]);
								q = q + p;
							}
							if (!((S.a == 0) && (S.b == 0)))
								bm.ToMap(l - 1, S);
							else
								bm.Erase(l - 1);

							if (S.a < 0)
								S.a = abs(S.a);
							if (S.b < 0)
								S.b = abs(S.b);

							if (S.b < S.a) {
								z = S.a;
								S.a = S.b;
								S.b = z;
							}
							if ((j1 < n2) && (S.b > MAX.a)) {
								MAX = S;
								jh = j1;
								lh = l;
							}

						}
					}
					//cout << "MAX= [" << MAX.a << " ; " << MAX.b << "]" << endl;
					if ((MAX.a == 0) && (MAX.b == 0))
						st = 5;
					else {
						tmpi = bm.FromMap(lh - 1);
						MAX = tmpi.Inverse();//(ione / tmpi);
						r[jh - 1] = k;
						for (int i = 1; i <= p; i++) {
							tmpi = bm.FromMap(i - 1);
							S = (MAX * tmpi);
							if (!((S.a == 0) && (S.b == 0)))
								bm.ToMap(i - 1, S);
							else
								bm.Erase(i - 1);
						}
						jh = 0;
						q = 0;
						for (j1 = 1; j1 <= kh; j1++) {
							S = this->X[q + lh - 1];
							for (int i = 1; i <= p; i++)
								if (i != lh) {
									jh = jh + 1;
									S1 = bm.FromMap(i - 1);
//									tmpi = S1.Opposite();
//									S1 = Hull(S1, S1.);
									this->X[jh - 1] = this->X[q + i - 1] - (S * S1);
								}
							q = q + p;
						}
						for (int i = 1; i <= p; i++) {
							if (i != lh) {
								jh = jh + 1;
								tmpi = bm.FromMap(i - 1);
								this->X[jh - 1] = tmpi;
								//filestr << jh-1 << ": X= [" << this->X[jh - 1].a << " ; " << this->X[jh - 1].b << "]" << endl;
							}
						}
						p = p - 1;
					}

					if (j == m - 1) {
						j = 0;
					}
				}

				bm.Clear();

			} while (!((k == n1) || (st == 5)));

			if (st == 0) {
				for (int k = 1; k <= n1; k++) {
					rh = r[k - 1];
					if (rh != k) {
						S = this->X[k - 1];
						this->X[k - 1] = this->X[rh - 1];
						i = r[rh - 1];
						while (i != k) {
							this->X[rh - 1] = this->X[i - 1];
							r[rh - 1] = rh;
							rh = i;
							i = r[rh - 1];
						}
						this->X[rh - 1] = S;
						r[rh - 1] = rh;
					}
				}
			}
		}
	//	filestr.close();
		return 0;
}

//The explicit instantiation part
template class GPESolver3C<long double> ;
//template class GPESolver3C<mpreal> ;

}
/* namespace intervalarth */
