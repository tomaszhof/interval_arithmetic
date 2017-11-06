/*
 * PoissonSolver4OrderAM.cpp
 *
 *  Created on: 4 wrz 2017
 *      Author: tomhof
 */

#include "PoissonSolver4OrderAM.h"

namespace interval_arithmetic {


template<typename T>
PoissonSolver4OrderAM<T>::PoissonSolver4OrderAM() {
	// TODO Auto-generated constructor stub

}

template<typename T>
PoissonSolver4OrderAM<T>::~PoissonSolver4OrderAM() {
	// TODO Auto-generated destructor stub
}

template<typename T>
int PoissonSolver4OrderAM<T>::SetExample(int eid) {
	switch (eid) {
	case 1:
		//bc = new Example01();
		break;
	case 2:
		bc = new Example02<T>();
		break;
	case 3:
		bc = new Example03<T>();
		break;
	case 4:
		bc = new Example04<T>();
		break;
	case 5:
		bc = new Example05<T>();
		break;
	case 6:
		bc = new Example06<T>();
		break;
	case 7:
		bc = new Example07<T>();
		break;
	case 8:
		bc = new Example08<T>();
		break;
	case 9:
		bc = new Example09<T>();
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
int PoissonSolver4OrderAM<T>::SolveFP() {
	int i, j, jh, j1, k, kh, l, lh, l1, l2, n1, n2, p, q, rh, st;
	long double af, cf, h1, k1, h2, k2, h4, k4,hh1, kk1, max, s, tmpM, tmpN, h2d12, k2d12,
			hhp1, hhm1, kkp1, kkm1;

	//coefficients
	long double ca, cb, cc, s1, s2;

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
	else if ((alpha < 0) || (beta < 0))
		st = 2;
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
		h4 = h2 * h2;
		k4 = k2 * k2;
		ca = h2 + k2;
		cb = 2*(5*h2-k2);
		cc = 2*(5*k2-h2);
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
//
//		ofstream rowFile;
//		rowFile.open("row.txt");

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
			hhp1 = alpha + (i + 1) * h1;
			kkp1 = beta + (j + 1) * k1;
			hhm1 = alpha + (i - 1) * h1;
			kkm1 = beta + (j - 1) * k1;

			//changed - in order to generalize to elliptic PDE

			af = 1.0 / h1; //bc->a(hh1, kk1) / h1; //1 / h1;
			cf = 1.0 / k1; //bc->c(hh1, kk1) / k1; //1 / k1;

			//fourth order method - option 6 (Zhang)
			if (i > 1) {
				a1[l1 - 1] = cc; //af / h1;

				if (j > 1) {
					a1[l1 - 2] = ca;
				}
				if (j < m - 1) {
					a1[l1] = ca;
				}

			}

			a1[l2 - 1] = (-20.0) * ca; //-2 * (af / h1 + cf / k1);

			if (j > 1) {
				a1[l2 - 2] = cb; //cf / k1;
			}
			if (j < m - 1)
				a1[l2] = cb; //cf / k1;

			l1 = l2 + m - 1;
			if (i < n - 1) {
				a1[l1 - 1] = cc;
				if (j > 1) {
					a1[l1 - 2] = ca;
				}
				if (j < m - 1) {
					a1[l1] = ca;
				}
			}

//		    Fortuna Macukow option 5
			s = 8.0 * bc->f(hh1, kk1);
			s = s
					+ (bc->f(hhp1, kk1) + bc->f(hhm1, kk1)
									+ bc->f(hh1, kkp1) + bc->f(hh1, kkm1));

//			cout << "k = " << k << "; s = " << s << endl;

			//Zhang- option 6
			if (i == 1) {
				s1 = bc->phi1(kkm1);
				s2 = bc->phi1(kkp1);
				s1 = ca * (s1 + s2);
				s2 = bc->phi1(kk1);
				s1 = s1 + ca*s2;
				s = s -s1;
//				if (j < m)
//				{
//					s = s - (af / h1) * bc->phi1(kkp1);
//				}
//				if (j > 0)
//					s = s - (af / h1) * bc->phi1(kkm1);
				if (j == 1) {
					s1 = bc->phi2(hhp1);
					s = s-ca*s1;
					s1 = bc->phi2(hh1);
					s = s-cb*s1;
				}
				if (j == m - 1) {
				    s1 = bc->phi4(hhp1);
				    s =s-ca*s1;
				    s1 =bc->phi4(hh1);
				    s =s-cb*s1;
				}
			} else if (i == n - 1) {
				s1 = bc->phi3(kk1);
				s2 = bc->phi3(kkp1);
				s1 = ca*(s1+s2);
				s2 = bc->phi3(kk1);
				s1 = s1+cc*s2;
				s = s - s1;
//				if (j < m)
//					s = s - (af / h1) * bc->phi3(kkp1);
//				if (j > 0)
//					s = s - (af / h1) * bc->phi3(kkm1);
				if (j == 1) {
					s1 = bc-> phi2(hhm1);
					s = s-ca*s1;
					s1 = bc->phi2(hh1);
					s = s-cb*s1;
				}
				if (j == m - 1) {
				  s1 = bc-> phi4(hhm1);
				  s=s-ca*s1;
				  s1=bc->phi4(hh1);
				  s=s-cb*s1;
				}
			} else {
				if (j == 1) {
					s1 =bc->phi2(hhm1);
					  s2=bc->phi2(hhp1);
					  s1=s1+s2;
					  s=s-ca*s1;
					  s1=bc->phi2(hh1);
					  s=s-cb*s1;
				}
				if (j == m - 1) {
				  s1 = bc->phi4(hhm1);
				  s2 = bc->phi4(hhp1);
				  s1 = s1+s2;
				  s = s-ca*s1;
				  s1 = bc->phi4(hh1);
				  s = s-cb*s1;
				}
			}

			//if (i < 10) cout << "s = " << s << endl;

			a1[n2 - 1] = s;

//			for (int i = 1; i <= n1; i++) {
//				rowFile << a1[i - 1] << ";";
//			}
//			rowFile << s << "\n";

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
						//cout << jh-1 << ": x= " << x[jh-1] << endl;
					}
				}
				p = p - 1;
			}
			if (j == m - 1)
				j = 0;
		} while ((k != n1) && (st != 5));

//		rowFile.close();

		delete[] a1;
		delete[] b1;
		if (st == 0) {
			u = new long double*[n + 1];
			for (int i = 0; i <= n; i++)
				u[i] = new long double[m + 1];
			for (k = 1; k <= n1; k++) {
				rh = r[k - 1];
				if (rh != k) {
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

			maxM = 0;
			maxN = 0;
			if (this->_estimateMN) {

				for (int i = 3; i <= n - 3; i++)
					for (int j = 3; j <= m - 3; j++) {
						tmpM = abs(
								-12 * u[i][j] + 6 * (u[i][j-1] + u[i][j+1])
										+ 8 * (u[i-1][j] + u[i+1][j])
										- 4 * (u[i-1][j-1] + u[i-1][j+1]
										+ u[i+1][j-1] + u[i+1][j+1])
										- 2 * (u[i-2][j] + u[i+2][j])
										+ u[i-2][j-1] + u[i-2][j+1]
										+ u[i+2][j-1] + u[i+2][j+1]);
						tmpM = (1 / (h1 * h1 * h1 * h1 * k1 * k1)) * tmpM;
						if (tmpM > maxM)
							maxM = tmpM;
						tmpN = abs(
								-12 * u[i][j] + 6 * (u[i-1][j] + u[i+1][j])
										+ 8 * (u[i][j - 1] + u[i][j + 1])
										- 4 * (u[i-1][j-1] + u[i+1][j-1]
										+ u[i-1][j+1] + u[i+1][j+1])
										- 2 * (u[i][j-2] + u[i][j+2])
										+ u[i-1][j-2] + u[i+1][j-2]
										+ u[i-1][j+2] + u[i+1][j+2]);
						tmpN = (1 / (h1 * h1 * h1 * h1 * k1 * k1)) * tmpN;
						if (tmpN > maxN)
							maxN = tmpN;
					}

				for (int i = 0; i <= n; i++)
					delete[] u[i];
				delete[] u;
			}

		}
		delete[] x;
		delete[] r;
	}

	return 0;
}

template<typename T>
int PoissonSolver4OrderAM<T>::SolvePIA() {
//	fstream filestr;
//	string fname = "tmpLog.txt";
//	filestr.open(fname.c_str(), fstream::out);

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

	const Interval<T> izero = { 0, 0 };
	const Interval<T> ione = { 1, 1 };
	const Interval<T> itwo = { 2, 2 };
	const Interval<T> ithree = { 3, 3 };
	const Interval<T> ifour = { 4, 4 };
	const Interval<T> ifive = { 5, 5 };
	const Interval<T> itwelve = { 12, 12 };
	const Interval<T> itwenty = { 20, 20 };
	const Interval<T> imtwenty = { -20, -20 };
	const Interval<T> ithreesixty = { 360, 360 };
	const Interval<T> ininety = { 90, 90 };

	Interval<T> tmpi = { 0, 0 };
	int i, j, jh, j1, k, kh, l, lh, l1, l2, n1, n2, n3, p, q, rh;
	int num;
	Interval<T> AF, BB0, BB1, CF, H1, HH, HH1, II, JJ, K1, KK, KK1, MAX, MM, AA,
			CC, PPconst, QQconst, NN, S, S1, S2, S3, S4, S5, H1POW2, K1POW2,
			IIP1, IIM1, JJP1, JJM1, KKP1, KKM1, HHP1, HHM1, H1POW2K1POW2;
	Interval<T> H2, H4, K2, K4, PQ,CA, CB;
	bool list_exists;
	Interval<T> aij;
	int* r;
	T z;
	THashMap<T> bm;
	THashMap<T>* B1M = new THashMap<T>();

	st = 0;
	NN.a = n;
	NN.b = n;
	MM.a = m;
	MM.b = m;
	PPconst.a = -bc->GetConstM();
	PPconst.b = bc->GetConstM();
	QQconst.a = -bc->GetConstN();
	QQconst.b = bc->GetConstN();

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
		H2 = H1 * H1;
		H4 = H2 * H2;
		K2 = K1 * K1;
		K4 = K2 * K2;
		HH.a = -H1.b;
		HH.b = H1.a;
		KK.a = -K1.b;
		KK.b = K1.a;
		H1POW2 = H1 * H1;
		K1POW2 = K1 * K1;
		H1POW2K1POW2 = H1POW2 * K1POW2;

        CA=H2+K2;     // interval h^2+k^2
        CB=(ifive*H2)- K2;
        CB =itwo*CB;   // interval 2*(5*h^2-k^2)
        CC=(ifive*K2)-H2;
        CC=itwo*CC;   // interval 2*(5*k^2-h^2);
        PQ=H4*PPconst;
        PQ=PQ+K4*QQconst;
        PQ=PQ/itwenty;
        PPconst= PPconst+ QQconst;
        PPconst=(H2*K2)*PPconst;
        PPconst=PPconst/itwelve;
        PQ=PQ+PPconst;
		 // PQ = (h^4*[-P,P]+k^4*[-Q,Q])/20
		 //       +h^2*k^2*([-P,P]+{-Q,Q])/12
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
		num = 1;

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
			IIP1.a = i + 1;
			IIP1.b = i + 1;
			JJP1.a = j + 1;
			JJP1.b = j + 1;
			IIM1.a = i - 1;
			IIM1.b = i - 1;
			JJM1.a = j - 1;
			JJM1.b = j - 1;

			HH1 = ALPHA + (II * H1);
			KK1 = BETA + (JJ * K1);

			HHM1 = ALPHA + (IIM1 * H1);
			KKM1 = BETA + (JJM1 * K1);

			HHP1 = ALPHA + (IIP1 * H1);
			KKP1 = BETA + (JJP1 * K1);
			AA = bc->A(HH1, KK1, st);
			CC = bc->C(HH1, KK1, st);
			AF = AA / H1;
			CF = CC / K1;
			S1 = imtwenty * AF / H1;

			//fourth order method - option 6 (Zhang)
			if (i > 1) {
				bm.ToMap(l1 - 1, ifour * AF / H1);

				if (j > 1) {
					bm.ToMap(l1 - 2, AF / H1);
				}
				if (j < m - 1) {
					bm.ToMap(l1, AF / H1);
				}

			}

			bm.ToMap(l2 - 1, S1);

			if (j > 1) {
				bm.ToMap(l2 - 2, ifour * (AF / H1));
			}
			if (j < m - 1)
				bm.ToMap(l2, ifour * (AF / H1));

			l1 = l2 + m - 1;

			if (i < n - 1) {
				bm.ToMap(l1 - 1, ifour * (AF / H1));
				if (j > 1) {
					bm.ToMap(l1 - 2, (AF / H1));
				}
				if (j < m - 1) {
					bm.ToMap(l1, (AF / H1));
				}
			}

			//right site

			//Zhang- option 6
			S = ifour * bc->F(HH1, KK1, st);
			S = S
					+ (ione / itwo)
							* (bc->F(HHP1, KK1, st) + bc->F(HHM1, KK1, st)
									+ bc->F(HH1, KKP1, st)
									+ bc->F(HH1, KKM1, st));
			cout << "k = " << k << "; S= [" << S.a << " ; " << S.b << "]" << endl;
//			filestr << k << " B: S= [" << S.a << " ; " << S.b << "]" << endl;

			if (st == 0) {
				//second order truncation error
//				S5 = ((AA * H1POW2) * H1POW2K1POW2) * MMconst;
//				S3 = ((CC * K1POW2) * H1POW2K1POW2) * NNconst;
//				S1 = S3 + S5;
//				S = S + (S1 / itwelve);

				//fourth order truncation error
				S5 = H1POW2 * H1POW2  / ithreesixty;
				S3 = (bc->PSI1(HH1, KK1, st) + bc->PSI2(HH1, KK1, st));
				S1 = S5 * S3;
				S = S + S1;
				S5 = H1POW2 * H1POW2 / ininety;
				S3 = bc->PSI3(HH1, KK1, st);
				S1 = S5 * S3;
				S = S + S1;
				S = S + S5 * (PPconst + QQconst);

				if (i == 1) {
					S1 = ifour * (AF / H1) * bc->PHI1(KK1, st);
					S = S - S1;

					if (j < m) {
						S1 = (AF / H1) * bc->PHI1(KKP1, st);
						S = S - S1;
					}
					if (j > 0) {
						S1 = (AF / H1) * bc->PHI1(KKM1, st);
						S = S - S1;
					}

					if (j == 1) {
						S1 = ifour * (AF / H1) * bc->PHI2(HH1, st);
						S = S - S1;

						S1 = (AF / H1) * bc->PHI2(HHP1, st);
						S = S - S1;
					}

					if (j == m - 1) {
						S1 = ifour * (AF / H1) * bc->PHI4(HH1, st);
						S = S - S1;

						S1 = (AF / H1) * bc->PHI4(HHP1, st);
						S = S - S1;
					}

				} else if (i == n - 1) {
					S1 = ifour * (AF / H1) * bc->PHI3(KK1, st);
					S = S - S1;
					if (j < m) {
						S1 = (AF / H1) * bc->PHI3(KKP1, st);
						S = S - S1;
					}
					if (j > 0) {
						S1 = (AF / H1) * bc->PHI3(KKM1, st);
						S = S - S1;
					}

					if (j == 1) {
						S1 = ifour * (AF / H1) * bc->PHI2(HH1, st);
						S = S - S1;
						S1 = (AF / H1) * bc->PHI2(HHM1, st);
						S = S - S1;
					}

					if (j == m - 1) {
						S1 = ifour * (AF / H1) * bc->PHI4(HH1, st);
						S = S - S1;
						S1 = (AF / H1) * bc->PHI4(HHM1, st);
						S = S - S1;
					}

				} else {

					if (j == 1) {
						S1 = ifour * (AF / H1) * bc->PHI2(HH1, st);
						S = S - S1;
						S1 = (AF / H1) * bc->PHI2(HHP1, st);
						S = S - S1;
						S1 = (AF / H1) * bc->PHI2(HHM1, st);
						S = S - S1;
					}

					if (j == m - 1) {
						S1 = ifour * (AF / H1) * bc->PHI4(HH1, st);
						S = S - S1;
						S1 = (AF / H1) * bc->PHI4(HHP1, st);
						S = S - S1;
						S1 = (AF / H1) * bc->PHI4(HHM1, st);
						S = S - S1;
					}
				}

			}

//			filestr << k << " E: S= [" << S.a << " ; " << S.b << "]" << endl;
			if (st == 0) {
				bm.ToMap(n2 - 1, S);
//				if (i < 10) cout << "S.Mid() = " << S.Mid() << endl;

				for (int i = 1; i <= n1; i++) {
					rh = r[i - 1];
					if (rh != 0)
					{
						B1M->ToMap(rh -1, bm.FromMap(i-1));
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
							S = S - (B1M->FromMap(i-1) * this->X[q - 1]);
							q = q + p;
						}
						if (!((S.a >= 0) && (S.b <= 0)))
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
							//cout << jh-1 << ": X= [" << this->X[jh - 1].a << " ; " << this->X[jh - 1].b << "]" << endl;
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
int PoissonSolver4OrderAM<T>::SolveDIA() {
	//	fstream filestr;
	//	string fname = "tmpLog.txt";
	//	filestr.open(fname.c_str(), fstream::out);

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

		const Interval<T> izero = { 0, 0 };
		const Interval<T> ione = { 1, 1 };
		const Interval<T> itwo = { 2, 2 };
		const Interval<T> ithree = { 3, 3 };
		const Interval<T> ifour = { 4, 4 };
		const Interval<T> itwelve = { 12, 12 };
		const Interval<T> imtwenty = { -20, -20 };
		const Interval<T> ithreesixty = { 360, 360 };
		const Interval<T> ininety = { 90, 90 };

		Interval<T> tmpi = { 0, 0 };
		int i, j, jh, j1, k, kh, l, lh, l1, l2, n1, n2, n3, p, q, rh;
		int num;
		Interval<T> AF, BB0, BB1, CF, H1, HH, HH1, II, JJ, K1, KK, KK1, MAX, MM, AA,
				CC, MMconst, NNconst, NN, S, S1, S2, S3, S4, S5, H1POW2, K1POW2,
				IIP1, IIM1, JJP1, JJM1, KKP1, KKM1, HHP1, HHM1, H1POW2K1POW2;
		bool list_exists;
		Interval<T> aij;
		int* r;
		T z;
		THashMap<T> bm;
		THashMap<T>* B1M = new THashMap<T>();

		st = 0;
		NN.a = n;
		NN.b = n;
		MM.a = m;
		MM.b = m;
		MMconst.a = bc->GetConstM();
		MMconst.b = -bc->GetConstM();
		NNconst.a = bc->GetConstN();
		NNconst.b = -bc->GetConstN();

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
			HH.b = H1.a;
			H1POW2 = H1 * H1;
			K1POW2 = K1 * K1;
			H1POW2K1POW2 = H1POW2 * K1POW2;
			KK.a = -K1.b;
			KK.b = K1.a;
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
			num = 1;

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
				IIP1.a = i + 1;
				IIP1.b = i + 1;
				JJP1.a = j + 1;
				JJP1.b = j + 1;
				IIM1.a = i - 1;
				IIM1.b = i - 1;
				JJM1.a = j - 1;
				JJM1.b = j - 1;

				HH1 = ALPHA + (II * H1);
				KK1 = BETA + (JJ * K1);

				HHM1 = ALPHA + (IIM1 * H1);
				KKM1 = BETA + (JJM1 * K1);

				HHP1 = ALPHA + (IIP1 * H1);
				KKP1 = BETA + (JJP1 * K1);
				AA = bc->A(HH1, KK1, st);
				CC = bc->C(HH1, KK1, st);
				AF = AA / H1;
				CF = CC / K1;
				S1 = imtwenty * AF / H1;

				//fourth order method - option 6 (Zhang)
				if (i > 1) {
					bm.ToMap(l1 - 1, ifour * AF / H1);

					if (j > 1) {
						bm.ToMap(l1 - 2, AF / H1);
					}
					if (j < m - 1) {
						bm.ToMap(l1, AF / H1);
					}

				}

				bm.ToMap(l2 - 1, S1);

				if (j > 1) {
					bm.ToMap(l2 - 2, ifour * (AF / H1));
				}
				if (j < m - 1)
					bm.ToMap(l2, ifour * (AF / H1));

				l1 = l2 + m - 1;

				if (i < n - 1) {
					bm.ToMap(l1 - 1, ifour * (AF / H1));
					if (j > 1) {
						bm.ToMap(l1 - 2, (AF / H1));
					}
					if (j < m - 1) {
						bm.ToMap(l1, (AF / H1));
					}
				}

				//right site

				//Zhang- option 6
				S = ifour * bc->F(HH1, KK1, st);
				S = S
						+ (ione / itwo)
								* (bc->F(HHP1, KK1, st) + bc->F(HHM1, KK1, st)
										+ bc->F(HH1, KKP1, st)
										+ bc->F(HH1, KKM1, st));
				cout << "k = " << k << "; S= [" << S.a << " ; " << S.b << "]" << endl;
	//			filestr << k << " B: S= [" << S.a << " ; " << S.b << "]" << endl;

				//second order truncation error
//				S5 = H1POW2 / itwelve;
//				S3 = S5 * bc->PSI(HH1, KK1, st) * bc->OMEGA(HH1, KK1, st) - itwo*MMconst;
//				S1 = S5 * S3;
//				S = S + S1;

				//fourth order truncation error
				S5 = H1POW2 * H1POW2  / ithreesixty;
				S3 = (bc->PSI1(HH1, KK1, st) + bc->PSI2(HH1, KK1, st));
				S1 = S5 * S3;
				S = S + S1;
				S5 = H1POW2 * H1POW2 / ininety;
				S3 = bc->PSI3(HH1, KK1, st);
				S1 = S5 * S3;
				S = S + S1;

				if (st == 0) {

					if (i == 1) {
						S1 = ifour * (AF / H1) * bc->PHI1(KK1, st);
						S = S - S1.Dual();

						if (j < m) {
							S1 = (AF / H1) * bc->PHI1(KKP1, st);
							S = S - S1.Dual();
						}
						if (j > 0) {
							S1 = (AF / H1) * bc->PHI1(KKM1, st);
							S = S - S1.Dual();
						}

						if (j == 1) {
							S1 = ifour * (AF / H1) * bc->PHI2(HH1, st);
							S = S - S1.Dual();

							S1 = (AF / H1) * bc->PHI2(HHP1, st);
							S = S - S1.Dual();
						}

						if (j == m - 1) {
							S1 = ifour * (AF / H1) * bc->PHI4(HH1, st);
							S = S - S1.Dual();

							S1 = (AF / H1) * bc->PHI4(HHP1, st);
							S = S - S1.Dual();
						}

					} else if (i == n - 1) {

						S1 = ifour * (AF / H1) * bc->PHI3(KK1, st);
						S = S - S1.Dual();
						if (j < m) {
							S1 = (AF / H1) * bc->PHI3(KKP1, st);
							S = S - S1.Dual();
						}
						if (j > 0) {
							S1 = (AF / H1) * bc->PHI3(KKM1, st);
							S = S - S1.Dual();
						}

						if (j == 1) {
							S1 = ifour * (AF / H1) * bc->PHI2(HH1, st);
							S = S - S1.Dual();
							S1 = (AF / H1) * bc->PHI2(HHM1, st);
							S = S - S1.Dual();
						}

						if (j == m - 1) {
							S1 = ifour * (AF / H1) * bc->PHI4(HH1, st);
							S = S - S1.Dual();
							S1 = (AF / H1) * bc->PHI4(HHM1, st);
							S = S - S1.Dual();
						}

					} else {
						if (j == 1) {
							S1 = ifour * (AF / H1) * bc->PHI2(HH1, st);
							S = S - S1.Dual();
							S1 = (AF / H1) * bc->PHI2(HHP1, st);
							S = S - S1.Dual();
							S1 = (AF / H1) * bc->PHI2(HHM1, st);
							S = S - S1.Dual();
						}

						if (j == m - 1) {
							S1 = ifour * (AF / H1) * bc->PHI4(HH1, st);
							S = S - S1.Dual();
							S1 = (AF / H1) * bc->PHI4(HHP1, st);
							S = S - S1.Dual();
							S1 = (AF / H1) * bc->PHI4(HHM1, st);
							S = S - S1.Dual();
						}
					}

				}

	//			filestr << k << " E: S= [" << S.a << " ; " << S.b << "]" << endl;
				if (st == 0) {
					S5 = ((AA * H1POW2) * H1POW2K1POW2) * MMconst;
					S3 = ((CC * K1POW2) * H1POW2K1POW2) * NNconst;
					S1 = S3 + S5;
					S = S + (S1 / itwelve);
					bm.ToMap(n2 - 1, S);
	//				if (i < 10) cout << "S.Mid() = " << S.Mid() << endl;

					for (int i = 1; i <= n1; i++) {
						rh = r[i - 1];
						if (rh != 0)
						{
							B1M->ToMap(rh -1, bm.FromMap(i-1));
						}
	//					if ((k > 1) && (rh == k - 1)) {
	//						BB1 = bm.FromMap(i - 1);
	//					}
	//					if ((k > m - 1) && (rh == k - m + 1)) {
	//						BB0 = bm.FromMap(i - 1);
	//					}
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
								S = S - (B1M->FromMap(i-1) * this->X[q - 1]);
								q = q + p;
							}
							if (!((S.a >= 0) && (S.b <= 0)))
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
						MAX = (ione * tmpi.Inverse());
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
								//cout << jh-1 << ": X= [" << this->X[jh - 1].a << " ; " << this->X[jh - 1].b << "]" << endl;
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
						this->X[rh - 1] = S.Projection();
						r[rh - 1] = rh;
					}
				}
			}
		}
	//	filestr.close();
		return 0;
}

//The explicit instantiation part
template class PoissonSolver4OrderAM<long double> ;
template class PoissonSolver4OrderAM<mpreal> ;

}
