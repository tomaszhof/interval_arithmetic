/*
 * NakaoExperiment2DApprox.cpp
 *
 *  Created on: Jun 10, 2020
 *      Author: numeric
 */

#include "NakaoExperiment2DApprox.h"
#include "Matrix2D.h"

int NakaoExperiment2DApprox::i = 0;
int NakaoExperiment2DApprox::j = 0;
long double NakaoExperiment2DApprox::h = 0.0;

NakaoExperiment2DApprox::NakaoExperiment2DApprox() {
	this->integrator = new GSLIntegrator();
}

NakaoExperiment2DApprox::~NakaoExperiment2DApprox() {
}


void NakaoExperiment2DApprox::execute() {

	double integral_err = 0.0;
	Interval<long double> intErr = {0.0, 0.0};
	long double tmp_int_res = 0.0;

	Interval<long double> iIntCij1 = {0.0, 0.0};
	Interval<long double> iIntCij2 = {0.0, 0.0};
	Interval<long double> iIntCij3 = {0.0, 0.0};
	Interval<long double> iIntCij4 = {0.0, 0.0};
	Interval<long double> iIntCij5 = {0.0, 0.0};
	Interval<long double> iIntCij6 = {0.0, 0.0};
	Interval<long double> iIntCij7 = {0.0, 0.0};


	//global setting
	cout.setf(std::ios_base::scientific);

	cout << "Solving the boundary value problem in general form:" << endl;
	cout
			<< "  -d^2u/dx^2 - d^2*u/dy^2 -c(x,y)*u = f(x,y),"
			<< endl;
	cout << "  u(x, 0) = u(0,y) = u(x,1) = u(1,y) = 0" << endl;
	cout << "by the Galerkin approximation and Nakao''s method." << endl;
	cout << endl;

	cout << "Enter the number of subinterval in [0,1]     n = ";
	cin >> n;

	cout << endl;
	cout << "Do you want to set integration precision (y - yes, n - no) answer = ";
	cin >> answer;

	if (answer == 'y'){
		cout << "integrator_calls = ";
		cin >> integrator_calls;
		integrator->setCalls(integrator_calls);
	}

	cout << endl;
	while ((output != 's') && (output != 'f')) {
		cout << "Select the output (s - screen, f - file) output = ";
		cin >> output;
	}

	if (output == 'f') {
		while (answer != 'y') {
			cout << "Enter the name of file" << endl;
			cout << "(.txt will be added automatically)     file_name = ";
			cin >> file_name;
			file_name = file_name + ".txt";
			cout << "The program will create a resulting file:       "
					<< file_name << endl;
			cout << "Do you accept this name (y - yes, n - no) answer = ";
			cin >> answer;
		}
		cout << endl;
		cout << "Please, wait ..." << endl;

		results.open(file_name.c_str(), fstream::out);
		results.setf(std::ios_base::scientific);

		results << "Solving the boundary value problem" << endl;
		results
				<< "  -d^2u/dx^2 - d^2*u/dy^2 -pi*u = (2*pi - 1) * sin(pi*x)*sin(pi*y),"
				<< endl;
		results << "  u(x, 0) = u(0,y) = u(x,1) = u(1,y) = 0" << endl;
		results << "by the Galerkin approximation and Nakao''s method." << endl;
		results << endl;
		results << "The number of subintervals in [0,1]  n = " << n << endl;
		results << endl;
		results << "Galerkin's approximation" << endl;
		results << endl;
	} else {
		cout << endl;
		cout << "Galerkin's approximation" << endl;
		cout << endl;
	}

	auto t1 = std::chrono::high_resolution_clock::now();
	h = 1.0 / n;
	n1 = (n - 1) * (n - 1);
	n2 = n1 + 1;
	p = n2;

	cout << "Preparing integrals coefficients... (it may take some time)" << endl;


	Matrix2D<Interval<long double>>* MIC_1 = new Matrix2D<Interval<long double>>(n2, n2);
	Matrix2D<Interval<long double>>* MIC_2 = new Matrix2D<Interval<long double>>(n2, n2);
	Matrix2D<Interval<long double>>* MIC_3 = new Matrix2D<Interval<long double>>(n2, n2);
	Matrix2D<Interval<long double>>* MIC_4 = new Matrix2D<Interval<long double>>(n2, n2);
	Matrix2D<Interval<long double>>* MIC_5 = new Matrix2D<Interval<long double>>(n2, n2);
	Matrix2D<Interval<long double>>* MIC_6 = new Matrix2D<Interval<long double>>(n2, n2);
	Matrix2D<Interval<long double>>* MIC_7 = new Matrix2D<Interval<long double>>(n2, n2);

	Matrix2D<long double>* MIF = new Matrix2D<long double>(n2, n2);
	Matrix2D<Interval<long double>>* MID = new Matrix2D<Interval<long double>>(n2, n2);

	for (i=0; i<n;i++){
		for (j=0; j<n;j++){
			tmp_int_res = integrator->integrate(&g_int_c_ij1, 0, 1, 0, 1,&integral_err)*(1.0 / (h*h));
			tmpstr = boost::lexical_cast<string>(tmp_int_res);
			iIntCij1 = IntRead<long double>(tmpstr);
			intErr.a = - integral_err;
			intErr.b = integral_err;
			iIntCij1 = iIntCij1 + im11 * intErr;
			MIC_1->item(i,j) = iIntCij1;

			tmp_int_res = integrator->integrate(&g_int_c_ij2, 0, 1, 0, 1,&integral_err)*(1.0 / (h*h));
			tmpstr = boost::lexical_cast<string>(tmp_int_res);
			iIntCij2 = IntRead<long double>(tmpstr);
			intErr.a = - integral_err;
			intErr.b = integral_err;
			iIntCij2 = iIntCij2 + im11 * intErr;
			MIC_2->item(i,j) = iIntCij2;

			tmp_int_res = integrator->integrate(&g_int_c_ij3, 0, 1, 0, 1,&integral_err)*(1.0 / (h*h));
			tmpstr = boost::lexical_cast<string>(tmp_int_res);
			iIntCij3 = IntRead<long double>(tmpstr);
			intErr.a = - integral_err;
			intErr.b = integral_err;
			iIntCij3 = iIntCij3 + im11 * intErr;
			MIC_3->item(i,j) = iIntCij3;

			tmp_int_res = integrator->integrate(&g_int_c_ij4, 0, 1, 0, 1,&integral_err)*(1.0 / (h*h));
			tmpstr = boost::lexical_cast<string>(tmp_int_res);
			iIntCij4 = IntRead<long double>(tmpstr);
			intErr.a = - integral_err;
			intErr.b = integral_err;
			iIntCij4 = iIntCij4 + im11 * intErr;
			MIC_4->item(i,j) = iIntCij4;

			tmp_int_res = integrator->integrate(&g_int_c_ij5, 0, 1, 0, 1,&integral_err)*(1.0 / (h*h));
			tmpstr = boost::lexical_cast<string>(tmp_int_res);
			iIntCij5 = IntRead<long double>(tmpstr);
			intErr.a = - integral_err;
			intErr.b = integral_err;
			iIntCij5 = iIntCij5 + im11 * intErr;
			MIC_5->item(i,j) = iIntCij5;


			tmp_int_res = integrator->integrate(&g_int_c_ij6, 0, 1, 0, 1,&integral_err)*(1.0 / (h*h));
			tmpstr = boost::lexical_cast<string>(tmp_int_res);
			iIntCij6 = IntRead<long double>(tmpstr);
			intErr.a = - integral_err;
			intErr.b = integral_err;
			iIntCij6 = iIntCij6 + im11 * intErr;
			MIC_6->item(i,j) = iIntCij6;


			tmp_int_res = integrator->integrate(&g_int_c_ij7, 0, 1, 0, 1,&integral_err)*(1.0 / (h*h));
			tmpstr = boost::lexical_cast<string>(tmp_int_res);
			iIntCij7 = IntRead<long double>(tmpstr);
			intErr.a = - integral_err;
			intErr.b = integral_err;
			iIntCij7 = iIntCij7 + im11 * intErr;
			MIC_7->item(i,j) = iIntCij7;

			tmp_int_res = integrator->integrate(&g_f_phi, 0, 1, 0, 1,&integral_err)*(-1.0 / (h*h));
			MIF->item(i,j) = tmp_int_res;
			intErr.a = - integral_err;
			intErr.b = integral_err;
			tmpstr = boost::lexical_cast<string>(tmp_int_res);
			MID->item(i,j) = IntRead<long double>(tmpstr) + intErr;

		}
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
	cout << "Integration finished! Start calculations." << endl;
	cout << endl << "Integration time: " << duration << " [ms]" << endl;
	results << endl << "Integration time: " << duration << " [ms]" << endl;


	r = new int[n1 + 1];
	for (int i = 1; i <= n2; ++i) {
		r[i - 1] = 0;
	}
	k = 0;
	j = 0;
	a1 = new long double[n1 + 1];
	b1 = new long double[n1 + 1];
	x = new long double[(n1 + 2) * (n1 + 2) / 4];
//	a = 4.0 / pow(h, 2.0) - PI / 2.0;
//	b = -1.0 / pow(h, 2.0) - PI / 12.0;
//	b_dash = (-1.0 * PI) / 12.0;
	d = (1.0 - 2.0 * PI) * (1.0 - cos(PI * h)) / pow(PI * h, 2.0);
	sph = sin(PI * h) / (PI * h);

	mhpm2 = -1.0/(h*h);

	finish = false;
	while (k != n1) {
		k = k + 1;
		for (i = 1; i <= n1; ++i) {
			a1[i - 1] = 0;
		}
		j = j + 1;
		i = (k - 1) / (n - 1) + 1;
		l1 = (i - 2) * (n - 1) + j;
		l2 = l1 + n - 1;


		aij = mhpm2*(4.0 - MIC_1->item(i,j).a);
		aij_l = mhpm2*(1.0 + MIC_2->item(i,j).a);
		aij_u = mhpm2*(1.0 + MIC_3->item(i,j).a);

		bij = mhpm2*(1.0 + MIC_5->item(i,j).a);
		bij_l = mhpm2*(MIC_7->item(i,j).a);

		b_dash_ij = mhpm2*(1.0 + MIC_4->item(i,j).a);
		b_dash_ij_u = mhpm2*(MIC_6->item(i,j).a);

		if (i > 1) {
			a1[l1 - 1] = b_dash_ij;
			if (j < n - 1) {
				a1[l1] = b_dash_ij_u;
			}
		}
		a1[l2 - 1] = aij;
		if (j > 1)
			a1[l2 - 2] = aij_l;
		if (j < n - 1)
			a1[l2] = aij_u;
		l1 = l2 + n - 1;
		if (i < n - 1) {
			a1[l1 - 1] = bij;
			if (j > 1)
				a1[l1 - 2] = bij_l;
		}


		a1[n2 - 1] = MIF->item(i,j);

		long double testInt = d
						* (cos(PI * (i + j) * h) - sph * cos(PI * (i - j) * h));
		long double test_f_ij = (1.0 - 2.0 * M_PI)/(M_PI*M_PI)*(1.0 - cos(M_PI*h))*(1.0 / (h*h));
		test_f_ij *= (cos(M_PI*(i+j)*h)-sph*cos((i-j)*h));
//		cout << "------------------------" << endl;
//		cout << "Integral(f_i=" << i << ", j=" << j << ")_exact = " << testInt << endl;
//		cout << "Integral(f_i=" << i << ", j=" << j << ")_approx = " << a1[n2 - 1] << endl;
//		cout << "Integral(f_i=" << i << ", j=" << j << ")_approx_err = " << integral_err << endl;
		for (i = 1; i <= n1; ++i) {
			rh = r[i - 1];
			if (rh != 0)
				b1[rh - 1] = a1[i - 1];
		}
		kh = k - 1;
		l = 0;
		max = 0;
		for (j1 = 1; j1 <= n2; ++j1) {
			if (r[j1 - 1] == 0) {
				s = a1[j1 - 1];
				l = l + 1;
				q =l;
				for (i = 1; i <= kh; ++i) {
					s = s - b1[i - 1] * x[q - 1];
					q = q + p;
				}
				a1[l - 1] = s;
				if (s < 0) {
					s = abs(s);
				}
				if ((j1 < n2) && (s > max)) {
					max = s;
					jh = j1;
					lh = l;
				}
			} //end if inside for
		} //end for loop (find max value)
		max = 1.0 / a1[lh - 1];
		r[jh - 1] = k;
		for (i = 1; i <= p; ++i) {
			a1[i - 1] = max * a1[i - 1];
		}
		jh = 0;
		q = 0;
		for (j1 = 1; j1 <= kh; ++j1) {
			s = x[q + lh - 1];
			for (i = 1; i <= p; ++i) {
				if (i != lh) {
					jh = jh + 1;
					x[jh - 1] = x[q + i - 1] - s * a1[i - 1];
				}
			}
			q = q + p;
		}
		for (i = 1; i <= p; ++i) {
			if (i != lh) {
				jh = jh + 1;
				x[jh - 1] = a1[i - 1];
			}
		}
		p = p - 1;
		if (j == n - 1)
			j = 0;
	}; //end of loop while(k!=n1)

	delete[] a1;
	delete[] b1;

	for (k = 1; k <= n1; ++k) {
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
	u = new long double*[n];
	for (i = 0; i <= n-1; ++i)
		u[i] = new long double[n];

	for (i = 1; i <= n - 1; ++i) {
		for (j = 1; j <= n - 1; ++j) {
			u[i][j] = x[(i - 1) * (n - 1) + j - 1];
		}
	}

	delete[] x;

	for (i = 1; i <= n - 1; ++i)
		for (j = 1; j <= n - 1; ++j) {
			exact = u_exact(i*h, j*h);
			if (output == 's') {
				cout << "u(" << i * h << "," << j * h << ") = " << u[i][j]
						<< endl;
				cout << "         exact = " << exact;
				cout << "   error = " << abs(u[i][j] - exact) << endl;
			} else {
				results << "u(" << i * h << "," << j * h << ") = " << u[i][j]
						<< endl;
				results << "         exact = " << exact;
				results << "   error = " << abs(u[i][j] - exact) << endl;
			}

		} //end of double for loop

	if (output == 's') {
		cout << "Press Enter to continue ..." << endl;
		std::getchar();
	}

	if (output == 's') {
		cout << endl;
		cout << "Nakao's method" << endl;
		cout << endl;
	} else {
		results << endl;
		results << "Nakao's method" << endl;
		results << endl;
	}

	ipi = Interval<long double>::IPi();
	tmpstr = boost::lexical_cast<string>(h);
	ih.a = LeftRead<long double>(tmpstr);
	ih.b = RightRead<long double>(tmpstr);

	ih2 = ISqr(ih, error);

	ia = i4/ih2;

	ib = im1 / ih2;

	id = ICos(ipi * ih);
	id = (i1 - id)/(ISqr(ipi * ih, error));
	id = (i1 - (i2*ipi)) * id;

	ic = im11 * (Interval<long double>::ISqr2()/i2);


	ipi12 = ipi/i12;

	ia1 = ISqr(ipi*ih, error) / i2;


	ib1 = i1 - ICos(ipi*ih);
	ib1 = ((i2/ipi) - i4) * ib1;

	ic1 = ISqr(((i2 * ipi) - i1), error) / i4;

	alpha_k = new long double*[n];
	alpha_km1 = new long double*[n];
	iu_k = new Interval<long double>*[n];
	iu_km1 = new Interval<long double>*[n];

	for (int i = 0; i < n; ++i) {
		alpha_k[i] = new long double[n];
		alpha_km1[i] = new long double[n];
		iu_k[i] = new Interval<long double>[n];
		iu_km1[i] = new Interval<long double>[n];
	}

	iu_km1[0][0] = i0;
	for (int i = 1; i <= n - 1; ++i) {
		for (int j = 1; j <= n - 1; ++j) {
			tmpstr = boost::lexical_cast<string>(u[i][j]);
			iu_km1[i][j].a = LeftRead<long double>(tmpstr);
			iu_km1[i][j].b = RightRead<long double>(tmpstr);
			alpha_km1[i][j] = 0.0;
		}
	}
	delete[] u;

	it = 0;
	finish = false;
	delta = 1e-8;
	epsilon = 1e-8;
	interval_a1 = new Interval<long double> [n1 + 1];
	interval_b1 = new Interval<long double> [n1 + 1];
	interval_x = new Interval<long double> [(n1 + 2) * (n1 + 2) / 4];


	while (!finish) {
		it = it + 1;
		cout << "it = " << it << endl;
		p = n2;
		for (i = 1; i <= n2; ++i) {
			r[i - 1] = 0;
		}

		k = 0;
		j = 0;

		while (k != n1) {
			k = k + 1;
			for (i = 1; i <= n1; ++i) {
				interval_a1[i - 1] = i0;
			}
			j = j + 1;

			i = (k - 1) / (n - 1) + 1;
			l1 = (i - 2) * (n - 1) + j;
			l2 = l1 + n - 1;

			//read all integrals for coefficients
			iIntCij1 = MIC_1->item(i,j);
			iIntCij2 = MIC_2->item(i,j);
			iIntCij3 = MIC_3->item(i,j);
			iIntCij4 = MIC_4->item(i,j);
			iIntCij5 = MIC_5->item(i,j);
			iIntCij6 = MIC_6->item(i,j);
			iIntCij7 = MIC_7->item(i,j);

			if (i > 1) {
				interval_a1[l1 - 1] = ib;
			}
			interval_a1[l2 - 1] = ia;
			if (j > 1) {
				interval_a1[l2 - 2] = ib;
			}
			if (j < n - 1) {
				interval_a1[l2] = ib;
			}

			l1 = l2 + n - 1;
			if (i < n - 1) {
				interval_a1[l1 - 1] = ib;
			}

			id1 = iIntCij1 * iu_km1[i][j]; //i6 * iu_km1[i][j];

			if (i == 1) {
				if (j == 1) {
					id1 = id1 + (iIntCij3*iu_km1[1][2] + iIntCij5*iu_km1[2][1]);
				} else {
					if (j == n - 1) {
						id1 = id1 + iIntCij2*iu_km1[1][n - 2] + iIntCij7*iu_km1[2][n - 2]
								+ iIntCij5*iu_km1[2][n - 1];
					} else {
						id1 = id1 + iIntCij2*iu_km1[1][j - 1] + iIntCij3*iu_km1[1][j + 1]
								+ iIntCij7*iu_km1[2][j - 1] + iIntCij5*iu_km1[2][j];
					}
				}
			} else if (i == n - 1) {
				if (j == 1) {
					id1 = id1 + iIntCij4*iu_km1[n - 2][1] + iIntCij6*iu_km1[n - 2][2]
							+ iIntCij3*iu_km1[n - 1][2];
				} else {
					if (j == n - 1) {
						id1 = id1 + iIntCij2*iu_km1[n - 1][n - 2] + iIntCij4*iu_km1[n - 2][n - 1];
					} else {
						id1 = id1 + iIntCij3*iu_km1[n - 1][j + 1] + iIntCij2*iu_km1[n - 1][j - 1]
								+ iIntCij6*iu_km1[n - 2][j + 1] + iIntCij4*iu_km1[n - 2][j];
					}
				}
			} else {
				if (j == 1) {
					id1 = id1 + iIntCij4*iu_km1[i - 1][1] + iIntCij6*iu_km1[i - 1][2]
							+ iIntCij3*iu_km1[i][2] + iIntCij5*iu_km1[i + 1][1];
				} else {
					if (j == n - 1) {
						id1 = id1 + iIntCij4*iu_km1[i - 1][n - 1] + iIntCij2*iu_km1[i][n - 2]
								+ iIntCij7*iu_km1[i + 1][n - 2] + iIntCij5*iu_km1[i + 1][n - 1];
					} else {
						id1 = id1 + iIntCij4*iu_km1[i - 1][j] + iIntCij6*iu_km1[i - 1][j + 1]
								+ iIntCij2*iu_km1[i][j - 1] + iIntCij3*iu_km1[i][j + 1]
								+ iIntCij7*iu_km1[i + 1][j - 1] + iIntCij5*iu_km1[i + 1][j];
					}
				}
			}

			iz = MID->item(i,j);

			id1 = id1 + iz + im11*intErr;

			tmpstr = boost::lexical_cast<string>(alpha_km1[i][j]);
			iz.a = LeftRead<long double>(tmpstr);
			iz.b = RightRead<long double>(tmpstr);

			interval_a1[n2 - 1] = id1 + (ic * iz);
			for (int i = 1; i <= n1; ++i) {
				rh = r[i - 1];
				if (rh != 0) {
					interval_b1[rh - 1] = interval_a1[i - 1];
				}
			}

			kh = k - 1;
			l = 0;
			imax = i0;


			for (j1 = 1; j1 <= n2; ++j1) {
				if (r[j1 - 1] == 0) {
					interval_s = interval_a1[j1 - 1];
					l = l + 1;
					q = l;
					for (int i = 1; i <= kh; ++i) {
						interval_s = interval_s
								- (interval_b1[i - 1] * interval_x[q - 1]);
						q = q + p;
					}
					interval_a1[l - 1] = interval_s;

					if (interval_s.a < 0) {
						interval_s.a = abs(interval_s.a);
					}
					if (interval_s.b < 0) {
						interval_s.b = abs(interval_s.b);
					}

					if (interval_s.b < interval_s.a) {
						s = interval_s.a;
						interval_s.a = interval_s.b;
						interval_s.b = s;
					}
					if ((j1 < n2) && (interval_s.b > imax.a)) {
						imax = interval_s;
						jh = j1;
						lh = l;
					}
				}
			}
			imax = i1 / interval_a1[lh - 1];
			r[jh - 1] = k;

			for (i = 1; i <= p; ++i) {
				interval_a1[i - 1] = imax * interval_a1[i - 1];
			}
			jh = 0;
			q = 0;
			for (j1 = 1; j1 <= kh; ++j1) {
				interval_s = interval_x[q + lh - 1];
				for (int i = 1; i <= p; ++i) {
					if (i != lh) {
						jh = jh + 1;
						interval_x[jh - 1] = interval_x[q + i - 1]
								- (interval_s * interval_a1[i - 1]);
					}
				}
				q = q + p;
			}

			for (int i = 1; i <= p; ++i) {
				if (i != lh) {
					jh = jh + 1;
					interval_x[jh - 1] = interval_a1[i - 1];
				}
			}
			p = p - 1;
			if (j == n - 1)
				j = 0;
		} //inner while loop

		for (k = 1; k <= n1; ++k) {
			rh = r[k - 1];
			if (rh != k) {
				interval_s = interval_x[k - 1];
				interval_x[k - 1] = interval_x[rh - 1];
				i = r[rh - 1];

				while (i != k) {
					interval_x[rh - 1] = interval_x[i - 1];
					r[rh - 1] = rh;
					rh = i;
					i = r[rh - 1];
				}
				interval_x[rh - 1] = interval_s;
				r[rh - 1] = rh;
			}
		}

		for (int i = 1; i <= n - 1; ++i)
			for (int j = 1; j <= n - 1; ++j) {
				iu_k[i][j] = interval_x[(i - 1) * (n - 1) + j - 1];
			}

		if (output == 's') {
			for (int i = 1; i <= n - 1; ++i) {
				for (int j = 1; j <= n - 1; ++j) {
					iu_k[i][j].IEndsToStrings(left, right);
					cout << "    u(" << i * h << ", " << j * h << ") = ["
							<< left << ", " << right << "]" << endl;
					exact = u_exact(i*h, j*h);
					cout << "          exact = " << exact << endl;
				}
			}
		} // if output = 's'

		for (int i = 1; i <= n - 1; ++i) {
			for (int j = 1; j <= n - 1; ++j) {
				d = cos(PI * (i + j) * h) - sph * cos(PI * (i - 1) * h);
				tmpstr = boost::lexical_cast<string>(d);
				ibeta.a = LeftRead<long double>(tmpstr);
				ibeta.b = RightRead<long double>(tmpstr);
				ibeta = ib1 * ibeta * iu_km1[i][j];
				ibeta = ibeta + (ia1 * iu_km1[i][j] * iu_km1[i][j]);
				ibeta = ibeta + ic1;
				ibeta = interval_arithmetic::ISqrt(ibeta, error);
				if (error == 0) {
					tmpstr = boost::lexical_cast<string>(
					PI * h * alpha_km1[i][j]);
					id1.a = LeftRead<long double>(tmpstr);
					id1.b = RightRead<long double>(tmpstr);
					ibeta = ih * (ibeta + id1);
					alpha_k[i][j] = ibeta.b;
				} else {
					cout << "Impossible to calculate Sqrt(ibeta) for i=" << i
							<< " and j=" << j << endl;
				}
			}
		} //end of double "for" loop (i,j)

//		if (output == 's') {
//			for (int i = 1; i <= n - 1; ++i) {
//				for (int j = 1; j <= n - 1; ++j) {
//					cout << "alpha(" << i * h << ", " << j * h << ") = "
//							<< alpha_k[i][j] << endl;
//				}
//			}
//		} // if output = 's'

		norm_u = 0;
		for (int i = 1; i <= n - 1; ++i)
			for (int j = 1; j <= n - 1; ++j) {
				norm_uij = abs(iu_k[i][j].a - iu_km1[i][j].a);
				if (abs(iu_k[i][j].b - iu_km1[i][j].b) > norm_uij)
					norm_uij = abs(iu_k[i][j].b - iu_km1[i][j].b);
				if (norm_uij > norm_u) {
					norm_u = norm_uij;
				}
			}

		if (norm_u < epsilon)
			u_OK = true;
		else
			u_OK = false;

		abs_alpha = 0;
		for (int i = 1; i <= n - 1; ++i)
			for (int j = 1; j <= n - 1; ++j) {
				abs_alphaij = abs(alpha_k[i][j] - alpha_km1[i][j]);
				if (abs_alphaij > abs_alpha) {
					abs_alpha = abs_alphaij;
				}
			}

		if (abs_alpha < epsilon)
			alpha_OK = true;
		else
			alpha_OK = false;

		if (u_OK && alpha_OK) {
			finish = true;
			if (output == 's') {
				cout << "Press Enter to continue..." << endl;
				std::getchar();
				cout << "Interval solution after-delta extension:" << endl
						<< endl;
			} else {
				results << "Number of iterations: " << it << endl;
				results << "Interval solution (after-delta extension):" << endl
						<< endl;
			}

			for (int i = 1; i <= n - 1; ++i)
				for (int j = 1; j <= n - 1; ++j) {
					iu_k[i][j].a = iu_k[i][j].a - delta;
					iu_k[i][j].b = iu_k[i][j].b + delta;
					iu_k[i][j].IEndsToStrings(left, right);

					if (output == 's') {
						cout << "    u(" << i * h << ", " << j * h << ") = ["
								<< left << ", " << right << "]" << endl;
						cout << "                  = [" << iu_k[i][j].a << ", "
								<< iu_k[i][j].b << "]" << endl;
					} else {
						results << "    u(" << i * h << ", " << j * h << ") = ["
								<< left << ", " << right << "]" << endl;
						results << "                  = [" << iu_k[i][j].a
								<< ", " << iu_k[i][j].b << "]" << endl;
					}

					exact =u_exact(i*h, j*h);
					width = abs(iu_k[i][j].b - iu_k[i][j].a);
					if (output == 's')
						cout << "           exact = " << exact << endl;
					else
						results << "        exact = " << exact << endl;

					if (output == 's')
						cout << "           width = " << width << endl;
					else
						results << "        width = " << width << endl;
				} // end for (i,j) loop
		}

		if (!finish) {
			for (int i = 1; i <= n - 1; ++i)
				for (int j = 1; j <= n - 1; ++j) {
					iu_km1[i][j] = iu_k[i][j];
					alpha_km1[i][j] = alpha_k[i][j];
				} // end of for (i,j) loop

			if (output == 's') {
				//cout << "Press Enter to continue..." << endl;
				//std::getchar();
			}
		}
	} //end of while(!finish) loop


	delete[] interval_b1;
	delete[] interval_x;
	delete[] r;

	if (output == 'f')
		results.close();

	cout << endl;
	cout << "END OF PROGRAM" << endl;
	std::getchar();

	return;
}

