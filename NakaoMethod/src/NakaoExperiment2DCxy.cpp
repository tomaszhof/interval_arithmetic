/*
 * NakaoExperiment2DCxy.cpp
 *
 *  Created on: Jul 21, 2018
 *      Author: numeric
 */

#include "NakaoExperiment2DCxy.h"

NakaoExperiment2DCxy::NakaoExperiment2DCxy() {
	// TODO Auto-generated constructor stub

}

NakaoExperiment2DCxy::~NakaoExperiment2DCxy() {
	// TODO Auto-generated destructor stub
}
void NakaoExperiment2DCxy::execute() {
	//global setting
	cout.setf(std::ios_base::scientific);

	cout << "Solving the boundary value problem" << endl;
	cout
			<< "  -d^2u/dx^2 - d^2*u/dy^2 -20*x*y*u = (2*pi - 1) * sin(pi*x)*sin(pi*y),"
			<< endl;
	cout << "  u(x, 0) = u(0,y) = u(x,1) = u(1,y) = 0" << endl;
	cout << "by the Galerkin approximation and Nakao''s method." << endl;
	cout << endl;

	cout << "Enter the number of subinterval in [0,1]     n = ";
	cin >> n;

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
				<< "  -d^2u/dx^2 - d^2*u/dy^2 -20*x*y*u = (2*pi - 1) * sin(pi*x)*sin(pi*y),"
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

	h = 1.0 / n;
	n1 = (n - 1) * (n - 1);
	n2 = n1 + 1;
	p = n2;
	r = new int[n1 + 1];
	for (int i = 1; i <= n2; ++i) {
		r[i - 1] = 0;
	}
	k = 0;
	j = 0;
	a1 = new long double[n1 + 1];
	b1 = new long double[n1 + 1];
	x = new long double[(n1 + 2) * (n1 + 2) / 4];
	a = 4.0L / (h * h);
	b = -1.0L / (h * h);
	d = (1.0L - 2.0L * PI) * (1.0L - cos(PI * h)) / pow(PI * h, 2.0);
	sph = sin(PI * h) / (PI * h);
	twenty_h2 = 20.0L * h * h;
	finish = false;
	while (k != n1) {
		k = k + 1;
		for (i = 1; i <= n1; ++i) {
			a1[i - 1] = 0.0;
		}
		j = j + 1;
		i = (k - 1) / (n - 1) + 1;
		l1 = (i - 2) * (n - 1) + j;
		l2 = l1 + n - 1;
		i_over_24 = i / 24.0L;
		j_over_24 = j / 24.0L;
		ij_over_12 = i * j / 12.0L;
		cout << "i=" << i << endl << "j=" << j << endl;
		if (i > 1) {
			a1[l1 - 1] = b
					- twenty_h2 * (ij_over_12 - j_over_24 - 1.0L / 360.0L);
			if (j < n - 1) {
				a1[l1] = -twenty_h2
						* (ij_over_12 + i_over_24 - j_over_24 - 1.0L / 45.0L);
			}
		}
		a1[l2 - 1] = a - twenty_h2 * (i * j / 2.0L - 1.0L / 36.0L);
		if (j > 1)
			a1[l2 - 2] = b
					- twenty_h2 * (ij_over_12 - i_over_24 - 1.0L / 360.0L);
		if (j < n - 1)
			a1[l2] = b - twenty_h2 * (ij_over_12 + i_over_24 - 1.0L / 360.0L);
		l1 = l2 + n - 1;
		if (i < n - 1) {
			a1[l1 - 1] = b
					- twenty_h2 * (ij_over_12 + j_over_24 - 1.0L / 360.0L);
			if (j > 1)
				a1[l1 - 2] = -twenty_h2
						* (ij_over_12 - i_over_24 + j_over_24 - 1.0L / 45.0L);
		}
		a1[n2 - 1] = d * (cos(PI * (i + j) * h) - sph * cos(PI * (i - j) * h));
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
				q = l;
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
		cout << max << endl;
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
	for (i = 0; i <= n - 1; ++i)
		u[i] = new long double[n];

	for (i = 1; i <= n - 1; ++i) {
		for (j = 1; j <= n - 1; ++j) {
			u[i][j] = x[(i - 1) * (n - 1) + j - 1];
		}
	}

	delete[] x;

	for (i = 1; i <= n - 1; ++i)
		for (j = 1; j <= n - 1; ++j) {
			if (output == 's') {
				cout << "u(" << i * h << "," << j * h << ") = " << u[i][j]
						<< endl;
			} else {
				results << "u(" << i * h << "," << j * h << ") = " << u[i][j]
						<< endl;
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

	ia = i4 / ih2;

	ib = im1 / ih2;

	id = ICos(ipi * ih);
	id = (i1 - id) / (ISqr(ipi * ih, error));
	id = (i1 - (i2 * ipi)) * id;

	ic = im11 * (Interval<long double>::ISqr2() / i2);

	i20h4 = (i20 * (ih2 * ih2));

	isph = ISin(ipi * ih);
	isph = (isph / (ipi * ih));

	alpha_k = new long double*[n];
	alpha_km1 = new long double*[n];
	iu_k = new Interval<long double>*[n];
	iu_km1 = new Interval<long double>*[n];

	for (int i = 0; i < n; ++i) {
		alpha_k[i] = new long double[n];
		alpha_km1[i] = new long double[n];
		iu_k[i] = new Interval<long double> [n];
		iu_km1[i] = new Interval<long double> [n];
	}

	//iu_km1[0][0] = i0; ??
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
	delta = 1e-3;
	epsilon = 1e-4;
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

			interval_i.a = i;
			interval_i.b = i;
			interval_j.a = j;
			interval_j.b = j;

			interval_i_over_24 = interval_i / i24;
			interval_j_over_24 = interval_j / i24;
			interval_ij_over_12 = (interval_i * interval_j) / i12;

			iz = ((interval_i * interval_j) / i2) - (i1 / i36);
			iz = i20h4 * iz;
			id1 = iz * iu_km1[i][j];

			if (i == 1) {
				if (j == 1) {
					iz = interval_ij_over_12 + interval_i_over_24;
					iz = i20h4 * (iz - (i1 / i360));
					id1 = id1 + (iz * iu_km1[1][2]);

					iz = interval_ij_over_12 + interval_j_over_24;
					iz = i20h4 * (iz - (i1 / i360));
					id1 = id1 + (iz * iu_km1[2][1]);
				} else {
					if (j == n - 1) {
						iz = interval_ij_over_12 - interval_i_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[1][n - 2];

						iz = interval_ij_over_12 - interval_i_over_24
								+ interval_j_over_24;
						iz = i20h4 * (iz - i1 / i45);
						id1 = id1 + iz * iu_km1[2][n - 2];

						iz = interval_ij_over_12 + interval_j_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[2][n - 1];
					} else {
						iz = interval_ij_over_12 - interval_i_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[1][j - 1];

						iz = interval_ij_over_12 + interval_i_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[1][j + 1];

						iz = interval_ij_over_12 - interval_i_over_24;
						iz = iz + interval_j_over_24;
						iz = i20h4 * (iz - i1 / i45);
						id1 = id1 + iz * iu_km1[2][j - 1];

						iz = interval_ij_over_12 + interval_j_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[2][j];
					}
				}
			} else if (i == n - 1) {
				if (j == 1) {
					iz = interval_ij_over_12 - interval_j_over_24;
					iz = i20h4 * (iz - (i1 / i360));
					id1 = id1 + iz * iu_km1[n - 2][1];

					iz = interval_ij_over_12 + interval_i_over_24;
					iz = iz - interval_j_over_24;
					iz = i20h4 * (iz - i1 / i45);
					id1 = id1 + iz * iu_km1[n - 2][2];

					iz = interval_ij_over_12 + interval_i_over_24;
					iz = i20h4 * (iz - (i1 / i360));
					id1 = id1 + iz * iu_km1[n - 1][2];
				} else {
					if (j == n - 1) {
						iz = interval_ij_over_12 - interval_j_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[n - 2][n - 1];

						iz = interval_ij_over_12 - interval_i_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[n - 1][n - 2];
					} else {
						iz = interval_ij_over_12 - interval_j_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[n - 2][j];

						iz = interval_ij_over_12 + interval_i_over_24;
						iz = iz - interval_j_over_24;
						iz = i20h4 * (iz - i1 / i45);
						id1 = id1 + iz * iu_km1[n - 2][j + 1];

						iz = interval_ij_over_12 - interval_i_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[n - 1][j - 1];

						iz = interval_ij_over_12 + interval_i_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[n - 1][j + 1];

					}
				}
			} else {
				if (j == 1) {
					iz = interval_ij_over_12 - interval_j_over_24;
					iz = i20h4 * (iz - (i1 / i360));
					id1 = id1 + iz * iu_km1[i - 1][1];

					iz = interval_ij_over_12 + interval_i_over_24;
					iz = iz - interval_j_over_24;
					iz = i20h4 * (iz - i1 / i45);
					id1 = id1 + iz * iu_km1[i - 1][2];

					iz = interval_ij_over_12 + interval_i_over_24;
					iz = i20h4 * (iz - (i1 / i360));
					id1 = id1 + iz * iu_km1[i][2];

					iz = interval_ij_over_12 + interval_j_over_24;
					iz = i20h4 * (iz - (i1 / i360));
					id1 = id1 + iz * iu_km1[i + 1][1];
				} else {
					if (j == n - 1) {
						iz = interval_ij_over_12 - interval_j_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[i - 1][n - 1];

						iz = interval_ij_over_12 - interval_i_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[i][n - 2];

						iz = interval_ij_over_12 - interval_i_over_24;
						iz = iz + interval_j_over_24;
						iz = i20h4 * (iz - i1 / i45);
						id1 = id1 + iz * iu_km1[i + 1][n - 2];

						iz = interval_ij_over_12 + interval_j_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[i + 1][n - 1];

					} else {
						iz = interval_ij_over_12 - interval_j_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[i - 1][j];

						iz = interval_ij_over_12 + interval_i_over_24;
						iz = iz - interval_j_over_24;
						iz = i20h4 * (iz - i1 / i45);
						id1 = id1 + iz * iu_km1[i - 1][j + 1];

						iz = interval_ij_over_12 - interval_i_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[i][j - 1];

						iz = interval_ij_over_12 + interval_i_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[i][j + 1];

						iz = interval_ij_over_12 - interval_i_over_24;
						iz = iz + interval_j_over_24;
						iz = i20h4 * (iz - i1 / i45);
						id1 = id1 + iz * iu_km1[i + 1][j - 1];

						iz = interval_ij_over_12 + interval_j_over_24;
						iz = i20h4 * (iz - (i1 / i360));
						id1 = id1 + iz * iu_km1[i + 1][j];

					}
				}
			}

			i_plus_j_ih = (interval_i + interval_j) * ih;
			i_minus_j_ih = (interval_i - interval_j) * ih;
			iz = isph * ICos(ipi * i_minus_j_ih);
			iz =  ICos(ipi * i_plus_j_ih) - iz;
			id1 = (id1/ih2) + (id * iz);

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
//					exact = (1 / PI) * sin(PI * i * h) * sin(PI * j * h);
//					cout << "          exact = " << exact << endl;
				}
			}
		} // if output = 's'

		for (int i = 1; i <= n - 1; ++i) {
			for (int j = 1; j <= n - 1; ++j) {
				interval_i.a = i;
				interval_i.b = i;
				interval_j.a = j;
				interval_j.b = j;
				interval_ij = interval_i * interval_j;
				iz = interval_ij * interval_ij;
				iz = (iz * (interval_i + interval_j)) / i2;
				ibeta1 = im1 * iz;
				iz = (interval_i * interval_i) + (interval_ij / i2);
				iz = iz + (interval_j * interval_j);
				iz = (interval_ij * iz) / i3;
				ibeta1 = ibeta1 - iz;
				iz = interval_i*interval_i*interval_i;
				iz = (iz / i3) + (interval_ij * (interval_i + interval_j));
				iz = iz + (ISqr(interval_j, error) * interval_j / i3);
				ibeta1 = ibeta1 - (iz / i4);
				iz = (ISqr(interval_i, error) / i3) + interval_ij;
				iz = iz + (ISqr(interval_j, error) / i3);
				ibeta1 = ibeta1 - (iz / i15);

				iz = ((interval_i + interval_j) / i180);
				ibeta1 = (ibeta1 - iz) + (i1 / i120);

				iz = ih2 * ih2 * ih2;
				ibeta1 = iz * ibeta1;
				ibeta1 = ISqr((i20 * iu_km1[i][j]), error) * ibeta1;
				i_plus_j_ih = (interval_i + interval_j) * ih;
				i_minus_j_ih = (interval_i - interval_j) * ih;
				ipih = ipi * ih;
				i2pih = i2 * ipih;
				iz1 = ih * (ICos(i2pih) + ICos(ipih));
				iz1 = iz1 - (i3 * ISin(i2pih)) / (i2 * ipi);
				iz = (interval_i - interval_j) * iz1;
				iz1 = (i3 * interval_i) - (i7 * interval_j) + i4;
				iz1 = (iz1 * ISin(ipih) / ipi);
				iz = iz + iz1;
				iz1 = ipi * i_minus_j_ih;
				ibeta2 = (iz * ISin(iz1)) / (i2 * ipi);
				iz = interval_ij * ih;
				iz1 = i5 / (i2 * (ISqr(ipi, error) * ih));
				iz = iz + iz1;
				iz1 = (i2 * ISin(ipih)) - ISin(i2pih);
				iz = (iz * iz1) + (ih * (ISin(i2pih)));
				iz1 = (i3 * ICos(i2pih)) - ICos(ipih);
				iz = iz + (iz1 / ipi);
				iz1 = (ipi * i_minus_j_ih);
				ibeta2 = ibeta2 + (iz * ICos(iz1)) / (i2 * ipi);
				iz1 = (i1 / ipi) - ((ih * ISin(ipih)) / i2);
				iz = (interval_i + interval_j) * iz1;
				iz1 = ((i2 * interval_i) * ICos(ipih)) / ipi;
				iz = iz - iz1;
				iz1 = ipi * i_plus_j_ih;
				ibeta2 = ibeta2 + ih * iz * ISin(iz1);
				iz1 = (i1 / ipi) - ih * ISin(ipih);
				iz = iz1 / ipi;
				iz1 = (interval_ij + (i1 / i6)) * ICos(ipih);
				iz1 = iz1 + (i1 / i3) - interval_ij;
				iz1 = iz + (ISqr(ih, error) * iz1);
				iz1 = ipi * i_plus_j_ih;
				ibeta2 = ibeta2 + iz * ICos(iz1);
				iz1 = i_minus_j_ih * ISin(ipi * interval_j * ih);
				iz1 = iz1 - (ICos(ipi * interval_j * ih) / ipi);
				iz = (i2 * iz1) / ipi;
				iz1 = ICos(ipi * ih) * ICos(ipi * interval_i * ih);
				ibeta2 = ibeta2 + (iz * iz1);
				iz = i20 * ((i2 * ipi) - i1);
				iz = (iz * iu_km1[i][j]) / ISqr(ipi, error);
				ibeta2 = iz * ibeta2;
				iz = ISqr(((i2 * ipi) - i1), error);
				ibeta = ibeta1 + (i2 * ibeta2) + (iz / i4);
				ibeta = ISqrt(ibeta, error);

				if (error == 0) {
					iz = (i20 * ih) * (alpha_km1[i][j]);
					ibeta = ih * (ibeta + iz);
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
				} // end for (i,j) loop
		}

		if (!finish) {
			for (int i = 1; i <= n - 1; ++i)
				for (int j = 1; j <= n - 1; ++j) {
					iu_km1[i][j] = iu_k[i][j];
					alpha_km1[i][j] = alpha_k[i][j];
				} // end of for (i,j) loop

			if (output == 's') {
				cout << "Press Enter to continue..." << endl;
				std::getchar();
			}
		}
	} //end of while(!finish) loop

//	for (int i = 0; i < n; ++i) {
//			delete[] alpha_k[i];
//			delete[] alpha_km1[i];
//			delete[] iu_k[i];
//			delete[] iu_km1[i];
//		}

	//delete[] interval_a1;
	delete[] interval_b1;
	delete[] interval_x;
	delete[] r;
	//delete[] iu_km1;
	//delete[] iu_k;

	if (output == 'f')
		results.close();

	cout << endl;
	cout << "END OF PROGRAM" << endl;
	std::getchar();

	return;
}

