/*
 * NakaoExperimentApprox.cpp
 *
 *  Created on: Feb 25, 2018
 *      Author: numeric
 */

#include "NakaoExperimentApprox.h"
#include "Integration.h"

using namespace integration;

NakaoExperimentApprox::NakaoExperimentApprox() {
	// TODO Auto-generated constructor stub
	output = NULL;
	answer = NULL;
}

NakaoExperimentApprox::~NakaoExperimentApprox() {
	// TODO Auto-generated destructor stub
}

long double phi(int i, long double x) {
	long double h = xv[1] - xv[0];
	if ((x>=xv[i-1])&&(x<=xv[i+1]))
		return (1.0L - (std::abs(x-xv[i])/h));

	return 0.0L;
}

long double phi_p(int i, long double x) {
	long double h = xv[1] - xv[0];
	if ((x>=xv[i-1])&&(x<=xv[i]))
		return (1.0L/h);

	if ((x>xv[i])&&(x<=xv[i+1]))
		return (-1.0L/h);

	return 0.0L;
}

void NakaoExperimentApprox::initialize() {
}

void NakaoExperimentApprox::initialize_vx(long double a, long double b, int n) {
	h = (b-a)/n;
	for (int i=0; i<=n; ++i)
		{
			xv.push_back(a + h*i);
		}
}

long double f(long double x) {
	return (M_PI - 1)*std::sin(M_PI*x);
}

void NakaoExperimentApprox::execute() {
	//global setting
	cout.setf(std::ios_base::scientific);

	cout << "Solving the boundary value problem" << endl;
	cout << "  -u'''' - pi*u = (pi - 1) * sin(pi*x), u(0) = u(1) = 0" << endl;
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
			cout << "Do you accept this name (y - yes, n - no) answetr = ";
			cin >> answer;
		}
		cout << endl;
		cout << "Please, wait ..." << endl;

		results.open(file_name.c_str(), fstream::out);
		results.setf(std::ios_base::scientific);

		results << "Solving the boundary value problem" << endl;
		results << "  -u'''' - pi*u = (pi - 1) * sin(pi*x), u(0) = u(1) = 0"
				<< endl;
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

	u = new long double[n];
	y = new long double[n];
	z = new long double[n];

	h = 1.0 / n;
	b = -1.0 / (h * h) - M_PIl / 6.0;
	a = 2.0 * (1.0 / (h * h) - M_PIl / 3.0);
	alpha = std::pow(2.0 * std::sin(M_PIl * h / 2.0) / (M_PIl * h), 2)
			* (M_PIl - 1.0);

	initialize_vx(0.0L, 1.0L, n);
//	c = -M_PI;
//	long double f_1 = GaussLegendreTH(f, 4, GL_X4, GL_A4, xv[0], xv[2], error);
//	long double d_1 = 1.0/h * f_1;
//	a = GaussLegendreTH(phi_p, 1, 1, 4, GL_X4, GL_A4, 0.0L, 1.0L, error);
//	a = a + c*GaussLegendreTH(phi, 1, 1, 4, GL_X4, GL_A4, 0.0L, 1.0L, error);
//	a = a/h;

	long double d_0 = GaussLegendreTH(f, phi, 1, 4, GL_X4, GL_A4, xv[0], xv[2], error) / h;
	y[0] = d_0 * sin(M_PIl * h) / a;
	z[0] = b / a;
	for (int i = 1; i < n - 2; ++i)
		z[i] = b / (a - b * z[i - 1]);

	for (int i = 0; i < n - 2; ++i) {
		long double d_i = GaussLegendreTH(f, phi, i+1, 4, GL_X4, GL_A4, xv[i], xv[i+2], error) / h;
		y[i + 1] = d_i * std::sin((i + 2.0) * M_PIl * h);
		y[i + 1] = (y[i + 1] - b * y[i]) / (a - b * z[i]);
	}
	u[n - 1 - 1] = y[n - 1 - 1];

	for (int i = n - 2 - 1; i >= 0; --i) {
		u[i] = y[i] - z[i] * u[i + 1];
	}

	for (int i = 0; i < n - 1; ++i) {
		exact = (1.0 / M_PIl) * std::sin(M_PIl * (i + 1) * h);
		if (output == 's') {
			cout << "u(" << i * h << ") = " << u[i];
			cout << endl;
			cout << "  exact = " << exact;
			cout << "  error = " << abs(u[i] - exact) << endl;
		} else {
			results << "u(" << i * h << ") = " << u[i];
			results << endl;
			results << "  exact = " << exact;
			results << "  error = " << abs(u[i] - exact) << endl;
		}
	}

	if (output == 's') {
		cout << endl;
		cout << "Press Enter to continue ..." << endl;
		std::getchar();
	}

	delete[] y;
	delete[] z;

	if (output == 's') {
		cout << endl;
		cout << "Nakao's method" << endl;
		cout << endl;
	} else {
		results << endl;
		results << "Nakao's method" << endl;
		results << endl;
	}

	a = 2.0 / (h * h);
	string astr = boost::lexical_cast<string>(a);
	ia.a = LeftRead<long double>(astr);
	ia.b = RightRead<long double>(astr);

	b = -1.0 / (h * h);
	string bstr = boost::lexical_cast<string>(b);
	ib.a = LeftRead<long double>(bstr);
	ib.b = RightRead<long double>(bstr);

	d = 2.0 * std::sin(M_PIl * h / 2.0) / (M_PIl * h);
	d = d * d * (M_PIl - 1.0);
	string dstr = boost::lexical_cast<string>(d);
	id.a = LeftRead<long double>(dstr);
	id.b = RightRead<long double>(dstr);

	c = std::sqrt(2.0 * h / 3.0);
	string cstr = boost::lexical_cast<string>(c);
	ic.a = LeftRead<long double>(cstr);
	ic.b = RightRead<long double>(cstr);

	string pistr = boost::lexical_cast<string>(M_PIl);
	iPi.a = LeftRead<long double>(pistr);
	iPi.b = RightRead<long double>(pistr);

	string hstr = boost::lexical_cast<string>(h);
	ih.a = LeftRead<long double>(hstr);
	ih.b = RightRead<long double>(hstr);

	string a1str = boost::lexical_cast<string>(2.0 * M_PIl * M_PIl * h / 3.0);
	ia1.a = LeftRead<long double>(a1str);
	ia1.b = RightRead<long double>(a1str);

	string b1str = boost::lexical_cast<string>(8.0 * (1.0 - 1.0/M_PIl)*std::pow(std::sin(M_PIl * h/2.0), 2.0));
	ib1.a = LeftRead<long double>(b1str);
	ib1.b = RightRead<long double>(b1str);

	string c1str = boost::lexical_cast<string>(std::pow(M_PIl - 1.0, 2) / 2.0);
	ic1.a = LeftRead<long double>(c1str);
	ic1.b = RightRead<long double>(c1str);

	alpha_k = new long double[n];
	alpha_km1 = new long double[n];

	iu_k = new Interval<long double> [n];
	iu_km1 = new Interval<long double> [n];
	iy = new Interval<long double> [n];
	iz = new Interval<long double> [n];

	for (int i = 0; i < n; ++i) {
		alpha_k[i] = 0.0;
		alpha_km1[i] = 0.0;

		iu_k[i].a = 0.0;
		iu_k[i].b = 0.0;

		iu_km1[i].a = 0.0;
		iu_km1[i].b = 0.0;

		iy[i].a = 0.0;
		iy[i].b = 0.0;

		iz[i].a = 0.0;
		iz[i].b = 0.0;
	}

	for (int i = 0; i < n - 1; ++i) {
		string uistr = boost::lexical_cast<string>(u[i]);
		iu_km1[i].a = LeftRead<long double>(uistr);
		iu_km1[i].b = RightRead<long double>(uistr);
	}

	k = 0;
	finish = false;
	delta = 1e-8;
	epsilon = 1e-8;

	while (!finish) {
		k++;
		if (output == 's') {
			cout << "k = " << k << endl;
			cout << endl;
		} else {
			results << "k = " << k << endl;
			results << endl;
		}

		iy[0] = (iPi * (i4 * iu_km1[0] + iu_km1[1])) / i6;

		string sinpihstr = boost::lexical_cast<string>(std::sin(M_PIl * h));
		iSin.a = LeftRead<long double>(sinpihstr);
		iSin.b = RightRead<long double>(sinpihstr);

		iy[0] = iy[0] + (id * iSin);

		string alphakm1str = boost::lexical_cast<string>(alpha_km1[0]);
		ialpha.a = LeftRead<long double>(alphakm1str);
		ialpha.b = RightRead<long double>(alphakm1str);

		iy[0] = (iy[0] + (ic * ialpha)) / ia;
		iz[0] = ib / ia;

		for (int i = 1; i < n - 2; ++i) {
			iz[i] = ib / (ia - (ib * iz[i - 1]));
		}

		for (i = 0; i < n - 2; ++i) {
			iy[i + 1] = iu_km1[i] + (i4 * iu_km1[i + 1]);
			if (i < n - 2 - 1) {
				iy[i + 1] = iy[i + 1] + iu_km1[i + 2];
			}

			sinpihstr = boost::lexical_cast<string>(
					std::sin((i + 1.0 + 1.0) * M_PIl * h));
			iSin.a = LeftRead<long double>(sinpihstr);
			iSin.b = RightRead<long double>(sinpihstr);

			iy[i + 1] = ((iPi * iy[i + 1]) / i6) + (id * iSin);

			string alphakm1str = boost::lexical_cast<string>(alpha_km1[i + 1]);
			ialpha.a = LeftRead<long double>(alphakm1str);
			ialpha.b = RightRead<long double>(alphakm1str);
			ialpha = im11 * ialpha;

			iy[i + 1] = iy[i + 1] + (ic * ialpha);
			iy[i + 1] = (iy[i + 1] - (ib * iy[i])) / (ia - (ib * iz[i]));
		}

		iu_k[n - 1 - 1] = iy[n - 1 - 1];

		for (int i = n - 2 - 1; i >= 0; --i) {
			iu_k[i] = iy[i] - iz[i] * iu_k[i + 1];
		}

		for (int i = 0; i < n - 1; ++i) {
			iu_k[i].IEndsToStrings(left, right);
			exact = (1.0 / M_PIl) * std::sin(M_PIl * (i + 1) * h);
			if (output == 's') {
				cout << std::setprecision(2) << "u(" << (i + 1) * h
						<< std::setprecision(7) << ") = [" << left << "; "
						<< right << "]" << endl;
				cout << "  exact = " << exact << endl;
				cout << endl;
			} else {
				results << "u(" << std::setprecision(2) << (i + 1) * h
						<< std::setprecision(7) << ") = [" << left << ", "
						<< right << "]" << endl;
				results << "  exact = " << exact << endl;
				results << endl;
			}
		}

		if (output == 's') {
			cout << endl;
		} else {
			results << endl;
		}

		for (int i = 0; i < n - 1; ++i) {
			beta = std::sin(i * M_PIl * h);
//			cout << "beta = " << beta << endl;

			string betastr = boost::lexical_cast<string>(beta);
			ibeta.a = LeftRead<long double>(betastr);
			ibeta.b = RightRead<long double>(betastr);

			ibeta = (ib1 * iu_km1[i]) * ibeta;

			ia1.IEndsToStrings(left, right);
//												cout << "ia1 = "<< std::setprecision(7) << "[" << left << "; "
//																		<< right << "]" << endl;

			ibeta = ia1 * (iu_km1[i] * iu_km1[i]) + ibeta;

			ibeta.IEndsToStrings(left, right);
//									cout << "ibeta = "<< std::setprecision(7) << "[" << left << "; "
//															<< right << "]" << endl;

			ibeta = ibeta + ic1;
			ibeta = ISqrt(ibeta, error);

			if (error == 0) {
				string id1str = boost::lexical_cast<string>(
						std::sin(M_PIl * h) * alpha_km1[i]);
				id1.a = LeftRead<long double>(id1str);
				id1.b = RightRead<long double>(id1str);
				ibeta = (ibeta + id1) * ih;
				alpha_k[i] = ibeta.b;
			} else {
				cout << "Impossible to calculate Sqrt(beta) for i = " << i
						<< endl;
				break;
			}
		}

		for (int i = 0; i < n - 1; ++i) {
			if (output == 's') {
				cout << "alpha(" << (i + 1) * h << ") = " << alpha_k[i] << endl;
			} else {
				results << "alpha(" << (i + 1) * h << ") = " << alpha_k[i]
						<< endl;
			}
		}

		if (output == 's') {
			cout << endl;
		} else
			results << endl;

		norm_u = 0.0;
		for (int i = 0; i < n - 1; ++i) {
			norm_ui = std::abs(iu_k[i].a - iu_km1[i].a);
			if (std::abs(iu_k[i].b - iu_km1[i].b) > norm_u) {
				norm_ui = std::abs(iu_k[i].b - iu_km1[i].b);
			}
			if (norm_ui > norm_u)
				norm_u = norm_ui;
		}
		if (norm_u < epsilon)
			u_OK = true;
		else
			u_OK = false;

		abs_alpha = 0.0;
		for (int i = 0; i < n - 1; ++i) {
			abs_alphai = std::abs(alpha_k[i] - alpha_km1[i]);
			if (abs_alphai > abs_alpha) {
				abs_alpha = abs_alphai;
			}
		}
		if (abs_alpha < epsilon)
			alpha_OK = true;
		else
			alpha_OK = false;

		if (u_OK && alpha_OK) {
			finish = true;

			if (output == 's') {
				cout << "Press Enter to continue ..." << endl;
				std::getchar();
				cout << "Interval solution after delta-extension:" << endl
						<< endl;
			} else {
				results << "Interval solution after delta-extension:" << endl
						<< endl;
			}

			for (int i = 0; i < n - 1; ++i) {
				iu_k[i].a = iu_k[i].a - delta;
				iu_k[i].b = iu_k[i].b + delta;
				iu_k[i].IEndsToStrings(left, right);

				exact = (1.0 / M_PIl) * std::sin(M_PIl * (i + 1) * h);

				if (output == 's') {
					cout << "u(" << std::setprecision(2) << (i + 1) * h
							<< std::setprecision(7) << ") = [" << left << "; "
							<< right << "]" << endl;
					cout << "     = [" << iu_k[i].a << "; " << iu_k[i].b << "]"
							<< endl;
					cout << "     exact = " << exact << endl;
				} else {
					results << "u(" << (i + 1) * h << ") = [" << left << ', '
							<< right << "]" << endl;
					results << "     = [" << iu_k[i].a << ', ' << iu_k[i].b
							<< "]" << endl;
					results << "     exact = " << exact << endl;
				}
			}
		}
		if (!finish) {
			for (int i = 0; i < n - 1; ++i) {
				iu_km1[i] = iu_k[i];
				alpha_km1[i] = alpha_k[i];
			}

			if (output == 's') {
				cout << "Press Enter to continue ..." << endl;
				std::getchar();
			}

		}
	} //end while (!finish)

	delete[] u;
	delete[] alpha_k;
	delete[] alpha_km1;
	delete[] iu_k;
	delete[] iu_km1;
	delete[] iy;
	delete[] iz;

	if (output == 'f')
		results.close();

	cout << endl;
	cout << "END OF PROGRAM" << endl;
	std::getchar();
}

