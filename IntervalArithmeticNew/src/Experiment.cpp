/*
 * Experiment.cpp
 *
 *  Created on: Jan 5, 2013
 *      Author: tomaszhof
 */

#include "Experiment.h"

namespace interval_arithmetic
{

template<typename T>
Experiment<T>::Experiment()
{
	// TODO Auto-generated constructor stub
	_param_initialized = false;
	_solver_initialized = false;
}

template<typename T>
Experiment<T>::~Experiment()
{
	// TODO Auto-generated destructor stub
}

template<typename T>
void Experiment<T>::IEndsToStrings(Interval<T> i, string& left, string& right)
{
	stringstream ss;

	fesetround(FE_DOWNWARD);
	int prec = std::numeric_limits<long double>::digits10;
	ss.setf(std::ios_base::scientific);
	ss << std::setprecision(prec) << i.a;
	left = ss.str();
	//left = boost::lexical_cast<string>(i.a);
	ss.str(std::string());

	fesetround(FE_UPWARD);
	ss << std::setprecision(prec) << i.b;
	right = ss.str();
	//right = boost::lexical_cast<string>(i.b);
	ss.clear();
	fesetround(FE_TONEAREST);
}

template<typename T>
string Experiment<T>::FloatToString(long double f)
{
	stringstream ss;
	int prec = std::numeric_limits<long double>::digits10;
	ss.setf(std::ios_base::scientific);
	ss << std::setprecision(prec) << f;
	string res = ss.str();
	ss.clear();
	return res;
}

template<typename T>
void Experiment<T>::SetExample(int eid, int arth_mode)
{
	switch (eid)
	{
	case 1:
		//_example = new Example01();
		break;
	case 2:
		//_example = new Example02();
		break;
	case 3:
		_example = new Example03<T>();
		break;
	case 4:
		_example = new Example04<T>();
		break;
	default:
		_example = NULL;
		break;
	}

//	if (_example != NULL)
//		_example->SetArithmeticMode(arth_mode);
}

template<typename T>
void Experiment<T>::DoExperiment(int eid)
{
	int i, imod, j, jmod, l, n, m, st, ui, uj;
	T alpha, beta, exact, h, k, w;
	bool OK, OK1, dint_mode;
	Interval<T> HH, KK, NN, MM, sol, UIH, UII, UJJ, UJK;
	char z, z1;
	string file_name, left, right, time;
	fstream filestr;
	clock_t time1, time2;
	double time_diff;
	int dprec = std::numeric_limits<long double>::digits10;
	std::setprecision(dprec);

	alpha = 1;
	beta = 1;
	Interval<T> intalpha =
	{ 1, 1 };
	Interval<T> intbeta =
	{ 1, 1 };
	Interval<T> intgamma =
	{ 2, 2 };
	Interval<T> intdelta =
	{ 2, 2 };
	alpha = 1;
	beta = 1;
	Interval<T>* X;

	//int prec = eps.getPrecision();

	do
	{
		cout << "Interval arithmetic mode (d - directed, p - proper): ";
		cin >> z;
	} while ((z != 'd') && (z != 'D') && (z != 'p') && (z != 'P'));
	dint_mode = ((z == 'd') || (z == 'D'));

	cout
			<< "SOLVING THE GENERALIZED POISSON EQUATION WITH GIVEN BOUNDARY CONDITIONS"
			<< endl;
	cout << "BY AN INTERVAL DIFFERENCE METHOD WITH "
			<< (dint_mode ? "DIRECTED" : "PROPER") << " INTERVAL ARITHMETIC"
			<< endl << endl;
	cout << "Insert data" << endl;

	do
	{
		OK = false;
		cout << "n = m = ";
		cin >> n;
		if ((n >= 20) && (n % 10 == 0))
		{
			OK = true;
		}
		else
			cout
					<< "Error: the value of n = m must be an integer number equal to 20, 30, 40, ... ."
					<< endl;
	} while (!OK);

	m = n;
	NN.a = n;
	NN.b = n;
	MM.a = m;
	MM.b = m;

	do
	{
		cout << "Output of results (f - file, s - screen): ";
		cin >> z;
	} while ((z != 'f') && (z != 'F') && (z != 's') && (z != 'S'));

	if ((z == 'f') || (z == 'F'))
	{
		cout << endl;
		do
		{
			OK = false;
			cout << "Insert the name of file" << endl;
			cout << "the extension .txt will be added automatically): ";
			cin >> file_name;
			file_name += ".txt";
			OK1 = false;

			ifstream ifile(file_name.c_str());
			if (ifile)
			{
				ifile.close();
				cout << "The file " + file_name + " already exists. " << endl;
				do
				{
					cout
							<< "Do you want to delete this file and create a new one in its place?"
							<< endl;
					cout << "(y - yes, n - no): ";
					cin >> z1;
				} while ((z1 != 'n') && (z1 != 'N') && (z1 != 'y')
						&& (z1 != 'Y'));

				if ((z1 == 'y') || (z1 == 'Y'))
					OK1 = true;
			}
			else
				OK1 = true;
			if (OK1)
			{
				filestr.open(file_name.c_str(), fstream::out);
				OK = true;
				filestr
						<< "SOLVING THE GENERALIZED POISSON EQUATION WITH GIVEN BOUNDARY CONDITIONS"
						<< endl;
				filestr << "BY AN INTERVAL DIFFERENCE METHOD WITH "
						<< (dint_mode ? "DIRECTED" : "PROPER")
						<< " INTERVAL ARITHMETIC" << endl;
			}
		} while (!OK);

		cout << endl;
		cout << "Calculating... Please, wait!" << endl;
		cout << endl;
		int mode = dint_mode ? DINT_MODE : PINT_MODE;
		this->SetExample(eid, mode);
		//DiffPoisson* dp = new DiffPoisson();
		//dp->SetBoundaryConditions(*_example);
		//IntervalArithmetic* ia = new IntervalArithmetic();

		const long double Mconst = _example->GetConstM(); //1627;
		const long double Nconst = _example->GetConstN();
		const long double eps = 1.0e-16;

		//dp->SetArithmeticMode(mode);

		time1 = clock();
		if (dint_mode)
			;//dp->DIntGPDEbyAM(n, m, intalpha, intbeta, intgamma, intdelta, eps,
			//		Mconst, X, st);
		//dp->DIntPoisson(n, m, intalpha, intbeta, eps, Mconst, X, st);
		else
			;//dp->IntPoisson(n, m, intalpha, intbeta, eps, Mconst, X, st);
			//dp->IntGPDE3(n, m, intalpha, intbeta, intgamma, intdelta, eps,
			//		Mconst, X, st);
		time2 = clock();
		time_diff = ((double) time2 - time1) / CLOCKS_PER_SEC;
		time = boost::lexical_cast<string>(time_diff);
		cout << "status = " << st << ", time = " << time_diff;
		if ((z == 'f') || (z == 'F'))
		{
			filestr << " " << endl;
			filestr << "status = " << st << ", time = " << time << " [s]"
					<< endl;
		}
		else
			return;

		filestr << " u - a solution obtained by interval method" << endl;
		filestr << " eu - the exact (approximate) solution" << endl;

		filestr.setf(std::ios::scientific);
		if (st != 0)
			return;

		h = alpha / n;
		k = beta / m;
		HH = intalpha / NN;
		KK = intbeta / MM;
		l = 0;
		imod = n / 10;
		jmod = m / 10;
		ui = n / 2;
		uj = m / 2;
		UII.a = ui;
		UII.b = ui;
		UJJ.a = uj;
		UJJ.b = uj;
		filestr << endl;
		for (j = 0; j <= m; j++)
		{
			if (j % jmod == 0)
			{
				exact = _example->ExactSol(alpha + ui * h, beta + j * k); //exp(3 * ui * h) * cos(3 * j * k); // exact solution
				if ((j != 0) && (j != m))
					sol = X[(ui - 1) * (m - 1) + j - 1];
				else
				{
					UIH = intalpha + (UII * HH);
					if (j == 0)
						sol = _example->PHI2(UIH, st);
					else
						sol = _example->PHI4(UIH, st);
				}
				sol = sol.Projection();
				w = sol.GetWidth();
				sol.IEndsToStrings(left, right);
				filestr << " " << endl;
				filestr << std::setprecision(2) << " u(" << alpha + ui * h
						<< "," << beta + j * k << ") = ";
				//fprintf(filestr, "[%e , %e] \n", left, right)
				filestr << "[" << left << "," << right << "]" << endl;
				filestr << "      width =  " << std::setprecision(dprec) << w
						<< endl;
				filestr << std::setprecision(2) << "eu(" << alpha + ui * h
						<< "," << beta + j * k << ") ="
						<< std::setprecision(dprec) << exact << endl;
				//fprintf(filestr, "u(")
			}
		}
		for (int i = 0; i <= n; i++)
		{
			if (i % imod == 0)
			{
				exact = _example->ExactSol(alpha + i * h, beta + uj * k); //exp(3 * i * h) * cos(3 * uj * k); // exact solution
				if ((i != 0) && (i != n))
					sol = X[(i - 1) * (m - 1) + uj - 1];
				else
				{
					UJK = intbeta + (UJJ * KK);
					if (i == 0)
						sol = _example->PHI1(UJK, st);
					else
						sol = _example->PHI3(UJK, st);
				}
				sol = sol.Projection();
				w = sol.GetWidth();
				sol.IEndsToStrings(left, right);
				filestr << " " << endl;
				filestr << std::setprecision(2) << " u(" << alpha + i * h << ","
						<< beta + uj * k << ") = ";
				filestr << "[" << left << "," << right << "]" << endl;
				filestr << "      width =  " << std::setprecision(dprec) << w
						<< endl;
				filestr << std::setprecision(2) << "eu(" << alpha + i * h << ","
						<< beta + uj * k << ") =" << std::setprecision(dprec)
						<< exact << endl;
			}
		}
		filestr.close();
	}
}

template<typename T>
void Experiment<T>::DoClassicallExperiment(int eid)
{
	int i, imod, j, jmod, l, n, m, st, ui, uj, nmin, nmax, mmax, step;
	long double alpha, beta, gamma, delta, exact, h, k, w;
	bool OK, OK1;
	//interval HH, KK, NN, MM, sol, UIH, UII, UJJ, UJK;
	long double sol;
	char z, z1;
	string file_name, left, right, time;
	fstream filestr;
	clock_t time1, time2;
	double time_diff;
	alpha = 1;
	beta = 1;
	gamma = 2;
	delta = 2;
	long double Mconst = 0;
	const long double eps = 1.0e-16; // LDBL_EPSILON; //1.0e-16;
	int dprec = std::numeric_limits<long double>::digits10;
	std::setprecision(dprec);
	cout.setf(std::ios_base::scientific);

	long double** u;

	this->SetExample(eid, PINT_MODE);
	//DiffPoisson* dp = new DiffPoisson();
	//dp->SetBoundaryConditions(*_example);
	cout << "SOLVING THE POISSON EQUATION WITH GIVEN BOUNDARY CONDITIONS"
			<< endl;
	cout << "BY AN INTERVAL DIFFERENCE METHOD WITH FLOATING-POINT ARITHMETIC"
			<< endl << endl;
	cout << "Insert data" << endl;

	do
	{
		OK = false;
		cout << "n = m = ";
		cin >> n;

		if ((n >= 20) && (n % 10 == 0))
		{
			OK = true;
		}
		else
			cout
					<< "Error: the value of n = m must be an integer number equal to 20, 30, 40, ... ."
					<< endl;
	} while (!OK);

	m = n;
	//NN.a = n;
	//NN.b = n;
	//MM.a = m;
	//MM.b = m;

	do
	{
		cout << "Output of results (f - file, s - screen): ";
		cin >> z;
	} while ((z != 'f') && (z != 'F') && (z != 's') && (z != 'S'));

	if ((z == 'f') || (z == 'F'))
	{
		cout << endl;
		do
		{
			OK = false;
			cout << "Insert the name of file" << endl;
			cout << "the extension .txt will be added automatically): ";
			cin >> file_name;
			file_name += ".txt";
			OK1 = false;

			ifstream ifile(file_name.c_str());
			if (ifile)
			{
				ifile.close();
				cout << "The file " + file_name + " already exists. " << endl;
				do
				{
					cout
							<< "Do you want to delete this file and create a new one in its place?"
							<< endl;
					cout << "(y - yes, n - no): ";
					cin >> z1;
				} while ((z1 != 'n') && (z1 != 'N') && (z1 != 'y')
						&& (z1 != 'Y'));

				if ((z1 == 'y') || (z1 == 'Y'))
					OK1 = true;
			}
			else
				OK1 = true;
			if (OK1)
			{
				filestr.open(file_name.c_str(), fstream::out);
				filestr.setf(std::ios_base::scientific);
				OK = true;
				filestr
						<< "SOLVING THE POISSON EQUATION WITH GIVEN BOUNDARY CONDITIONS"
						<< endl;
				filestr
						<< "BY AN INTERVAL DIFFERENCE METHOD WITH FLOATING-POINT ARITHMETIC"
						<< endl;
			}
		} while (!OK);

		cout << endl;
		cout << "Calculating... Please, wait!" << endl;
		cout << endl;

		m = n;
		u = new long double*[n + 1];
		for (int i = 0; i <= n; i++)
			u[i] = new long double[m + 1];
		time1 = clock();
		//dp->solveClassicallPoisson2(n, m, alpha, beta, gamma, delta, eps, u,
		//		st);
		time2 = clock();
		time_diff = ((double) time2 - time1) / CLOCKS_PER_SEC;
		time = boost::lexical_cast<string>(time_diff);
		cout << "m=n=" << m;
		cout << "; status = " << st << "; time = " << time << " [s]"
				<< "; Mconst = " << boost::lexical_cast<string>(Mconst) << endl;

		filestr << "m=n=" << m;
		filestr << "; status = " << st << "; time = " << time << " [s]"
				<< "; Mconst = " << boost::lexical_cast<string>(Mconst) << endl;
		filestr << " u - a solution obtained by interval method" << endl;
		filestr << " eu - the exact (approximate) solution" << endl;

		if (st != 0)
			return;

		h = alpha / n;
		k = beta / m;
		l = 0;
		imod = n / 10;
		jmod = m / 10;
		ui = n / 2;
		uj = m / 2;

		filestr << endl;

		for (j = 0; j <= m; j++)
		{
			if (j % jmod == 0)
			{
				exact = _example->ExactSol(alpha + ui * h, beta + j * k);//exp(3 * ui * h) * cos(3 * j * k); // exact solution
				if ((j != 0) && (j != m))
					sol = u[ui][j];
				else
				{
					//UIH = ia->IMul(UII, HH);
					if (j == 0)
						sol = _example->phi2(alpha + ui * h);
					else
						sol = _example->phi4(alpha + ui * h);
				}
				//sol = ia->Projection(sol);
				w = abs(sol - exact); //ia->DIntWidth(sol);
				//ia->IEndsToStrings(sol, left, right);
				filestr << " " << endl;
				filestr << " u(" << std::setprecision(2) << alpha + ui * h
						<< "," << beta + j * k << ") = ";
				filestr << std::setprecision(dprec) << sol << endl;
				filestr << "      err =  " << w << endl;
				filestr << "eu(" << std::setprecision(2) << alpha + ui * h
						<< "," << beta + j * k << ") ="
						<< std::setprecision(dprec) << exact << endl;
			}
		}
		for (int i = 0; i <= n; i++)
		{
			if (i % imod == 0)
			{
				exact = _example->ExactSol(alpha + i * h, beta + uj * k); //exp(3 * i * h) * cos(3 * uj * k); // exact solution
				if ((i != 0) && (i != n))
					sol = u[i][uj];
				else
				{
					//UJK = ia->IMul(UJJ, KK);
					if (i == 0)
						sol = _example->phi1(beta + uj * k);
					else
						sol = _example->phi3(beta + uj * k);
				}
				//sol = ia->Projection(sol);
				w = abs(sol - exact); //ia->DIntWidth(sol);
				//IEndsToStrings(sol, left, right);
				filestr << " " << endl;
				filestr << " u(" << std::setprecision(2) << alpha + i * h << ","
						<< beta + uj * k << ") = ";
				filestr << std::setprecision(dprec) << sol << endl;
				filestr << "      err =  " << w << endl;
				filestr << "eu(" << std::setprecision(2) << alpha + i * h << ","
						<< beta + uj * k << ") =" << std::setprecision(dprec)
						<< exact << endl;
			}
		}
		//filestr.close();
		stringstream ss;
		int prec = std::numeric_limits<long double>::digits10;
		//ss.setf( std::ios_base::scientific);
		ss << file_name << "_m_" << m;

		fstream datafilestr;
		datafilestr.open(ss.str().c_str(), fstream::out);
		datafilestr.setf(std::ios_base::scientific);
		for (int i = 0; i <= n; i++)
		{
			for (int j = 0; j <= m; j++)
			{

				if (i % imod == 0)
				{
					if (j % jmod == 0)
					{
						if ((i != 0) && (i != n))
							sol = u[i][j];
						else
						{
							if (i == 0)
								sol = _example->phi1(beta + j * k);
							else
								sol = _example->phi3(beta + j * k);
						}

						if ((j != 0) && (j != m))
							sol = u[i][j];
						else
						{
							if (j == 0)
								sol = _example->phi2(alpha + i * h);
							else
								sol = _example->phi4(alpha + i * h);
						}
						//datafilestr << "u[" << i*h << " , " << j*k << "]="<< std::setprecision(prec) << sol << endl;
						datafilestr << std::setprecision(dprec) << sol << ";";
					}

				}
			}
			datafilestr << endl;
		}

		datafilestr.close();
		for (int i = 0; i <= n; i++)
			delete[] u[i];
		delete[] u;

		filestr.close();
	}

}

template<typename T>
void Experiment<T>::DoConstMExperiment(int eid)
{
	int i, imod, j, jmod, l, n, m, st, ui, uj, nmin, nmax, mmax, step;
	long double alpha, beta, gamma, delta, exact, h, k, w;
	bool OK, OK1;
	//interval HH, KK, NN, MM, sol, UIH, UII, UJJ, UJK;
	long double sol;
	char z, z1;
	string file_name, left, right, time;
	fstream filestr;
	clock_t time1, time2;
	double time_diff;
	alpha = 1;
	beta = 1;
	gamma = 2;
	delta = 2;
	long double Mconst = 0;
	long double Nconst = 0;
	const long double eps = 1.0e-16; // LDBL_EPSILON; //1.0e-16;
	//intervalvector X;
	long double** u;
	this->SetExample(eid, PINT_MODE);
	//DiffPoisson* dp = new DiffPoisson();
	//dp->SetBoundaryConditions(*_example);
	//int dprec = std::numeric_limits<long double>::digits10;
	//std::setprecision(dprec);

	cout << "TTHH SOLVING THE POISSON EQUATION WITH GIVEN BOUNDARY CONDITIONS"
			<< endl;
	cout << "BY AN INTERVAL DIFFERENCE METHOD WITH DIRECTED INTERVAL ARITHMETIC"
			<< endl << endl;
	cout << "Insert data" << endl;

	do
	{
		OK = false;
		cout << "nmin = mmin = ";
		cin >> nmin;
		cout << "nmax = mmax = ";
		cin >> nmax;
		cout << "step = ";
		cin >> step;
		if ((nmin >= 20) && (nmin % 10 == 0) && (nmax >= nmin)
				&& (nmax % 10 == 0))
		{
			OK = true;
		}
		else
			cout
					<< "Error: the value of n = m must be an integer number equal to 20, 30, 40, ... ."
					<< endl;
	} while (!OK);

	m = n;
	//NN.a = n;
	//NN.b = n;
	//MM.a = m;
	//MM.b = m;

	do
	{
		cout << "Output of results (f - file, s - screen): ";
		cin >> z;
	} while ((z != 'f') && (z != 'F') && (z != 's') && (z != 'S'));

	if ((z == 'f') || (z == 'F'))
	{
		cout << endl;
		do
		{
			OK = false;
			cout << "Insert the name of file" << endl;
			cout << "the extension .txt will be added automatically): ";
			cin >> file_name;
			file_name += ".txt";
			OK1 = false;

			ifstream ifile(file_name.c_str());
			if (ifile)
			{
				ifile.close();
				cout << "The file " + file_name + " already exists. " << endl;
				do
				{
					cout
							<< "Do you want to delete this file and create a new one in its place?"
							<< endl;
					cout << "(y - yes, n - no): ";
					cin >> z1;
				} while ((z1 != 'n') && (z1 != 'N') && (z1 != 'y')
						&& (z1 != 'Y'));

				if ((z1 == 'y') || (z1 == 'Y'))
					OK1 = true;
			}
			else
				OK1 = true;
			if (OK1)
			{
				filestr.open(file_name.c_str(), fstream::out);
				OK = true;
				filestr
						<< "SOLVING THE POISSON EQUATION WITH GIVEN BOUNDARY CONDITIONS"
						<< endl;
				filestr
						<< "BY AN INTERVAL DIFFERENCE METHOD WITH DIRECTED INTERVAL ARITHMETIC"
						<< endl;
			}
		} while (!OK);

		cout << endl;
		cout << "Calculating... Please, wait!" << endl;
		cout << endl;
		for (n = nmin; n <= nmax; n = n + step)
		{
			m = n;
			u = new long double*[n + 1];
			for (int i = 0; i <= n; i++)
				u[i] = new long double[m + 1];
			time1 = clock();
			//dp->diffPoisson(n, m, alpha, beta, eps, u, Mconst, st);
			//dp->solveClassicallPoisson(n, m, alpha, beta, gamma, delta, eps, u,
			//		st, Mconst, Nconst);
			time2 = clock();
			time_diff = ((double) time2 - time1) / CLOCKS_PER_SEC;
			time = boost::lexical_cast<string>(time_diff);
			cout << "m=n=" << m;
			cout << "; status = " << st << "; time = " << time << " [s]"
					<< "; Mconst = " << boost::lexical_cast<string>(Mconst)
					<< "; Nconst = " << boost::lexical_cast<string>(Nconst)
					<< endl;

			filestr << "m=n=" << m;
			filestr << "; status = " << st << "; time = " << time << " [s]"
					<< "; Mconst = " << boost::lexical_cast<string>(Mconst)
					<< "; Nconst = " << boost::lexical_cast<string>(Nconst)
					<< endl;
			filestr << " u - a solution obtained by interval method" << endl;
			filestr << " eu - the exact (approximate) solution" << endl;

			if (st != 0)
				return;

			h = alpha / n;
			k = beta / m;
			l = 0;
			imod = n / 10;
			jmod = m / 10;
			ui = n / 2;
			uj = m / 2;

			filestr << endl;

			for (j = 0; j <= m; j++)
			{
				if (j % jmod == 0)
				{
					exact = exp(3 * ui * h) * cos(3 * j * k); // exact solution
					if ((j != 0) && (j != m))
						sol = u[ui][j];
					else
					{
						if (j == 0)
							;//sol = dp->phi2(ui * h);
						else
							;//sol = dp->phi4(ui * h);
					}
					//sol = ia->Projection(sol);
					w = abs(sol - exact); //ia->DIntWidth(sol);
					//ia->IEndsToStrings(sol, left, right);
					filestr << " " << endl;
					filestr << " u(" << ui * h << "," << j * k << ") = ";
					filestr << boost::lexical_cast<string>(sol) << endl;
					filestr << "      err =  " << boost::lexical_cast<string>(w)
							<< endl;
					filestr << "eu(" << ui * h << "," << j * k << ") ="
							<< boost::lexical_cast<string>(exact) << endl;
				}
			}
			for (int i = 0; i <= n; i++)
			{
				if (i % imod == 0)
				{
					exact = exp(3 * i * h) * cos(3 * uj * k); // exact solution
					if ((i != 0) && (i != n))
						sol = u[i][uj];
					else
					{
						if (i == 0)
							;//sol = dp->phi1(uj * k);
						else
							;//sol = dp->phi3(uj * k);
					}
					//sol = ia->Projection(sol);
					w = abs(sol - exact); //ia->DIntWidth(sol);
					//IEndsToStrings(sol, left, right);
					filestr << " " << endl;
					filestr << " u(" << i * h << "," << uj * k << ") = ";
					filestr << boost::lexical_cast<string>(sol) << endl;
					filestr << "      err =  " << boost::lexical_cast<string>(w)
							<< endl;
					filestr << "eu(" << i * h << "," << uj * k << ") ="
							<< boost::lexical_cast<string>(exact) << endl;
				}
			}
			//filestr.close();
			stringstream ss;
			int prec = std::numeric_limits<long double>::digits10;
			//ss.setf( std::ios_base::scientific);
			ss << file_name << "_m_" << m;

			fstream datafilestr;
			datafilestr.open(ss.str().c_str(), fstream::out);
			datafilestr.setf(std::ios_base::scientific);
			for (int i = 0; i <= n; i++)
			{
				for (int j = 0; j <= m; j++)
				{

					if (i % imod == 0)
					{
						if (j % jmod == 0)
						{
							if ((i != 0) && (i != n))
								sol = u[i][j];
							else
							{
								if (i == 0)
									;//sol = dp->phi1(j * k);
								else
									;//sol = dp->phi3(j * k);
							}

							if ((j != 0) && (j != m))
								sol = u[i][j];
							else
							{
								if (j == 0)
									;//sol = dp->phi2(i * h);
								else
									;//sol = dp->phi4(i * h);
							}
							//datafilestr << "u[" << i*h << " , " << j*k << "]="<< std::setprecision(prec) << sol << endl;
							datafilestr << std::setprecision(prec) << sol
									<< ";";
						}

					}
				}
				datafilestr << endl;
			}

			datafilestr.close();
			for (int i = 0; i <= n; i++)
				delete[] u[i];
			delete[] u;
		}
		filestr.close();
	}

}

template<typename T>
Experiment<T>::Experiment(int ac, char* av[])
{
	_example = NULL;
	solver = NULL;
	_param_initialized = false;
	_solver_initialized = false;
	ConfigReader<T>* cf = new ConfigReader<T>(ac, av);
	this->SetParameters(cf->GetParameters());
	delete cf;
}

template<typename T>
void Experiment<T>::SetSolver(Parameters<T> p)
{
	switch (p.selected_solver)
	{
	case GPDE_SOLVER:
		solver = new GPDESolver<T>();
		break;
	default:
		solver = NULL;
		break;
	}
	_solver_initialized = (solver != NULL);

	if (_solver_initialized)
	{
		solver->SetParameters(parameters);
	}
}

template<typename T>
void Experiment<T>::SetParameters(Parameters<T> p)
{
	this->parameters = p;
	_param_initialized = true;
}

template<typename T>
void Experiment<T>::Initialize()
{
	if (!_param_initialized)
		return;
	this->SetSolver(parameters);
}

template<typename T>
void Experiment<T>::SetExampleForSolver(int eid)
{
	if (_solver_initialized)
		solver->SetExample(eid);
}

template<typename T>
void Experiment<T>::Execute()
{
	if (_solver_initialized)
	{
		solver->Execute();
		solver->WriteResults();
	}
}

//The explicit instantiation part
template class Experiment<long double>;
}

