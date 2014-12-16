/*
 * Solver.cpp
 *
 *  Created on: 25-01-2014
 *      Author: thof
 */

#include "Solver.h"

using namespace interval_arithmetic;
using namespace std;

namespace interval_arithmetic {

template<typename T>
Solver<T>::Solver() {
	this->_initparams = false;
	this->_estimateMN = false;
	this->maxM = 0;
	this->maxN = 0;
	this->X = NULL;
	this->bc = NULL;
	this->u = NULL;
}

template<typename T>
Solver<T>::~Solver() {
	// TODO Auto-generated destructor stub
}

template<typename T>
void Solver<T>::WriteFPResultsToFile() {
	int i, imod, j, jmod, l, ui, uj, nmin, nmax, mmax, step;
	long double exact, h, k, w;
	bool OK, OK1;
	long double sol;
	char z, z1;
	string file_name, left, right, time;
	fstream filestr;
	clock_t time1, time2;
	double time_diff;
	long double Mconst = 0;

	if (!_initparams)
		throw runtime_error("Parameters not initialized!");

	int st = 0;
	int m = params.m;
	int n = params.n;
	T alpha = params.alpha;
	T beta = params.beta;
	T delta = params.delta;
	T gamma = params.gamma;
	T eps = params.eps;

	int dprec = std::numeric_limits<T>::digits10;
	std::setprecision(dprec);
	cout.setf(std::ios_base::scientific);
	filestr.open(params.file_name.c_str(), fstream::out);
	filestr.setf(std::ios_base::scientific);
	filestr << "SOLVING THE POISSON EQUATION WITH GIVEN BOUNDARY CONDITIONS"
			<< endl;
	filestr << "IN FLOATING-POINT ARITHMETIC" << endl;
	filestr << "m=n=" << m;
	filestr << "; status = " << st << "; time = " << time << " [s]" << endl;
	filestr << " u - a solution obtained by interval method" << endl;
	filestr << " eu - the exact (approximate) solution" << endl;

	if (st != 0)
		return;

	h = (delta - alpha) / n;
	k = (gamma - beta) / m;
	l = 0;
	imod = n / 10;
	jmod = m / 10;
	ui = n / 2;
	uj = m / 2;

	filestr << endl;

	for (j = 0; j <= m; j++) {
		if (j % jmod == 0) {
			exact = bc->ExactSol(alpha + ui * h, beta + j * k); // exact solution
			if ((j != 0) && (j != m))
				sol = u[ui][j];
			else {
				if (j == 0)
					sol = bc->phi2(alpha + ui * h);
				else
					sol = bc->phi4(alpha + ui * h);
			}
			w = abs(sol - exact);
			filestr << " " << endl;
			filestr << " u(" << std::setprecision(2) << alpha + ui * h << ","
					<< beta + j * k << ") = ";
			filestr << std::setprecision(dprec) << sol << endl;
			filestr << "      err =  " << w << endl;
			filestr << "eu(" << std::setprecision(2) << alpha + ui * h << ","
					<< beta + j * k << ") =" << std::setprecision(dprec)
					<< exact << endl;
		}
	}
	for (int i = 0; i <= n; i++) {
		if (i % imod == 0) {
			exact = bc->ExactSol(alpha + i * h, beta + uj * k); // exact solution
			if ((i != 0) && (i != n))
				sol = u[i][uj];
			else {
				if (i == 0)
					sol = bc->phi1(beta + uj * k);
				else
					sol = bc->phi3(beta + uj * k);
			}
			w = abs(sol - exact);
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
	stringstream ss;
	int prec = std::numeric_limits<T>::digits10;
	ss << file_name << "_m_" << m;
}

template<typename T>
void Solver<T>::WriteIntervalResultsToFile() {
	int i, imod, j, jmod, l, ui, uj;
	long double exact, h, k, w;
	bool OK, OK1, dint_mode;
	Interval<T> HH, KK, NN, MM, sol, UIH, UII, UJJ, UJK;
	char z, z1;
	string left, right, time;
	fstream filestr;
	clock_t time1, time2;
	double time_diff;
	int dprec = std::numeric_limits<T>::digits10;
	std::setprecision(dprec);

	if (!_initparams)
		throw runtime_error("Parameters not initialized!");

	int st = 0;
	int n = params.n;
	int m = params.m;
	T alpha = params.alpha;
	T beta = params.beta;
	T delta = params.delta;
	T gamma = params.gamma;
	Interval<T> intalpha = { params.alpha, params.alpha };
	Interval<T> intbeta = { params.beta, params.beta };
	Interval<T> GAMMA = { params.gamma, params.gamma };
	Interval<T> DELTA = { params.delta, params.delta };
	T eps = params.eps;
	NN.a = n;
	NN.b = n;
	MM.a = m;
	MM.b = m;

	filestr.open(params.file_name.c_str(), fstream::out);
	OK = true;
	filestr
			<< "SOLVING THE GENERALIZED POISSON EQUATION WITH GIVEN BOUNDARY CONDITIONS"
			<< endl;
	filestr << "BY AN INTERVAL DIFFERENCE METHOD WITH "
			<< (dint_mode ? "DIRECTED" : "PROPER") << " INTERVAL ARITHMETIC"
			<< endl;

	cout << "status = " << st << ", time = " << time_diff;
	filestr << " " << endl;
	filestr << "status = " << st << ", time = " << time << " [s]" << endl;
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
	for (j = 0; j <= m; j++) {
		if (j % jmod == 0) {
			exact = bc->ExactSol(alpha + ui * h, beta + j * k); // exact solution
			if ((j != 0) && (j != m))
				sol = X[(ui - 1) * (m - 1) + j - 1];
			else {
				UIH = intalpha + (UII * HH);
				if (j == 0)
					sol = bc->PHI2(UIH, st);
				else
					sol = bc->PHI4(UIH, st);
			}
			sol = sol.Projection();
			w = sol.GetWidth();
			sol.IEndsToStrings(left, right);
			filestr << " " << endl;
			filestr << std::setprecision(2) << " u(" << alpha + ui * h << ","
					<< beta + j * k << ") = ";
			filestr << "[" << left << "," << right << "]" << endl;
			filestr << "      width =  " << std::setprecision(dprec) << w
					<< endl;
			filestr << std::setprecision(2) << "eu(" << alpha + ui * h << ","
					<< beta + j * k << ") =" << std::setprecision(dprec)
					<< exact << endl;
		}
	}
	for (int i = 0; i <= n; i++) {
		if (i % imod == 0) {
			exact = bc->ExactSol(alpha + i * h, beta + uj * k); // exact solution
			if ((i != 0) && (i != n))
				sol = X[(i - 1) * (m - 1) + uj - 1];
			else {
				UJK = intbeta + (UJJ * KK);
				if (i == 0)
					sol = bc->PHI1(UJK, st);
				else
					sol = bc->PHI3(UJK, st);
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

template<typename T>
void Solver<T>::WriteResults() {
	switch (params.exp_mode) {
	case CONST_M_EXP:
		this->WriteConstMResults();
		break;
	case CLASSICAL_EXP:
		this->WriteFPResultsToFile();
		this->WriteFPResultsToCsv();
		break;
	case INTERVAL_EXP:
		this->WriteIntervalResultsToFile();
		break;
	}
}

template<typename T>
void Solver<T>::InitializeX(int m, int n) {
	Interval<T> izero = { 0, 0 };
	int n3 = (n * m - n - m + 4) * (n * m - n - m + 4) / 4;
	X = new Interval<T> [n3];
	for (int i = 1; i <= n3; i++)
		X[i - 1] = izero;
}

template<typename T>
int Solver<T>::ConstMExperiment() {
	this->SetEstimateMN(true);
	for (int i = 2; i < 6; ++i) {
		this->params.m = (i + 1) * 10;
		this->params.n = this->params.m;
		this->SolveFP();
		this->vecConstM.push_back(this->GetMaxM());
		this->vecConstN.push_back(this->GetMaxN());
	}
	return 0;
}

template<typename T>
void Solver<T>::WriteConstMResults() {
	fstream filestr;
	filestr.open(params.file_name.c_str(), fstream::out);
	filestr << "CONST M EXPERIMENT RESULTS" << endl;
	if (vecConstM.size() != vecConstN.size())
		return;

	long double constM = 0;
	long double constN = 0;
	for (int i = 0; i < vecConstM.size(); ++i) {
		constM = this->vecConstM.at(i);
		constN = this->vecConstN.at(i);
		filestr << constM << ";" << constN << endl;
	}
	filestr.close();
}

template<typename T>
int Solver<T>::SetExample(int eid) {
	return -1;
}

template<typename T>
int Solver<T>::SolveFP() {
	return 0;
}

template<typename T>
int Solver<T>::SolvePIA() {
	return 0;
}

template<typename T>
int Solver<T>::SolveDIA() {
	return 0;
}

template<typename T>
void Solver<T>::SetParameters(Parameters<T>& p) {
	this->params = p;
	this->SetExample(p.example_id);
	_initparams = true;
}

template<typename T>
void Solver<T>::Execute() {
	switch (params.exp_mode) {
	case CONST_M_EXP:
		this->ConstMExperiment();
		break;
	case CLASSICAL_EXP:
		this->SolveFP();
		break;
	case INTERVAL_EXP:
		this->SolveInterval();
		break;
	}
}

template<typename T>
int Solver<T>::SolveInterval() {
	int result = -1;
	switch (params.ia_mode) {
	case DINT_MODE:
		result = this->SolveDIA();
		break;
	case PINT_MODE:
		result = this->SolvePIA();
		break;
	default:
		result = this->SolvePIA();
		break;
	}
	return result;
}

template<typename T>
bool Solver<T>::SetEstimateMN(bool b) {
	this->_estimateMN = b;
	return this->_estimateMN;
}

template<typename T>
bool Solver<T>::GetEstimateMN() {
	return this->_estimateMN;
}

template<typename T>
long double Solver<T>::GetMaxM() {
	return this->maxM;
}

template<typename T>
long double Solver<T>::GetMaxN() {
	return this->maxN;
}

template<typename T>
inline void Solver<T>::WriteFPResultsToCsv() {
	int i, imod, j, jmod, l, ui, uj, nmin, nmax, mmax, step;
	long double exact, h, k, w;
	bool OK, OK1;
	long double sol;
	char z, z1;
	string file_name, left, right, time;
	fstream fp_filestr, exact_filestr;
	clock_t time1, time2;
	double time_diff;
	long double Mconst = 0;

	if (!_initparams)
		throw runtime_error("Parameters not initialized!");

	int st = 0;
	int m = params.m;
	int n = params.n;
	T alpha = params.alpha;
	T beta = params.beta;
	T delta = params.delta;
	T gamma = params.gamma;
	T eps = params.eps;

	int dprec = std::numeric_limits<T>::digits10;
	std::setprecision(dprec);
	cout.setf(std::ios_base::scientific);
	fp_filestr.open("res_fp.csv", fstream::out);
	exact_filestr.open("res_exact.csv", fstream::out);

	if (st != 0)
		return;

	h = (delta - alpha) / n;
	k = (gamma - beta) / m;
	l = 0;
	imod = n / 10;
	jmod = m / 10;
	ui = n / 2;
	uj = m / 2;

	fp_filestr << endl;

	for (int i = 0; i <= n; i++) {
		for (j = 0; j <= m; j++) {
			exact = bc->ExactSol(alpha + i * h, beta + j * k); // exact solution
			if ((j != 0) && (j != m) && (i != 0) && (i != n))
				sol = u[i][j];
			else {
				if (j == 0)
					sol = bc->phi2(alpha + i * h);
				else if (j == m)
					sol = bc->phi4(alpha + i * h);
				else if (i == 0)
					sol = bc->phi1(beta + j * k);
				else
					sol = bc->phi3(beta + j * k);
			}
			fp_filestr << std::setprecision(dprec) << sol << ";";
			exact_filestr << std::setprecision(dprec) << exact << ";";
		}
		fp_filestr << std::setprecision(dprec) << endl;
		exact_filestr << endl;
	}
	fp_filestr.close();
	exact_filestr.close();
}

template<typename T>
inline void Solver<T>::WriteIntervalResultsToCsv() {

}

//The explicit instantiation part
template class Solver<long double> ;

}
/* namespace intervalarth */
