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
	this->maxP = 0;
	this->maxQ = 0;
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
	int imod, j, jmod, ui, uj; // i,l, nmin, nmax; // mmax;
	long double exact, h, k, w;
//	bool OK;
	long double sol;
//	char z;
	string file_name, left, right, time;
	fstream filestr;
//	clock_t time1;
//	double time_diff;
//	T Mconst = 0;

	if (!_initparams)
		throw runtime_error("Parameters not initialized!");

	int st = 0;
	int m = params.m;
	int n = params.n;
	long double alpha = params.alpha;
	long double beta = params.beta;
	long double delta = params.delta;
	long double gamma = params.gamma;
//	long double eps = params.eps;

	int dprec = Interval<T>::GetOutDigits();

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
//	l = 0;
	imod = n / 20;
	jmod = m / 20;
	ui = n / 2;
	uj = m / 2;

	filestr << endl;

	for (int i = 0; i <= n; i++) {
		for (j = 0; j <= m; j++) {
			if ((i % imod == 0) && (j % jmod == 0)) {
				exact = bc->ExactSol(alpha + i * h, beta + j * k); // exact solution
				if ((j != 0) && (j != m))
					sol = u[i][j];
				else {
					if (j == 0)
						sol = bc->phi2(alpha + i * h);
					else
						sol = bc->phi4(alpha + i * h);
				}
				w = abs(sol - exact);
				filestr << " " << endl;
				filestr << " u(" << std::setprecision(2) << alpha + i * h << ","
						<< beta + j * k << ") = ";
				filestr << std::setprecision(dprec) << sol << endl;
				filestr << "      err =  " << w << endl;
				filestr << "eu(" << std::setprecision(2) << alpha + i * h << ","
						<< beta + j * k << ") =" << std::setprecision(dprec)
						<< exact << endl;
			}
		}
	}
//	for (int i = 0; i <= n; i++) {
//		if (i % imod == 0) {
//			exact = bc->ExactSol(alpha + i * h, beta + uj * k); // exact solution
//			if ((i != 0) && (i != n))
//				sol = u[i][uj];
//			else {
//				if (i == 0)
//					sol = bc->phi1(beta + uj * k);
//				else
//					sol = bc->phi3(beta + uj * k);
//			}
//			w = abs(sol - exact);
//			filestr << " " << endl;
//			filestr << " u(" << std::setprecision(2) << alpha + i * h << ","
//					<< beta + uj * k << ") = ";
//			filestr << std::setprecision(dprec) << sol << endl;
//			filestr << "      err =  " << w << endl;
//			filestr << "eu(" << std::setprecision(2) << alpha + i * h << ","
//					<< beta + uj * k << ") =" << std::setprecision(dprec)
//					<< exact << endl;
//		}
//	}
	stringstream ss;
//	int prec = Interval<T>::GetOutDigits();
	ss << file_name << "_m_" << m;
}

template<typename T>
void Solver<T>::WriteIntervalResultsToFile() {
	int imod, j, jmod, ui, uj, i; //i,l;
	long double exact, h, k;
	T w;
	bool dint_mode;
	Interval<T> HH, KK, NN, MM, sol, UIH, UII, UJJ, UJK;
//	char z;// z1;
	string left, right;// time;
	fstream filestr;
//	clock_t time1, time2;
//	double time_diff;
	int dprec = Interval<T>::GetOutDigits();
	std::setprecision(dprec);

	if (!_initparams)
		throw runtime_error("Parameters not initialized!");

	int st = 0;
	int n = params.n;
	int m = params.m;
	long double alpha = params.alpha;
	long double beta = params.beta;
	long double delta = params.delta;
	long double gamma = params.gamma;
	Interval<T> intalpha = { params.alpha, params.alpha };
	Interval<T> intbeta = { params.beta, params.beta };
	Interval<T> GAMMA = { params.gamma, params.gamma };
	Interval<T> DELTA = { params.delta, params.delta };
//	T eps = params.eps;
	dint_mode = (params.ia_mode == IAMode::DINT_MODE);
	NN.a = n;
	NN.b = n;
	MM.a = m;
	MM.b = m;

	filestr.open(params.file_name.c_str(), fstream::out);
//	OK = true;
	filestr
			<< "SOLVING THE GENERALIZED POISSON EQUATION WITH GIVEN BOUNDARY CONDITIONS"
			<< endl;
	filestr << "BY AN INTERVAL DIFFERENCE METHOD WITH "
			<< (dint_mode ? "DIRECTED" : "PROPER") << " INTERVAL ARITHMETIC"
			<< endl;

//	time_diff = time2 - time1;
//	cout << "status = " << st << ", time = " << time_diff;
	filestr << " " << endl;
//	filestr << "status = " << st << ", time = " << time << " [s]" << endl;
	filestr << " u - a solution obtained by interval method" << endl;
	filestr << " eu - the exact (approximate) solution" << endl;

	filestr.setf(std::ios::scientific);
	if (st != 0)
		return;

	h = (delta - alpha) / n;
	k = (gamma - beta) / m;
	HH = (DELTA - intalpha) / NN;
	KK = (GAMMA - intbeta) / MM;
//	l = 0;
	imod = n / 10;
	jmod = m / 10;
//	ui = n / 4;
//	uj = m / 4;
//	UII.a = ui;
//	UII.b = ui;
//	UJJ.a = uj;
//	UJJ.b = uj;
	filestr << endl;

	for (i = 0; i <= n; i++) {
		for (j = 0; j <= m; j++) {
				ui = i;
				uj = j;
				UII.a = ui;
				UII.b = ui;
				UJJ.a = uj;
				UJJ.b = uj;
		//		if (j % jmod == 0) {
					exact = bc->ExactSol(alpha + i * h, beta + j * k); // exact solution
					if (((i != 0) && (i != n))&&((j != 0) && (j != m)))
						sol = X[(i - 1) * (m - 1) + uj - 1];
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
					filestr << std::setprecision(2) << " u(" << alpha + i * h << ","
							<< beta + j * k << ") = ";
					filestr << "[" << left << ", " << right << "]" << endl;
					filestr << "      width =  " << std::setprecision(dprec) << w
							<< endl;
					filestr << std::setprecision(2) << "eu(" << alpha + i * h << ","
							<< beta + j * k << ") = " << std::setprecision(dprec)
							<< exact << endl;
		//		}
			}
		}

//	for (j = 0; j <= m; j++) {
////		if (j % jmod == 0) {
//			exact = bc->ExactSol(alpha + ui * h, beta + j * k); // exact solution
//			if ((j != 0) && (j != m))
//				sol = X[(ui - 1) * (m - 1) + j - 1];
//			else {
//				UIH = intalpha + (UII * HH);
//				if (j == 0)
//					sol = bc->PHI2(UIH, st);
//				else
//					sol = bc->PHI4(UIH, st);
//			}
//			sol = sol.Projection();
//			w = sol.GetWidth();
//			sol.IEndsToStrings(left, right);
//			filestr << " " << endl;
//			filestr << std::setprecision(2) << " u(" << alpha + ui * h << ","
//					<< beta + j * k << ") = ";
//			filestr << "[" << left << ", " << right << "]" << endl;
//			filestr << "      width =  " << std::setprecision(dprec) << w
//					<< endl;
//			filestr << std::setprecision(2) << "eu(" << alpha + ui * h << ","
//					<< beta + j * k << ") = " << std::setprecision(dprec)
//					<< exact << endl;
////		}
//	}
//	for (int i = 0; i <= n; i++) {
////		if (i % imod == 0) {
//			exact = bc->ExactSol(alpha + i * h, beta + uj * k); // exact solution
//			if ((i != 0) && (i != n))
//				sol = X[(i - 1) * (m - 1) + uj - 1];
//			else {
//				UJK = intbeta + (UJJ * KK);
//				if (i == 0)
//					sol = bc->PHI1(UJK, st);
//				else
//					sol = bc->PHI3(UJK, st);
//			}
//			sol = sol.Projection();
//			w = sol.GetWidth();
//			sol.IEndsToStrings(left, right);
//			filestr << " " << endl;
//			filestr << std::setprecision(2) << " u(" << alpha + i * h << ","
//					<< beta + uj * k << ") = ";
//			filestr << "[" << left << ", " << right << "]" << endl;
//			filestr << "      width =  " << std::setprecision(dprec) << w
//					<< endl;
//			filestr << std::setprecision(2) << "eu(" << alpha + i * h << ","
//					<< beta + uj * k << ") = " << std::setprecision(dprec)
//					<< exact << endl;
////		}
//	}
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
		if (params.print_csv)
		{
			this->WriteFPResultsToCsv();
		}
		break;
	case INTERVAL_EXP:
		this->WriteIntervalResultsToFile();
		if (params.print_csv)
		{
			this->WriteIntervalResultsToCsv();
		}
		else cout << "Print csv is False!" << endl;
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
	for (int i = 0; i < 12; ++i) {
		this->params.m = (i + 1) * 10;
		this->params.n = this->params.m;
		this->SolveFP();
		this->vecConstM.push_back(this->GetMaxM());
		this->vecConstN.push_back(this->GetMaxN());
		this->vecConstP.push_back(this->GetMaxP());
		this->vecConstQ.push_back(this->GetMaxQ());
		cout << "M = " << this->GetMaxM() << endl;
		cout << "N = " << this->GetMaxN() << endl;
	}
	return 0;
}

template<typename T>
void Solver<T>::WriteConstMResults() {
	fstream filestr;
	filestr.open(params.file_name.c_str(), fstream::out);
	filestr << "CONST M EXPERIMENT RESULTS" << endl;
	if (vecConstP.size() != vecConstQ.size())
		return;

	long double constM = 0;
	long double constN = 0;
	long double constP = 0;
	long double constQ = 0;
	for (size_t i = 0; i < vecConstP.size(); ++i) {
		constM = this->vecConstM.at(i);
		constN = this->vecConstN.at(i);
		constP = this->vecConstP.at(i);
		constQ = this->vecConstQ.at(i);
		filestr << constM << ";" << constN << ";" << constP << ";" << constQ << endl;
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
void Solver<T>::SetParameters(Parameters<long double>& p) {
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
long double Solver<T>::GetMaxP() {
	return this->maxP;
}

template<typename T>
long double Solver<T>::GetMaxQ() {
	return this->maxQ;
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
	int j;
	long double exact, h, k;
	T sol;
	string file_name, left, right, time;
	fstream fp_filestr, exact_filestr;

	if (!_initparams)
		throw runtime_error("Parameters not initialized!");

	int st = 0;
	int m = params.m;
	int n = params.n;
	long double alpha = params.alpha;
	long double beta = params.beta;
	long double delta = params.delta;
	long double gamma = params.gamma;
	string sep = ";";
	int dprec = Interval<T>::GetOutDigits();
	std::setprecision(dprec);
	cout.setf(std::ios_base::scientific);
	fs::path p(params.file_name);
	string fname = fs::basename(params.file_name);

	string dir = p.parent_path().string();
	fp_filestr.open((dir + "/" + fname +"_f.csv").c_str(), fstream::out);
	exact_filestr.open((dir + "/" + fname+"_e.csv").c_str(), fstream::out);

	if (st != 0)
		return;

	h = (delta - alpha) / n;
	k = (gamma - beta) / m;

//	imod = n / 10;
//	jmod = m / 10;

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
			sep = (j==m) ? "\n" : ";";
			fp_filestr << std::setprecision(dprec) << sol << sep;
			exact_filestr << std::setprecision(dprec) << exact << sep;
		}
	}
	fp_filestr.close();
	exact_filestr.close();
}

template<typename T>
inline void Solver<T>::WriteIntervalResultsToCsv() {
	string file_name, left, right, time;
	fstream l_filestr, r_filestr, w_filestr;

	if (!_initparams)
		throw runtime_error("Parameters not initialized!");

	int st = 0;
	int m = params.m;
	int n = params.n;

	fs::path p(params.file_name);
	string fname = fs::basename(params.file_name);
	string dir = p.parent_path().string();
	if (dir.length() > 0)
		{
		  dir = dir +"/";
		}
	string sep = ";";
	T w = 0.0;

	Interval<T> HH, KK, NN, MM, sol, UIH, UII, UJJ, UJK;

	int dprec = Interval<T>::GetOutDigits();
	std::setprecision(dprec);

	Interval<T> intalpha = { params.alpha, params.alpha };
	Interval<T> intbeta = { params.beta, params.beta };
	Interval<T> GAMMA = { params.gamma, params.gamma };
	Interval<T> DELTA = { params.delta, params.delta };

	NN.a = n;
	NN.b = n;
	MM.a = m;
	MM.b = m;

	std::setprecision(dprec);
	cout.setf(std::ios_base::scientific);
	w_filestr.open((dir + fname+"_w.csv").c_str(), fstream::out);
	l_filestr.open((dir + fname+"_l.csv").c_str(), fstream::out);
	r_filestr.open((dir + fname+"_r.csv").c_str(), fstream::out);

	if (st != 0)
		return;

	HH = intalpha / NN;
	KK = intbeta / MM;

	for (int i = 0; i <= n; i++) {
		for (int j = 0; j <= m; j++) {
			if ((j != 0) && (j != m) && (i != 0) && (i != n))
				sol = X[(i - 1) * (m - 1) + j - 1];
			else {
				if (j == 0)
				{
					Interval<T> II = {(long double)i, (long double)i};
					sol = bc->PHI2(intalpha + II * HH, st);
				}
				else if (j == m)
				{
					Interval<T> II = {(long double) i, (long double)i};
					sol = bc->PHI4(intalpha + II * HH, st);
				}
				else if (i == 0)
				{
					Interval<T> JJ = {(long double)j, (long double)j};
					sol = bc->PHI1(intbeta + JJ * KK, st);
				}
				else
				{
					Interval<T> JJ = {(long double)j, (long double)j};
					sol = bc->PHI3(intbeta + JJ * KK, st);
				}
			}
			sol = sol.Projection();
			w = sol.GetWidth();

			sol.IEndsToStrings(left, right);
			sep = (j==m) ? "\n" : ";";
			l_filestr << std::setprecision(dprec) << left << sep;
			r_filestr << std::setprecision(dprec) << right << sep;
			w_filestr << std::setprecision(dprec) << w << sep;

		}
	}
	l_filestr.close();
	r_filestr.close();
	w_filestr.close();
}

//The explicit instantiation part
template class Solver<long double> ;
template class Solver<mpreal> ;

}
/* namespace intervalarth */
