/*
 * Config.h
 *
 *  Created on: 12-04-2014
 *      Author: thof
 */

#ifndef CONFIGREADER_H_
#define CONFIGREADER_H_

#include "Utils.h"
#include <string>
#include <mpfr.h>
#include <mpreal.h>


using namespace std;
using namespace mpfr;

namespace interval_arithmetic {

template<typename T>
class ConfigReader {
private:
	IAMode arth_mode; //experiment mode
	ExperimentMode exp_mode; //experiment mode
	Solvers solver_id;
	int m;     //grid size m
	int n;     //grid size n
	int e;     //example id
	T alpha1;
	T alpha2;
	T beta1;
	T beta2;
	string file_name;
	bool print_csv;
	void parseCommandArgs(int ac, char *av[]);
	void readConfigFile(string fileName);
	void readParametersFromConsole();
public:
	ConfigReader();
	ConfigReader(int ac, char *av[]);
	virtual ~ConfigReader();
	ExperimentMode SetExperimentMode(ExperimentMode mode);
	IAMode SetArithmeticMode(IAMode mode);
	Solvers SetSolverId(Solvers sol_id);
	int SetM(int m);
	int SetN(int n);
	T SetAlpha1(T alpha1);
	T SetAlpha2(T alpha2);
	T SetBeta1(T beta1);
	T SetBeta2(T beta2);
	int SetExampleId(int e);
	bool SetPrintCSV(bool b);
	void SetFileName(string fname);
	ExperimentMode GetExperimentMode();
	IAMode GetArithmeticMode();
	Solvers GetSolverId();
	int GetM();
	int GetN();
	int GetExampleId();
	bool GetPrintCSV();
	T GetAlpha1();
	T GetAlpha2();
	T GetBeta1();
	T GetBeta2();
	string GetFileName();
	Parameters<T> GetParameters();
	void SetDefaultParameters();
};

} /* namespace intervalarth */
#endif /* CONFIGREADER_H_ */
