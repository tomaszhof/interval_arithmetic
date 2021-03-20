/*
 * Config.cpp
 *
 *  Created on: 12-04-2014
 *      Author: thof
 */

#include "ConfigReader.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>

namespace po = boost::program_options;
using namespace std;

namespace interval_arithmetic
{

template<typename T>
ConfigReader<T>::ConfigReader()
{
	//default configuration
	this->exp_mode =  INTERVAL_EXP; //CLASSICAL_EXP; //CONST_M_EXP; //INTERVAL_EXP; //CLASSICAL_EXP; //CONST_M_EXP
	this->arth_mode = PINT_MODE;
	this->m = 20; //default grid size m=n=20
	this->n = 20;
	this->e = 4; //default example_id = 1
	this->alpha1 = 1;
	this->alpha2 = 2;
	this->beta1 = 1;
	this->beta2 = 2;
	this->file_name = "ex04pint.txt";
	this->print_csv = false;
	this->solver_id = GPE_SOLVER3C;
}

template<typename T>
ConfigReader<T>::ConfigReader(int ac, char* av[])
{
	//default configuration
	this->exp_mode =  INTERVAL_EXP; //CLASSICAL_EXP; //CONST_M_EXP; //INTERVAL_EXP; //CLASSICAL_EXP; //CONST_M_EXP
	this->arth_mode = PINT_MODE;
	this->m = 20; //default grid size m=n=20
	this->n = 20;
	this->e = 8; //default example_id = 1
	this->alpha1 = 1;
	this->alpha2 = 2;
	this->beta1 = 1;
	this->beta2 = 2;
	this->file_name = "ex04pint.txt";
	this->print_csv = false;
	this->solver_id = GPE_SOLVER3C;
	this->parseCommandArgs(ac, av);
}

template<typename T>
ConfigReader<T>::~ConfigReader()
{
	// TODO Auto-generated destructor stub
}

template<typename T>
ExperimentMode ConfigReader<T>::GetExperimentMode()
{
	return this->exp_mode;
}

template<typename T>
IAMode ConfigReader<T>::GetArithmeticMode()
{
	return this->arth_mode;
}

template<typename T>
int ConfigReader<T>::GetM()
{
	return this->m;
}

template<typename T>
int ConfigReader<T>::GetN()
{
	return this->n;
}



template<typename T>
ExperimentMode ConfigReader<T>::SetExperimentMode(ExperimentMode mode)
{
	this->exp_mode = mode;
	return this->exp_mode;
}

template<typename T>
IAMode ConfigReader<T>::SetArithmeticMode(IAMode mode)
{
	this->arth_mode = mode;
	return this->arth_mode;
}

template<typename T>
int ConfigReader<T>::SetM(int m)
{
	this->m = m;
	return this->m;
}

template<typename T>
int ConfigReader<T>::SetN(int n)
{
	this->n = n;
	return this->n;
}

template<typename T>
int ConfigReader<T>::SetExampleId(int e)
{
	this->e = e;
	return this->e;
}

template<typename T>
int ConfigReader<T>::GetExampleId()
{
	return this->e;
}

template<typename T>
bool ConfigReader<T>::SetPrintCSV(bool b)
{
	this->print_csv = b;
	return this->print_csv;
}

template<typename T>
bool ConfigReader<T>::GetPrintCSV()
{
	return this->print_csv;
}

template<typename T>
void ConfigReader<T>::readConfigFile(string fileName)
{
	if (!boost::filesystem::exists(fileName.c_str()))
	{
		std::cout << "Can't find config file!" << std::endl;
		return;
	}

	//default values of config's parameters
	int m, exId;
	T alpha1, alpha2, beta1, beta2;
	string ia_mode, exp_mode, solver_name, output_file;

	//setup options
	po::options_description desc("Program options");
	desc.add_options()("arithmetic.mode", po::value<string>(&ia_mode),
			"interval arithmetic mode")("experiment.mode",
			po::value<string>(&exp_mode), "experiment mode")("experiment.m",
			po::value<int>(&m), "grid size m")("experiment.exId",
			po::value<int>(&exId), "example ID")("experiment.alpha1",
			po::value<T>(&alpha1))("experiment.alpha2",
			po::value<T>(&alpha2))("experiment.beta1",
			po::value<T>(&beta1))("experiment.beta2",
			po::value<T>(&beta2))("solver.name",
			po::value<string>(&solver_name))("solver.output",
			po::value<string>(&output_file));
	po::variables_map vm;

	std::ifstream settings_file(fileName.c_str());

	//clear the map.
	vm = po::variables_map();

	po::store(po::parse_config_file(settings_file, desc), vm);
	po::notify(vm);

	//setting parameters values
	if (ia_mode == "dint")
		this->SetArithmeticMode(DINT_MODE);
	else
		this->SetArithmeticMode(PINT_MODE);

	if (exp_mode == "interval_exp")
		this->SetExperimentMode(INTERVAL_EXP);
	else if (exp_mode == "const_m_exp")
		this->SetExperimentMode(CONST_M_EXP);
	else
		this->SetExperimentMode(CLASSICAL_EXP);
	this->SetExampleId(exId);

	if (solver_name == "gpe5c")
		this->SetSolverId(GPE_SOLVER5C);
	if (solver_name == "gpe3c")
			this->SetSolverId(GPE_SOLVER3C);


	this->SetFileName(output_file);
	this->SetAlpha1(alpha1);
	this->SetAlpha2(alpha2);
	this->SetBeta1(beta1);
	this->SetBeta2(beta2);
	this->SetM(m);
	this->SetN(m);
}

template<typename T>
T ConfigReader<T>::SetAlpha1(T alpha1)
{
	return this->alpha1 = alpha1;
}

template<typename T>
T ConfigReader<T>::SetAlpha2(T alpha2)
{
	return this->alpha2 = alpha2;
}

template<typename T>
T ConfigReader<T>::SetBeta1(T beta1)
{
	return this->beta1 = beta1;
}

template<typename T>
T ConfigReader<T>::SetBeta2(T beta2)
{
	return this->beta2 = beta2;
}

template<typename T>
T ConfigReader<T>::GetAlpha1()
{
	return this->alpha1;
}

template<typename T>
T ConfigReader<T>::GetAlpha2()
{
	return this->alpha2;
}

template<typename T>
T ConfigReader<T>::GetBeta1()
{
	return this->beta1;
}

template<typename T>
T ConfigReader<T>::GetBeta2()
{
	return this->beta2;
}

template<typename T>
void ConfigReader<T>::parseCommandArgs(int ac, char* av[])
{
	string ia_mode, ex_mode, solver;

	//default options values
	string configFileName, outFileName, programMode;
	int m, e;
	T alpha1, alpha2, beta1, beta2;
	bool print_csv;

	//parsing values
	po::options_description desc("Allowed options");

	desc.add_options()("help", "print help")("m", po::value<int>(),
			"set grid size m")("mode", po::value<string>(), "set program mode")(
			"arth_mode", po::value<string>(),
			"set interval arithmetic mode [pint-proper, dint-directed]")
			("mode", po::value<string>(), "set program mode")(
						"print_csv", po::value<bool>(),
						"print all (warining!) results to csv")(
			"exp_mode", po::value<string>(),
			"set experiment mode [c-classic, m-constM]")
			("solver", po::value<string>(),
						"set solver [gpe5c, gpe3c]")("e", po::value<int>(),
			"set example [e=1,2,3,4]")("file", po::value<string>(),
			"set configuration file name")("out_file", po::value<string>(),
			"set output file name")("alpha1", po::value<T>(),
			"set alpha1")("alpha2", po::value<T>(), "set alpha2")(
			"beta1", po::value<T>(), "set beta1")("beta2",
			po::value<T>(), "set beta2");

	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);

	if (vm.count("help"))
	{
		cout << desc << "\n";
		return;
	}

	if (vm.count("mode"))
	{
		programMode = vm["mode"].as<string>();
		if (programMode == "console")
			this->readParametersFromConsole();
		return;
	}

	if (vm.count("file"))
	{
		configFileName = vm["file"].as<string>();
		this->readConfigFile(configFileName);
		cout << "Configuration file was set to " << configFileName << ".\n";
		return;
	}

	if (vm.count("print_csv"))
		{
			print_csv = vm["print_csv"].as<bool>();
			this->SetPrintCSV(print_csv);
			cout << "Print_csv was set to " << print_csv << ".\n";
		}

	if (vm.count("out_file"))
	{
		outFileName = vm["out_file"].as<string>();
		this->SetFileName(outFileName);
		cout << "Results output file was set to " << outFileName << ".\n";
	}
	else
	{
		cout << "Results output file was not set.\n";
	}

	if (vm.count("arth_mode"))
	{
		ia_mode = vm["arth_mode"].as<string>();
		if (ia_mode == "dint")
			this->SetArithmeticMode(DINT_MODE);
		else
			this->SetArithmeticMode(PINT_MODE);
		cout << "Arithmetic mode was set to " << ia_mode << ".\n";
	}
	else
	{
		cout << "Arithmetic mode was not set.\n";
	}

	if (vm.count("solver"))
		{
			solver = vm["solver"].as<string>();
			if (solver == "gpe5c")
				this->SetSolverId(GPE_SOLVER5C);
			else
				this->SetSolverId(GPE_SOLVER3C);
			cout << "Solver set to " << solver << ".\n";
		}
		else
		{
			cout << "Solver was not set.\n";
		}


	if (vm.count("exp_mode"))
	{
		ex_mode = vm["exp_mode"].as<string>();
		if (ex_mode == "interval_exp")
			this->SetExperimentMode(INTERVAL_EXP);
		else if (ex_mode == "const_m_exp")
			this->SetExperimentMode(CONST_M_EXP);
		else if (ex_mode == "classical_exp")
					this->SetExperimentMode(CLASSICAL_EXP);
		else
			this->SetExperimentMode(CLASSICAL_EXP);
		cout << "Experiment mode was set to " << ex_mode << ".\n";
	}
	else
	{
		cout << "Experiment mode was not set.\n";
	}

	if (vm.count("m"))
	{
		m = vm["m"].as<int>();
		this->SetM(m);
		this->SetN(m);
		cout << "Grid size m=n was set to " << m << ".\n";
	}
	else
	{
		cout << "Grid size was not set.\n";
	}

	if (vm.count("e"))
	{
		e = vm["e"].as<int>();
		this->SetExampleId(e);
		cout << "Example was set to " << e << ".\n";
	}
	else
	{
		cout << "Example was not set.\n";
	}

	if (vm.count("alpha1"))
	{
		alpha1 = vm["alpha1"].as<T>();
		this->SetAlpha1(alpha1);
		cout << "Alpha1 was set to " << alpha1 << ".\n";
	}
	else
	{
		cout << "Alpha1 was not set, used default value.\n";
	}

	if (vm.count("alpha2"))
	{
		alpha2 = vm["alpha2"].as<T>();
		this->SetAlpha2(alpha2);
		cout << "Alpha2 was set to " << alpha2 << ".\n";
	}
	else
	{
		cout << "Alpha2 was not set, used default value.\n";
	}

	if (vm.count("beta1"))
	{
		beta1 = vm["beta1"].as<T>();
		this->SetBeta1(beta1);
		cout << "Beta1 was set to " << beta1 << ".\n";
	}
	else
	{
		cout << "Beta1 was not set, used default value.\n";
	}

	if (vm.count("beta2"))
	{
		beta2 = vm["beta2"].as<T>();
		this->SetBeta2(beta2);
		cout << "Beta2 was set to " << beta2 << ".\n";
	}
	else
	{
		cout << "Beta2 was not set, used default value.\n";
	}
}

template<typename T>
void ConfigReader<T>::readParametersFromConsole()
{
	bool OK, OK1, dint_mode;
	char z, z1;

	cout
			<< "CONFIG READER"
			<< endl;
	cout << "Please insert necessary parameters." << endl;

	//read arithmetic mode
	do
	{
		cout << "Interval arithmetic mode (d - directed, p - proper): ";
		cin >> z;
	} while ((z != 'd') && (z != 'D') && (z != 'p') && (z != 'P'));
	dint_mode = ((z == 'd') || (z == 'D'));
	if (z == 'd')
		this->SetArithmeticMode(DINT_MODE);
	else
		this->SetArithmeticMode(PINT_MODE);

	//read grid size
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

	//read output file name
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
			OK = OK1;
		} while (!OK);
	}
	cout << endl;

	//read experiment mode
	do
	{
		cout << "Experiment mode [c-classic,m-constM, i - interval]: ";
		cin >> z;
	} while ((z != 'c') && (z != 'm') && (z != 'i'));
	if (z == 'c')
		this->SetExperimentMode(CLASSICAL_EXP);
	if (z == 'm')
		this->SetExperimentMode(CONST_M_EXP);
	if (z == 'i')
		this->SetExperimentMode(INTERVAL_EXP);

	//read example id
	cout << "Example id [1,2,3,4]: ";
	int exId = 0;
	while (!(cin >> exId) || (exId < 0) || (exId > 4))
	{
		cin.clear();
		cout << "Allowed values: 1,2,3,4. Select one.";
	}
	this->SetExampleId(exId);

	//read boundary region
	cout << "Alpha1 = ";
	T alpha1;
	while (!(cin >> alpha1))
	{
		cin.clear();
		cout << "Only floating point numbers.";
	}
	cout << "Alpha2 = ";
	T alpha2;
	while (!(cin >> alpha2))
	{
		cin.clear();
		cout << "Only floating point numbers.";
	}
	cout << "Beta1 = ";
	T beta1;
	while (!(cin >> beta1))
	{
		cin.clear();
		cout << "Only floating point numbers.";
	}
	cout << "Beta2 = ";
	T beta2;
	while (!(cin >> beta2))
	{
		cin.clear();
		cout << "Only floating point numbers.";
	}
	this->SetAlpha1(alpha1);
	this->SetAlpha2(alpha2);
	this->SetBeta1(beta1);
	this->SetBeta2(beta2);
}

template<typename T>
void ConfigReader<T>::SetFileName(string fname)
{
	this->file_name = fname;
}

template<typename T>
string ConfigReader<T>::GetFileName()
{
	return this->file_name;
}

template<typename T>
Parameters<T> ConfigReader<T>::GetParameters()
{
	Parameters<T> currParameters;
	currParameters.eps = 1e-16; //don't change!
	currParameters.m = this->GetM();
	currParameters.n = this->GetN();
	currParameters.alpha = this->GetAlpha1();
	currParameters.beta = this->GetBeta1();
	currParameters.delta = this->GetAlpha2();
	currParameters.gamma = this->GetBeta2();
	currParameters.file_name = this->GetFileName();
	currParameters.ia_mode = this->GetArithmeticMode();
	currParameters.exp_mode = this->GetExperimentMode();
	currParameters.selected_solver = this->GetSolverId();
	currParameters.example_id = this->GetExampleId();
	currParameters.print_csv = this->GetPrintCSV();
	return currParameters;
}

template<typename T>
Solvers ConfigReader<T>::SetSolverId(Solvers sol_id)
{
	this->solver_id = sol_id;
	return this->solver_id;
}

template<typename T>
Solvers ConfigReader<T>::GetSolverId()
{
	return this->solver_id;
}

template<typename T>
void ConfigReader<T>::SetDefaultParameters()
{
	//default configuration
	this->exp_mode =  INTERVAL_EXP; //CLASSICAL_EXP; //CONST_M_EXP; //INTERVAL_EXP; //CLASSICAL_EXP; //CONST_M_EXP
	this->arth_mode = DINT_MODE;
	this->m = 20; //default grid size m=n=20
	this->n = 20;
	this->e = 3; //default example_id = 1
	this->alpha1 = 1;
	this->alpha2 = 2;
	this->beta1 = 1;
	this->beta2 = 2;
	this->file_name = "ex03dint.txt";
	this->print_csv = false;
	this->solver_id = GPE_SOLVER3C;
}

//The explicit instantiation part
template class ConfigReader<long double>;

}
/* namespace intervalarth */
