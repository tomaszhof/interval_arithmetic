/*
 * Configuration.cpp
 *
 *  Created on: 7 lis 2020
 *      Author: tomhof
 */

#include "Configuration.h"


Configuration::Configuration() {
	// TODO Auto-generated constructor stub

}

Configuration::~Configuration() {
	// TODO Auto-generated destructor stub
}

void Configuration::parseCommandArgs(int ac, char *av[]) {
	string ia_mode;
	string ex_mode;

	//default options values
	string configFileName, outFileName, programMode;
	int m, e;

	//parsing values
	po::options_description desc("Allowed options");

	desc.add_options()("help", "print help")("n", po::value<int>(),
			"set grid size n")("e", po::value<int>(),
			"set example by id [1-4]")("mode", po::value<string>(), "set program mode")
			("f", po::value<string>(), "set results file name")("arth_mode", po::value<string>(),
			"set interval arithmetic mode [pint-proper, dint-directed]");

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
		if (programMode == "boost")
			this->use_boost = true;
	}

	if (vm.count("f"))
	{
		this->fname = vm["f"].as<string>();
	}

	if (vm.count("n"))
	{
		this->n = vm["n"].as<int>();
	}

	if (vm.count("e"))
	{
		this->example = vm["e"].as<int>();
	}
}
