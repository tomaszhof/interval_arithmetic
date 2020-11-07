/*
 * Configuration.h
 *
 *  Created on: 7 lis 2020
 *      Author: tomhof
 */

#ifndef SRC_CONFIGURATION_H_
#define SRC_CONFIGURATION_H_

#include <string>
#include <iostream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace std;

class Configuration {

public:
	string mode;
	bool use_boost = false;
	int example = 1;
	int n = 10;

	Configuration();
	virtual ~Configuration();

	void parseCommandArgs(int ac, char* av[]);
};

#endif /* SRC_CONFIGURATION_H_ */
