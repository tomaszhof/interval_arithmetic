//============================================================================
// Name        : IntervalAnalyzer.cpp
// Author      : Tomasz Hoffmann
// Version     :
// Copyright   : TH
// Description : Analyzer for interval results
//============================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <string>
#include <vector>
#include <limits.h>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/exception/diagnostic_information.hpp>

using namespace std;

namespace alg = boost::algorithm;
namespace fs = boost::filesystem;

int main(int ac, char *av[])
{
	long double eps = std::numeric_limits<long double>::epsilon();
	try
	{

		if (ac != 7)
		{
			cout << "Interval Analyzer" << endl;
			cout << "Usage: \n" << fs::basename(av[0]) << " <l_filename>"
					<< " <r_filename>" << " <e_filename>" << " <f_filename>"
					<< " <oe_filename" << " <of_filename>" << endl << endl;
			cout << "where: \n" << "l_ - left ends of interval solution \n"
					<< "r_ - right ends of interval solution \n"
					<< "e_ - exact solution (if known) \n"
					<< "f_ - fp arithmetic solutions \n"
					<< "oe_ - exact distribution output file name \n"
					<< "of_ - fp distribution output file name \n" << endl;
			return 0;
		}
		ifstream l_file(av[1]);
		ifstream r_file(av[2]);
		ifstream e_file(av[3]);
		ifstream f_file(av[4]);
		fstream oe_file(av[5], fstream::out);
		fstream of_file(av[6], fstream::out);

		//string fname = av[4];
		//		o_file.open(fname.c_str(), fstream::out);

		int dprec = std::numeric_limits<long double>::digits10;

		bool allopen = true;
		allopen &= l_file.is_open();
		allopen &= r_file.is_open();
		allopen &= e_file.is_open();
		allopen &= f_file.is_open();

		if (!allopen)
		{
			cout << "One or more files cannot be open! \n";
			return 0;
		}
		else
		{
			string line;
			vector<string> v;
			vector<long double> lv;
			vector<long double> rv;
			vector<long double> ev;
			vector<long double> fv;

			while ((!l_file.eof()) && (!r_file.eof()) && (!e_file.eof())
					&& (!f_file.eof()))
			{
				getline(l_file, line);
				alg::split(v, line, alg::is_any_of(";"));
				vector<string>::iterator it;
				if (v.size() < 2)
					break;
//				cout << "Start left..." << endl;
				for (it = v.begin(); it != v.end(); ++it)
				{
					lv.push_back(boost::lexical_cast<long double>(*it));
				}
//				cout << "Read left done." << endl;

				getline(r_file, line);
//				cout << "Start right..." << endl;
				alg::split(v, line, alg::is_any_of(";"));
				if (v.size() < 2)
					break;
				for (it = v.begin(); it != v.end(); ++it)
				{
					rv.push_back(boost::lexical_cast<long double>(*it));
				}
//				cout << "Read right done." << endl;

				getline(e_file, line);
//				cout << "Start exact..." << endl;
				if (v.size() < 2)
					break;
				alg::split(v, line, alg::is_any_of(";"));
				for (it = v.begin(); it != v.end(); ++it)
				{
//					cout << (*it);
					ev.push_back(boost::lexical_cast<long double>(*it));
				}
//				cout << "Read exact done." << endl;

				getline(f_file, line);
				alg::split(v, line, alg::is_any_of(";"));
				for (it = v.begin(); it != v.end(); ++it)
				{
					fv.push_back(boost::lexical_cast<long double>(*it));
				}
//				cout << "Read floating-point done." << endl;

				if ((lv.size() == rv.size()) && (rv.size() == ev.size())
						&& (ev.size() == fv.size()))
				{
					long double l, r, e, f, w1, w, m;

					vector<long double> res;
					//caluclate exact solutions distribution
					for (size_t i = 0; i < lv.size(); ++i)
					{
						l = lv.at(i);
						r = rv.at(i);
						e = ev.at(i);
						w = r - l;
						m = (r + l) / 2;
						if ((w < eps) || (e == r) || (e == l))
						{
							res.push_back(0.5);
							continue;
						}
						w1 = abs(e - m);
						if ((w1 / w) > 0.5)
						{
							cout << " \n Wrong element: "
									<< std::setprecision(dprec) << e;
							cout << std::setprecision(dprec) << "\n [ " << l
									<< " ; " << r << " ] \n";
							cout << std::setprecision(dprec) << "m = " << m
									<< " \n";
						}
						res.push_back(w1 / w);
					}

					stringstream ss;
					vector<long double>::iterator it;
					for (it = res.begin(); it != res.end(); ++it)
					{
						ss << std::setprecision(dprec) << (*it) << ";";
					}
					ss << endl;
					oe_file << ss.str();

					//caluclate floating-pint solutions distribution
					res.clear();
					for (size_t i = 0; i < lv.size(); ++i)
					{
						l = lv.at(i);
						r = rv.at(i);
						f = fv.at(i);
						w = r - l;
						m = (r + l) / 2;
						if ((w < eps) || (f == r) || (f == l))
						{
							res.push_back(0.5);
							continue;
						}
						w1 = abs(f - m);
						if ((w1 / w) > 0.5)
						{
							cout << " \n Wrong element: "
									<< std::setprecision(dprec) << e;
							cout << std::setprecision(dprec) << "\n [ " << l
									<< " ; " << r << " ] \n";
							cout << std::setprecision(dprec) << "m = " << m
									<< " \n";
						}
						res.push_back(w1 / w);
					}

					ss.str(string());
					for (it = res.begin(); it != res.end(); ++it)
					{
						ss << std::setprecision(dprec) << (*it) << ";";
					}
					ss << endl;
					of_file << ss.str();

					//clear input vectors for current row
					lv.clear();
					rv.clear();
					ev.clear();
					fv.clear();
				}
			}

		}
		l_file.close();
		r_file.close();
		e_file.close();
		oe_file.close();
	} catch (boost::exception const& ex)
	{
		cout
				<< "Execution problem "
						+ boost::current_exception_diagnostic_information();
	}
	return 0;
}
