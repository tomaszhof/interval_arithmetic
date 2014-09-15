/*
 * Tester.cpp
 *
 *  Created on: 12-04-2014
 *      Author: thof
 */

#include "Tester.h"

namespace intervalarth
{

Tester::Tester()
{
	// TODO Auto-generated constructor stub

}

Tester::~Tester()
{
	// TODO Auto-generated destructor stub
}

void Tester::PrintBinary(long double x)
{
	//string sa = "0.1";
//	mpfr_t rop, op;
//	mpfr_init2(rop, 80);
//	mpfr_set_str(rop, sa.c_str(), 10, MPFR_RNDN);

//cout << "Internal representation of the interval:" << endl <<
// "i.a =" << endl;
	long double le = x; //mpfr_get_ld(rop, MPFR_RNDN);
	ieee854_long_double ie;
	ie.d = le;
	std::bitset<32> bs_m0(ie.ieee.mantissa0);
	std::bitset<32> bs_m1(ie.ieee.mantissa1);
	std::bitset<15> bs_exp(ie.ieee.exponent);
	std::bitset<1> bs_neg(ie.ieee.negative);
	std::bitset<16> bs_empty(ie.ieee.empty);
	cout << bs_neg << bs_exp << bs_m0 << bs_m1 << endl;
	cout
			<< "|<---exponent--><---------------------------mantissa--------------------------->"
			<< endl << "sign" << endl;

}

void Tester::IEndsToStrings(interval i, string& left, string& right)
{
	stringstream ss;
	int prec = std::numeric_limits<long double>::digits10;
	ss.setf(std::ios_base::scientific);
	ss << std::setprecision(prec) << i.a;
	left = ss.str();
	ss.str(std::string());

	ss << std::setprecision(prec) << i.b;
	right = ss.str();
	ss.clear();
}

void Tester::ArithmeticTest()
{
		string x1, x2, y1, y2, z1, z2;
		interval x, y, z;

		IntervalArithmetic ia;
		string xa, xb;
		x = ia.IntRead("3.14");
		interval three = {3,3};
		x = ia.IDiv(x, three);
		ia.IEndsToStrings(x,xa,xb);
		cout << xa << " ; " << xb;

		cout << endl;
		cout << "Test for real numbers:" << endl;
		cout << "x = ";
		cin >> x1;
		cout << endl;
		x = ia.IntRead(x1);
		ia.IEndsToStrings(x, z1, z2);
		cout << "x converted" << endl;
		cout << "to machine interval X = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of X:" << endl << "X.a =" << endl;
		PrintBinary(x.a);
		cout << "X.b =" << endl;
		PrintBinary(x.b);

		cout << "y = ";
		cin >> y1;
		cout << endl;
		y = ia.IntRead(y1);
		ia.IEndsToStrings(y, z1, z2);
		cout << "y converted" << endl;
		cout << "to machine interval Y = [" << z1 << "," << z2 << "]" << endl;
		cout << endl;

		cout << "Internal representation of Y:" << endl << "Y.a =" << endl;
		PrintBinary(y.a);
		cout << "Y.b =" << endl;
		PrintBinary(y.b);

		z = ia.IAdd(x, y);
		ia.IEndsToStrings(z, z1, z2);
		cout << x1 << " + " << y1 << " = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);

		z = ia.ISub(x, y);
		ia.IEndsToStrings(z, z1, z2);
		cout << x1 << " - " << y1 << " = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);

		z = ia.IMul(x, y);
		ia.IEndsToStrings(z, z1, z2);
		cout << x1 << " * " << y1 << " = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);

		z = ia.IDiv(x, y);
		ia.IEndsToStrings(z, z1, z2);
		cout << x1 << " - " << y1 << " = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);
		cout << endl;

		cout << "Test for proper intervals:" << endl;
		cout << "X.a = ";
		cin >> x1;
		x.a = ia.LeftRead(x1);
		cout << "X.b (X.b>=X.a) = ";
		cin >> x2;
		x.b = ia.RightRead(x2);
		ia.IEndsToStrings(x, z1, z2);
		cout << "X converted" << endl;
		cout << "to machine interval X = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of X:" << endl << "X.a =" << endl;
		PrintBinary(x.a);
		cout << "X.b =" << endl;
		PrintBinary(x.b);

		cout << "Y.a = ";
		cin >> y1;
		y.a = ia.LeftRead(y1);
		cout << "Y.b (Y.b>=Y.a) = " << endl;
		cin >> y2;
		y.b = ia.RightRead(y2);
		ia.IEndsToStrings(y, z1, z2);
		cout << "Y converted" << endl;
		cout << "to machine interval Y = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Y:" << endl << "Y.a =" << endl;
		PrintBinary(y.a);
		cout << "Y.b =" << endl;
		PrintBinary(y.b);

		cout << endl;
		z = ia.IAdd(x, y);
		ia.IEndsToStrings(z, z1, z2);
		cout << "X + Y = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);

		z = ia.ISub(x, y);
		ia.IEndsToStrings(z, z1, z2);
		cout << "X - Y = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);
		z = ia.IMul(x, y);
		ia.IEndsToStrings(z, z1, z2);
		cout << "X * Y = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);
		z = ia.IDiv(x, y);
		ia.IEndsToStrings(z, z1, z2);
		cout << "X / Y = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);
}

void Tester::ArithmeticTestNew()
{
		string x1, x2, y1, y2, z1, z2;
		Interval<long double> x, y, z;
		string xa, xb;
		x = x.IntRead("3.14");
		Interval<long double> three (3,3);
		x = x / three;

		x.IEndsToStrings(xa,xb);
		cout << xa << " ; " << xb;

		cout << endl;
		cout << "Test for real numbers:" << endl;
		cout << "x = ";
		cin >> x1;
		cout << endl;
		x = x.IntRead(x1);
		x.IEndsToStrings(z1, z2);
		cout << "x converted" << endl;
		cout << "to machine interval X = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of X:" << endl << "X.a =" << endl;
		PrintBinary(x.a);
		cout << "X.b =" << endl;
		PrintBinary(x.b);

		cout << "y = ";
		cin >> y1;
		cout << endl;
		y = y.IntRead(y1);
		y.IEndsToStrings(z1, z2);
		cout << "y converted" << endl;
		cout << "to machine interval Y = [" << z1 << "," << z2 << "]" << endl;
		cout << endl;

		cout << "Internal representation of Y:" << endl << "Y.a =" << endl;
		PrintBinary(y.a);
		cout << "Y.b =" << endl;
		PrintBinary(y.b);

		z = x + y;
		z.IEndsToStrings(z1, z2);
		cout << x1 << " + " << y1 << " = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);

		z = x - y;
		z.IEndsToStrings(z1, z2);
		cout << x1 << " - " << y1 << " = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);

		z = x * y;
		z.IEndsToStrings(z1, z2);
		cout << x1 << " * " << y1 << " = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);

		z = x / y;
		z.IEndsToStrings(z1, z2);
		cout << x1 << " - " << y1 << " = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);
		cout << endl;

		cout << "Test for proper intervals:" << endl;
		cout << "X.a = ";
		cin >> x1;
		x.a = x.LeftRead(x1);
		cout << "X.b (X.b>=X.a) = ";
		cin >> x2;
		x.b = x.RightRead(x2);
		x.IEndsToStrings(z1, z2);
		cout << "X converted" << endl;
		cout << "to machine interval X = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of X:" << endl << "X.a =" << endl;
		PrintBinary(x.a);
		cout << "X.b =" << endl;
		PrintBinary(x.b);

		cout << "Y.a = ";
		cin >> y1;
		y.a = y.LeftRead(y1);
		cout << "Y.b (Y.b>=Y.a) = " << endl;
		cin >> y2;
		y.b = y.RightRead(y2);
		y.IEndsToStrings(z1, z2);
		cout << "Y converted" << endl;
		cout << "to machine interval Y = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Y:" << endl << "Y.a =" << endl;
		PrintBinary(y.a);
		cout << "Y.b =" << endl;
		PrintBinary(y.b);

		cout << endl;
		z = x + y;
		z.IEndsToStrings(z1, z2);
		cout << "X + Y = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);

		z = x - y;
		z.IEndsToStrings(z1, z2);
		cout << "X - Y = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);
		z = x * y;
		z.IEndsToStrings(z1, z2);
		cout << "X * Y = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);
		z = x / y;
		z.IEndsToStrings(z1, z2);
		cout << "X / Y = [" << z1 << "," << z2 << "]" << endl;
		cout << "Internal representation of Z:" << endl << "Z.a =" << endl;
		PrintBinary(z.a);
		cout << "Z.b =" << endl;
		PrintBinary(z.b);
}
void Tester::StuffTest()
{
	//
	//		string sa = "0.1";
	//		mpfr_t rop, op;
	//		mpfr_init2(rop, 80);
	//		mpfr_set_str(rop, sa.c_str(), 10, MPFR_RNDN);
	//		long double le = mpfr_get_ld(rop, MPFR_RNDN);
	//		ieee854_long_double ie;
	//		ie.d = le;
	//		std::bitset<32> bs_m0(ie.ieee.mantissa0);
	//		std::bitset<32> bs_m1(ie.ieee.mantissa1);
	//		std::bitset<15> bs_exp(ie.ieee.exponent);
	//		std::bitset<1> bs_neg(ie.ieee.negative);
	//		std::bitset<16> bs_empty(ie.ieee.empty);
	//		cout << "mantisa = " << bs_m0 << bs_m1  << endl;
	//		cout << "exponent = " << bs_exp  << endl;
	//		cout << "negative = " << bs_neg  << endl;
	//		cout << "empty = " << bs_empty  << endl;
	//		|<---exponent--><---------------------------mantissa--------------------------->
	//		sign
	/*
	 interval i1 = {1,1};
	 interval i2 = {2,2};
	 THashMap* tm = new THashMap();
	 tm->ToMap(3, i1);
	 interval fm = tm->FromMap(3);
	 tm->ToMap(3, i2);
	 fm = tm->FromMap(3);
	 int n,m,st;
	 n = m = 20;
	 //test
	 interval intalpha = {1,1};
	 interval intbeta ={1,1};
	 const long double Mconst = 1627;
	 const long double eps = 1.0e-16;
	 intervalvector X;
	 DiffPoisson* dp = new DiffPoisson();
	 dp->DIntPoisson(n,m, intalpha, intbeta, eps, Mconst, X, st);

	 IntervalArithmetic ia;
	 interval sqr2 = ia.ISqr3();
	 cout << "sqr2 left endpoint  = " << setprecision(30) << sqr2.a << endl;
	 cout << "sqr2 right endpoint  = " << setprecision(30) << sqr2.b << endl;
	 cout << "sqr2 interval width  = " << setprecision(30) << sqr2.b - sqr2.a << endl;
	 long double sqr2d = boost::lexical_cast<long double>("1.414213562373095049");
	 cout << "sqr2d  = " << setprecision(30) << sqr2d << endl;
	 mpfr_t x,y;
	 mpfr_init2(x,80);

	 string ds = "1.414213562373095049";
	 mpfr_set_str(x,ds.c_str(), 10, MPFR_RNDD);
	 mpfr_printf ("%.128R*f", MPFR_RNDN, x);
	 */
	/* list_elem* h = NULL;
	 list_elem* e = NULL;

	 for (int i = 0; i < 10; i++)
	 {
	 e = new list_elem();
	 e->value = i;
	 add_elem_to_list(e,h);
	 }
	 print_list(h); */

	/*
	 IntervalArithmetic i; //= new IntervalArithmetic();
	 string sv = "2.449489742783178098";
	 long double ldt = boost::lexical_cast<long double>(sv);
	 long double ld = strtold(sv.c_str(), NULL);
	 ld = 2.449489742783178098L;
	 int e = 0;
	 long double ld_m = frexp(ld, &e);
	 double dd = strtod(sv.c_str(), NULL);
	 double dd_m = frexp(dd, &e);
	 char c1 = '0';
	 int c1i = c1 - '0';
	 cout << "c1i=" << c1i << endl;
	 //for (int i = 15; i > 0; i--)
	 //	cout << i << endl;
	 //int m,n;
	 m = n = 172;
	 int mn = (n*m-n-m+4)*(n*m-n-m+4) / 4;
	 long double * tabd = new long double[mn];
	 for (int i=0; i < mn; i++)
	 {
	 tabd[i] = 1.5*i;
	 }
	 tabd[28900] = 11.2;
	 cout << "tabd[28900]=" << tabd[28900] << endl;

	 cout << "double dd = " << setprecision(25) << dd << endl;
	 cout << "long double ld = " << setprecision(25) << ld << endl;
	 cout << "long double ldt = " << setprecision(25) << ldt << endl;

	 cout << "size(long double)=" <<  sizeof(long double) << endl;

	 cout << "size(double)=" <<  sizeof(double) << endl;
	 cout << "size(float)=" <<  sizeof(float) << endl;
	 cout << "numeric_limits(long double)=" <<  std::numeric_limits<long double>::digits10 << endl;
	 cout << "numeric_limits(double)=" <<  std::numeric_limits<double>::digits10 << endl;
	 */
}
} /* namespace intervalarth */
