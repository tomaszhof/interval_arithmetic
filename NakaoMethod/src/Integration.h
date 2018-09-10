/*
 * Integration.h
 *
 *  Created on: 02.09.2018
 *      Author: numeric
 */

#ifndef INTEGRATION_H_
#define INTEGRATION_H_

#include <cmath>
#include <vector>

namespace integration {

std::vector<long double> xv;

/*! Implementation of Gauss-Legendre quadrature
 *  http://en.wikipedia.org/wiki/Gaussian_quadrature
 *  http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
 *
 */
template<int N>
class GaussLegendreQuadrature {
public:
	enum {
		eDEGREE = N
	};

	/*! Compute the integral of a functor
	 *
	 *   @param a    lower limit of integration
	 *   @param b    upper limit of integration
	 *   @param f    the function to integrate
	 *   @param err  callback in case of problems
	 */
	template<typename Function>
	double integrate(double a, double b, Function f) {
		double p = (b - a) / 2;
		double q = (b + a) / 2;
		const LegendrePolynomial& legpoly = s_LegendrePolynomial;

		double sum = 0;
		for (int i = 1; i <= eDEGREE; ++i) {
			sum += legpoly.weight(i) * f(p * legpoly.root(i) + q);
		}

		return p * sum;
	}

	/*! Print out roots and weights for information
	 */
	void print_roots_and_weights(std::ostream& out) const {
		const LegendrePolynomial& legpoly = s_LegendrePolynomial;
		out << "Roots:  ";
		for (int i = 0; i <= eDEGREE; ++i) {
			out << ' ' << legpoly.root(i);
		}
		out << '\n';
		out << "Weights:";
		for (int i = 0; i <= eDEGREE; ++i) {
			out << ' ' << legpoly.weight(i);
		}
		out << '\n';
	}
private:
	/*! Implementation of the Legendre polynomials that form
	 *   the basis of this quadrature
	 */
	class LegendrePolynomial {
	public:
		LegendrePolynomial() {
			// Solve roots and weights
			for (int i = 0; i <= eDEGREE; ++i) {
				double dr = 1;

				// Find zero
				Evaluation eval(cos(M_PI * (i - 0.25) / (eDEGREE + 0.5)));
				do {
					dr = eval.v() / eval.d();
					eval.evaluate(eval.x() - dr);
				} while (fabs(dr) > 2e-16);

				this->_r[i] = eval.x();
				this->_w[i] = 2
						/ ((1 - eval.x() * eval.x()) * eval.d() * eval.d());
			}
		}

		double root(int i) const {
			return this->_r[i];
		}
		double weight(int i) const {
			return this->_w[i];
		}
	private:
		double _r[eDEGREE + 1];
		double _w[eDEGREE + 1];

		/*! Evaluate the value *and* derivative of the
		 *   Legendre polynomial
		 */
		class Evaluation {
		public:
			explicit Evaluation(double x) :
					_x(x), _v(1), _d(0) {
				this->evaluate(x);
			}

			void evaluate(double x) {
				this->_x = x;

				double vsub1 = x;
				double vsub2 = 1;
				double f = 1 / (x * x - 1);

				for (int i = 2; i <= eDEGREE; ++i) {
					this->_v = ((2 * i - 1) * x * vsub1 - (i - 1) * vsub2) / i;
					this->_d = i * f * (x * this->_v - vsub1);

					vsub2 = vsub1;
					vsub1 = this->_v;
				}
			}

			double v() const {
				return this->_v;
			}
			double d() const {
				return this->_d;
			}
			double x() const {
				return this->_x;
			}

		private:
			double _x;
			double _v;
			double _d;
		};
	};

	/*! Pre-compute the weights and abscissae of the Legendre polynomials
	 */
	static LegendrePolynomial s_LegendrePolynomial;
};

template<int N>
typename GaussLegendreQuadrature<N>::LegendrePolynomial GaussLegendreQuadrature<
		N>::s_LegendrePolynomial;

// This to avoid issues with exp being a templated function
//double RosettaExp(double x) {
//    return exp(x);
//}

//usage example
//int main() {
//    Rosetta::GaussLegendreQuadrature<5> gl5;
//
//    std::cout << std::setprecision(10);
//
//    gl5.print_roots_and_weights(std::cout);
//    std::cout << "Integrating Exp(X) over [-3, 3]: " << gl5.integrate(-3., 3., RosettaExp) << '\n';
//    std::cout << "Actual value:                    " << RosettaExp(3) - RosettaExp(-3) << '\n';
//}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

//    {---------------------------------------------------------------------------}
//    {                                                                           }
//    {  The function GaussLegendre calculates an approximate value of the        }
//    {  integral from a function f(x) in the interval [-1,1] by the Gauss-       }
//    {  -Legendre quadrature formula.                                            }
//    {  Data:                                                                    }
//    {    f - a Turbo Pascal function which for the given x calculates the value }
//    {        of integrand at x,                                                 }
//    {    n - number of nodes minus 1 (the nodes are numbered from 0 to n),      }
//    {    x - an array containing roots of the Legendre polynomial of the degree }
//    {        n+1 (the element x[k] should contain the value of the k-th root;   }
//    {        k=0,1,...,n div 2).                                                }
//    {  Result:                                                                  }
//    {    GaussLegendre(f,n,x,st) - approximate value of the integral.           }
//    {  Other parameters:                                                        }
//    {    st - a variable which within the function GaussLegendre is assigned    }
//    {         the value of:                                                     }
//    {           1, if n<1,                                                      }
//    {           0, otherwise.                                                   }
//    {         Note: If st=1, then GaussLegendre(f,n,x,st) is not calculated.    }
//    {  Unlocal identifiers:                                                     }
//    {    vector - a type identifier of extended array [q1..qn2], where q1<=1    }
//    {             and qn2>=n div 2,                                             }
//    {    fx     - a procedural-type identifier defined as follows               }
//    {               type fx = function (x : Extended): Extended;                }
//    {  Note: A function passed as a parameter should be declared with a far     }
//    {        directive or compiled in the $F+ state.                            }
//    {                                                                           }
//    {---------------------------------------------------------------------------}
template<typename F>
long double GaussLegendre(F f, int n, std::vector<long double> x, int st) {
	int k, l, m;
	long double p0, p1, p2, s, x1;
	if (n > 0) {
		st = 0;
		s = 0;
		l = n / 2;
		for (k = 0; k <= n; ++k) {
			p0 = 1;
			if (k <= l)
				x1 = x[k];
			else
				x1 = -x[n - k];
			p1 = x1;
			for (m = 1; m <= n - 1; ++m) {
				p2 = ((2 * m + 1) * x1 * p1 - m * p0) / (m + 1);
				p0 = p1;
				p1 = p2;
			}
			s = s - 2 * f(x1) * (x1 * x1 - 1) / ((n + 1) * (n + 1) * p1 * p1);
		}
		return s;
	} else
		st = 1;

	return -1;
}

//long double e(long double x){
//	return std::exp(x);
//}

template<typename F>
long double GaussLegendreTH(F f, int n, std::vector<long double> X, std::vector<long double> A, long double a, long double b, int st) {
	int k, l, m;
	long double p0, p1, p2, s, x1;
	long double abd2 = (a + b) / 2.0;
	long double bmad2 = (b - a) / 2.0;
	long double tk;

	if (n > 0) {
		st = 0;
		s = 0;
		l = n / 2;
		for (k = 0; k <= n; ++k) {
			tk = abd2 + bmad2*X[k];
			s = s + A[k]* f(tk);
		}
		s = bmad2 * s;
		return s;
	} else
		st = 1;
	return -1;
}

template<typename F, typename PHI>
long double GaussLegendreTH(F f, PHI phi, int i, int n, std::vector<long double> X, std::vector<long double> A, long double a, long double b, int st) {
	int k, l, m;
	long double p0, p1, p2, s, x1;
	long double abd2 = (a + b) / 2.0;
	long double bmad2 = (b - a) / 2.0;
	long double tk;

	if (n > 0) {
		st = 0;
		s = 0;
		l = n / 2;
		for (k = 0; k <= n; ++k) {
			tk = abd2 + bmad2*X[k];
			s = s + A[k]* f(tk)*phi(i, tk);
		}
		s = bmad2 * s;
		return s;
	} else
		st = 1;
	return -1;
}

template<typename PHI>
long double GaussLegendreTH(PHI phi, int i, int j, int n, std::vector<long double> X, std::vector<long double> A, long double a, long double b, int st) {
	int k, l, m;
	long double p0, p1, p2, s, x1;
	long double abd2 = (a + b) / 2.0;
	long double bmad2 = (b - a) / 2.0;
	long double tk;

	if (n > 0) {
		st = 0;
		s = 0;
		l = n / 2;
		for (k = 0; k <= n; ++k) {
			tk = abd2 + bmad2*X[k];
			s = s + A[k]* phi(tk, i)*phi(j, tk);
		}
		s = bmad2 * s;
		return s;
	} else
		st = 1;
	return -1;
}

}

#endif /* INTEGRATION_H_ */
