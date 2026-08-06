// Evaluation of the integrals arising in the bubble diagrams
// Amina Alic, May 2017
// later updates by RZ

// Calculates integrals over frequency in the bubble formula for transport properties
// using analytical expressions for integrals over energy that were evaluated in Mathematica.
// 
// The full description is given in the specification (bubble.pdf)
// Syntax: bubble <m> <n> <o> <T> <mu> <resigma> <imsigma>
// The input files must contain a table of space-separated (frequency, self-energy) pairs
//
// In integration over frequencies, the Gauss-Kronrod rules with 15 points are used. 
// Interpolation is implemented using the GSL, with possibility to choose between linear, 
// cubic and Akima splines.
// In case that there is no analytical expression for the kernel function Phi(epsilon), 
// the user has to provide the file "Phi.dat", which contains a table of space-separated 
// (epsilon, Phi(epsilon)) pairs.

// CHANGE LOG
// 7.6.2017 - code cleanup (rz)
// 10.10.2017 - support for the -f switch
// 22.11.2017 - support for the -e switch
// 23.12.2017 - compute the spectral function (debugging aid) (rz)
//            - clipping of Im Sigma
// 6.8.2026 - merge support for the -d, -e, and -f switches

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <cassert>
#include <string>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <complex>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <getopt.h>
#include <unistd.h>

#include "Faddeeva.hh"

using namespace std;

// The Mathematica-generated expressions use integer coefficients with
// complex<double>. Normalize those mixed operations for standard libraries
// that require both operands of the complex scalar overloads to have the same
// value type.
complex<double> operator+(const complex<double>& lhs, int rhs) { return lhs + static_cast<double>(rhs); }
complex<double> operator+(int lhs, const complex<double>& rhs) { return static_cast<double>(lhs) + rhs; }
complex<double> operator-(const complex<double>& lhs, int rhs) { return lhs - static_cast<double>(rhs); }
complex<double> operator-(int lhs, const complex<double>& rhs) { return static_cast<double>(lhs) - rhs; }
complex<double> operator*(const complex<double>& lhs, int rhs) { return lhs * static_cast<double>(rhs); }
complex<double> operator*(int lhs, const complex<double>& rhs) { return static_cast<double>(lhs) * rhs; }
complex<double> operator/(const complex<double>& lhs, int rhs) { return lhs / static_cast<double>(rhs); }
complex<double> operator/(int lhs, const complex<double>& rhs) { return static_cast<double>(lhs) / rhs; }

// Interpolation object for re and im parts of Sigma(omega)
gsl_interp_accel *acc_reSigma;
gsl_spline *spline_reSigma;
gsl_interp_accel *acc_imSigma;
gsl_spline *spline_imSigma;
double omega_min, omega_max;
double reSigma_asymp_neg, reSigma_asymp_pos; // asymptotic values

// Interpolation object for generic Phi(epsilon) function
gsl_interp_accel *acc_Phi;
gsl_spline *spline_Phi;
double eps_min, eps_max; // Interval boundaries

const string VERSION = "1.4";

// Mandatory parameters
int m, n, o;
double T;
double mu;
string fnReSigma, fnImSigma;

// Optional parameters (with defaults)
int b = 1;
bool verbose = false;
int key = GSL_INTEG_GAUSS15;
double abs_error = 1.0e-7;
double rel_error = 1.0e-8;
double cutoff = 15.0;
string fnPhi = "Phi.dat";
bool calcdos = false; // compute the spectral function
enum ff { f_f, f_derivative };
ff f_type = f_derivative;
int e = 0; // power of epsilon

const double EPSILON = 1e-10; // some very small value...
// For clipping Im Sigma. 
const double EPSILON2 = 1e-8; // some small value...

void about()
{
   cout << "bubble version " << VERSION << endl;
}

void usage()
{
   about();
   cout << "Usage: bubble <m> <n> <o> <T> <mu> <resigma> <imsigma>" << endl;
   cout << "Options:" << endl;
   cout << "-v : increase verbosity" << endl;
   cout << "-i : interpolation (default = 1, 1=>linear, 2=>cspline, 3=>Akima spline)" << endl;
   cout << "-k : integration rule (default = 1, 1 => 15, 2 => 21, etc.)" << endl;
   cout << "-a : absolute error (default = 1e-7)" << endl;
   cout << "-r : relative error (default = 1e-8)" << endl;
   cout << "-c : frequency interval cutoff in units of T (default = 15)" << endl;
   cout << "-p : filename for Phi tables (default = Phi.dat)" << endl;
   cout << "-d : compute the epsilon integrals only, for m=0 (output = dos.dat)" << endl;
   cout << "-e : additional power of epsilon when using the m=0 code (default=0)" << endl;
   cout << "-f : switch (-df/dw) to f in the w integration (incompatible with -d)" << endl;
}

void cmd_line(int argc, char *argv[])
{
   char c;
    
   while ((c = getopt(argc, argv, "vi:k:a:r:c:p:dfe:")) != -1) {
      switch (c) {
      case 'v':
	 verbose = true;
	 break;
      case 'i':
	 b = atoi(optarg);
	 break;
      case 'k':
         key = atoi(optarg);
         break;
      case 'a':
	 abs_error = atof(optarg);
	 break;
      case 'r':
	 rel_error = atof(optarg);
	 break;
      case 'c':
	 cutoff = atof(optarg);
	 break;
      case 'p':
	 fnPhi = string(optarg);
	 break;
      case 'd':
	 calcdos = true;
	 break;
      case 'f':
	 f_type = f_f;
	 break;
      case 'e':
	 e = atoi(optarg);
	 break;
      default:
	 cerr << "Option not implemented." << endl;
	 usage();
	 abort();
      }
   }
    
   int remaining = argc-optind;
    
   if (remaining != 7) {
      usage();
      abort();
   }
   
   m = atoi(argv[optind]);
   n = atoi(argv[optind+1]);
   o = atoi(argv[optind+2]);
   T = atof(argv[optind+3]);
   mu = atof(argv[optind+4]);
   fnReSigma = string(argv[optind+5]);
   fnImSigma = string(argv[optind+6]);

   if (calcdos && m != 0) {
      cerr << "Option -d requires m=0." << endl;
      exit(EXIT_FAILURE);
   }
   if (calcdos && f_type == f_f) {
      cerr << "Options -d and -f cannot be used together." << endl;
      exit(EXIT_FAILURE);
   }
    
   if (verbose) {
      about();
      cout << "m=" << m << "(";
      switch (m) {
      case 0:
	 cout << "generic";
	 break;
      case 1:
	 cout << "flat";
	 break;
      case 2:
	 cout << "eps";
	 break;
      case 3:
	 cout << "semicirc";
	 break;
      case 4:
	 cout << "eps*semicirc";
	 break;
      case 5:
	 cout << "sqrt^(3/2)";
	 break;
      case 6:
	 cout << "eps*sqrt^(3/2)";
	 break;
      case 7:
	 cout << "gaussian";
	 break;
      case 8:
	 cout << "eps*gaussian";
	 break;
      }
      cout << ")";
      cout << " n=" << n << " o=" << o;
      cout << " T=" << T << " mu=" << mu;
      cout << " Sigma: re=" << fnReSigma << ",im=" << fnImSigma << endl;
      cout << "i=" << b << " key=" << key << endl;
      cout << "abs_error=" << abs_error << " rel_error=" << rel_error << " cutoff=" << cutoff << endl;
      cout << "Phi=" << fnPhi << " e=" << e << endl;
      if (calcdos) {
	 cout << "epsilon integrals only (dos.dat)" << endl;
      } else {
	 switch (f_type) {
	 case f_f:
	    cout << "f(T)" << endl;
	    break;
	 case f_derivative:
	    cout << "-df/dw" << endl;
	    break;
	 default:
	    cerr << "Not implemented." << endl;
	    break;
	 }
      }
   }
}

// Faddeeva function
complex<double> Erfi(complex<double> z) {
   double relerr = 0;
   return Faddeeva::erfi(z,relerr);
}

const complex<double> I(0,1);

// Functions Jmn, for m = 1,...,8 and n = 0,1,2,3
double J_10(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return 2.0;
}

double J_11(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return (atan((1.0 - x)/y) + atan((1.0 + x)/y))/M_PI;
}


double J_12(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return ((2*y*(1 + pow(y,2) - pow(x,2)) + (pow(y,2) + pow(1 - x,2)) * (pow(y,2) + pow(1 + x,2))*(atan((1 - x)/y) + atan((1 + x)/y))))/(2.*y*pow(M_PI,2)*(pow(y,2) + pow(1 - x,2))*(pow(y,2) + pow(1 + x,2)));
}

double J_13(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return ((2 * y * (pow(1 + pow(y,2),2) * (3 + 5 * pow(y,2)) - (9 + 2 * pow(y,2) + pow(y,4)) * pow(x,2) + 9 * (1 - y) * (1 + y) * pow(x,4) - 3 * pow(x,6)) + 3 * pow(1 + 2 * (y - x) * (y + x) + pow(pow(y,2) + pow(x,2),2),2) * (atan((1 - x)/y) + atan((1 + x)/y))))/(8. * pow(y,2) * pow(M_PI,3) * pow(1+ 2 * (y - x) * (y + x) + pow(pow(y,2) + pow(x,2),2),2));
}

double J_20(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return 0.0;
}

double J_21(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return (2 * x * atan((1 - x)/y) + 2 * x * atan((1 + x)/y) + y * log(1 - (4*x)/(pow(y,2) + pow(1 + x,2))))/(2.*M_PI);
}


double J_22(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return  (x*(-2*y*(-1 + pow(x,2) + pow(y,2)) - (pow(-1 + x,2) + pow(y,2))*(pow(1 + x,2) + pow(y,2))*atan((-1 + x)/y) +(pow(-1 + x,2) + pow(y,2))*(pow(1 + x,2) + pow(y,2))*atan((1 + x)/y)))/(2.*pow(M_PI,2)*y*(pow(-1 + x,2) + pow(y,2))*(pow(1 + x,2) + pow(y,2)));
}

double J_23(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return -(x*(6*pow(-1 + pow(x,2),3)*y + 2*(-11 + 2*pow(x,2) + 9*pow(x,4))*pow(y,3) + 2*(-5 + 9*pow(x,2))*pow(y,5) + 6*pow(y,7) + 3*pow(pow(x,4) + 2*pow(x,2)*(-1 + pow(y,2)) + pow(1 + pow(y,2),2),2)*atan((-1 + x)/y) - 3*pow(pow(x,4) + 2*pow(x,2)*(-1 + pow(y,2)) + pow(1 + pow(y,2),2),2)*atan((1 + x)/y)))/(8.*pow(M_PI,3)*pow(y,2)*pow(pow(x,4) + 2*pow(x,2)*(-1 + pow(y,2)) + pow(1 + pow(y,2),2),2));
}

double J_30(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return M_PI/2;
}

double J_31(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = (1./2.)*(I*x - I*xbar + (sqrt(1 - x*x) + sqrt(1 - xbar*xbar))*s);
    
   return z.real();
}

double J_32(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = -(-2 - (I * x * s)/sqrt(1 - x*x) + (I * xbar * s)/sqrt(1 - xbar*xbar)- (2 * (- x + xbar + I * (sqrt(1 - x*x) + sqrt(1 - xbar*xbar)) * s)/(x-xbar)))/(4.*M_PI);
   return z.real();
}

double J_33(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = ((pow(1 - x*x,-1.5) + pow(1 - xbar*xbar,-1.5) - (6*(1 - x*xbar + sqrt(1 - x*x)*sqrt(1 - xbar*xbar)))/(sqrt(1 - x*x)*(x - xbar)*(x - xbar)) - (6*(1 - x*xbar + sqrt(1 - x*x)*sqrt(1 - xbar*xbar)))/((x - xbar)*(x - xbar)*sqrt(1 - xbar*xbar)))*s)/(16.*pow(M_PI,2));
    
   return z.real();
}

double J_40(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return 0.0;
}

double J_41(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = (-I*xbar*xbar + xbar*sqrt(1 - xbar*xbar)*s + I*x*x + x*sqrt(1 - x*x)*s)/2.;
   
   return z.real();
}

double J_42(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = -((I*(2*x*sqrt(1 - x*x)*xbar*xbar - x*(sqrt(1-x*x) + sqrt(1 - xbar*xbar)) - xbar*(sqrt(1-x*x) + sqrt(1-xbar*xbar) - 2*x*x*sqrt(1-xbar*xbar)))*s)/(4*M_PI*sqrt(1 - x*x)*(x-xbar)*sqrt(1 - xbar*xbar)));
    
   return z.real();
}

double J_43(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z =  -I*(2 - I*x*(-3 + 2*x*x)*s/pow(1 - x*x,1.5) - (2 - 2*xbar*xbar - 3*I*xbar*s/sqrt(1 - xbar*xbar) + 2*I*xbar*xbar*xbar*s/sqrt(1 - xbar*xbar))/(1 - xbar*xbar) + 6*I*(-I - (x - 2*x*xbar*xbar + pow(xbar,3) + x*sqrt(1 - x*x)*sqrt(1 - xbar*xbar))*s/((x - xbar)*(x - xbar)*sqrt(1 - xbar*xbar))) - 6*I*(-I + (x*x*x + xbar*(1 - 2*x*x + sqrt(1 - x*x)*sqrt(1- xbar*xbar)))*s/(sqrt(1 - x*x)*(x - xbar)*(x - xbar))))/(16.*M_PI*M_PI);
   
   return z.real();
}

double J_50(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return 3*M_PI/8;
}

double J_51(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = (1./4.)*I*(3*x - 2*pow(x,3) - 3*xbar + 2*pow(xbar,3) - 2*I*pow(1 - x*x,1.5)*s - 2*I*pow(1 - xbar*xbar,1.5)*s);
   
   return z.real();
}

double J_52(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = (-pow(x,3) + pow(xbar,3) + I*(x*x*sqrt(1 - x*x) + 2*(sqrt(1 - x*x) + sqrt(1 - xbar*xbar)))*s + xbar*xbar*(-3*x + I*sqrt(1 - xbar*xbar)*s) + 3*x*xbar*(x - I*(sqrt(1 - x*x) + sqrt(1 - xbar*xbar))*s))/(4.*M_PI*(x - xbar));
   
   return z.real();
}

double J_53(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z =  -(3*(2*x*sqrt(1 - x*x)*pow(xbar,3) + 4*(sqrt(1 - x*x) + sqrt(1 - xbar*xbar)) - x*x*(sqrt(1 - x*x)+ 3*sqrt(1 - xbar*xbar)) + xbar*(2*pow(x,3)*sqrt(1 - xbar*xbar) - 4*x*(sqrt(1 - x*x) + sqrt(1 - xbar*xbar))) + xbar*xbar*(-3*sqrt(1 - x*x) - sqrt(1 - xbar*xbar) + 2*x*x*(sqrt(1 - x*x) + sqrt(1 - xbar*xbar))))*s) /(16.*sqrt(1 - x*x)*M_PI*M_PI*(x - xbar)*(x - xbar)*sqrt(1 - xbar*xbar));
   
   return z.real();
}

double J_60(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return 0.0;
}

double J_61(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z =  (-3*I*xbar*xbar + 2*I*pow(xbar,4) + 2*xbar*sqrt(1 - xbar*xbar)*s - 2*pow(xbar,3)*sqrt(1 - xbar*xbar)*s + x*(I*x*(3 - 2*x*x) + 2*pow(1 - x*x,1.5)*s))/4.;
   
   return z.real();
}

double J_62(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z =  -(x*(-3 + 4*x*x) - 3*xbar + 4*pow(xbar,3) - I*sqrt(1 - x*x)*(-1 + 4*x*x)*s - I*sqrt(1 - xbar*xbar)*s + 4*I*xbar*xbar*sqrt(1 - xbar*xbar)*s - (I*(3*I*x*x - 2*I*pow(x,4) - 3*I*xbar*xbar + 2*I*pow(xbar,4) + 2*x*pow(1 - x*x,1.5)*s + 2*xbar*pow(1 - xbar*xbar,1.5)*s))/(x - xbar))/(4.*M_PI);
   
   return z.real();
}

double J_63(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z =  -(3*(sqrt(1 - x*x)*(-1 + 4*x*x)*pow(xbar,3) + xbar*xbar*(-4*x*sqrt(1 - x*x) - 3*x*sqrt(1 - xbar*xbar) + 4*pow(x,3)*sqrt(1 - xbar*xbar)) + x*(-x*x*sqrt(1 - xbar*xbar) + 2*(sqrt(1 - x*x) + sqrt(1 - xbar*xbar))) + xbar*(2*(sqrt(1 - x*x) + sqrt(1 - xbar*xbar)) - x*x*(3*sqrt(1 - x*x) + 4*sqrt(1 - xbar*xbar))))*s)/(16.*sqrt(1 - x*x)*M_PI*M_PI*(x - xbar)*(x - xbar)*sqrt(1 - xbar*xbar));
   
   return z.real();
}

double J_70(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return sqrt(M_PI/2);
}

double J_71(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z =  I*(exp(2*xbar*xbar)*(M_PI*Erfi(sqrt(2)*x) - log(-(1/x)) - log(x)) + exp(2*x*x)*(-(M_PI*Erfi(sqrt(2)*xbar)) + log(-(1/xbar)) + log(xbar)))/(2.*exp(2*(x*x + xbar*xbar))*M_PI);
   
   return z.real();
}

double J_72(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = (2*exp(2*(x*x + xbar*xbar))*x*sqrt(2*M_PI) - exp(2*xbar*xbar)*(1 + 2*x*x)*M_PI*Erfi(sqrt(2)*x) + exp(2*x*x)*M_PI*Erfi(sqrt(2)*xbar) + exp(2*xbar*xbar)*log(-(1/x)) + 2*exp(2*xbar*xbar)*x*x*log(-(1/x)) - exp(2*xbar*xbar)*log(1/x) - 2*exp(2*xbar*xbar)*x*x*log(1/x) - exp(2*x*x)*log(-(1/xbar)) + exp(2*x*x)*log(1/xbar) + 2*exp(2*x*x)*xbar*xbar*(M_PI*Erfi(sqrt(2)*xbar) - log(-(1/xbar)) + log(1/xbar)) - 2*xbar*(exp(2*(x*x + xbar*xbar))*sqrt(2*M_PI) - exp(2*xbar*xbar)*x*M_PI*Erfi(sqrt(2)*x) + exp(2*x*x)*x*M_PI*Erfi(sqrt(2)*xbar) + exp(2*xbar*xbar)*x*log(-(1/x)) - exp(2*xbar*xbar)*x*log(1/x) - exp(2*x*x)*x*log(-(1/xbar)) + exp(2*x*x)*x*log(1/xbar)))/(2.*exp(2*(x*x + xbar*xbar))*M_PI*M_PI*(x - xbar));
   
   return z.real();
}

double J_73(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = (I*(6*exp(2*(x*x + xbar*xbar))*x*sqrt(2*M_PI) + 2*exp(2*(x*x + xbar*xbar))*pow(x,3)*sqrt(2*M_PI) - exp(2*xbar*xbar)*(3 + 5*x*x + 4*pow(x,4))*M_PI*Erfi(sqrt(2)*x) + 3*exp(2*x*x)*M_PI*Erfi(sqrt(2)*xbar) - exp(2*x*x)*x*x*M_PI*Erfi(sqrt(2)*xbar) + 3*exp(2*xbar*xbar)*log(-(1/x)) + 5*exp(2*xbar*xbar)*x*x*log(-(1/x)) + 4*exp(2*xbar*xbar)*pow(x,4)*log(-(1/x)) + 3*exp(2*xbar*xbar)*log(x) + 5*exp(2*xbar*xbar)*x*x*log(x) + 4*exp(2*xbar*xbar)*pow(x,4)*log(x) - 3*exp(2*x*x)*log(-(1/xbar)) + exp(2*x*x)*x*x*log(-(1/xbar)) + 4*exp(2*x*x)*pow(xbar,4)*(M_PI*Erfi(sqrt(2)*xbar) - log(-(1/xbar)) - log(xbar)) - 3*exp(2*x*x)*log(xbar) + exp(2*x*x)*x*x*log(xbar) - 2*exp(2*x*x)*pow(xbar,3)*(exp(2*xbar*xbar)*sqrt(2*M_PI) + 4*x*M_PI*Erfi(sqrt(2)*xbar) - 4*x*log(-(1/xbar)) - 4*x*log(xbar)) - 2*xbar*(3*exp(2*(x*x + xbar*xbar))*sqrt(2*M_PI) + 3*exp(2*(x*x + xbar*xbar))*x*x*sqrt(2*M_PI) - 2*exp(2*xbar*xbar)*x*(1 + 2*x*x)*M_PI*Erfi(sqrt(2)*x) + 2*exp(2*x*x)*x*M_PI*Erfi(sqrt(2)*xbar) + 2*exp(2*xbar*xbar)*x*log(-(1/x)) + 4*exp(2*xbar*xbar)*pow(x,3)*log(-(1/x)) + 2*exp(2*xbar*xbar)*x*log(x) + 4*exp(2*xbar*xbar)*pow(x,3)*log(x) - 2*exp(2*x*x)*x*log(-(1/xbar)) - 2*exp(2*x*x)*x*log(xbar)) + xbar*xbar*(6*exp(2*(x*x + xbar*xbar))*x*sqrt(2*M_PI) - exp(2*xbar*xbar)*(-1 + 4*x*x)*M_PI*Erfi(sqrt(2)*x) + exp(2*x*x)*(5 + 4*x*x)*M_PI*Erfi(sqrt(2)*xbar) - exp(2*xbar*xbar)*log(-(1/x)) + 4*exp(2*xbar*xbar)*x*x*log(-(1/x)) - exp(2*xbar*xbar)*log(x) + 4*exp(2*xbar*xbar)*x*x*log(x) - 5*exp(2*x*x)*log(-(1/xbar)) - 4*exp(2*x*x)*x*x*log(-(1/xbar)) - 5*exp(2*x*x)*log(xbar) - 4*exp(2*x*x)*x*x*log(xbar))))/(4.*exp(2*(x*x + xbar*xbar))*pow(M_PI,3)*pow(x - xbar,2));
   return z.real();
}

double J_80(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return 0.0;
}

double J_81(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z =  (I*(exp(2*xbar*xbar)*x*(M_PI*Erfi(sqrt(2)*x) - log(-(1/x)) - log(x)) - exp(2*x*x)*xbar*(M_PI*Erfi(sqrt(2)*xbar) - log(-(1/xbar)) - log(xbar))))/(2.*exp(2*(x*x + xbar*xbar))*M_PI);
   
   return z.real();
}

double J_82(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = -(-4*exp(2*x*x)*pow(xbar,3)*(M_PI*Erfi(sqrt(2)*xbar) - log(-(1/xbar)) - log(xbar)) + x*(-2*exp(2*(x*x + xbar*xbar))*x*sqrt(2*M_PI) + exp(2*xbar*xbar)*(1 + 4*x*x)*M_PI*Erfi(sqrt(2)*x) - exp(2*x*x)*M_PI*Erfi(sqrt(2)*xbar) - exp(2*xbar*xbar)*log(-(1/x)) - 4*exp(2*xbar*xbar)*x*x*log(-(1/x)) - exp(2*xbar*xbar)*log(x) - 4*exp(2*xbar*xbar)*x*x*log(x) + exp(2*x*x)*log(-(1/xbar)) + exp(2*x*x)*log(xbar)) + xbar*(-(exp(2*xbar*xbar)*(-1 + 4*x*x)*M_PI*Erfi(sqrt(2)*x)) - exp(2*x*x)*M_PI*Erfi(sqrt(2)*xbar) - exp(2*xbar*xbar)*log(-(1/x)) + 4*exp(2*xbar*xbar)*x*x*log(-(1/x)) - exp(2*xbar*xbar)*log(x) + 4*exp(2*xbar*xbar)*x*x*log(x) + exp(2*x*x)*log(-(1/xbar)) + exp(2*x*x)*log(xbar)) + 2*exp(2*x*x)*xbar*xbar*(exp(2*xbar*xbar)*sqrt(2*M_PI) + 2*x*M_PI*Erfi(sqrt(2)*xbar) - 2*x*log(-(1/xbar)) - 2*x*log(xbar)))/(4.*exp(2*(x*x + xbar*xbar))*M_PI*M_PI*(x - xbar));
   
   return z.real();
}

double J_83(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = (-I*(2*exp(2*xbar*xbar)*xbar*xbar*(3*exp(2*x*x)*sqrt(2*M_PI) + x*(-3 + 4*x*x)*M_PI*Erfi(sqrt(2)*x) + (3*x - 4*pow(x,3))*log(-(1/x)) + 3*x*log(x) - 4*pow(x,3)*log(x)) - 8*exp(2*x*x)*pow(xbar,5)*(M_PI*Erfi(sqrt(2)*xbar) - log(-(1/xbar)) - log(xbar)) - x*(6*exp(2*(x*x + xbar*xbar))*x*sqrt(2*M_PI) + 4*exp(2*(x*x + xbar*xbar))*pow(x,3)*sqrt(2*M_PI) - exp(2*xbar*xbar)*(3 + 6*x*x + 8*pow(x,4))*M_PI*Erfi(sqrt(2)*x) + 3*exp(2*x*x)*M_PI*Erfi(sqrt(2)*xbar) + 3*exp(2*xbar*xbar)*log(-(1/x)) + 6*exp(2*xbar*xbar)*x*x*log(-(1/x)) + 8*exp(2*xbar*xbar)*pow(x,4)*log(-(1/x)) + 3*exp(2*xbar*xbar)*log(x) + 6*exp(2*xbar*xbar)*x*x*log(x) + 8*exp(2*xbar*xbar)*pow(x,4)*log(x) - 3*exp(2*x*x)*log(-(1/xbar)) - 3*exp(2*x*x)*log(xbar)) + 4*exp(2*x*x)*pow(xbar,4)*(exp(2*xbar*xbar)*sqrt(2*M_PI) + 4*x*M_PI*Erfi(sqrt(2)*xbar) - 4*x*log(-(1/xbar)) - 4*x*log(xbar)) - 2*exp(2*x*x)*pow(xbar,3)*(4*exp(2*xbar*xbar)*x*sqrt(2*M_PI) + (3 + 4*x*x)*M_PI*Erfi(sqrt(2)*xbar) - (3 + 4*x*x)*log(-(1/xbar)) - 3*log(xbar) - 4*x*x*log(xbar)) + xbar*(8*exp(2*(x*x + xbar*xbar))*pow(x,3)*sqrt(2*M_PI) - exp(2*xbar*xbar)*(-3 + 16*pow(x,4))*M_PI*Erfi(sqrt(2)*x) + 3*exp(2*x*x)*(-1 + 2*x*x)*M_PI*Erfi(sqrt(2)*xbar) - 3*exp(2*xbar*xbar)*log(-(1/x)) + 16*exp(2*xbar*xbar)*pow(x,4)*log(-(1/x)) - 3*exp(2*xbar*xbar)*log(x) + 16*exp(2*xbar*xbar)*pow(x,4)*log(x) + 3*exp(2*x*x)*log(-(1/xbar)) - 6*exp(2*x*x)*x*x*log(-(1/xbar)) + 3*exp(2*x*x)*log(xbar) - 6*exp(2*x*x)*x*x*log(xbar))))/(8.*exp(2*(x*x + xbar*xbar))*pow(M_PI,3)*(x - xbar)*(x - xbar));
    
   return z.real();
}

// Gsl function as an integrand for general J integral 
// (when there is no analytical expression for it)
double f_gsl (double epsilon, void * params) 
{
   const complex<double> OMEGA = *(complex<double> *) params;
   double phi_interp = gsl_spline_eval(spline_Phi,epsilon,acc_Phi);
   complex<double> G = 1/(OMEGA - epsilon);
    
   return phi_interp*pow(((-1/M_PI)*G.imag()),n);
}

void load_Phi()
{
   // Opening and reading file that contains epsilon, Phi(epsilon) pairs
   ifstream phi_file(fnPhi.c_str());
   if (!phi_file.is_open()) {
      cerr << "Error opening file " << fnPhi << endl;
      exit(EXIT_FAILURE);
   }
   vector< vector <double> > data;
    
   if(phi_file.is_open()) {
      double vals;
      vector<double> phi_vals;
      while(phi_file >> vals){
	 phi_vals.push_back(vals);
	 if (phi_vals.size() == 2){
	    data.push_back(phi_vals);
	    phi_vals.clear();
	 }
      }
      phi_file.close();
   }
   
   // Making epsilon and phi vectors out of columns from file Phi.dat
   vector<double> eps, phi;

   int N = data.size();
   for (int i=0; i<N; ++i) {
      double x = data[i][0];
      double y = data[i][1];
      y *= pow(x, e); // multiply by epsilon^e
      eps.push_back(x);
      phi.push_back(y);
   }
   
   eps_min = eps[0];
   eps_max = eps[N - 1];
    
   // GSL interpolation of Phi(epsilon) (currently this is hard
   // coded, but akima spline interpolation may be replaced by an
   // arbitrary interpolation type)
   const gsl_interp_type * Interp_type = gsl_interp_akima;
   acc_Phi = gsl_interp_accel_alloc();
   spline_Phi = gsl_spline_alloc(Interp_type, N);
   gsl_spline_init(spline_Phi, &eps[0], &phi[0], N);
}

// General J integral
double J_0n (complex<double> OMEGA) 
{
   // GSL integration to get J(OMEGA)
   gsl_function F;
   F.function = &f_gsl;
   void *params_ptr = &OMEGA;
   F.params = params_ptr;

   const size_t ws_size = 1000;
   
   gsl_integration_workspace *w = gsl_integration_workspace_alloc (ws_size);
   gsl_set_error_handler_off();

   double result_Phi, error_Phi;

   gsl_integration_qag(&F,
		       eps_min,
		       eps_max,
		       abs_error,
		       rel_error,
		       ws_size,
		       key,
		       w,
		       &result_Phi,
		       &error_Phi);
   
   return result_Phi;
}

// Switches between J for different m and n indices.
double J_mn (complex<double> OMEGA)
{
   // Handle the special case first
   if (m == 0) {
      return J_0n(OMEGA);
   }

   switch(n) {
   case 0:
      switch (m) {
      case 1:
	 return J_10(OMEGA);
      case 2:
	 return J_20(OMEGA);
      case 3:
	 return J_30(OMEGA);
      case 4:
	 return J_40(OMEGA);
      case 5:
	 return J_50(OMEGA);
      case 6:
	 return J_60(OMEGA);
      case 7:
	 return J_70(OMEGA);
      case 8:
	 return J_80(OMEGA);
      }

   case 1:
      switch (m) {
      case 1:
	 return J_11(OMEGA);
      case 2:
	 return J_21(OMEGA);
      case 3:
	 return J_31(OMEGA);
      case 4:
	 return J_41(OMEGA);
      case 5:
	 return J_51(OMEGA);
      case 6:
	 return J_61(OMEGA);
      case 7:
	 return J_71(OMEGA);
      case 8:
	 return J_81(OMEGA);
      }
            
   case 2:
      switch (m) {
      case 1:
	 return J_12(OMEGA);
      case 2:
	 return J_22(OMEGA);
      case 3:
	 return J_32(OMEGA);
      case 4:
	 return J_42(OMEGA);
      case 5:
	 return J_52(OMEGA);
      case 6:
	 return J_62(OMEGA);
      case 7:
	 return J_72(OMEGA);
      case 8:
	 return J_82(OMEGA);
      }
      
   case 3:
      switch (m) {
      case 1:
	 return J_13(OMEGA);
      case 2:
	 return J_23(OMEGA);
      case 3:
	 return J_33(OMEGA);
      case 4:
	 return J_43(OMEGA);
      case 5:
	 return J_53(OMEGA);
      case 6:
	 return J_63(OMEGA);
      case 7:
	 return J_73(OMEGA);
      case 8:
	 return J_83(OMEGA);
      }
      
   default:
      cerr << "Jmn not implemented for m=" << m << ", n=" << n << endl;
      abort();
   }
}

// Integrand of Imno integral
double integrand (double omega, void * params) 
{
   double T = *(double *) params;
   double sigma_re, sigma_im;
   if (omega_min <= omega && omega <= omega_max) {
      sigma_re = gsl_spline_eval(spline_reSigma, omega, acc_reSigma);
      sigma_im = gsl_spline_eval(spline_imSigma, omega, acc_imSigma);
   } else {
      sigma_im = -EPSILON;
      if (omega < omega_min)
	 sigma_re = reSigma_asymp_neg;
      if (omega > omega_max)
	 sigma_re = reSigma_asymp_pos;
   }
	 
   complex<double> OMEGA = omega + mu - sigma_re - sigma_im * I;

   double f_factor;
   switch (f_type) {
   case f_derivative:
      f_factor = 1/(T*(2 + 2*cosh(omega/T)));
      break;
   case f_f:
      f_factor = 1/(1+exp(omega/T));
      break;
   default:
      cerr << "Not implemented." << endl;
      exit(EXIT_FAILURE);
   }

   return f_factor * J_mn(OMEGA) * pow(omega,o);
}

// Calculate the epsilon integral only. For m=0 case only.
// For n=1 this is the density of states, which is the case of main
// interest.
void calc_DOS()
{
   assert(m == 0);
   if (n != 1) {
      cout << "Warning: maybe you want n=1 here?" << endl;
   }
   
   ofstream F("dos.dat");
   if (!F.is_open()) {
      cerr << "Error opening file dos.dat" << endl;
      exit(EXIT_FAILURE);
   }
   
   const double range = omega_max - omega_min;
   const int nr_points = 10000;
   const double step = range/nr_points;
   assert(step > 0);
   
   for (int i = 0; i <= nr_points; ++i) {
      const double omega = (i == nr_points ? omega_max : omega_min + i*step);
      double sigma_re, sigma_im;
      if (omega_min <= omega && omega <= omega_max) {
	 sigma_re = gsl_spline_eval(spline_reSigma, omega, acc_reSigma);
	 sigma_im = gsl_spline_eval(spline_imSigma, omega, acc_imSigma);
      } else {
	 sigma_im = -EPSILON;
	 if (omega < omega_min)
	    sigma_re = reSigma_asymp_neg;
	 if (omega > omega_max)
	    sigma_re = reSigma_asymp_pos;
      }
      
      complex<double> OMEGA = omega + mu - sigma_re - sigma_im * I;
      
      double result = J_mn(OMEGA); // -1/Pi Im[] already included in f_gsl() integrand
      
      F << omega << " " << result << endl;
   }
}

void load_Sigma()
{
   const int COLUMNS = 2;
   vector< vector <double> > data1, data2;
    
   ifstream FileSigmaRe (fnReSigma.c_str());
   ifstream FileSigmaIm (fnImSigma.c_str());
    
   if (FileSigmaRe.is_open()) {
      double num_re;
      vector <double> numbers_re;
      while (FileSigmaRe >> num_re) {
	 numbers_re.push_back(num_re);
	 if (numbers_re.size() == COLUMNS) {
	    data1.push_back(numbers_re);
	    numbers_re.clear();
	 }
      }
      FileSigmaRe.close();
   } else {
      cerr << "Error opening file " << fnReSigma << endl;
      abort();
   }
    
   if (FileSigmaIm.is_open()) {
      double num_im;
      vector <double> numbers_im;
      while (FileSigmaIm >> num_im) {
	 numbers_im.push_back(num_im);
	 if (numbers_im.size() == COLUMNS) {
	    data2.push_back(numbers_im);
	    numbers_im.clear();
	 }
      }
      FileSigmaIm.close();
   } else {
      cerr << "Error opening file " << fnImSigma << endl;
      abort();
   }
    
   assert(data1.size() == data2.size()); // sanity check

   int N = data1.size();
   vector <double> omega, reSigma, imSigma;
    
   for(int i=0; i<N; ++i) {
      omega.push_back(data1[i][0]);
      reSigma.push_back(data1[i][1]);
   }
    
   omega_min = data1[0][0];
   omega_max = data1[N-1][0];
   assert(omega_min < omega_max);

   reSigma_asymp_neg = data1[0][1]; // use this for omega<omega_min
   reSigma_asymp_pos = data1[N-1][1]; // use this for omega>omega_max
   // imSigma assumed to be zero outside the [omega_min:omega_max] interval
      
   for(int i=0; i<N; ++i) {
      double val = data2[i][1];
      // IMPORTANT: Perform clipping of Im Sigma !!
      if (val > -EPSILON2)
	 val = -EPSILON2;
      imSigma.push_back(val);
   }
   
   assert(omega_min == data2[0][0]);
   assert(omega_max == data2[N-1][0]);
    
   if (verbose) {
      cout << "omega_min=" << omega_min << " omega_max=" << omega_max << endl;
      cout << "reSigma asymptotic values neg=" << reSigma_asymp_neg << " pos=" << reSigma_asymp_pos << endl;
   }
   
   // Create interpolation objects
   const gsl_interp_type * Interp_type;
   switch(b) {
   case 1:
      Interp_type = gsl_interp_linear;
      break;
   case 2:
      Interp_type = gsl_interp_cspline;
      break;
   case 3:
      Interp_type = gsl_interp_akima;
      break;
   default:
      cerr << "Interpolation " << b << " is not implemented." << endl;
      abort();
   }
   
   // GSL interpolation with linear interpolation as default
   acc_reSigma = gsl_interp_accel_alloc();
   spline_reSigma  = gsl_spline_alloc(Interp_type, N);
   gsl_spline_init(spline_reSigma, &omega[0], &reSigma[0], N);
    
   acc_imSigma = gsl_interp_accel_alloc();
   spline_imSigma = gsl_spline_alloc(Interp_type, N);
   gsl_spline_init(spline_imSigma, &omega[0], &imSigma[0], N);
}

void calc()
{
   // GSL integration
   const size_t ws_size = 1000;
   gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (ws_size);

   double lower_limit;
   switch (f_type) {
   case f_derivative:
      lower_limit = -cutoff*T;
      break;
   case f_f:
      lower_limit = omega_min;
      break;
   default:
      cerr << "Not implemented." << endl;
      exit(EXIT_FAILURE);
   }

   double upper_limit = cutoff*T;
   if (!(lower_limit < upper_limit)) {
      cerr << "Invalid integration interval [" << lower_limit << ", " << upper_limit << "]." << endl;
      exit(EXIT_FAILURE);
   }
   double result;
   double error;
    
   gsl_function My_function;
   My_function.function = &integrand;
   void *params_ptr = &T;
   My_function.params = params_ptr;
    
   gsl_set_error_handler_off();
   gsl_integration_qag (&My_function,       // integrand function
			lower_limit,        // lower integration boundary
			upper_limit,        // upper integration boundary
			abs_error,          // preferred absolute error
			rel_error,          // preferred relative error
			ws_size,               // size of workspace
			key,                // Gauss-Kronrod rule
			work_ptr,           // integration workspace
			&result,            // final approximation
			&error);            // estimate of absolute error
    
   const int width = 20;
   const int precision_result = 16;
   const int precision_error = 4;
   
   if (verbose) {
      cout << "Result = " << setprecision(precision_result) << setw(width) << result << endl;
      cout << "Error  = " << setprecision(precision_error) << setw(width) << error << endl;
   } else {
      cout << setprecision(precision_result) << result << endl;
   }
}

int main (int argc, char *argv[]) 
{
   cmd_line(argc, argv);

   load_Sigma();
   if (m == 0)
      load_Phi();
   
   if (calcdos)
      calc_DOS();
   else
      calc();
    
   return 0;
}
