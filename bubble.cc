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
// 7.8.2026 - support finite-frequency optical conductivity through -O

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
#include <cerrno>
#include <cmath>
#include <complex>
#include <limits>

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

const string VERSION = "1.5";

// Mandatory parameters
int m, n, o;
double T;
double mu;
string fnReSigma, fnImSigma;

// Optional parameters (with defaults)
int b = 1;
bool verbose = false;
bool quiet_warnings = false;
int key = GSL_INTEG_GAUSS15;
double abs_error = 1.0e-7;
double rel_error = 1.0e-8;
double cutoff = 15.0;
double sigma_clip = 1.0e-8;
string fnPhi = "Phi.dat";
bool calcdos = false; // compute the spectral function
enum ff { f_f, f_derivative };
ff f_type = f_derivative;
int e = 0; // power of epsilon
bool optical_mode = false;
double optical_frequency = 0.0;

const double EPSILON = 1e-10; // some very small value...

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
   cout << "-q : suppress non-fatal warnings" << endl;
   cout << "-i : interpolation (default = 1, 1=>linear, 2=>cspline, 3=>Akima spline)" << endl;
   cout << "-k : integration rule (default = 1, 1 => 15, 2 => 21, etc.)" << endl;
   cout << "-a : absolute error (default = 1e-7)" << endl;
   cout << "-r : relative error (default = 1e-8)" << endl;
   cout << "-c : frequency interval cutoff in units of T (default = 15)" << endl;
   cout << "-s FLOOR : positive ImSigma clipping floor (default = 1e-8)" << endl;
   cout << "-p : filename for Phi tables (default = Phi.dat)" << endl;
   cout << "-d : compute the epsilon integrals only, for m=0 (output = dos.dat)" << endl;
   cout << "-e : additional power of epsilon when using the m=0 code (default=0)" << endl;
   cout << "-f : switch (-df/dw) to f in the w integration (incompatible with -d)" << endl;
   cout << "-O OMEGA : external optical frequency (requires OMEGA>=0 and n=2)" << endl;
}

double parse_optical_frequency(const char *value)
{
   char *end = NULL;
   errno = 0;
   const double result = strtod(value, &end);

   if (errno == ERANGE || end == value || *end != '\0' || !isfinite(result) || result < 0.0) {
      cerr << "Invalid optical frequency: " << value << endl;
      exit(EXIT_FAILURE);
   }

   return result;
}

double parse_sigma_clip(const char *value)
{
   char *end = NULL;
   errno = 0;
   const double result = strtod(value, &end);

   if (errno == ERANGE || end == value || *end != '\0' || !isfinite(result) || result <= 0.0) {
      cerr << "Invalid ImSigma clipping floor: " << value << endl;
      exit(EXIT_FAILURE);
   }

   return result;
}

void cmd_line(int argc, char *argv[])
{
   int c;
    
   while ((c = getopt(argc, argv, "vqi:k:a:r:c:s:p:dfe:O:")) != -1) {
      switch (c) {
      case 'v':
	 verbose = true;
	 break;
      case 'q':
	 quiet_warnings = true;
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
      case 's':
	 sigma_clip = parse_sigma_clip(optarg);
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
      case 'O':
	 optical_mode = true;
	 optical_frequency = parse_optical_frequency(optarg);
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
   if (optical_mode && n != 2) {
      cerr << "Option -O requires n=2." << endl;
      exit(EXIT_FAILURE);
   }
   if (optical_mode && calcdos) {
      cerr << "Options -O and -d cannot be used together." << endl;
      exit(EXIT_FAILURE);
   }
   if (optical_mode && f_type == f_f) {
      cerr << "Options -O and -f cannot be used together." << endl;
      exit(EXIT_FAILURE);
   }
   if (optical_mode && (!isfinite(T) || T <= 0.0 || !isfinite(cutoff) || cutoff <= 0.0)) {
      cerr << "Optical conductivity requires positive finite T and cutoff." << endl;
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
      cout << "ImSigma clipping floor=" << sigma_clip << endl;
      cout << "Phi=" << fnPhi << " e=" << e << endl;
      if (optical_mode)
	 cout << "external Omega=" << optical_frequency << endl;
      if (calcdos) {
	 cout << "epsilon integrals only (dos.dat)" << endl;
      } else if (optical_mode && optical_frequency > 0.0) {
	 cout << "[f(omega)-f(omega+Omega)]/Omega" << endl;
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

double J_builtin(int kernel, int power, complex<double> z);

// Functions Jmn, for m = 1,...,8 and n = 0,1,2,3
// The historical generated n=2,3 expressions are retained with a _generated
// suffix for provenance; production dispatch uses the stable functions below.
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


double J_12_generated(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return ((2*y*(1 + pow(y,2) - pow(x,2)) + (pow(y,2) + pow(1 - x,2)) * (pow(y,2) + pow(1 + x,2))*(atan((1 - x)/y) + atan((1 + x)/y))))/(2.*y*pow(M_PI,2)*(pow(y,2) + pow(1 - x,2))*(pow(y,2) + pow(1 + x,2)));
}

double J_13_generated(complex<double> OMEGA) {
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


double J_22_generated(complex<double> OMEGA) {
   double x = OMEGA.real();
   double y = OMEGA.imag();
   return  (x*(-2*y*(-1 + pow(x,2) + pow(y,2)) - (pow(-1 + x,2) + pow(y,2))*(pow(1 + x,2) + pow(y,2))*atan((-1 + x)/y) +(pow(-1 + x,2) + pow(y,2))*(pow(1 + x,2) + pow(y,2))*atan((1 + x)/y)))/(2.*pow(M_PI,2)*y*(pow(-1 + x,2) + pow(y,2))*(pow(1 + x,2) + pow(y,2)));
}

double J_23_generated(complex<double> OMEGA) {
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

double J_32_generated(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = -(-2 - (I * x * s)/sqrt(1 - x*x) + (I * xbar * s)/sqrt(1 - xbar*xbar)- (2 * (- x + xbar + I * (sqrt(1 - x*x) + sqrt(1 - xbar*xbar)) * s)/(x-xbar)))/(4.*M_PI);
   return z.real();
}

double J_33_generated(complex<double> OMEGA) {
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

double J_42_generated(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = -((I*(2*x*sqrt(1 - x*x)*xbar*xbar - x*(sqrt(1-x*x) + sqrt(1 - xbar*xbar)) - xbar*(sqrt(1-x*x) + sqrt(1-xbar*xbar) - 2*x*x*sqrt(1-xbar*xbar)))*s)/(4*M_PI*sqrt(1 - x*x)*(x-xbar)*sqrt(1 - xbar*xbar)));
    
   return z.real();
}

double J_43_generated(complex<double> OMEGA) {
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

double J_52_generated(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = (-pow(x,3) + pow(xbar,3) + I*(x*x*sqrt(1 - x*x) + 2*(sqrt(1 - x*x) + sqrt(1 - xbar*xbar)))*s + xbar*xbar*(-3*x + I*sqrt(1 - xbar*xbar)*s) + 3*x*xbar*(x - I*(sqrt(1 - x*x) + sqrt(1 - xbar*xbar))*s))/(4.*M_PI*(x - xbar));
   
   return z.real();
}

double J_53_generated(complex<double> OMEGA) {
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

double J_62_generated(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z =  -(x*(-3 + 4*x*x) - 3*xbar + 4*pow(xbar,3) - I*sqrt(1 - x*x)*(-1 + 4*x*x)*s - I*sqrt(1 - xbar*xbar)*s + 4*I*xbar*xbar*sqrt(1 - xbar*xbar)*s - (I*(3*I*x*x - 2*I*pow(x,4) - 3*I*xbar*xbar + 2*I*pow(xbar,4) + 2*x*pow(1 - x*x,1.5)*s + 2*xbar*pow(1 - xbar*xbar,1.5)*s))/(x - xbar))/(4.*M_PI);
   
   return z.real();
}

double J_63_generated(complex<double> OMEGA) {
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

double J_72_generated(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = (2*exp(2*(x*x + xbar*xbar))*x*sqrt(2*M_PI) - exp(2*xbar*xbar)*(1 + 2*x*x)*M_PI*Erfi(sqrt(2)*x) + exp(2*x*x)*M_PI*Erfi(sqrt(2)*xbar) + exp(2*xbar*xbar)*log(-(1/x)) + 2*exp(2*xbar*xbar)*x*x*log(-(1/x)) - exp(2*xbar*xbar)*log(1/x) - 2*exp(2*xbar*xbar)*x*x*log(1/x) - exp(2*x*x)*log(-(1/xbar)) + exp(2*x*x)*log(1/xbar) + 2*exp(2*x*x)*xbar*xbar*(M_PI*Erfi(sqrt(2)*xbar) - log(-(1/xbar)) + log(1/xbar)) - 2*xbar*(exp(2*(x*x + xbar*xbar))*sqrt(2*M_PI) - exp(2*xbar*xbar)*x*M_PI*Erfi(sqrt(2)*x) + exp(2*x*x)*x*M_PI*Erfi(sqrt(2)*xbar) + exp(2*xbar*xbar)*x*log(-(1/x)) - exp(2*xbar*xbar)*x*log(1/x) - exp(2*x*x)*x*log(-(1/xbar)) + exp(2*x*x)*x*log(1/xbar)))/(2.*exp(2*(x*x + xbar*xbar))*M_PI*M_PI*(x - xbar));
   
   return z.real();
}

double J_73_generated(complex<double> OMEGA) {
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

double J_82_generated(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = -(-4*exp(2*x*x)*pow(xbar,3)*(M_PI*Erfi(sqrt(2)*xbar) - log(-(1/xbar)) - log(xbar)) + x*(-2*exp(2*(x*x + xbar*xbar))*x*sqrt(2*M_PI) + exp(2*xbar*xbar)*(1 + 4*x*x)*M_PI*Erfi(sqrt(2)*x) - exp(2*x*x)*M_PI*Erfi(sqrt(2)*xbar) - exp(2*xbar*xbar)*log(-(1/x)) - 4*exp(2*xbar*xbar)*x*x*log(-(1/x)) - exp(2*xbar*xbar)*log(x) - 4*exp(2*xbar*xbar)*x*x*log(x) + exp(2*x*x)*log(-(1/xbar)) + exp(2*x*x)*log(xbar)) + xbar*(-(exp(2*xbar*xbar)*(-1 + 4*x*x)*M_PI*Erfi(sqrt(2)*x)) - exp(2*x*x)*M_PI*Erfi(sqrt(2)*xbar) - exp(2*xbar*xbar)*log(-(1/x)) + 4*exp(2*xbar*xbar)*x*x*log(-(1/x)) - exp(2*xbar*xbar)*log(x) + 4*exp(2*xbar*xbar)*x*x*log(x) + exp(2*x*x)*log(-(1/xbar)) + exp(2*x*x)*log(xbar)) + 2*exp(2*x*x)*xbar*xbar*(exp(2*xbar*xbar)*sqrt(2*M_PI) + 2*x*M_PI*Erfi(sqrt(2)*xbar) - 2*x*log(-(1/xbar)) - 2*x*log(xbar)))/(4.*exp(2*(x*x + xbar*xbar))*M_PI*M_PI*(x - xbar));
   
   return z.real();
}

double J_83_generated(complex<double> OMEGA) {
   complex<double> x = OMEGA;
   complex<double> xbar = conj(OMEGA);
   double s = (OMEGA.imag() > 0.0 ? 1.0 : -1.0);
   complex<double> z = (-I*(2*exp(2*xbar*xbar)*xbar*xbar*(3*exp(2*x*x)*sqrt(2*M_PI) + x*(-3 + 4*x*x)*M_PI*Erfi(sqrt(2)*x) + (3*x - 4*pow(x,3))*log(-(1/x)) + 3*x*log(x) - 4*pow(x,3)*log(x)) - 8*exp(2*x*x)*pow(xbar,5)*(M_PI*Erfi(sqrt(2)*xbar) - log(-(1/xbar)) - log(xbar)) - x*(6*exp(2*(x*x + xbar*xbar))*x*sqrt(2*M_PI) + 4*exp(2*(x*x + xbar*xbar))*pow(x,3)*sqrt(2*M_PI) - exp(2*xbar*xbar)*(3 + 6*x*x + 8*pow(x,4))*M_PI*Erfi(sqrt(2)*x) + 3*exp(2*x*x)*M_PI*Erfi(sqrt(2)*xbar) + 3*exp(2*xbar*xbar)*log(-(1/x)) + 6*exp(2*xbar*xbar)*x*x*log(-(1/x)) + 8*exp(2*xbar*xbar)*pow(x,4)*log(-(1/x)) + 3*exp(2*xbar*xbar)*log(x) + 6*exp(2*xbar*xbar)*x*x*log(x) + 8*exp(2*xbar*xbar)*pow(x,4)*log(x) - 3*exp(2*x*x)*log(-(1/xbar)) - 3*exp(2*x*x)*log(xbar)) + 4*exp(2*x*x)*pow(xbar,4)*(exp(2*xbar*xbar)*sqrt(2*M_PI) + 4*x*M_PI*Erfi(sqrt(2)*xbar) - 4*x*log(-(1/xbar)) - 4*x*log(xbar)) - 2*exp(2*x*x)*pow(xbar,3)*(4*exp(2*xbar*xbar)*x*sqrt(2*M_PI) + (3 + 4*x*x)*M_PI*Erfi(sqrt(2)*xbar) - (3 + 4*x*x)*log(-(1/xbar)) - 3*log(xbar) - 4*x*x*log(xbar)) + xbar*(8*exp(2*(x*x + xbar*xbar))*pow(x,3)*sqrt(2*M_PI) - exp(2*xbar*xbar)*(-3 + 16*pow(x,4))*M_PI*Erfi(sqrt(2)*x) + 3*exp(2*x*x)*(-1 + 2*x*x)*M_PI*Erfi(sqrt(2)*xbar) - 3*exp(2*xbar*xbar)*log(-(1/x)) + 16*exp(2*xbar*xbar)*pow(x,4)*log(-(1/x)) - 3*exp(2*xbar*xbar)*log(x) + 16*exp(2*xbar*xbar)*pow(x,4)*log(x) + 3*exp(2*x*x)*log(-(1/xbar)) - 6*exp(2*x*x)*x*x*log(-(1/xbar)) + 3*exp(2*x*x)*log(xbar) - 6*exp(2*x*x)*x*x*log(xbar))))/(8.*exp(2*(x*x + xbar*xbar))*pow(M_PI,3)*(x - xbar)*(x - xbar));
    
   return z.real();
}

double J_12(complex<double> z) { return J_builtin(1, 2, z); }
double J_13(complex<double> z) { return J_builtin(1, 3, z); }
double J_22(complex<double> z) { return J_builtin(2, 2, z); }
double J_23(complex<double> z) { return J_builtin(2, 3, z); }
double J_32(complex<double> z) { return J_builtin(3, 2, z); }
double J_33(complex<double> z) { return J_builtin(3, 3, z); }
double J_42(complex<double> z) { return J_builtin(4, 2, z); }
double J_43(complex<double> z) { return J_builtin(4, 3, z); }
double J_52(complex<double> z) { return J_builtin(5, 2, z); }
double J_53(complex<double> z) { return J_builtin(5, 3, z); }
double J_62(complex<double> z) { return J_builtin(6, 2, z); }
double J_63(complex<double> z) { return J_builtin(6, 3, z); }
double J_72(complex<double> z) { return J_builtin(7, 2, z); }
double J_73(complex<double> z) { return J_builtin(7, 3, z); }
double J_82(complex<double> z) { return J_builtin(8, 2, z); }
double J_83(complex<double> z) { return J_builtin(8, 3, z); }

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

struct optical_kernel_params {
   complex<double> z1;
   complex<double> z2;
};

double f_gsl_optical(double epsilon, void *params)
{
   const optical_kernel_params &p = *(optical_kernel_params *)params;
   const double phi_interp = gsl_spline_eval(spline_Phi, epsilon, acc_Phi);
   const complex<double> G1 = 1.0/(p.z1 - epsilon);
   const complex<double> G2 = 1.0/(p.z2 - epsilon);
   const double A1 = -G1.imag()/M_PI;
   const double A2 = -G2.imag()/M_PI;

   return phi_interp*A1*A2;
}

void add_optical_peak_breakpoints(vector<double> &points, complex<double> z)
{
   const double center = z.real();
   const double width = abs(z.imag());
   const double range = eps_max - eps_min;
   if (!isfinite(center) || !isfinite(width) || width == 0.0 ||
       (center < eps_min && eps_min - center > width) ||
       (center > eps_max && center - eps_max > width))
      return;

   if (eps_min < center && center < eps_max)
      points.push_back(center);
   double distance = width;
   for (int scale = 0; scale < 24; ++scale) {
      const double lower_point = center - distance;
      const double upper_point = center + distance;
      if (eps_min < lower_point && lower_point < eps_max)
	 points.push_back(lower_point);
      if (eps_min < upper_point && upper_point < eps_max)
	 points.push_back(upper_point);
      if (distance >= range)
	 break;
      distance *= 8.0;
   }
}

double J_0_optical(complex<double> z1, complex<double> z2)
{
   optical_kernel_params params = {z1, z2};
   gsl_function F;
   F.function = &f_gsl_optical;
   F.params = &params;

   const size_t ws_size = 1000;
   gsl_integration_workspace *work = gsl_integration_workspace_alloc(ws_size);
   if (work == NULL) {
      cerr << "Unable to allocate optical kernel integration workspace." << endl;
      exit(EXIT_FAILURE);
   }

   double result;
   double error;
   vector<double> points;
   points.push_back(eps_min);
   points.push_back(eps_max);
   add_optical_peak_breakpoints(points, z1);
   add_optical_peak_breakpoints(points, z2);
   sort(points.begin(), points.end());
   const double spacing = min(0.25*(eps_max - eps_min),
      64.0*numeric_limits<double>::epsilon()*
      max(1.0, max(abs(eps_min), abs(eps_max))));
   vector<double> filtered_points;
   filtered_points.push_back(eps_min);
   for (vector<double>::const_iterator point = points.begin(); point != points.end(); ++point) {
      if (eps_min < *point && *point < eps_max &&
	  *point - filtered_points.back() > spacing && eps_max - *point > spacing)
	 filtered_points.push_back(*point);
   }
   filtered_points.push_back(eps_max);
   points.swap(filtered_points);

   const int status = gsl_integration_qagp(&F,
					   &points[0],
					   points.size(),
					   abs_error,
					   rel_error,
					   ws_size,
					   work,
					   &result,
					   &error);
   gsl_integration_workspace_free(work);

   if ((status != GSL_SUCCESS && status != GSL_EROUND) || !isfinite(result)) {
      cerr << "Optical epsilon integration failed for z1=" << z1 << ", z2=" << z2
	   << ": " << gsl_strerror(status) << endl;
      exit(EXIT_FAILURE);
   }
   if (status == GSL_EROUND && !quiet_warnings) {
      static bool warning_emitted = false;
      if (!warning_emitted) {
	 cerr << "Warning: optical epsilon integration was limited by roundoff." << endl;
	 warning_emitted = true;
      }
   }

   return result;
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
	 break;

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
	 break;
            
   case 2:
   case 3:
      if (1 <= m && m <= 8)
	 return J_builtin(m, n, OMEGA);
      break;
      
   default:
      break;
   }

   cerr << "Jmn not implemented for m=" << m << ", n=" << n << endl;
   abort();
}

struct hilbert_pair {
   complex<double> value;
   complex<double> derivative;
   complex<double> second_derivative;
};

// H_m(z) = integral Phi_m(epsilon)/(z-epsilon) d epsilon.

hilbert_pair flat_hilbert(int kernel, complex<double> z)
{
   if (abs(z) > 2.0) {
      const complex<double> u = 1.0/z;
      const complex<double> u2 = u*u;

      if (kernel == 1) {
	 complex<double> value = 0.0;
	 complex<double> term = u;
	 for (int k = 0; k < 64; ++k, term *= u2)
	    value += 2.0*term/(2*k + 1.0);

	 return hilbert_pair{value,
	    -2.0*u2/(1.0 - u2),
	    4.0*u*u2/((1.0 - u2)*(1.0 - u2))};
      }

      complex<double> value = 0.0;
      complex<double> derivative = 0.0;
      complex<double> second_derivative = 0.0;
      complex<double> term = u2;
      for (int k = 1; k < 64; ++k, term *= u2) {
	 value += 2.0*term/(2*k + 1.0);
	 derivative -= 4.0*k*term*u/(2*k + 1.0);
	 second_derivative += 4.0*k*term*u2;
      }
      return hilbert_pair{value, derivative, second_derivative};
   }

   const complex<double> value1 = log(z + 1.0) - log(z - 1.0);
   const complex<double> derivative1 = 2.0/(1.0 - z*z);
   const complex<double> second_derivative1 =
      4.0*z/((1.0 - z*z)*(1.0 - z*z));
   if (kernel == 1)
      return hilbert_pair{value1, derivative1, second_derivative1};

   return hilbert_pair{z*value1 - 2.0,
		       value1 + z*derivative1,
		       2.0*derivative1 + z*second_derivative1};
}

hilbert_pair gaussian_hilbert_asymptotic(int kernel, complex<double> z)
{
   const double mass = sqrt(M_PI/2.0);
   const complex<double> u = 1.0/z;
   const complex<double> u2 = u*u;
   complex<double> power = 1.0;
   complex<double> value7 = 0.0;
   complex<double> derivative7 = 0.0;
   complex<double> second_derivative7 = 0.0;
   complex<double> value8 = 0.0;
   complex<double> derivative8 = 0.0;
   complex<double> second_derivative8 = 0.0;
   double moment = mass;

   for (int k = 0; k < 32; ++k) {
      const complex<double> term = moment*power;
      value7 += term*u;
      derivative7 -= (2*k + 1.0)*term*u*u;
      second_derivative7 += (2*k + 1.0)*(2*k + 2.0)*term*u*u*u;

      if (k != 0) {
	 value8 += term;
	 derivative8 -= 2.0*k*term*u;
	 second_derivative8 += 2.0*k*(2*k + 1.0)*term*u*u;
      }

      power *= u2;
      moment *= (2*k + 1.0)/4.0;
   }

   if (kernel == 7)
      return hilbert_pair{value7, derivative7, second_derivative7};
   return hilbert_pair{value8, derivative8, second_derivative8};
}

bool lower_half_plane(complex<double> z)
{
   return z.imag() < 0.0 || (z.imag() == 0.0 && signbit(z.imag()));
}

hilbert_pair kernel_hilbert(int kernel, complex<double> z)
{
   if (kernel == 1 || kernel == 2)
      return flat_hilbert(kernel, z);

   if (3 <= kernel && kernel <= 6) {
      const complex<double> root = sqrt(z - 1.0)*sqrt(z + 1.0);
      const complex<double> t = 1.0/(z + root);

      if (kernel == 3)
	 return hilbert_pair{M_PI*t, -M_PI*t/root, M_PI/(root*root*root)};
      if (kernel == 4)
	 return hilbert_pair{0.5*M_PI*t*t,
			     -M_PI*t*t/root,
			     M_PI*t*t*(z + 2.0*root)/(root*root*root)};
      if (kernel == 5)
	 return hilbert_pair{0.25*M_PI*t*(3.0 - t*t),
			     -1.5*M_PI*t*t,
			     3.0*M_PI*t*t/root};
      return hilbert_pair{0.125*M_PI*t*t*(2.0 - t*t),
			  -M_PI*t*t*t,
			  3.0*M_PI*t*t*t/root};
   }

   if (kernel == 7 || kernel == 8) {
      if (abs(z) > 8.0)
	 return gaussian_hilbert_asymptotic(kernel, z);

      const double sqrt2 = sqrt(2.0);
      const double mass = sqrt(M_PI/2.0);
      const complex<double> value7 = lower_half_plane(z)
	 ? I*M_PI*Faddeeva::w(-sqrt2*z)
	 : -I*M_PI*Faddeeva::w(sqrt2*z);
      const complex<double> derivative7 = 4.0*mass - 4.0*z*value7;
      const complex<double> second_derivative7 = -4.0*value7 - 4.0*z*derivative7;

      if (kernel == 7)
	 return hilbert_pair{value7, derivative7, second_derivative7};
      return hilbert_pair{z*value7 - mass,
			  value7 + z*derivative7,
			  2.0*derivative7 + z*second_derivative7};
   }

   cerr << "Optical kernel not implemented for m=" << kernel << endl;
   exit(EXIT_FAILURE);
}

namespace {

const long double PI_L = acosl(-1.0L);

void compensated_add(long double term, long double &sum, long double &correction)
{
   const long double adjusted = term - correction;
   const long double updated = sum + adjusted;
   correction = (updated - sum) - adjusted;
   sum = updated;
}

long double finite_moment_series(int kernel, int power, long double x, long double y)
{
   // Expanding J itself avoids cancellation between Hilbert-transform terms
   // when the Lorentzian center is far from a finite band.
   const long double radius = hypotl(x, y);
   const long double inverse_radius = 1.0L/radius;
   const long double cosine = x/radius;
   long double gegenbauer_previous = 0.0L;
   long double gegenbauer = 1.0L;
   long double inverse_power = 1.0L;
   long double moment = 0.0L;

   if (kernel == 3)
      moment = PI_L/2.0L;
   else if (kernel == 4)
      moment = PI_L/8.0L;
   else if (kernel == 5)
      moment = 3.0L*PI_L/8.0L;
   else if (kernel == 6)
      moment = PI_L/16.0L;

   long double sum = 0.0L;
   long double correction = 0.0L;
   int negligible_terms = 0;
   for (int k = 0; k <= 160; ++k) {
      long double kernel_moment = 0.0L;
      if (kernel == 1 && k % 2 == 0) {
         kernel_moment = 2.0L/(k + 1.0L);
      } else if (kernel == 2 && k % 2 == 1) {
         kernel_moment = 2.0L/(k + 2.0L);
      } else if (kernel == 3 && k % 2 == 0) {
         kernel_moment = moment;
         moment *= (k + 1.0L)/(k + 4.0L);
      } else if (kernel == 4 && k % 2 == 1) {
         kernel_moment = moment;
         moment *= (k + 2.0L)/(k + 5.0L);
      } else if (kernel == 5 && k % 2 == 0) {
         kernel_moment = moment;
         moment *= (k + 1.0L)/(k + 6.0L);
      } else if (kernel == 6 && k % 2 == 1) {
         kernel_moment = moment;
         moment *= (k + 2.0L)/(k + 7.0L);
      }

      if (kernel_moment != 0.0L) {
         const long double term = kernel_moment*gegenbauer*inverse_power;
         compensated_add(term, sum, correction);
         if (k > 12 && fabsl(term) <=
             4.0L*numeric_limits<double>::epsilon()*fabsl(sum))
            ++negligible_terms;
         else
            negligible_terms = 0;
         if (negligible_terms == 3)
            break;
      }

      const long double next_gegenbauer =
         (2.0L*(k + power)*cosine*gegenbauer
          - (k + 2.0L*power - 1.0L)*gegenbauer_previous)/(k + 1.0L);
      gegenbauer_previous = gegenbauer;
      gegenbauer = next_gegenbauer;
      inverse_power *= inverse_radius;
   }

   const long double scaled_y = y/(PI_L*radius);
   return powl(scaled_y, power)*sum/powl(radius, power);
}

long double flat_exterior_kernel(int kernel, int power, long double x, long double y)
{
   const long double distance = x - 1.0L;
   const long double eta = y/distance;
   const long double eta2 = eta*eta;
   const long double ratio = distance/(x + 1.0L);
   const long double ratio2 = ratio*ratio;
   int coefficient_index = 2*power - 1;
   long double ratio_power = powl(ratio, coefficient_index);
   long double ratio_previous = ratio_power/ratio;
   long double binomial = 1.0L;
   long double eta_power = 1.0L;
   long double sign = 1.0L;
   long double sum = 0.0L;
   long double correction = 0.0L;

   for (int j = 0; j < 96; ++j) {
      const long double h1 = (ratio_power - 1.0L)/coefficient_index;
      long double h = h1;
      if (kernel == 2) {
         const long double h1_previous =
            (1.0L - ratio_previous)/(coefficient_index - 1.0L);
         h = x*h1 + distance*h1_previous;
      }

      const long double term = sign*binomial*(-h)*eta_power;
      compensated_add(term, sum, correction);
      if (j > 8 && fabsl(term) <=
          4.0L*numeric_limits<long double>::epsilon()*fabsl(sum))
         break;

      binomial *= (power + j)/(j + 1.0L);
      eta_power *= eta2;
      sign = -sign;
      coefficient_index += 2;
      ratio_power *= ratio2;
      ratio_previous *= ratio2;
   }

   return powl(eta, power)*sum/
      (powl(PI_L, power)*powl(distance, power - 1));
}

long double flat_kernel(int kernel, int power, long double x, long double y)
{
   const long double y2 = y*y;
   const long double lower_distance = 1.0L - x;
   const long double upper_distance = 1.0L + x;
   const long double lower_denominator = lower_distance*lower_distance + y2;
   const long double upper_denominator = upper_distance*upper_distance + y2;
   const long double a1 = lower_distance/lower_denominator
      + upper_distance/upper_denominator;
   const long double theta = atan2l(2.0L*y, (x - 1.0L)*(x + 1.0L) + y2);
   const long double j12 = (theta + y*a1)/(2.0L*PI_L*PI_L*y);

   if (power == 2) {
      if (kernel == 1)
         return j12;
      return x*(j12 - 2.0L*y2/
         (PI_L*PI_L*lower_denominator*upper_denominator));
   }

   const long double a2 = lower_distance/(lower_denominator*lower_denominator)
      + upper_distance/(upper_denominator*upper_denominator);
   const long double j13 = (3.0L*theta + 3.0L*y*a1 + 2.0L*y*y2*a2)/
      (8.0L*PI_L*PI_L*PI_L*y2);
   if (kernel == 1)
      return j13;
   return x*(j13 - y*y2*(lower_denominator + upper_denominator)/
      (PI_L*PI_L*PI_L*lower_denominator*lower_denominator*
       upper_denominator*upper_denominator));
}

long double square_root_kernel(int kernel, int power, long double x, long double y)
{
   // These branch variables keep every numerator nonnegative for x >= 0 and
   // avoid subtracting nearly equal square roots at and outside the band edge.
   const long double y2 = y*y;
   const long double band_coordinate = (x - 1.0L)*(x + 1.0L);
   const long double radius = hypotl(band_coordinate - y2, 2.0L*x*y);
   const long double sum_coordinate = band_coordinate + y2;
   const long double a = sum_coordinate >= 0.0L
      ? 2.0L*y2/(radius + sum_coordinate)
      : (radius - sum_coordinate)/2.0L;
   const long double q = sqrtl(y2 + a);
   const long double g = a/(q + y);
   const long double scale = a*a + y2;
   long double numerator;

   if (power == 2) {
      if (kernel == 3)
         numerator = q*a*a/scale;
      else if (kernel == 4)
         numerator = x*a*a*a/(q*scale);
      else if (kernel == 5)
         numerator = g*g*(q + 2.0L*y);
      else
         numerator = x*g*g*g*(q + 3.0L*y)/q;
      return numerator/(2.0L*PI_L*y);
   }

   if (kernel == 3) {
      numerator = q*a*a*a*(3.0L*a*a*a + 7.0L*a*y2 + 4.0L*y2*y2)/
         (scale*scale*scale);
   } else if (kernel == 4) {
      numerator = x*a*a*a*a*
         (3.0L*a*a*a + 2.0L*a*a*y2 + 7.0L*a*y2 + 6.0L*y2*y2)/
         (q*scale*scale*scale);
   } else if (kernel == 5) {
      numerator = 3.0L*q*a*a*a/scale;
   } else {
      numerator = 3.0L*x*a*a*a*a/(q*scale);
   }
   return numerator/(8.0L*PI_L*PI_L*y2);
}

bool gaussian_tail_kernel(int kernel, int power, long double x, long double y,
                          long double &result)
{
   // The algebraic moment series misses an exponentially small near-real pole;
   // in the narrow exterior sector both contributions are needed.
   const long double radius2 = x*x + y*y;
   const long double alpha = -2.0L*x/radius2;
   const long double beta = 1.0L/radius2;
   long double coefficient_previous = 0.0L;
   long double coefficient = 1.0L;
   long double moment = sqrtl(PI_L/2.0L);
   if (kernel == 8)
      moment /= 4.0L;

   long double sum = 0.0L;
   long double correction = 0.0L;
   long double block_maximum = 0.0L;
   long double previous_block_maximum = numeric_limits<long double>::infinity();
   int block_terms = 0;
   int rising_blocks = 0;
   bool converged = false;

   for (int k = 0; k <= 255; ++k) {
      if ((kernel == 7 && k % 2 == 0) || (kernel == 8 && k % 2 == 1)) {
         const long double term = moment*coefficient;
         compensated_add(term, sum, correction);
         block_maximum = max(block_maximum, fabsl(term));
         ++block_terms;

         if (kernel == 7)
            moment *= (k + 1.0L)/4.0L;
         else
            moment *= (k + 2.0L)/4.0L;

         if (block_terms == 4) {
            if (block_maximum <= 8.0L*numeric_limits<double>::epsilon()*fabsl(sum)) {
               converged = true;
               break;
            }
            if (block_maximum > previous_block_maximum)
               ++rising_blocks;
            else
               rising_blocks = 0;
            if (rising_blocks >= 2)
               return false;
            previous_block_maximum = block_maximum;
            block_maximum = 0.0L;
            block_terms = 0;
         }
      }

      const long double next_coefficient =
         -(alpha*(k + power)*coefficient
           + beta*(k + 2.0L*power - 1.0L)*coefficient_previous)/(k + 1.0L);
      coefficient_previous = coefficient;
      coefficient = next_coefficient;
   }

   if (!converged)
      return false;

   const long double radius = sqrtl(radius2);
   result = powl(y/(PI_L*radius), power)*sum/powl(radius, power);

   if (y < x && x*y <= 0.25L) {
      const complex<long double> z(x, y);
      const complex<long double> gaussian = exp(-2.0L*z*z);
      complex<long double> psi;
      complex<long double> psi_derivative;
      complex<long double> psi_second_derivative;
      if (kernel == 7) {
         psi = gaussian;
         psi_derivative = -4.0L*z*gaussian;
         psi_second_derivative = (16.0L*z*z - 4.0L)*gaussian;
      } else {
         psi = z*gaussian;
         psi_derivative = (1.0L - 4.0L*z*z)*gaussian;
         psi_second_derivative = (16.0L*z*z*z - 12.0L*z)*gaussian;
      }

      if (power == 2) {
         result += (psi.real()/y + psi_derivative.imag())/(2.0L*PI_L);
      } else {
         result += (3.0L*psi.real()/(y*y) + 3.0L*psi_derivative.imag()/y
                    - psi_second_derivative.real())/(8.0L*PI_L*PI_L);
      }
   }

   return isfinite(result) && result >= 0.0L;
}

bool gaussian_hilbert_kernel(int kernel, int power, double x, double y,
                             long double &result)
{
   const complex<double> z(x, y);
   const double mass = sqrt(M_PI/2.0);
   const complex<double> h7 = -I*M_PI*Faddeeva::w(sqrt(2.0)*z);
   const complex<double> h7_derivative = 4.0*mass - 4.0*z*h7;
   const complex<double> h7_second_derivative = -4.0*h7 - 4.0*z*h7_derivative;
   complex<double> h = h7;
   complex<double> derivative = h7_derivative;
   complex<double> second_derivative = h7_second_derivative;

   const double error0 = abs(h7);
   const double error1 = 4.0*mass + 4.0*abs(z)*error0;
   const double error2 = 4.0*error0 + 4.0*abs(z)*error1;
   double propagated0 = error0;
   double propagated1 = error1;
   double propagated2 = error2;

   if (kernel == 8) {
      h = z*h7 - mass;
      derivative = h7 + z*h7_derivative;
      second_derivative = 2.0*h7_derivative + z*h7_second_derivative;
      propagated0 = abs(z)*error0 + mass;
      propagated1 = error0 + abs(z)*error1;
      propagated2 = 2.0*error1 + abs(z)*error2;
   }

   long double numerator;
   long double propagated_error;
   if (power == 2) {
      numerator = static_cast<long double>(derivative.real())
         - static_cast<long double>(h.imag())/y;
      propagated_error = propagated1 + propagated0/y;
      result = numerator/(2.0L*PI_L*PI_L);
   } else {
      numerator = static_cast<long double>(second_derivative.imag())
         + 3.0L*static_cast<long double>(derivative.real())/y
         - 3.0L*static_cast<long double>(h.imag())/(y*y);
      propagated_error = propagated2 + 3.0L*propagated1/y
         + 3.0L*propagated0/(y*y);
      result = numerator/(8.0L*PI_L*PI_L*PI_L);
   }

   if (!isfinite(result) || result < 0.0L || numerator == 0.0L)
      return false;
   const long double relative_error = numeric_limits<double>::epsilon()*
      propagated_error/fabsl(numerator);
   return relative_error <= 1.0e-11L;
}

struct gaussian_quadrature_params {
   int kernel;
   int power;
   long double x;
   long double y;
   bool local_coordinate;
};

double gaussian_quadrature_integrand(double argument, void *raw_params)
{
   const gaussian_quadrature_params &params =
      *static_cast<gaussian_quadrature_params *>(raw_params);
   const long double x = params.x;
   const long double y = params.y;
   const long double epsilon = params.local_coordinate
      ? x + y*argument : argument;
   if (epsilon < 0.0L || epsilon > 10.0L)
      return 0.0;

   const long double qminus = params.local_coordinate
      ? y*y*(1.0L + static_cast<long double>(argument)*argument)
      : (x - epsilon)*(x - epsilon) + y*y;
   const long double qplus = (x + epsilon)*(x + epsilon) + y*y;
   const long double gaussian = expl(-2.0L*epsilon*epsilon);
   long double rational;
   if (params.kernel == 7) {
      rational = 1.0L/powl(qminus, params.power)
         + 1.0L/powl(qplus, params.power);
   } else if (params.power == 2) {
      rational = epsilon*4.0L*x*epsilon*(qplus + qminus)/
         (qminus*qminus*qplus*qplus);
   } else {
      rational = epsilon*4.0L*x*epsilon*
         (qplus*qplus + qplus*qminus + qminus*qminus)/
         (qminus*qminus*qminus*qplus*qplus*qplus);
   }

   long double value = gaussian*rational*powl(y/PI_L, params.power);
   if (params.local_coordinate)
      value *= y;
   return static_cast<double>(value);
}

bool integrate_gaussian_segment(gsl_function &function,
                                double lower, double upper,
                                gsl_integration_workspace *workspace,
                                long double &sum)
{
   if (!(lower < upper))
      return true;
   double value;
   double error;
   const int status = gsl_integration_qag(&function, lower, upper, 0.0, 2.0e-13,
      300, GSL_INTEG_GAUSS61, workspace, &value, &error);
   if (!isfinite(value) || !isfinite(error))
      return false;
   if (status != GSL_SUCCESS && status != GSL_EROUND) {
      if (status != GSL_ESING || sum == 0.0L ||
          fabsl(value) + error > 1.0e-12L*fabsl(sum))
         return false;
   }
   sum += value;
   return true;
}

bool gaussian_quadrature_kernel(int kernel, int power,
                                long double x, long double y,
                                long double &result)
{
   // This path is only used when cancellation diagnostics reject the exact
   // Faddeeva identities and the asymptotic series is not applicable.
   gsl_set_error_handler_off();
   gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(300);
   if (workspace == NULL)
      return false;

   const long double local_scale = y < 1.0e-7L ? 1.0e6L : 1024.0L;
   const long double local_lower = min(10.0L, max(0.0L, x - local_scale*y));
   const long double local_upper = max(0.0L, min(10.0L, x + local_scale*y));
   const double landmarks[] = {0.0, 0.25, 0.5, 1.0, 2.0, 3.0,
                               4.0, 5.0, 6.0, 8.0, 10.0};
   vector<double> direct_points(landmarks, landmarks + sizeof(landmarks)/sizeof(*landmarks));
   direct_points.push_back(static_cast<double>(local_lower));
   direct_points.push_back(static_cast<double>(local_upper));
   sort(direct_points.begin(), direct_points.end());
   direct_points.erase(unique(direct_points.begin(), direct_points.end()), direct_points.end());

   gaussian_quadrature_params params = {kernel, power, x, y, false};
   gsl_function function;
   function.function = &gaussian_quadrature_integrand;
   function.params = &params;
   long double sum = 0.0L;
   bool success = true;
   if (local_lower < local_upper) {
      const double local_landmarks[] = {-1.0e5, -1.0e4, -1.0e3, -256.0, -64.0,
         -16.0, -4.0, -1.0, 0.0, 1.0, 4.0, 16.0, 64.0, 256.0,
         1.0e3, 1.0e4, 1.0e5};
      vector<double> local_points;
      const double lower = static_cast<double>((local_lower - x)/y);
      const double upper = static_cast<double>((local_upper - x)/y);
      local_points.push_back(lower);
      for (size_t i = 0; i < sizeof(local_landmarks)/sizeof(*local_landmarks); ++i)
         if (lower < local_landmarks[i] && local_landmarks[i] < upper)
            local_points.push_back(local_landmarks[i]);
      local_points.push_back(upper);

      params.local_coordinate = true;
      for (size_t i = 1; success && i < local_points.size(); ++i)
         success = integrate_gaussian_segment(function, local_points[i - 1],
            local_points[i], workspace, sum);
   }

   params.local_coordinate = false;
   for (size_t i = 1; success && i < direct_points.size(); ++i) {
      const double lower = direct_points[i - 1];
      const double upper = direct_points[i];
      const long double midpoint = 0.5L*(lower + upper);
      if (local_lower < midpoint && midpoint < local_upper)
         continue;
      success = integrate_gaussian_segment(function, lower, upper, workspace, sum);
   }

   gsl_integration_workspace_free(workspace);
   result = sum;
   return success && isfinite(result) && result >= 0.0L;
}

long double gaussian_kernel(int kernel, int power, long double x, long double y)
{
   long double result;
   if ((x >= 5.5L || hypotl(x, y) >= 8.0L) &&
       gaussian_tail_kernel(kernel, power, x, y, result))
      return result;

   if (gaussian_hilbert_kernel(kernel, power,
                               static_cast<double>(x), static_cast<double>(y), result))
      return result;

   if (gaussian_quadrature_kernel(kernel, power, x, y, result))
      return result;

   cerr << "Gaussian clean-limit kernel integration failed for m=" << kernel
        << ", n=" << power << ", z=(" << static_cast<double>(x)
        << "," << static_cast<double>(y) << ")" << endl;
   exit(EXIT_FAILURE);
}

} // namespace

double J_builtin(int kernel, int power, complex<double> z)
{
   if (kernel < 1 || kernel > 8 || (power != 2 && power != 3) ||
       !isfinite(z.real()) || !isfinite(z.imag()) || z.imag() == 0.0) {
      cerr << "Invalid built-in kernel request: m=" << kernel
           << ", n=" << power << ", z=" << z << endl;
      exit(EXIT_FAILURE);
   }

   const long double x = fabsl(static_cast<long double>(z.real()));
   const long double y = fabsl(static_cast<long double>(z.imag()));
   if (kernel % 2 == 0 && x == 0.0L)
      return 0.0;

   long double value;
   if (kernel <= 6) {
      if (hypotl(x, y) >= 2.0L) {
         value = finite_moment_series(kernel, power, x, y);
      } else if (kernel <= 2 && x > 1.0L && y < 0.75L*(x - 1.0L)) {
         value = flat_exterior_kernel(kernel, power, x, y);
      } else if (kernel <= 2) {
         value = flat_kernel(kernel, power, x, y);
      } else {
         value = square_root_kernel(kernel, power, x, y);
      }
   } else {
      value = gaussian_kernel(kernel, power, x, y);
   }

   if (kernel % 2 == 0 && z.real() < 0.0)
      value = -value;
   if (power % 2 == 1 && z.imag() < 0.0)
      value = -value;
   return static_cast<double>(value);
}

bool same_analytic_region(int kernel, complex<double> a, complex<double> b)
{
   if ((a.imag() > 0.0 && b.imag() > 0.0) ||
       (a.imag() < 0.0 && b.imag() < 0.0))
      return true;

   if (kernel <= 6 && ((a.real() > 1.0 && b.real() > 1.0) ||
		       (a.real() < -1.0 && b.real() < -1.0)))
      return true;

   return false;
}

complex<double> hilbert_divided_difference(int kernel,
					   complex<double> a,
					   complex<double> b)
{
   if (a == b)
      return -kernel_hilbert(kernel, a).derivative;

   if (3 <= kernel && kernel <= 6) {
      // Factor H(a)-H(b) in t=1/(z+sqrt(z-1)*sqrt(z+1)).
      const complex<double> root_a = sqrt(a - 1.0)*sqrt(a + 1.0);
      const complex<double> root_b = sqrt(b - 1.0)*sqrt(b + 1.0);
      const complex<double> ua = a + root_a;
      const complex<double> ub = b + root_b;
      const complex<double> ta = 1.0/ua;
      const complex<double> tb = 1.0/ub;
      const complex<double> root_sum = root_a + root_b;
      const complex<double> factor_a = a - b + root_sum;
      const complex<double> factor_b = b - a + root_sum;
      const complex<double> denominator = abs(factor_a) >= abs(factor_b)
	 ? ub*factor_a
	 : ua*factor_b;
      const complex<double> ratio = 2.0/denominator;

      if (kernel == 3)
	 return M_PI*ratio;
      if (kernel == 4)
	 return 0.5*M_PI*ratio*(ta + tb);
      if (kernel == 5)
	 return 0.25*M_PI*ratio*(3.0 - ta*ta - ta*tb - tb*tb);
      return 0.125*M_PI*ratio*(ta + tb)*(2.0 - ta*ta - tb*tb);
   }

   const complex<double> difference = b - a;
   const complex<double> center = 0.5*(a + b);
   const double local_scale = kernel <= 6
      ? min(abs(center - 1.0), abs(center + 1.0))
      : max(1.0, abs(center));
   static const double threshold = pow(numeric_limits<double>::epsilon(), 1.0/7.0);

   if (same_analytic_region(kernel, a, b) &&
       abs(difference) <= threshold*local_scale) {
      const complex<double> offset = 0.5*sqrt(3.0/5.0)*difference;
      return -(5.0/18.0)*(kernel_hilbert(kernel, center - offset).derivative +
			  kernel_hilbert(kernel, center + offset).derivative)
	     -(4.0/9.0)*kernel_hilbert(kernel, center).derivative;
   }

   return (kernel_hilbert(kernel, a).value - kernel_hilbert(kernel, b).value)/difference;
}

double J_m2_optical(complex<double> z1, complex<double> z2)
{
   if (m == 0)
      return J_0_optical(z1, z2);
   if (z1 == z2)
      return J_builtin(m, 2, z1);

   const complex<double> cross = hilbert_divided_difference(m, conj(z1), z2);
   const complex<double> same = hilbert_divided_difference(m, z1, z2);
   return (cross.real() - same.real())/(2.0*M_PI*M_PI);
}

complex<double> effective_frequency(double omega)
{
   double sigma_re;
   double sigma_im;
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

   return omega + mu - sigma_re - sigma_im*I;
}

double optical_fermi_factor(double omega)
{
   // Algebraic branches avoid both exponential overflow and close subtraction.
   const double x = omega/T;
   const double shifted_x = (omega + optical_frequency)/T;
   const double shift = optical_frequency/T;
   const double fermi_difference_scale = -expm1(-shift);

   if (x >= 0.0) {
      const double exp_minus_x = exp(-x);
      const double exp_minus_shifted_x = exp(-shifted_x);
      return exp_minus_x*fermi_difference_scale/
	 (optical_frequency*(1.0 + exp_minus_x)*(1.0 + exp_minus_shifted_x));
   }

   if (shifted_x <= 0.0) {
      const double exp_x = exp(x);
      const double exp_shifted_x = exp(shifted_x);
      return exp_shifted_x*fermi_difference_scale/
	 (optical_frequency*(1.0 + exp_x)*(1.0 + exp_shifted_x));
   }

   return fermi_difference_scale/
      (optical_frequency*(1.0 + exp(x))*(1.0 + exp(-shifted_x)));
}

// Integrand of Imno integral
double integrand (double omega, void * params)
{
   double T = *(double *) params;
   const complex<double> OMEGA = effective_frequency(omega);

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

double optical_integrand(double omega, void *)
{
   const complex<double> z1 = effective_frequency(omega);
   const complex<double> z2 = effective_frequency(omega + optical_frequency);
   return optical_fermi_factor(omega)*J_m2_optical(z1, z2)*pow(omega, o);
}

// Calculate the epsilon integral only. For m=0 case only.
// For n=1 this is the density of states, which is the case of main
// interest.
void calc_DOS()
{
   assert(m == 0);
   if (n != 1 && !quiet_warnings) {
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
      const complex<double> OMEGA = effective_frequency(omega);
      
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
      if (val > -sigma_clip)
	 val = -sigma_clip;
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

void calc_optical()
{
   assert(optical_mode && optical_frequency > 0.0);

   const size_t ws_size = 1000;
   gsl_integration_workspace *work = gsl_integration_workspace_alloc(ws_size);
   if (work == NULL) {
      cerr << "Unable to allocate optical frequency integration workspace." << endl;
      exit(EXIT_FAILURE);
   }

   const double lower_limit = -cutoff*T - optical_frequency;
   const double upper_limit = cutoff*T;
   const double boundaries[] = {lower_limit, -optical_frequency, 0.0, upper_limit};

   if (!quiet_warnings &&
       (lower_limit < omega_min || upper_limit + optical_frequency > omega_max)) {
      cerr << "Warning: optical calculation evaluates Sigma outside ["
	   << omega_min << ", " << omega_max << "]." << endl;
   }

   gsl_function F;
   F.function = &optical_integrand;
   F.params = NULL;

   gsl_set_error_handler_off();
   double result = 0.0;
   double error = 0.0;
   bool roundoff_limited = false;
   for (int segment = 0; segment < 3; ++segment) {
      double segment_result;
      double segment_error;
      const int status = gsl_integration_qag(&F,
					     boundaries[segment],
					     boundaries[segment + 1],
					     abs_error/3.0,
					     rel_error,
					     ws_size,
					     key,
					     work,
					     &segment_result,
					     &segment_error);
      if ((status != GSL_SUCCESS && status != GSL_EROUND) || !isfinite(segment_result)) {
	 gsl_integration_workspace_free(work);
	 cerr << "Optical frequency integration failed on ["
	      << boundaries[segment] << ", " << boundaries[segment + 1]
	      << "]: " << gsl_strerror(status) << endl;
	 exit(EXIT_FAILURE);
      }
      if (status == GSL_EROUND)
	 roundoff_limited = true;
      result += segment_result;
      error += segment_error;
   }
   gsl_integration_workspace_free(work);

   if (roundoff_limited && !quiet_warnings)
      cerr << "Warning: optical frequency integration was limited by roundoff; estimated error="
	   << error << "." << endl;

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

#ifndef BUBBLE_NO_MAIN
int main (int argc, char *argv[])
{
   cmd_line(argc, argv);

   load_Sigma();
   if (m == 0)
      load_Phi();
   
   if (calcdos)
      calc_DOS();
   else if (optical_mode && optical_frequency > 0.0)
      calc_optical();
   else
      calc();
    
   return 0;
}
#endif
