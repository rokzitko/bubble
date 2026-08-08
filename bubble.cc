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
//            - support Fermi-level epsilon windows through -M
// 8.8.2026 - strict integration error handling with -E compatibility mode
//            - strict numeric parsing and negative positional arguments
//            - enforce the ImSigma floor after interpolation

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <set>
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
vector<double> sigma_omega_knots;
bool particle_hole_symmetric_sigma = false;

// Interpolation object for generic Phi(epsilon) function
gsl_interp_accel *acc_Phi;
gsl_spline *spline_Phi;
double eps_min, eps_max; // Interval boundaries

const string VERSION = "1.9";

// Mandatory parameters
int m, n, o;
double T;
double mu;
string fnReSigma, fnImSigma;

// Optional parameters (with defaults)
int b = 1;
bool verbose = false;
bool quiet_warnings = false;
bool ignore_integration_errors = false;
int key = GSL_INTEG_GAUSS15;
double abs_error = 1.0e-7;
double rel_error = 1.0e-8;
double cutoff = 15.0;
double sigma_clip = 1.0e-8;
double epsilon_window_limit = 0.0;
double epsilon_window_center = 0.0;
double epsilon_window_lower = 0.0;
double epsilon_window_upper = 0.0;
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
   cout << "Usage: bubble [options] <m> <n> <o> <T> <mu> <resigma> <imsigma>" << endl;
   cout << "Options must precede the positional arguments." << endl;
   cout << "Options:" << endl;
   cout << "-v : increase verbosity" << endl;
   cout << "-q : suppress non-fatal warnings" << endl;
   cout << "-E, --ignore-integration-errors : continue with finite partial integration results" << endl;
   cout << "-i I : interpolation (default = 1, 1=>linear, 2=>cspline, 3=>Akima spline)" << endl;
   cout << "-k K : integration rule (default = 1, 1 => 15, 2 => 21, etc.)" << endl;
   cout << "-a A : nonnegative absolute error (default = 1e-7)" << endl;
   cout << "-r R : nonnegative relative error (default = 1e-8)" << endl;
   cout << "-c C : positive frequency cutoff in units of T (default = 15)" << endl;
   cout << "-s FLOOR : positive in-table ImSigma clipping floor (default = 1e-8)" << endl;
   cout << "-M LIMIT : epsilon-window half-width around the Fermi level (default = 0, unrestricted)" << endl;
   cout << "-p : filename for Phi tables (default = Phi.dat)" << endl;
   cout << "-d : compute the epsilon integrals only, for m=0 (output = dos.dat)" << endl;
   cout << "-e E : nonnegative power of epsilon when using the m=0 code (default=0)" << endl;
   cout << "-f : switch (-df/dw) to f in the w integration (incompatible with -d)" << endl;
   cout << "-O OMEGA : external optical frequency (requires OMEGA>=0 and n=2)" << endl;
}

void invalid_numeric_value(const char *description, const char *value)
{
   cerr << "Invalid " << description << ": " << value << endl;
   exit(EXIT_FAILURE);
}

int parse_integer(const char *value, const char *description)
{
   char *end = NULL;
   errno = 0;
   const long result = strtol(value, &end, 10);

   if (errno == ERANGE || end == value || *end != '\0' ||
       result < numeric_limits<int>::min() || result > numeric_limits<int>::max())
      invalid_numeric_value(description, value);

   return static_cast<int>(result);
}

double parse_finite_double(const char *value, const char *description)
{
   char *end = NULL;
   errno = 0;
   const double result = strtod(value, &end);

   if (errno == ERANGE || end == value || *end != '\0' || !isfinite(result))
      invalid_numeric_value(description, value);

   return result;
}

double parse_optical_frequency(const char *value)
{
   const double result = parse_finite_double(value, "optical frequency");

   if (result < 0.0)
      invalid_numeric_value("optical frequency", value);

   return result;
}

double parse_sigma_clip(const char *value)
{
   const double result = parse_finite_double(value, "ImSigma clipping floor");

   if (result <= 0.0)
      invalid_numeric_value("ImSigma clipping floor", value);

   return result;
}

double parse_epsilon_window_limit(const char *value)
{
   const double result = parse_finite_double(value, "epsilon-window limit");

   if (result < 0.0)
      invalid_numeric_value("epsilon-window limit", value);

   return result;
}

void cmd_line(int argc, char *argv[])
{
   const int positional_count = 7;
   if (argc < positional_count + 1) {
      usage();
      exit(EXIT_FAILURE);
   }

   const int first_positional = argc - positional_count;
   vector<char *> option_argv(argv, argv + first_positional);
   option_argv.push_back(NULL);

   int c;
   static const struct option long_options[] = {
      {"ignore-integration-errors", no_argument, NULL, 'E'},
      {NULL, 0, NULL, 0}
   };

   optind = 1;
   while ((c = getopt_long(first_positional, &option_argv[0],
                           "+vqEi:k:a:r:c:s:M:p:dfe:O:",
                           long_options, NULL)) != -1) {
      switch (c) {
      case 'v':
	 verbose = true;
	 break;
      case 'q':
	 quiet_warnings = true;
	 break;
      case 'E':
	 ignore_integration_errors = true;
	 break;
      case 'i':
	 b = parse_integer(optarg, "interpolation selector");
	 break;
      case 'k':
         key = parse_integer(optarg, "integration rule");
         break;
      case 'a':
	 abs_error = parse_finite_double(optarg, "absolute integration tolerance");
	 break;
      case 'r':
	 rel_error = parse_finite_double(optarg, "relative integration tolerance");
	 break;
      case 'c':
	 cutoff = parse_finite_double(optarg, "frequency cutoff");
	 break;
      case 's':
	 sigma_clip = parse_sigma_clip(optarg);
	 break;
      case 'M':
	 epsilon_window_limit = parse_epsilon_window_limit(optarg);
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
	 e = parse_integer(optarg, "epsilon power");
	 break;
      case 'O':
	 optical_mode = true;
	 optical_frequency = parse_optical_frequency(optarg);
	 break;
      default:
	 cerr << "Invalid or incomplete command-line option." << endl;
	 usage();
	 exit(EXIT_FAILURE);
      }
   }

   if (optind != first_positional) {
      cerr << "Options must precede the seven positional arguments." << endl;
      usage();
      exit(EXIT_FAILURE);
   }

   m = parse_integer(argv[first_positional], "kernel index m");
   n = parse_integer(argv[first_positional + 1], "spectral power n");
   o = parse_integer(argv[first_positional + 2], "frequency power o");
   T = parse_finite_double(argv[first_positional + 3], "temperature T");
   mu = parse_finite_double(argv[first_positional + 4], "chemical potential mu");
   fnReSigma = string(argv[first_positional + 5]);
   fnImSigma = string(argv[first_positional + 6]);

   if (m < 0 || m > 8) {
      cerr << "Unsupported kernel index m=" << m << "." << endl;
      exit(EXIT_FAILURE);
   }
   if (n < 0 || (m != 0 && n > 3)) {
      cerr << "Unsupported spectral power n=" << n << " for m=" << m << "." << endl;
      exit(EXIT_FAILURE);
   }
   if (b < 1 || b > 3) {
      cerr << "Unsupported interpolation selector i=" << b << "." << endl;
      exit(EXIT_FAILURE);
   }
   if (key < GSL_INTEG_GAUSS15 || key > GSL_INTEG_GAUSS61) {
      cerr << "Unsupported integration rule k=" << key << "." << endl;
      exit(EXIT_FAILURE);
   }
   if (abs_error < 0.0) {
      cerr << "Absolute integration tolerance must be nonnegative." << endl;
      exit(EXIT_FAILURE);
   }
   if (rel_error < 0.0) {
      cerr << "Relative integration tolerance must be nonnegative." << endl;
      exit(EXIT_FAILURE);
   }
   if (m == 0 && e < 0) {
      cerr << "Tabulated epsilon power e must be nonnegative." << endl;
      exit(EXIT_FAILURE);
   }

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
   if (!calcdos) {
      if (o < 0) {
         cerr << "Frequency power o must be nonnegative during frequency integration." << endl;
         exit(EXIT_FAILURE);
      }
      if (T <= 0.0) {
         cerr << "Frequency integration requires positive finite T." << endl;
         exit(EXIT_FAILURE);
      }
      if (cutoff <= 0.0) {
         cerr << "Frequency integration requires a positive finite cutoff." << endl;
         exit(EXIT_FAILURE);
      }
      const double thermal_bound = cutoff*T;
      if (!isfinite(thermal_bound) || thermal_bound <= 0.0 ||
          (optical_mode && !isfinite(thermal_bound + optical_frequency))) {
         cerr << "The frequency integration bounds are not representable." << endl;
         exit(EXIT_FAILURE);
      }
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
      if (ignore_integration_errors)
	 cout << "finite partial integration results are accepted" << endl;
      cout << "ImSigma clipping floor=" << sigma_clip << endl;
      if (epsilon_window_limit > 0.0)
	 cout << "epsilon-window half-width M=" << epsilon_window_limit << endl;
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

bool integration_result_acceptable(const string &operation,
                                   int status, double result, double error)
{
   const bool finite_result = isfinite(result);
   const bool valid_error = isfinite(error) && error >= 0.0;
   if (status == GSL_SUCCESS && finite_result && valid_error)
      return true;

   ostringstream message;
   message << setprecision(17) << operation;
   if (status != GSL_SUCCESS)
      message << ": " << gsl_strerror(status) << " (status=" << status << ")";
   else if (!finite_result)
      message << ": non-finite result";
   else
      message << ": invalid error estimate";
   message << ", result=" << result << ", error=" << error;

   if (!finite_result) {
      cerr << "Error: " << message.str() << endl;
      return false;
   }

   if (ignore_integration_errors) {
      if (!quiet_warnings) {
         static set<string> warned_operations;
         if (warned_operations.insert(operation).second)
            cerr << "Warning: ignoring integration error: " << message.str() << endl;
      }
      return true;
   }

   cerr << "Error: " << message.str() << endl;
   return false;
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

   gsl_set_error_handler_off();
   gsl_integration_workspace *w = gsl_integration_workspace_alloc (ws_size);
   if (w == NULL) {
      cerr << "Unable to allocate custom-kernel integration workspace." << endl;
      exit(EXIT_FAILURE);
   }

   double result_Phi = numeric_limits<double>::quiet_NaN();
   double error_Phi = numeric_limits<double>::quiet_NaN();

   const int status = gsl_integration_qag(&F,
				          eps_min,
				          eps_max,
				          abs_error,
				          rel_error,
				          ws_size,
				          key,
				          w,
				          &result_Phi,
				          &error_Phi);
   gsl_integration_workspace_free(w);

   if (!integration_result_acceptable("Custom epsilon integration",
                                      status, result_Phi, error_Phi))
      exit(EXIT_FAILURE);

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

   double result = numeric_limits<double>::quiet_NaN();
   double error = numeric_limits<double>::quiet_NaN();
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

   if (!integration_result_acceptable("Optical epsilon integration",
                                      status, result, error)) {
      cerr << "Optical arguments were z1=" << z1 << ", z2=" << z2 << "." << endl;
      exit(EXIT_FAILURE);
   }

   return result;
}

struct epsilon_interval {
   double lower;
   double upper;
   bool truncated;
   bool empty;
};

epsilon_interval restricted_epsilon_interval(int kernel);
double J_restricted_single(int kernel, int power, complex<double> z,
                           double lower, double upper);
double J_restricted_optical(int kernel, complex<double> z1, complex<double> z2,
                            double lower, double upper);

// Switches between J for different m and n indices.
double J_mn (complex<double> OMEGA)
{
   if (epsilon_window_limit > 0.0) {
      const epsilon_interval interval = restricted_epsilon_interval(m);
      if (interval.empty)
	 return 0.0;
      if (interval.truncated) {
	 if (n == 0) {
	    static bool cached = false;
	    static double cached_value = 0.0;
	    static int cached_kernel = -1;
	    static double cached_lower = 0.0;
	    static double cached_upper = 0.0;
	    if (!cached || cached_kernel != m || cached_lower != interval.lower ||
		cached_upper != interval.upper) {
	       cached_value = J_restricted_single(m, n, OMEGA,
		  interval.lower, interval.upper);
	       cached_kernel = m;
	       cached_lower = interval.lower;
	       cached_upper = interval.upper;
	       cached = true;
	    }
	    return cached_value;
	 }
	 return J_restricted_single(m, n, OMEGA, interval.lower, interval.upper);
      }
   }

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
   double value = numeric_limits<double>::quiet_NaN();
   double error = numeric_limits<double>::quiet_NaN();
   const int status = gsl_integration_qag(&function, lower, upper, 0.0, 2.0e-13,
      300, GSL_INTEG_GAUSS61, workspace, &value, &error);
   const bool negligible_failure = (status == GSL_ESING || status == GSL_EROUND) &&
      sum != 0.0L &&
      isfinite(value) && isfinite(error) && error >= 0.0 &&
      fabsl(value) + error <= 1.0e-12L*fabsl(sum);
   if ((status != GSL_SUCCESS || !isfinite(value) || !isfinite(error) || error < 0.0) &&
       !negligible_failure) {
      if (!integration_result_acceptable("Gaussian epsilon integration",
                                         status, value, error))
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

epsilon_interval restricted_epsilon_interval(int kernel)
{
   if (epsilon_window_limit == 0.0)
      return epsilon_interval{0.0, 0.0, false, false};

   if (kernel == 7 || kernel == 8)
      return epsilon_interval{epsilon_window_lower, epsilon_window_upper, true,
                              !(epsilon_window_lower < epsilon_window_upper)};

   double natural_lower;
   double natural_upper;
   if (kernel == 0) {
      natural_lower = eps_min;
      natural_upper = eps_max;
   } else if (1 <= kernel && kernel <= 6) {
      natural_lower = -1.0;
      natural_upper = 1.0;
   } else {
      cerr << "Epsilon window not implemented for m=" << kernel << endl;
      exit(EXIT_FAILURE);
   }

   if (epsilon_window_lower <= natural_lower && natural_upper <= epsilon_window_upper)
      return epsilon_interval{natural_lower, natural_upper, false, false};

   const double lower = max(natural_lower, epsilon_window_lower);
   const double upper = min(natural_upper, epsilon_window_upper);
   return epsilon_interval{lower, upper, true, !(lower < upper)};
}

namespace {

const size_t RESTRICTED_WORKSPACE_SIZE = 1000;

gsl_integration_workspace *restricted_workspace()
{
   gsl_set_error_handler_off();
   static gsl_integration_workspace *workspace =
      gsl_integration_workspace_alloc(RESTRICTED_WORKSPACE_SIZE);
   if (workspace == NULL) {
      cerr << "Unable to allocate restricted epsilon integration workspace." << endl;
      exit(EXIT_FAILURE);
   }
   return workspace;
}

double restricted_relative_error()
{
   return max(5.0e-13, min(1.0e-12, rel_error));
}

double restricted_acceptable_relative_error()
{
   return max(rel_error, 10.0*restricted_relative_error());
}

double restricted_acceptable_absolute_error()
{
   return max(1.0e-300, min(1.0e-12, 1.0e-4*abs_error));
}

long double restricted_error_bound(double result, long double result_scale = 1.0L)
{
   return max(static_cast<long double>(restricted_acceptable_absolute_error())/
                 fabsl(result_scale),
              static_cast<long double>(restricted_acceptable_relative_error())*
                 abs(result));
}

long double restricted_integrand_scale(long double full_scale)
{
   const long double magnitude = fabsl(full_scale);
   return min(1.0e200L, max(1.0e-200L, magnitude));
}

void warn_restricted_roundoff(const char *description)
{
   if (quiet_warnings)
      return;
   static bool warning_emitted = false;
   if (!warning_emitted) {
      cerr << "Warning: " << description
           << " was limited by roundoff within the accepted error bound." << endl;
      warning_emitted = true;
   }
}

long double restricted_phi_value(int kernel, long double epsilon)
{
   if (kernel == 0) {
      double point = static_cast<double>(epsilon);
      point = max(eps_min, min(eps_max, point));
      return gsl_spline_eval(spline_Phi, point, acc_Phi);
   }

   if (kernel <= 6 && fabsl(epsilon) > 1.0L)
      return 0.0L;

   const long double one_minus_square = max(0.0L, 1.0L - epsilon*epsilon);
   switch (kernel) {
   case 1:
      return 1.0L;
   case 2:
      return epsilon;
   case 3:
      return sqrtl(one_minus_square);
   case 4:
      return epsilon*sqrtl(one_minus_square);
   case 5:
      return one_minus_square*sqrtl(one_minus_square);
   case 6:
      return epsilon*one_minus_square*sqrtl(one_minus_square);
   case 7:
      return expl(-2.0L*epsilon*epsilon);
   case 8:
      return epsilon*expl(-2.0L*epsilon*epsilon);
   default:
      cerr << "Restricted kernel not implemented for m=" << kernel << endl;
      exit(EXIT_FAILURE);
   }
}

long double restricted_phi_value_at_offset(int kernel, long double center,
                                           long double offset)
{
   const long double epsilon = center + offset;
   if (kernel < 3 || kernel > 6)
      return restricted_phi_value(kernel, epsilon);

   // Keep distances from a finite-band edge before adding a sub-ULP offset
   // to a center near +/-1.
   const long double one_minus_square = max(0.0L,
      (1.0L - center - offset)*(1.0L + center + offset));
   const long double root = sqrtl(one_minus_square);
   switch (kernel) {
   case 3:
      return root;
   case 4:
      return epsilon*root;
   case 5:
      return one_minus_square*root;
   case 6:
      return epsilon*one_minus_square*root;
   default:
      return 0.0L;
   }
}

long double restricted_angle_span(long double lower_offset,
                                  long double upper_offset,
                                  long double width)
{
   const long double scale = max(width, max(fabsl(lower_offset), fabsl(upper_offset)));
   const long double scaled_width = width/scale;
   const long double scaled_lower = lower_offset/scale;
   const long double scaled_upper = upper_offset/scale;
   return atan2l(scaled_width*(scaled_upper - scaled_lower),
                 scaled_width*scaled_width + scaled_lower*scaled_upper);
}

struct restricted_single_params {
   int kernel;
   int power;
   long double center;
   long double width;
   long double signed_imaginary;
   bool transformed;
   bool logarithmic;
   bool scaled;
   long double transformed_origin;
   long double transformed_near;
   int transformed_direction;
   long double integrand_scale;
};

double restricted_single_integrand(double argument, void *raw_params)
{
   const restricted_single_params &params =
      *static_cast<restricted_single_params *>(raw_params);
   if (params.scaled) {
      const long double distance = static_cast<long double>(argument)*argument;
      const long double offset = params.transformed_direction*
         params.width*distance;
      return static_cast<double>(restricted_phi_value_at_offset(params.kernel,
         params.center, offset)*2.0L*argument/
         powl(1.0L + distance*distance, params.power)*params.integrand_scale);
   }
   if (params.logarithmic) {
      const long double distance = params.transformed_near*
         expl(static_cast<long double>(argument));
      const long double inverse = 1.0L/distance;
      const long double weight = distance <= 1.0L
         ? distance/powl(1.0L + distance*distance, params.power)
         : powl(inverse, 2*params.power - 1)/
           powl(1.0L + inverse*inverse, params.power);
      const long double offset = params.transformed_direction*
         params.width*distance;
      return static_cast<double>(restricted_phi_value_at_offset(params.kernel,
         params.center, offset)*weight*params.integrand_scale);
   }
   if (params.transformed) {
      const long double theta = params.transformed_origin + argument;
      const long double cosine = cosl(theta);
      const long double offset = params.width*tanl(theta);
      return static_cast<double>(restricted_phi_value_at_offset(params.kernel,
         params.center, offset)*
         powl(cosine, 2*params.power - 2)*params.integrand_scale);
   }

   const long double epsilon = argument;
   long double value = restricted_phi_value(params.kernel, epsilon);
   if (params.power != 0) {
      const long double difference = params.center - epsilon;
      const long double spectral = params.signed_imaginary/
         (PI_L*(difference*difference + params.width*params.width));
      value *= powl(spectral, params.power);
   }
   return static_cast<double>(value);
}

struct restricted_odd_pair_params {
   int kernel;
   int power;
   long double center;
   long double width;
   long double angle_origin;
   long double integrand_scale;
   bool transformed;
};

double restricted_odd_pair_integrand(double argument, void *raw_params)
{
   const restricted_odd_pair_params &params =
      *static_cast<restricted_odd_pair_params *>(raw_params);
   const long double theta = params.angle_origin + argument;
   const long double cosine = params.transformed ? cosl(theta) : 1.0L;
   const long double epsilon = params.transformed
      ? params.center + params.width*tanl(theta)
      : static_cast<long double>(argument);
   const long double scaled_center = params.center/params.width;
   const long double scaled_epsilon = epsilon/params.width;
   const long double minus_denominator =
      1.0L + (scaled_center - scaled_epsilon)*(scaled_center - scaled_epsilon);
   const long double plus_denominator =
      1.0L + (scaled_center + scaled_epsilon)*(scaled_center + scaled_epsilon);
   const long double minus_value = 1.0L/minus_denominator;
   const long double plus_value = 1.0L/plus_denominator;
   const long double difference = 4.0L*scaled_center*scaled_epsilon/
      (minus_denominator*plus_denominator);
   long double power_sum = 1.0L;
   if (params.power == 2)
      power_sum = minus_value + plus_value;
   else if (params.power == 3)
      power_sum = minus_value*minus_value + minus_value*plus_value +
                  plus_value*plus_value;
   const long double paired_spectral = difference*power_sum;
   return static_cast<double>(restricted_phi_value(params.kernel, epsilon)*
      paired_spectral/(cosine*cosine)*params.integrand_scale);
}

struct restricted_fixed_result {
   double value;
   double difference;
};

restricted_fixed_result integrate_restricted_fixed_check(gsl_function &function,
                                                          double lower, double upper)
{
   static gsl_integration_glfixed_table *coarse_table =
      gsl_integration_glfixed_table_alloc(64);
   static gsl_integration_glfixed_table *fine_table =
      gsl_integration_glfixed_table_alloc(128);
   if (coarse_table == NULL || fine_table == NULL) {
      cerr << "Unable to allocate restricted integration check table." << endl;
      exit(EXIT_FAILURE);
   }
   const double coarse = gsl_integration_glfixed(&function, lower, upper, coarse_table);
   const double fine = gsl_integration_glfixed(&function, lower, upper, fine_table);
   return restricted_fixed_result{fine, abs(fine - coarse)};
}

double integrate_restricted_qag(gsl_function &function, double lower, double upper,
                                const char *description,
                                long double result_scale = 1.0L)
{
   double result = numeric_limits<double>::quiet_NaN();
   double error = numeric_limits<double>::quiet_NaN();
   const int status = gsl_integration_qag(&function, lower, upper, 0.0,
      restricted_relative_error(), RESTRICTED_WORKSPACE_SIZE, GSL_INTEG_GAUSS61,
      restricted_workspace(), &result, &error);
   double verified_error = error;
   if (status == GSL_EROUND) {
      // GSL can miss the tighter internal target after already satisfying the
      // requested tolerance; otherwise require convergence of an independent rule.
      const restricted_fixed_result check =
         integrate_restricted_fixed_check(function, lower, upper);
      const double independent_error = max(abs(result - check.value), check.difference);
      if (isfinite(independent_error) && independent_error >= 0.0)
         verified_error = isfinite(error) && error >= 0.0
            ? min(error, independent_error) : independent_error;
      else
         verified_error = numeric_limits<double>::infinity();
   }
   const bool roundoff_acceptable = status == GSL_EROUND &&
       isfinite(result) && isfinite(verified_error) &&
       verified_error <= restricted_error_bound(result, result_scale);
   const bool successful = status == GSL_SUCCESS && isfinite(result) &&
      isfinite(error) && error >= 0.0;
   if (!successful && !roundoff_acceptable &&
       !integration_result_acceptable(description, status, result, error))
      exit(EXIT_FAILURE);
   if (roundoff_acceptable)
      warn_restricted_roundoff(description);
   return result;
}

struct restricted_optical_params {
   int kernel;
   complex<double> z1;
   complex<double> z2;
};

double restricted_optical_direct_integrand(double epsilon, void *raw_params)
{
   const restricted_optical_params &params =
      *static_cast<restricted_optical_params *>(raw_params);
   const complex<double> g1 = 1.0/(params.z1 - epsilon);
   const complex<double> g2 = 1.0/(params.z2 - epsilon);
   return static_cast<double>(restricted_phi_value(params.kernel, epsilon))
      *(-g1.imag()/M_PI)*(-g2.imag()/M_PI);
}

void add_restricted_peak_points(vector<double> &points, complex<double> z,
                                double lower, double upper)
{
   const double center = z.real();
   const double width = abs(z.imag());
   const double range = upper - lower;
   if (!isfinite(center) || !isfinite(width) || width == 0.0)
      return;
   const double interval_distance = center < lower ? lower - center
      : (center > upper ? center - upper : 0.0);
   if (interval_distance > 8192.0*width)
      return;

   if (lower < center && center < upper)
      points.push_back(center);
   double distance = width;
   for (int scale = 0; scale < 24; ++scale) {
      if (lower < center - distance && center - distance < upper)
         points.push_back(center - distance);
      if (lower < center + distance && center + distance < upper)
         points.push_back(center + distance);
      if (distance >= range || distance >= 4096.0*width)
         break;
      distance *= 8.0;
   }
}

vector<double> filter_integration_points(vector<double> points,
                                         double lower, double upper)
{
   sort(points.begin(), points.end());
   vector<double> filtered;
   filtered.push_back(lower);
   for (vector<double>::const_iterator point = points.begin(); point != points.end(); ++point) {
      if (lower < *point && *point < upper && filtered.back() < *point)
         filtered.push_back(*point);
   }
   filtered.push_back(upper);
   return filtered;
}

double integrate_restricted_qagp(gsl_function &function, vector<double> points,
                                 const char *description,
                                 long double result_scale = 1.0L)
{
   double result = numeric_limits<double>::quiet_NaN();
   double error = numeric_limits<double>::quiet_NaN();
   const int status = gsl_integration_qagp(&function, &points[0], points.size(),
      0.0, restricted_relative_error(), RESTRICTED_WORKSPACE_SIZE,
      restricted_workspace(), &result, &error);
   double verified_error = error;
   if (status == GSL_EROUND) {
      double fine = 0.0;
      double check_difference = 0.0;
      for (size_t interval = 1; interval < points.size(); ++interval) {
         const restricted_fixed_result check = integrate_restricted_fixed_check(function,
            points[interval - 1], points[interval]);
         fine += check.value;
         check_difference += check.difference;
      }
      const double independent_error = max(abs(result - fine), check_difference);
      if (isfinite(independent_error) && independent_error >= 0.0)
         verified_error = isfinite(error) && error >= 0.0
            ? min(error, independent_error) : independent_error;
      else
         verified_error = numeric_limits<double>::infinity();
   }
   const bool roundoff_acceptable = status == GSL_EROUND &&
       isfinite(result) && isfinite(verified_error) &&
       verified_error <= restricted_error_bound(result, result_scale);
   const bool successful = status == GSL_SUCCESS && isfinite(result) &&
      isfinite(error) && error >= 0.0;
   if (!successful && !roundoff_acceptable &&
       !integration_result_acceptable(description, status, result, error))
      exit(EXIT_FAILURE);
   if (roundoff_acceptable)
      warn_restricted_roundoff(description);
   return result;
}

double integrate_restricted_piecewise_qag(gsl_function &function,
                                           const vector<double> &points,
                                           const char *description,
                                           long double result_scale = 1.0L)
{
   double total = 0.0;
   double compensation = 0.0;
   double questionable_error = 0.0;
   int questionable_status = GSL_SUCCESS;
   for (size_t interval = 1; interval < points.size(); ++interval) {
      double result = numeric_limits<double>::quiet_NaN();
      double error = numeric_limits<double>::quiet_NaN();
      const int status = gsl_integration_qag(&function, points[interval - 1],
         points[interval], 0.0, restricted_relative_error(),
         RESTRICTED_WORKSPACE_SIZE, GSL_INTEG_GAUSS61,
         restricted_workspace(), &result, &error);
      const bool valid_output = isfinite(result) && isfinite(error) && error >= 0.0;
      if (status == GSL_EROUND && valid_output) {
         const restricted_fixed_result check = integrate_restricted_fixed_check(function,
             points[interval - 1], points[interval]);
         if (!isfinite(check.value) || !isfinite(check.difference)) {
            if (!integration_result_acceptable(description, status, result, error))
               exit(EXIT_FAILURE);
         } else {
            questionable_error += min(error,
               max(abs(result - check.value), check.difference));
            questionable_status = status;
         }
      } else if (status != GSL_SUCCESS || !valid_output) {
         if (!integration_result_acceptable(description, status, result, error))
            exit(EXIT_FAILURE);
      }
      const double adjusted = result - compensation;
      const double updated = total + adjusted;
      compensation = (updated - total) - adjusted;
      total = updated;
   }
   if (questionable_status != GSL_SUCCESS &&
         questionable_error > restricted_error_bound(total, result_scale)) {
      if (!integration_result_acceptable(description, questionable_status,
                                         total, questionable_error))
         exit(EXIT_FAILURE);
      questionable_status = GSL_SUCCESS;
   }
   if (questionable_status != GSL_SUCCESS)
      warn_restricted_roundoff(description);
   if (!isfinite(total)) {
      if (!integration_result_acceptable(description, GSL_SUCCESS, total, 0.0))
         exit(EXIT_FAILURE);
   }
   return total;
}

double integrate_restricted_single_direct(gsl_function &function, int kernel,
                                          double lower, double upper,
                                          const char *description)
{
   if (kernel != 7 && kernel != 8)
      return integrate_restricted_qag(function, lower, upper, description);

   static const double gaussian_landmarks[] = {
      -10.0, -8.0, -5.0, -3.0, -2.0, -1.0, 0.0,
      1.0, 2.0, 3.0, 5.0, 8.0, 10.0
   };
   vector<double> points;
   points.push_back(lower);
   points.push_back(upper);
   for (size_t i = 0; i < sizeof(gaussian_landmarks)/sizeof(*gaussian_landmarks); ++i)
      if (lower < gaussian_landmarks[i] && gaussian_landmarks[i] < upper)
         points.push_back(gaussian_landmarks[i]);
   return integrate_restricted_qagp(function,
      filter_integration_points(points, lower, upper), description);
}

struct restricted_optical_transformed_params {
   int kernel;
   long double anchor_center;
   long double anchor_width;
   long double angle_origin;
   complex<double> other;
   long double integrand_scale;
};

double restricted_optical_transformed_integrand(double argument, void *raw_params)
{
   const restricted_optical_transformed_params &params =
      *static_cast<restricted_optical_transformed_params *>(raw_params);
   const long double tangent = tanl(params.angle_origin + argument);
   const long double offset = params.anchor_width*tangent;
   const long double other_width = abs(params.other.imag());
   const long double normalized_difference =
      (static_cast<long double>(params.other.real()) - params.anchor_center)/other_width
      - params.anchor_width/other_width*tangent;
   const long double normalized_other =
      1.0L/(normalized_difference*normalized_difference + 1.0L);
   return static_cast<double>(restricted_phi_value_at_offset(params.kernel,
      params.anchor_center, offset)*normalized_other*params.integrand_scale);
}

struct restricted_optical_log_params {
   int kernel;
   long double anchor_center;
   long double anchor_width;
   long double near_distance;
   int direction;
   complex<double> other;
   long double integrand_scale;
};

double restricted_optical_log_integrand(double argument, void *raw_params)
{
   const restricted_optical_log_params &params =
      *static_cast<restricted_optical_log_params *>(raw_params);
   const long double distance = params.near_distance*
      expl(static_cast<long double>(argument));
   const long double offset = params.direction*params.anchor_width*distance;
   const long double other_width = abs(params.other.imag());
   const long double normalized_difference =
      (static_cast<long double>(params.other.real()) - params.anchor_center)/other_width
      - params.direction*params.anchor_width/other_width*distance;
   const long double normalized_other =
      1.0L/(normalized_difference*normalized_difference + 1.0L);
   const long double weight = distance <= 1.0L
      ? distance/(1.0L + distance*distance)
      : 1.0L/(distance + 1.0L/distance);
   return static_cast<double>(restricted_phi_value_at_offset(params.kernel,
      params.anchor_center, offset)*normalized_other*weight*params.integrand_scale);
}

double integrate_restricted_optical_direct(int kernel, complex<double> z1,
                                           complex<double> z2,
                                           double lower, double upper)
{
   restricted_optical_params params = {kernel, z1, z2};
   gsl_function function;
   function.function = &restricted_optical_direct_integrand;
   function.params = &params;
   vector<double> points;
   points.push_back(lower);
   points.push_back(upper);
   add_restricted_peak_points(points, z1, lower, upper);
   add_restricted_peak_points(points, z2, lower, upper);
   return integrate_restricted_qagp(function,
      filter_integration_points(points, lower, upper),
      "Restricted optical epsilon integration");
}

double integrate_restricted_optical_exterior_log(int kernel,
                                                  complex<double> anchor,
                                                  complex<double> other,
                                                  long double lower,
                                                  long double upper)
{
   const long double center = anchor.real();
   const long double width = abs(anchor.imag());
   const int direction = center < lower ? 1 : -1;
   const long double near_distance = direction > 0
      ? (lower - center)/width : (center - upper)/width;
   const long double far_distance = direction > 0
      ? (upper - center)/width : (center - lower)/width;
   const double log_lower = 0.0;
   const double log_upper = static_cast<double>(logl(far_distance/near_distance));
   if (!(log_lower < log_upper)) {
      const double direct_lower = static_cast<double>(lower);
      const double direct_upper = static_cast<double>(upper);
      if (direct_lower < direct_upper)
         return integrate_restricted_optical_direct(kernel, anchor, other,
                                                     direct_lower, direct_upper);
      return 0.0;
   }

   const long double full_scale =
      (anchor.imag() < 0.0 ? -1.0L : 1.0L)*
      (other.imag() < 0.0 ? -1.0L : 1.0L)/
      (PI_L*PI_L*static_cast<long double>(abs(other.imag())));
   const long double integrand_scale = restricted_integrand_scale(full_scale);
   const long double result_scale = full_scale/integrand_scale;
   restricted_optical_log_params params = {
      kernel, center, width, near_distance, direction, other, integrand_scale
   };
   gsl_function function;
   function.function = &restricted_optical_log_integrand;
   function.params = &params;

   vector<double> energy_points;
   const double energy_lower = static_cast<double>(lower);
   const double energy_upper = static_cast<double>(upper);
   if (energy_lower < energy_upper) {
      energy_points.push_back(energy_lower);
      energy_points.push_back(energy_upper);
      add_restricted_peak_points(energy_points, other, energy_lower, energy_upper);
   }
   vector<double> log_points;
   log_points.push_back(log_lower);
   log_points.push_back(log_upper);
   for (vector<double>::const_iterator point = energy_points.begin();
        point != energy_points.end(); ++point) {
      if (lower < *point && *point < upper) {
         const long double distance = direction > 0
             ? (*point - center)/width : (center - *point)/width;
         log_points.push_back(static_cast<double>(logl(distance/near_distance)));
      }
   }
   log_points = filter_integration_points(log_points, log_lower, log_upper);
   const long double normalized = integrate_restricted_piecewise_qag(function,
      log_points, "Log-transformed restricted optical epsilon integration", result_scale);
   return static_cast<double>(result_scale*normalized);
}

double integrate_restricted_optical_transformed(int kernel,
                                                 complex<double> anchor,
                                                 complex<double> other,
                                                 long double lower,
                                                 long double upper)
{
   const long double width = abs(anchor.imag());
   const long double center = anchor.real();
   if (center < lower || center > upper)
      return integrate_restricted_optical_exterior_log(kernel, anchor, other,
                                                         lower, upper);
   const long double lower_offset = lower - center;
   const long double upper_offset = upper - center;
   const long double theta_lower = atan2l(lower_offset, width);
   const double transformed_lower = 0.0;
   const double transformed_upper = static_cast<double>(
      restricted_angle_span(lower_offset, upper_offset, width));
   if (!(transformed_lower < transformed_upper)) {
      const double direct_lower = static_cast<double>(lower);
      const double direct_upper = static_cast<double>(upper);
      if (direct_lower < direct_upper)
         return integrate_restricted_optical_direct(kernel, anchor, other,
                                                     direct_lower, direct_upper);
      return 0.0;
   }

   const long double full_scale =
      (anchor.imag() < 0.0 ? -1.0L : 1.0L)*
      (other.imag() < 0.0 ? -1.0L : 1.0L)/
      (PI_L*PI_L*static_cast<long double>(abs(other.imag())));
   const long double integrand_scale = restricted_integrand_scale(full_scale);
   const long double result_scale = full_scale/integrand_scale;
   restricted_optical_transformed_params params = {
      kernel, center, width, theta_lower, other, integrand_scale
   };
   gsl_function function;
   function.function = &restricted_optical_transformed_integrand;
   function.params = &params;

   vector<double> energy_points;
   const double energy_lower = static_cast<double>(lower);
   const double energy_upper = static_cast<double>(upper);
   if (energy_lower < energy_upper) {
      energy_points.push_back(energy_lower);
      energy_points.push_back(energy_upper);
   }
   const double other_distance = other.real() < lower ? lower - other.real()
      : (other.real() > upper ? other.real() - upper : 0.0);
   if (energy_lower < energy_upper && other_distance <= 16.0*abs(other.imag()))
      add_restricted_peak_points(energy_points, other, energy_lower, energy_upper);
   vector<double> theta_points;
   theta_points.push_back(transformed_lower);
   theta_points.push_back(transformed_upper);
   for (vector<double>::const_iterator point = energy_points.begin();
         point != energy_points.end(); ++point) {
      if (lower < *point && *point < upper)
         theta_points.push_back(static_cast<double>(restricted_angle_span(
            lower_offset, *point - center, width)));
   }
   theta_points = filter_integration_points(theta_points,
      transformed_lower, transformed_upper);
   const long double normalized = integrate_restricted_piecewise_qag(function,
      theta_points, "Transformed restricted optical epsilon integration", result_scale);
   return static_cast<double>(result_scale*normalized);
}

bool peak_relevant_to_interval(complex<double> z, double lower, double upper)
{
   if (lower <= z.real() && z.real() <= upper)
      return true;
   const double distance = z.real() < lower ? lower - z.real() : z.real() - upper;
   return distance <= 8192.0*abs(z.imag());
}

double peak_distance_to_interval(complex<double> z, double lower, double upper)
{
   if (z.real() < lower)
      return lower - z.real();
   if (z.real() > upper)
      return z.real() - upper;
   return 0.0;
}

long double restricted_single_scale(long double width, int power,
                                     long double sign)
{
   const long double spectral_scale = 1.0L/(PI_L*width);
   return sign*width*powl(spectral_scale, power);
}

long double configure_restricted_single_scale(restricted_single_params &params)
{
   const long double sign = params.signed_imaginary < 0.0L && params.power % 2 == 1
      ? -1.0L : 1.0L;
   const long double full_scale =
      restricted_single_scale(params.width, params.power, sign);
   params.integrand_scale = restricted_integrand_scale(full_scale);
   return full_scale/params.integrand_scale;
}

double integrate_restricted_odd_pair(int kernel, int power,
                                     complex<double> z, double upper)
{
   const long double center = abs(z.real());
   const long double width = abs(z.imag());
   const bool direct = center >= upper && center - upper == center;
   const long double lower_offset = -center;
   const long double upper_offset = upper - center;
   const long double angle_origin = atan2l(lower_offset, width);
   const double angle_span = direct ? upper : static_cast<double>(
      restricted_angle_span(lower_offset, upper_offset, width));
   if (!(angle_span > 0.0))
      return 0.0;
   const long double imaginary_sign = z.imag() < 0.0 && power % 2 == 1
      ? -1.0L : 1.0L;
   const long double center_sign = z.real() < 0.0 ? -1.0L : 1.0L;
   const long double full_scale = direct
      ? center_sign*imaginary_sign*powl(1.0L/(PI_L*width), power)
      : center_sign*restricted_single_scale(width, power, imaginary_sign);
   const long double integrand_scale = restricted_integrand_scale(full_scale);
   const long double result_scale = full_scale/integrand_scale;
   restricted_odd_pair_params params = {
      kernel, power, center, width, angle_origin, integrand_scale, !direct
   };
   gsl_function function;
   function.function = &restricted_odd_pair_integrand;
   function.params = &params;
   const long double normalized = integrate_restricted_qag(function, 0.0, angle_span,
      "Symmetry-paired restricted epsilon integration", result_scale);
   return static_cast<double>(result_scale*normalized);
}

double integrate_restricted_single_exterior_log(gsl_function &function,
                                                 restricted_single_params &params,
                                                 double lower, double upper)
{
   const int direction = params.center < lower ? 1 : -1;
   const long double near_distance = direction > 0
      ? (lower - params.center)/params.width
      : (params.center - upper)/params.width;
   const long double far_distance = direction > 0
      ? (upper - params.center)/params.width
      : (params.center - lower)/params.width;
   const double span = static_cast<double>(logl(far_distance/near_distance));
   if (!(span > 0.0))
      return integrate_restricted_qag(function, lower, upper,
         "Restricted exterior epsilon integration");

   params.logarithmic = true;
   params.transformed_near = near_distance;
   params.transformed_direction = direction;
   const long double result_scale = configure_restricted_single_scale(params);
   vector<double> points;
   points.push_back(0.0);
   points.push_back(span);
   const double logarithmic_step = log(8.0);
   for (double point = logarithmic_step; point < span; point += logarithmic_step)
      points.push_back(point);
   points = filter_integration_points(points, 0.0, span);
   const long double normalized = integrate_restricted_piecewise_qag(function, points,
      "Log-transformed restricted epsilon integration", result_scale);
   return static_cast<double>(result_scale*normalized);
}

double integrate_restricted_single_band_edge(gsl_function &function,
                                             restricted_single_params &params,
                                             double lower, double upper)
{
   const int direction = params.center < 0.0L ? 1 : -1;
   const long double far_distance = (upper - lower)/params.width;
   params.transformed_direction = direction;
   params.scaled = true;
   const long double result_scale = configure_restricted_single_scale(params);
   const double scaled_upper = static_cast<double>(sqrtl(min(1.0L, far_distance)));
   long double normalized = integrate_restricted_qag(function, 0.0, scaled_upper,
      "Band-edge transformed restricted epsilon integration", result_scale);

   if (far_distance > 1.0L) {
      params.scaled = false;
      params.logarithmic = true;
      params.transformed_near = 1.0L;
      const double span = static_cast<double>(logl(far_distance));
      vector<double> points;
      points.push_back(0.0);
      points.push_back(span);
      const double logarithmic_step = log(8.0);
      for (double point = logarithmic_step; point < span; point += logarithmic_step)
         points.push_back(point);
      points = filter_integration_points(points, 0.0, span);
      normalized += integrate_restricted_piecewise_qag(function, points,
         "Log-transformed band-edge restricted epsilon integration", result_scale);
   }

   return static_cast<double>(result_scale*normalized);
}

} // namespace

double J_restricted_single(int kernel, int power, complex<double> z,
                           double lower, double upper)
{
   if (!(lower < upper))
      return 0.0;
   if (power < 0 || !isfinite(z.real()) || !isfinite(z.imag())) {
      cerr << "Invalid restricted kernel request: m=" << kernel
           << ", n=" << power << ", z=" << z << endl;
      exit(EXIT_FAILURE);
   }
   if (kernel != 0 && kernel % 2 == 0 && lower == -upper &&
       (power == 0 || z.real() == 0.0))
      return 0.0;

   const long double width = abs(z.imag());
   restricted_single_params params = {kernel, power, z.real(), width,
                                       z.imag(), false, false, false,
                                       0.0L, 0.0L, 0, 1.0L};
   gsl_function function;
   function.function = &restricted_single_integrand;
   function.params = &params;
   if (power == 0)
      return integrate_restricted_single_direct(function, kernel, lower, upper,
					 "Restricted epsilon integration");
   if (width == 0.0L) {
      cerr << "Restricted spectral kernel requires a nonzero linewidth." << endl;
      exit(EXIT_FAILURE);
   }
   const long double absolute_center = abs(z.real());
   const bool unresolved_exterior = absolute_center >= upper &&
      absolute_center - upper == absolute_center;
   if (kernel != 0 && kernel % 2 == 0 && lower == -upper &&
       (absolute_center <= 16.0L*width || unresolved_exterior))
      return integrate_restricted_odd_pair(kernel, power, z, upper);

   const long double distance = z.real() < lower ? lower - z.real()
      : (z.real() > upper ? z.real() - upper : 0.0L);
   if (distance > 0.0L)
      return integrate_restricted_single_exterior_log(function, params, lower, upper);
   if (3 <= kernel && kernel <= 6 &&
       ((params.center == -1.0L && lower == -1.0) ||
        (params.center == 1.0L && upper == 1.0)))
      return integrate_restricted_single_band_edge(function, params, lower, upper);

   const long double lower_offset = lower - params.center;
   const long double upper_offset = upper - params.center;
   const long double theta_lower = atan2l(lower_offset, width);
   const double transformed_lower = 0.0;
   const double transformed_upper = static_cast<double>(
      restricted_angle_span(lower_offset, upper_offset, width));
   if (!(transformed_lower < transformed_upper))
      return integrate_restricted_single_direct(function, kernel, lower, upper,
					 "Restricted exterior epsilon integration");

   params.transformed = true;
   params.transformed_origin = theta_lower;
   const long double result_scale = configure_restricted_single_scale(params);
   const long double normalized = integrate_restricted_qag(function,
      transformed_lower, transformed_upper, "Transformed restricted epsilon integration",
      result_scale);
   return static_cast<double>(result_scale*normalized);
}

double J_restricted_optical(int kernel, complex<double> z1, complex<double> z2,
                            double lower, double upper)
{
   if (!(lower < upper))
      return 0.0;
   if (kernel != 0 && kernel % 2 == 0 && lower == -upper &&
       z1.real() == -z2.real() && abs(z1.imag()) == abs(z2.imag()))
      return 0.0;
   if (z1 == z2)
      return J_restricted_single(kernel, 2, z1, lower, upper);
   if (z1.imag() == 0.0 || z2.imag() == 0.0) {
      cerr << "Restricted optical kernel requires nonzero linewidths." << endl;
      exit(EXIT_FAILURE);
   }

   const bool first_relevant = peak_relevant_to_interval(z1, lower, upper);
   const bool second_relevant = peak_relevant_to_interval(z2, lower, upper);
   const bool first_much_narrower = abs(z1.imag()) < abs(z2.imag())/16.0;
   const bool second_much_narrower = abs(z2.imag()) < abs(z1.imag())/16.0;

   if (first_relevant && second_relevant &&
       abs(z1.real() - z2.real()) > 8.0*max(abs(z1.imag()), abs(z2.imag()))) {
      const bool first_is_left = z1.real() <= z2.real();
      const complex<double> left = first_is_left ? z1 : z2;
      const complex<double> right = first_is_left ? z2 : z1;
      const long double split = 0.5L*(static_cast<long double>(left.real()) +
                                     static_cast<long double>(right.real()));
      if (lower < split && split < upper && left.real() < split && split < right.real()) {
	 return integrate_restricted_optical_transformed(kernel, left, right,
	                                                 static_cast<long double>(lower), split)
	    + integrate_restricted_optical_transformed(kernel, right, left, split,
	                                                static_cast<long double>(upper));
      }
   }

   complex<double> anchor;
   complex<double> other;
   const double first_distance = peak_distance_to_interval(z1, lower, upper);
   const double second_distance = peak_distance_to_interval(z2, lower, upper);
   if (first_much_narrower || (!second_much_narrower &&
       (first_distance < second_distance ||
	(first_distance == second_distance && abs(z1.imag()) <= abs(z2.imag()))))) {
      anchor = z1;
      other = z2;
   } else {
      anchor = z2;
      other = z1;
   }
   return integrate_restricted_optical_transformed(kernel, anchor, other, lower, upper);
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
   if (epsilon_window_limit > 0.0) {
      const epsilon_interval interval = restricted_epsilon_interval(m);
      if (interval.empty)
	 return 0.0;
      if (interval.truncated)
	 return J_restricted_optical(m, z1, z2, interval.lower, interval.upper);
   }

   if (m == 0)
      return J_0_optical(z1, z2);
   if (z1 == z2)
      return J_builtin(m, 2, z1);

   const complex<double> cross = hilbert_divided_difference(m, conj(z1), z2);
   const complex<double> same = hilbert_divided_difference(m, z1, z2);
   return (cross.real() - same.real())/(2.0*M_PI*M_PI);
}

double clipped_im_sigma(double value)
{
   return value > -sigma_clip ? -sigma_clip : value;
}

complex<double> effective_frequency(double omega)
{
   double sigma_re;
   double sigma_im;
   if (omega_min <= omega && omega <= omega_max) {
      sigma_re = gsl_spline_eval(spline_reSigma, omega, acc_reSigma);
      // Higher-order splines can overshoot causal input knots.
      sigma_im = clipped_im_sigma(
         gsl_spline_eval(spline_imSigma, omega, acc_imSigma));
   } else {
      sigma_im = -EPSILON;
      if (omega < omega_min)
	 sigma_re = reSigma_asymp_neg;
      if (omega > omega_max)
	 sigma_re = reSigma_asymp_pos;
   }

   return omega + mu - sigma_re - sigma_im*I;
}

namespace {

double real_effective_frequency(double omega)
{
   double sigma_re;
   if (omega_min <= omega && omega <= omega_max)
      sigma_re = gsl_spline_eval(spline_reSigma, omega, acc_reSigma);
   else
      sigma_re = omega < omega_min ? reSigma_asymp_neg : reSigma_asymp_pos;
   return omega + mu - sigma_re;
}

double restricted_crossing_function(double omega, double shift, double edge)
{
   return real_effective_frequency(omega + shift) - edge;
}

double real_effective_derivative(double omega)
{
   if (omega_min <= omega && omega <= omega_max)
      return 1.0 - gsl_spline_eval_deriv(spline_reSigma, omega, acc_reSigma);
   return 1.0;
}

double real_effective_second_derivative(double omega)
{
   if (omega_min <= omega && omega <= omega_max)
      return -gsl_spline_eval_deriv2(spline_reSigma, omega, acc_reSigma);
   return 0.0;
}

double bisect_effective_derivative(double lower, double upper)
{
   double lower_value = real_effective_derivative(lower);
   for (int iteration = 0; iteration < 80; ++iteration) {
      const double middle = 0.5*(lower + upper);
      if (middle == lower || middle == upper)
         break;
      const double middle_value = real_effective_derivative(middle);
      if (middle_value == 0.0)
         return middle;
      if (signbit(lower_value) != signbit(middle_value)) {
         upper = middle;
      } else {
         lower = middle;
         lower_value = middle_value;
      }
   }
   return 0.5*(lower + upper);
}

double bisect_effective_second_derivative(double lower, double upper)
{
   double lower_value = real_effective_second_derivative(lower);
   for (int iteration = 0; iteration < 80; ++iteration) {
      const double middle = 0.5*(lower + upper);
      if (middle == lower || middle == upper)
         break;
      const double middle_value = real_effective_second_derivative(middle);
      if (middle_value == 0.0)
         return middle;
      if (signbit(lower_value) != signbit(middle_value)) {
         upper = middle;
      } else {
         lower = middle;
         lower_value = middle_value;
      }
   }
   return 0.5*(lower + upper);
}

vector<double> add_effective_frequency_extrema(vector<double> &partitions, double shift)
{
   vector<double> extrema;
   for (size_t interval = 1; interval < partitions.size(); ++interval) {
      const double outer_lower = partitions[interval - 1];
      const double outer_upper = partitions[interval];
      const double actual_middle = 0.5*(outer_lower + outer_upper) + shift;
      if (!(omega_min < actual_middle && actual_middle < omega_max))
         continue;

      const double actual_lower = nextafter(outer_lower + shift, outer_upper + shift);
      const double actual_upper = nextafter(outer_upper + shift, outer_lower + shift);
      vector<double> derivative_partitions;
      derivative_partitions.push_back(actual_lower);
      const double second_lower = real_effective_second_derivative(actual_lower);
      const double second_upper = real_effective_second_derivative(actual_upper);
      if (second_lower == 0.0) {
         derivative_partitions.push_back(actual_lower);
      } else if (second_upper == 0.0) {
         derivative_partitions.push_back(actual_upper);
      } else if (signbit(second_lower) != signbit(second_upper)) {
         derivative_partitions.push_back(
            bisect_effective_second_derivative(actual_lower, actual_upper));
      }
      derivative_partitions.push_back(actual_upper);
      sort(derivative_partitions.begin(), derivative_partitions.end());
      derivative_partitions.erase(
         unique(derivative_partitions.begin(), derivative_partitions.end()),
         derivative_partitions.end());

      for (size_t part = 1; part < derivative_partitions.size(); ++part) {
         const double lower = derivative_partitions[part - 1];
         const double upper = derivative_partitions[part];
         const double lower_value = real_effective_derivative(lower);
         const double upper_value = real_effective_derivative(upper);
         if (lower_value == 0.0)
            extrema.push_back(lower - shift);
         if (upper_value == 0.0)
            extrema.push_back(upper - shift);
         if (lower_value != 0.0 && upper_value != 0.0 &&
             signbit(lower_value) != signbit(upper_value))
            extrema.push_back(bisect_effective_derivative(lower, upper) - shift);
      }
   }
   partitions.insert(partitions.end(), extrema.begin(), extrema.end());
   sort(partitions.begin(), partitions.end());
   partitions.erase(unique(partitions.begin(), partitions.end()), partitions.end());
   return extrema;
}

double bisect_restricted_crossing(double lower, double upper,
                                  double shift, double edge)
{
   double lower_value = restricted_crossing_function(lower, shift, edge);
   for (int iteration = 0; iteration < 80; ++iteration) {
      const double middle = 0.5*(lower + upper);
      if (middle == lower || middle == upper)
         break;
      const double middle_value = restricted_crossing_function(middle, shift, edge);
      if (middle_value == 0.0)
         return middle;
      if (signbit(lower_value) != signbit(middle_value)) {
         upper = middle;
      } else {
         lower = middle;
         lower_value = middle_value;
      }
   }
   return 0.5*(lower + upper);
}

void add_restricted_crossing_neighborhood(vector<double> &points,
                                          double crossing,
                                          double lower, double upper,
                                          double shift)
{
   points.push_back(crossing);
   const double frequency = crossing + shift;
   const double linewidth = abs(effective_frequency(frequency).imag());
   const double derivative_step = max(linewidth,
      1.0e-6*max(1.0, abs(frequency)));
   const double slope = (real_effective_frequency(frequency + derivative_step)
      - real_effective_frequency(frequency - derivative_step))/(2.0*derivative_step);
   double crossing_width = linewidth/max(1.0e-3, abs(slope));
   const double curvature = abs(real_effective_second_derivative(frequency));
   if (curvature > 0.0 &&
       abs(slope) < sqrt(2.0*linewidth*curvature))
      crossing_width = max(crossing_width, sqrt(2.0*linewidth/curvature));
   double distance = crossing_width;
   for (int scale = 0; scale < 7; ++scale) {
      if (lower < crossing - distance && crossing - distance < upper)
         points.push_back(crossing - distance);
      if (lower < crossing + distance && crossing + distance < upper)
         points.push_back(crossing + distance);
      distance *= 4.0;
   }
}

void add_restricted_crossings(vector<double> &points,
                              double lower, double upper,
                              double shift, double edge)
{
   vector<double> partitions;
   partitions.push_back(lower);
   partitions.push_back(upper);
   for (vector<double>::const_iterator knot = sigma_omega_knots.begin();
        knot != sigma_omega_knots.end(); ++knot) {
      const double shifted_knot = *knot - shift;
      if (lower < shifted_knot && shifted_knot < upper) {
         partitions.push_back(shifted_knot);
         if (omega_min < *knot && *knot < omega_max) {
            const double left_slope = real_effective_derivative(
               nextafter(*knot, -numeric_limits<double>::infinity()));
            const double right_slope = real_effective_derivative(
               nextafter(*knot, numeric_limits<double>::infinity()));
            if (left_slope == 0.0 || right_slope == 0.0 ||
                signbit(left_slope) != signbit(right_slope))
               add_restricted_crossing_neighborhood(points, shifted_knot,
                                                     lower, upper, shift);
         }
      }
   }
   sort(partitions.begin(), partitions.end());
   partitions.erase(unique(partitions.begin(), partitions.end()), partitions.end());
   const vector<double> extrema = add_effective_frequency_extrema(partitions, shift);
   for (vector<double>::const_iterator extremum = extrema.begin();
        extremum != extrema.end(); ++extremum)
      add_restricted_crossing_neighborhood(points, *extremum,
                                           lower, upper, shift);

   for (size_t interval = 1; interval < partitions.size(); ++interval) {
      const double interval_lower = partitions[interval - 1];
      const double interval_upper = partitions[interval];
      const double lower_value = restricted_crossing_function(interval_lower, shift, edge);
      const double upper_value = restricted_crossing_function(interval_upper, shift, edge);
      if (lower_value == 0.0)
         add_restricted_crossing_neighborhood(points, interval_lower,
                                              lower, upper, shift);
      if (upper_value == 0.0)
         add_restricted_crossing_neighborhood(points, interval_upper,
                                              lower, upper, shift);
      if (lower_value != 0.0 && upper_value != 0.0 &&
          signbit(lower_value) != signbit(upper_value))
         add_restricted_crossing_neighborhood(points,
            bisect_restricted_crossing(interval_lower, interval_upper, shift, edge),
            lower, upper, shift);
   }
}

vector<double> restricted_outer_points(double lower, double upper,
                                        const epsilon_interval &interval,
                                        bool include_shifted)
{
   vector<double> points;
   points.push_back(lower);
   points.push_back(upper);
   add_restricted_crossings(points, lower, upper, 0.0, interval.lower);
   add_restricted_crossings(points, lower, upper, 0.0, interval.upper);
   if (include_shifted) {
      add_restricted_crossings(points, lower, upper, optical_frequency, interval.lower);
      add_restricted_crossings(points, lower, upper, optical_frequency, interval.upper);
   }
   return filter_integration_points(points, lower, upper);
}

bool restricted_outer_failure_acceptable(int status, gsl_function &function,
                                         const vector<double> &points,
                                         double result, double error)
{
   if (status != GSL_EROUND && status != GSL_ESING)
      return false;
   if (!isfinite(result))
      return false;
   double fine = 0.0;
   double check_difference = 0.0;
   for (size_t interval = 1; interval < points.size(); ++interval) {
      const restricted_fixed_result check = integrate_restricted_fixed_check(function,
         points[interval - 1], points[interval]);
      if (!isfinite(check.value) || !isfinite(check.difference))
         return false;
      fine += check.value;
      check_difference += check.difference;
   }
   const double independent_error = max(abs(result - fine), check_difference);
   if (!isfinite(independent_error) || independent_error < 0.0)
      return false;
   const double verified_error = status == GSL_EROUND &&
      isfinite(error) && error >= 0.0
      ? min(error, independent_error) : independent_error;
   return isfinite(verified_error) &&
      verified_error <= abs_error + rel_error*abs(result);
}

struct symmetric_outer_params {
   double (*function)(double, void *);
   void *params;
};

double symmetric_outer_integrand(double omega, void *raw_params)
{
   const symmetric_outer_params &params =
      *static_cast<symmetric_outer_params *>(raw_params);
   return params.function(omega, params.params) +
          params.function(-omega, params.params);
}

vector<double> symmetric_outer_points(const vector<double> &points, double upper)
{
   vector<double> symmetric_points;
   symmetric_points.push_back(0.0);
   symmetric_points.push_back(upper);
   for (vector<double>::const_iterator point = points.begin(); point != points.end(); ++point) {
      const double magnitude = abs(*point);
      if (0.0 < magnitude && magnitude < upper)
         symmetric_points.push_back(magnitude);
   }
   return filter_integration_points(symmetric_points, 0.0, upper);
}

} // namespace

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
   gsl_set_error_handler_off();
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
   sigma_omega_knots = omega;
   assert(omega_min < omega_max);

   reSigma_asymp_neg = data1[0][1]; // use this for omega<omega_min
   reSigma_asymp_pos = data1[N-1][1]; // use this for omega>omega_max
   // imSigma assumed to be zero outside the [omega_min:omega_max] interval
      
   for(int i=0; i<N; ++i) {
      imSigma.push_back(clipped_im_sigma(data2[i][1]));
   }

   particle_hole_symmetric_sigma = true;
   for (int i = 0; i < N; ++i) {
      const int mirrored = N - 1 - i;
      if (omega[i] != -omega[mirrored] ||
          reSigma[i] != -reSigma[mirrored] ||
          imSigma[i] != imSigma[mirrored]) {
         particle_hole_symmetric_sigma = false;
         break;
      }
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

void initialize_epsilon_window()
{
   if (epsilon_window_limit == 0.0)
      return;

   if (!(omega_min <= 0.0 && 0.0 <= omega_max)) {
      cerr << "Option -M requires the self-energy frequency grid to contain omega=0."
           << endl;
      exit(EXIT_FAILURE);
   }

   const double re_sigma_zero = gsl_spline_eval(spline_reSigma, 0.0, acc_reSigma);
   epsilon_window_center = mu - re_sigma_zero;
   const long double lower = static_cast<long double>(epsilon_window_center)
      - epsilon_window_limit;
   const long double upper = static_cast<long double>(epsilon_window_center)
      + epsilon_window_limit;
   const long double maximum = numeric_limits<double>::max();
   if (!isfinite(epsilon_window_center) || lower < -maximum || upper > maximum) {
      cerr << "The epsilon-window bounds are not representable." << endl;
      exit(EXIT_FAILURE);
   }
   epsilon_window_lower = static_cast<double>(lower);
   epsilon_window_upper = static_cast<double>(upper);
   if (!(epsilon_window_lower < epsilon_window_center &&
         epsilon_window_center < epsilon_window_upper)) {
      cerr << "The positive epsilon-window limit is too small at the resolved center."
           << endl;
      exit(EXIT_FAILURE);
   }
   if (epsilon_window_center != 0.0 &&
       epsilon_window_lower == -epsilon_window_upper) {
      cerr << "The epsilon-window bounds do not preserve the resolved center." << endl;
      exit(EXIT_FAILURE);
   }

   if (verbose) {
      cout << "ReSigma(0)=" << re_sigma_zero << endl;
      cout << "epsilon window=[" << epsilon_window_lower << ","
           << epsilon_window_upper << "]" << endl;
   }
}

void calc()
{
   // GSL integration
   const size_t ws_size = 1000;

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
   gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(ws_size);
   if (work_ptr == NULL) {
      cerr << "Unable to allocate frequency integration workspace." << endl;
      exit(EXIT_FAILURE);
   }
   double result = numeric_limits<double>::quiet_NaN();
   double error = numeric_limits<double>::quiet_NaN();
    
   gsl_function My_function;
   My_function.function = &integrand;
   void *params_ptr = &T;
   My_function.params = params_ptr;
    
   gsl_set_error_handler_off();
   const epsilon_interval interval = epsilon_window_limit > 0.0
      ? restricted_epsilon_interval(m)
      : epsilon_interval{0.0, 0.0, false, false};
   const bool parity_zero = epsilon_window_limit > 0.0 && interval.truncated &&
      f_type == f_derivative && particle_hole_symmetric_sigma &&
      epsilon_window_center == 0.0 && m != 0 && m % 2 == 0 && o % 2 == 0 &&
      interval.lower == -interval.upper;
   if (parity_zero || (epsilon_window_limit > 0.0 && interval.empty)) {
      result = 0.0;
      error = 0.0;
   } else if (epsilon_window_limit > 0.0 && interval.truncated) {
      vector<double> points;
      points.push_back(lower_limit);
      points.push_back(upper_limit);
      if (n != 0)
         points = restricted_outer_points(lower_limit, upper_limit, interval, false);
      gsl_function *integration_function = &My_function;
      symmetric_outer_params paired_params;
      gsl_function paired_function;
      if (f_type == f_derivative && particle_hole_symmetric_sigma &&
          m != 0 && m % 2 == 0 && o % 2 == 0 &&
          lower_limit == -upper_limit && interval.lower == -interval.upper) {
         paired_params.function = My_function.function;
         paired_params.params = My_function.params;
         paired_function.function = &symmetric_outer_integrand;
         paired_function.params = &paired_params;
         integration_function = &paired_function;
         points = symmetric_outer_points(points, upper_limit);
      }
      const size_t restricted_size = max(ws_size, 32*points.size());
      gsl_integration_workspace *restricted_work = restricted_size == ws_size
         ? work_ptr : gsl_integration_workspace_alloc(restricted_size);
      if (restricted_work == NULL) {
	 gsl_integration_workspace_free(work_ptr);
	 cerr << "Unable to allocate restricted frequency integration workspace." << endl;
	 exit(EXIT_FAILURE);
      }
      const int status = gsl_integration_qagp(integration_function,
                                             &points[0], points.size(),
					     abs_error, rel_error, restricted_size,
					     restricted_work, &result, &error);
      if (restricted_work != work_ptr)
	 gsl_integration_workspace_free(restricted_work);
      const bool acceptable_failure = restricted_outer_failure_acceptable(
         status, *integration_function, points, result, error);
      const bool successful = status == GSL_SUCCESS && isfinite(result) &&
         isfinite(error) && error >= 0.0;
      if (!successful && !acceptable_failure &&
          !integration_result_acceptable("Restricted frequency integration",
                                         status, result, error)) {
	 gsl_integration_workspace_free(work_ptr);
	 exit(EXIT_FAILURE);
      }
      if (acceptable_failure && !quiet_warnings)
	 cerr << "Warning: restricted frequency integration returned "
	      << gsl_strerror(status) << "; independent quadrature agreed within tolerance."
	      << endl;
   } else {
      const int status = gsl_integration_qag(
	 &My_function,       // integrand function
	 lower_limit,        // lower integration boundary
	 upper_limit,        // upper integration boundary
	 abs_error,          // preferred absolute error
	 rel_error,          // preferred relative error
	 ws_size,            // size of workspace
	 key,                // Gauss-Kronrod rule
	 work_ptr,           // integration workspace
	 &result,            // final approximation
	 &error);            // estimate of absolute error
      if (!integration_result_acceptable("Frequency integration",
                                         status, result, error)) {
	 gsl_integration_workspace_free(work_ptr);
	 exit(EXIT_FAILURE);
      }
   }
   gsl_integration_workspace_free(work_ptr);
    
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
   bool independently_verified = false;
   const epsilon_interval interval = epsilon_window_limit > 0.0
      ? restricted_epsilon_interval(m)
      : epsilon_interval{0.0, 0.0, false, false};
   if (epsilon_window_limit > 0.0 && interval.empty) {
      result = 0.0;
      error = 0.0;
   } else if (epsilon_window_limit > 0.0 && interval.truncated) {
      vector<double> points = restricted_outer_points(lower_limit, upper_limit,
                                                      interval, true);
      points.push_back(-optical_frequency);
      points.push_back(0.0);
      points = filter_integration_points(points, lower_limit, upper_limit);
      const size_t restricted_size = max(ws_size, 32*points.size());
      gsl_integration_workspace *restricted_work = restricted_size == ws_size
         ? work : gsl_integration_workspace_alloc(restricted_size);
      if (restricted_work == NULL) {
	 gsl_integration_workspace_free(work);
	 cerr << "Unable to allocate restricted optical frequency integration workspace."
	      << endl;
	 exit(EXIT_FAILURE);
      }
      const int status = gsl_integration_qagp(&F, &points[0], points.size(),
					     abs_error, rel_error, restricted_size,
					     restricted_work, &result, &error);
      if (restricted_work != work)
	 gsl_integration_workspace_free(restricted_work);
      const bool acceptable_failure = restricted_outer_failure_acceptable(
         status, F, points, result, error);
      const bool successful = status == GSL_SUCCESS && isfinite(result) &&
         isfinite(error) && error >= 0.0;
      if (!successful && !acceptable_failure &&
          !integration_result_acceptable("Restricted optical frequency integration",
                                         status, result, error)) {
	 gsl_integration_workspace_free(work);
	 exit(EXIT_FAILURE);
      }
      independently_verified = acceptable_failure;
   } else {
      for (int segment = 0; segment < 3; ++segment) {
	 double segment_result = numeric_limits<double>::quiet_NaN();
	 double segment_error = numeric_limits<double>::quiet_NaN();
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
	 if (!integration_result_acceptable("Optical frequency integration",
	                                    status, segment_result, segment_error)) {
	    gsl_integration_workspace_free(work);
	    exit(EXIT_FAILURE);
	 }
	 result += segment_result;
	 error += segment_error;
      }
   }
   gsl_integration_workspace_free(work);

   if (!integration_result_acceptable("Optical frequency integration total",
                                      GSL_SUCCESS, result, error))
      exit(EXIT_FAILURE);

   if (independently_verified && !quiet_warnings) {
      cerr << "Warning: restricted optical frequency integration required independent "
           << "verification; estimated GSL error=" << error << "." << endl;
   }

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
   initialize_epsilon_window();
   
   if (calcdos)
      calc_DOS();
   else if (optical_mode && optical_frequency > 0.0)
      calc_optical();
   else
      calc();
    
   return 0;
}
#endif
