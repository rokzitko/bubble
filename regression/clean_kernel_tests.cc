#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

double J_builtin(int kernel, int power, std::complex<double> z);
double J_restricted_single(int kernel, int power, std::complex<double> z,
                           double lower, double upper);
double J_restricted_optical(int kernel, std::complex<double> z1,
                            std::complex<double> z2,
                            double lower, double upper);
void add_sigma_transition_points(std::vector<double> &points,
                                 double lower, double upper,
                                 bool include_shifted);
extern bool allow_legacy_self_energy_extrapolation;
extern double omega_min;
extern double omega_max;
extern double optical_frequency;

int main(int argc, char **argv)
{
   if (argc > 4) {
      std::cerr << "Usage: clean_kernel_tests [clean-references [restricted-references [restricted-optical-references]]]\n";
      return EXIT_FAILURE;
   }
   const char *filename = argc >= 2 ? argv[1] : "clean_kernel_references.dat";
   std::ifstream input(filename);
   if (!input.is_open()) {
      std::cerr << "Cannot open " << filename << "\n";
      return EXIT_FAILURE;
   }

   int passed = 0;
   int failed = 0;
   int reference_rows = 0;
   std::string line;
   while (std::getline(input, line)) {
      if (line.empty() || line[0] == '#')
         continue;

      int kernel;
      int power;
      double x;
      double y;
      double expected;
      double relative_tolerance;
      double absolute_tolerance;
      std::istringstream values(line);
      if (!(values >> kernel >> power >> x >> y >> expected
                   >> relative_tolerance >> absolute_tolerance)) {
         std::cerr << "Invalid reference row: " << line << "\n";
         return EXIT_FAILURE;
      }
      ++reference_rows;

      const double value = J_builtin(kernel, power, std::complex<double>(x, y));
      const double tolerance = std::max(absolute_tolerance,
         relative_tolerance*std::abs(expected));
      if (std::isfinite(value) && std::abs(value - expected) <= tolerance) {
         ++passed;
      } else {
         ++failed;
         std::cerr << std::setprecision(17)
                   << "FAILED: m=" << kernel << " n=" << power
                   << " z=(" << x << "," << y << ")"
                   << " expected=" << expected << " value=" << value
                   << " tolerance=" << tolerance << "\n";
      }

      const double mirrored = J_builtin(kernel, power, std::complex<double>(-x, y));
      const double parity_value = kernel % 2 == 0 ? -value : value;
      if (mirrored == parity_value) {
         ++passed;
      } else {
         ++failed;
         std::cerr << std::setprecision(17)
                   << "FAILED parity: m=" << kernel << " n=" << power
                   << " x=" << x << " positive=" << value
                   << " negative=" << mirrored << "\n";
      }
   }

   if (reference_rows == 0) {
      std::cerr << "No clean-kernel reference rows found in " << filename << "\n";
      return EXIT_FAILURE;
   }

   for (int kernel = 2; kernel <= 8; kernel += 2) {
      for (int power = 2; power <= 3; ++power) {
         const double value = J_builtin(kernel, power, std::complex<double>(0.0, 1e-12));
         if (value == 0.0) {
            ++passed;
         } else {
            ++failed;
            std::cerr << "FAILED exact odd-kernel zero: m=" << kernel
                      << " n=" << power << " value=" << value << "\n";
         }
      }
   }

   const char *restricted_filename = argc >= 3
      ? argv[2] : "restricted_kernel_references.dat";
   std::ifstream restricted_input(restricted_filename);
   if (!restricted_input.is_open()) {
      std::cerr << "Cannot open " << restricted_filename << "\n";
      return EXIT_FAILURE;
   }
   int restricted_rows = 0;
   while (std::getline(restricted_input, line)) {
      if (line.empty() || line[0] == '#')
         continue;

      int kernel;
      int power;
      double x;
      double y;
      double lower;
      double upper;
      double expected;
      double relative_tolerance;
      double absolute_tolerance;
      std::istringstream values(line);
      if (!(values >> kernel >> power >> x >> y >> lower >> upper >> expected
                   >> relative_tolerance >> absolute_tolerance)) {
         std::cerr << "Invalid restricted reference row: " << line << "\n";
         return EXIT_FAILURE;
      }
      ++restricted_rows;

      const double value = J_restricted_single(kernel, power,
         std::complex<double>(x, y), lower, upper);
      const double tolerance = std::max(absolute_tolerance,
         relative_tolerance*std::abs(expected));
      if (std::isfinite(value) && std::abs(value - expected) <= tolerance) {
         ++passed;
      } else {
         ++failed;
         std::cerr << std::setprecision(17)
                   << "FAILED restricted: m=" << kernel << " n=" << power
                   << " z=(" << x << "," << y << ") interval=["
                   << lower << "," << upper << "] expected=" << expected
                   << " value=" << value << " tolerance=" << tolerance << "\n";
      }

      const double mirrored = J_restricted_single(kernel, power,
         std::complex<double>(-x, y), -upper, -lower);
      const double parity_value = kernel % 2 == 0 ? -value : value;
      if (std::abs(mirrored - parity_value) <= tolerance) {
         ++passed;
      } else {
         ++failed;
         std::cerr << std::setprecision(17)
                   << "FAILED restricted parity: m=" << kernel << " n=" << power
                   << " value=" << value << " mirrored=" << mirrored << "\n";
      }
   }
   if (restricted_rows == 0) {
      std::cerr << "No restricted-kernel reference rows found in "
                << restricted_filename << "\n";
      return EXIT_FAILURE;
   }

   const char *optical_filename = argc == 4
      ? argv[3] : "restricted_optical_references.dat";
   std::ifstream optical_input(optical_filename);
   if (!optical_input.is_open()) {
      std::cerr << "Cannot open " << optical_filename << "\n";
      return EXIT_FAILURE;
   }
   int optical_rows = 0;
   while (std::getline(optical_input, line)) {
      if (line.empty() || line[0] == '#')
         continue;

      int kernel;
      double x1;
      double y1;
      double x2;
      double y2;
      double lower;
      double upper;
      double expected;
      double relative_tolerance;
      double absolute_tolerance;
      std::istringstream values(line);
      if (!(values >> kernel >> x1 >> y1 >> x2 >> y2 >> lower >> upper
                   >> expected >> relative_tolerance >> absolute_tolerance)) {
         std::cerr << "Invalid restricted optical reference row: " << line << "\n";
         return EXIT_FAILURE;
      }
      ++optical_rows;

      const double value = J_restricted_optical(kernel,
         std::complex<double>(x1, y1), std::complex<double>(x2, y2), lower, upper);
      const double tolerance = std::max(absolute_tolerance,
         relative_tolerance*std::abs(expected));
      if (std::isfinite(value) && std::abs(value - expected) <= tolerance) {
         ++passed;
      } else {
         ++failed;
         std::cerr << std::setprecision(17)
                   << "FAILED restricted optical: m=" << kernel
                   << " z1=(" << x1 << "," << y1 << ") z2=("
                   << x2 << "," << y2 << ") interval=[" << lower << ","
                   << upper << "] expected=" << expected << " value=" << value
                   << " tolerance=" << tolerance << "\n";
      }

      const double mirrored = J_restricted_optical(kernel,
         std::complex<double>(-x1, y1), std::complex<double>(-x2, y2),
         -upper, -lower);
      const double parity_value = kernel % 2 == 0 ? -value : value;
      if (std::abs(mirrored - parity_value) <= tolerance) {
         ++passed;
      } else {
         ++failed;
         std::cerr << std::setprecision(17)
                   << "FAILED restricted optical parity: m=" << kernel
                   << " value=" << value << " mirrored=" << mirrored << "\n";
      }
   }
   if (optical_rows == 0) {
      std::cerr << "No restricted optical reference rows found in "
                << optical_filename << "\n";
      return EXIT_FAILURE;
   }

   const bool saved_legacy_extrapolation = allow_legacy_self_energy_extrapolation;
   const double saved_omega_min = omega_min;
   const double saved_omega_max = omega_max;
   const double saved_optical_frequency = optical_frequency;
   allow_legacy_self_energy_extrapolation = true;
   omega_min = -1.0;
   omega_max = 1.0;
   optical_frequency = 0.3;

   std::vector<double> transition_points;
   add_sigma_transition_points(transition_points, -1.8, 1.5, true);
   std::sort(transition_points.begin(), transition_points.end());
   const double shifted_min = static_cast<double>(
      static_cast<long double>(omega_min) - optical_frequency);
   const double shifted_max = static_cast<double>(
      static_cast<long double>(omega_max) - optical_frequency);
   const std::vector<double> expected_transitions = {
      shifted_min, omega_min, shifted_max, omega_max
   };
   if (transition_points == expected_transitions) {
      ++passed;
   } else {
      ++failed;
      std::cerr << "FAILED shifted self-energy transition points\n";
   }

   transition_points.clear();
   add_sigma_transition_points(transition_points, -1.8, 1.5, false);
   std::sort(transition_points.begin(), transition_points.end());
   const std::vector<double> expected_unshifted = {omega_min, omega_max};
   if (transition_points == expected_unshifted) {
      ++passed;
   } else {
      ++failed;
      std::cerr << "FAILED unshifted self-energy transition points\n";
   }

   allow_legacy_self_energy_extrapolation = false;
   transition_points.clear();
   add_sigma_transition_points(transition_points, -1.8, 1.5, true);
   if (transition_points.empty()) {
      ++passed;
   } else {
      ++failed;
      std::cerr << "FAILED strict self-energy transition points\n";
   }

   allow_legacy_self_energy_extrapolation = saved_legacy_extrapolation;
   omega_min = saved_omega_min;
   omega_max = saved_omega_max;
   optical_frequency = saved_optical_frequency;

   std::cout << "CLEAN KERNEL OK=" << passed << " FAILED=" << failed << "\n";
   return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
