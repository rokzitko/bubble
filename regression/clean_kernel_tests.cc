#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

double J_builtin(int kernel, int power, std::complex<double> z);
double K_builtin(int kernel, std::complex<double> z1,
                 std::complex<double> z2);
double J_restricted_single(int kernel, int power, std::complex<double> z,
                           double lower, double upper);
double J_restricted_optical(int kernel, std::complex<double> z1,
                            std::complex<double> z2,
                            double lower, double upper);
void add_sigma_transition_points(std::vector<double> &points,
                                 double lower, double upper,
                                 bool include_shifted);
void add_spectral_frequency_transition_points(std::vector<double> &points,
                                               double lower, double upper,
                                               bool include_shifted);
extern bool allow_legacy_self_energy_extrapolation;
extern double omega_min;
extern double omega_max;
extern double optical_frequency;
extern double spectral_frequency_window_half_width;

int main(int argc, char **argv)
{
   if (argc > 5) {
      std::cerr << "Usage: clean_kernel_tests [clean-references [restricted-references [restricted-optical-references [clean-optical-references]]]]\n";
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
   std::set<std::tuple<int, int, double, double> > clean_keys;
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
      std::string extra;
      std::istringstream values(line);
      if (!(values >> kernel >> power >> x >> y >> expected
                   >> relative_tolerance >> absolute_tolerance) || values >> extra) {
         std::cerr << "Invalid reference row: " << line << "\n";
         return EXIT_FAILURE;
      }
      if (!clean_keys.insert(std::make_tuple(kernel, power, x, y)).second) {
         std::cerr << "Duplicate clean-kernel reference row: " << line << "\n";
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

      const double conjugate_side = J_builtin(kernel, power,
         std::complex<double>(x, -y));
      const double spectral_parity_value = power % 2 == 0 ? value : -value;
      if (conjugate_side == spectral_parity_value) {
         ++passed;
      } else {
         ++failed;
         std::cerr << std::setprecision(17)
                   << "FAILED linewidth sign: m=" << kernel << " n=" << power
                   << " positive=" << value << " negative=" << conjugate_side
                   << "\n";
      }

   }

   if (reference_rows != 370) {
      std::cerr << "Expected 370 clean-kernel reference rows in " << filename
                << ", found " << reference_rows << "\n";
      return EXIT_FAILURE;
   }

   for (int kernel = 2; kernel <= 8; kernel += 2) {
      for (int power = 1; power <= 3; ++power) {
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

   const double extreme_linewidth = 1e308;
   const double scaled_extreme_values[] = {
      2.0/std::acos(-1.0), 0.5, 0.375, 1.0/std::sqrt(2.0*std::acos(-1.0))
   };
   for (int index = 0; index < 4; ++index) {
      const int kernel = 2*index + 1;
      const double value = J_builtin(kernel, 1,
         std::complex<double>(0.37, extreme_linewidth));
      const double scaled = value*extreme_linewidth;
      const double expected = scaled_extreme_values[index];
      if (std::isfinite(scaled) && std::abs(scaled - expected) <= 1e-12) {
         ++passed;
      } else {
         ++failed;
         std::cerr << std::setprecision(17)
                   << "FAILED extreme linewidth: m=" << kernel
                   << " scaled=" << scaled << " expected=" << expected << "\n";
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
   std::set<std::tuple<int, int, double, double, double, double> > restricted_keys;
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
      std::string extra;
      std::istringstream values(line);
      if (!(values >> kernel >> power >> x >> y >> lower >> upper >> expected
                   >> relative_tolerance >> absolute_tolerance) || values >> extra) {
         std::cerr << "Invalid restricted reference row: " << line << "\n";
         return EXIT_FAILURE;
      }
      if (!restricted_keys.insert(
            std::make_tuple(kernel, power, x, y, lower, upper)).second) {
         std::cerr << "Duplicate restricted reference row: " << line << "\n";
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
   if (restricted_rows != 128) {
      std::cerr << "Expected 128 restricted-kernel reference rows in "
                << restricted_filename << ", found " << restricted_rows << "\n";
      return EXIT_FAILURE;
   }

   const char *optical_filename = argc >= 4
      ? argv[3] : "restricted_optical_references.dat";
   std::ifstream optical_input(optical_filename);
   if (!optical_input.is_open()) {
      std::cerr << "Cannot open " << optical_filename << "\n";
      return EXIT_FAILURE;
   }
   int optical_rows = 0;
   std::set<std::tuple<int, double, double, double, double, double, double> >
      optical_keys;
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
      std::string extra;
      std::istringstream values(line);
      if (!(values >> kernel >> x1 >> y1 >> x2 >> y2 >> lower >> upper
                   >> expected >> relative_tolerance >> absolute_tolerance) ||
          values >> extra) {
         std::cerr << "Invalid restricted optical reference row: " << line << "\n";
         return EXIT_FAILURE;
      }
      if (!optical_keys.insert(
            std::make_tuple(kernel, x1, y1, x2, y2, lower, upper)).second) {
         std::cerr << "Duplicate restricted optical reference row: " << line << "\n";
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
   if (optical_rows != 136) {
      std::cerr << "Expected 136 restricted optical reference rows in "
                << optical_filename << ", found " << optical_rows << "\n";
      return EXIT_FAILURE;
   }

   const char *clean_optical_filename = argc == 5
      ? argv[4] : "clean_optical_references.dat";
   std::ifstream clean_optical_input(clean_optical_filename);
   if (!clean_optical_input.is_open()) {
      std::cerr << "Cannot open " << clean_optical_filename << "\n";
      return EXIT_FAILURE;
   }
   int clean_optical_rows = 0;
   std::set<std::tuple<int, double, double, double, double> > clean_optical_keys;
   while (std::getline(clean_optical_input, line)) {
      if (line.empty() || line[0] == '#')
         continue;

      int kernel;
      double x1;
      double y1;
      double x2;
      double y2;
      double expected;
      double relative_tolerance;
      double absolute_tolerance;
      std::string extra;
      std::istringstream values(line);
      if (!(values >> kernel >> x1 >> y1 >> x2 >> y2 >> expected
                   >> relative_tolerance >> absolute_tolerance) || values >> extra) {
         std::cerr << "Invalid clean optical reference row: " << line << "\n";
         return EXIT_FAILURE;
      }
      if (!clean_optical_keys.insert(
            std::make_tuple(kernel, x1, y1, x2, y2)).second) {
         std::cerr << "Duplicate clean optical reference row: " << line << "\n";
         return EXIT_FAILURE;
      }
      ++clean_optical_rows;

      const std::complex<double> z1(x1, y1);
      const std::complex<double> z2(x2, y2);
      const double value = K_builtin(kernel, z1, z2);
      const double tolerance = std::max(absolute_tolerance,
         relative_tolerance*std::abs(expected));
      if (std::isfinite(value) && std::abs(value - expected) <= tolerance) {
         ++passed;
      } else {
         ++failed;
         std::cerr << std::setprecision(17)
                   << "FAILED clean optical: m=" << kernel
                   << " z1=(" << x1 << "," << y1 << ") z2=("
                   << x2 << "," << y2 << ") expected=" << expected
                   << " value=" << value << " tolerance=" << tolerance << "\n";
      }

      const double exchanged = K_builtin(kernel, z2, z1);
      if (std::abs(exchanged - value) <= tolerance) {
         ++passed;
      } else {
         ++failed;
         std::cerr << std::setprecision(17)
                   << "FAILED clean optical exchange symmetry: m=" << kernel
                   << " forward=" << value << " exchanged=" << exchanged << "\n";
      }

      const double mirrored = K_builtin(kernel,
         std::complex<double>(-x1, y1), std::complex<double>(-x2, y2));
      const double parity_value = kernel % 2 == 0 ? -value : value;
      if (std::abs(mirrored - parity_value) <= tolerance) {
         ++passed;
      } else {
         ++failed;
         std::cerr << std::setprecision(17)
                   << "FAILED clean optical parity: m=" << kernel
                   << " value=" << value << " mirrored=" << mirrored << "\n";
      }

      if (kernel % 2 == 1 && y1*y2 > 0.0) {
         if (value >= 0.0) {
            ++passed;
         } else {
            ++failed;
            std::cerr << "FAILED clean optical positivity: m=" << kernel
                      << " value=" << value << "\n";
         }
      }
   }
   if (clean_optical_rows != 58) {
      std::cerr << "Expected 58 clean optical reference rows in "
                << clean_optical_filename << ", found " << clean_optical_rows << "\n";
      return EXIT_FAILURE;
   }

   for (int kernel = 1; kernel <= 8; ++kernel) {
      const std::complex<double> z(0.37, 0.2);
      const double optical = K_builtin(kernel, z, z);
      const double dc = J_builtin(kernel, 2, z);
      if (optical == dc) {
         ++passed;
      } else {
         ++failed;
         std::cerr << std::setprecision(17)
                   << "FAILED coincident clean optical limit: m=" << kernel
                   << " optical=" << optical << " dc=" << dc << "\n";
      }
   }

   const bool saved_legacy_extrapolation = allow_legacy_self_energy_extrapolation;
   const double saved_omega_min = omega_min;
   const double saved_omega_max = omega_max;
   const double saved_optical_frequency = optical_frequency;
   const double saved_frequency_window = spectral_frequency_window_half_width;
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

   spectral_frequency_window_half_width = 0.4;
   transition_points.clear();
   add_spectral_frequency_transition_points(transition_points, -1.0, 1.0, true);
   std::sort(transition_points.begin(), transition_points.end());
   const double shifted_frequency_lower = static_cast<double>(
      -static_cast<long double>(spectral_frequency_window_half_width)
      - optical_frequency);
   const double shifted_frequency_upper = static_cast<double>(
      static_cast<long double>(spectral_frequency_window_half_width)
      - optical_frequency);
   const std::vector<double> expected_frequency_transitions = {
      shifted_frequency_lower, -spectral_frequency_window_half_width,
      shifted_frequency_upper, spectral_frequency_window_half_width
   };
   if (transition_points == expected_frequency_transitions) {
      ++passed;
   } else {
      ++failed;
      std::cerr << "FAILED optical spectral-frequency transition points\n";
   }

   transition_points.clear();
   add_spectral_frequency_transition_points(transition_points, -1.0, 1.0, false);
   std::sort(transition_points.begin(), transition_points.end());
   const std::vector<double> expected_frequency_unshifted = {
      -spectral_frequency_window_half_width, spectral_frequency_window_half_width
   };
   if (transition_points == expected_frequency_unshifted) {
      ++passed;
   } else {
      ++failed;
      std::cerr << "FAILED DC spectral-frequency transition points\n";
   }

   spectral_frequency_window_half_width = 0.0;
   transition_points.clear();
   add_spectral_frequency_transition_points(transition_points, -1.0, 1.0, true);
   if (transition_points.empty()) {
      ++passed;
   } else {
      ++failed;
      std::cerr << "FAILED unrestricted spectral-frequency transition points\n";
   }

   allow_legacy_self_energy_extrapolation = saved_legacy_extrapolation;
   omega_min = saved_omega_min;
   omega_max = saved_omega_max;
   optical_frequency = saved_optical_frequency;
   spectral_frequency_window_half_width = saved_frequency_window;

   std::cout << "CLEAN KERNEL OK=" << passed << " FAILED=" << failed << "\n";
   return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
