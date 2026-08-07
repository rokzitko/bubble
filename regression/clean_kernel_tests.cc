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

double J_builtin(int kernel, int power, std::complex<double> z);

int main(int argc, char **argv)
{
   const char *filename = argc == 2 ? argv[1] : "clean_kernel_references.dat";
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

   std::cout << "CLEAN KERNEL OK=" << passed << " FAILED=" << failed << "\n";
   return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
