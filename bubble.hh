#ifndef BUBBLE_HH
#define BUBBLE_HH

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace bubble {

enum class Interpolation { linear = 1, cspline, akima, steffen };
enum class PhiInterpolation { automatic, linear, cspline, akima, steffen };
enum class ErrorPolicy { ignore, warn, fail };

struct Config {
   double temperature = 0.0; // Positive for frequency integrals; unused by spectral().
   double chemical_potential = 0.0;
   Interpolation interpolation = Interpolation::linear;
   PhiInterpolation phi_interpolation = PhiInterpolation::automatic;
   ErrorPolicy error_policy = ErrorPolicy::fail;
   int integration_rule = 1;
   double absolute_error = 1e-7;
   double relative_error = 1e-8;
   std::size_t workspace_limit = 1000;
   double cutoff = 30.0;
   double sigma_clip_floor = 1e-8;
};

struct SelfEnergy {
   std::vector<double> omega, real, imaginary;
};

struct PhiTable {
   std::vector<double> epsilon, value;
   int epsilon_power = 0; // Applied to knots before interpolation, as with CLI -e.
};

struct Result {
   double value;
   double outer_error; // Not an end-to-end error bound; see README.
   std::string diagnostics; // Warnings, never written to stdout/stderr by the API.
};

// Owns the configuration, interpolants and numerical workspaces. Inputs are copied.
// Invalid input/numerical failures throw std::runtime_error (or std::bad_alloc).
// Calls are serialized internally because GSL's error handler is process-global.
// Do not use other GSL code concurrently with these calls. No callbacks into users.
class Context {
public:
   explicit Context(const Config &, const SelfEnergy &, const PhiTable & = PhiTable{});
   ~Context();
   Context(const Context &) = delete;
   Context &operator=(const Context &) = delete;

   // CLI m,n,o conventions and raw normalization, including for optical/energy.
   Result dc(int kernel = 5, int frequency_power = 0, int spectral_power = 2,
             bool mesh = false);
   Result optical(double frequency, int kernel = 5, int frequency_power = 0);
   Result occupied(int kernel, int spectral_power = 1, int frequency_power = 0);
   // Equivalent to occupied(kernel,1,0). Use m=4 for the unnormalized Bethe
   // energy kernel, or m=0 with a DOS PhiTable and epsilon_power=1.
   Result kinetic_energy(int kernel = 4);
   // J_mn at arbitrary supplied frequencies; n=1 is the lattice spectral value
   // when the selected kernel is a DOS. No implicit DOS normalization.
   std::vector<double> spectral(const std::vector<double> &frequencies,
                                int kernel = 3, int spectral_power = 1);
   std::string diagnostics() const; // Most recent call (also after a failure).

private:
   struct Impl;
   std::unique_ptr<Impl> impl_;
};

} // namespace bubble
#endif
