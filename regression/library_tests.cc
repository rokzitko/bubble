#include "bubble.hh"

#include <gsl/gsl_errno.h>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

using namespace bubble;

namespace {
int checks = 0;
void check(bool condition, const std::string &message) {
   ++checks;
   if (!condition) throw std::runtime_error(message);
}
void close(double actual, double expected, double relative = 2e-10, double absolute = 1e-12) {
   std::ostringstream message;
   message << std::setprecision(17) << "actual=" << actual << " expected=" << expected;
   check(std::isfinite(actual) && std::abs(actual - expected) <=
         absolute + relative*std::abs(expected), message.str());
}
void handler(const char *, const char *, int, int) {}
void handler_restored() {
   gsl_error_handler_t *previous = gsl_set_error_handler(&handler);
   check(previous == &handler, "GSL handler was not restored");
}
void fails(const std::function<void()> &call, const std::string &text) {
   bool caught = false;
   try { call(); }
   catch (const std::runtime_error &error) {
      caught = std::string(error.what()).find(text) != std::string::npos;
      check(caught, std::string("Wrong exception: ") + error.what());
   }
   check(caught, "Expected exception containing " + text);
   handler_restored();
}
SelfEnergy load_sigma(const std::string &root, const std::string &re, const std::string &im) {
   std::ifstream real(root + "/" + re), imaginary(root + "/" + im);
   check(real.good() && imaginary.good(), "Cannot read Sigma fixtures");
   SelfEnergy result;
   double w, r, wi, i;
   while (real >> w >> r) {
      check(bool(imaginary >> wi >> i) && w == wi, "Invalid Sigma fixture");
      result.omega.push_back(w); result.real.push_back(r); result.imaginary.push_back(i);
   }
   check(real.eof() && !(imaginary >> wi), "Mismatched Sigma fixtures");
   return result;
}
std::vector<double> values(Context &context) {
   std::vector<double> result;
   for (int moment = 0; moment <= 2; ++moment) result.push_back(context.dc(5, moment).value);
   result.push_back(context.optical(0.3).value);
   result.push_back(context.kinetic_energy().value);
   const auto spectral = context.spectral({-0.3, 0.0, 0.4});
   result.insert(result.end(), spectral.begin(), spectral.end());
   handler_restored();
   return result;
}
void emit(const std::string &name, double value) {
   std::cout << "VALUE " << name << " " << std::setprecision(17) << value << "\n";
}
}

int main(int argc, char **argv) {
   const std::string root = argc > 1 ? argv[1] : ".";
   gsl_error_handler_t *original = gsl_set_error_handler(&handler);
   try {
      Config a_config;
      a_config.temperature = 0.1;
      SelfEnergy a_sigma{{-10, -5, 0, 5, 10}, {0, 0, 0, 0, 0}, {-0.1, -0.1, -0.1, -0.1, -0.1}};
      Config b_config = a_config;
      b_config.temperature = 0.07;
      b_config.chemical_potential = 0.23;
      b_config.interpolation = Interpolation::steffen;
      b_config.integration_rule = 6;
      b_config.workspace_limit = 2000;
      SelfEnergy b_sigma = a_sigma;
      b_sigma.real = {-0.3, -0.2, 0.12, 0.2, 0.3};
      b_sigma.imaginary = {-0.4, -0.3, -0.2, -0.3, -0.4};
      Context a(a_config, a_sigma), b(b_config, b_sigma);
      handler_restored();
      const auto av = values(a);
      const auto bv = values(b);
      check(av != bv, "A and B fixtures must differ");
      check(values(a) == av, "A/B/A changed A results");
      check(values(b) == bv, "A/B/A/B changed B results");
      // Input storage is not retained by reference.
      a_sigma.real.assign(5, 999.0);
      a_config.chemical_potential = 42.0;
      check(values(a) == av, "Context did not own inputs/configuration");
      a_sigma.real.assign(5, 0.0);
      a_config.chemical_potential = 0.0;
      for (int repeat = 0; repeat < 3; ++repeat) {
         Context fresh_a(a_config, a_sigma), fresh_b(b_config, b_sigma);
         check(values(fresh_a) == av && values(fresh_b) == bv,
               "Destruction/reconstruction leaked dataset state");
      }
      for (int moment = 0; moment <= 2; ++moment) {
         emit("dc" + std::to_string(moment), av[moment]);
         close(a.dc(5, moment, 2, true).value, av[moment], 1e-8);
      }
      emit("optical", av[3]);
      emit("kinetic", av[4]);
      check(a.optical(0.0).value == av[0], "Omega=0 is not exact DC");
      close(a.optical(1e-8).value, av[0], 1e-9);

      Config tight = a_config;
      tight.cutoff = 15; tight.absolute_error = tight.relative_error = 1e-12;
      Context optical(tight, a_sigma);
      const double references[] = {0.4885903945334426, 0.4708182948628991,
                                   0.4420559614703037, 0.4343363161009058};
      for (int k = 0; k < 4; ++k)
         close(optical.optical(0.3, 2*k+1).value, references[k], 1e-10, 5e-12);

      // Identical Phi domain, different values: catches the old static mass cache.
      PhiTable phi;
      phi.epsilon = {-1, -0.5, 0, 0.5, 1};
      phi.value = {1, 1, 1, 1, 1};
      PhiTable phi_b = phi; phi_b.value.assign(5, 3.0);
      Context tab_a(a_config, a_sigma, phi), tab_b(b_config, b_sigma, phi_b);
      close(tab_a.spectral({0}, 0, 0)[0], 2.0, 0, 0);
      close(tab_b.spectral({0}, 0, 0)[0], 6.0, 0, 0);
      close(tab_a.spectral({0}, 0, 0)[0], 2.0, 0, 0);
      const double tab_dc = tab_a.dc(0).value;
      emit("tab_dc", tab_dc);
      emit("tab_optical", tab_a.optical(0.3, 0).value);
      tab_b.optical(0.3, 0);
      check(tab_a.dc(0).value == tab_dc, "Tabulated A/B/A changed DC");
      close(tab_dc, a.dc(1).value, 1e-8);
      close(tab_a.optical(0.3, 0).value, a.optical(0.3, 1).value, 1e-8);
      const auto spectral = tab_a.spectral({-0.3, 0.0, 0.4}, 0);
      for (std::size_t i = 0; i < spectral.size(); ++i)
         close(spectral[i], a.spectral({-0.3, 0.0, 0.4}, 1)[i]);
      emit("tab_spectral_zero", spectral[1]);
      phi.epsilon_power = 1;
      Context energy(a_config, a_sigma, phi);
      close(energy.kinetic_energy(0).value, a.kinetic_energy(2).value, 1e-8);
      emit("tab_kinetic", energy.kinetic_energy(0).value);
      for (int method = 0; method <= 4; ++method) {
         Config config = a_config;
         config.phi_interpolation = static_cast<PhiInterpolation>(method);
         config.interpolation = static_cast<Interpolation>(method ? method : 1);
         Context context(config, a_sigma, phi);
         close(context.kinetic_energy(0).value, a.kinetic_energy(2).value, 1e-8);
      }
      Config clip_a = a_config, clip_b = a_config;
      clip_a.sigma_clip_floor = 1e-3; clip_b.sigma_clip_floor = 0.2;
      SelfEnergy positive = a_sigma; positive.imaginary.assign(5, 0.1);
      Context clipped_a(clip_a, positive), clipped_b(clip_b, positive);
      const double clipped = clipped_a.spectral({0})[0];
      close(clipped, std::sqrt(1+1e-6)-1e-3);
      close(clipped_b.spectral({0})[0], std::sqrt(1+0.04)-0.2);
      check(clipped_a.spectral({0})[0] == clipped, "Clipping floor leaked across datasets");

      // Independent hidden-band and remote-knot references from run_dc_mesh_tests.
      Config mesh;
      mesh.temperature = 1; mesh.sigma_clip_floor = 1e-12; mesh.absolute_error = 1e-14;
      SelfEnergy affine{{-30, 30}, {-30-16384*(-30-2.75), 30-16384*(30-2.75)}, {-0.125, -0.125}};
      const double weight = std::exp(-2.75)/std::pow(1+std::exp(-2.75), 2);
      for (int rule : {1, 6}) {
         mesh.integration_rule = rule;
         Context hidden(mesh, affine);
         for (int moment = 0; moment <= 2; ++moment)
            close(hidden.dc(5, moment, 2, true).value,
                  3/(16*0.125*16384)*weight*std::pow(2.75, moment), 1e-8, 0);
      }
      mesh.temperature = 1e8; mesh.chemical_potential = -2.75; mesh.absolute_error = 1e-20;
      Context remote(mesh, SelfEnergy{{-1e100, 1e100}, {0, 0}, {-0.125, -0.125}});
      close(remote.dc(5, 0, 2, true).value, 3.75e-9, 1e-8, 0);
      Config hot;
      hot.temperature = 1e308; hot.cutoff = hot.sigma_clip_floor = 1e-300;
      Context thermal_overflow(hot, SelfEnergy{{-1e8, 1e8}, {-1e8, 1e8}, {-1e-300, -1e-300}});
      fails([&] { thermal_overflow.dc(5, 0, 2, true); }, "unsupported thermal denominator");

      const auto dmft = load_sigma(root, "Bethe_lattice_test/resigma.dat", "Bethe_lattice_test/imsigma.dat");
      Config bethe;
      bethe.temperature = 0.27; bethe.interpolation = Interpolation::akima; bethe.integration_rule = 6;
      Context dmft_context(bethe, dmft);
      const double pi = std::acos(-1.0);
      const double d0 = dmft_context.dc().value;
      const double d1 = dmft_context.dc(5, 1).value;
      const double d2 = dmft_context.dc(5, 2).value;
      close(2*pi*d0, 9.777905785035401, 2e-6, 1e-10);
      close(-d1/(bethe.temperature*d0), -0.282541905816773, 2e-6, 1e-10);
      close(d2/(bethe.temperature*bethe.temperature*d0), 2.34072492836283, 2e-6, 1e-10);
      close(2/pi*dmft_context.kinetic_energy().value, -0.14512837773776, 2e-6, 1e-10);
      emit("dmft_dc", d0);
      emit("dmft_optical", dmft_context.optical(0.1).value);
      emit("dmft_kinetic", dmft_context.kinetic_energy().value);
      check(values(a) == av, "DMFT evaluation changed A");

      // Validation failures never terminate the host and never poison later calls.
      fails([&] { Context bad(a_config, SelfEnergy{}); }, "empty");
      fails([&] { Context bad(a_config, SelfEnergy{{0}, {0}, {-0.1}}); }, "at least 2 rows");
      PhiTable invalid_phi = phi; invalid_phi.value.pop_back();
      fails([&] { Context bad(a_config, a_sigma, invalid_phi); }, "mismatched");
      invalid_phi = phi; invalid_phi.epsilon[2] = invalid_phi.epsilon[1];
      fails([&] { Context bad(a_config, a_sigma, invalid_phi); }, "strictly increasing");
      SelfEnergy invalid = a_sigma; invalid.imaginary.pop_back();
      fails([&] { Context bad(a_config, invalid); }, "mismatched");
      invalid = a_sigma; invalid.omega[2] = invalid.omega[1];
      fails([&] { Context bad(a_config, invalid); }, "strictly increasing");
      invalid = a_sigma; invalid.real[2] = std::numeric_limits<double>::quiet_NaN();
      fails([&] { Context bad(a_config, invalid); }, "non-finite");
      Config bad_config = a_config; bad_config.workspace_limit = 0;
      fails([&] { Context bad(bad_config, a_sigma); }, "workspace");
      bad_config = a_config; bad_config.interpolation = static_cast<Interpolation>(9);
      fails([&] { Context bad(bad_config, a_sigma); }, "selector");
      fails([&] { a.dc(0); }, "Phi table");
      fails([&] { a.dc(9); }, "kernel");
      fails([&] { a.dc(5, -1); }, "Frequency power");
      fails([&] { a.optical(-1); }, "nonnegative");
      fails([&] { a.optical(8); }, "closed interval");
      fails([&] { a.spectral({11}); }, "outside");
      fails([&] { b.dc(5, 0, 2, true); }, "--dc-mesh");
      check(values(a) == av && values(b) == bv, "Failure poisoned context");
      check(a.diagnostics().empty(), "Failed call diagnostics leaked into successful call");
      Config spectral_only;
      Context zero_temperature(spectral_only, a_sigma);
      close(zero_temperature.spectral({0})[0], a.spectral({0})[0], 0, 0);
      fails([&] { zero_temperature.dc(); }, "positive finite T");

      // Workspace/error-policy state must not survive A/B/A or suppress warnings.
      Config limited = a_config; limited.workspace_limit = 1;
      Context strict(limited, a_sigma, phi);
      fails([&] { strict.optical(0.3, 1); }, "iterations");
      fails([&] { strict.dc(0); }, "iterations"); // Nested callback failure.
      close(strict.spectral({0}, 0, 0)[0], 0.0, 0, 1e-15); // Signed Phi mass.
      limited.error_policy = ErrorPolicy::warn;
      Context warn(limited, a_sigma), warn_b(limited, b_sigma);
      const auto warned = warn.optical(0.3, 1);
      check(!warned.diagnostics.empty(), "Missing integration warning");
      warn_b.optical(0.3, 1);
      check(warn.optical(0.3, 1).diagnostics == warned.diagnostics, "Warning leaked across calls");
      limited.error_policy = ErrorPolicy::ignore;
      Context ignore(limited, a_sigma);
      const auto ignored = ignore.optical(0.3, 1);
      check(ignored.value == warned.value && ignored.diagnostics.empty(), "Ignore policy changed result");
      check(values(a) == av, "Workspace/error-policy state leaked to A");
      handler_restored();
      gsl_set_error_handler(original);
      std::cout << "LIBRARY OK=" << checks << " FAILED=0\n";
      return 0;
   } catch (const std::exception &error) {
      gsl_set_error_handler(original);
      std::cerr << "LIBRARY FAILED after " << checks << " checks: " << error.what() << "\n";
      return 1;
   }
}
