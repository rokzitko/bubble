#include <bubble.hh>

#ifdef BUBBLE_LIBRARY
#error BUBBLE_LIBRARY must not propagate to consumers
#endif

#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
   try {
      bubble::Config config;
      config.temperature = 0.1;
      bubble::SelfEnergy sigma{{-10.0, 10.0}, {0.0, 0.0}, {-0.1, -0.1}};
      bubble::Context context(config, sigma);
      const double dc = context.dc().value;
      const double optical = context.optical(0.3).value;
      const double spectral = context.spectral({0.0})[0];
      if (!std::isfinite(dc) || !std::isfinite(optical) || !std::isfinite(spectral) ||
          std::abs(dc - 1.4962863230896153) > 1e-10 ||
          std::abs(optical - 0.44205596151879128) > 1e-10 ||
          std::abs(spectral - (std::sqrt(1.01) - 0.1)) > 1e-12)
         throw std::runtime_error("Installed Bubble numerical smoke test failed");
      std::cout << "Installed Bubble DC/optical/spectral checks passed\n";
   } catch (const std::exception &error) {
      std::cerr << error.what() << '\n';
      return 1;
   }
}
