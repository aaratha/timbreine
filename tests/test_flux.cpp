#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>

#include "dsp.hpp"

TEST_CASE("computeFlux compares power spectra") {
  const int fftSize = 8;
  const int nSpec = fftSize / 2 + 1;

  std::vector<float> psCurr(nSpec, 1.0f);
  std::vector<float> psPrev(nSpec, 0.0f);

  float flux = 0.0f;
  DSP::computeFlux(psCurr, psPrev, fftSize, flux);

  const float expected = std::sqrt(static_cast<float>(nSpec));
  REQUIRE(flux == Catch::Approx(expected));
}
