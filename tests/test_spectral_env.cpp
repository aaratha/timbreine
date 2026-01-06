#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>

#include "dsp.hpp"

TEST_CASE("computeSpectralEnv produces nonnegative mel bands") {
  const int fftSize = 8;
  const int sampleRate = 8000;
  const int nSpec = fftSize / 2 + 1;
  const int nMelBands = 4;

  std::vector<float> powerSpec(nSpec, 1.0f);
  std::vector<float> melEnv;

  DSP::computeSpectralEnv(powerSpec, melEnv, sampleRate, fftSize, nMelBands);

  REQUIRE(melEnv.size() == static_cast<size_t>(nMelBands));

  float sum = 0.0f;
  for (float value : melEnv) {
    REQUIRE(std::isfinite(value));
    REQUIRE(value >= 0.0f);
    sum += value;
  }
  REQUIRE(sum > 0.0f);
}
