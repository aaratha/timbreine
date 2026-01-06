#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "dsp.hpp"

TEST_CASE("computeRolloff finds threshold bin") {
  const int fftSize = 8;
  const int sampleRate = 8000;
  const int nSpec = fftSize / 2 + 1;

  std::vector<float> powerSpec(nSpec, 1.0f);

  float rolloff = 0.0f;
  DSP::computeRolloff(powerSpec, fftSize, sampleRate, 0.6f, rolloff);

  const float expected = 2.0f * sampleRate / fftSize;
  REQUIRE(rolloff == Catch::Approx(expected));
}
