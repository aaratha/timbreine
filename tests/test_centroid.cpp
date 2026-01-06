#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "dsp.hpp"

TEST_CASE("computeCentroid uses spectral energy distribution") {
  const int fftSize = 8;
  const int sampleRate = 8000;
  const int nSpec = fftSize / 2 + 1;

  std::vector<float> powerSpec(nSpec, 0.0f);
  powerSpec[2] = 1.0f;

  float centroid = 0.0f;
  DSP::computeCentroid(powerSpec, fftSize, sampleRate, centroid);

  const float expected = 2.0f * sampleRate / fftSize;
  REQUIRE(centroid == Catch::Approx(expected));
}
