#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "dsp.hpp"

TEST_CASE("computeMFCC handles flat input") {
  const int melBands = 8;
  const int coeffs = 6;

  std::vector<float> melEnv(melBands, 1.0f);
  std::vector<float> mfccs;

  DSP::computeMFCC(melEnv, mfccs, coeffs);

  REQUIRE(mfccs.size() == static_cast<size_t>(coeffs));
  for (float value : mfccs) {
    REQUIRE(value == Catch::Approx(0.0f).margin(1e-4f));
  }
}
