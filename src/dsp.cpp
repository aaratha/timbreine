#include "dsp.hpp"
#include "globals.hpp"

#include "pffft.hpp"

#include <algorithm>

void DSP::window(std::vector<float> &input, int size) {
  // Apply Hann window
  for (int n = 0; n < size; ++n) {
    float w = 0.5f * (1 - cos(2 * M_PI * n / (size - 1)));
    input[n] *= w; // Note: input is const, so this function does not modify it
  }
}

void DSP::computeFFT(const std::vector<float> &input,
                     std::vector<std::complex<float>> &output, int size) {

  // Create PFFFT object
  pffft::detail::PFFFT_Setup *pffft =
      pffft_new_setup(size, pffft::detail::PFFFT_REAL);

  if (!pffft) {
    throw std::runtime_error("Failed to initialize PFFFT");
  }

  float *pffftOutput = new float[size]; // PFFFT output in interleaved format

  // Perform real FFT
  pffft_transform(pffft, input.data(), pffftOutput, nullptr,
                  pffft::detail::PFFFT_FORWARD);

  // Convert interleaved floats -> std::complex<float>
  for (size_t k = 0; k < size; k++) {
    float re = pffftOutput[2 * k + 0];
    float im = pffftOutput[2 * k + 1];
    output[k] = std::complex<float>(re, im);
  }

  // Clean up
  pffft_destroy_setup(pffft);
}

void DSP::computeIFFT(const std::vector<std::complex<float>> &input,
                      std::vector<float> &output, int size) {

  // Prepare buffer for interleaved format (Re, Im)
  std::vector<float> fftInterleaved(2 * size);
  for (size_t k = 0; k < size; ++k) {
    fftInterleaved[2 * k + 0] = input[k].real();
    fftInterleaved[2 * k + 1] = input[k].imag();
  }

  // 4. Allocate output buffer
  output.resize(size);

  // Perform inverse FFT
  pffft::detail::PFFFT_Setup *setup =
      pffft_new_setup(size, pffft::detail::PFFFT_REAL);

  pffft_transform_ordered(setup,
                          fftInterleaved.data(), // input interleaved Re+Im
                          output.data(),         // output time-domain
                          nullptr, pffft::detail::PFFFT_BACKWARD);

  pffft_destroy_setup(setup);
}

void DSP::computePowerSpectrum(const std::vector<std::complex<float>> &input,
                               std::vector<float> &output, int size) {
  output.resize(size);
  for (int k = 0; k < size; ++k) {
    float re = input[k].real();
    float im = input[k].imag();
    output[k] = re * re + im * im; // Power = Re^2
  }
}

void DSP::computeSpectralEnv(
    const std::vector<float> &powerSpec, // |FFT|^2, size = fftSize/2 + 1
    std::vector<float> &melEnv, int sampleRate, int fftSize, int nMelBands) {
  const int nSpec = fftSize / 2 + 1;
  melEnv.assign(nMelBands, 0.0f);

  // --- Hz <-> Mel ---
  auto hzToMel = [](float f) {
    return 2595.0f * std::log10(1.0f + f / 700.0f);
  };
  auto melToHz = [](float m) {
    return 700.0f * (std::pow(10.0f, m / 2595.0f) - 1.0f);
  };

  const float fMin = 0.0f;
  const float fMax = sampleRate * 0.5f;

  const float melMin = hzToMel(fMin);
  const float melMax = hzToMel(fMax);

  // --- Mel band edges (Hz) ---
  std::vector<float> melHz(nMelBands + 2);
  for (int i = 0; i < nMelBands + 2; ++i) {
    float mel = melMin + (melMax - melMin) * i / (nMelBands + 1);
    melHz[i] = melToHz(mel);
  }

  // --- Convert Hz -> FFT bin indices ---
  std::vector<int> bins(nMelBands + 2);
  for (int i = 0; i < nMelBands + 2; ++i) {
    int b = static_cast<int>(std::floor(melHz[i] * fftSize / sampleRate));
    bins[i] = std::clamp(b, 0, nSpec - 1);
  }

  // --- Apply triangular mel filters ---
  for (int m = 0; m < nMelBands; ++m) {
    int left = bins[m];
    int center = bins[m + 1];
    int right = bins[m + 2];

    if (center <= left || right <= center)
      continue; // avoid zero-width filters

    float energy = 0.0f;
    float norm = 0.0f;

    // Rising slope
    for (int k = left; k < center; ++k) {
      float w = (k - left) / float(center - left);
      energy += powerSpec[k] * w;
      norm += w;
    }

    // Falling slope
    for (int k = center; k < right; ++k) {
      float w = (right - k) / float(right - center);
      energy += powerSpec[k] * w;
      norm += w;
    }

    // Normalize for equal-area filters
    if (norm > 0.0f)
      energy /= norm;

    melEnv[m] = energy;
  }
}

void DSP::computeMFCC(const std::vector<float> &input,
                      std::vector<float> &output,
                      int numCoeffs) // number of MFCCs to compute
{
  const int M = input.size(); // number of mel bands
  output.resize(numCoeffs);

  // --- Step 1: log compression ---
  constexpr float eps = 1e-10f; // prevent log(0)
  std::vector<float> logMel(M);
  for (int m = 0; m < M; ++m)
    logMel[m] = std::log(input[m] + eps);

  // --- Step 2: precompute cosine table ---
  // Optional: can also store as a class member to reuse every frame
  std::vector<std::vector<float>> cosTable(numCoeffs, std::vector<float>(M));
  for (int n = 0; n < numCoeffs; ++n)
    for (int m = 0; m < M; ++m)
      cosTable[n][m] = std::cos(M_PI * n * (m + 0.5f) / M);

  // --- Step 3: DCT-II ---
  for (int n = 0; n < numCoeffs; ++n) {
    float sum = 0.0f;
    for (int m = 0; m < M; ++m)
      sum += logMel[m] * cosTable[n][m];

    output[n] = sum;
  }

  // --- Step 4: optional orthonormal scaling ---
  // Many MFCC implementations scale like this:
  // Useful for comparisons with other MFCCs
  output[0] *= std::sqrt(1.0f / M);
  for (int n = 1; n < numCoeffs; ++n)
    output[n] *= std::sqrt(2.0f / M);
}

void DSP::computeCentroid(const std::vector<float> &input, int size,
                          float &centroid) {
  const int nSpec = size / 2 + 1;
  float num = 0.0f;
  float denom = 0.0f;
  for (int k = 0; k < nSpec; k++) {
    num += (k * DEVICE_SAMPLE_RATE / size) * input[k];
    denom += input[k];
  }
  if (denom > 0.0f) {
    centroid = num / denom;
  } else {
    centroid = 0.0f;
  }
}

void DSP::computeFlux(const std::vector<float> &psCurr,
                      const std::vector<float> &psPrev, int size, float &flux) {
  const int nSpec = size / 2 + 1;
  flux = 0.0f;
  for (int k = 0; k < nSpec; k++) {
    float diff = psCurr[k] - psPrev[k];
    flux += diff * diff;
  }
  flux = std::sqrt(flux);
}

// Computes spectral roll-off frequency in Hz
// ps: power or magnitude spectrum (length = size/2 + 1)
// size: original FFT size
// sampleRate: audio sample rate
// threshold: fraction of total spectral energy to consider (0.85 = 85%)
// rolloff: output frequency in Hz
void DSP::computeRolloff(const std::vector<float> &ps, int size,
                         float &rolloff) {
  const int nSpec = size / 2 + 1; // unique FFT bins
  float totalEnergy = 0.0f;

  // Compute total spectral energy
  for (int k = 0; k < nSpec; ++k)
    totalEnergy += ps[k];

  // Compute cumulative energy and find bin exceeding threshold
  float cumulative = 0.0f;
  int rollBin = 0;
  float target = ROLLOFF_THRESHOLD * totalEnergy;

  for (int k = 0; k < nSpec; ++k) {
    cumulative += ps[k];
    if (cumulative >= target) {
      rollBin = k;
      break;
    }
  }

  // Convert bin to frequency in Hz
  rolloff = (rollBin * float(DEVICE_SAMPLE_RATE)) / size;
}
