#include "dsp.hpp"
#include "globals.hpp"

#include "pffft.hpp"

#include <algorithm>

void DSP::computePS(const std::vector<float> &input,
                    std::vector<std::vector<float>> &frameOutputs,
                    std::vector<float> &psOutput, int frameSize) {
  size_t hop = frameSize / 4; // 75% overlap
  std::vector<float> window(frameSize);

  // Create Hann window
  for (size_t n = 0; n < frameSize; ++n)
    window[n] = 0.5f * (1 - cos(2 * M_PI * n / (frameSize - 1)));

  // Create PFFFT object
  pffft::detail::PFFFT_Setup *pffft =
      pffft_new_setup(frameSize, pffft::detail::PFFFT_REAL);

  if (!pffft) {
    throw std::runtime_error("Failed to initialize PFFFT");
  }

  // Temporary buffers
  std::vector<float> fftIn(frameSize);
  std::vector<float> fftOut(
      frameSize); // PFFFT real output in interleaved format

  // Prepare total power spectrum
  psOutput.assign(frameSize / 2 + 1, 0.0f);
  size_t numFrames = 0;

  // Process overlapping frames
  for (size_t start = 0; start + frameSize <= input.size(); start += hop) {
    // Copy frame and apply window
    for (size_t n = 0; n < frameSize; ++n)
      fftIn[n] = input[start + n] * window[n];

    // Perform real FFT
    pffft_transform(pffft, fftIn.data(), fftOut.data(), nullptr,
                    pffft::detail::PFFFT_FORWARD);

    // Compute magnitude spectrum (first N/2+1 bins)
    std::vector<float> mag(frameSize / 2 + 1);
    mag[0] = std::abs(fftOut[0]); // DC
    for (size_t k = 1; k < frameSize / 2; ++k) {
      float re = fftOut[2 * k];
      float im = fftOut[2 * k + 1];
      mag[k] = std::sqrt(re * re + im * im);
    }
    mag[frameSize / 2] = std::abs(fftOut[1]); // Nyquist

    // Store magnitude spectrum for this frame
    frameOutputs.push_back(mag);

    // Accumulate for total PS (sum of magnitudes squared)
    for (size_t k = 0; k < mag.size(); ++k) {
      psOutput[k] += mag[k] * mag[k];
    }

    numFrames++;
  }

  // Take average or sqrt to get amplitude spectrum per bin
  for (size_t k = 0; k < psOutput.size(); ++k) {
    psOutput[k] = std::sqrt(psOutput[k] / numFrames);
  }

  // Clean up
  pffft_destroy_setup(pffft);
}

void DSP::computeSignificantFreqs(const std::vector<float> &input,
                                  std::vector<float> &output, int sampleRate,
                                  int maxCount) {
  output.clear();
  if (input.size() < 3 || maxCount <= 0)
    return;

  struct Peak {
    float value;
    size_t index;
  };

  std::vector<Peak> peaks;
  peaks.reserve(input.size());

  for (size_t n = 1; n + 1 < input.size(); ++n) {
    if (input[n] > input[n - 1] && input[n] > input[n + 1]) {
      peaks.push_back({input[n], n});
    }
  }

  std::sort(peaks.begin(), peaks.end(),
            [](const Peak &a, const Peak &b) { return a.value > b.value; });

  size_t fftSize = (input.size() - 1) * 2;
  size_t count = std::min(static_cast<size_t>(maxCount), peaks.size());
  output.reserve(count);
  for (size_t i = 0; i < count; ++i) {
    float freqHz =
        static_cast<float>(peaks[i].index) * sampleRate / fftSize;
    output.push_back(freqHz);
  }
}

// input: power spectrum (length = N/2+1)
// output: Mel spectral envelope (length = nMelBands)
void DSP::computeSpectralEnv(const std::vector<float> &input,
                             std::vector<float> &output, int sampleRate,
                             int nMelBands) {
  int nSpec = input.size();
  output.assign(nMelBands, 0.0f);

  // 1. Frequency range
  float fMin = 0.0f;
  float fMax = sampleRate / 2.0f;

  // 2. Convert fMin/fMax to Mel scale
  auto hzToMel = [](float f) { return 2595.0f * log10f(1.0f + f / 700.0f); };
  auto melToHz = [](float m) {
    return 700.0f * (powf(10.0f, m / 2595.0f) - 1.0f);
  };

  float melMin = hzToMel(fMin);
  float melMax = hzToMel(fMax);

  // 3. Compute Mel bin edges
  std::vector<float> melEdges(nMelBands + 2);
  for (int i = 0; i < melEdges.size(); ++i)
    melEdges[i] = melToHz(melMin + i * (melMax - melMin) / (nMelBands + 1));

  // 4. Convert Mel edges to FFT bin indices
  std::vector<int> binEdges(nMelBands + 2);
  for (int i = 0; i < binEdges.size(); ++i)
    binEdges[i] = std::min(
        nSpec - 1, static_cast<int>(roundf(melEdges[i] / fMax * (nSpec - 1))));

  // 5. Compute Mel band energies (triangular filter)
  for (int m = 0; m < nMelBands; ++m) {
    int left = binEdges[m];
    int center = binEdges[m + 1];
    int right = binEdges[m + 2];

    float energy = 0.0f;

    // Rising slope
    for (int k = left; k < center; ++k)
      energy += input[k] * (k - left) / float(center - left);

    // Falling slope
    for (int k = center; k < right; ++k)
      energy += input[k] * (right - k) / float(right - center);

    output[m] = energy;
  }

  // Optional: log scale
  for (auto &v : output)
    v = logf(v + 1e-8f); // avoid log(0)
}
