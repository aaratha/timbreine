#pragma once

#include <complex>
#include <vector>

// TODO: Refactor to take globals as parameters instead of using macros

namespace DSP {

void window(std::vector<float> &input, int size);

// Compute FFT of input samples (re + im)
void computeFFT(const std::vector<float> &input,
                std::vector<std::complex<float>> &output, int size);

void computeIFFT(const std::vector<std::complex<float>> &input,
                 std::vector<float> &output, int size);

// Compute power spectrum from FFT input
void computePowerSpectrum(const std::vector<std::complex<float>> &input,
                          std::vector<float> &output, int size);

// input: Power spectrum (length = nSpec)
// output: Mel-scale spectral envelope (length = nMelBands) - NOT LOG-SCALED YET
void computeSpectralEnv(const std::vector<float> &powerSpec,
                        std::vector<float> &melEnv, int sampleRate, int fftSize,
                        int nMelBands);

// Compute MFCCs from power spectrum
// input: Mel-scale spectral envelope (length = nMelBands)
void computeMFCC(const std::vector<float> &input, std::vector<float> &output,
                 int size);

// Compute centroid from input power spectrum
void computeCentroid(const std::vector<float> &input, int size,
                     float &centroid);

// Compute flux from current and previous power spectra
void computeFlux(const std::vector<float> &psCurr,
                 const std::vector<float> &psPrev, int size, float &flux);

// Compute rolloff from input power spectrum
void computeRolloff(const std::vector<float> &ps, int size, float &rolloff);

// Compute UMAP from MFCCs, centroid, flux, and rolloff
void computeUMAP(const std::vector<std::vector<float>> &input,
                 std::vector<std::vector<float>> &output, int nComponents);
}; // namespace DSP
