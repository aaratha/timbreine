#pragma once

#include <complex>
#include <vector>

namespace DSP {

void window(std::vector<float> &input, int size);

// Compute FFT of input samples (re + im)
void computeFFT(const std::vector<float> &input,
                std::vector<std::complex<float>> &output, int size);

void computeIFFT(const std::vector<std::complex<float>> &input,
                 std::vector<float> &output, int size);

void computeRandomPhase(float *phaseOutput, int size);

// Use STFT to compute power spectrum
void computePS(const std::vector<float> &input,
               std::vector<std::vector<float>> &frameOutputs,
               std::vector<float> &psOutput, int frameSize);

// Compute peak frequencies from a spectrum (Hz output)
void computePeakFrequenciesHz(const std::vector<float> &input,
                              std::vector<float> &output, int sampleRate,
                              int maxCount = 10);

// Compute spectral envelope (mel-scaled)
void computeSpectralEnv(const std::vector<float> &input,
                        std::vector<float> &output, int sampleRate,
                        int nMelBands = 40);

// Compute MFCCs from power spectrum
void computeMFCC(const float *input, std::vector<float> &output, int size);

// Compute centroid from input sample
void computeCentroid(const float *input, int size, float &centroid);

// Compute flux from input sample
void computeFlux(const float *input, const float *prevInput, int size,
                 float &flux);

// Compute rolloff from input sample
void computeRolloff(const float *input, int size, float &rolloff);

// Compute UMAP from MFCCs, centroid, flux, and rolloff
void computeUMAP(const std::vector<std::vector<float>> &input,
                 std::vector<std::vector<float>> &output, int nComponents);
}; // namespace DSP
