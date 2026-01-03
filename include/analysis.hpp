#pragma once

#include <vector>

class AnalysisCore {
  std::vector<float> inputRaw;
  size_t inputBinSize;
  // store binned input samples and later windowed bins
  std::vector<std::vector<float>> binnedInput;
  // store PS for each frame (STFT) to use for flux
  std::vector<std::vector<float>> framePS;
  // store overall PS for each bin
  std::vector<std::vector<float>> binPS;
  // store significant frequency components for each bin
  std::vector<std::vector<float>> binFreqComponents;
  // store spectral envelopes for each bin (mel-scaled)
  std::vector<std::vector<float>> binSpectralEnv;

public:
  AnalysisCore();
  ~AnalysisCore();

  void readFile(const std::string &filename);

  // Bin input file into 24 equal time segments
  // Input bin size is determined by total samples / 24
  // Includes 1/24th of samples, + overlapping 75% (25% in history + look ahead = hop)
  void binInput();

  // Apply han window to each bin
  void windowBins();

  // Decompose each bin into frequency components (power spectrum) 
  // Splits bins into subbins for STFT analysis
  void decomposeBins();

  // Find significant frequency components and spectral envelopes for synthesis
  void findSynthesisFeatures();

  std::vector<float> &getInputRaw();
  std::vector<std::vector<float>> &getBinFreqComponents();
  std::vector<std::vector<float>> &getBinSpectralEnvs();
};
