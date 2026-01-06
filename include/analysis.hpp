#pragma once

#include <vector>
#include <complex>

class AnalysisCore {
  std::vector<float> inputRaw;
  size_t inputBinSize;
  // store binned input samples and later windowed bins
  std::vector<std::vector<float>> bins;
  // store fft for each bin
  std::vector<std::vector<std::complex<float>>> binFFTs;
  // per-bin: mfccs (12), centroid, flux, rolloff
  std::vector<std::vector<float>> umapFeatures;
  // store per-bin phase accumulators for deterministic resynthesis
  std::vector<std::vector<float>> binPhaseAccumulators;

  size_t binIndex{0};

public:
  AnalysisCore();
  ~AnalysisCore();

  void readFile(const std::string &filename);

  // Bin input file into 4096 sample bins with 75% overlap
  void binInput();

  // Decompose each bin into frequency components (re + im) 
  void decomposeBins();

  void resynthesizeBin(size_t binIndex, std::vector<float> &output);

  // compute and store MFCCs, centroid, flux, rolloff for each bin
  void computeUmapFeatures();

  size_t getBinIndex();
  std::vector<float> &getInputRaw();
  size_t getBinCount() const;
  const std::vector<float> &getBin(size_t index) const;

  void setBinIndex(size_t index);
};
