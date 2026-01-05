#pragma once

#include <vector>
#include <complex>

struct BinFeatures {
  std::vector<float> peakFrequenciesHz;
  std::vector<float> spectralEnvelope;
};

class AnalysisCore {
  std::vector<float> inputRaw;
  size_t inputBinSize;
  // store binned input samples and later windowed bins
  std::vector<std::vector<float>> bins;
  // store fft for each bin
  std::vector<std::vector<std::complex<float>>> binFFTs;
  // store spectral features per bin
  std::vector<BinFeatures> binFeatures;

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

  // Find significant frequency components and spectral envelopes for synthesis
  void findSynthesisFeatures();

  size_t getBinIndex();
  std::vector<float> &getInputRaw();
  const std::vector<BinFeatures> &getBinFeatures() const;
  size_t getBinCount() const;
  const std::vector<float> &getBin(size_t index) const;

  void setBinIndex(size_t index);
};
