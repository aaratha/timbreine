#pragma once

#include <complex>
#include <vector>

struct VoronoiEdge {
  float x1;
  float y1;
  float x2;
  float y2;
  int site;
};

struct VoronoiCell {
  int site;
  std::vector<float> coords;
};

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
  // Final UMAP timbre coordinates
  std::vector<float> binTimbreX;
  std::vector<float> binTimbreY;
  // Voronoi diagram data for timbre map
  std::vector<VoronoiEdge> voronoiEdges;
  std::vector<VoronoiCell> voronoiCells;

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

  // compute UMAP coordinates from features
  void computeUmapCoordinates();

  // compute Voronoi diagram edges from UMAP coordinates
  void computeVoronoiEdges();

  size_t getBinIndex();
  std::vector<float> &getInputRaw();
  size_t getBinCount() const;
  const std::vector<float> &getBin(size_t index) const;
  const std::vector<float> &getBinTimbreX() const;
  const std::vector<float> &getBinTimbreY() const;
  const std::vector<float> &getVoronoiCoords() const;
  const std::vector<VoronoiEdge> &getVoronoiEdges() const;
  const std::vector<VoronoiCell> &getVoronoiCells() const;

  void setBinIndex(size_t index);
};
