#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#include "analysis.hpp"
#include "dsp.hpp"
#include "globals.hpp"
#include "miniaudio.h"

#define JC_VORONOI_IMPLEMENTATION
#include "jc_voronoi.h"

AnalysisCore::AnalysisCore() {}

AnalysisCore::~AnalysisCore() {}

void AnalysisCore::readFile(const std::string &filename) {
  ma_decoder decoder;

  ma_decoder_config config =
      ma_decoder_config_init(ma_format_f32, 2, DEVICE_SAMPLE_RATE);

  // 0 = use file's channels and sample rate
  ma_result result = ma_decoder_init_file(filename.c_str(), &config, &decoder);
  if (result != MA_SUCCESS) {
    std::cerr << "Failed to open file\n";
    return;
  }

  // Get total frame count
  ma_uint64 totalFrameCount = 0;
  result = ma_decoder_get_length_in_pcm_frames(&decoder, &totalFrameCount);
  if (result != MA_SUCCESS) {
    std::cerr << "Failed to read file length\n";
    ma_decoder_uninit(&decoder);
    return;
  }

  // Prepare buffer: channels * frames
  std::vector<float> buffer(totalFrameCount * decoder.outputChannels);

  // Read frames into buffer
  ma_uint64 framesRead = 0;
  result = ma_decoder_read_pcm_frames(&decoder, buffer.data(), totalFrameCount,
                                      &framesRead);
  if (result != MA_SUCCESS || framesRead != totalFrameCount) {
    std::cerr << "Warning: did not read all frames\n";
  }

  inputBinSize = INPUT_BIN_SIZE;

  // Downmix to mono for analysis
  if (decoder.outputChannels > 1) {
    std::vector<float> mono(totalFrameCount);
    for (ma_uint64 frame = 0; frame < totalFrameCount; ++frame) {
      float sum = 0.0f;
      for (ma_uint32 ch = 0; ch < decoder.outputChannels; ++ch) {
        sum += buffer[frame * decoder.outputChannels + ch];
      }
      mono[frame] = sum / decoder.outputChannels;
    }
    inputRaw = std::move(mono);
  } else {
    inputRaw = buffer;
  }

  std::cout << "Read " << framesRead << " frames from " << filename << "\n";

  ma_decoder_uninit(&decoder);
}

void AnalysisCore::binInput() {
  size_t hopSize = inputBinSize / BIN_DECIMATION_FACTOR; // 75% overlap
  size_t startIdx = 0;

  // Goes until the end of inputRaw, avoids out-of-bounds
  while (startIdx < inputRaw.size()) {
    size_t endIdx = startIdx + inputBinSize;
    if (endIdx > inputRaw.size())
      endIdx = inputRaw.size();

    std::vector<float> bin(inputRaw.begin() + startIdx,
                           inputRaw.begin() + endIdx);

    // Optional: zero-pad to full inputBinSize
    if (bin.size() < inputBinSize) {
      bin.resize(inputBinSize, 0.0f);
    }

    bins.push_back(bin);
    startIdx += hopSize;
  }
}

void AnalysisCore::decomposeBins() {
  for (auto &bin : bins) {
    // 1. Apply window
    DSP::hannWindow(bin, bin.size());

    // 2. Prepare FFT output buffer
    std::vector<std::complex<float>> fftOutput(bin.size());

    // 3. Compute FFT
    DSP::computeFFT(bin, fftOutput, bin.size());

    // 5. Store
    binFFTs.push_back(std::move(fftOutput));
  }
}

void AnalysisCore::resynthesizeBin(size_t binIndex,
                                   std::vector<float> &output) {
  if (bins.empty()) {
    output.clear();
    return;
  }

  if (binIndex >= bins.size()) {
    binIndex = 0;
  }

  if (binFFTs.empty() || binIndex >= binFFTs.size()) {
    output = bins[binIndex];
    return;
  }

  const auto &spectrum = binFFTs[binIndex];
  size_t N = spectrum.size();

  // Extract magnitudes
  std::vector<float> mag(N);
  for (size_t k = 0; k < N; ++k) {
    mag[k] = std::abs(spectrum[k]);
  }

  std::vector<std::complex<float>> newSpectrum(N);
  constexpr float twoPi = 2.0f * static_cast<float>(M_PI);

  if (binPhaseAccumulators.size() < binFFTs.size()) {
    binPhaseAccumulators.resize(binFFTs.size());
  }
  if (binPhaseAccumulators[binIndex].size() != N) {
    binPhaseAccumulators[binIndex].assign(N, 0.0f);
  }

#if RANDOM_PHASE
  static std::mt19937 rng(static_cast<unsigned>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count()));
  static std::uniform_real_distribution<float> phaseDist(0.f, twoPi);
  for (size_t k = 0; k < N; ++k) {
    float phi = phaseDist(rng);
    newSpectrum[k] =
        std::complex<float>(mag[k] * std::cos(phi), mag[k] * std::sin(phi));
  }
#else
  for (size_t k = 0; k < N; ++k) {
    float phi = binPhaseAccumulators[binIndex][k];
    newSpectrum[k] =
        std::complex<float>(mag[k] * std::cos(phi), mag[k] * std::sin(phi));
    float phaseStep = twoPi * static_cast<float>(k) / static_cast<float>(N);
    phi += phaseStep;
    if (phi >= twoPi) {
      phi -= twoPi;
    }
    binPhaseAccumulators[binIndex][k] = phi;
  }
#endif

  DSP::computeIFFT(newSpectrum, output, N);

  // Normalize IFFT output
  for (size_t i = 0; i < N; ++i)
    output[i] /= static_cast<float>(N);

  // apply synthesis window (Hann)
  DSP::hannWindow(output, output.size());
}

void AnalysisCore::computeUmapFeatures() {
  // 1. compute power spectra
  // 2. compute stft (2 frames, decimation factor 2)
  // 3. compute mfccs (12), centroid, flux, rolloff
  for (size_t i = 0; i < binFFTs.size(); ++i) {
    const auto &spectrum = binFFTs[i];
    size_t N = spectrum.size();

    // Compute power spectrum
    std::vector<float> powerSpec;
    DSP::computePowerSpectrum(spectrum, powerSpec, N);

    // Compute mel-scale spectral envelope
    std::vector<float> melEnv;
    DSP::computeSpectralEnv(powerSpec, melEnv, DEVICE_SAMPLE_RATE, N, 40);

    //  Compute MFCCs
    std::vector<float> mfccs;
    DSP::computeMFCC(melEnv, mfccs, MFCC_COEFF_COUNT);

    // Compute centroid
    float centroid = 0.0f;
    DSP::computeCentroid(powerSpec, N, DEVICE_SAMPLE_RATE, centroid);

    // Compute STFT
    std::vector<std::vector<std::complex<float>>> stftFrames;
    DSP::computeSTFT(bins[i], stftFrames, N, 2);

    // Compute Flux
    float flux = 0.0f;
    if (stftFrames.size() >= 2) {
      std::vector<float> psCurr;
      std::vector<float> psPrev;
      DSP::computePowerSpectrum(stftFrames[0], psCurr, N);
      DSP::computePowerSpectrum(stftFrames[1], psPrev, N);
      DSP::computeFlux(psCurr, psPrev, N, flux);
    }

    // 6. Compute rolloff
    float rolloff = 0.0f;
    DSP::computeRolloff(powerSpec, N, DEVICE_SAMPLE_RATE, ROLLOFF_THRESHOLD,
                        rolloff);

    // Store features
    std::vector<float> features;
    features.insert(features.end(), mfccs.begin(), mfccs.end());
    features.push_back(centroid);
    features.push_back(flux);
    features.push_back(rolloff);

    umapFeatures.push_back(std::move(features));
  }
}

void AnalysisCore::computeUmapCoordinates() {
  std::vector<float> xCoords;
  std::vector<float> yCoords;
  DSP::computeUMAP(umapFeatures, xCoords, yCoords);
  binTimbreX = std::move(xCoords);
  binTimbreY = std::move(yCoords);
}

void AnalysisCore::computeVoronoiEdges() {
  const size_t count = binTimbreX.size();
  if (count == 0 || binTimbreY.size() != count) {
    return;
  }

  // ------------------------------------------------------------
  // 1. Convert points
  // ------------------------------------------------------------
  std::vector<jcv_point> points(count);
  for (size_t i = 0; i < count; ++i) {
    points[i].x = binTimbreX[i];
    points[i].y = binTimbreY[i];
  }

  // ------------------------------------------------------------
  // 2. Setup bounds (with padding so points don't sit on edges)
  // ------------------------------------------------------------
  float minX = binTimbreX[0];
  float maxX = binTimbreX[0];
  float minY = binTimbreY[0];
  float maxY = binTimbreY[0];
  for (size_t i = 1; i < count; ++i) {
    minX = std::min(minX, binTimbreX[i]);
    maxX = std::max(maxX, binTimbreX[i]);
    minY = std::min(minY, binTimbreY[i]);
    maxY = std::max(maxY, binTimbreY[i]);
  }

  jcv_rect bounds;
  bounds.min.x = minX;
  bounds.min.y = minY;
  bounds.max.x = maxX;
  bounds.max.y = maxY;

  // ------------------------------------------------------------
  // 3. Generate Voronoi diagram
  // ------------------------------------------------------------
  jcv_diagram diagram;
  memset(&diagram, 0, sizeof(jcv_diagram));

  jcv_diagram_generate((int)count, points.data(), &bounds, nullptr, &diagram);

  // ------------------------------------------------------------
  // 4. Extract edges (graph)
  // ------------------------------------------------------------
  voronoiEdges.clear();
  voronoiCells.clear();

  const jcv_site *sites = jcv_diagram_get_sites(&diagram);
  for (int i = 0; i < diagram.numsites; ++i) {
    const jcv_site *site = &sites[i];

    std::vector<jcv_point> cellPoints;
    for (const jcv_graphedge *e = site->edges; e; e = e->next) {
      // Line segment
      VoronoiEdge edge = VoronoiEdge{e->pos[0].x, e->pos[0].y, e->pos[1].x,
                                     e->pos[1].y, site->index};

      voronoiEdges.emplace_back(edge);

      cellPoints.push_back(e->pos[0]);
      cellPoints.push_back(e->pos[1]);
    }

    if (cellPoints.empty()) {
      continue;
    }

    std::vector<jcv_point> uniquePoints;
    uniquePoints.reserve(cellPoints.size());
    const float epsilon = 1e-5f;
    for (const auto &p : cellPoints) {
      bool exists = false;
      for (const auto &q : uniquePoints) {
        if (std::abs(p.x - q.x) < epsilon && std::abs(p.y - q.y) < epsilon) {
          exists = true;
          break;
        }
      }
      if (!exists) {
        uniquePoints.push_back(p);
      }
    }

    if (uniquePoints.size() < 3) {
      continue;
    }

    std::sort(uniquePoints.begin(), uniquePoints.end(),
              [&](const jcv_point &a, const jcv_point &b) {
                const float angleA =
                    std::atan2(a.y - site->p.y, a.x - site->p.x);
                const float angleB =
                    std::atan2(b.y - site->p.y, b.x - site->p.x);
                return angleA < angleB;
              });

    VoronoiCell cell;
    cell.site = site->index;
    cell.coords.reserve(uniquePoints.size() * 2);
    for (const auto &p : uniquePoints) {
      cell.coords.push_back(static_cast<float>(p.x));
      cell.coords.push_back(static_cast<float>(p.y));
    }
    voronoiCells.push_back(std::move(cell));
  }

  // ------------------------------------------------------------
  // 5. Cleanup
  // ------------------------------------------------------------
  jcv_diagram_free(&diagram);
}

std::vector<float> &AnalysisCore::getInputRaw() { return inputRaw; }

size_t AnalysisCore::getBinIndex() { return binIndex; }

void AnalysisCore::setBinIndex(size_t index) { binIndex = index; }

size_t AnalysisCore::getBinCount() const { return bins.size(); }

const std::vector<float> &AnalysisCore::getBin(size_t index) const {
  static const std::vector<float> empty;
  if (index >= bins.size()) {
    return empty;
  }
  return bins[index];
}

const std::vector<float> &AnalysisCore::getBinTimbreX() const {
  return binTimbreX;
}

const std::vector<float> &AnalysisCore::getBinTimbreY() const {
  return binTimbreY;
}

const std::vector<VoronoiEdge> &AnalysisCore::getVoronoiEdges() const {
  return voronoiEdges;
}

const std::vector<VoronoiCell> &AnalysisCore::getVoronoiCells() const {
  return voronoiCells;
}
