#pragma once

#include <atomic>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "miniaudio.h"
#include "analysis.hpp"

struct Oscillator {
  float phase;
  float frequency;
  float sampleRate;

  Oscillator() : phase(0.0f), frequency(440.0f), sampleRate(48000.0f) {}

  float process() {
    float value = sinf(phase * 2.0f * 3.14159265f);
    phase += frequency / sampleRate;
    if (phase >= 1.0f)
      phase -= 1.0f;
    return value;
  }
};

class AudioCore {
  AnalysisCore &analysisCore;

  Oscillator osc;

  ma_device_config deviceConfig{};
  ma_device device{};
  std::atomic<bool> running{false};
  bool audioInitialized{false};
  std::atomic<bool> callbackSeen{false};

  std::vector<float> inputRaw;
  std::atomic<size_t> playhead{0};

  std::atomic<size_t> binIndex{0};
  std::vector<float> phaseAccumulators;
  std::vector<float> noisePhaseAccumulators;
  std::vector<float> noiseFrequencies;
  size_t lastBinIndex{std::numeric_limits<size_t>::max()};

  static void dataCallback(ma_device *pDevice, void *pOutput,
                           const void * /*pInput*/, ma_uint32 frameCount);

public:
  AudioCore(AnalysisCore &analysisCore);
  ~AudioCore();

  void setBinIndex(size_t index);

  // Called per audio buffer
  void processAudio(float *out, ma_uint32 frameCount);

  // Optional: load samples
  void loadSample(const std::string &id, const std::vector<float> &data);
};
