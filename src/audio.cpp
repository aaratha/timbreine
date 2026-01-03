#include "audio.hpp"
#include "globals.hpp"

#include <iostream>

AudioCore::AudioCore(AnalysisCore &analysisCore) : analysisCore(analysisCore) {
  if (audioInitialized)
    return;

  deviceConfig = ma_device_config_init(ma_device_type_playback);
  deviceConfig.playback.format = DEVICE_FORMAT;
  deviceConfig.playback.channels = DEVICE_CHANNELS;
  deviceConfig.sampleRate = static_cast<ma_uint32>(DEVICE_SAMPLE_RATE);
  deviceConfig.dataCallback = AudioCore::dataCallback;
  deviceConfig.pUserData = this;

  osc.sampleRate = DEVICE_SAMPLE_RATE;

  if (ma_device_init(NULL, &deviceConfig, &device) != MA_SUCCESS) {
    std::cerr << "ma_device_init failed\n";
    return;
  }

  if (ma_device_start(&device) != MA_SUCCESS) {
    std::cerr << "ma_device_start failed\n";
    ma_device_uninit(&device);
    return;
  }

  playhead.store(0);
  audioInitialized = true;
  running.store(true);
}

AudioCore::~AudioCore() {
  if (audioInitialized) {
    ma_device_uninit(&device);
    audioInitialized = false;
    running.store(false);
  }
}

void AudioCore::setBinIndex(size_t index) {
  const auto &bins = analysisCore.getBinFreqComponents();
  if (bins.empty()) {
    binIndex.store(0);
    return;
  }

  if (index >= bins.size()) {
    index = bins.size() - 1;
  }

  binIndex.store(index);
}

void AudioCore::dataCallback(ma_device *pDevice, void *pOutput,
                             const void * /*pInput*/, ma_uint32 frameCount) {
  auto *core = static_cast<AudioCore *>(pDevice->pUserData);
  float *out = static_cast<float *>(pOutput);

  if (!core->callbackSeen.exchange(true)) {
    std::cerr << "audio callback active\n";
  }

  core->processAudio(out, frameCount);
}

void AudioCore::processAudio(float *out, ma_uint32 frameCount) {
  // for (ma_uint32 frame = 0; frame < frameCount; ++frame) {
  // // Sample playback functionality
  // if (!analysisCore.getInputRaw().empty()) {
  //   // use playhead to track position in inputRaw
  //   size_t idx = playhead.load();
  //   size_t channels = static_cast<size_t>(DEVICE_CHANNELS);
  //   size_t total = analysisCore.getInputRaw().size();
  //   for (ma_uint32 ch = 0; ch < DEVICE_CHANNELS; ++ch) {
  //     size_t sampleIndex = (idx + ch) % total;
  //     *out++ = analysisCore.getInputRaw()[sampleIndex];
  //   }
  //   idx = (idx + channels) % total;
  //   playhead.store(idx);
  // } else {
  //   // Generate a simple sine wave if no input file is loaded
  //   float sample = 0.2f * osc.process();
  //   for (ma_uint32 ch = 0; ch < DEVICE_CHANNELS; ++ch) {
  //     *out++ = sample;
  //   }
  // }
  // Reproduce sample bin based on rectangle clicked
  // }
  const auto &binFreqs = analysisCore.getBinFreqComponents();
  const auto &binEnvs = analysisCore.getBinSpectralEnvs();
  if (binFreqs.empty() || binEnvs.empty()) {
    for (ma_uint32 frame = 0; frame < frameCount; ++frame) {
      for (ma_uint32 ch = 0; ch < DEVICE_CHANNELS; ++ch) {
        *out++ = 0.0f;
      }
    }
    return;
  }

  size_t currentBin = binIndex.load();
  if (currentBin >= binFreqs.size()) {
    currentBin = 0;
  }
  if (currentBin >= binEnvs.size()) {
    currentBin = 0;
  }

  const auto &freqs = binFreqs[currentBin]; // dominant freqs
  const auto &env = binEnvs[currentBin];    // single envelope per bin
  size_t channels = DEVICE_CHANNELS;

  if (currentBin != lastBinIndex || phaseAccumulators.size() != freqs.size()) {
    phaseAccumulators.assign(freqs.size(), 0.0f);
    noisePhaseAccumulators.assign(env.size(), 0.0f);
    noiseFrequencies.resize(env.size());
    if (!env.empty()) {
      for (size_t i = 0; i < env.size(); ++i) {
        float t = (static_cast<float>(i) + 0.5f) / env.size();
        noiseFrequencies[i] = t * (DEVICE_SAMPLE_RATE / 2.0f);
      }
    }
    lastBinIndex = currentBin;
  }

  for (ma_uint32 frame = 0; frame < frameCount; ++frame) {
    float sample = 0.0f;
    const float harmonicGain = 0.1f;
    const float noiseGain = 0.02f;

    // For each significant frequency, generate sine wave weighted by envelope
    for (size_t i = 0; i < freqs.size(); ++i) {
      float f = freqs[i]; // frequency in Hz

      // Map frequency to envelope bin (assuming envelope = N/2+1 bins)
      size_t envBin =
        static_cast<size_t>(f / (DEVICE_SAMPLE_RATE / 2.0f) * env.size());
      if (envBin >= env.size())
        envBin = env.size() - 1;
      float A = expf(env[envBin]); // amplitude from spectral envelope

      // Track phase for continuity
      float &phase = phaseAccumulators[i];
      sample += harmonicGain * A * sinf(phase);
      phase += 2.0f * M_PI * f / DEVICE_SAMPLE_RATE;
      if (phase > 2.0f * M_PI)
        phase -= 2.0f * M_PI;
    }

    // Add noise shaped by the spectral envelope
    for (size_t i = 0; i < env.size(); ++i) {
      float A = expf(env[i]);
      float &phase = noisePhaseAccumulators[i];
      float f = noiseFrequencies[i];
      sample += noiseGain * A * sinf(phase);
      phase += 2.0f * M_PI * f / DEVICE_SAMPLE_RATE;
      if (phase > 2.0f * M_PI)
        phase -= 2.0f * M_PI;
    }

    // Write the sample to all output channels
    for (ma_uint32 ch = 0; ch < channels; ++ch)
      *out++ = sample;
  }
}
