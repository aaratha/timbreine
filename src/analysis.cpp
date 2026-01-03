#include <iostream>

#include "analysis.hpp"
#include "dsp.hpp"
#include "globals.hpp"
#include "miniaudio.h"

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

  inputBinSize = INPUT_BUFFER_SIZE; // buffer.size() / 24;

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
  size_t hopSize = inputBinSize / 4; // 75% overlap
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

    binnedInput.push_back(bin);
    startIdx += hopSize;
  }
}

void AnalysisCore::windowBins() {
  for (auto &bin : binnedInput) {
    size_t N = bin.size();
    for (size_t n = 0; n < N; ++n) {
      // Hann window
      float w = 0.5f * (1.0f - cosf((2.0f * 3.14159265f * n) / (N - 1)));
      bin[n] *= w;
    }
  }
}

void AnalysisCore::decomposeBins() {
  for (const auto &bin : binnedInput) {
    std::vector<std::vector<float>> frameSpectrums;
    std::vector<float> powerSpectrum;
    DSP::computePS(bin, frameSpectrums, powerSpectrum,
                   static_cast<int>(bin.size()));

    // store overall bin PS and per-frame PS for flux later
    binPS.push_back(powerSpectrum);
    framePS.insert(framePS.end(), frameSpectrums.begin(), frameSpectrums.end());
  }
}

void AnalysisCore::findSynthesisFeatures() {
  for (const auto &binSpectrum : binPS) {
    std::vector<float> freqComponents;
    std::vector<float> spectralEnv;

    DSP::computeSignificantFreqs(
        binSpectrum, freqComponents, static_cast<int>(DEVICE_SAMPLE_RATE));
    DSP::computeSpectralEnv(binSpectrum, spectralEnv,
                            static_cast<int>(DEVICE_SAMPLE_RATE));

    binFreqComponents.push_back(freqComponents);
    binSpectralEnv.push_back(spectralEnv);
  }
}

std::vector<float> &AnalysisCore::getInputRaw() { return inputRaw; }
std::vector<std::vector<float>> &AnalysisCore::getBinFreqComponents() {
  return binFreqComponents;
}
std::vector<std::vector<float>> &AnalysisCore::getBinSpectralEnvs() {
  return binSpectralEnv;
}
