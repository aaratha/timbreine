#include <iostream>
#include <chrono>
#include <random>

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

    bins.push_back(bin);
    startIdx += hopSize;
  }
}

void AnalysisCore::decomposeBins() {
  for (auto &bin : bins) {
    // 1. Apply window
    DSP::window(bin, bin.size());

    // 2. Prepare FFT output buffer
    std::vector<std::complex<float>> fftOutput(bin.size()); 

    // 3. Compute FFT
    DSP::computeFFT(bin, fftOutput, bin.size());


    // 5. Store
    binFFTs.push_back(std::move(fftOutput));
  }
}

void AnalysisCore::resynthesizeBin(size_t binIndex, std::vector<float> &output) {
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

    // Generate random phase
    static std::mt19937 rng(
        static_cast<unsigned>(std::chrono::high_resolution_clock::now()
                                  .time_since_epoch()
                                  .count()));
    static std::uniform_real_distribution<float> phaseDist(0.f, 2.f * M_PI);

    std::vector<std::complex<float>> newSpectrum(N);
    for (size_t k = 0; k < N; ++k) {
        float phi = phaseDist(rng);
        newSpectrum[k] = std::complex<float>(mag[k] * std::cos(phi),
                                             mag[k] * std::sin(phi));
    }

    DSP::computeIFFT(newSpectrum, output, N);
    
    // Normalize IFFT output
    for (size_t i = 0; i < N; ++i)
        output[i] /= static_cast<float>(N);

    // apply synthesis window (Hann)
    DSP::window(output, output.size());
}


std::vector<float> &AnalysisCore::getInputRaw() { return inputRaw; }
const std::vector<BinFeatures> &AnalysisCore::getBinFeatures() const {
  return binFeatures;
}

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
