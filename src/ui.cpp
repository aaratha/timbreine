#include "ui.hpp"

#include <iostream>

#include "analysis.hpp"
#include "audio.hpp"

UiController::UiController(AudioCore *audioCore, AnalysisCore *analysisCore,
                           QObject *parent)
    : QObject(parent), audioCore(audioCore), analysisCore(analysisCore) {}

void UiController::rectangleClicked(int index) {
  lastClickedIndex = index;
  std::cerr << "rectangle clicked: " << index << "\n";
  emit lastClickedChanged();

  if (audioCore) {
    audioCore->setBinIndex(static_cast<size_t>(index));
  }
}

int UiController::lastClicked() const { return lastClickedIndex; }
