#pragma once

#include <QObject>

class AnalysisCore;
class AudioCore;

class UiController : public QObject {
  Q_OBJECT
  Q_PROPERTY(int lastClicked READ lastClicked NOTIFY lastClickedChanged)

public:
  explicit UiController(AudioCore *audioCore, AnalysisCore *analysisCore,
                        QObject *parent = nullptr);

  Q_INVOKABLE void rectangleClicked(int index);
  int lastClicked() const;

signals:
  void lastClickedChanged();

private:
  AudioCore *audioCore;
  AnalysisCore *analysisCore;
  int lastClickedIndex{-1};
};
