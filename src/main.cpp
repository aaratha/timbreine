#include <QGuiApplication>
#include <QQmlApplicationEngine>
#include <QQmlContext>

#include "audio.hpp"
#include "ui.hpp"
#include "timbre.hpp"

int main(int argc, char *argv[]) {
    QGuiApplication app(argc, argv);

    QQmlApplicationEngine engine;

    AnalysisCore analysisCore;
    if (argc > 1) {
        analysisCore.readFile(argv[1]); // Load audio file from command line argument
        analysisCore.binInput();
        analysisCore.windowBins();
        analysisCore.decomposeBins();
        analysisCore.findSynthesisFeatures();
    }

    AudioCore audioCore(analysisCore);
    UiController uiController(&audioCore, &analysisCore);
    engine.rootContext()->setContextProperty("ui", &uiController);
    
    engine.loadFromModule("MyApp", "MainView");

    if (engine.rootObjects().isEmpty())
        return -1;

    return app.exec();
}
