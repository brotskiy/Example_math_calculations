#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QObject>
#include <QMainWindow>
#include <QMenuBar>
#include <QToolBar>
#include <QAction>
#include <QComboBox>

#include <optional>

#include "engine/mainengine.h"
#include "engine/RK4ParticleEngine.h"
#include "interface/widget.h"

class MainWindow : public QMainWindow
{
  Q_OBJECT

// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

private:
    RK4Engine* p_rk4_engine_ = nullptr;    // является thread-placeable!
    SceneWidget* p_scene_    = nullptr;    // за удаление отвечает объект-родитель!

    QMenuBar* p_menuBar_   = nullptr;
    QToolBar* p_toolBar_   = nullptr;
    QAction* p_txtEdit_    = nullptr;
    QComboBox* p_comboBox_ = nullptr;

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    void createMenuBar();

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

private slots:
    void onStartFull();
    void onStartParticlesOnly();
    void onPlaneSection();
    void onStreamFunctionCoefficients();
    void onStreamFunctionValues();

    void onComputeIntersectionTime();

// :::::::::::::::::::::::::::::::::::::::: Служебное. ::::::::::::::::::::::::::::::::::::::::

private:
    static std::optional<RK4Engine::Settings> tryGetRK4EngineSettings(const QString& engineConfigPath);
    static std::optional<RK4ParticleEngine::Settings> tryGetRK4ParticleEngineSettings(const QString &filePath_Settings,
                                                                                      const QString &filePath_Thetas);

// :::::::::::::::::::::::::::::::::::::::: Ограничения. ::::::::::::::::::::::::::::::::::::::::

private:
    MainWindow()                             = delete;
    MainWindow(const MainWindow&)            = delete;
    MainWindow(MainWindow&&)                 = delete;
    MainWindow& operator=(const MainWindow&) = delete;
    MainWindow& operator=(MainWindow&&)      = delete;
};

#endif // MAINWINDOW_H
