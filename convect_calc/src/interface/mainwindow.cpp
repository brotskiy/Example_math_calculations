#include "mainwindow.h"

#include <QMenu>
#include <QCoreApplication>
#include <QFileInfo>
#include <QFile>
#include <QTimer>
#include <QThreadPool>
#include <QString>
#include <QStringList>
#include <QDebug>

#include <vector>
#include <numeric>
#include <cmath>
#include <utility>
#include <exception>

#include <iostream>
#include <iomanip>

#include "helper/EngineConfiger.h"
#include "helper/IThreadPlaceable.h"

#include "helper/Configer.h"
#include "helper/DataTextReader.h"
#include "helper/DataTextWriter2.h"
#include "helper/miscellanea.h"

#include "engine/RHSCreator.h"
#include "math/PlaneEquation.h"
#include "math/NewtonComputer.h"
#include "engine/RK4Computer.h"

namespace
{

QAction* GET_ACTION_OF_MENU(QMenuBar* menuBar, int menuIndex, int actionIndex)
{
    return menuBar->actions().at(menuIndex)->menu()->actions().at(actionIndex);
}

bool IS_GOOD_CONFIG(const QString &configPath)
{
    const QFileInfo configInfo(configPath);
    return configInfo.isFile() && configInfo.exists();
};

QString RK4ENGINE_CONFIG_PATH()
{
    return "config/config.ini";
}

void CHECK(bool isConnected)
{
    if (not isConnected)
        throw std::runtime_error("Connection can't be created!");
}

}

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

MainWindow::MainWindow(QWidget* parent) :
    QMainWindow(parent)
{
    setFixedSize(750, 750);
    createMenuBar();

//    onComputeIntersectionTime();
}

void MainWindow::createMenuBar()
{
    p_toolBar_ = new QToolBar(this);
    p_toolBar_->setMovable(false);
    addToolBar(Qt::TopToolBarArea, p_toolBar_);

    p_comboBox_ = new QComboBox(p_toolBar_);
    p_comboBox_->addItem("MASTER");
    p_comboBox_->setEditable(false);
    p_comboBox_->setEnabled(false);

    p_toolBar_->addWidget(p_comboBox_);

    p_txtEdit_ = new QAction(p_toolBar_);
    p_txtEdit_->setEnabled(false);

    p_toolBar_->addAction(p_txtEdit_);

    p_menuBar_ = new QMenuBar(this);
    setMenuBar(p_menuBar_);

    QMenu*const menuFile = new QMenu("Work", p_menuBar_);
    menuFile->addAction("Start full simulation");        // 0.
    menuFile->addAction("Start particle simulation");    // 1.
    menuFile->addSeparator();                            // 2.
    menuFile->addAction("Exit");                         // 3.

    p_menuBar_->addAction(menuFile->menuAction());       // menu 0.

    QAction*const buttonStartFull = GET_ACTION_OF_MENU(p_menuBar_, 0, 0);         // Start full.
    CHECK(connect(buttonStartFull, &QAction::triggered,
                  this,            &MainWindow::onStartFull));

    QAction*const buttonStartParticles = GET_ACTION_OF_MENU(p_menuBar_, 0, 1);    // Start particle.
    CHECK(connect(buttonStartParticles, &QAction::triggered,
                  this,                &MainWindow::onStartParticlesOnly));

    QAction*const buttonExit = GET_ACTION_OF_MENU(p_menuBar_, 0, 3);              // Exit.
    CHECK(connect(buttonExit,                   &QAction::triggered,
                  QCoreApplication::instance(), &QCoreApplication::quit));

// ----------------------------------------------------------------------------------------------
    QMenu*const menuPoincare = new QMenu("Poincare", p_menuBar_);
    menuPoincare->addAction("Plane section");    // 0.

    p_menuBar_->addAction(menuPoincare->menuAction());    // menu 1.

    QAction*const buttonPlane = GET_ACTION_OF_MENU(p_menuBar_, 1, 0);    // Plane section.
    CHECK(connect(buttonPlane, &QAction::triggered,
                  this,        &MainWindow::onPlaneSection));

// ----------------------------------------------------------------------------------------------
    QMenu*const menuStreamFunction = new QMenu("Stream function", p_menuBar_);
    p_menuBar_->addAction(menuStreamFunction->menuAction());                      // menu 2.

    menuStreamFunction->addAction("Stream function coefficients");    // 0.
    menuStreamFunction->addAction("Stream function values");          // 1.

    QAction*const buttonCoeffs = GET_ACTION_OF_MENU(p_menuBar_, 2, 0);
    CHECK(connect(buttonCoeffs, &QAction::triggered,
                  this,         &MainWindow::onStreamFunctionCoefficients));

    QAction*const buttonValues = GET_ACTION_OF_MENU(p_menuBar_, 2, 1);
    CHECK(connect(buttonValues, &QAction::triggered,
                  this,         &MainWindow::onStreamFunctionValues));
}

MainWindow::~MainWindow()
{
    qDebug() << "RK4Engine will be terminated!";
}

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

void MainWindow::onStartParticlesOnly()
{
    using tech::helper::readEngineSettings;
    using tech::helper::writeEngineSettings;
    using tech::Configer;

    const QString filePath_Settings = "config/prtcl_sim_config.ini";
    const QString filePath_Thetas   = "config/prtcl_sim_thetas.txt";

    const std::optional<RK4ParticleEngine::Settings> settings = tryGetRK4ParticleEngineSettings(filePath_Settings, filePath_Thetas);

    if (not settings.has_value())    // если не смогли считать файл настроек.
    {
        Configer<RK4ParticleEngine::Settings>(filePath_Settings, readEngineSettings, writeEngineSettings)
            .writeSettings(RK4ParticleEngine::Settings());

        return;    // создадим настройки по умолчанию и выходим!
    }

    for (size_t index = 1; index <= settings.value().initial_particles_.size(); index++)    // создаем вкладку для каждой частицы в выпадающем списке.
        p_comboBox_->addItem(QString::number(index));

    SceneWidget::Settings sceneSettings;
    sceneSettings.outputDirectory_ = "output/images";
    sceneSettings.a_               = settings.value().a_;
    sceneSettings.b_               = settings.value().b_;
    sceneSettings.particleCount_   = settings.value().initial_particles_.size();

    p_scene_ = new SceneWidget(sceneSettings, this);                 // создаем виджет с изображением контейнера.
    setCentralWidget(p_scene_);                                      // устанавливаем его в окно.

    CHECK(connect(p_comboBox_, &QComboBox::currentTextChanged,
                  p_scene_,    &SceneWidget::onSetCurrentImage));    // можно выбирать, за траекторией какой частицы будем следить.

    p_toolBar_->setEnabled(true);          // активируем выпадающий список частиц, за траекториями которых будем следить.
    p_menuBar_->setEnabled(false);         // деактивируем меню, чтобы не спамили кнопки.
    p_txtEdit_->setText("starting ..");

    QThreadPool::globalInstance()->start(
    [this, sttngs = settings.value()]() -> void
    {
        using math::NaturalConvection2D::Particle;
        using ParticlesWriter = tech::DataTextWriter2<double, std::vector<Particle>>;
        using tech::helper::timeParticlesToString;
        using tech::helper::makeUniqueFileName;

        const QString particlesFilePath = makeUniqueFileName("E:/EVENING/output/sim_particles.txt");
        ParticlesWriter particlesWriter(particlesFilePath, ParticlesWriter::OpenMode::kTruncate, timeParticlesToString);

        const auto storeParticles = [&particlesWriter](const double& time, const std::vector<Particle>& particles) -> void
        {
            particlesWriter.writeLine(time, particles);
        };

        const auto onStepCallback = [this](const int& step, const double& elapsed,
                                           const double& time, const std::vector<Particle>& particles) -> void    // в конце каждого шага алгоритма.
        {
            Q_UNUSED(time)

            if (step % 150 != 0)    // выводим информацию только для одного шага из заданного количества.
                return;

            QTimer::singleShot(1, this,
            [this, particles]() -> void        // испукая сообщение в очередь из другого потока.
            {
                this->p_scene_->onDrawParticles(particles);    // рисуем частицы на экране.
            });

            QTimer::singleShot(1, this,
            [this, step, elapsed]() -> void    // испукая сообщение в очередь из другого потока.
            {
                if (not this->p_comboBox_->isEnabled())
                    this->p_comboBox_->setEnabled(true);

                this->p_txtEdit_->setText(QString("step %1, %2 ms.")
                    .arg(QString::number(step), QString::number(elapsed)));    // выводим информацию о произведенном шаге.
            });
        };

        RK4ParticleEngine().simulateParticles(sttngs, storeParticles, onStepCallback);    // симулируем движение частиц методом Рунге-Кутты.

        QTimer::singleShot(1, this, [this](){ this->p_menuBar_->setEnabled(true); });     // ВКЛЮЧАЕМ кнопки меню через очередь, испустив сообщение из другого потока.
    });
}

void MainWindow::onStartFull()
{
    const QString configPath = RK4ENGINE_CONFIG_PATH();
    const std::optional<RK4Engine::Settings> engineSettings = tryGetRK4EngineSettings(configPath);

    if (not engineSettings.has_value())
        return;

    const auto factoryEngine = [&settings = engineSettings.value()]() { return new RK4Engine(settings, nullptr); };
    const auto onEngineStop  = []() { qDebug() << "RK4Engine is stopping!"; };

    p_rk4_engine_ = tech::heavy::create<RK4Engine>(this, factoryEngine, onEngineStop);    // создался в другом потоке.
    if (p_rk4_engine_)
    {
        GET_ACTION_OF_MENU(p_menuBar_, 0, 0)->setEnabled(false);    // отключаем кнопку Start.
        p_txtEdit_->setText("starting computations...");

        for (size_t index = 1; index <= engineSettings.value().particles_.size(); index++)
            p_comboBox_->addItem(QString::number(index));

        SceneWidget::Settings sceneSettings;
        sceneSettings.outputDirectory_ = "output/images";
        sceneSettings.a_               = engineSettings.value().a_;
        sceneSettings.b_               = engineSettings.value().b_;
        sceneSettings.particleCount_   = engineSettings.value().particles_.size();

        p_scene_ = new SceneWidget(sceneSettings, this);
        setCentralWidget(p_scene_);

        CHECK(connect(p_comboBox_, &QComboBox::currentTextChanged,
                      p_scene_,    &SceneWidget::onSetCurrentImage));

        CHECK(connect(p_rk4_engine_, &RK4Engine::particlesCoordinates,
                      p_scene_,      &SceneWidget::onDrawParticles));

        CHECK(connect(p_rk4_engine_, &RK4Engine::stepInfo, this,
        [this](int step, int elapsed_ms) -> void
        {
            if (not p_comboBox_->isEnabled())
                p_comboBox_->setEnabled(true);

            p_txtEdit_->setText(QString("step %1, %2 ms.")
                .arg(QString::number(step), QString::number(elapsed_ms)));
        }));

        QTimer::singleShot(500, p_rk4_engine_, [this](){ p_rk4_engine_->onStart(); });
    }
}

void MainWindow::onPlaneSection()
{
    const QString configPath = RK4ENGINE_CONFIG_PATH();
    const std::optional<RK4Engine::Settings> engineSettings = tryGetRK4EngineSettings(configPath);

    if (not engineSettings.has_value())
        return;

    const auto factoryEngine = [&settings = engineSettings.value()]() { return new RK4Engine(settings, nullptr); };
    const auto onEngineStop  = []() { qDebug() << "RK4Engine is stopping!"; };

    p_rk4_engine_ = tech::heavy::create<RK4Engine>(this, factoryEngine, onEngineStop);    // создался в другом потоке.
    if (p_rk4_engine_)
    {
        GET_ACTION_OF_MENU(p_menuBar_, 1, 0)->setEnabled(false);    // отключаем кнопку Plane section.
        p_txtEdit_->setText("starting reading file...");

        CHECK(connect(p_rk4_engine_, &RK4Engine::stepInfo, this,
        [this](int step, int elapsed_ms) -> void
        {
            if (not p_comboBox_->isEnabled())
                p_comboBox_->setEnabled(true);

            p_txtEdit_->setText(QString("step %1, %2 ms.")
                .arg(QString::number(step), QString::number(elapsed_ms)));
        }));

        QTimer::singleShot(500, p_rk4_engine_, [this](){ p_rk4_engine_->onPoincare(); });
    }
}

void MainWindow::onStreamFunctionCoefficients()
{
    const QString configPath = RK4ENGINE_CONFIG_PATH();
    const std::optional<RK4Engine::Settings> engineSettings = tryGetRK4EngineSettings(configPath);

    if (not engineSettings.has_value())
        return;

    const auto factoryEngine = [&settings = engineSettings.value()]() { return new RK4Engine(settings, nullptr); };
    const auto onEngineStop  = []() { qDebug() << "RK4Engine is stopping!"; };

    p_rk4_engine_ = tech::heavy::create<RK4Engine>(this, factoryEngine, onEngineStop);    // создался в другом потоке.
    if (p_rk4_engine_)
    {
        p_menuBar_->setEnabled(false);                          // отключаем все кнопки меню.
        p_txtEdit_->setText("starting reading theta file...");

        CHECK(connect(p_rk4_engine_, &RK4Engine::stepInfo, this,
        [this](int step, int elapsed_ms) -> void
        {
            p_txtEdit_->setText(QString("step %1, %2 ms.")
                .arg(QString::number(step), QString::number(elapsed_ms)));
        }));

        QTimer::singleShot(500, p_rk4_engine_, [this](){ p_rk4_engine_->onComputeStreamFunctionCoefficients(); });
    }
}

void MainWindow::onStreamFunctionValues()
{
    const QString configPath = RK4ENGINE_CONFIG_PATH();
    const std::optional<RK4Engine::Settings> engineSettings = tryGetRK4EngineSettings(configPath);

    if (not engineSettings.has_value())
        return;

    const auto factoryEngine = [&settings = engineSettings.value()]() { return new RK4Engine(settings, nullptr); };
    const auto onEngineStop  = []() { qDebug() << "RK4Engine is stopping!"; };

    p_rk4_engine_ = tech::heavy::create<RK4Engine>(this, factoryEngine, onEngineStop);    // создался в другом потоке.
    if (p_rk4_engine_)
    {
        p_menuBar_->setEnabled(false);                          // отключаем все кнопки меню.
        p_txtEdit_->setText("starting reading psi file...");

        CHECK(connect(p_rk4_engine_, &RK4Engine::stepInfo, this,
        [this](int step, int elapsed_ms) -> void
        {
            p_txtEdit_->setText(QString("step %1, %2 ms.")
                .arg(QString::number(step), QString::number(elapsed_ms)));
        }));

        QTimer::singleShot(500, p_rk4_engine_, [this](){ p_rk4_engine_->onComputeStreamFunction(); });
    }
}

void MainWindow::onComputeIntersectionTime()
{
    using math::NaturalConvection2D::RHSCreator;
    using math::RK4Computer;
    using math::NewtonComputer;
    using ThetasReader = tech::DataTextReader<double, std::vector<double>>;
    using ThetasWriter2 = tech::DataTextWriter2<double, std::vector<double>>;

    const double RAYLEIGH = 0.23;
    const double A = 20;
    const double B = 70;
    const int MAX_I = 10, MAX_J = 10;

    const double EPSILON = 1e-10;
    const int STEP_MAX   = 100;

    const auto stringToTimeThetas = [](const QString &str, double &time_out, std::vector<double> &thetas_out) -> void
    {
        const QStringList list_time_values = str.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);

        if (list_time_values.size() <= 1)
            throw std::runtime_error("<stringToPairTimeValues()>: invalid line has been read!");

        time_out = list_time_values.front().toDouble();
        thetas_out.resize(list_time_values.size() - 1);

        std::transform(std::next(list_time_values.cbegin()), list_time_values.cend(), thetas_out.begin(),
        [](const QString& value_string) -> double
        {
            return value_string.toDouble();
        });
    };

    const QString path_Start  = "config/thetas_start.txt";
    const QString path_Before = "config/thetas_before.txt";
    const QString path_After  = "config/thetas_after.txt";
    const QString path_finish = "output/thetas_finish.txt";

    double T_Start;
    std::vector<double> THETA_Start;
    ThetasReader(path_Start, stringToTimeThetas).readLineInto(T_Start, THETA_Start);

    double T_Before;
    std::vector<double> THETA_Before;
    ThetasReader(path_Before, stringToTimeThetas).readLineInto(T_Before, THETA_Before);

    double T_After;
    std::vector<double> THETA_After;
    ThetasReader(path_After, stringToTimeThetas).readLineInto(T_After, THETA_After);

    const math::PlaneEquation Plane(std::vector<double>(THETA_Start.size(), 1), THETA_Start);

    const auto Diff = [&RAYLEIGH, rhs = RHSCreator::createRHS(A, B, MAX_I, MAX_J)->diffs_]
                      (const double &time, const std::vector<double> &thetas, std::vector<double> &thetasDiff_out) -> void
    {
        Q_UNUSED(time);

        if (rhs.size() != thetas.size() || thetasDiff_out.size() != thetas.size())
            throw std::runtime_error("<MainWindow::onComputeIntersectionTime>:dimension missmatch!");

        const auto thetaIndexToArgument = [&thetas](const int index) { return thetas.at(index); };

        #pragma omp parallel for default(shared)
        for (size_t eqIndex = 0; eqIndex < thetasDiff_out.size(); eqIndex++)
        {
            thetasDiff_out.at(eqIndex) = RK4Engine::RK4_computeDiffEquation(rhs.at(eqIndex), RAYLEIGH, thetaIndexToArgument);
        }
    };

    const auto Estim = [&T_Before, &THETA_Before, &Plane, &Diff](const double &t) -> double
    {
        std::cout << "T_Finish = " << std::setprecision(20) << t << std::endl;

        const double h = t - T_Before;
        if (h < 0)
            throw std::runtime_error("<MainWindow::onComputeIntersectionTime>: rk4 step is less than zero!");

        std::vector<double> THETA_Finish(THETA_Before.size());
        RK4Computer().RK4(Diff, h, T_Before, THETA_Before, THETA_Finish);

        const double distance = Plane.substitute(THETA_Finish);
        return distance;
    };

    const auto EstimEuc = [&T_Before, &THETA_Before, &THETA_Start, &Diff](const double &t) -> double
    {
        std::cout << "T_Finish = " << std::setprecision(20) << t << std::endl;
        const double h = t - T_Before;
        if (h < 0)
            throw std::runtime_error("<MainWindow::onComputeIntersectionTime()>: rk4 step is less than zero!");
        std::vector<double> THETA_Finish(THETA_Before.size());
        RK4Computer().RK4(Diff, h, T_Before, THETA_Before, THETA_Finish);

        double distance = 0.0;
        for (size_t i = 0; i < THETA_Start.size(); i++)
        {
            distance += (THETA_Start.at(i) - THETA_Finish.at(i)) * (THETA_Start.at(i) - THETA_Finish.at(i));
        }
        distance = std::sqrt(distance);

        return distance;
    };

    const auto dEstim = [&T_Before, &THETA_Before, &Diff, normal = Plane.getNormal()](const double &t) -> double
    {
        if (normal.size() != THETA_Before.size())
            throw std::runtime_error("<MainWindow::onComputeIntersectionTime>: dimension missmatch!");

        const double h = t - T_Before;
        if (h < 0)
            throw std::runtime_error("<MainWindow::onComputeIntersectionTime>: rk4 step is less than zero!");

        std::vector<double> THETA_Finish_(THETA_Before.size());
        RK4Computer().RK4(Diff, h, T_Before, THETA_Before, THETA_Finish_);

        std::vector<double> diff_THETA_Finish_(THETA_Finish_.size());
        Diff(t, THETA_Finish_, diff_THETA_Finish_);

        const double distance = std::inner_product(normal.cbegin(), normal.cend(), diff_THETA_Finish_.cbegin(), 0.0);
        return distance;
    };

    const auto dEstimEuc = [&T_Before, &THETA_Before, &THETA_Start, &Diff](const double &t) -> double
    {
        const double h = t - T_Before;
        if (h < 0)
            throw std::runtime_error("<MainWindow::onComputeIntersectionTime>: rk4 step is less than zero!");

        std::vector<double> THETA_Finish(THETA_Before.size());
        RK4Computer().RK4(Diff, h, T_Before, THETA_Before, THETA_Finish);

        std::vector<double> d_THETA_Finish(THETA_Finish.size());
        Diff(t, THETA_Finish, d_THETA_Finish);

        std::vector<double> d_THETA_Start(THETA_Start.size());
        Diff(t, THETA_Start, d_THETA_Start);

        double numerator   = 0.0;
        double denominator = 0.0;
        for (size_t i = 0; i < THETA_Start.size(); i++)
        {
            const double first  = THETA_Start.at(i) - THETA_Finish.at(i);
            const double second = d_THETA_Start.at(i) - d_THETA_Finish.at(i);
            numerator   += first * second;
            denominator += (THETA_Start.at(i) - THETA_Finish.at(i)) * (THETA_Start.at(i) - THETA_Finish.at(i));
        }
        denominator = std::sqrt(denominator);

        const double distance = numerator / denominator;

        return distance;
    };

    const double T_Finish_real = NewtonComputer::computeRoot(T_After, Estim, dEstim, EPSILON, STEP_MAX);
    //const double T_Finish_real = NewtonComputer::computeRoot(T_Before, EstimEuc, dEstimEuc, EPSILON, STEP_MAX);

    const double Step = T_Finish_real - T_Before;
    std::vector<double> THETA_Finish_real(THETA_Before.size());
    RK4Computer().RK4(Diff, Step, T_Before, THETA_Before, THETA_Finish_real);

    ThetasWriter2(path_finish, ThetasWriter2::OpenMode::kTruncate,
    [](const double &time, const std::vector<double> &thetas) -> QString
    {
        QString line = QString::number(time, 'e', 20);

        for (const double &theta : thetas)
            line.append(" " + QString::number(theta, 'e', 20));

        return line;
    })
    .writeLine(T_Finish_real, THETA_Finish_real);

    std::cout  << "FINAL time = " << std::setprecision(20) << T_Finish_real << " | FINAL step = "<< Step << std::endl;
}

// :::::::::::::::::::::::::::::::::::::::: Служебное. ::::::::::::::::::::::::::::::::::::::::

std::optional<RK4Engine::Settings> MainWindow::tryGetRK4EngineSettings(const QString& engineConfigPath)
{
    if (engineConfigPath.isEmpty())
        throw std::invalid_argument("<MainWindow>: Runge-Kutta engine config path is empty!");

    if (IS_GOOD_CONFIG(engineConfigPath))    // если файл конфига существует.
    {
        RK4Engine::Settings settings;

        if (tech::EngineConfiger(engineConfigPath).readEngineConfig(settings))    // и мы смогли его прочитать.
            return settings;                                                      // то возвратим настройки.
    }

    QFile::remove(engineConfigPath);    // если по какой-то причине не удалось прочитать существующий файл (или его нет), то сначала (возможно) удаляем его.

    tech::EngineConfiger(engineConfigPath).writeEngineConfig(RK4Engine::Settings());    // а затем создаем конфиг по умолчанию.

    return std::nullopt;
}

std::optional<RK4ParticleEngine::Settings> MainWindow::tryGetRK4ParticleEngineSettings(const QString &filePath_Settings,
                                                                                       const QString &filePath_Thetas)
{
    using tech::helper::readEngineSettings;
    using tech::helper::writeEngineSettings;
    using tech::helper::stringToPairTimeValues;
    using tech::Configer;
    using tech::DataTextReader;

    if (not QFileInfo::exists(filePath_Settings))    // если нет файла настроек, то даже не будем пытаться что-либо делать.
        return std::nullopt;

    if (not QFileInfo::exists(filePath_Thetas))
        throw std::runtime_error("<MainWindow::tryGetRK4ParticleEngineSettings()>: thetas file doesn't exist!");

    std::optional<RK4ParticleEngine::Settings> settings = RK4ParticleEngine::Settings();

    Configer<RK4ParticleEngine::Settings> configReader(filePath_Settings, readEngineSettings, writeEngineSettings);
    configReader.readSettings(settings.value());

    DataTextReader<std::pair<double, std::vector<double>>> thetasReader(filePath_Thetas, stringToPairTimeValues);
    std::pair<double, std::vector<double>> time_values;

    while (thetasReader.readLineInto(time_values))
    {
        auto &times_thetas = settings.value().times_thetas_;    // ссылка-псевдоним!

        if (not times_thetas.empty())    // если в наборе заданных на периоде THETA уже содержится хоть одно значение.
        {
            const size_t prevSize = times_thetas.back().second.size();
            const size_t curSize  = time_values.second.size();

            if (prevSize != curSize)
                throw std::runtime_error("<MainWindow::tryGetRK4ParticleEngineSettings()>: dimension missmatch!");    // сравним длину последних считанных THETA.
        }

        times_thetas.push_back(std::move(time_values));
    }

    return settings;
}
