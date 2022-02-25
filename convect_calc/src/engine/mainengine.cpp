#include "mainengine.h"

#include <QElapsedTimer>
#include <QDebug>

#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <exception>
#include <optional>
#include <utility>
#include <iterator>

#include "RHSCreator.h"
#include "helper/DataWriter.h"
#include "helper/miscellanea.h"

#include "helper/DataReaderWriterPoincare.h"
#include "helper/DataTextWriter.h"

Q_DECLARE_METATYPE(std::vector<math::NaturalConvection2D::Particle>)

namespace
{

const std::vector<double> INITIALS =
{
-4.63284e-19
 ,38.6988
 ,-2.98969e-19
 ,23.5746
 ,-3.42062e-19
 ,19.537
 ,4.22416e-20
 ,9.85456
 ,-1.81796e-19
 ,11.4874
 ,-3.22674
 ,-2.77009e-19
 ,7.61948
 ,6.68391e-19
 ,-5.18037
 ,-5.03296e-19
 ,1.23528
 ,2.21027e-19
 ,-0.264903
 ,-3.94445e-20
 ,4.16586e-19
 ,-20.3969
 ,-1.4133e-19
 ,4.404
 ,6.09555e-20
 ,-1.6562
 ,-7.10071e-20
 ,2.33034
 ,-3.85437e-20
 ,-2.42196
 ,-0.702664
 ,1.11383e-19
 ,-1.23792
 ,-2.37579e-19
 ,0.948973
 ,1.84944e-19
 ,0.200405
 ,-4.19336e-20
 ,-0.142047
 ,-2.06295e-20
 ,3.05559e-20
 ,-1.63442
 ,-4.90378e-21
 ,-1.88966
 ,1.46869e-20
 ,-0.518299
 ,-1.32844e-20
 ,-1.26727
 ,2.30745e-20
 ,-0.471605
 ,0.191256
 ,-1.37501e-20
 ,-0.193361
 ,1.29185e-20
 ,0.254859
 ,-4.18143e-21
 ,-0.138631
 ,-9.97758e-21
 ,-0.00366688
 ,3.62473e-21
 ,-6.16507e-22
 ,-0.862342
 ,1.6401e-20
 ,-0.602381
 ,1.78969e-21
 ,-0.499712
 ,5.37459e-21
 ,-0.237328
 ,3.39759e-21
 ,-0.0806903
 ,-0.0100364
 ,4.33619e-21
 ,-0.0489198
 ,-9.77955e-21
 ,0.0409536
 ,4.28166e-21
 ,-0.0313764
 ,-1.55355e-21
 ,0.0199315
 ,1.01582e-21
 ,5.94397e-21
 ,-0.334496
 ,3.62851e-21
 ,-0.315815
 ,4.42177e-21
 ,-0.1874
 ,1.51181e-21
 ,-0.118578
 ,2.20906e-21
 ,-0.0886138
 ,0.0106497
 ,1.3683e-22
 ,-0.0376783
 ,-3.08722e-21
 ,0.0316891
 ,2.77231e-21
 ,-0.00698332
 ,-1.42789e-21
 ,0.00325285
 ,5.12599e-22
};

}

// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

RK4Engine::Settings::Settings()
{
    a_ = 20.0;
    b_ = 46.0;

    t_beg_ = 0.0;
    t_end_ = 500.0;

    rayleigh_ = 0.119;

    eq_count_   = 100;
    step_count_ = 50000;

    repeat_count_ = 1;

    file_path_thetas_   = "output/thetas/thetas.txt";
    file_path_particle_ = "output/particles/particle.txt";

    particles_ = { {5, 12}, {10, 12}, {15, 12},
                   {5, 24}, {10, 24}, {15, 24},
                   {5, 36}, {10, 36}, {15, 36} };
}

void RK4Engine::printSettings(const RK4Engine::Settings &sttngs)
{
    qDebug()
        << "a            = " << sttngs.a_           << Qt::endl
        << "b            = " << sttngs.b_           << Qt::endl
        << "t_beg        = " << sttngs.t_beg_       << Qt::endl
        << "t_end        = " << sttngs.t_end_       << Qt::endl
        << "Rayleigh     = " << sttngs.rayleigh_    << Qt::endl
        << "eq_count     = " << sttngs.eq_count_    << Qt::endl
        << "step_count   = " << sttngs.step_count_  << Qt::endl
        << "repeat_count = " << sttngs.repeat_count_<< Qt::endl;
}

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

RK4Engine::RK4Engine(const RK4Engine::Settings &engnSttngs, QObject* parent) :
    QObject(parent),
    settings_(engnSttngs)
{
    qRegisterMetaType<std::vector<math::NaturalConvection2D::Particle>>();
}

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

void RK4Engine::onStart()
{
    using math::NaturalConvection2D::SystemRHS;
    using math::NaturalConvection2D::SystemData;
    using math::NaturalConvection2D::RHSCreator;

    printSettings(settings_);

    const int max_i = static_cast<int>(std::sqrt(settings_.eq_count_));
    const int max_j = max_i;

    const SystemRHS rhs = *RHSCreator::createRHS(settings_.a_, settings_.b_, max_i, max_j);
    SystemData data     = *createData(settings_);

    std::vector<std::vector<double>> KThetaBuffers = *create_K_coefficients(settings_.eq_count_);    // по 4 коэффициента метода Рунге-Кутты для каждой THETA в системе.

    const double h = (settings_.t_end_ - settings_.t_beg_) / settings_.step_count_;

    qDebug() << "RK4 start!";

    for (int repeat = 0; repeat < settings_.repeat_count_; repeat++)
    {
        maybeEmitInfo(repeat * settings_.step_count_, 0, data.coords_.at(0));    // отображаем начальные шаги каждого повтора.

        for (int step = 1; step <= settings_.step_count_; step++)
        {
            QElapsedTimer timer;
            timer.start();

            RK4_step(rhs, settings_.rayleigh_, h, step, KThetaBuffers, data);

            const int realStep   = repeat * settings_.step_count_ + step;
            const int elapsed_ms = timer.elapsed();
            maybeEmitInfo(realStep, elapsed_ms, data.coords_.at(step));
        }

        const size_t endIndex = data.time_.size() - 1;    // на самом последнем повторе не будет записан самый последний шаг - ничего страшного.
        writeData(data, endIndex);
        cycleData(data);
    }

    qDebug() << "RK4 finish!";
}

void RK4Engine::onPoincare()
{
    using math::NaturalConvection2D::SystemRHS;
    using math::NaturalConvection2D::RHSCreator;
    using math::PoincareSection::PlaneEquation;
    using math::NaturalConvection2D::Particle;

    printSettings(settings_);

    const QString INPUT_FILE  = "E:/TORUS/LAST/output/thetas/thetas (0.230000).txt.4768";    // TODO сделать нормально!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    const QString OUTPUT_FILE = "E:/TORUS/LAST/output/poincare/poincare_thetas.txt";             // TODO сделать нормально!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    const QString INPUT_PRT_FILE  = "E:/TORUS/LAST/output/particles/particles (0.230000).txt.4768";    // TODO сделать нормально!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    const QString OUTPUT_PRT_FILE = "E:/TORUS/LAST/output/poincare/poincare_particles.txt";    // TODO сделать нормально!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    const int max_i = static_cast<int>(std::sqrt(settings_.eq_count_));
    const int max_j = max_i;
    const SystemRHS rhs = *RHSCreator::createRHS(settings_.a_, settings_.b_, max_i, max_j);    // создаем систему дифференциальных уравнений.
    const PlaneEquation planeEquation = *PN_createPlane(settings_.eq_count_);                  // создаем уравнение секущей плоскости.

    tech::PoincareSection::DataReaderWriterPoincare thetaProvider(INPUT_FILE, OUTPUT_FILE);
    tech::PoincareSection::DataReaderWriterPoincare particleProvider(INPUT_PRT_FILE, OUTPUT_PRT_FILE);

    const auto getData = [&thetaProvider, &particleProvider](double& time_out,
                                                             std::vector<double>& thetas_out,
                                                             std::vector<Particle>& particles_out) -> bool
    {
        const bool okThetas = thetaProvider.readLineInto(time_out, thetas_out);

        double time;
        const bool okParticles = particleProvider.readLineInto(time,particles_out);

        if (std::abs(time_out - time) >= 0.1)
            throw std::runtime_error("INVALID PARTICLES TIME!");

        return okThetas && okParticles;
    };

    const auto storeData = [&thetaProvider, &particleProvider](const double& time,
                                                               const std::vector<double> &thetas,
                                                               const std::vector<Particle> &particles) -> void
    {
        thetaProvider.writeLine(time, thetas);
        particleProvider.writeLine(time, particles);
    };

    qDebug() << "Poincare start!";

    bool canProceed = true;   // значение, чтобы войти в цикл.
    int step = 0;
    while (canProceed)        // пока мы можем продолжать выполнение шагов.
    {
        QElapsedTimer timer;
        timer.start();

        canProceed = PN_step(getData, storeData, settings_.rayleigh_, planeEquation, rhs);

        if (step % 500 == 0)
            emit stepInfo(step, timer.elapsed());

        step++;
    }

    qDebug() << "Poincare finish!";
}

void RK4Engine::onComputeStreamFunctionCoefficients()
{
    using math::NaturalConvection2D::Equation1Term;
    using math::NaturalConvection2D::RHSCreator;
    using tech::PoincareSection::DataReaderWriterPoincare;

    printSettings(settings_);

    const QString INPUT_FILE  = "E:/TORUS/LAST/output/thetas/thetas (0.230000).txt.4768";         // TODO сделать нормально!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    const QString OUTPUT_FILE = "E:/TORUS/LAST/output/psis/psi_coeffs.txt";    // TODO сделать нормально!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DataReaderWriterPoincare dataProvider(INPUT_FILE, OUTPUT_FILE);

    const auto getData = [&dataProvider](double& time_out, std::vector<double>& thetas_out) -> bool
    {
        return dataProvider.readLineInto(time_out, thetas_out);
    };
    const auto storeData = [&dataProvider](const double& time, const std::vector<double>& psi) -> void
    {
        dataProvider.writeLine(time, psi);
    };

    const int max_i = static_cast<int>(std::lround(std::sqrt(settings_.eq_count_)));
    const int max_j = max_i;

    const auto eq1 = *RHSCreator::createEquation1(settings_.a_, settings_.b_, *RHSCreator::createBasis(max_i, max_j));

    qDebug() << "Stream function coefficients computing start!";

    bool canProceed = true;    // значение, чтобы войти в цикл.
    int step = 1;
    while (canProceed)         // пока мы можем продолжать выполнение шагов.
    {
        QElapsedTimer timer;
        timer.start();

        canProceed = SF_step_computeCoefficients(getData, storeData, eq1);

        if (step % 500 == 0)
            emit stepInfo(step, timer.elapsed());

        step++;
    }

    qDebug() << "Stream function coefficients computing finish!";
}

void RK4Engine::onComputeStreamFunction()
{
    using tech::PoincareSection::DataReaderWriterPoincare;
    using tech::DataTextWriter;
    using math::NaturalConvection2D::Particle;
    using math::NaturalConvection2D::RHSCreator;

    printSettings(settings_);

    const QString INPUT_FILE  = "E:/TORUS/LAST/output/psis/psi_coeffs.6768.txt";    // TODO сделать нормально!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    const QString STUB_FILE   = "E:/TORUS/LAST/output/stream_function/stub.txt";
    const QString OUTPUT_FILE = "E:/TORUS/LAST/output/stream_function/stream_function.txt";        // TODO сделать нормально!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DataReaderWriterPoincare dataProvider(INPUT_FILE, STUB_FILE);
    DataTextWriter dataWriter(tech::helper::makeUniqueFileName(OUTPUT_FILE), DataTextWriter::OpenMode::kAppend);

    const auto getPsis = [&dataProvider](double& time_out, std::vector<double>& psis_out) -> bool
    {
        return dataProvider.readLineInto(time_out, psis_out);
    };

    const auto storeStreamFunction = [&dataWriter](const double &time,
                                                   const std::vector<Particle> &points,
                                                   const std::vector<double> &streamFunction) -> void
    {
        if (points.size() != streamFunction.size())
            throw std::runtime_error("<RK4Engine::onComputeStreamFunction()>: dimension missmatch!");

        const auto writeAction = [&time, &points, &streamFunction]() -> QString
        {
            QString result = QString("%1").arg(time);

            for (size_t index = 0; index < points.size(); index++)
                result += QString(" %1 %2 %3")
                    .arg(points.at(index).x_).arg(points.at(index).y_).arg(streamFunction.at(index));

            return result;
        };

        dataWriter.writeData(writeAction);
    };



    std::vector<Particle> PARTICLES;    // TODO временно добавил сетку для расчета функции тока!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for (int i = 1; i <= 69; i++)
        for (int j = 1; j <= 19; j++)
        {
            Particle prt;
            prt.x_ = 1.0 * j;
            prt.y_ = 1.0 * i;
            PARTICLES.push_back(prt);
        }


    for (const auto & part : PARTICLES)
        qDebug() << part.x_ << part.y_;


    qDebug() << "Stream function computing start!";

    bool canProceed = true;    // значение, чтобы войти в цикл.
    int step = 0;
    while (canProceed)         // пока мы можем продолжать выполнение шагов.
    {
        QElapsedTimer timer;
        timer.start();

        const int max_i = static_cast<int>(std::lround(std::sqrt(settings_.eq_count_)));
        const int max_j = max_i;
        canProceed = SF_step_computeStreamFunction(settings_.a_,
                                                   settings_.b_,
                                                   PARTICLES, //settings_.particles_,
                                                   *RHSCreator::createBasis(max_i, max_j),
                                                   getPsis,
                                                   storeStreamFunction);

        if (step % 500 == 0)
            emit stepInfo(step, timer.elapsed());

        step++;
    }

    qDebug() << "Stream function computing finish!";
}

// :::::::::::::::::::::::::::::::::::::::: Функция тока. ::::::::::::::::::::::::::::::::::::::::

bool RK4Engine::SF_step_computeCoefficients(const std::function<bool(double&, std::vector<double>&)> &getThetas,
                                            const std::function<void(const double&, const std::vector<double>&)> &storePsis,
                                            const std::vector<std::vector<math::NaturalConvection2D::Equation1Term>> &eq1)
{
    double currentTime;
    std::vector<double> currentThetas;
    if (not getThetas(currentTime, currentThetas))    // если невозможно получить текущие THETA.
        return false;                                 // то мы не можем продолжить выполнение шагов.

    const std::vector<double> currentPsis = SF_computeEquation1(currentThetas, eq1);    // поэлементное копирование производиться не будет!
    storePsis(currentTime, currentPsis);                                                // сохраняем посчитанные коэффициенты разложения функции тока в ряд.

    return true;     // мы успешно посчитали коэффициенты разложения функции тока - попробуем повторить это на следующем шаге.
}

auto RK4Engine::SF_computeEquation1(const std::vector<double> &thetas,
                                    const std::vector<std::vector<math::NaturalConvection2D::Equation1Term>> &eq1)
    -> std::vector<double>
{
    using math::NaturalConvection2D::Equation1Term;

    if (thetas.size() != eq1.size())
        throw std::runtime_error("<RK4Engine::SF_applyEquation1()>: dimension missmatch!");    // проекций и коэффициентов theta столько, сколько базисных функций!

    const auto thetaIndexToValue = [&thetas](const int thetaIndex) -> double
    {
        return thetas.at(thetaIndex);
    };

    std::vector<double> psis(eq1.size());                           // набор коэффициентов psi_k в разложении функции тока в ряд.
    std::transform(eq1.cbegin(), eq1.cend(), psis.begin(),
    [&thetaIndexToValue](const std::vector<Equation1Term> &psiEq) -> double
    {
        const double psi = SF_computePsiEquation(psiEq, thetaIndexToValue);
        return psi;
    });

    return psis;    // копирования не будет вследствие "copy elision" либо наличия у вектора конструктора перемещения!
}

double RK4Engine::SF_computePsiEquation(const std::vector<math::NaturalConvection2D::Equation1Term> &psiEq,
                                        const std::function<double(const int)> &thetaIndexToValue)
{
    using math::NaturalConvection2D::Equation1Term;

    double psi = 0;
    for (const Equation1Term &term : psiEq)
        psi += SF_computePsiTerm(term, thetaIndexToValue);

    return psi;
}

double RK4Engine::SF_computePsiTerm(const math::NaturalConvection2D::Equation1Term &term,
                                    const std::function<double(const int)> &thetaIndexToValue)
{
    return term.coeff_ * thetaIndexToValue(term.theta_index_);
}

bool RK4Engine::SF_step_computeStreamFunction(const double &a,
                                              const double &b,
                                              const std::vector<math::NaturalConvection2D::Particle> &points,
                                              const std::vector<math::NaturalConvection2D::BasisFunction> &basis,
                                              const std::function<bool(double&, std::vector<double>&)> &getPsis,
                                              const std::function<void(const double&,
                                                                       const std::vector<math::NaturalConvection2D::Particle>&,
                                                                       const std::vector<double>&)> &storeStreamFunction)
{
    using math::NaturalConvection2D::Particle;

    double currentTime;
    std::vector<double> currentPsis;
    if (not getPsis(currentTime, currentPsis))    // если невозможно получить текущие коэффициенты PSI.
        return false;                             // то мы не можем продолжить выполнение шагов.

    std::vector<double> streamFunction(points.size());
    std::transform(points.cbegin(), points.cend(), streamFunction.begin(),
    [&a, &b, &currentPsis, &basis](const Particle &point) -> double
    {
        return SF_computeStreamFunctionEquation(a, b, currentPsis, basis, point);
    });

    storeStreamFunction(currentTime, points, streamFunction);

    return true;     // мы успешно посчитали функцию тока в заданных точках - попробуем повторить это на следующем шаге.
}

double RK4Engine::SF_computeStreamFunctionEquation(const double &a,
                                                   const double &b,
                                                   const std::vector<double> &psis,
                                                   const std::vector<math::NaturalConvection2D::BasisFunction> &basis,
                                                   const math::NaturalConvection2D::Particle &point)
{
    if (psis.size() != basis.size())
        throw std::runtime_error("<RK4Engine::SF_computeStreamFunction()>: dimension missmatch!");

    double accum = 0;
    for (size_t index = 0; index < psis.size(); index++)
        accum += SF_computeStreamFunctionTerm(a, b, psis.at(index), basis.at(index), point);    // у меня нет std::transform_reduce()!

    return accum;
}

double RK4Engine::SF_computeStreamFunctionTerm(const double &a,
                                               const double &b,
                                               const double &psiCoeff,
                                               const math::NaturalConvection2D::BasisFunction &basisFunction,
                                               const math::NaturalConvection2D::Particle &point)
{
    const int ik = basisFunction.i_;
    const int jk = basisFunction.j_;

    const double x = point.x_;
    const double y = point.y_;

    const double basisFuncVal = 2.0 / std::sqrt(a*b) * std::sin(ik * M_PI*x/a) * std::sin(jk * M_PI*y/b);
    const double result       = psiCoeff * basisFuncVal;

    return result;
}

// :::::::::::::::::::::::::::::::::::::::: Сечение Пуанкаре. ::::::::::::::::::::::::::::::::::::::::

bool RK4Engine::PN_step(const std::function<bool(double&,
                                                 std::vector<double>&,
                                                 std::vector<math::NaturalConvection2D::Particle>&)> &getData,
                        const std::function<void(const double&,
                                                 const std::vector<double>&,
                                                 const std::vector<math::NaturalConvection2D::Particle>&)> &storeData,
                        const double &Rayleigh,
                        const math::PoincareSection::PlaneEquation &planeEquation,
                        const math::NaturalConvection2D::SystemRHS &rhs)
{
    using math::NaturalConvection2D::Particle;

    const auto isOnPlane = [&planeEquation](const std::vector<double> &thetas) -> bool
    {
        const double ZERO_EPSILON = 1e-4;    // TODO подставить нормальное значение!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        return std::abs(PN_applyPlane(thetas, planeEquation)) <= ZERO_EPSILON;
    };

    const auto isAtOppositeSides = [&planeEquation](const std::vector<double> &currentThetas,
                                                    const std::vector<double> &nextThetas) -> bool
    {
        const double currentP = PN_applyPlane(currentThetas, planeEquation);
        const double nextP    = PN_applyPlane(nextThetas, planeEquation);
        return  currentP * nextP < 0;
    };

    double currentTime;
    std::vector<double> currentThetas;
    std::vector<Particle> currentParticles;
    if (not getData(currentTime, currentThetas, currentParticles))    // если невозможно получить текущие THETA и частицы.
        return false;                                                 // то мы не можем продолжить выполнение шагов.

    if (isOnPlane(currentThetas))    // если текущие THETA лежат на секущей плоскости.
    {
        storeData(currentTime, currentThetas, currentParticles);      // сохраняем текущие THETA и частицы как искомые.
        return true;                                                  // сообщаем, что можно переходить на следующий шаг.
    }

    double nextTime;
    std::vector<double> nextThetas;
    std::vector<Particle> nextParticles;
    if (not getData(nextTime, nextThetas, nextParticles))     // пытаемся получить следующие THETA и частицы.
        return false;                                         // если это невозможно, то мы не можем продолжить выполнение шагов.

    if (isOnPlane(nextThetas))    // если следующие THETA лежат на секущей плоскости.
    {
        storeData(nextTime, nextThetas, nextParticles);       // сохраняем следующие THETA и частицы как искомые.
        return true;                                          // сообщаем, что можно переходить на следующий шаг.
    }

    if (isAtOppositeSides(currentThetas, nextThetas))         // если текущие и следующие THETA находятся с РАЗНЫХ сторон секущей плоскти.
    {
        const std::optional<double> middleTime = PN_tryComputeIntersectionTime(Rayleigh, currentTime, nextTime, nextThetas, planeEquation, rhs.diffs_);

        if (not middleTime.has_value())                       // если нам не удалось вычислить время THETA, лежащих на секущей плоскости.
            return true;                                      // то нам ничего не остается, кроме как перейти на следующий шаг.

        const double h = middleTime.value() - currentTime;    // получаем новую длину шага для метода Рунге-Кутты.

        std::vector<std::vector<double>> K_BUFFERS(4, std::vector<double>(rhs.diffs_.size()));    // TODO выделить заранее!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        std::vector<double> middleThetas(rhs.diffs_.size());
        RK4Engine::RK4_computeDiff(rhs.diffs_, Rayleigh, h, currentThetas, middleThetas, K_BUFFERS);                        // делаем один шаг метода Рунге-Кутты для THETA.
        std::vector<Particle> middleParticles(currentParticles.size());
        RK4Engine::RK4_computeCoord(rhs.eqX_, rhs.eqY_, h, currentThetas, K_BUFFERS, currentParticles, middleParticles);    // делаем один шаг метода Рунге-Кутты для частиц.

        if (isOnPlane(middleThetas))    // если новые THETA лежат на секущей плоскости.
        {
            storeData(middleTime.value(), middleThetas, middleParticles);    // сохраняем их и соответствующие им частицы как искомые.

            qDebug() << QString("time = %1 in [%2; %3]: thetas are on plane.")
                .arg(middleTime.value()).arg(currentTime).arg(nextTime);
        }
        else
        {
            qDebug() << QString("thetas with time=%1 don't lay on plane!").arg(middleTime.value());
        }
    }

    return true;    // в любом случае мы успешно считывали данные на текущем шаге, поэтому попробуем перейти на следующий.
}

auto RK4Engine::PN_tryComputeIntersectionTime(const double &Rayleigh,
                                              const double &t_prev,
                                              const double &t_next,
                                              const std::vector<double> &thetas_next,
                                              const math::PoincareSection::PlaneEquation &planeEquation,
                                              const std::vector<std::vector<math::NaturalConvection2D::EquationTerm>> &diffRHS)
    -> std::optional<double>
{
    const auto isInsideTimeInterval = [&t_prev, &t_next](const double &t_mid) -> bool
    {
        return (t_prev < t_mid) && (t_mid < t_next);
    };

    const double t_mid = PN_newtonStep(Rayleigh, t_next, thetas_next, planeEquation, diffRHS);

    if (isInsideTimeInterval(t_mid))
    {
        return t_mid;
    }
    else
    {
        qDebug() << QString("%1 NOT in ( %2 ; %3 ).").arg(t_mid).arg(t_prev).arg(t_next);
        return std::nullopt;
    }
}

double RK4Engine::PN_newtonStep(const double &Rayleigh,
                                const double &time,
                                const std::vector<double> &thetas,
                                const math::PoincareSection::PlaneEquation &planeEquation,
                                const std::vector<std::vector<math::NaturalConvection2D::EquationTerm>> &diffRHS)
{
    const double P  = PN_applyPlane(thetas, planeEquation);
    const double dP = PN_applyPlaneDerivative(Rayleigh, thetas, planeEquation, diffRHS);

    return time - P/dP;
}

double RK4Engine::PN_applyPlane(const std::vector<double> &thetas,
                                const math::PoincareSection::PlaneEquation &planeEquation)
{
    if (thetas.size() != planeEquation.coeffs_.size())
        throw std::runtime_error("<RK4Engine>: dimension missmath!");

    const double initial = planeEquation.coeffFree_;

    return std::inner_product(thetas.cbegin(), thetas.cend(), planeEquation.coeffs_.cbegin(), initial);
}

double RK4Engine::PN_applyPlaneDerivative(const double &Rayleigh,
                                          const std::vector<double> &thetas,
                                          const math::PoincareSection::PlaneEquation &planeEquation,
                                          const std::vector<std::vector<math::NaturalConvection2D::EquationTerm>> &diffRHS)
{
    if (thetas.size() != planeEquation.coeffs_.size() || thetas.size() != diffRHS.size())
        throw std::runtime_error("<RK4Engine>: dimension missmath!");

    const auto thetaIndexToArgument = [&thetas](const int thetaIndex) -> double
    {
        return thetas.at(thetaIndex);
    };

    std::vector<double> DATA_BUFFER(thetas.size());    // TODO нужно выделить место заранее!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#pragma omp parallel for default(shared)
    for (size_t eqIndex = 0; eqIndex < DATA_BUFFER.size(); eqIndex++)
    {
        DATA_BUFFER.at(eqIndex) = RK4Engine::RK4_computeDiffEquation(diffRHS.at(eqIndex), Rayleigh, thetaIndexToArgument);
    }

    const double initial = 0;    // у производной уравнения плоскости нету свободного коэффициента!

    return std::inner_product(DATA_BUFFER.cbegin(), DATA_BUFFER.cend(), planeEquation.coeffs_.cbegin(), initial);
}

auto RK4Engine::PN_createPlane(const size_t dimension)
    -> std::unique_ptr<math::PoincareSection::PlaneEquation>
{
    using math::PoincareSection::PlaneEquation;

    const auto GENERATOR = [](const size_t dimension) -> std::unique_ptr<PlaneEquation>
    {
        if (dimension == 0)
            throw std::invalid_argument("<RK4Engine>: Poincare plane section equation is empty!");

        const double SHIFT            = -32.0;
        const size_t THETA_AXIS_INDEX = 1;

        auto p_plane = std::make_unique<PlaneEquation>();

        p_plane->coeffs_.assign(dimension, 0.0);
        p_plane->coeffs_.at(THETA_AXIS_INDEX) = 1.0;

        p_plane->coeffFree_ = SHIFT;

        return p_plane;
    };

    return GENERATOR(dimension);
}

// :::::::::::::::::::::::::::::::::::::::: Рунге-Кутта. ::::::::::::::::::::::::::::::::::::::::

void RK4Engine::RK4_step(const math::NaturalConvection2D::SystemRHS &rhs,
                         const double &Rayleigh,
                         const double &h,
                         const int step,
                         std::vector<std::vector<double>> &KThetaBuffers,
                         math::NaturalConvection2D::SystemData &data_in_out)
{
    const int prevStep = step - 1;

    data_in_out.time_.at(step) = data_in_out.time_.at(prevStep) + h;                           // считаем значение времени на текущем шаге.

    RK4_computeDiff(rhs.diffs_, Rayleigh, h, data_in_out.diffs_.at(prevStep),
        data_in_out.diffs_.at(step), KThetaBuffers);                                           // вычисляем все THETA для текущего шага.

    RK4_computeCoord(rhs.eqX_, rhs.eqY_, h, data_in_out.diffs_.at(prevStep), KThetaBuffers,
        data_in_out.coords_.at(prevStep), data_in_out.coords_.at(step));                       // вычисляем координаты всех точек для текущего шага.
}

void RK4Engine::RK4_computeDiff(const std::vector<std::vector<math::NaturalConvection2D::EquationTerm>> &diffRHS,
                                const double &Rayleigh,
                                const double &h,
                                const std::vector<double> &dataPrevStep,
                                std::vector<double> &dataCurStep_out,
                                std::vector<std::vector<double>> &KThetaVectors_out)
{
    if (diffRHS.size() != dataPrevStep.size() || dataCurStep_out.size() != dataPrevStep.size())
        throw std::runtime_error("<RK4Engine>: dimension missmath!");

    if (KThetaVectors_out.size() != 4)
        throw std::runtime_error("<RK4Engine>: wrong Runge-Kutta coefficient count!");

    std::vector<double> &K1 = KThetaVectors_out.at(0);
    RK4_computeDiffKVector(diffRHS, Rayleigh, dataPrevStep, std::nullopt, std::vector<double>(), K1);

    std::vector<double> &K2 = KThetaVectors_out.at(1);
    RK4_computeDiffKVector(diffRHS, Rayleigh, dataPrevStep, h/2, K1, K2);

    std::vector<double> &K3 = KThetaVectors_out.at(2);
    RK4_computeDiffKVector(diffRHS, Rayleigh, dataPrevStep, h/2, K2, K3);

    std::vector<double> &K4 = KThetaVectors_out.at(3);
    RK4_computeDiffKVector(diffRHS, Rayleigh, dataPrevStep, h, K3, K4);

#pragma omp parallel for default(shared)
    for (size_t i = 0; i < dataCurStep_out.size(); i++)
    {
        dataCurStep_out.at(i) = RK4_formula(dataPrevStep.at(i), h, K1.at(i), K2.at(i), K3.at(i), K4.at(i));
    }
}

void RK4Engine::RK4_computeDiffKVector(const std::vector<std::vector<math::NaturalConvection2D::EquationTerm>> &diffRHS,
                                       const double &Rayleigh,
                                       const std::vector<double> &dataPrevStep,
                                       const std::optional<double> &factorK,
                                       const std::vector<double> &KVector_in,
                                       std::vector<double> &KVector_out)
{
    if (diffRHS.size() != dataPrevStep.size() || KVector_out.size() != dataPrevStep.size())
        throw std::runtime_error("<RK4Engine>: dimension missmath!");

    if (factorK.has_value() && (KVector_in.size() != dataPrevStep.size()))
        throw std::runtime_error("<RK4Engine>: dimension missmath!");

    const auto thetaIndexToArgument = [&dataPrevStep, &factorK, &KVector_in](const int thetaIndex) -> double
    {
        const double addend = factorK.has_value()
            ? (factorK.value() * KVector_in.at(thetaIndex))
            : 0.0;

        return dataPrevStep.at(thetaIndex) + addend;
    };

#pragma omp parallel for default(shared)
    for (size_t eqIndex = 0; eqIndex <diffRHS.size(); eqIndex++)
    {
        KVector_out.at(eqIndex) = RK4_computeDiffEquation(diffRHS.at(eqIndex), Rayleigh, thetaIndexToArgument);
    }
}

double RK4Engine::RK4_computeDiffEquation(const std::vector<math::NaturalConvection2D::EquationTerm> &diffEquation,
                                          const double &Rayleigh,
                                          const std::function<double(const int)> &thetaIndexToArgument)
{
    using math::NaturalConvection2D::EquationTerm;

    if (diffEquation.empty())
        throw std::runtime_error("<RK4Engine>: empty right-hand side of the differential equation given!");

    double accumulator = 0;

    for (const EquationTerm &term : diffEquation)
        accumulator += RK4_computeDiffTerm(term, Rayleigh, thetaIndexToArgument);

    return accumulator;
}

double RK4Engine::RK4_computeDiffTerm(const math::NaturalConvection2D::EquationTerm &term,
                                      const double &Rayleigh,
                                      const std::function<double(const int)> &thetaIndexToArgument)
{
    double accumulator = term.coeff_;                       // берем числовой коэффициент перед произведением заданного количества THETA.

    if (term.has_Rayleigh_)
        accumulator *= Rayleigh;                            // домножаем на число Рэлея при его наличии в данном члене-слагаемом.

    for (const int thetaIndex : term.theta_indices_)
        accumulator *= thetaIndexToArgument(thetaIndex);    // вычисляем произведение заданного количества THETA.

    return accumulator;                                     // возвращаем посчитанное числовое значение члена-слагаемого.
}

void RK4Engine::RK4_computeCoord(const std::vector<math::NaturalConvection2D::CoordinateTerm> &equationX,
                                 const std::vector<math::NaturalConvection2D::CoordinateTerm> &equationY,
                                 const double &h,
                                 const std::vector<double> &thetasPrevStep,
                                 const std::vector<std::vector<double>> &KThetaVectors,
                                 const std::vector<math::NaturalConvection2D::Particle> &particlesPrevStep,
                                 std::vector<math::NaturalConvection2D::Particle> &particlesCurStep_out)
{
    using math::NaturalConvection2D::Particle;

    if (particlesPrevStep.size() != particlesCurStep_out.size())
        throw std::runtime_error("<RK4Engine>: dimension missmath!");

    const std::vector<double> &K1ThetaVector = KThetaVectors.at(0);
    const std::vector<double> &K2ThetaVector = KThetaVectors.at(1);
    const std::vector<double> &K3ThetaVector = KThetaVectors.at(2);

#pragma omp parallel for default(shared)
    for (size_t index = 0; index < particlesCurStep_out.size(); index++)
    {
        const Particle &prtcl_in = particlesPrevStep.at(index);
        Particle &prtcl_out      = particlesCurStep_out.at(index);

        KCoordCoeff K1;
        RK4_computeCoordK(equationX, equationY, thetasPrevStep, prtcl_in, std::nullopt, std::vector<double>(), {0, 0}, K1);

        KCoordCoeff K2;
        RK4_computeCoordK(equationX, equationY, thetasPrevStep, prtcl_in, h/2, K1ThetaVector, K1, K2);

        KCoordCoeff K3;
        RK4_computeCoordK(equationX, equationY, thetasPrevStep, prtcl_in, h/2, K2ThetaVector, K2, K3);

        KCoordCoeff K4;
        RK4_computeCoordK(equationX, equationY, thetasPrevStep, prtcl_in, h, K3ThetaVector, K3, K4);

        prtcl_out.x_ = RK4_formula(prtcl_in.x_, h, K1.x_, K2.x_, K3.x_, K4.x_);
        prtcl_out.y_ = RK4_formula(prtcl_in.y_, h, K1.y_, K2.y_, K3.y_, K4.y_);
    }
}

void RK4Engine::RK4_computeCoordK(const std::vector<math::NaturalConvection2D::CoordinateTerm> &equationX,
                                  const std::vector<math::NaturalConvection2D::CoordinateTerm> &equationY,
                                  const std::vector<double> &thetasPrevStep,
                                  const math::NaturalConvection2D::Particle &particlePrevStep,
                                  const std::optional<double> &factorK,
                                  const std::vector<double> &KThetaVector_in,
                                  const KCoordCoeff &K_in,
                                  KCoordCoeff &K_out)
{
    if (factorK.has_value() && (thetasPrevStep.size() != KThetaVector_in.size()))
        throw std::runtime_error("<RK4Engine>: dimension missmath!");

    const auto thetaIndexToArgument = [&thetasPrevStep, &factorK, &KThetaVector_in](const int thetaIndex) -> double
    {
        const double addend = factorK.has_value()
            ? (factorK.value() * KThetaVector_in.at(thetaIndex))
            : 0.0;

        return thetasPrevStep.at(thetaIndex) + addend;
    };

    const auto wrap = [&factorK](const double &var, const double &KCoeff) -> double
    {
        const double addend = factorK.has_value()
            ? (factorK.value() * KCoeff)
            : 0.0;

        return var + addend;
    };

    K_out.x_ = RK4_computeCoordEquation(equationX, wrap(particlePrevStep.x_, K_in.x_), wrap(particlePrevStep.y_, K_in.y_), thetaIndexToArgument);
    K_out.y_ = RK4_computeCoordEquation(equationY, wrap(particlePrevStep.y_, K_in.y_), wrap(particlePrevStep.x_, K_in.x_), thetaIndexToArgument);
}

double RK4Engine::RK4_computeCoordEquation(const std::vector<math::NaturalConvection2D::CoordinateTerm> &coordinateEquation,
                                           const double &arg_sin, const double &arg_cos,
                                           const std::function<double(const int)> &thetaIndexToArgument)
{
    using math::NaturalConvection2D::CoordinateTerm;

    if (coordinateEquation.empty())
        throw std::runtime_error("<RK4Engine>: empty right-hand side of the coordinate equation given!");

    double accumulator = 0;

    for (const CoordinateTerm &term : coordinateEquation)
        accumulator += RK4_computeCoordTerm(term, arg_sin, arg_cos, thetaIndexToArgument);

    return accumulator;
}

double RK4Engine::RK4_computeCoordTerm(const math::NaturalConvection2D::CoordinateTerm &term,
                                       const double &arg_sin, const double &arg_cos,
                                       const std::function<double(const int)> &thetaIndexToArgument)
{
    return term.coeff_
        * std::sin(term.coeff_sin_*arg_sin)
        * std::cos(term.coeff_cos_*arg_cos)
        * thetaIndexToArgument(term.theta_index_);
}

double RK4Engine::RK4_formula(const double &scalarPrevStep, const double &h, const double &K1, const double &K2, const double &K3, const double &K4)
{
    return scalarPrevStep + h/6 * (K1 + 2*K2 + 2*K3 + K4);
};

// :::::::::::::::::::::::::::::::::::::::: Полезности. ::::::::::::::::::::::::::::::::::::::::

auto RK4Engine::createData(const RK4Engine::Settings &settings)
    -> std::unique_ptr<math::NaturalConvection2D::SystemData>
{
    using math::NaturalConvection2D::SystemData;
    using math::NaturalConvection2D::Particle;

    std::unique_ptr<math::NaturalConvection2D::SystemData> data = std::make_unique<SystemData>(settings.step_count_+1, settings.eq_count_, settings.particles_.size());



    /*
    std::random_device seed;
    std::mt19937 engine(seed());
    std::normal_distribution<double> gauss(3.0, 2.0);
    for (double &initialTheta : data->diffs_.front())
        initialTheta = gauss(engine);
    */
    if (data->diffs_.front().size() != INITIALS.size())
        throw std::runtime_error("INVALID INITIALS!");
    else
        data->diffs_.front() = INITIALS;



    std::transform(settings.particles_.cbegin(), settings.particles_.cend(), data->coords_.front().begin(),
    [](const Particle &particle) -> Particle
    {
        return particle;
    });

    data->time_.front() = settings.t_beg_;

    return data;
}

auto RK4Engine::create_K_coefficients(const int eq_count)
    -> std::unique_ptr<std::vector<std::vector<double>>>
{
    auto result = std::make_unique<std::vector<std::vector<double>>>(4);

    for (std::vector<double> &vec : *result)
        vec.resize(eq_count);

    return result;
}

void RK4Engine::writeData(const math::NaturalConvection2D::SystemData &data,  const size_t endIndex) const
{
    using tech::helper::pair_time_thetas;
    using tech::helper::pair_time_particles;
    using math::NaturalConvection2D::Particle;

    if (data.time_.size() != data.diffs_.size() ||
        data.time_.size() != data.coords_.size())
        throw std::runtime_error("<RK4Engine::writeData()>: data dimensions missmatch!");

    const auto addRayleigh = [rayleigh = settings_.rayleigh_](const QString &oldBase) -> QString
    {
        return QString("%1 (%2)").arg(oldBase, QString::number(rayleigh, 'f'));
    };

    std::vector<pair_time_thetas> time_thetas;
    time_thetas.reserve(endIndex);
    std::transform(data.time_.cbegin(), std::next(data.time_.cbegin(), endIndex), data.diffs_.cbegin(), std::back_inserter(time_thetas),
    [](const double &time, const std::vector<double> &thetas) -> pair_time_thetas
    {
        return pair_time_thetas(std::cref(time), std::cref(thetas));
    });
    tech::DataWriter(tech::helper::makeUniqueFileName(settings_.file_path_thetas_, addRayleigh))
        .writeData<pair_time_thetas>(time_thetas.cbegin(), time_thetas.cend(), tech::helper::serializeThetas);

    std::vector<pair_time_particles> time_particles;
    time_particles.reserve(endIndex);
    std::transform(data.time_.cbegin(), std::next(data.time_.cbegin(), endIndex), data.coords_.cbegin(), std::back_inserter(time_particles),
    [](const double &time, const std::vector<Particle> &particles) -> pair_time_particles
    {
        return pair_time_particles(std::cref(time), std::cref(particles));
    });
    tech::DataWriter(tech::helper::makeUniqueFileName(settings_.file_path_particle_, addRayleigh))
        .writeData<pair_time_particles>(time_particles.cbegin(), time_particles.cend(), tech::helper::serializeParticles);
}

void RK4Engine::cycleData(math::NaturalConvection2D::SystemData &data)
{
    const auto cycle = [](auto &data_to, const auto &data_from) -> void
    {
        if (data_to.size() != data_from.size())
            throw std::runtime_error("<RK4Engine>: incompatible data elements count!");

        for (size_t index = 0; index < data_to.size(); index++)
            data_to.at(index) = data_from.at(index);
    };

    cycle(data.diffs_.front(), data.diffs_.back());
    cycle(data.coords_.front(), data.coords_.back());
    data.time_.front() = data.time_.back();
}

void RK4Engine::maybeEmitInfo(const int step, const int elapsed_ms, const std::vector<math::NaturalConvection2D::Particle> &currentCoordinates)
{
    const auto mayDisplay = [](const int step) -> bool
    {
        const int STEP_TO_DISPLAY = 500;        // отображаем каждый пятисотый шаг.
        return step % STEP_TO_DISPLAY == 0;
    };

    if (mayDisplay(step))
    {
        emit stepInfo(step, elapsed_ms);
        emit particlesCoordinates(currentCoordinates);
    }
}
