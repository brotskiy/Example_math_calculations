#include "RK4ParticleEngine.h"

#include <QElapsedTimer>

#include <cmath>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <exception>

#include "engine/RHSCreator.h"
#include "engine/RHSComputer.h"
#include "engine/RK4Computer.h"
#include "helper/miscellanea.h"

// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

RK4ParticleEngine::Settings::Settings()
{
    a_ = 20;
    b_ = 70;
    t_beg_ = 0;
    step_  = 0.01;
    step_count_ = 400000;
    repeat_count_ = 10;
    eq_count_ = 100;

    initial_particles_ = { {5, 10}, {10, 10}, {15, 10},
                           {5, 32}, {10, 32}, {15, 32},
                           {5, 63}, {10, 63}, {15, 63} };
}

void RK4ParticleEngine::printSettings(const RK4ParticleEngine::Settings &settings)
{
    std::cout << "a_ ="            << settings.a_            << std::endl
              << "b_ ="            << settings.b_            << std::endl
              << "t_beg_ ="        << settings.t_beg_        << std::endl
              << "step_ ="         << settings.step_         << std::endl
              << "step_count_ ="   << settings.step_count_   << std::endl
              << "repeat_count_ =" << settings.repeat_count_ << std::endl
              << "eq_count_ ="     << settings.eq_count_     << std::endl;
}

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

void RK4ParticleEngine::simulateParticles(const RK4ParticleEngine::Settings &settings,
                                          const std::function<void(const double&,
                                                                   const std::vector<math::NaturalConvection2D::Particle>&)> &storeParticles,
                                          const std::function<void(const int&,
                                                                   const double&,
                                                                   const double&,
                                                                   const std::vector<math::NaturalConvection2D::Particle>&)> &onStepCallback)
{
    using math::NaturalConvection2D::CoordinateTerm;
    using math::NaturalConvection2D::Particle;

    printSettings(settings);

    const auto system = createDiffSystemForCoordinates(settings.times_thetas_, settings.eq_count_, settings.a_, settings.b_);    // системы дифф. уравнений для координат частицы.

    std::vector<std::pair<double, std::vector<Particle>>> particlesAtTimes;                                           // выделяем заранее память под местоположения
    prepareParticleData(1 + settings.step_count_, settings.t_beg_, settings.initial_particles_, particlesAtTimes);    // всех частиц для всех шагов.

    // const double h = (settings.t_end_ - settings.t_beg_) / settings.step_count_;    // вычисляем постоянный шаг алгоритма по времени.
    const double h = settings.step_;

    std::cout << "Particles RK4 start!" << std::endl;

    for (int repeat = 0; repeat < settings.repeat_count_; repeat++)
    {
        for (int step = 1; step <= settings.step_count_; step++)
        {
            QElapsedTimer timer;
            timer.start();

            simulateParticles_step(system, h, particlesAtTimes.at(step-1), particlesAtTimes.at(step));

            const double elapsed  = timer.elapsed();
            const double realStep = repeat * settings.step_count_ + step;
            onStepCallback(realStep, elapsed, particlesAtTimes.at(step).first, particlesAtTimes.at(step).second);
        }

        storeParticleData(particlesAtTimes, storeParticles);
        cycleParticleData(particlesAtTimes);
    }

    std::cout << "Particles RK4 finish!" << std::endl;
}

void RK4ParticleEngine::simulateParticles_step(const std::function<void(const double&, const std::vector<double>&, std::vector<double>&)> &system,
                                               const double &h,
                                               const std::pair<double, std::vector<math::NaturalConvection2D::Particle>> &timeParticlesPrevStep,
                                               std::pair<double, std::vector<math::NaturalConvection2D::Particle>> &timeParticlesCurStep_out)
{
    using math::NaturalConvection2D::Particle;
    using math::RK4Computer;

    if (timeParticlesPrevStep.second.size() != timeParticlesCurStep_out.second.size())
        throw std::runtime_error("<RK4ParticleEngine::simulateParticles_step()>: dimension missmatch!");

    timeParticlesCurStep_out.first = timeParticlesPrevStep.first + h;                   // вычисляем новое значение времени!

    #pragma omp parallel for default(shared)
    for (size_t index = 0; index < timeParticlesCurStep_out.second.size(); index++)     // перебираем все частицы (параллельно).
    {
        const Particle &particlePrev = timeParticlesPrevStep.second.at(index);
        Particle &particleCur_out    = timeParticlesCurStep_out.second.at(index);

        const std::vector<double> input = { particlePrev.x_, particlePrev.y_ };
        std::vector<double> output(2);

        const double timeOfInput = timeParticlesPrevStep.first;

        RK4Computer().RK4(system, h, timeOfInput, input, output);    // вычисляем новое местоположение частицы методом Рунге-Кутты.

        particleCur_out.x_ = output.at(0);
        particleCur_out.y_ = output.at(1);
    }
}

// :::::::::::::::::::::::::::::::::::::::: Полезности. ::::::::::::::::::::::::::::::::::::::::

auto RK4ParticleEngine::createDiffSystemForCoordinates(const std::vector<std::pair<double, std::vector<double>>> &times_thetas,
                                                       const int &basisSize,
                                                       const double &a,
                                                       const double &b)
    -> std::function<void(const double&, const std::vector<double>&, std::vector<double>&)>
{
    using math::NaturalConvection2D::RHSCreator;
    using math::RHSComputer;

    const int max_i = static_cast<int>(std::lround(std::sqrt(basisSize)));
    const int max_j = max_i;

    const auto basis    = *RHSCreator::createBasis(max_i, max_j);
    const auto eq1      = *RHSCreator::createEquation1(a, b, basis);
    const auto coordEqs = *RHSCreator::createCoords(basis, eq1, a, b);

    // создаём обёртку для правой части дифференциальных уравнений координат частицы.

    return [eqX = coordEqs.first,
            eqY = coordEqs.second,
            &times_thetas](const double &time, const std::vector<double> &input, std::vector<double> &output) -> void
    {
        std::vector<double> thetasAtTime;
        computeThetasAtTime(times_thetas, time, thetasAtTime);                        // находим THETA в интересующий момент, используя заданные на периоде значения.

        const auto thetaIndexToValue = [&thetasAtTime](const int &index) -> double    // функция получения значения THETA по индексу.
        {
            return thetasAtTime.at(index);
        };

        const double &inputX = input.at(0);
        const double &inputY = input.at(1);

        output.at(0) = RHSComputer::computeCoordinateEquation(eqX, inputX, inputY, thetaIndexToValue);    // считаем новую X координату.
        output.at(1) = RHSComputer::computeCoordinateEquation(eqY, inputY, inputX, thetaIndexToValue);    // считаем новую Y координату.
    };
}

void RK4ParticleEngine::computeThetasAtTime(const std::vector<std::pair<double, std::vector<double>>> &times_thetas,
                                            const double &time,
                                            std::vector<double> &thetas_out)
{
    if (times_thetas.size() == 1)
        throw std::runtime_error("<RK4ParticleEngine::computeThetasAtTime()>: set of THETA contains only initial values!");

    const double T0   = times_thetas.front().first;    // время начала периода THETA.
    const double TMax = times_thetas.back().first;     // время КОНЦА периода THETA!

    if (tech::helper::isSame(T0, TMax))
        throw std::runtime_error("<RK4ParticleEngine::computeThetasAtTime()>: initial and final THETA time values are binary equal!");

    double normalizedTime = math::normalizeValue(T0, TMax - T0, time);    // сносим интересующее нас значение времени в заданные рамки периода!
    if (normalizedTime > TMax)
    {
        std::cout << "Rounding error! normT > TMax : " << normalizedTime << " > " << TMax << std::endl;

        normalizedTime = TMax;    // если вследствие ошибки округления мы перескочили за конец периода, то делаем значение снова валидным.
    }

    auto t1_It = std::upper_bound(times_thetas.cbegin(), times_thetas.cend(), normalizedTime,    // находим набор THETA с временем БОЛЬШИМ, чем интересующее нас.
    [](const double &time, const std::pair<double, std::vector<double>> &time_thetas) -> bool
    {
        return time < time_thetas.first;    // переданное значение времени СТРОГО МЕНЬШЕ временной отметки текущего набора THETA.
    });

    if (t1_It == times_thetas.cbegin())
        throw std::runtime_error("<RK4ParticleEngine::computeThetasAtTime()>: normalized time can't be less than T0 by the calculation method!");

    if (t1_It == times_thetas.cend())                // если нет ни одной отметки времени, значение которой больше, чем нормализованное время.
        t1_It = std::prev(times_thetas.cend());      // значит нормализованное время попало ровно в TMax, поэтому корректируем итератор.

    const auto t0_It = std::prev(t1_It);             // получаем набор THETA с предыдущим значением времени.

    const double &t0              = t0_It->first;    // ссылка-псевдоним для удобства!
    const std::vector<double> &v0 = t0_It->second;   // ссылка-псевдоним для удобства!

    const double &t1              = t1_It->first;    // ссылка-псевдоним для удобства!
    const std::vector<double> &v1 = t1_It->second;   // ссылка-псевдоним для удобства!

    thetas_out.resize(v0.size());
    math::interpLin(t0, v0, t1, v1, normalizedTime, thetas_out);
}

void RK4ParticleEngine::prepareParticleData(const int &timestampCount,
                                            const double &initialTime,
                                            const std::vector<math::NaturalConvection2D::Particle> &initialParticles,
                                            std::vector<std::pair<double, std::vector<math::NaturalConvection2D::Particle>>> &particlesAtTimes_out)
{
    if (timestampCount < 1)
        throw std::runtime_error("<RK4ParticleEngine::prepareParticleData()>: zero timestamp count!");

    if (initialParticles.empty())
        throw std::runtime_error("<RK4ParticleEngine::prepareParticleData()>: zero initial particles count!");

    particlesAtTimes_out.resize(timestampCount);

    for (auto &time_particles : particlesAtTimes_out)
        time_particles.second.resize(initialParticles.size());

    particlesAtTimes_out.front().first  = initialTime;
    particlesAtTimes_out.front().second = initialParticles;
}

void RK4ParticleEngine::storeParticleData(const std::vector<std::pair<double, std::vector<math::NaturalConvection2D::Particle>>> &particlesAtTimes,
                                          const std::function<void(const double&, const std::vector<math::NaturalConvection2D::Particle>&)> &storeParticles)
{
    const auto beg = particlesAtTimes.cbegin();
    const auto end = std::prev(particlesAtTimes.cend());    // конец зацикливается в начало, поэтому он не должен быть продублирован.

    for (auto it = beg; it != end; it++)
        storeParticles(it->first, it->second);
}

void RK4ParticleEngine::cycleParticleData(std::vector<std::pair<double, std::vector<math::NaturalConvection2D::Particle>>> &particlesAtTimes_in_out)
{
    particlesAtTimes_in_out.front() = particlesAtTimes_in_out.back();
}
