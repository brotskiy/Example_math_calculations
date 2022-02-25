#ifndef RK4PARTICLEENGINE_H
#define RK4PARTICLEENGINE_H

#include <vector>
#include <utility>
#include <functional>

#include "engine/NaturalConvection2DStructs.h"

class RK4ParticleEngine
{
// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

public:
    struct Settings
    {
        // T = 146,61531329952732144

        double a_, b_;
        double t_beg_;
        double step_;
        int step_count_;
        int repeat_count_;
        int eq_count_;

        std::vector<std::pair<double, std::vector<double>>> times_thetas_;

        std::vector<math::NaturalConvection2D::Particle> initial_particles_;

        Settings();
    };
    static void printSettings(const RK4ParticleEngine::Settings &settings);

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public:
    static void simulateParticles(const RK4ParticleEngine::Settings &settings,
                                  const std::function<void(const double&,
                                                           const std::vector<math::NaturalConvection2D::Particle>&)> &storeParticles,
                                  const std::function<void(const int&,
                                                           const double&,
                                                           const double&,
                                                           const std::vector<math::NaturalConvection2D::Particle>&)> &onStepCallback);

private:
    static void simulateParticles_step(const std::function<void(const double&, const std::vector<double>&, std::vector<double>&)> &system,
                                       const double &h,
                                       const std::pair<double, std::vector<math::NaturalConvection2D::Particle>> &timeParticlesPrevStep,
                                       std::pair<double, std::vector<math::NaturalConvection2D::Particle>> &timeParticlesCurStep_out);

// :::::::::::::::::::::::::::::::::::::::: Полезности. ::::::::::::::::::::::::::::::::::::::::

private:
    static auto createDiffSystemForCoordinates(const std::vector<std::pair<double, std::vector<double>>> &times_thetas,
                                               const int &basisSize,
                                               const double &a,
                                               const double &b)
        -> std::function<void(const double&, const std::vector<double>&, std::vector<double>&)>;

    static void computeThetasAtTime(const std::vector<std::pair<double, std::vector<double>>> &times_thetas,
                                    const double &time,
                                    std::vector<double> &thetas_out);

    static void prepareParticleData(const int &timestampCount,
                                    const double &initialTime,
                                    const std::vector<math::NaturalConvection2D::Particle> &initialParticles,
                                    std::vector<std::pair<double, std::vector<math::NaturalConvection2D::Particle>>> &particlesAtTimes_out);

    static void storeParticleData(const std::vector<std::pair<double, std::vector<math::NaturalConvection2D::Particle>>> &particlesAtTimes,
                                  const std::function<void(const double&, const std::vector<math::NaturalConvection2D::Particle>&)> &storeParticles);
    static void cycleParticleData(std::vector<std::pair<double, std::vector<math::NaturalConvection2D::Particle>>> &particlesAtTimes_in_out);

};

#endif // RK4PARTICLEENGINE_H
