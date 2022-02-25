#include "EngineConfiger.h"

namespace
{
    const QString BASICS_    = "BASICS";
    const QString PARTICLES_ = "PARTICLES";
}

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

tech::EngineConfiger::EngineConfiger(const QString& configFilePath) :
    tech::ConfigerBase(configFilePath)
{
}

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

void tech::EngineConfiger::writeEngineConfig(const RK4Engine::Settings& engnSttngs)
{
    using math::NaturalConvection2D::Particle;

    const auto write = [&engnSttngs](QSettings& config) -> void
    {
        config.beginGroup(BASICS_);

        config.setValue("a_",                  engnSttngs.a_);
        config.setValue("b_",                  engnSttngs.b_);
        config.setValue("t_beg_",              engnSttngs.t_beg_);
        config.setValue("t_end_",              engnSttngs.t_end_);
        config.setValue("rayleigh_",           engnSttngs.rayleigh_);
        config.setValue("eq_count_",           engnSttngs.eq_count_);
        config.setValue("step_count_",         engnSttngs.step_count_);
        config.setValue("repeat_count_",       engnSttngs.repeat_count_);
        config.setValue("file_path_thetas_",   engnSttngs.file_path_thetas_);
        config.setValue("file_path_particle_", engnSttngs.file_path_particle_);

        const int particleCount = static_cast<int>(engnSttngs.particles_.size());

        config.beginWriteArray(PARTICLES_, particleCount);
        for (int index = 0; index < particleCount; index++)
        {
            const Particle &particle = engnSttngs.particles_.at(index);

            config.setArrayIndex(index);
            config.setValue("x_", particle.x_);
            config.setValue("y_", particle.y_);
        }
        config.endArray();

        config.endGroup();
    };

    writeWrapper(write);
}

bool tech::EngineConfiger::readEngineConfig(RK4Engine::Settings &engnSttngs_out)
{
    using math::NaturalConvection2D::Particle;

    const auto read = [&engnSttngs_out](QSettings& config) -> bool
    {
        config.beginGroup(BASICS_);

        engnSttngs_out.a_                  = config.value("a_",                  engnSttngs_out.a_).toDouble();
        engnSttngs_out.b_                  = config.value("b_",                  engnSttngs_out.b_).toDouble();
        engnSttngs_out.t_beg_              = config.value("t_beg_",              engnSttngs_out.t_beg_).toDouble();
        engnSttngs_out.t_end_              = config.value("t_end_",              engnSttngs_out.t_end_).toDouble();
        engnSttngs_out.rayleigh_           = config.value("rayleigh_",           engnSttngs_out.rayleigh_).toDouble();
        engnSttngs_out.eq_count_           = config.value("eq_count_",           engnSttngs_out.eq_count_).toInt();
        engnSttngs_out.step_count_         = config.value("step_count_",         engnSttngs_out.step_count_).toInt();
        engnSttngs_out.repeat_count_       = config.value("repeat_count_",       engnSttngs_out.repeat_count_).toInt();
        engnSttngs_out.file_path_thetas_   = config.value("file_path_thetas_",   engnSttngs_out.file_path_thetas_).toString();
        engnSttngs_out.file_path_particle_ = config.value("file_path_particle_", engnSttngs_out.file_path_particle_).toString();

        const int particleCount = config.beginReadArray(PARTICLES_);
        engnSttngs_out.particles_.resize(static_cast<size_t>(particleCount));
        for (int index = 0; index < particleCount; index++)
        {
            Particle &particle = engnSttngs_out.particles_.at(index);

            config.setArrayIndex(index);
            particle.x_ = config.value("x_", particle.x_).toDouble();
            particle.y_ = config.value("y_", particle.y_).toDouble();
        }
        config.endArray();

        config.endGroup();

        return true;
    };

    return readWrapper(read);
}
