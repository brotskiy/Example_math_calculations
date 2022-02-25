#include "miscellanea.h"

#include <QFileInfo>
#include <QCoreApplication>
#include <QStringList>
#include <QRegularExpression>

#include <cstring>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <exception>

QString tech::helper::makeUniqueFileName(QString filePath)
{
    filePath.append("." + QString::number(QCoreApplication::applicationPid()));

    return filePath;
}

QString tech::helper::makeUniqueFileName(QString filePath, const std::function<QString(const QString&)> &transformBaseName)
{
    const QString oldBase = QFileInfo(filePath).baseName();
    const QString newBase = transformBaseName(oldBase);

    const int from = filePath.lastIndexOf(oldBase);
    const int size = oldBase.size();
    filePath.replace(from, size, newBase);

    return makeUniqueFileName(filePath);
}

QString tech::helper::serializeThetas(const tech::helper::pair_time_thetas &time_thetas)
{
    QString line = QString::number(time_thetas.first.get(), 'e', 20);

    for (const double &theta : time_thetas.second.get())
        line.append(" " + QString::number(theta, 'e', 20));

    return line;
};

QString tech::helper::serializeParticles(const tech::helper::pair_time_particles &time_particles)
{
    using math::NaturalConvection2D::Particle;

    QString line = QString::number(time_particles.first.get(), 'e', 20);

    for (const Particle& prt : time_particles.second.get())
        line.append(" " + QString::number(prt.x_, 'e', 20) + " " + QString::number(prt.y_, 'e', 20));

    return line;
};

bool tech::helper::isSame(const double &first, const double &second)
{
    return 0 == std::memcmp(&first, &second, sizeof(double));
}

void tech::helper::readEngineSettings(QSettings &file, RK4ParticleEngine::Settings &settings_out)
{
    using math::NaturalConvection2D::Particle;

    const QString grp_BASICS    = "BASICS";
    const QString grp_PARTICLES = "PARTICLES";

    file.beginGroup(grp_BASICS);
    settings_out.a_            = file.value("a_",            settings_out.a_).toDouble();
    settings_out.b_            = file.value("b_",            settings_out.b_).toDouble();
    settings_out.t_beg_        = file.value("t_beg_",        settings_out.t_beg_).toDouble();
    settings_out.step_         = file.value("step_",        settings_out.step_).toDouble();
    settings_out.step_count_   = file.value("step_count_",   settings_out.step_count_).toDouble();
    settings_out.repeat_count_ = file.value("repeat_count_", settings_out.repeat_count_).toDouble();
    settings_out.eq_count_     = file.value("eq_count_",     settings_out.eq_count_).toDouble();
    file.endGroup();

    const int particleCount = file.beginReadArray(grp_PARTICLES);
    settings_out.initial_particles_.resize(particleCount);

    for (int index = 0; index < particleCount; index++)
    {
        Particle &particle = settings_out.initial_particles_.at(index);

        file.setArrayIndex(index);
        particle.x_ = file.value("x_", particle.x_).toDouble();
        particle.y_ = file.value("y_", particle.y_).toDouble();
    }

    file.endArray();
}

void tech::helper::writeEngineSettings(QSettings &file, const RK4ParticleEngine::Settings &settings)
{
    using math::NaturalConvection2D::Particle;

    const QString grp_BASICS    = "BASICS";
    const QString grp_PARTICLES = "PARTICLES";

    file.beginGroup(grp_BASICS);
    file.setValue("a_",            settings.a_);
    file.setValue("b_",            settings.b_);
    file.setValue("t_beg_",        settings.t_beg_);
    file.setValue("step_",         settings.step_);
    file.setValue("step_count_",   settings.step_count_);
    file.setValue("repeat_count_", settings.repeat_count_);
    file.setValue("eq_count_",     settings.eq_count_);
    file.endGroup();

    const int particleCount = settings.initial_particles_.size();
    file.beginWriteArray(grp_PARTICLES, particleCount);

    for (int index = 0; index < particleCount; index++)
    {
        const Particle &particle = settings.initial_particles_.at(index);

        file.setArrayIndex(index);
        file.setValue("x_", particle.x_);
        file.setValue("y_", particle.y_);
    }

    file.endArray();
}

void tech::helper::stringToPairTimeValues(const QString &str, std::pair<double, std::vector<double>> &time_values_out)
{
    const QStringList list_time_values = str.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);

    if (list_time_values.size() <= 1)
        throw std::runtime_error("<stringToPairTimeValues()>: invalid line has been read!");

    time_values_out.first = list_time_values.front().toDouble();
    time_values_out.second.resize(list_time_values.size() - 1);

    std::transform(std::next(list_time_values.cbegin()), list_time_values.cend(), time_values_out.second.begin(),
    [](const QString& value_string) -> double
    {
        return value_string.toDouble();
    });
}

QString tech::helper::timeParticlesToString(const double &time, const std::vector<math::NaturalConvection2D::Particle> &particles)
{
    using math::NaturalConvection2D::Particle;

    QString str = QString::number(time, 'e', 20);

    for (const Particle& prt : particles)
        str.append(" " + QString::number(prt.x_, 'e', 20) + " " + QString::number(prt.y_, 'e', 20));

    return str;
}

double math::normalizeValue(const double &T0, const double &periodT, const double &t)
{
    if (not (periodT > 0))
        throw std::runtime_error("<math::normalizeValue()>: period T must be greater than zero!");

    double distance = std::fmod(t - T0, periodT);
    if (distance < 0)
        distance += periodT;                // находим расстояние по модулю евклида между t и начальным значением T0.

    const double result = T0 + distance;    // сносим t в заданный промежуток.

    return result;
}

double math::interpLin(const double &x0, const double &y0,
                       const double &x1, const double &y1,
                       const double &x)
{
    using tech::helper::isSame;

    if (isSame(x0, x1))
        throw std::runtime_error("<math::linInterp()>: x0 and x1 are binary equal!");

    if (x0 > x1)
        throw std::runtime_error("<math::linInterp()>: x0 is greater than x1!");

    if (x < x0 || x1 < x)
        throw std::runtime_error("<math::linInterp()>: x is out of the interval [x0; x1]!");

    if (isSame(x, x0))    // если попали точно в левую границу.
        return y0;

    if (isSame(x, x1))    // если попали точно в правую границу.
        return y1;

    return y0 + (y1 - y0)/(x1 - x0) * (x - x0);
}

void math::interpLin(const double &t0, const std::vector<double> &v0,
                     const double &t1, const std::vector<double> &v1,
                     const double &t, std::vector<double> &v_out)
{
    if (v0.size() != v1.size() || v0.size() != v_out.size())
        throw std::runtime_error("<math::linInterpVec()>: dimension missmatch!");

    for (size_t i = 0; i < v_out.size(); i++)
        v_out.at(i) = interpLin(t0, v0.at(i), t1, v1.at(i), t);
}
