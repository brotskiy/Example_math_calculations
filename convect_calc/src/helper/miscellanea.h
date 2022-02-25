#ifndef MISCELLANEA_H
#define MISCELLANEA_H

#include <QString>
#include <QSettings>

#include <vector>
#include <utility>
#include <functional>

#include "engine/NaturalConvection2DStructs.h"
#include "engine/RK4ParticleEngine.h"

namespace tech
{
namespace helper
{

typedef std::pair<std::reference_wrapper<const double>, std::reference_wrapper<const std::vector<double>>> pair_time_thetas;
typedef std::pair<std::reference_wrapper<const double>, std::reference_wrapper<const std::vector<math::NaturalConvection2D::Particle>>> pair_time_particles;

QString makeUniqueFileName(QString filePath);
QString makeUniqueFileName(QString filePath, const std::function<QString(const QString&)> &transformBaseName);

QString serializeThetas(const tech::helper::pair_time_thetas &time_thetas);
QString serializeParticles(const tech::helper::pair_time_particles &time_particles);

bool isSame(const double &first, const double &second);

void readEngineSettings(QSettings &file, RK4ParticleEngine::Settings &settings_out);
void writeEngineSettings(QSettings &file, const RK4ParticleEngine::Settings &settings);

void stringToPairTimeValues(const QString &str, std::pair<double, std::vector<double>> &time_values_out);
QString timeParticlesToString(const double &time, const std::vector<math::NaturalConvection2D::Particle> &particles);

}
}

namespace math
{

double normalizeValue(const double &T0, const double &periodT, const double &t);

double interpLin(const double &x0, const double &y0,
                 const double &x1, const double &y1,
                 const double &x);

void interpLin(const double &t0, const std::vector<double> &v0,
               const double &t1, const std::vector<double> &v1,
               const double &t, std::vector<double> &v_out);
}

#endif // MISCELLANEA_H
