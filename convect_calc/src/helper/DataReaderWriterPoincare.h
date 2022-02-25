#ifndef DATAREADERWRITERPOINCARE_H
#define DATAREADERWRITERPOINCARE_H

#include <QString>
#include <QTextStream>
#include <QFile>

#include <vector>
#include <memory>

#include "engine/NaturalConvection2DStructs.h"

namespace tech
{
namespace PoincareSection
{

class DataReaderWriterPoincare
{

// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

private:
    std::unique_ptr<QFile> p_file_input_, p_file_output_;
    std::unique_ptr<QTextStream> p_stream_input_, p_stream_output_;

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

public:
    DataReaderWriterPoincare(const QString &inputThetasPath, const QString &outputThetasPath);
    ~DataReaderWriterPoincare();

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public:
    bool readLineInto(double &time_out, std::vector<double> &thetas_out);
    void writeLine(const double &time, const std::vector<double> &thetas);

    bool readLineInto(double &time_out, std::vector<math::NaturalConvection2D::Particle> &particles_out);
    void writeLine(const double &time, const std::vector<math::NaturalConvection2D::Particle> &particles);

};

}
}

#endif // DATAREADERWRITERPOINCARE_H
