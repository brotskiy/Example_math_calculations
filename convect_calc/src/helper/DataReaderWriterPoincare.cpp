#include "DataReaderWriterPoincare.h"

#include <QFileInfo>
#include <QDir>
#include <QCoreApplication>
#include <QStringList>
#include <QRegularExpression>

#include <exception>
#include <algorithm>
#include <iterator>

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

tech::PoincareSection::DataReaderWriterPoincare::DataReaderWriterPoincare(const QString &inputThetasPath, const QString &outputThetasPath)
{
    if (inputThetasPath.isEmpty() ||
        not QFileInfo::exists(inputThetasPath) ||
        not QFileInfo(inputThetasPath).isFile())
        throw std::runtime_error("<DataReaderWriterPoincare>: invalid input thetas file!");

    p_file_input_.reset( new QFile(inputThetasPath) );
    if (not p_file_input_->open(QIODevice::ReadOnly | QIODevice::Text))
        throw std::runtime_error("<DataReaderWriterPoincare>: input thetas file can not be opened!");

    p_stream_input_.reset( new QTextStream(p_file_input_.get()) );

    const auto insertPid = [](QString str) -> QString
    {
        const QString baseName    = QFileInfo(str).baseName();
        const QString newBaseName = QString("%1%2%3").arg(baseName, ".", QString::number(QCoreApplication::applicationPid()));
        str.replace(str.lastIndexOf(baseName), baseName.size(), newBaseName);

        return str;
    };

    QDir().mkpath(QFileInfo(outputThetasPath).absolutePath());

    p_file_output_.reset( new QFile(insertPid(outputThetasPath)) );
    if (not p_file_output_->open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate))
        throw std::runtime_error("<DataReaderWriterPoincare>: output thetas file can not be opened!");

    p_stream_output_.reset( new QTextStream(p_file_output_.get()) );
}

tech::PoincareSection::DataReaderWriterPoincare::~DataReaderWriterPoincare()
{
    p_file_input_->close();
    p_file_output_->close();
}

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

bool tech::PoincareSection::DataReaderWriterPoincare::readLineInto(double &time_out, std::vector<double> &thetas_out)
{
    QString str;
    if (not p_stream_input_->readLineInto(&str))
        return false;

    const QStringList time_thetas = str.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);

    if (time_thetas.size() <= 1)
        throw std::runtime_error("<DataReaderWriterPoincare>: invalid line has been read!");

    time_out = time_thetas.front().toDouble();

    thetas_out.resize(time_thetas.size() - 1);
    std::transform(std::next(time_thetas.cbegin()), time_thetas.cend(), thetas_out.begin(),
    [](const QString &theta) -> double
    {
        return theta.toDouble();
    });

    return true;
}

void tech::PoincareSection::DataReaderWriterPoincare::writeLine(const double &time, const std::vector<double> &thetas)
{
    *p_stream_output_ << time << " ";

    for (const double &theta : thetas)
        *p_stream_output_ << theta << " ";

    *p_stream_output_ << Qt::endl;
}

bool tech::PoincareSection::DataReaderWriterPoincare::readLineInto(double &time_out, std::vector<math::NaturalConvection2D::Particle> &particles_out)
{
    using math::NaturalConvection2D::Particle;

    QString str;
    if (not p_stream_input_->readLineInto(&str))
        return false;

    const QStringList time_particles = str.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);

    if (time_particles.size() <= 1)
        throw std::runtime_error("<DataReaderWriterPoincare>: invalid line has been read!");

    time_out = time_particles.front().toDouble();

    particles_out.clear();
    auto pIt = time_particles.cbegin();
    while (pIt != time_particles.cend())
    {
        const auto itx = std::next(pIt, 1);
        if (itx == time_particles.cend())
            return true;

        const auto ity = std::next(pIt, 2);
        if (ity == time_particles.cend())
            throw std::runtime_error("<DataReaderWriterPoincare>: invalid particles line has been read!");

        Particle prt;
        prt.x_ = QString(*itx).toDouble();
        prt.y_ = QString(*ity).toDouble();
        particles_out.push_back(prt);

        pIt = ity;
    }

    return true;
}


void tech::PoincareSection::DataReaderWriterPoincare::writeLine(const double &time, const std::vector<math::NaturalConvection2D::Particle> &particles)
{
    using math::NaturalConvection2D::Particle;

    *p_stream_output_ << time;

    for (const Particle &prt : particles)
        *p_stream_output_ << " " << prt.x_ << " " << prt.y_;

    *p_stream_output_ << Qt::endl;
}
