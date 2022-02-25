#include "widget.h"

#include <QPainter>
#include <QPen>
#include <QRectF>
#include <QColor>
#include <QBrush>
#include <QPointF>

#include <random>

#include <cmath>
#include <algorithm>
#include <iterator>


// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

const int SceneWidget::HEIGHT_ = 750;
const int SceneWidget::WIDTH_  = 750;

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

SceneWidget::SceneWidget(const SceneWidget::Settings &settings, QWidget *parentWidget ) :
    QWidget(parentWidget),
    settings_(settings)
{
    std::generate_n(std::back_inserter(images_), settings_.particleCount_ + 1,
    []() -> std::unique_ptr<QImage>
    {
        return std::make_unique<QImage>(SceneWidget::WIDTH_, SceneWidget::HEIGHT_, QImage::Format_ARGB32);
    });

    p_currentImage_ = images_.front().get();

    const double scale = computeScale(settings_.a_, settings_.b_);

    for (std::unique_ptr<QImage> &p_image : images_)
    {
        p_image->fill(Qt::white);

        std::unique_ptr<QPainter> p_painter(createPainter(p_image.get(), settings_.a_, settings_.b_));
        p_painter->setPen(QPen(Qt::black, 1, Qt::SolidLine, Qt::FlatCap));
        p_painter->drawRect(QRectF(0.0, 0.0, scale*settings_.a_, scale*settings_.b_));
    }

    update();
}

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

void SceneWidget::paintEvent(QPaintEvent* event)
{
  QPainter(this).drawImage(0, 0, *p_currentImage_);

  QWidget::paintEvent(event);
}

void SceneWidget::onDrawParticles(const std::vector<math::NaturalConvection2D::Particle> &particles)
{
    using math::NaturalConvection2D::Particle;

    auto colour = [rndm = std::default_random_engine()]() mutable -> int    // mutable лямбда!
    {
        return 40 + rndm() % 180;
    };

    const auto pen = [&colour]() -> QPen
    {
        return QPen(QBrush(QColor(colour(), colour(), colour())), 3, Qt::SolidLine, Qt::RoundCap);
    };

    std::unique_ptr<QPainter> painterMaster(createPainter(images_.front().get(), settings_.a_, settings_.b_ ));

    for (size_t index = 0; index < particles.size(); index++)
    {
        const double scale = computeScale(settings_.a_, settings_.b_);

        std::unique_ptr<QPainter> painterParticle(createPainter(images_.at(index + 1).get(), settings_.a_, settings_.b_));
        painterMaster->setPen(pen());
        painterParticle->setPen(pen());

        const QPointF particle(scale*particles.at(index).x_, scale*particles.at(index).y_);

        painterMaster->drawPoint(particle);
        painterParticle->drawPoint(particle);
    }

    update();
}

void SceneWidget::onSetCurrentImage(const QString& newImage)
{
    if (newImage.contains("MASTER"))
        p_currentImage_ = images_.front().get();
    else
        p_currentImage_ = images_.at(newImage.toInt()).get();

    update();
}

// :::::::::::::::::::::::::::::::::::::::: Полезности. ::::::::::::::::::::::::::::::::::::::::

QPainter* SceneWidget::createPainter(QPaintDevice *p_paintDevice, double a, double b)
{
    const double scale = computeScale(a, b);

    QPainter* p_painter = new QPainter(p_paintDevice);
    p_painter->setRenderHint(QPainter::Antialiasing);
    p_painter->translate(SceneWidget::WIDTH_/2 - a*scale/2, SceneWidget::HEIGHT_/2 + b*scale/2);
    p_painter->scale(1.0, -1.0);

    return p_painter;
}

double SceneWidget::computeScale(double a, double b)
{
    return std::abs(a) > std::abs(b)
        ? (SceneWidget::WIDTH_ - 70.0) / a
        : (SceneWidget::HEIGHT_ - 70.0) / b;
}
