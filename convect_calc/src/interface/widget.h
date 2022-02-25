#ifndef WIDGET_H
#define WIDGET_H

#include <QObject>
#include <QString>
#include <QWidget>
#include <QImage>
#include <QPainter>
#include <QPaintDevice>
#include <QPaintEvent>

#include <vector>
#include <memory>

#include "engine/NaturalConvection2DStructs.h"

class SceneWidget : public QWidget
{
    Q_OBJECT

// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

public:
    struct Settings
    {
        QString outputDirectory_;

        int particleCount_;
        double a_, b_;
    };

private:
    static const int HEIGHT_;
    static const int WIDTH_;
    const Settings settings_;

    std::vector<std::unique_ptr<QImage>> images_;
    const QImage* p_currentImage_ = nullptr;

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

public:
    SceneWidget(const SceneWidget::Settings &settings, QWidget *parentWidget = nullptr);
    ~SceneWidget() override = default;

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

protected:
    void paintEvent(QPaintEvent* event) override;

public slots:
    void onDrawParticles(const std::vector<math::NaturalConvection2D::Particle> &particles);
    void onSetCurrentImage(const QString& newImage);

// :::::::::::::::::::::::::::::::::::::::: Полезности. ::::::::::::::::::::::::::::::::::::::::

private:
    static QPainter* createPainter(QPaintDevice *p_paintDevice, double a, double b);
    static double computeScale(double a, double b);

// :::::::::::::::::::::::::::::::::::::::: Ограничения. ::::::::::::::::::::::::::::::::::::::::

private:
    SceneWidget()                              = delete;
    SceneWidget(const SceneWidget&)            = delete;
    SceneWidget(SceneWidget&&)                 = delete;
    SceneWidget& operator=(const SceneWidget&) = delete;
    SceneWidget& operator=(SceneWidget&&)      = delete;
};

#endif // WIDGET_H
