#include "ConfigerBase.h"

#include <QDir>
#include <QFileInfo>

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::


tech::ConfigerBase::ConfigerBase(const QString& configFilePath) :
    configFilePath_(configFilePath)
{
    p_readWriteMutex_.reset( new QMutex() );
    makeDir(configFilePath_);
}

tech::ConfigerBase::~ConfigerBase()    // деструктор абстрактного класса все равно должен быть определен!
{
}

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

void tech::ConfigerBase::writeWrapper(std::function<void(QSettings&)> write)
{
    if (isInvalidPath(configFilePath_))
        return;

    QMutexLocker locker(p_readWriteMutex_.get());
    makeDir(configFilePath_);
    QSettings config(configFilePath_, FORMAT_);

    write(config);
}

bool tech::ConfigerBase::readWrapper(std::function<bool(QSettings&)> read)
{
    if (isInvalidPath(configFilePath_))
        return false;

    QMutexLocker locker(p_readWriteMutex_.get());
    QSettings config(configFilePath_, FORMAT_);

    return read(config);
}

// :::::::::::::::::::::::::::::::::::::::: Полезности. ::::::::::::::::::::::::::::::::::::::::

bool tech::ConfigerBase::isInvalidPath(const QString& filePath)
{
    return filePath.isEmpty();
}

void tech::ConfigerBase::makeDir(const QString& filePath)
{
    if (isInvalidPath(filePath))
        return;

    QDir().mkpath(QFileInfo(filePath).absolutePath());
}
