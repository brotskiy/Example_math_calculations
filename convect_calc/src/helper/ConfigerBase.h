#ifndef CONFIGERBASE_H
#define CONFIGERBASE_H

#include <QString>
#include <QSettings>
#include <QMutex>
#include <QMutexLocker>

#include <memory>
#include <functional>

namespace tech
{

class ConfigerBase
{

// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

private:
    static const QSettings::Format FORMAT_ = QSettings::IniFormat;

    const QString configFilePath_;
    std::unique_ptr<QMutex> p_readWriteMutex_ = nullptr;

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

public:
    explicit ConfigerBase(const QString& configFilePath);
    virtual ~ConfigerBase() = 0;                            // абстрактный класс!

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public:
    void writeWrapper(std::function<void(QSettings&)> write);
    bool readWrapper(std::function<bool(QSettings&)> read);

// :::::::::::::::::::::::::::::::::::::::: Полезности. ::::::::::::::::::::::::::::::::::::::::

public:
    static bool isInvalidPath(const QString& filePath);
    static void makeDir(const QString& filePath);

// :::::::::::::::::::::::::::::::::::::::: Ограничения. ::::::::::::::::::::::::::::::::::::::::

private:
    ConfigerBase()                               = delete;
    ConfigerBase(const ConfigerBase&)            = delete;
    ConfigerBase& operator=(const ConfigerBase&) = delete;
    ConfigerBase(ConfigerBase&&)                 = delete;
    ConfigerBase& operator=(ConfigerBase&&)      = delete;
};

}

#endif // CONFIGERBASE_H
