#ifndef CONFIGER_H
#define CONFIGER_H

#include <QSettings>
#include <QString>
#include <QMutex>
#include <QMutexLocker>
#include <QDir>
#include <QFileInfo>

#include <memory>
#include <functional>
#include <exception>

namespace tech
{

template<class SettingsType>
class Configer
{
// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

private:
    static const QSettings::Format FORMAT_ = QSettings::IniFormat;

    const QString configFilePath_;
    const std::function<void(QSettings&, SettingsType&)> readFromIni_;
    const std::function<void(QSettings&, const SettingsType&)> writeIntoIni_;

    std::unique_ptr<QMutex> p_readWriteMutex_ = nullptr;

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

public:
    Configer(const QString& configFilePath,
             const std::function<void(QSettings&, SettingsType&)> &readFromIni,
             const std::function<void(QSettings&, const SettingsType&)> &writeIntoIni);
    ~Configer() = default;

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public:
    void readSettings(SettingsType &settings_out);
    void writeSettings(const SettingsType &settings);

// :::::::::::::::::::::::::::::::::::::::: Ограничения. ::::::::::::::::::::::::::::::::::::::::

private:
    Configer()                           = delete;
    Configer(const Configer&)            = delete;
    Configer& operator=(const Configer&) = delete;
    Configer(Configer&&)                 = delete;
    Configer& operator=(Configer&&)      = delete;
};

}



template<class SettingsType>
tech::Configer<SettingsType>::Configer(const QString &configFilePath,
                                       const std::function<void(QSettings&, SettingsType&)> &readFromIni,
                                       const std::function<void(QSettings&, const SettingsType&)> &writeIntoIni) :
    configFilePath_(configFilePath),
    readFromIni_(readFromIni),
    writeIntoIni_(writeIntoIni)
{
    if (configFilePath_.isEmpty())
        throw std::runtime_error("<Configer::Configer()>: config file path is empty!");

    QDir().mkpath(QFileInfo(configFilePath).absolutePath());
    p_readWriteMutex_.reset( new QMutex() );
}

template<class SettingsType>
void tech::Configer<SettingsType>::readSettings(SettingsType &settings_out)
{
    QMutexLocker locker(p_readWriteMutex_.get());
    QSettings settingsFile(configFilePath_, FORMAT_);

    readFromIni_(settingsFile, settings_out);
}

template<class SettingsType>
void tech::Configer<SettingsType>::writeSettings(const SettingsType &settings)
{
    QMutexLocker locker(p_readWriteMutex_.get());
    QSettings settingsFile(configFilePath_, FORMAT_);

    writeIntoIni_(settingsFile, settings);
}

#endif // CONFIGER_H
