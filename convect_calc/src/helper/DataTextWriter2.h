#ifndef DATATEXTWRITER2_H
#define DATATEXTWRITER2_H

#include <QString>
#include <QFile>
#include <QTextStream>
#include <QDir>
#include <QFileInfo>
#include <QIODevice>

#include <functional>
#include <memory>

namespace tech
{

template<class DatumType, class... OtherDatumTypes>
class DataTextWriter2
{
// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

public:
    enum class OpenMode : int
    {
        kNormal,      // просто открываем существующий файл сначала или создаем новый.
        kTruncate,    // открываем существующий файл, очищая его содержимое, или создаем новый.
        kAppend       // открываем существующий файл на дозапись или создаем новый.
    };

private:
    const std::function<QString(const DatumType&, const OtherDatumTypes&...)> datumToString_;

    std::unique_ptr<QFile> p_file_         = nullptr;
    std::unique_ptr<QTextStream> p_stream_ = nullptr;

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

public:
    DataTextWriter2(const QString &filePath,
                    const OpenMode &openMode,
                    const std::function<QString(const DatumType&, const OtherDatumTypes&...)> &datumToString);
    ~DataTextWriter2();

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public:
    void writeLine(const DatumType& datum, const OtherDatumTypes&... otherData);

// :::::::::::::::::::::::::::::::::::::::: Служебное. ::::::::::::::::::::::::::::::::::::::::

private:
    static QIODevice::OpenMode toQtOpenMode(const OpenMode &openMode);

// :::::::::::::::::::::::::::::::::::::::: Ограничения. ::::::::::::::::::::::::::::::::::::::::

private:
    DataTextWriter2()                                  = delete;
    DataTextWriter2(const DataTextWriter2&)            = delete;
    DataTextWriter2& operator=(const DataTextWriter2&) = delete;
    DataTextWriter2(DataTextWriter2&&)                 = delete;
    DataTextWriter2& operator=(DataTextWriter2&&)      = delete;
};

}



// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

template<class DatumType, class... OtherDatumTypes>
tech::DataTextWriter2<DatumType, OtherDatumTypes...>::DataTextWriter2(const QString &filePath,
                                                                      const OpenMode &openMode,
                                                                      const std::function<QString(const DatumType&, const OtherDatumTypes&...)> &datumToString) :
    datumToString_(datumToString)
{
    if (filePath.isEmpty())
        throw std::runtime_error("<DataTextWriter2::DataTextWriter2()>: file path is empty!");

    QDir().mkpath(QFileInfo(filePath).absolutePath());    // в любом случае пытаемся создать папку для файла.
    p_file_.reset(new QFile(filePath));

    const QIODevice::OpenMode fileOpenModes = QIODevice::WriteOnly | QIODevice::Text | toQtOpenMode(openMode);

    if (not p_file_->open(fileOpenModes))
        throw std::runtime_error("<DataTextWriter2::DataTextWriter2()>: can't open file!");

    p_stream_.reset(new QTextStream(p_file_.get()));
}

template<class DatumType, class... OtherDatumTypes>
tech::DataTextWriter2<DatumType, OtherDatumTypes...>::~DataTextWriter2()
{
    p_file_->close();
}

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

template<class DatumType, class... OtherDatumTypes>
void tech::DataTextWriter2<DatumType, OtherDatumTypes...>::writeLine(const DatumType& datum, const OtherDatumTypes&... otherData)
{
    *p_stream_ << datumToString_(datum, otherData...) << Qt::endl;
}

// :::::::::::::::::::::::::::::::::::::::: Служебное. ::::::::::::::::::::::::::::::::::::::::

template<class DatumType, class... OtherDatumTypes>
QIODevice::OpenMode tech::DataTextWriter2<DatumType, OtherDatumTypes...>::toQtOpenMode(const OpenMode &openMode)
{
    switch (openMode)
    {
        case OpenMode::kNormal:
            return QIODevice::OpenMode();    // просто открываем существующий файл сначала или создаем новый.
        case OpenMode::kTruncate:
            return QIODevice::Truncate;      // открываем существующий файл, очищая его содержимое, или создаем новый.
        case OpenMode::kAppend:
            return QIODevice::Append;        // открываем существующий файл на дозапись или создаем новый.
        default:
            throw std::runtime_error("<DataTextWriter2::toQtOpenMode()>: file open mode is invalid!");
    }
}

#endif // DATATEXTWRITER2_H
