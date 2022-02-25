#ifndef DATATEXTREADER_H
#define DATATEXTREADER_H

#include <QString>
#include <QFile>
#include <QTextStream>
#include <QFileInfo>

#include <memory>
#include <functional>
#include <exception>

namespace tech
{

template<class DatumType, class... OtherDatumTypes>
class DataTextReader
{
// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

private:
    const std::function<void(const QString&, DatumType&, OtherDatumTypes&...)> stringToDatum_;

    std::unique_ptr<QFile> p_file_         = nullptr;
    std::unique_ptr<QTextStream> p_stream_ = nullptr;

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

public:
    DataTextReader(const QString &filePath,
                   const std::function<void(const QString&, DatumType&, OtherDatumTypes&...)> &stringToDatum);
    ~DataTextReader();

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public:
    bool readLineInto(DatumType& datum, OtherDatumTypes&... otherData);

// :::::::::::::::::::::::::::::::::::::::: Ограничения. ::::::::::::::::::::::::::::::::::::::::

private:
    DataTextReader()                                 = delete;
    DataTextReader(const DataTextReader&)            = delete;
    DataTextReader& operator=(const DataTextReader&) = delete;
    DataTextReader(DataTextReader&&)                 = delete;
    DataTextReader& operator=(DataTextReader&&)      = delete;
};

}



// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

template<class DatumType, class... OtherDatumTypes>
tech::DataTextReader<DatumType, OtherDatumTypes...>::DataTextReader(const QString &filePath,
                                                                    const std::function<void(const QString&, DatumType&, OtherDatumTypes&...)> &stringToDatum) :
    stringToDatum_(stringToDatum)
{
    if (filePath.isEmpty())
        throw std::runtime_error("<DataTextReader::DataTextReader()>: file path is empty!");

    const QFileInfo fileInfo(filePath);

    if (not (fileInfo.exists() && fileInfo.isFile()))
        throw std::runtime_error("<DataTextReader::DataTextReader()>: invalid file!");

    p_file_.reset( new QFile(filePath) );

    if (not p_file_->open(QIODevice::ReadOnly | QIODevice::Text))
        throw std::runtime_error("<DataTextReader::DataTextReader()>: can't open file!");

    p_stream_.reset( new QTextStream(p_file_.get()) );
}

template<class DatumType, class... OtherDatumTypes>
tech::DataTextReader<DatumType, OtherDatumTypes...>::~DataTextReader()
{
    p_file_->close();
}

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

template<class DatumType, class... OtherDatumTypes>
bool tech::DataTextReader<DatumType, OtherDatumTypes...>::readLineInto(DatumType& datum, OtherDatumTypes&... otherData)
{
    QString str;
    if (not p_stream_->readLineInto(&str))
        return false;

    stringToDatum_(str, datum, otherData...);

    return true;
}

#endif // DATATEXTREADER_H
