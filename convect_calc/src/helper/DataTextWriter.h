#ifndef DATATEXTWRITER_H
#define DATATEXTWRITER_H

#include <QString>
#include <QFile>
#include <QTextStream>

#include <memory>
#include <functional>

namespace tech
{

class DataTextWriter
{
// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

private:
    std::unique_ptr<QFile> p_file_         = nullptr;
    std::unique_ptr<QTextStream> p_stream_ = nullptr;

public:
    enum class OpenMode : int
    {
        kNormal,
        kTruncate,
        kAppend
    };

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

public:
    DataTextWriter(const QString &filePath, const tech::DataTextWriter::OpenMode &openMode);
    ~DataTextWriter();

    DataTextWriter(tech::DataTextWriter &&other);
    DataTextWriter& operator=(tech::DataTextWriter &&other);

    DataTextWriter(const tech::DataTextWriter &other)      = delete;
    DataTextWriter& operator=(tech::DataTextWriter &other) = delete;

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public:
    void writeData(const std::function<QString()> &writeAction);

};

}


#endif // DATATEXTWRITER_H
