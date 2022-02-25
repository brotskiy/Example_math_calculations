#include "DataWriter.h"

#include <QDir>

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

tech::DataWriter::DataWriter(const QString &dataFilePath)
    : dataFilePath_(dataFilePath)
{
    if (dataFilePath_.isEmpty())
        throw std::invalid_argument("<DataWriter>: empty data file path!");

    QDir().mkpath(QFileInfo(dataFilePath_).absolutePath());
}
