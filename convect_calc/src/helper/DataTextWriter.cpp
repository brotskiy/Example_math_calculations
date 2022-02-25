#include "DataTextWriter.h"

#include <QDir>
#include <QFileInfo>
#include <QIODevice>

#include <exception>
#include <optional>
#include <utility>

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

tech::DataTextWriter::DataTextWriter(const QString &filePath, const tech::DataTextWriter::OpenMode &openMode)
{
    if (filePath.isEmpty())
        throw std::runtime_error("<DataTextWriter::DataTextWriter()>: file path is empty!");

    const auto toQtOpenMode = [](const OpenMode &openMode) -> std::optional<QIODevice::OpenMode>
    {
        switch (openMode)
        {
            case OpenMode::kNormal:
                return std::nullopt;
            case OpenMode::kTruncate:
                return QIODevice::Truncate;
            case OpenMode::kAppend:
                return QIODevice::Append;
            default:
                throw std::runtime_error("<DataTextWriter::DataTextWriter()>: file open mode is invalid!");
        }
    };

    QDir().mkpath(QFileInfo(filePath).absolutePath());
    p_file_.reset(new QFile(filePath));

    const QIODevice::OpenMode fileOpenModes = toQtOpenMode(openMode).has_value()
        ? QIODevice::WriteOnly | QIODevice::Text | toQtOpenMode(openMode).value()
        : QIODevice::WriteOnly | QIODevice::Text;

    if (not p_file_->open(fileOpenModes))
        throw std::runtime_error("<DataTextWriter::DataTextWriter()>: file can not be opened!");

    p_stream_.reset(new QTextStream(p_file_.get()));
}

tech::DataTextWriter::~DataTextWriter()
{
    p_file_->close();
}

tech::DataTextWriter::DataTextWriter(tech::DataTextWriter &&other) :
    p_file_(std::move(other.p_file_)),
    p_stream_(std::move(other.p_stream_))
{
}

tech::DataTextWriter& tech::DataTextWriter::operator=(tech::DataTextWriter &&other)
{
    p_file_   = std::move(other.p_file_);
    p_stream_ = std::move(other.p_stream_);

    return *this;
}

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

void tech::DataTextWriter::writeData(const std::function<QString()> &writeAction)
{
    *p_stream_ << writeAction() << Qt::endl;
}
