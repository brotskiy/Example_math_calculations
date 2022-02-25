#ifndef DATAWRITER_H
#define DATAWRITER_H

#include <QString>
#include <QFile>
#include <QTextStream>
#include <QDebug>

#include <functional>
#include <exception>
#include <type_traits>

namespace tech
{

class DataWriter
{
// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

private:
    const QString dataFilePath_;    // TODO заменить на QFile и QTextStream!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

public:
    DataWriter(const QString &dataFilePath);

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public:
    template<class DataType, class InputIt>
    void writeData(const InputIt &it_beg, const InputIt &it_end,
                   const std::function<QString(const DataType&)> &serialize) const;
};

}



template<class DataType, class InputIt>
void tech::DataWriter::writeData(const InputIt &it_beg, const InputIt &it_end,
                                 const std::function<QString(const DataType&)> &serialize) const
{
    static_assert(std::is_same_v<std::decay_t<DataType>, std::decay_t<decltype(*it_beg)>>,
        "<DataWriter::writeData()>: iterator must point to element of the specified type!");

    if (it_beg == it_end)
        throw std::runtime_error("<DataWriter::writeData()>: can't open data file!");

    QFile dataFile(dataFilePath_);

    const bool isOpened = dataFile.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Text);
    if (not isOpened)
        throw std::runtime_error("<DataWriter::writeData()>: can't open data file!");

    QTextStream stream(&dataFile);

    for (InputIt it = it_beg; it != it_end; it++)
        stream << serialize(*it) << Qt::endl;

    dataFile.close();
}

#endif // DATAWRITER_H
