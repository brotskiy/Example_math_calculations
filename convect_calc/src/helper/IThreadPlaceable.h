#ifndef IWORKER_H
#define IWORKER_H

#include <QObject>
#include <QThread>
#include <QDebug>

#include <functional>
#include <type_traits>

namespace tech
{
namespace heavy
{

template<class Derived_t>
Derived_t* create(QObject* threadParent, std::function<Derived_t*()> factory, std::function<void()> onThreadStop);

template<class Derived_t>
void destroy(Derived_t* me);

}
}

template<class Derived_t>
Derived_t* tech::heavy::create(QObject* threadParent, std::function<Derived_t*()> factory, std::function<void()> onThreadStop)
{
    static_assert(std::is_base_of<QObject, Derived_t>::value, "Only QObject subclass can be placed inside a thread!");

    Derived_t*const worker = factory();
    if (worker)
    {
        QThread*const workerThread = new QThread(threadParent);
        worker->moveToThread(workerThread);

        QObject::connect(workerThread, &QThread::finished, static_cast<QObject*>(worker), onThreadStop);

        workerThread->start();
    }

    return worker;
}

template<class Derived_t>
void tech::heavy::destroy(Derived_t* me)
{
    static_assert(std::is_base_of<QObject, Derived_t>::value, "Only QObject subclass can be placed inside a thread!");

    if (me)
    {
        QThread*const workerThread = me->thread();
        if (workerThread == nullptr)
        {
            delete me;    // если мы по какой-то причине не можем получить поток обекта, то просто удаляем объект.
            return;
        }

        workerThread->quit();
        workerThread->wait();

        delete me;    // обязательно сначала нужно удалить сам объект, а только потом его поток.

        if (workerThread->parent() == nullptr)    // если потоку не был передан родитель, то мы сами отвечаем за его удаление.
        {
            delete workerThread;
            qDebug() << "Thread is destroyed manually!";
        }
    }
}

#endif // IWORKER_H
