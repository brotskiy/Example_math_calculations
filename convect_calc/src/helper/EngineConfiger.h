#ifndef ENGINECONFIGER_H
#define ENGINECONFIGER_H

#include <QString>

#include "ConfigerBase.h"
#include "engine/mainengine.h"

namespace tech
{

class EngineConfiger : private tech::ConfigerBase
{
// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

public:
    explicit EngineConfiger(const QString& configFilePath);
    ~EngineConfiger() override = default;

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public:
    void writeEngineConfig(const RK4Engine::Settings &engnSttngs);
    bool readEngineConfig(RK4Engine::Settings &engnSttngs_out);

// :::::::::::::::::::::::::::::::::::::::: Ограничения. ::::::::::::::::::::::::::::::::::::::::

private:
    EngineConfiger()                                 = delete;
    EngineConfiger(const EngineConfiger&)            = delete;
    EngineConfiger(EngineConfiger&&)                 = delete;
    EngineConfiger& operator=(const EngineConfiger&) = delete;
    EngineConfiger& operator=(EngineConfiger&&)      = delete;
};

}

#endif // ENGINECONFIGER_H
