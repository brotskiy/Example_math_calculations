QT += core gui widgets

MY_ROOT_DIRECTORY = $$PWD/..

#CONFIG (release, debug|release) {
    MY_TARGET_NAME   = convec
    #MY_TMP_DIRECTORY = &&{MY_ROOT_DIRECTORY}/tmp/release
#} else {
#    MY_TARGET_NAME   = convec_debug
#    MY_TMP_DIRECTORY = &&{MY_ROOT_DIRECTORY}/tmp/debug
#}

TEMPLATE    = app
TARGET      = $${MY_TARGET_NAME}
#DESTDIR     = $${MY_ROOT_DIRECTORY}/bin
#OBJECTS_DIR = $${MY_TMP_DIRECTORY}/objs
#MOC_DIR     = $${MY_TMP_DIRECTORY}/mocs
#RCC_DIR     = $${MY_TMP_DIRECTORY}/qrcs
#UI_DIR      = $${MY_TMP_DIRECTORY}/uis

CONFIG += c++17

INCLUDEPATH += $${MY_ROOT_DIRECTORY}/src

SOURCES += $$files($${MY_ROOT_DIRECTORY}/src/*.cpp)           \
           $$files($${MY_ROOT_DIRECTORY}/src/helper/*.cpp)    \
           $$files($${MY_ROOT_DIRECTORY}/src/engine/*.cpp)    \
           $$files($${MY_ROOT_DIRECTORY}/src/math/*.cpp)      \
           $$files($${MY_ROOT_DIRECTORY}/src/interface/*.cpp)

HEADERS += $$files($${MY_ROOT_DIRECTORY}/src/*.h)             \
           $$files($${MY_ROOT_DIRECTORY}/src/helper/*.h)      \
           $$files($${MY_ROOT_DIRECTORY}/src/engine/*.h)      \
           $$files($${MY_ROOT_DIRECTORY}/src/math/*.h)        \
           $$files($${MY_ROOT_DIRECTORY}/src/interface/*.h)
QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS += -fopenmp

LIBS += -lgomp -lpthread
