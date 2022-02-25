#include <QApplication>

#include "interface/mainwindow.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w(nullptr);
    w.show();

    return a.exec();
}
