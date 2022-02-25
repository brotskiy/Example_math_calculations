#ifndef PLANEEQUATION_H
#define PLANEEQUATION_H

#include <vector>

namespace math
{

class PlaneEquation
{
// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

private:
    std::vector<double> NORMAL_;
    double D_;

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

public:
    PlaneEquation(const std::vector<double> &NORMAL, const double &D);
    PlaneEquation(const std::vector<double> &NORMAL, const std::vector<double> &point);

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public:
    auto getNormal() const -> std::vector<double>;
    double substitute(const std::vector<double> &point) const;
};

}

#endif // PLANEEQUATION_H
