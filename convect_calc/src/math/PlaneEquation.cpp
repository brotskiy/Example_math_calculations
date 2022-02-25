#include "PlaneEquation.h"

#include <stdexcept>
#include <numeric>
#include <functional>

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

math::PlaneEquation::PlaneEquation(const std::vector<double> &NORMAL, const double &D) :
    NORMAL_(NORMAL),
    D_(D)
{
    if (NORMAL_.size() <= 1)
        throw std::runtime_error("<PlaneEquation::PlaneEquation()>: invalid dimension!");
}

math::PlaneEquation::PlaneEquation(const std::vector<double> &NORMAL, const std::vector<double> &point) :
    NORMAL_(NORMAL)
{
    if (NORMAL_.size() <= 1)
        throw std::runtime_error("<PlaneEquation::PlaneEquation()>: invalid dimension!");

    if (NORMAL_.size() != point.size())
        throw std::runtime_error("<PlaneEquation::PlaneEquation()>: dimension missmatch!");

    D_ = std::inner_product(NORMAL_.cbegin(), NORMAL_.cend(), point.cbegin(), 0.0, std::minus<>(), std::multiplies<>());
}

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

auto math::PlaneEquation::getNormal() const -> std::vector<double>
{
    return NORMAL_;
}

double math::PlaneEquation::substitute(const std::vector<double> &point) const
{
    if (NORMAL_.size() != point.size())
        throw std::runtime_error("<PlaneEquation::PlaneEquation()>: dimension missmatch!");

    return std::inner_product(NORMAL_.cbegin(), NORMAL_.cend(), point.cbegin(), D_);
}
