#include "RHSComputer.h"

#include <cmath>
#include <exception>

// :::::::::::::::::::::::::::::::::::::::: Дифференциальное уравнение координаты частицы. ::::::::::::::::::::::::::::::::::::::::

double math::RHSComputer::computeCoordinateEquation(const std::vector<math::NaturalConvection2D::CoordinateTerm> &coordEq,
                                                    const double &var_sin,
                                                    const double &var_cos,
                                                    const std::function<double(const int&)> &thetaIndexToValue)
{
    using math::NaturalConvection2D::CoordinateTerm;

    if (coordEq.empty())
        throw std::runtime_error("<RHSComputer::computeCoordinateEquation()>: empty right-hand side of the coordinate equation!");

    double accumulator = 0;

    for (const CoordinateTerm &term : coordEq)
        accumulator += computeCoordinateTerm(term, var_sin, var_cos, thetaIndexToValue);

    return accumulator;
}

double math::RHSComputer::computeCoordinateTerm(const math::NaturalConvection2D::CoordinateTerm &term,
                                                const double &var_sin,
                                                const double &var_cos,
                                                const std::function<double(const int&)> &thetaIndexToValue)
{
    return term.coeff_
        * std::sin(term.coeff_sin_*var_sin)
        * std::cos(term.coeff_cos_*var_cos)
        * thetaIndexToValue(term.theta_index_);
}
