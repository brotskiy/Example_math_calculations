#ifndef RHSCOMPUTER_H
#define RHSCOMPUTER_H

#include <vector>
#include <functional>

#include "engine/NaturalConvection2DStructs.h"

namespace math
{

class RHSComputer
{

// :::::::::::::::::::::::::::::::::::::::: Дифференциальное уравнение координаты частицы. ::::::::::::::::::::::::::::::::::::::::

public:
    static double computeCoordinateEquation(const std::vector<math::NaturalConvection2D::CoordinateTerm> &coordEq,
                                            const double &var_sin,
                                            const double &var_cos,
                                            const std::function<double(const int&)> &thetaIndexToValue);

private:
    static double computeCoordinateTerm(const math::NaturalConvection2D::CoordinateTerm &term,
                                        const double &var_sin,
                                        const double &var_cos,
                                        const std::function<double(const int&)> &thetaIndexToValue);
};

}

#endif // RHSCOMPUTER_H
