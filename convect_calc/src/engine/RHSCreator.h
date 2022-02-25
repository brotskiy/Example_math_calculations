#ifndef RHSCREATOR_H
#define RHSCREATOR_H

#include <vector>
#include <memory>
#include <optional>
#include <utility>

#include "NaturalConvection2DStructs.h"

namespace math
{
namespace NaturalConvection2D
{

class RHSCreator
{

// :::::::::::::::::::::::::::::::::::::::: Служебное. ::::::::::::::::::::::::::::::::::::::::

private:
    static double intSin2D(const int ik, const double bound);
    static double intCosSin(const int i, const int ik, const double bound);
    static double intSinSin(const int i, const int ik, const double bound);
    static double intSub(const int i, const int ii, const double bound);
    static double intSinCosSin(const int i, const int ii, const int ik, const double bound);

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public:
    static auto createRHS(const double a, const double b, const int max_i, const int max_j)
        -> std::unique_ptr<SystemRHS>;

    static auto createBasis(const int max_i, const int max_j)
        -> std::unique_ptr<std::vector<BasisFunction>>;

    static auto createEquation1(const double a, const double b, const std::vector<BasisFunction> &basis)
        -> std::unique_ptr<std::vector<std::vector<Equation1Term>>>;

    static auto createCoords(const std::vector<BasisFunction> &basis,
                             const std::vector<std::vector<Equation1Term>> &equation1,
                             const double a, const double b)
        -> std::unique_ptr<std::pair<std::vector<CoordinateTerm>, std::vector<CoordinateTerm>>>;

private:
    static auto createEquation2(const double a, const double b,
                                const std::vector<BasisFunction> &basis,
                                const std::vector<std::vector<Equation1Term>> &equation1)
        -> std::unique_ptr<std::vector<std::vector<EquationTerm>>>;

    static auto createEquation2Sub4(const double xpr1, const double a, const double b,
                                    const int ik, const int jk, const int thetaIndex)
        -> std::optional<EquationTerm>;

    static auto createEquation2Sub5(const std::vector<Equation1Term> &equationPSI,
                                    const double xpr1, const double a, const double b,
                                    const int ik, const int jk, const int i, const int j)
        -> std::optional<std::vector<EquationTerm>>;

    static auto createEquation2Sub32(const std::vector<Equation1Term> &equationPSI,
                                     const double xpr1, const double a, const double b,
                                     const int ik, const int jk, const int i, const int j, const int ii, const int jj, const int n)
        -> std::optional<std::vector<EquationTerm>>;

};

}
}

#endif // RHSCREATOR_H
