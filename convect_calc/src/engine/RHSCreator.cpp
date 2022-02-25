#include "RHSCreator.h"

#include <cmath>
#include <stdexcept>

namespace
{
    const double ZERO_EPSILON = 1e-9;
}

// :::::::::::::::::::::::::::::::::::::::: Служебное. ::::::::::::::::::::::::::::::::::::::::

double math::NaturalConvection2D::RHSCreator::intSin2D(const int ik, const double bound)
{
    return bound/2;
}

double math::NaturalConvection2D::RHSCreator::intCosSin(const int i, const int ik, const double bound)
{
    if (i == ik)
    {
        return 0;
    }
    else
    {
        const double pi  = i * M_PI;
        const double pik = ik * M_PI;

        return bound * ik * (std::cos(pi) * std::cos(pik) - 1) / (M_PI * (i*i - ik*ik));
    }
}

double math::NaturalConvection2D::RHSCreator::intSinSin(const int i, const int ik, const double bound)
{
    if (i == ik)
        return bound/2;
    else
        return 0;
}

double math::NaturalConvection2D::RHSCreator::intSub(const int i, const int ii, const double bound)
{
    return -1 * bound/4;
}

double math::NaturalConvection2D::RHSCreator::intSinCosSin(const int i, const int ii, const int ik, const double bound)
{
    // существует 4 корня знаменателя:
    // 1) ik = i + ii;
    // 2) ik = ii - i;
    // 3) ik = i - ii;
    // 4) ik = -i - ii.
    // т.к. по способу задания индексы строго больше 0, то:
    // a) корень #4 не возможен, т.к. он был бы отрицательным;
    // b) случай ii==i в корнях #3 и #4 не возможен, т.к. они были бы нулевыми.

    if ((i <= 0) || (ii <= 0) || (ik <= 0))
        throw std::invalid_argument("index can't be zero or less.");

    if (ik == i + ii)
        return bound/4;

    if (ik == ii - i)
        return intSub(i, ii, bound);

    if (ik == i - ii)
        return -1 * intSub(i, ii, bound);

    return 0;
}

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

auto math::NaturalConvection2D::RHSCreator::createRHS(const double a, const double b, const int max_i, const int max_j)
    -> std::unique_ptr<SystemRHS>
{
    std::unique_ptr<SystemRHS> rhs = std::make_unique<SystemRHS>();

    const auto basis = createBasis(max_i, max_j);
    const auto eq1   = createEquation1(a, b, *basis);
    auto diffs       = createEquation2(a, b, *basis, *eq1);
    auto coords      = createCoords(*basis, *eq1, a, b);

    rhs->diffs_ = std::move(*diffs);
    rhs->eqX_   = std::move(coords->first);
    rhs->eqY_   = std::move(coords->second);

    return rhs;
}

auto math::NaturalConvection2D::RHSCreator::createBasis(const int max_i, const int max_j)
    -> std::unique_ptr<std::vector<BasisFunction>>
{
    auto basis = std::make_unique<std::vector<BasisFunction>>(max_i * max_j);
    auto funcIt = basis->begin();

    for (int i = 1; i <= max_i; i++)
    {
        for (int j = 1; j <= max_j; j++)
        {
            funcIt->i_ = i;
            funcIt->j_ = j;

            funcIt++;
        }
    }

    return basis;
}

auto math::NaturalConvection2D::RHSCreator::createEquation1(const double a, const double b, const std::vector<BasisFunction> &basis)
    -> std::unique_ptr<std::vector<std::vector<Equation1Term>>>
{
    auto equation1System = std::make_unique<std::vector<std::vector<Equation1Term>>>(basis.size());    // проекций столько, сколько базисных функций.
    auto eqIt = equation1System->begin();

    for (const BasisFunction &projFunc : basis)    // считаем проекцию первого уравнения на каждую базисную функцию.
    {
        const int ik = projFunc.i_;
        const int jk = projFunc.j_;

        const double xpr1 = (-1*ik*ik/(a*a) - jk*jk/(b*b)) * M_PI*M_PI * 4/(a*b) * intSin2D(ik, a) * intSin2D(jk, b);

        for (size_t thetaIndex = 0; thetaIndex < basis.size(); thetaIndex++)    // каждое PSI выражается через все THETA.
        {
            const int i = basis[thetaIndex].i_;
            const int j = basis[thetaIndex].j_;

            const double xpr2  = -1*i/a * M_PI * 4/(a*b) * intCosSin(i, ik, a) * intSinSin(j, jk, b);
            const double coeff = -1 * xpr2/xpr1;

            if (std::abs(coeff) > ZERO_EPSILON)
                eqIt->push_back(Equation1Term(coeff, thetaIndex));
        }

        eqIt++;    // переходим к заполнению следующего уравнения.
    }

    return equation1System;
}

auto math::NaturalConvection2D::RHSCreator::createEquation2(const double a, const double b,
                                                            const std::vector<BasisFunction> &basis,
                                                            const std::vector<std::vector<Equation1Term>> &equation1)
    -> std::unique_ptr<std::vector<std::vector<EquationTerm>>>
{
    const int basisSize = static_cast<int>(basis.size());
    auto diffs = std::make_unique<std::vector<std::vector<EquationTerm>>>(basisSize);

    for (int thetaIndex = 0; thetaIndex < basisSize; thetaIndex++)    // будем вычислять проекцию на каждую базисную функцию.
    {
        std::vector<EquationTerm> &currentEquation = diffs->at(thetaIndex);

        const int ik = basis[thetaIndex].i_;
        const int jk = basis[thetaIndex].j_;

        const double xpr1 = intSin2D(ik, a) * intSin2D(jk, b);

        const std::optional<EquationTerm> sub4 = createEquation2Sub4(xpr1, a, b, ik, jk, thetaIndex);
        if (sub4.has_value())
            currentEquation.push_back(sub4.value());

        for (int k = 0; k < basisSize; k++)
        {
            const int i = basis[k].i_;
            const int j = basis[k].j_;

            const std::optional<std::vector<EquationTerm>> sub5 = createEquation2Sub5(equation1.at(k), xpr1, a, b, ik, jk, i, j);
            if (sub5.has_value())
                currentEquation.insert(currentEquation.end(), sub5.value().cbegin(), sub5.value().cend());

            for (int n = 0; n < basisSize; n++)
            {
                const int ii = basis[n].i_;
                const int jj = basis[n].j_;

                const std::optional<std::vector<EquationTerm>> sub32 = createEquation2Sub32(equation1.at(k), xpr1, a, b, ik, jk, i, j, ii, jj, n);
                if (sub32.has_value())
                    currentEquation.insert(currentEquation.end(), sub32.value().cbegin(), sub32.value().cend());
            }
        }
    }

    return diffs;
}

auto math::NaturalConvection2D::RHSCreator::createEquation2Sub4(const double xpr1, const double a, const double b,
                                                                const int ik, const int jk, const int thetaIndex)
    -> std::optional<EquationTerm>
{
    const double xpr4  = M_PI*M_PI * (ik*ik/(a*a) + jk*jk/(b*b)) * intSin2D(ik, a) * intSin2D(jk, b);
    const double coeff = -1 * xpr4/xpr1;

    if (std::abs(coeff) > ZERO_EPSILON)
    {
        EquationTerm equation2Term;

        equation2Term.coeff_ = coeff;
        equation2Term.theta_indices_.push_back(thetaIndex);
        equation2Term.has_Rayleigh_ = false;

        return std::optional<EquationTerm>(equation2Term);
    }
    else
    {
        return std::nullopt;
    }
}

auto math::NaturalConvection2D::RHSCreator::createEquation2Sub5(const std::vector<Equation1Term> &equationPSI,
                                                                const double xpr1, const double a, const double b,
                                                                const int ik, const int jk, const int i, const int j)
    -> std::optional<std::vector<EquationTerm>>
{
    const double xpr5 = i * M_PI / a * intCosSin(i, ik, a) * intSinSin(j, jk, b);
    const double coeff = xpr5 / xpr1;

    if (std::abs(coeff) > ZERO_EPSILON)
    {
        std::vector<EquationTerm> result;

        for (const Equation1Term &termPSI : equationPSI)
        {
            EquationTerm termResult;

            termResult.has_Rayleigh_ = true;
            termResult.coeff_ = coeff * termPSI.coeff_;
            termResult.theta_indices_.push_back(termPSI.theta_index_);

            result.push_back(termResult);
        }

        return std::make_optional<std::vector<EquationTerm>>(std::move(result));
    }
    else
    {
        return std::nullopt;
    }
}

auto math::NaturalConvection2D::RHSCreator::createEquation2Sub32(const std::vector<Equation1Term> &equationPSI,
                                                                 const double xpr1, const double a, const double b,
                                                                 const int ik, const int jk, const int i, const int j, const int ii, const int jj, const int n)
    -> std::optional<std::vector<EquationTerm>>
{
    const double xpr2 = 2/(std::sqrt(a*b)) * M_PI*M_PI * j/b * ii/a * intSinCosSin(i, ii, ik, a) * intSinCosSin(jj, j, jk, b);
    const double xpr3 = 2/(std::sqrt(a*b)) * M_PI*M_PI * i/a * jj/b * intSinCosSin(ii, i, ik, a) * intSinCosSin(j, jj, jk, b);
    const double coeff = (xpr3-xpr2) / xpr1;

    if (std::abs(coeff) > ZERO_EPSILON)
    {
        std::vector<EquationTerm> result;

        for (const Equation1Term &termPSI : equationPSI)
        {
            EquationTerm termResult;

            termResult.has_Rayleigh_ = false;
            termResult.coeff_ = coeff * termPSI.coeff_;
            termResult.theta_indices_.push_back(n);
            termResult.theta_indices_.push_back(termPSI.theta_index_);

            result.push_back(termResult);
        }

        return std::make_optional<std::vector<EquationTerm>>(std::move(result));
    }
    else
    {
        return std::nullopt;
    }
}

auto math::NaturalConvection2D::RHSCreator::createCoords(const std::vector<BasisFunction> &basis,
                                                         const std::vector<std::vector<Equation1Term>> &equation1,
                                                         const double a, const double b)
    -> std::unique_ptr<std::pair<std::vector<CoordinateTerm>, std::vector<CoordinateTerm>>>
{
    auto coordEquation = std::make_unique<std::pair<std::vector<CoordinateTerm>, std::vector<CoordinateTerm>>>();

    for (size_t basisIndex = 0; basisIndex < basis.size(); basisIndex++)
    {
        const BasisFunction &currentBasis            = basis.at(basisIndex);
        const std::vector<Equation1Term> &currentPSI = equation1.at(basisIndex);

        const int i = currentBasis.i_;
        const int j = currentBasis.j_;

        for (const Equation1Term &termPSI : currentPSI)
        {
            const double coeffX = -1 * 2/(std::sqrt(a*b)) * j*M_PI/b * termPSI.coeff_;
            if (std::abs(coeffX) > ZERO_EPSILON)
            {
                CoordinateTerm termX;
                termX.coeff_     = coeffX;
                termX.coeff_sin_ = i * M_PI / a;
                termX.coeff_cos_ = j * M_PI / b;
                termX.theta_index_ = termPSI.theta_index_;

                coordEquation->first.push_back(termX);
            }

            const double coeffY = 2/(std::sqrt(a*b)) * i*M_PI/a * termPSI.coeff_;
            if (std::abs(coeffY) > ZERO_EPSILON)
            {
                CoordinateTerm termY;
                termY.coeff_     = coeffY;
                termY.coeff_cos_ = i * M_PI / a;
                termY.coeff_sin_ = j * M_PI / b;
                termY.theta_index_ = termPSI.theta_index_;

                coordEquation->second.push_back(termY);
            }
        }
    }

    return coordEquation;
}
