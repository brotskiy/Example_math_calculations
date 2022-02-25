#ifndef POINCARESECTIONSTRUCTS_H
#define POINCARESECTIONSTRUCTS_H

#include <vector>
#include <utility>

namespace math
{
namespace PoincareSection
{

struct PlaneEquation
{
    std::vector<double> coeffs_;    // коэффициенты уравнения плоскости перед переменными.
    double coeffFree_;              // свободный коэффициент уравнения плоскости.

    PlaneEquation()                                = default;
    PlaneEquation(const PlaneEquation&)            = default;
    PlaneEquation& operator=(const PlaneEquation&) = default;
    ~PlaneEquation()                               = default;

    PlaneEquation(PlaneEquation&&);
    PlaneEquation& operator=(PlaneEquation&&);
};

}
}



inline math::PoincareSection::PlaneEquation::PlaneEquation(PlaneEquation&& other) :
    coeffs_(std::move(other.coeffs_)),
    coeffFree_(other.coeffFree_)
{
}

inline auto math::PoincareSection::PlaneEquation::operator=(PlaneEquation&& other) -> PlaneEquation&
{
    coeffs_    = std::move(other.coeffs_);
    coeffFree_ = other.coeffFree_;
    return *this;
}

#endif // POINCARESECTIONSTRUCTS_H
