#ifndef NATURALCONVECTION2DSTRUCTS_H
#define NATURALCONVECTION2DSTRUCTS_H

#include <vector>
#include <utility>

namespace math
{
namespace NaturalConvection2D
{

struct Particle
{
    double x_, y_;
};

struct BasisFunction    // базисная функция вида PHI(x,y) = 2/sqrt(a*b) * sin(i * PI*x/a) * sin(j * PI*y/b).
{
    int i_, j_;
};

struct Equation1Term    // слагаемое в векторном уравнении 1, выражающем функцию тока PSI через отклонение температуры THETA.
{
    double coeff_;
    int theta_index_;   // THETA нумеруются от 0!

    Equation1Term() = default;
    Equation1Term(double coeff, double theta_index) : coeff_(coeff), theta_index_(theta_index) {}
};

struct EquationTerm    // слагаемое в результирующем векторном уравнении, выражающем производную THETA по времени через значения THETA.
{
    double coeff_;
    std::vector<int> theta_indices_;    // индексы всех THETA, из произведения которых и состоит слагаемое. THETA нумеруются от 0!
    bool has_Rayleigh_;                 // входит ли в данное слагаемое число Рэлея.
};

struct CoordinateTerm    // слагаемое в уравнении, выражающем скорость частицы через отклонения температуры THETA.
{
    double coeff_;
    int theta_index_;                 // идекс THETA, входящей в данное слагаемое.
    double coeff_sin_, coeff_cos_;    // числовые коэффициенты в аргументах синуса и косинуса.
};

struct SystemRHS    // результирующее векторное уравнение, выражающее производную THETA по времени через значения THETA.
{
    std::vector<std::vector<EquationTerm>> diffs_;    // правые части всех уравнений для THETA.
    std::vector<CoordinateTerm> eqX_;                 // выражение для нахождения X-координаты частицы.
    std::vector<CoordinateTerm> eqY_;                 // выражение для нахождения Y-координаты частицы.

    SystemRHS() = default;
    SystemRHS(const SystemRHS& other) = default;
    SystemRHS(SystemRHS &&other);
};

struct SystemData
{
    std::vector<std::vector<double>> diffs_;       // значения всех THETA_i на каждом шаге. в 0-ом столбце записаны начальные условия.
    std::vector<std::vector<Particle>> coords_;    // координаты всех частиц на каждом шаге. в 0-ом столбце записаны начальные положения.
    std::vector<double> time_;                     // значение времени на каждом шаге. в 0-ом элементе записано время старта.

    SystemData(const size_t stepCount,
               const size_t equationCount,
               const size_t particleCount);

    SystemData(const SystemData &other) = default;
    SystemData(SystemData &&other);
};

}
}

inline math::NaturalConvection2D::SystemRHS::SystemRHS(SystemRHS &&other) :
    diffs_(std::move(other.diffs_)),
    eqX_(std::move(other.eqX_)),
    eqY_(std::move(other.eqY_))
{
}

inline math::NaturalConvection2D::SystemData::SystemData(const size_t stepCount,
                                                         const size_t equationCount,
                                                         const size_t particleCount)
{
    diffs_.resize(stepCount);
    coords_.resize(stepCount);
    time_.resize(stepCount);

    for (std::vector<double> &thetasAtStep : diffs_)
        thetasAtStep.resize(equationCount);

    for (std::vector<Particle> &coordsAtStep : coords_)
        coordsAtStep.resize(particleCount);
}

inline math::NaturalConvection2D::SystemData::SystemData(SystemData &&other) :
    diffs_(std::move(other.diffs_)),
    coords_(std::move(other.coords_)),
    time_(std::move(other.time_))
{
}

#endif // NATURALCONVECTION2DSTRUCTS_H
