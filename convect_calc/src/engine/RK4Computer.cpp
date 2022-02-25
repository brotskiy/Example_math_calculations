#include "RK4Computer.h"

#include <exception>
#include <algorithm>

void math::RK4Computer::RK4(const std::function<void(const double&, const std::vector<double>&, std::vector<double>&)> &system,
                            const double &h,
                            const double &timePrev,
                            const std::vector<double> &inputAtTimePrev,
                            std::vector<double> &outputAtTimeNext) const
{
    if (inputAtTimePrev.size() != outputAtTimeNext.size())
        throw std::runtime_error("<RK4Computer::RK4()>: data vectors dimension missmatch!");

    K1_.resize(inputAtTimePrev.size());    // размер НЕ уменьшается, если уже выделено достаточное количество элементов!
    K2_.resize(inputAtTimePrev.size());
    K3_.resize(inputAtTimePrev.size());
    K4_.resize(inputAtTimePrev.size());
    argument_.resize(inputAtTimePrev.size());

    system(timePrev, inputAtTimePrev, K1_);                  // вычислить значение K1 для каждой компоненты решения системы.

    createArgument(inputAtTimePrev, K1_, h/2, argument_);    // вычислить аргумент, передаваемый в правую часть системы для нахождения K2.
    system(timePrev + h/2, argument_, K2_);                  // вычислить значение K2 для каждой компоненты решения системы.

    createArgument(inputAtTimePrev, K2_, h/2, argument_);    // вычислить аргумент, передаваемый в правую часть системы для нахождения K3.
    system(timePrev + h/2, argument_, K3_);                  // вычислить значение K3 для каждой компоненты решения системы.

    createArgument(inputAtTimePrev, K3_, h, argument_);      // вычислить аргумент, передаваемый в правую часть системы для нахождения K4.
    system(timePrev + h, argument_, K4_);                    // вычислить значение K4 для каждой компоненты решения системы.

    for (size_t i = 0; i < outputAtTimeNext.size(); i++)
    {
        outputAtTimeNext.at(i) = RK4Formula(inputAtTimePrev.at(i), h, K1_.at(i), K2_.at(i), K3_.at(i), K4_.at(i));
    }
}

// :::::::::::::::::::::::::::::::::::::::: Полезности. ::::::::::::::::::::::::::::::::::::::::

void math::RK4Computer::createArgument(const std::vector<double> &values,
                                       const std::vector<double> &K,
                                       const double &coeff,
                                       std::vector<double> &output)
{
    if (values.size() != K.size() || values.size() != output.size())
        throw std::runtime_error("<RK4Computer::RK4()>: dimension missmatch!");

    std::transform(values.cbegin(), values.cend(), K.cbegin(), output.begin(),
    [&coeff](const double &val, const double& k) -> double
    {
        return val + coeff * k;
    });
}

double math::RK4Computer::RK4Formula(const double &val,
                                     const double &h,
                                     const double &k1,
                                     const double &k2,
                                     const double &k3,
                                     const double &k4)
{
    const double a = 1.0 / 6.0;
    const double b = 2.0 / 6.0;

    return val + h * (a*k1 + b*k2 + b*k3 + a*k4);
}
