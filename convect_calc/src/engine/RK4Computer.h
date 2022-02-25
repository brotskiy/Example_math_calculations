#ifndef RK4_H
#define RK4_H

#include <vector>
#include <functional>

namespace math
{

class RK4Computer
{
// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

private:
    mutable std::vector<double> K1_, K2_, K3_, K4_, argument_;

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public:
    void RK4(const std::function<void(const double&, const std::vector<double>&, std::vector<double>&)> &system,
             const double &h,
             const double &timePrev,
             const std::vector<double> &inputAtTimePrev,
             std::vector<double> &outputAtTimeNext) const;

// :::::::::::::::::::::::::::::::::::::::: Полезности. ::::::::::::::::::::::::::::::::::::::::

private:
    static void createArgument(const std::vector<double> &values,
                               const std::vector<double> &K,
                               const double &coeff,
                               std::vector<double> &output);

    static double RK4Formula(const double &val,
                             const double &h,
                             const double &k1,
                             const double &k2,
                             const double &k3,
                             const double &k4);
};

}

#endif // RK4_H
