#ifndef NEWTONCOMPUTER_H
#define NEWTONCOMPUTER_H

#include <functional>

namespace math
{

class NewtonComputer
{
public:
    static double computeRoot(const double &t0,
                              const std::function<double(const double&)> &F,
                              const std::function<double(const double&)> &dF,
                              const double &EPSILON,
                              const int &MAX_STEP);
};

}

#endif // NEWTONCOMPUTER_H
