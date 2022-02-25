#include "NewtonComputer.h"

#include <cmath>
#include <stdexcept>

#include <iostream>
#include <iomanip>

double math::NewtonComputer::computeRoot(const double &t0,
                                         const std::function<double(const double&)> &F,
                                         const std::function<double(const double&)> &dF,
                                         const double &EPSILON,
                                         const int &MAX_STEP)
{
    int step_ = 1;

    double t     = t0;                  // значение свободной переменной равно начальному.
    double value = F(t0);               // значение функции равно начальному.

    while (std::abs(value) > EPSILON)   // продолжаем, пока значение функции от текущей свободной переменной не станет равным нулю.
    {
        t     = t - value / dF(t);      // вычисляем следующее приближение значения свободной переменной.
        value = F(t);                   // вычисляем значение функции от полученного значения свободной переменной.

        step_++;
        if (step_ > MAX_STEP)
            throw std::runtime_error("NewtonComputer::computeRoot(): step number exceeds MAX_STEP!");
    }

    return t;    // возвращаем корень - значение свободной переменной, при котором значение функции равно нулю.
}
