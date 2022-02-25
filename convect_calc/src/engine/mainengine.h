#ifndef MAINENGINE_H
#define MAINENGINE_H

#include <QObject>
#include <QString>

#include <vector>
#include <memory>
#include <functional>

#include "NaturalConvection2DStructs.h"
#include "PoincareSectionStructs.h"

class RK4Engine : public QObject
{
    Q_OBJECT

// :::::::::::::::::::::::::::::::::::::::: Содержимое. ::::::::::::::::::::::::::::::::::::::::

public:
    struct Settings
    {
        double a_, b_;
        double t_beg_, t_end_;
        double rayleigh_;
        int eq_count_;
        int step_count_;

        int repeat_count_;

        QString file_path_thetas_;
        QString file_path_particle_;

        std::vector<math::NaturalConvection2D::Particle> particles_;

        Settings();
    };
    static void printSettings(const RK4Engine::Settings &sttngs);

private:
    struct KCoordCoeff
    {
        double x_, y_;
    };

private:
    const Settings settings_;

// :::::::::::::::::::::::::::::::::::::::: Создание. ::::::::::::::::::::::::::::::::::::::::

public:
    RK4Engine(const RK4Engine::Settings &engnSttngs, QObject* parent);
    ~RK4Engine() = default;

// :::::::::::::::::::::::::::::::::::::::: Работа. ::::::::::::::::::::::::::::::::::::::::

public slots:
    void onStart();
    void onPoincare();
    void onComputeStreamFunctionCoefficients();
    void onComputeStreamFunction();

// :::::::::::::::::::::::::::::::::::::::: Функция тока. ::::::::::::::::::::::::::::::::::::::::

private:
    static bool SF_step_computeCoefficients(const std::function<bool(double&, std::vector<double>&)> &getThetas,
                                            const std::function<void(const double&, const std::vector<double>&)> &storePsis,
                                            const std::vector<std::vector<math::NaturalConvection2D::Equation1Term>> &eq1);

    static auto SF_computeEquation1(const std::vector<double> &thetas,
                                    const std::vector<std::vector<math::NaturalConvection2D::Equation1Term>> &eq1)
        -> std::vector<double>;

    static double SF_computePsiEquation(const std::vector<math::NaturalConvection2D::Equation1Term> &psiEq,
                                        const std::function<double(const int)> &thetaIndexToValue);

    static double SF_computePsiTerm(const math::NaturalConvection2D::Equation1Term &term,
                                    const std::function<double(const int)> &thetaIndexToValue);

    static bool SF_step_computeStreamFunction(const double &a,
                                              const double &b,
                                              const std::vector<math::NaturalConvection2D::Particle> &points,
                                              const std::vector<math::NaturalConvection2D::BasisFunction> &basis,
                                              const std::function<bool(double&, std::vector<double>&)> &getPsis,
                                              const std::function<void(const double&,
                                                                       const std::vector<math::NaturalConvection2D::Particle>&,
                                                                       const std::vector<double>&)> &storeStreamFunction);

    static double SF_computeStreamFunctionEquation(const double &a,
                                                   const double &b,
                                                   const std::vector<double> &psis,
                                                   const std::vector<math::NaturalConvection2D::BasisFunction> &basis,
                                                   const math::NaturalConvection2D::Particle &point);

    static double SF_computeStreamFunctionTerm(const double &a,
                                               const double &b,
                                               const double &psiCoeff,
                                               const math::NaturalConvection2D::BasisFunction &basisFunction,
                                               const math::NaturalConvection2D::Particle &point);



// :::::::::::::::::::::::::::::::::::::::: Сечение Пуанкаре. ::::::::::::::::::::::::::::::::::::::::

private:
    static bool PN_step(const std::function<bool(double&,
                                                 std::vector<double>&,
                                                 std::vector<math::NaturalConvection2D::Particle>&)> &getData,
                        const std::function<void(const double&,
                                                 const std::vector<double>&,
                                                 const std::vector<math::NaturalConvection2D::Particle>&)> &storeData,
                        const double &Rayleigh,
                        const math::PoincareSection::PlaneEquation &planeEquation,
                        const math::NaturalConvection2D::SystemRHS &rhs);

    static auto PN_tryComputeIntersectionTime(const double &Rayleigh,
                                              const double &t_prev,
                                              const double &t_next,
                                              const std::vector<double> &thetas_next,
                                              const math::PoincareSection::PlaneEquation &planeEquation,
                                              const std::vector<std::vector<math::NaturalConvection2D::EquationTerm>> &diffRHS)
        -> std::optional<double>;

    static double PN_newtonStep(const double &Rayleigh,
                                const double &time,
                                const std::vector<double> &thetas,
                                const math::PoincareSection::PlaneEquation &planeEquation,
                                const std::vector<std::vector<math::NaturalConvection2D::EquationTerm>> &diffRHS);

    static double PN_applyPlane(const std::vector<double> &thetas,
                                const math::PoincareSection::PlaneEquation &planeEquation);

    static double PN_applyPlaneDerivative(const double &Rayleigh,
                                          const std::vector<double> &thetas,
                                          const math::PoincareSection::PlaneEquation &planeEquation,
                                          const std::vector<std::vector<math::NaturalConvection2D::EquationTerm>> &diffRHS);

    static auto PN_createPlane(const size_t dimension)
        -> std::unique_ptr<math::PoincareSection::PlaneEquation>;

// :::::::::::::::::::::::::::::::::::::::: Рунге-Кутта. ::::::::::::::::::::::::::::::::::::::::

private:
    static void RK4_step(const math::NaturalConvection2D::SystemRHS &rhs,
                         const double &Rayleigh,
                         const double &h,
                         const int step,
                         std::vector<std::vector<double>> &KThetaBuffers,
                         math::NaturalConvection2D::SystemData &data_in_out);

public:
    static void RK4_computeDiff(const std::vector<std::vector<math::NaturalConvection2D::EquationTerm>> &diffRHS,
                                const double &Rayleigh,
                                const double &h,
                                const std::vector<double> &dataPrevStep,
                                std::vector<double> &dataCurStep_out,
                                std::vector<std::vector<double>> &KThetaVectors_out);

    static void RK4_computeDiffKVector(const std::vector<std::vector<math::NaturalConvection2D::EquationTerm>> &diffRHS,
                                       const double &Rayleigh,
                                       const std::vector<double> &dataPrevStep,
                                       const std::optional<double> &factorK,
                                       const std::vector<double> &KVector_in,
                                       std::vector<double> &KVector_out);

    static double RK4_computeDiffEquation(const std::vector<math::NaturalConvection2D::EquationTerm> &diffEquation,
                                          const double &Rayleigh,
                                          const std::function<double(const int)> &thetaIndexToArgument);

private:
    static double RK4_computeDiffTerm(const math::NaturalConvection2D::EquationTerm &term,
                                      const double &Rayleigh,
                                      const std::function<double(const int)> &thetaIndexToArgument);

    static void RK4_computeCoord(const std::vector<math::NaturalConvection2D::CoordinateTerm> &equationX,
                                 const std::vector<math::NaturalConvection2D::CoordinateTerm> &equationY,
                                 const double &h,
                                 const std::vector<double> &thetasPrevStep,
                                 const std::vector<std::vector<double>> &KThetaVectors,
                                 const std::vector<math::NaturalConvection2D::Particle> &particlesPrevStep,
                                 std::vector<math::NaturalConvection2D::Particle> &particlesCurStep_out);

    static void RK4_computeCoordK(const std::vector<math::NaturalConvection2D::CoordinateTerm> &equationX,
                                  const std::vector<math::NaturalConvection2D::CoordinateTerm> &equationY,
                                  const std::vector<double> &thetasPrevStep,
                                  const math::NaturalConvection2D::Particle &particlePrevStep,
                                  const std::optional<double> &factorK,
                                  const std::vector<double> &KThetaVector_in,
                                  const KCoordCoeff &K_in,
                                  KCoordCoeff &K_out);

    static double RK4_computeCoordEquation(const std::vector<math::NaturalConvection2D::CoordinateTerm> &coordinateEquation,
                                           const double &arg_sin, const double &arg_cos,
                                           const std::function<double(const int)> &thetaIndexToArgument);

    static double RK4_computeCoordTerm(const math::NaturalConvection2D::CoordinateTerm &term,
                                       const double &arg_sin, const double &arg_cos,
                                       const std::function<double(const int)> &thetaIndexToArgument);

    static double RK4_formula(const double &scalarPrevStep, const double &h, const double &K1, const double &K2, const double &K3, const double &K4);

// :::::::::::::::::::::::::::::::::::::::: Полезности. ::::::::::::::::::::::::::::::::::::::::

private:
    static auto createData(const RK4Engine::Settings &settings)
        -> std::unique_ptr<math::NaturalConvection2D::SystemData>;

    static auto create_K_coefficients(const int eq_count)
        -> std::unique_ptr<std::vector<std::vector<double>>>;

    void writeData(const math::NaturalConvection2D::SystemData &data, const size_t endIndex) const;

    static void cycleData(math::NaturalConvection2D::SystemData &data);

    void maybeEmitInfo(const int step, const int elapsed_ms, const std::vector<math::NaturalConvection2D::Particle> &currentCoordinates);

signals:
    void stepInfo(int step, int elapsed_ms);
    void particlesCoordinates(const std::vector<math::NaturalConvection2D::Particle> &coordinates);

// :::::::::::::::::::::::::::::::::::::::: Ограничения. ::::::::::::::::::::::::::::::::::::::::

private:
    RK4Engine()                            = delete;
    RK4Engine(const RK4Engine&)            = delete;
    RK4Engine(RK4Engine&&)                 = delete;
    RK4Engine& operator=(const RK4Engine&) = delete;
    RK4Engine& operator=(RK4Engine&&)      = delete;
};

#endif // MAINENGINE_H
