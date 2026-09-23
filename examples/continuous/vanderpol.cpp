/***************************************************************************
 *            vanderpol.cpp
 *
 *  Copyright  2017-20  Luca Geretti
 *
 ****************************************************************************/

#include "utility/stopwatch.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "dynamics/enclosure.hpp"
#include "ariadne_main.hpp"

void ariadne_main()
{
    CONCLOG_PRINTLN("van der Pol oscillator");

    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");
    VectorField dynamics({dot(x)=y, dot(y)=mu*y*(1-sqr(x))-x});

    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-12);

    Real x0=1.40_dec;
    Real y0=2.30_dec;
    Real eps_x0=0.15_dec;
    Real eps_y0=0.05_dec;
    RealExpressionBoundedConstraintSet initial_set({
        x0-eps_x0<=x<=x0+eps_x0,
        y0-eps_y0<=y<=y0+eps_y0
    });

    auto configure_evolver = [](VectorFieldEvolver& evolver, ExactDouble max_step) {
        evolver.configuration().set_maximum_enclosure_radius(1.0);
        evolver.configuration().set_maximum_step_size(max_step);
        evolver.configuration().set_maximum_spacial_error(1e-6);
        evolver.configuration().set_enable_reconditioning(false);
    };

    PreconditionedGradedTaylorSeriesIntegrator gronwall_integrator(
        StepMaximumError(1e-6),sweeper,lipschitz_tolerance=0.5_x,
        minimum_spacial_order=5,minimum_temporal_order=5,
        maximum_spacial_order=5,maximum_temporal_order=5);
    gronwall_integrator.set_preconditioning(TaylorSeriesPreconditioning::QR);
    gronwall_integrator.set_diagnostics(false);

    VectorFieldEvolver gronwall_evolver(dynamics,gronwall_integrator);
    configure_evolver(gronwall_evolver,0.02_x);

    Stopwatch<Milliseconds> gronwall_stopwatch;
    auto gronwall_orbit=
        gronwall_evolver.orbit(initial_set,Real(1.00_dec),Semantics::UPPER);
    gronwall_stopwatch.click();

    std::cerr << "[IntegratorBenchmark]"
              << " method=GRONWALL max_step=0.02"
              << " elapsed_seconds=" << gronwall_stopwatch.elapsed_seconds()
              << " reach_sets=" << gronwall_orbit.reach().size()
              << " intermediate_sets=" << gronwall_orbit.intermediate().size()
              << std::endl;

    GradedTaylorSeriesIntegrator graded_integrator(
        StepMaximumError(1e-6),sweeper,lipschitz_tolerance=0.5_x,
        minimum_spacial_order=5,minimum_temporal_order=5,
        maximum_spacial_order=5,maximum_temporal_order=5);

    VectorFieldEvolver graded_evolver(dynamics,graded_integrator);
    configure_evolver(graded_evolver,0.02_x);

    Stopwatch<Milliseconds> graded_stopwatch;
    auto graded_orbit=
        graded_evolver.orbit(initial_set,Real(1.00_dec),Semantics::UPPER);
    graded_stopwatch.click();

    std::cerr << "[IntegratorBenchmark]"
              << " method=GRADED max_step=0.02"
              << " elapsed_seconds=" << graded_stopwatch.elapsed_seconds()
              << " reach_sets=" << graded_orbit.reach().size()
              << " intermediate_sets=" << graded_orbit.intermediate().size()
              << std::endl;

    VectorFieldEvolver gronwall_evolver_004(dynamics,gronwall_integrator);
    configure_evolver(gronwall_evolver_004,0.04_x);

    Stopwatch<Milliseconds> gronwall_stopwatch_004;
    auto gronwall_orbit_004=
        gronwall_evolver_004.orbit(initial_set,Real(1.00_dec),Semantics::UPPER);
    gronwall_stopwatch_004.click();

    std::cerr << "[IntegratorBenchmark]"
              << " method=GRONWALL max_step=0.04"
              << " elapsed_seconds=" << gronwall_stopwatch_004.elapsed_seconds()
              << " reach_sets=" << gronwall_orbit_004.reach().size()
              << " intermediate_sets=" << gronwall_orbit_004.intermediate().size()
              << std::endl;

}
