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

    PreconditionedGradedTaylorSeriesIntegrator residual_probe(
        StepMaximumError(1e-6),sweeper,lipschitz_tolerance=0.5_x,
        minimum_spacial_order=5,minimum_temporal_order=5,
        maximum_spacial_order=5,maximum_temporal_order=5);
    residual_probe.set_preconditioning(TaylorSeriesPreconditioning::QR);
    residual_probe.set_diagnostics(true);
    VectorFieldEvolver residual_probe_evolver(dynamics,residual_probe);
    residual_probe_evolver.configuration().set_maximum_enclosure_radius(1.0);
    residual_probe_evolver.configuration().set_maximum_step_size(0.02);
    residual_probe_evolver.configuration().set_maximum_spacial_error(1e-6);
    residual_probe_evolver.configuration().set_enable_reconditioning(false);
    auto residual_probe_orbit=
        residual_probe_evolver.orbit(initial_set,Real(0.02_dec),Semantics::UPPER);
    std::cerr << "[RecurrenceResidualProbe]"
              << " reach_sets=" << residual_probe_orbit.reach().size()
              << std::endl;

    auto configure_evolver = [](VectorFieldEvolver& evolver, ExactDouble max_step) {
        evolver.configuration().set_maximum_enclosure_radius(1.0);
        evolver.configuration().set_maximum_step_size(max_step);
        evolver.configuration().set_maximum_spacial_error(1e-6);
        evolver.configuration().set_enable_reconditioning(false);
    };

    auto run_long_benchmark =
        [&](String const& method,
            IntegratorInterface const& integrator,
            ExactDouble tolerance,
            ExactDouble max_step) {
            VectorFieldEvolver evolver(dynamics,integrator);
            configure_evolver(evolver,max_step);

            Stopwatch<Milliseconds> stopwatch;
            auto orbit=evolver.orbit(
                initial_set,Real(5.00_dec),Semantics::UPPER);
            stopwatch.click();

            std::cerr << "[IntegratorLongBenchmark]"
                      << " method=" << method
                      << " tolerance=" << tolerance
                      << " max_step=" << max_step
                      << " horizon=5.0"
                      << " elapsed_seconds=" << stopwatch.elapsed_seconds()
                      << " reach_sets=" << orbit.reach().size()
                      << " intermediate_sets=" << orbit.intermediate().size()
                      << std::endl;
        };

    PreconditionedGradedTaylorSeriesIntegrator gronwall_1e6(
        StepMaximumError(1e-6),sweeper,lipschitz_tolerance=0.5_x,
        minimum_spacial_order=5,minimum_temporal_order=5,
        maximum_spacial_order=5,maximum_temporal_order=5);
    gronwall_1e6.set_preconditioning(TaylorSeriesPreconditioning::QR);
    gronwall_1e6.set_diagnostics(false);

    GradedTaylorSeriesIntegrator graded_1e6(
        StepMaximumError(1e-6),sweeper,lipschitz_tolerance=0.5_x,
        minimum_spacial_order=5,minimum_temporal_order=5,
        maximum_spacial_order=5,maximum_temporal_order=5);

    run_long_benchmark("GRONWALL",gronwall_1e6,1e-6_x,0.04_x);
    run_long_benchmark("GRADED",graded_1e6,1e-6_x,0.04_x);

    PreconditionedGradedTaylorSeriesIntegrator gronwall_1e8(
        StepMaximumError(1e-8),sweeper,lipschitz_tolerance=0.5_x,
        minimum_spacial_order=5,minimum_temporal_order=5,
        maximum_spacial_order=5,maximum_temporal_order=5);
    gronwall_1e8.set_preconditioning(TaylorSeriesPreconditioning::QR);
    gronwall_1e8.set_diagnostics(false);

    GradedTaylorSeriesIntegrator graded_1e8(
        StepMaximumError(1e-8),sweeper,lipschitz_tolerance=0.5_x,
        minimum_spacial_order=5,minimum_temporal_order=5,
        maximum_spacial_order=5,maximum_temporal_order=5);

    run_long_benchmark("GRONWALL",gronwall_1e8,1e-8_x,0.04_x);
    run_long_benchmark("GRADED",graded_1e8,1e-8_x,0.04_x);

}
