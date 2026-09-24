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

#include <array>
#include <memory>

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

    // Snapshot the actually carried Taylor expansion rather than cumulative
    // sweeper events. The evolver diagnostic reports the input and output
    // normalised mappings at steps 500, 1000, 1500 and 2000.
    const ExactDouble loose_tolerance=1e-2_x;
    const ExactDouble plateau_step=0.0025_x;

    auto report_orbit =
        [&](String const& method, Stopwatch<Milliseconds> const& stopwatch,
            auto const& orbit) {
            ARIADNE_ASSERT(!orbit.final().empty());
            auto achieved_error=
                orbit.final()[0u].state_function().get(0u).error();
            auto final_radius=orbit.final()[0u].radius();
            for(auto const& enclosure : orbit.final()) {
                final_radius=max(final_radius,enclosure.radius());
                for(SizeType i=0u;
                    i!=enclosure.state_function().result_size(); ++i) {
                    achieved_error=max(
                        achieved_error,
                        enclosure.state_function().get(i).error());
                }
            }
            std::cerr << "[IntegratorArchitectureBenchmark]"
                      << " method=" << method
                      << " elapsed_seconds=" << stopwatch.elapsed_seconds()
                      << " achieved_final_error=" << achieved_error
                      << " final_radius=" << final_radius
                      << " reach_sets=" << orbit.reach().size()
                      << std::endl;
        };

    auto configure_taylor_kernel =
        [&](Bool dense) {
            set_taylor_model_product_profile_enabled(false);
            set_taylor_model_early_discard_enabled(false);
            set_taylor_model_incremental_sweep_enabled(!dense);
            set_taylor_model_product_accumulator_enabled(dense);
            set_taylor_model_dense_accumulator_enabled(dense);
        };

    auto reset_taylor_kernel = [&]() {
        set_taylor_model_dense_accumulator_enabled(false);
        set_taylor_model_product_accumulator_enabled(false);
        set_taylor_model_incremental_sweep_enabled(true);
    };

    auto run_graded =
        [&](String const& method, double threshold, Bool dense) {
            ThresholdSweeper<FloatDP> probe_sweeper(DoublePrecision(),threshold);
            GradedTaylorSeriesIntegrator integrator(
                StepMaximumError(loose_tolerance),probe_sweeper,
                lipschitz_tolerance=0.5_x,
                minimum_spacial_order=5,minimum_temporal_order=5,
                maximum_spacial_order=5,maximum_temporal_order=5);
            VectorFieldEvolver evolver(dynamics,integrator);
            configure_evolver(evolver,plateau_step);

            configure_taylor_kernel(dense);
            Stopwatch<Milliseconds> stopwatch;
            auto orbit=evolver.orbit(
                initial_set,Real(5.00_dec),Semantics::UPPER);
            stopwatch.click();
            reset_taylor_kernel();
            report_orbit(method,stopwatch,orbit);
        };

    auto run_preconditioned_dense = [&]() {
            ThresholdSweeper<FloatDP> probe_sweeper(DoublePrecision(),3e-14);
            PreconditionedGradedTaylorSeriesIntegrator integrator(
                StepMaximumError(loose_tolerance),probe_sweeper,
                lipschitz_tolerance=0.5_x,
                minimum_spacial_order=5,minimum_temporal_order=5,
                maximum_spacial_order=5,maximum_temporal_order=5);
            integrator.set_preconditioning(TaylorSeriesPreconditioning::QR);
            integrator.set_diagnostics(false);
            integrator.set_carried_expansion_diagnostics(false);

            VectorFieldEvolver evolver(dynamics,integrator);
            configure_evolver(evolver,plateau_step);

            configure_taylor_kernel(true);
            Stopwatch<Milliseconds> stopwatch;
            auto orbit=evolver.orbit(
                initial_set,Real(5.00_dec),Semantics::UPPER);
            stopwatch.click();
            reset_taylor_kernel();
            report_orbit("preconditioned_dense_3e-14",stopwatch,orbit);
        };

    // Map the graded+dense accuracy/runtime frontier below 3e-14 and
    // compare against the established preconditioned+dense 3e-14 target
    // error (~8.54e-8). The 3e-14 graded point is retained as an anchor.
    run_graded("graded_dense_3e-14",3e-14,true);
    run_graded("graded_dense_1e-14",1e-14,true);
    run_graded("graded_dense_3e-15",3e-15,true);
    run_graded("graded_dense_1e-15",1e-15,true);

}
