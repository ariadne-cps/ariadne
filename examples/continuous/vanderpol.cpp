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

            ARIADNE_ASSERT(!orbit.final().empty());
            auto achieved_error=
                orbit.final()[0u].state_function().get(0u).error();
            for(auto const& enclosure : orbit.final()) {
                for(SizeType i=0u;
                    i!=enclosure.state_function().result_size(); ++i) {
                    achieved_error=max(
                        achieved_error,
                        enclosure.state_function().get(i).error());
                }
            }

            std::cerr << "[IntegratorCarriedErrorBenchmark]"
                      << " method=" << method
                      << " tolerance=" << tolerance
                      << " max_step=" << max_step
                      << " horizon=5.0"
                      << " elapsed_seconds=" << stopwatch.elapsed_seconds()
                      << " achieved_final_error=" << achieved_error
                      << " reach_sets=" << orbit.reach().size()
                      << " intermediate_sets=" << orbit.intermediate().size()
                      << " final_sets=" << orbit.final().size()
                      << std::endl;
        };

    // Spatial-order plateau diagnostic.  Keep time horizon, max step,
    // temporal order, sweeper and loose local tolerance fixed; vary only the
    // spatial order of the flow/state representation.
    const ExactDouble loose_tolerance=1e-2_x;
    const ExactDouble plateau_step=0.0025_x;

    auto run_spatial_order_probe = [&](DegreeType spatial_order) {
        PreconditionedGradedTaylorSeriesIntegrator gronwall(
            StepMaximumError(loose_tolerance),sweeper,
            lipschitz_tolerance=0.5_x,
            minimum_spacial_order=spatial_order,
            minimum_temporal_order=5,
            maximum_spacial_order=spatial_order,
            maximum_temporal_order=5);
        gronwall.set_preconditioning(TaylorSeriesPreconditioning::QR);
        gronwall.set_diagnostics(false);

        VectorFieldEvolver evolver(dynamics,gronwall);
        configure_evolver(evolver,plateau_step);

        Stopwatch<Milliseconds> stopwatch;
        auto orbit=evolver.orbit(
            initial_set,Real(5.00_dec),Semantics::UPPER);
        stopwatch.click();

        ARIADNE_ASSERT(!orbit.final().empty());
        auto achieved_error=
            orbit.final()[0u].state_function().get(0u).error();
        for(auto const& enclosure : orbit.final()) {
            for(SizeType i=0u;
                i!=enclosure.state_function().result_size(); ++i) {
                achieved_error=max(
                    achieved_error,
                    enclosure.state_function().get(i).error());
            }
        }

        std::cerr << "[IntegratorSpatialOrderBenchmark]"
                  << " method=GRONWALL"
                  << " spatial_order=" << spatial_order
                  << " temporal_order=5"
                  << " tolerance=" << loose_tolerance
                  << " max_step=" << plateau_step
                  << " horizon=5.0"
                  << " elapsed_seconds=" << stopwatch.elapsed_seconds()
                  << " achieved_final_error=" << achieved_error
                  << " reach_sets=" << orbit.reach().size()
                  << " intermediate_sets=" << orbit.intermediate().size()
                  << " final_sets=" << orbit.final().size()
                  << std::endl;
    };

    run_spatial_order_probe(4u);
    run_spatial_order_probe(5u);
    run_spatial_order_probe(6u);

}
