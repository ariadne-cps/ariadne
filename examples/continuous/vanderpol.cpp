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
    Real evolution_time=7;

    auto print_result =
        [&](const char* family, const char* setting,
            Stopwatch<Milliseconds> const& sw, auto const& orbit) {
            std::cerr << "[EqualAccuracyComparison]"
                      << " family=" << family
                      << " setting=" << setting
                      << " time_s=" << sw.elapsed_seconds()
                      << " reach_sets=" << orbit.reach().size()
                      << " intermediate_sets=" << orbit.intermediate().size()
                      << " final_sets=" << orbit.final().size();
            if(!orbit.final().empty()) {
                auto const& final_set=orbit.final()[0];
                std::cerr << " final_params=" << final_set.number_of_parameters()
                          << " final_radius=" << final_set.radius()
                          << " final_error=" << final_set.state_function().error()
                          << " final_box=" << final_set.euclidean_set().bounding_box();
            }
            std::cerr << std::endl;
        };

    auto run_graded =
        [&](StepMaximumError max_err, const char* setting) {
            GradedTaylorSeriesIntegrator integrator(
                max_err,sweeper,lipschitz_tolerance=0.5_x,
                minimum_spacial_order=5,minimum_temporal_order=5,
                maximum_spacial_order=5,maximum_temporal_order=5);

            VectorFieldEvolver evolver(dynamics,integrator);
            evolver.configuration().set_maximum_enclosure_radius(1.0);
            evolver.configuration().set_maximum_step_size(0.02);
            evolver.configuration().set_maximum_spacial_error(1e-6);
            evolver.configuration().set_enable_reconditioning(false);

            Stopwatch<Milliseconds> sw;
            sw.restart();
            auto orbit=evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
            sw.click();
            print_result("GradedTaylorSeries",setting,sw,orbit);
        };

    auto run_preconditioned_qr =
        [&](StepMaximumError max_err, const char* setting) {
            PreconditionedGradedTaylorSeriesIntegrator integrator(
                max_err,sweeper,lipschitz_tolerance=0.5_x,
                minimum_spacial_order=5,minimum_temporal_order=5,
                maximum_spacial_order=5,maximum_temporal_order=5);
            integrator.set_preconditioning(TaylorSeriesPreconditioning::QR);

            VectorFieldEvolver evolver(dynamics,integrator);
            evolver.configuration().set_maximum_enclosure_radius(1.0);
            evolver.configuration().set_maximum_step_size(0.02);
            evolver.configuration().set_maximum_spacial_error(1e-6);
            evolver.configuration().set_enable_reconditioning(false);

            Stopwatch<Milliseconds> sw;
            sw.restart();
            auto orbit=evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
            sw.click();
            print_result("PreconditionedGradedTaylorSeries",setting,sw,orbit);
        };

    run_graded(StepMaximumError(1e-6),"1e-6");
    run_graded(StepMaximumError(3e-7),"3e-7");
    run_graded(StepMaximumError(1e-7),"1e-7");
    run_graded(StepMaximumError(3e-8),"3e-8");
    run_preconditioned_qr(StepMaximumError(1e-6),"QR-1e-6");
}
