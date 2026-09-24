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

    // Sweeper-policy benchmark at the observed small-step floor.
    // Preserve all terms up to the configured spatial degree with a
    // GradedSweeper and compare against threshold sweeping.  This tests
    // whether we can retain long-horizon symbolic dependence without using
    // a globally tiny coefficient threshold.
    const ExactDouble loose_tolerance=1e-2_x;
    const ExactDouble plateau_step=0.0025_x;

    auto run_sweeper_policy_probe =
        [&](String const& policy, Sweeper<FloatDP> const& probe_sweeper) {
            PreconditionedGradedTaylorSeriesIntegrator gronwall(
                StepMaximumError(loose_tolerance),probe_sweeper,
                lipschitz_tolerance=0.5_x,
                minimum_spacial_order=5,
                minimum_temporal_order=5,
                maximum_spacial_order=5,
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

            std::cerr << "[IntegratorSweeperPolicyBenchmark]"
                      << " method=GRONWALL"
                      << " policy=" << policy
                      << " spatial_order=5"
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

    run_sweeper_policy_probe(
        "threshold_1e-12",
        Sweeper<FloatDP>(ThresholdSweeper<FloatDP>(DoublePrecision(),1e-12)));
    run_sweeper_policy_probe(
        "threshold_1e-14",
        Sweeper<FloatDP>(ThresholdSweeper<FloatDP>(DoublePrecision(),1e-14)));
    run_sweeper_policy_probe(
        "graded_degree_5",
        Sweeper<FloatDP>(GradedSweeper<FloatDP>(DoublePrecision(),5u)));
    run_sweeper_policy_probe(
        "graded_threshold_degree_5_1e-14",
        Sweeper<FloatDP>(
            GradedThresholdSweeper<FloatDP>(
                DoublePrecision(),5u,FloatDP(1e-14_x,DoublePrecision()))));

}
