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

    auto run_carried_expansion_probe =
        [&](String const& policy, double threshold) {
            ThresholdSweeper<FloatDP> probe_sweeper(DoublePrecision(),threshold);
            PreconditionedGradedTaylorSeriesIntegrator gronwall(
                StepMaximumError(loose_tolerance),probe_sweeper,
                lipschitz_tolerance=0.5_x,
                minimum_spacial_order=5,minimum_temporal_order=5,
                maximum_spacial_order=5,maximum_temporal_order=5);
            gronwall.set_preconditioning(TaylorSeriesPreconditioning::QR);
            gronwall.set_diagnostics(false);
            gronwall.set_carried_expansion_diagnostics(true);

            VectorFieldEvolver evolver(dynamics,gronwall);
            configure_evolver(evolver,plateau_step);

            reset_taylor_model_product_profile();
            set_taylor_model_product_profile_enabled(true);
            Stopwatch<Milliseconds> stopwatch;
            auto orbit=evolver.orbit(
                initial_set,Real(5.00_dec),Semantics::UPPER);
            stopwatch.click();
            set_taylor_model_product_profile_enabled(false);
            auto const product_profile=taylor_model_product_profile_snapshot();

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

            std::cerr << "[IntegratorCarriedExpansionBenchmark]"
                      << " policy=" << policy
                      << " elapsed_seconds=" << stopwatch.elapsed_seconds()
                      << " achieved_final_error=" << achieved_error
                      << " reach_sets=" << orbit.reach().size()
                      << std::endl;

            auto print_product_profile =
                [&](String const& context,
                    TaylorModelProductProfileCounters const& profile) {
                    ARIADNE_ASSERT(
                        profile.individual_products_below_threshold
                        + profile.individual_products_above_threshold
                        == profile.product_pairs);
                    ARIADNE_ASSERT(
                        profile.individual_products_below_threshold_collision
                        + profile.individual_products_below_threshold_new_term
                        + profile.individual_products_below_threshold_trailing
                        == profile.individual_products_below_threshold);
                    std::cerr << "[TaylorProductGenerationProfile]"
                              << " policy=" << policy
                              << " context=" << context
                              << " calls=" << profile.calls
                              << " product_pairs=" << profile.product_pairs
                              << " individual_below_threshold="
                              << profile.individual_products_below_threshold
                              << " individual_below_collision="
                              << profile.individual_products_below_threshold_collision
                              << " individual_below_new_term="
                              << profile.individual_products_below_threshold_new_term
                              << " individual_below_trailing="
                              << profile.individual_products_below_threshold_trailing
                              << " individual_above_threshold="
                              << profile.individual_products_above_threshold
                              << " individual_below_threshold_abs_mass="
                              << profile.individual_products_below_threshold_abs_mass
                              << " individual_above_threshold_abs_mass="
                              << profile.individual_products_above_threshold_abs_mass
                              << " sweep_passes=" << profile.sweep_passes
                              << " sweep_input_terms=" << profile.sweep_input_terms
                              << " sweep_output_terms=" << profile.sweep_output_terms
                              << " swept_terms=" << profile.swept_terms
                              << " max_sweep_input_terms="
                              << profile.maximum_sweep_input_terms
                              << " max_sweep_output_terms="
                              << profile.maximum_sweep_output_terms
                              << " elapsed_seconds=" << profile.elapsed_seconds
                              << std::endl;
                };
            print_product_profile("general",product_profile.general);
            print_product_profile("compose",product_profile.compose);
        };

    run_carried_expansion_probe("absolute_1e-12",1e-12);
    run_carried_expansion_probe("absolute_3e-14",3e-14);

}
