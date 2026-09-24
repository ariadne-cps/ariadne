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

    auto run_product_accumulator_probe =
        [&](String const& policy, double threshold, Bool accumulator) {
            ThresholdSweeper<FloatDP> probe_sweeper(DoublePrecision(),threshold);
            PreconditionedGradedTaylorSeriesIntegrator gronwall(
                StepMaximumError(loose_tolerance),probe_sweeper,
                lipschitz_tolerance=0.5_x,
                minimum_spacial_order=5,minimum_temporal_order=5,
                maximum_spacial_order=5,maximum_temporal_order=5);
            gronwall.set_preconditioning(TaylorSeriesPreconditioning::QR);
            gronwall.set_diagnostics(false);
            gronwall.set_carried_expansion_diagnostics(false);

            VectorFieldEvolver evolver(dynamics,gronwall);
            configure_evolver(evolver,plateau_step);

            set_taylor_model_product_profile_enabled(false);
            set_taylor_model_early_discard_enabled(false);
            set_taylor_model_incremental_sweep_enabled(false);
            set_taylor_model_product_accumulator_enabled(accumulator);
            if(accumulator) {
                reset_taylor_model_accumulator_profile();
            }
            Stopwatch<Milliseconds> stopwatch;
            auto orbit=evolver.orbit(
                initial_set,Real(5.00_dec),Semantics::UPPER);
            stopwatch.click();
            set_taylor_model_product_accumulator_enabled(false);
            set_taylor_model_incremental_sweep_enabled(true);

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

            std::cerr << "[IntegratorProductAccumulatorBenchmark]"
                      << " policy=" << policy
                      << " accumulator=" << accumulator
                      << " elapsed_seconds=" << stopwatch.elapsed_seconds()
                      << " achieved_final_error=" << achieved_error
                      << " reach_sets=" << orbit.reach().size()
                      << std::endl;
            if(accumulator) {
                auto const ap=taylor_model_accumulator_profile();
                const double duplication_ratio=
                    ap.unique_entries
                        ? static_cast<double>(ap.temporary_entries)
                            / static_cast<double>(ap.unique_entries)
                        : 0.0;
                std::cerr << "[TaylorProductAccumulatorProfile]"
                          << " policy=" << policy
                          << " calls=" << ap.calls
                          << " product_pairs=" << ap.product_pairs
                          << " temporary_entries=" << ap.temporary_entries
                          << " unique_entries=" << ap.unique_entries
                          << " duplication_ratio=" << duplication_ratio
                          << " max_temporary_entries="
                          << ap.maximum_temporary_entries
                          << " max_unique_entries="
                          << ap.maximum_unique_entries
                          << " max_argument_size="
                          << ap.maximum_argument_size
                          << " max_x_degree="
                          << ap.maximum_x_degree
                          << " max_y_degree="
                          << ap.maximum_y_degree
                          << " max_product_degree="
                          << ap.maximum_product_degree
                          << " max_dense_slots="
                          << ap.maximum_dense_slots
                          << std::endl;
            }
        };

    // Profile accumulator compression at the tighter cutoff, where
    // product volume is largest and the next data-structure choice matters most.
    run_product_accumulator_probe("accumulator_profile_3e-14",3e-14,true);

}
