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

    // Carried-expansion structure probe.
    // Compare the two useful absolute-threshold configurations and sample the
    // carried state expansion late in the run.  We need to understand which
    // degrees and coefficient magnitudes account for the accuracy gain of
    // 1e-14 over 1e-12 before designing a new selective sweeper.
    const ExactDouble loose_tolerance=1e-2_x;
    const ExactDouble plateau_step=0.0025_x;

    auto run_structure_probe =
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

            auto const& final_state=orbit.final()[0u].state_function();
            std::cerr << "[IntegratorExpansionStructureBenchmark]"
                      << " method=GRONWALL"
                      << " policy=" << policy
                      << " elapsed_seconds=" << stopwatch.elapsed_seconds()
                      << " achieved_final_error=" << achieved_error
                      << " final_errors=" << final_state.errors()
                      << std::endl;

            for(SizeType i=0u; i!=final_state.result_size(); ++i) {
                auto taylor_model=
                    std::dynamic_pointer_cast<
                        const ValidatedScalarMultivariateTaylorFunctionModelDP>(
                            final_state.get(i).managed_pointer());
                ARIADNE_ASSERT(taylor_model!=nullptr);
                auto const& expansion=taylor_model->model().expansion();
                SizeType max_degree=0u;
                for(auto const& term : expansion) {
                    max_degree=max(max_degree,term.index().degree());
                }
                for(SizeType degree=0u; degree<=max_degree; ++degree) {
                    SizeType count=0u;
                    ApproximateDouble max_abs=0.0;
                    ApproximateDouble sum_abs=0.0;
                    SizeType band_lt_1e14=0u;
                    SizeType band_1e14_1e12=0u;
                    SizeType band_1e12_1e10=0u;
                    SizeType band_ge_1e10=0u;
                    for(auto const& term : expansion) {
                        if(term.index().degree()!=degree) { continue; }
                        auto a=cast_exact(abs(term.coefficient())).get_d();
                        ++count;
                        sum_abs+=a;
                        if(a>max_abs) { max_abs=a; }
                        if(a<1e-14) { ++band_lt_1e14; }
                        else if(a<1e-12) { ++band_1e14_1e12; }
                        else if(a<1e-10) { ++band_1e12_1e10; }
                        else { ++band_ge_1e10; }
                    }
                    if(count!=0u) {
                        std::cerr << "[CarriedExpansionDegreeProfile]"
                                  << " policy=" << policy
                                  << " component=" << i
                                  << " degree=" << degree
                                  << " count=" << count
                                  << " max_abs=" << max_abs
                                  << " sum_abs=" << sum_abs
                                  << " lt_1e14=" << band_lt_1e14
                                  << " ge_1e14_lt_1e12=" << band_1e14_1e12
                                  << " ge_1e12_lt_1e10=" << band_1e12_1e10
                                  << " ge_1e10=" << band_ge_1e10
                                  << std::endl;
                    }
                }
            }
        };

    run_structure_probe(
        "absolute_1e-12",
        Sweeper<FloatDP>(ThresholdSweeper<FloatDP>(DoublePrecision(),1e-12)));
    run_structure_probe(
        "absolute_1e-14",
        Sweeper<FloatDP>(ThresholdSweeper<FloatDP>(DoublePrecision(),1e-14)));

}
