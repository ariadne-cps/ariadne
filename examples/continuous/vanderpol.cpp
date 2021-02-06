/***************************************************************************
 *            vanderpol.cpp
 *
 *  Copyright  2017-20  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "ariadne.hpp"

using namespace Ariadne;

int main(int argc, const char* argv[])
{
    ARIADNE_LOG_SET_VERBOSITY(get_verbosity(argc,argv));
    Logger::instance().configuration().set_theme(TT_THEME_DARK);
    Logger::instance().configuration().set_thread_name_printing_policy(ThreadNamePrintingPolicy::BEFORE);
    ConcurrencyManager::instance().set_concurrency(4);
    Logger::instance().use_blocking_scheduler();

    ARIADNE_LOG_PRINTLN("van der Pol oscillator");

    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");

    VectorField system({dot(x)=y, dot(y)=mu*y*(1-sqr(x))-x});

    double max_err = 1e-8;
    auto sweeper1 = ThresholdSweeper<FloatDP>(DoublePrecision(),max_err/10);
    auto sweeper2 = ThresholdSweeper<FloatDP>(DoublePrecision(),max_err/100);
    auto sweeper3 = ThresholdSweeper<FloatDP>(DoublePrecision(),max_err/1000);

    TaylorPicardIntegrator integrator(Configuration<TaylorPicardIntegrator>()
                                          .set_step_maximum_error(1e-8,1e-5)
                                          .set_maximum_temporal_order(8,15)
                                          .set_lipschitz_tolerance(1e-2,0.5)
                                          .set_sweeper({sweeper1,sweeper2,sweeper3})
                                          );

    typedef VectorFieldEvolver E; typedef TaskInput<E> I; typedef TaskOutput<E> O; typedef TaskObjective<E> OBJ;

    E evolver(system,Configuration<E>().set_integrator(integrator));
    ARIADNE_LOG_PRINTLN_VAR_AT(1,evolver.configuration());
    ARIADNE_LOG_PRINTLN_VAR_AT(1,evolver.configuration().search_space());

    double factor = 0.1;
    OBJ y_p275(y,PositiveFloatDPUpperBound(FloatDP(cast_exact(0.07*factor),DoublePrecision())),Dyadic(6.48_x));
    OBJ y_m275(y,PositiveFloatDPUpperBound(FloatDP(cast_exact(0.075*factor),DoublePrecision())),Dyadic(3.15_x));
    auto verification_parameter_p275 = ScalarRankingParameter<E>(y.name(), OptimisationCriterion::MINIMISE, [y](I const& i, O const& o, DurationType const& d) {
        return o.reach.bounding_box()[y].upper_bound().get_d(); });
    auto verification_parameter_m271 = ScalarRankingParameter<E>(y.name(), OptimisationCriterion::MAXIMISE, [y](I const& i, O const& o, DurationType const& d) {
        return o.reach.bounding_box()[y].lower_bound().get_d(); });
    auto verification_constraint_p275 = TaskRankingConstraint<E>(verification_parameter_p275, 2.75, RankingConstraintSeverity::CRITICAL);
    auto verification_constraint_m271 = TaskRankingConstraint<E>(verification_parameter_m271, -2.75, RankingConstraintSeverity::CRITICAL);
    auto refinement_target = ConfigurationPropertyRefinement<E>(
            ConfigurationPropertyPath("integrator").append("step_maximum_error"),{y_p275,y_m275},SaturateIfPositiveProportionalRefiner(-1e-3));
    VerificationManager::instance().add_safety_specification(evolver,
                                                             {verification_constraint_m271,verification_constraint_p275},
                                                             {refinement_target});

    Real x0 = 1.4_dec;
    Real y0 = 2.4_dec;
    Real eps_x0 = 0.15_dec;
    Real eps_y0 = 0.05_dec;

    auto function_factory = TaylorFunctionFactory(sweeper1);
    EnclosureConfiguration enclosure_config(function_factory);
    enclosure_config.set_reconditioning_num_blocks(3);
    auto initial_set = evolver.enclosure({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0},enclosure_config);
    ARIADNE_LOG_PRINTLN_VAR_AT(1,initial_set);

    Real evolution_time = 7;

    auto start = std::chrono::high_resolution_clock::now();
    ARIADNE_LOG_PRINTLN("Computing orbit... ");
    try {
        auto orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
        auto end = std::chrono::high_resolution_clock::now();
        ARIADNE_LOG_PRINTLN_AT(1,"Done in " << ((double)std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count())/1000 << " seconds.");

        ARIADNE_LOG_PRINTLN_AT(1,"Optimal point: " << ConcurrencyManager::instance().optimal_point());

        ARIADNE_LOG_PRINTLN("Plotting...");
        LabelledFigure fig({-2.5<=x<=2.5,-3<=y<=3});
        fig << fill_colour(1.0,0.75,0.5);
        fig.draw(orbit.reach());
        fig.write("vanderpol");
    } catch (CriticalRankingFailureException<VectorFieldEvolver>& ex) {
        ARIADNE_LOG_PRINTLN("Safety verification failure: " << ex.what());
    }

    ConcurrencyManager::instance().print_best_rankings();
    ConcurrencyManager::instance().print_refinement_values();
}
