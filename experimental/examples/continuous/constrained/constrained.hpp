/***************************************************************************
 *            constrained_vanderpol.cpp
 *
 *  Copyright  2023  Luca Geretti
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

#include "ariadne_main.hpp"
#include "helper/stopwatch.hpp"

using namespace std;
using namespace pExplore;

LabelledUpperBoxType bounding_box(ListSet<LabelledEnclosure> const& es) {
    HELPER_PRECONDITION(not es.empty())
    UpperBoxType result = es[0].euclidean_set().bounding_box();
    for (auto const& e : es) {
        auto const& bx = e.euclidean_set().bounding_box();
        for (size_t j=0; j<result.dimension(); ++j)
            result[j] = hull(result[j],bx[j]);
    }
    return {es[0].state_space(),result};
}

void constrained_execution(pExplore::String const& name, VectorField const& dynamics, List<RealExpression> const& constraints, RealExpressionBoundedConstraintSet const& initial_set, Real const& evolution_time) {
    auto configuration = Configuration<VectorFieldEvolver>().
            set_integrator(TaylorPicardIntegrator(Configuration<TaylorPicardIntegrator>()
                    .set_step_maximum_error(1e-6,1e-4)
                    .set_sweeper(ThresholdSweeperDP(dp,Configuration<ThresholdSweeperDP>().set_threshold(1e-8,1e-6)))
                    .set_minimum_temporal_order(0,5)
                    .set_maximum_temporal_order(5,10)
                    .set_bounder(EulerBounder(Configuration<EulerBounder>().set_lipschitz_tolerance(0.0675,0.5)))
                    ));

    CONCLOG_PRINTLN("Dynamics: " << dynamics)
    CONCLOG_PRINTLN("Constraints: " << constraints)
    CONCLOG_PRINTLN_VAR(configuration)
    CONCLOG_PRINTLN_VAR_AT(1,configuration.search_space())
    CONCLOG_PRINTLN_AT(1,"Total points: " << configuration.search_space().total_points())

    auto result = constrained_evolution(dynamics,initial_set,evolution_time,constraints,configuration);

    CONCLOG_PRINTLN("Constraint satisfaction result:" << result.satisfaction)

    auto best_scores = pExplore::TaskManager::instance().best_scores();
    CONCLOG_PRINTLN_AT(2,"Best scores:")
    for (auto const& b : best_scores) {
        CONCLOG_PRINTLN_AT(2,b)
    }

    CONCLOG_PRINTLN_AT(1,"Optimal point: " << pExplore::TaskManager::instance().optimal_point())

    auto variable_names = dynamics.state_space().variable_names();
    auto drawing_box = bounding_box(result.rigorous.reach());

    CONCLOG_PRINTLN("Plotting...")
    for (size_t i=0; i<dynamics.dimension()-1; i++) {
        auto xi = RealVariable(variable_names.at(i));
        for (size_t j=1; j<dynamics.dimension(); j++) {
            auto xj = RealVariable(variable_names.at(j));
            LabelledFigure fig=LabelledFigure({drawing_box[xi].lower_bound().get_d()<=xi<=drawing_box[xi].upper_bound().get_d(),drawing_box[xj].lower_bound().get_d()<=xj<=drawing_box[xj].upper_bound().get_d()});
            fig.clear();
            fig.draw(result.approximate);
            char var_char[64] = "";
            if (dynamics.dimension() > 2) snprintf(var_char,64,"[%s,%s]",xi.name().c_str(),xj.name().c_str());
            CONCLOG_RUN_MUTED(fig.write((name+"_approximate"+var_char).c_str()))
            fig.clear();
            fig.draw(result.rigorous);
            CONCLOG_RUN_MUTED(fig.write((name+"_rigorous"+var_char).c_str()))
            //auto approximator = NonlinearCandidateValidationInnerApproximator(ParallelLinearisationContractor(NativeSimplex(),2,1));
            //auto inner = approximator.compute_from(result.constrained.intermediate());
            //fig << fill_colour(red);
            //fig.draw(inner);
            fig.clear();
            fig.draw(result.constrained);
            CONCLOG_RUN_MUTED(fig.write((name+"_constrained"+var_char).c_str()))
        }
    }
}
