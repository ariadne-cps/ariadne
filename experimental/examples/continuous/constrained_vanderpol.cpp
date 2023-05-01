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

void ariadne_main()
{
    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    RealConstant ymax("ymax",2.8_x);
    RealConstant ymin("ymin",-3.0_x);
    RealConstant xmin("xmin",-1.0_x);
    RealConstant xmax("xmax",2.1_x);
    RealConstant rsqr("r^2",2.0_x);
    List<RealExpression> constraints = {y - ymin, x - xmin, ymax - y, xmax - x, sqr(x) + sqr(y) - rsqr};

    auto configuration = Configuration<VectorFieldEvolver>().
            set_maximum_enclosure_radius(1.0).
            set_both_enable_reconditioning().
            set_integrator(TaylorPicardIntegrator(Configuration<TaylorPicardIntegrator>().set_step_maximum_error(1e-6,1e-5))).
            set_maximum_step_size(0.005,1.0);
    CONCLOG_PRINTLN_VAR(configuration)
    CONCLOG_PRINTLN_VAR_AT(1,configuration.search_space())

    Real x0 = 1.40_dec;
    Real y0 = 2.40_dec;
    Real eps_x0 = 0.15_dec;
    Real eps_y0 = 0.05_dec;
    RealExpressionBoundedConstraintSet initial_set({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0});

    Real evolution_time = 7;

    auto result = constrained_evolution(dynamics,initial_set,evolution_time,constraints,configuration);

    auto const& approximate_orbit = get<0>(result);
    auto const& rigorous_orbit = get<1>(result);
    auto const& outcomes = get<2>(result);

    CONCLOG_PRINTLN("Constraint checking outcomes:")
    for (size_t m=0; m<constraints.size(); ++m) {
        CONCLOG_PRINTLN(constraints.at(m) << " -> " << outcomes.at(m))
    }

    auto best_scores = pExplore::TaskManager::instance().best_scores();
    CONCLOG_PRINTLN_AT(2,"Best scores:")
    for (auto const& b : best_scores) {
        CONCLOG_PRINTLN_AT(2,b)
    }

    LabelledFigure fig=LabelledFigure({-3<=x<=3,-3<=y<=3});
    CONCLOG_PRINTLN("Plotting...")
    fig.clear();
    fig.draw(approximate_orbit);
    CONCLOG_RUN_MUTED(fig.write("vanderpol_approximate"))
    fig.draw(rigorous_orbit);
    CONCLOG_RUN_MUTED(fig.write("vanderpol_rigorous"))
}
