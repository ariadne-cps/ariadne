/***************************************************************************
 *            attractor.cpp
 *
 *  Copyright  2017-20  Luca Geretti, Pieter Collins
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

void ariadne_main()
{
    RealVariable x("x"), y("y");
    VectorField system = {{dot(x)=2*x-x*y,dot(y)=2*x*x-y}};
    RealExpressionBoundedConstraintSet initial_set = {{0.9_dec<=x<=1,-2.2_dec<=y<=-2},{sqr(x)+sqr(y+2)<=1}};
    RealExpressionBoundedConstraintSet safe_set = {{-1<=x<=4,-4<=y<=6},{sqr(x-2)+sqr(y-1)<=22}};
    CONCLOG_PRINTLN_VAR(system);
    CONCLOG_PRINTLN_VAR(initial_set);
    CONCLOG_PRINTLN_VAR(safe_set);

    auto initial_constraint_set = initial_set.euclidean_set(system.state_space());
    auto safe_constraint_set = safe_set.euclidean_set(system.state_space());
    CONCLOG_PRINTLN_VAR(initial_constraint_set);
    CONCLOG_PRINTLN_VAR(safe_constraint_set);

    Real evolution_time = 50;

    VectorFieldSimulator simulator(system);
    simulator.configuration().set_step_size(0.1);
    CONCLOG_PRINTLN("Simulating...");
    auto orbit = simulator.orbit(initial_set,evolution_time);
    LabelledFigure g(Axes2d{{-2<=x<=5},{-4<=y<=6}});
    g << orbit;
    g.write("attractor_simulation");

    TaylorPicardIntegrator integrator(0.01);
    CONCLOG_PRINTLN("Evolving...");
    VectorFieldEvolver evolver(system,integrator);
    evolver.configuration().set_maximum_step_size(0.1);
    auto evolver_orbit = evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    g.clear();
    g << evolver_orbit;
    g.write("attractor_evolution");

    ContinuousReachabilityAnalyser analyser(evolver);
    analyser.configuration().set_transient_time(0.75_dec);
    analyser.configuration().set_lock_to_grid_time(0.75_dec);
    analyser.configuration().set_maximum_grid_extent(5);

    CONCLOG_PRINTLN("Computing safety...");
    auto safety = analyser.verify_safety(initial_constraint_set,safe_constraint_set);
    CONCLOG_PRINTLN_VAR(safety.is_safe);
    g.clear();
    g << fill_colour(lightgrey) << safety.safe_set
      << fill_colour(orange) << safety.chain_reach_set;
    g.write("attractor_chain_reach");
}
