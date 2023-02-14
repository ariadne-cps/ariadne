/***************************************************************************
 *            heading_control.cpp
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
#include "reach_avoid.hpp"
#include "utility/stopwatch.hpp"

void ariadne_main()
{
    Real deltat=0.1_dec, v=3;
    RealVariable x("x"), y("y"), theta("theta"), u("u");
    IteratedMap heading({next(x)=x+deltat*v*cos(theta),next(y)=y+deltat*v*sin(theta),next(theta)= theta+u,next(u)=u});

    auto dynamics = heading.function().zeros(3,4);
    for (SizeType i=0; i<3; ++i)
        dynamics[i] = heading.function().get(i);
    CONCLOG_PRINTLN_VAR(dynamics);

    double pi_ = pi.get_d();

    Grid state_grid({0.5,0.5,2*pi/8});
    BoundType theta_domain = {-2*pi_,2*pi_};
    BoundsBoxType state_domain({{0,5},{0,5},theta_domain});
    Grid control_grid({pi/4});
    BoundsBoxType control_domain({{-pi_,pi_}});
    SizeType depth = 0;

    ReachAvoid scs("heading", dynamics, state_grid, state_domain, control_grid, control_domain, depth, 1e-10_x);

    CONCLOG_PRINTLN_VAR_AT(1,scs.state_size())
    CONCLOG_PRINTLN_VAR_AT(1,scs.controller_size())

    scs.add_obstacle({{1,3.5},{4.5,5},theta_domain});
    scs.add_obstacle({{0,1},{2,3},theta_domain});
    scs.add_obstacle({{2.5,5},{2,3},theta_domain});
    scs.add_obstacle({{0,5},{0,0.5},theta_domain});
    CONCLOG_RUN_AT(2,scs.print_obstacles())
    CONCLOG_PRINTLN_VAR_AT(1,scs.obstacles_size())

    scs.add_goal({{4,5},{4.5,5},theta_domain});

    CONCLOG_RUN_AT(2,scs.print_goals())

    CONCLOG_PRINTLN_VAR_AT(1,scs.goals_size())
    CONCLOG_PRINTLN_VAR_AT(1,scs.unverified_size())

    Stopwatch<Milliseconds> sw;
    scs.compute_reachability_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of constructing reachability graph: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1,scs.num_transitions())

    sw.restart();
    scs.refine_to_safety_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of reducing graph to safe one: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1,scs.num_transitions())

    CONCLOG_PRINTLN_AT(1,"Safe abstract states: " << scs.num_sources())

    sw.restart();
    scs.refine_to_goal_reachable_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of reducing graph to goal-reachable one: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1,scs.num_transitions())

    CONCLOG_PRINTLN_AT(1,"Safe goal-reachable abstract states: " << scs.num_sources())

    scs.update_unverified();
    CONCLOG_PRINTLN_AT(1,"Unverified abstract states: " << scs.unverified_size() << " (" << scs.unverified_percentage() << "% left)")

    scs.plot({{0,5},{0,5},{-6.28_x,6.28_x}},0,1);
    scs.plot({{0,5},{0,5},{-6.28_x,6.28_x}},0,2);
    scs.plot({{0,5},{0,5},{-6.28_x,6.28_x}},1,2);

    CONCLOG_RUN_AT(3,scs.print_graph())
}