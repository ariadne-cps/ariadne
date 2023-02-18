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

Tuple<IteratedMap,Grid,BoundsBoxType> u_control() {
    Real deltat=0.1_dec, v=3;
    RealVariable x("x"), y("y"), theta("theta"), u("u");
    IteratedMap heading({next(x)=x+deltat*v*cos(theta),next(y)=y+deltat*v*sin(theta),next(theta)= theta+u,next(u)=u});

    Grid control_grid({pi/4});
    double pi_ = pi.get_d();
    BoundsBoxType control_domain({{-pi_-pi_/4,pi_+pi_/4}});

    return Tuple<IteratedMap,Grid,BoundsBoxType>(heading,control_grid,control_domain);
}

Tuple<IteratedMap,Grid,BoundsBoxType> cpwa_control() {
    Real deltat = 0.1_dec, v = 3;
    RealVariable x("x"), y("y"), theta("theta"), K1("K1"), K2("K2"), K3("K3"), b("b");
    IteratedMap heading({next(x)=x+deltat*v*cos(theta),next(y)=y+deltat*v*sin(theta),next(theta)= theta+deltat*(K1*x+K2*y+K3*theta+b),
                         next(K1)=K1,next(K2)=K2,next(K3)=K3,next(b)=b});

    double pi_ = pi.get_d();
    Grid control_grid({1,1,1,24*pi_/20});
    BoundsBoxType control_domain({{-1,1},{-1,1},{-1,1},{-12*pi_,12*pi_}});

    return Tuple<IteratedMap,Grid,BoundsBoxType>(heading,control_grid,control_domain);
}

ReachAvoid set_workspace_1(EffectiveVectorMultivariateFunction const& dynamics, Grid const& control_grid, BoundsBoxType const& control_domain) {
    double pi_ = pi.get_d();

    Grid state_grid({0.5,0.5,2*pi/8});
    BoundType theta_domain = {0,2*pi_};
    BoundsBoxType state_domain({{0,5},{0,5},theta_domain});
    SizeType depth = 0;

    ReachAvoid ra("heading", dynamics, state_grid, state_domain, control_grid, control_domain, depth, 1e-10_x);

    ra.add_obstacle({{1,3.5},{4.5,5},theta_domain});
    ra.add_obstacle({{0,1},{2,3},theta_domain});
    ra.add_obstacle({{2.5,5},{2,3},theta_domain});
    ra.add_obstacle({{0,5},{0,0.5},theta_domain});

    ra.add_goal({{4,5},{4.5,5},theta_domain});

    return ra;
}

ReachAvoid set_workspace_2(EffectiveVectorMultivariateFunction const& dynamics, Grid const& control_grid, BoundsBoxType const& control_domain) {

    double pi_ = pi.get_d();

    Grid state_grid({0.25,0.25,2*pi/8});
    BoundType theta_domain = {0,2*pi_};
    BoundsBoxType state_domain({{0,5},{0,5},theta_domain});
    SizeType depth = 0;

    ReachAvoid ra("heading", dynamics, state_grid, state_domain, control_grid, control_domain, depth, 1e-10_x);

    ra.add_obstacle({{0,0.5},{0,5},theta_domain});
    ra.add_obstacle({{1.5,2},{0,3},theta_domain});
    ra.add_obstacle({{1.5,2},{4,5},theta_domain});
    ra.add_obstacle({{3,3.5},{0,1},theta_domain});
    ra.add_obstacle({{3,3.5},{2,5},theta_domain});
    ra.add_obstacle({{4.5,5},{0,5},theta_domain});

    ra.add_goal({{3.75,4.25},{4.5,5},theta_domain});

    return ra;
}

void ariadne_main()
{
    auto sys = u_control();

    auto dynamics = get<0>(sys).function().zeros(3,get<0>(sys).dimension());
    for (SizeType i=0; i<3; ++i)
        dynamics[i] = get<0>(sys).function().get(i);
    CONCLOG_PRINTLN_VAR(dynamics);

    auto ra = set_workspace_1(dynamics,get<1>(sys), get<2>(sys));

    CONCLOG_PRINTLN_VAR_AT(1,ra.state_size())
    CONCLOG_PRINTLN_VAR_AT(1,ra.controller_size())

    CONCLOG_RUN_AT(2,ra.print_obstacles())
    CONCLOG_PRINTLN_VAR_AT(1,ra.obstacles_size())

    CONCLOG_RUN_AT(2,ra.print_goals())
    CONCLOG_PRINTLN_VAR_AT(1,ra.goals_size())

    CONCLOG_PRINTLN_VAR_AT(1,ra.unverified_size())

    Stopwatch<Milliseconds> sw;
    ra.compute_reachability_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of constructing reachability graph: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1,ra.num_transitions())

    sw.restart();
    ra.refine_to_safety_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of reducing graph to safe one: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1,ra.num_transitions())

    CONCLOG_PRINTLN_AT(1,"Safe abstract states: " << ra.num_sources())

    sw.restart();
    ra.refine_to_goal_reachable_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of reducing graph to goal-reachable one: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1,ra.num_transitions())

    CONCLOG_PRINTLN_AT(1,"Safe goal-reachable abstract states: " << ra.num_sources())

    ra.update_unverified();
    CONCLOG_PRINTLN_AT(1,"Unverified abstract states: " << ra.unverified_size() << " (" << ra.unverified_percentage() << "% left)")

    ra.plot({{0,5},{0,5},{-6.28_x,6.28_x}},0,1);
    ra.plot({{0,5},{0,5},{-6.28_x,6.28_x}},0,2);
    ra.plot({{0,5},{0,5},{-6.28_x,6.28_x}},1,2);

    CONCLOG_RUN_AT(3,ra.print_graph())
}