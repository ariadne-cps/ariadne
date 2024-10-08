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
#include "verification/reach_avoid.hpp"
#include "verification/reach_avoid_strategy.hpp"
#include "utility/stopwatch.hpp"
#include "utility/randomiser.hpp"

Tuple<IteratedMap,Grid,RealBox> u_control() {
    Real deltat=0.1_dec, v=3;
    RealVariable x("x"), y("y"), theta("theta"), u("u");
    IteratedMap heading({next(x)=x+deltat*v*cos(theta),next(y)=y+deltat*v*sin(theta),next(theta)= theta+u,next(u)=u});

    Grid control_grid({pi/4});
    RealBox control_domain({{-pi-pi/4,pi+pi/4}});

    return {heading,control_grid,control_domain};
}

Tuple<IteratedMap,Grid,RealBox> cpwa_control() {
    Real deltat = 0.1_dec, v = 3;
    RealVariable x("x"), y("y"), theta("theta"), K1("K1"), K2("K2"), K3("K3"), b("b");
    IteratedMap heading({next(x)=x+deltat*v*cos(theta),next(y)=y+deltat*v*sin(theta),next(theta)= theta+deltat*(K1*x+K2*y+K3*theta+b),
                         next(K1)=K1,next(K2)=K2,next(K3)=K3,next(b)=b});

    Grid control_grid({1,1,1,24*pi/20});
    RealBox control_domain({{-1,1},{-1,1},{-1,1},{-12*pi,12*pi}});

    return {heading,control_grid,control_domain};
}

ReachAvoid set_workspace_1(EffectiveVectorMultivariateFunction const& dynamics, Grid const& control_grid, RealBox const& control_domain) {

    Grid state_grid({0.5,0.5,2*pi/8});
    RealInterval theta_domain = {0, 2*pi};
    RealBox state_domain({{0,5},{0,5},theta_domain});
    SizeType depth = 0;

    ReachAvoid ra("heading", dynamics, state_grid, state_domain, control_grid, control_domain, depth, 1e-10_x, 1e-4);

    ra.add_obstacle({{1,3.5_x},{4.5_x,5},theta_domain});
    ra.add_obstacle({{0,1},{2,3},theta_domain});
    ra.add_obstacle({{2.5_x,5},{2,3},theta_domain});
    ra.add_obstacle({{0,5},{0,0.5_x},theta_domain});

    ra.add_goal({{4,5},{4.5_x,5},theta_domain});

    return ra;
}

ReachAvoid set_workspace_2(EffectiveVectorMultivariateFunction const& dynamics, Grid const& control_grid, RealBox const& control_domain) {

    Grid state_grid({0.25,0.25,2*pi/8});
    RealInterval theta_domain = {0, 2 * pi};
    RealBox state_domain({{0,5},{0,5},theta_domain});
    SizeType depth = 0;

    ReachAvoid ra("heading", dynamics, state_grid, state_domain, control_grid, control_domain, depth, 1e-10_x, 1e-4);

    ra.add_obstacle({{0,0.5_x},{0,5},theta_domain});
    ra.add_obstacle({{1.5_x,2},{0,3},theta_domain});
    ra.add_obstacle({{1.5_x,2},{4,5},theta_domain});
    ra.add_obstacle({{3,3.5_x},{0,1},theta_domain});
    ra.add_obstacle({{3,3.5_x},{2,5},theta_domain});
    ra.add_obstacle({{4.5_x,5},{0,5},theta_domain});

    ra.add_goal({{3.75_x,4.25_x},{4.5_x,5},theta_domain});

    return ra;
}

void check_scalabilities(SizeType n) {

    EffectiveVectorMultivariateFunction u_dynamics = EffectiveVectorMultivariateFunction::zeros(n,n+1);
    for (SizeType i=0; i<n; ++i)
        u_dynamics[i] = u_dynamics.coordinate(n+1,i) + u_dynamics.coordinate(n+1,n);
    CONCLOG_PRINTLN_VAR(u_dynamics)

    EffectiveVectorMultivariateFunction cpwa_dynamics = EffectiveVectorMultivariateFunction::zeros(n,2*n+1);
    for (SizeType i=0; i<n; ++i) {
        cpwa_dynamics[i] = cpwa_dynamics.coordinate(2*n+1,i);
    }
    for (SizeType j=n; j<2*n+1; ++j)
        cpwa_dynamics[0] = cpwa_dynamics[0] + cpwa_dynamics.coordinate(2*n+1,j)/(n+1);
    CONCLOG_PRINTLN_VAR(cpwa_dynamics)

    Grid state_grid(n,1);
    RealBox state_domain(n,{0,2});

    Grid u_control_grid(1,1);
    RealBox u_control_domain(1,{0,2});

    ReachAvoid u_ra("u_ra",u_dynamics,state_grid,state_domain,u_control_grid,u_control_domain,0,1e-10_x,1e-4);

    CONCLOG_PRINTLN_AT(1,"State size: " << u_ra.state_size())
    CONCLOG_PRINTLN_VAR_AT(1,u_ra.control_size())

    Stopwatch<Milliseconds> sw;
    u_ra.compute_free_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"(u) Time cost of constructing free graph: " << sw.elapsed_seconds() << " seconds")
    sw.restart();

    Grid cpwa_control_grid(n+1,1);
    RealBox cpwa_control_domain(n+1,{0,2});
    ReachAvoid cpwa_ra("cpwa_ra",cpwa_dynamics,state_grid,state_domain,cpwa_control_grid,cpwa_control_domain,0,1e-10_x,1e-4);
    CONCLOG_PRINTLN_VAR_AT(1,cpwa_ra.control_size())

    cpwa_ra.compute_free_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"(CPWA) Time cost of constructing free graph: " << sw.elapsed_seconds() << " seconds")
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
    CONCLOG_PRINTLN_VAR_AT(1,ra.control_size())

    CONCLOG_RUN_AT(2,ra.print_obstacles())
    CONCLOG_PRINTLN_VAR_AT(1,ra.obstacles_size())

    CONCLOG_RUN_AT(2,ra.print_goals())
    CONCLOG_PRINTLN_VAR_AT(1,ra.goals_size())

    CONCLOG_PRINTLN_VAR_AT(1,ra.unverified_size())

    Stopwatch<Milliseconds> sw;
    ra.compute_free_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of constructing free graph: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1, ra.unconstrained_num_transitions())

    sw.restart();
    ra.compute_avoid_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of constructing avoid graph: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1, ra.avoiding_num_transitions())

    CONCLOG_PRINTLN_AT(1,"Safe abstract states: " << ra.num_sources())

    sw.restart();
    ra.compute_possibly_reaching_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of constructing possibly reaching graph: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1,ra.possibly_reaching_num_transitions())

    CONCLOG_PRINTLN_AT(1,"Safe goal-reachable abstract states: " << ra.num_sources())

    ra.update_unverified();
    CONCLOG_PRINTLN_AT(1,"Unverified abstract states: " << ra.unverified_size() << " (" << ra.unverified_percentage() << "% left)")

    ra.plot();

    ReachAvoidStrategyBuilder strategy_builder(ra.dynamics(),ra.possibly_reaching_graph());
    auto assignments = strategy_builder.build().assignments();

    Randomiser<DegreeType> random_uint;
    Randomiser<ExactDouble> random_double;

    SizeType idx = static_cast<SizeType>(random_uint.get(0,assignments.size()-1));
    auto a_it = assignments.begin();
    auto dim = a_it->first.cell().dimension();
    for (SizeType i=0;i<idx;++i) ++a_it;
    auto initial_box = a_it->first.cell().box();

    CONCLOG_PRINTLN_AT(1,"Starting from box " << initial_box << " using control " << a_it->second.control().cell().box())

    PointType current(dim);
    for (SizeType i=0; i<dim; ++i) {
        current.at(i) = random_double.get(cast_exact(initial_box[i].lower_bound().get_d()),cast_exact(initial_box[i].upper_bound().get_d())).get_d();
    }

    auto const& vertex_factory = ra.vertex_factory();

    SizeType step = 0;

    List<PointType> sequence;

    while(true) {

        sequence.append(current);

        ExactBoxType boxed_current(dim);
        for (SizeType i=0; i<dim; ++i) {
            boxed_current[i] = ExactIntervalType(cast_exact(current[i]),cast_exact(current[i]+1e-10));
        }

        SPaving point_paving(ra.state_grid());
        point_paving.adjoin_over_approximation(boxed_current,ra.grid_depth());

        auto point_cell = *point_paving.begin();
        auto icell = vertex_factory.create(point_cell);

        auto target = assignments.get(icell).target_point();

        double delta_modulus = 0;
        for (SizeType i=0; i<dim; ++i) {
            delta_modulus += (target.at(i)-current.at(i))*(target.at(i)-current.at(i));
        }
        delta_modulus = std::sqrt(delta_modulus);

        PointType next(dim);
        for (SizeType i=0; i<dim; ++i) {
            next.at(i) = current.at(i) + (target.at(i)-current.at(i))/delta_modulus*0.3;
        }

        CONCLOG_PRINTLN_AT(1, step << ": from " << current << " (" << icell.id() << ") to " << next << " under control " << assignments.get(icell).control().cell().box())

        current = next;

        ++step;

        if (step > 50) break;
    }

    /*
    std::this_thread::sleep_for(std::chrono::milliseconds(500));

    std::cout << "sequence = [";
    for (auto const& p : sequence) {
        std::cout << p.at(0) << " " << p.at(1) << " " << p.at(2) << ";" << std::endl;
    }
    std::cout << "];";
    */
}