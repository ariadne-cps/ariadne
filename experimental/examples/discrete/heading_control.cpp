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

#include <iostream>
#include <fstream>
#include "ariadne_main.hpp"
#include "verification/reach_avoid.hpp"
#include "verification/reach_avoid_strategy.hpp"
#include "utility/stopwatch.hpp"
#include "utility/randomiser.hpp"

Tuple<IteratedMap,Grid,RealBox> u_control() {
    Real deltat=0.1_dec, v=3;
    RealVariable x("x"), y("y"), theta("theta"), u("u");
    IteratedMap heading({next(x)=x+deltat*v*cos(theta),next(y)=y+deltat*v*sin(theta),next(theta)= theta+deltat*u,next(u)=u});

    Grid control_grid({pi/4*10});
    RealBox control_domain({{-pi*10,pi*10}});

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

    ReachAvoid ra("heading", dynamics, state_grid, state_domain, control_grid, control_domain, depth, 1e-10_x);

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

    ReachAvoid ra("heading", dynamics, state_grid, state_domain, control_grid, control_domain, depth, 1e-10_x);

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

    ReachAvoid u_ra("u_ra",u_dynamics,state_grid,state_domain,u_control_grid,u_control_domain,0,1e-10_x);

    CONCLOG_PRINTLN_AT(1,"State size: " << u_ra.state_size())
    CONCLOG_PRINTLN_VAR_AT(1,u_ra.control_size())

    Stopwatch<Milliseconds> sw;
    u_ra.compute_bounded_domain_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"(u) Time cost of constructing free graph: " << sw.elapsed_seconds() << " seconds")
    sw.restart();

    Grid cpwa_control_grid(n+1,1);
    RealBox cpwa_control_domain(n+1,{0,2});
    ReachAvoid cpwa_ra("cpwa_ra",cpwa_dynamics,state_grid,state_domain,cpwa_control_grid,cpwa_control_domain,0,1e-10_x);
    CONCLOG_PRINTLN_VAR_AT(1,cpwa_ra.control_size())

    cpwa_ra.compute_bounded_domain_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"(CPWA) Time cost of constructing free graph: " << sw.elapsed_seconds() << " seconds")
}

IdentifiedCell point_to_cell(PointType const& pt, Grid const& grid, SizeType grid_depth, IdentifiedCellFactory const& vertex_factory) {

    auto dim = pt.size();

    ExactBoxType boxed_current(dim);
    for (SizeType i=0; i<dim; ++i) {
        boxed_current[i] = ExactIntervalType(cast_exact(pt[i]),cast_exact(pt[i]+1e-10));
    }

    SPaving point_paving(grid);
    point_paving.adjoin_over_approximation(boxed_current,grid_depth);

    auto point_cell = *point_paving.begin();

    return vertex_factory.create(point_cell);
}

void create_workspace_1_matlab_file(SizeType idx, List<PointType> const& sequence) {

    std::stringstream strstr;
    strstr << "scatter" << idx << ".m";
    std::ofstream file(strstr.str());

    file << "clear all; close all;" << std::endl;

    file << "x = [";
    for (auto pt : sequence) {
        file << pt[0] << " ";
    }
    file << "];" << std::endl;

    file << "y = [";
    for (auto pt : sequence) {
        file << pt[1] << " ";
    }
    file << "];" << std::endl;

    file << "theta = [";
    for (auto pt : sequence) {
        file << pt[2] << " ";
    }
    file << "];" << std::endl;
    file << "axis equal;" << std::endl;
    file << "xlim([0 5]); ylim([0 5]); zlim([0 2*pi]);" << std::endl;
    file << "xlabel('X'); ylabel('Y'); zlabel('Theta');" << std::endl;
    file << "figure(1);" << std::endl;
    file << "hold on;" << std::endl;
    file << "plot3(x(1),y(1),theta(1),'o-');" << std::endl;
    file << "plot3(x,y,theta,'r-');" << std::endl;
    file.close();
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

    CONCLOG_PRINTLN_VAR_AT(1,ra.num_sources())

    Stopwatch<Milliseconds> sw;
    ra.compute_bounded_domain_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of constructing bounded domain graph: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1,ra.num_sources())
    CONCLOG_PRINTLN_VAR_AT(1,ra.bounded_domain_graph().num_transitions())
    CONCLOG_PRINTLN_VAR_AT(1,ra.bounded_domain_graph().is_consistent())

    sw.restart();
    ra.compute_avoiding_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of constructing avoiding graph: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1,ra.num_sources())
    CONCLOG_PRINTLN_VAR_AT(1,ra.safe_goals().size())
    CONCLOG_PRINTLN_VAR_AT(1,ra.avoiding_graph().num_transitions())
    CONCLOG_PRINTLN_VAR_AT(1,ra.avoiding_graph().is_consistent())

    sw.restart();
    ra.compute_possibly_reaching_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of constructing possibly reaching graph: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1,ra.num_sources())
    CONCLOG_PRINTLN_VAR_AT(1,ra.possibly_reaching_graph().num_transitions())
    CONCLOG_PRINTLN_VAR_AT(1,ra.possibly_reaching_graph().is_consistent())

    CONCLOG_PRINTLN_AT(1,"Unverified abstract states: " << ra.unverified_size() << " (" << ra.unverified_percentage() << "% left)")

    ra.plot();

    ReachAvoidStrategyBuilder strategy_builder(ra.dynamics(),ra.possibly_reaching_graph());

    auto assignments = strategy_builder.build().assignments();

    Randomiser<DegreeType> random_uint;
    Randomiser<ExactDouble> random_double;

    SizeType num_points = 5;
    SizeType max_steps = 100;

    return;

    for (SizeType p=0; p < num_points; ++p) {

        SizeType idx = static_cast<SizeType>(random_uint.get(0,assignments.size()-1));
        auto a_it = assignments.begin();
        auto dim = a_it->first.cell().dimension();
        for (SizeType i=0;i<idx;++i) ++a_it;
        auto current_cell_box = a_it->first.cell().box();

        PointType current(dim);
        for (SizeType i=0; i<dim; ++i) {
            auto rand_value = random_double.get(cast_exact(current_cell_box[i].lower_bound().get_d()),cast_exact(current_cell_box[i].upper_bound().get_d()));
            current[i] = rand_value.get_d();
        }

        /*
        current[0] = 0.406792;
        current[1] = 3.94616;
        current[2] = 2.45582;
        */

        CONCLOG_PRINTLN("Point " << p << ": " << current)

        SizeType last_state_id = 10000000;

        List<PointType> sequence;

        SizeType s = 0;
        for(; s < max_steps; ++s) {

            sequence.append(current);

            auto current_icell = point_to_cell(current,ra.state_grid(),ra.grid_depth(),ra.vertex_factory());

            auto const& assignment = assignments.get(current_icell);

            auto control_icell = assignment.control_cell();
            auto target_icell = assignment.target_cell();

            auto control_box = control_icell.cell().box();
            auto control_midpoint = control_box.midpoint();
            auto control = control_box;
            for (SizeType i=0; i<control_box.dimension(); ++i) {
                control[i] = ExactIntervalType(control_midpoint[i], cast_exact(control_midpoint[i].get_d() + 1e-10));
            }

            if (current_icell.id() != last_state_id) {
                CONCLOG_PRINTLN_AT(1,"Now in " << current_cell_box << " targeting " << assignments.get(current_icell).target_cell() << " with control " << control.midpoint())
                last_state_id = current_icell.id();
            }

            auto combined = product(current_icell.cell().box(), control_box);
            auto next_box = cast_exact_box(apply(ra.dynamics(), combined).bounding_box());

            auto next_midpoint = next_box.midpoint();
            PointType next(dim);
            for (SizeType i=0; i<dim; ++i) {
                next[i] = next_midpoint.at(i).get_d();
            }

            CONCLOG_PRINTLN_AT(1,"Next point: " << next)

            auto next_icell = point_to_cell(next,ra.state_grid(),ra.grid_depth(),ra.vertex_factory());
            CONCLOG_PRINTLN_AT(1, s << ": from " << current << " (" << current_icell.id() << ") to " << next << " (" << next_icell.id() << ")")

            if (not ra.state_paving().superset(next_icell.cell())) {
                CONCLOG_PRINTLN("The next cell is outside of the state paving, terminating with failure.")
                break;
            }

            if (ra.goals().superset(next_icell.cell())) {
                CONCLOG_PRINTLN("The next cell is a goal, terminating with success in " << s+1 << " steps.")
                break;
            }

            if (ra.obstacles().superset(next_icell.cell())) {
                CONCLOG_PRINTLN("The next cell is an obstacle, terminating with failure.")
                break;
            }

            if (ra.unverified().superset(next_icell.cell())) {
                CONCLOG_PRINTLN("The next cell is an unverified (hence unsafe) cell, terminating with failure.")
                auto const& ts = ra.possibly_reaching_graph().internal().forward_transitions(current_icell);
                CONCLOG_PRINTLN("Cell " << current_icell.id() << " of desired target " << target_icell.cell().box() << " ( " << assignments.get(current_icell).target_cell().id() << ") would reach:")
                for (auto const& t : ts) {
                    List<SizeType> target_ids;
                    for (auto const& tgt : t.second) {
                        target_ids.append(tgt.first.id());
                    }
                    CONCLOG_PRINTLN(t.first.id() << ": " << target_ids)
                }
                break;
            }

            current = next;
        }
        if (s == max_steps) {
            CONCLOG_PRINTLN("Not terminated after maximum steps, failure.")
        }

        create_workspace_1_matlab_file(p,sequence);
    }
}