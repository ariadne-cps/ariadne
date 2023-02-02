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
#include "utility/stopwatch.hpp"

void ariadne_main()
{
    Real deltat=0.1_dec;
    Real v=0.5_dec;
    RealVariable x("x"), y("y"), theta("theta"), Kx("Kx"), Ky("Ky"), Kt("Kt"), b("b");
    IteratedMap heading({next(x)=x+deltat*v*cos(theta),next(y)=y+deltat*v*sin(theta),next(theta)=theta+deltat*(Kx*x+Ky*y+Kt*theta+b),
                         next(Kx)=Kx,next(Ky)=Ky,next(Kt)=Kt,next(b)=b});
    CONCLOG_PRINTLN_VAR(heading);

    Stopwatch<Milliseconds> sw;
    CONCLOG_PRINTLN("Computing evolution...");

    auto function = heading.function();

    typedef BinaryWord SWord;
    typedef BinaryWord CWord;
    typedef GridTreePaving SPaving;
    typedef GridTreePaving CPaving;

    FloatDP eps(1e-10_x,DoublePrecision());

    Grid sgrid({0.5,0.5,2*pi/8});
    ExactBoxType sdomain({{0+eps,5-eps},{0+eps,5-eps},{0+eps,2*pi-eps}});
    SPaving sdomain_paving(sgrid);
    sdomain_paving.adjoin_outer_approximation(sdomain,0);
    sdomain_paving.mince(0);
    auto num_state_domain_cells = sdomain_paving.size();
    CONCLOG_PRINTLN_VAR_AT(1,num_state_domain_cells)

    SPaving obstacle_paving(sgrid);
    ExactBoxType obstacle1({{1+eps,3.5_x-eps},{4.5_x+eps,5-eps},{0+eps,2*pi-eps}});
    ExactBoxType obstacle2({{0+eps,1-eps},{2+eps,3-eps},{0+eps,2*pi-eps}});
    ExactBoxType obstacle3({{2.5_x+eps,5-eps},{2+eps,3-eps},{0+eps,2*pi-eps}});
    ExactBoxType obstacle4({{0+eps,5-eps},{0+eps,0.5_x-eps},{0+eps,2*pi-eps}});
    obstacle_paving.adjoin_outer_approximation(obstacle1,0);
    obstacle_paving.adjoin_outer_approximation(obstacle2,0);
    obstacle_paving.adjoin_outer_approximation(obstacle3,0);
    obstacle_paving.adjoin_outer_approximation(obstacle4,0);
    auto num_obstacle_cells = obstacle_paving.size();
    CONCLOG_PRINTLN_VAR_AT(1,num_obstacle_cells)

    SPaving goal_paving(sgrid);
    ExactBoxType goal({{4+eps,5-eps},{4.5_x+eps,5-eps},{0+eps,2*pi-eps}});
    goal_paving.adjoin_outer_approximation(goal,0);
    auto num_goal_cells = goal_paving.size();
    CONCLOG_PRINTLN_VAR_AT(1,num_goal_cells)

    SPaving scandidate_paving = sdomain_paving;
    scandidate_paving.remove(obstacle_paving);
    scandidate_paving.remove(goal_paving);
    auto num_state_candidate_cells = scandidate_paving.size();
    CONCLOG_PRINTLN_VAR_AT(1,num_state_candidate_cells)

    Grid cgrid({0,0,0,0.5_x},{1,1,1,1});
    ExactBoxType cdomain({{-2+eps,2-eps},{-2+eps,2-eps},{-1+eps,1-eps},{-2.5_x+eps,2.5_x-eps}});
    CPaving cdomain_paving(cgrid);
    cdomain_paving.adjoin_outer_approximation(cdomain,0);
    cdomain_paving.mince(0);
    auto num_controller_cells = cdomain_paving.size();
    CONCLOG_PRINTLN_VAR_AT(1,num_controller_cells)

    ExactBoxType graphics_box({{-1,6},{-1,6},{-1,8}});
    Figure fig(graphics_box,0,1);
    fig << fill_colour(white) << sdomain_paving << fill_colour(green) << goal_paving << fill_colour(blue) << obstacle_paving;
    fig.write("heading_control");

    SPaving targets_paving;
    SPaving obstacles_paving;
    Map<SWord,Map<CWord,SPaving>> forward_graph;
    Map<SWord,Map<CWord,SPaving>> backward_graph;

    for (auto const& state_cell : scandidate_paving) {
        auto const& state_word = state_cell.word();
        Map<CWord,SPaving> targets;
        for (auto const& controller_cell : cdomain_paving) {
            auto const& controller_word = controller_cell.word();
            auto combined = product(state_cell.box(),controller_cell.box());
            SPaving target_cells(sgrid);
            target_cells.adjoin_outer_approximation(project(apply(function, combined),Range(0,3)),0);
            targets.insert(make_pair(controller_word,target_cells));
        }
        forward_graph.insert(make_pair(state_word,targets));
    }
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of constructing post graph: " << sw.elapsed_seconds() << " seconds (per state: " << sw.elapsed_seconds()/num_state_candidate_cells/num_controller_cells << ")")
}