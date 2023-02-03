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

typedef GridCell SCell;
typedef GridCell CCell;
typedef GridTreePaving SPaving;
typedef GridTreePaving CPaving;

typedef Map<SCell,Map<CCell,SPaving>> DirectedGraph;

std::string to_identifier(BinaryWord const& w, std::map<int,int>& hashed) {
    int expanded_id = 0;
    for (SizeType i=0; i<w.size(); ++i) expanded_id += (2<<i)*w.at(i);

    auto ref = hashed.find(expanded_id);
    int result = hashed.size();
    if (ref != hashed.end()) result = ref->second;
    else hashed.insert(make_pair(expanded_id,hashed.size()));

    return std::to_string(result);
}

void print_directed_graph(DirectedGraph const& g) {

    std::stringstream sstr;
    sstr << "{";

    std::map<int,int> hashed_space, hashed_controller;

    for (auto const& src : g) {
        sstr << to_identifier(src.first.word(),hashed_space) << ":[";
        for (auto const& ctrl : src.second) {
            sstr << to_identifier(ctrl.first.word(),hashed_controller) << "->(";
            for (auto const& tgt : ctrl.second) {
                sstr << to_identifier(tgt.word(),hashed_space) << ",";
            }
            sstr << "),";
        }
        sstr << "],";
    }
    sstr << "}";

    CONCLOG_PRINTLN_AT(1,sstr.str())
}

void ariadne_main()
{
    Real deltat=0.1_dec;
    Real v=3.3_dec;
    RealVariable x("x"), y("y"), theta("theta"), Kx("Kx"), Ky("Ky"), Kt("Kt"), b("b");
    IteratedMap heading({next(x)=x+deltat*v*cos(theta),next(y)=y+deltat*v*sin(theta),next(theta)=theta+deltat*(Kx*x+Ky*y+Kt*theta+b),
                         next(Kx)=Kx,next(Ky)=Ky,next(Kt)=Kt,next(b)=b});

    auto dynamics = heading.function().zeros(3,7);
    for (SizeType i=0; i<3; ++i)
        dynamics[i] = heading.function().get(i);
    CONCLOG_PRINTLN_VAR(dynamics);

    FloatDP eps(1e-10_x,DoublePrecision());

    Grid sgrid({0.5,0.5,2*pi/8});
    //Grid sgrid({1,1,2*pi});
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

    Grid cgrid({0,0,0,0.25_x},{0.5_x,0.5_x,0.5_x,0.5_x});
    //Grid cgrid({0,0,0,0.25_x},{2,2,1,2});
    ExactBoxType cdomain({{-1+eps,1-eps},{-1+eps,1-eps},{-0.5_x+eps,0.5_x-eps},{-1.25_x+eps,1.25_x-eps}});
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
    DirectedGraph forward_graph;
    DirectedGraph backward_graph;

    Stopwatch<Milliseconds> sw;

    for (auto const& state_cell : scandidate_paving) {
        Map<CCell,SPaving> targets;
        for (auto const& controller_cell : cdomain_paving) {
            auto combined = product(state_cell.box(),controller_cell.box());
            SPaving target_cells(sgrid);
            target_cells.adjoin_outer_approximation(apply(dynamics, combined),0);
            target_cells.mince(0);
            target_cells.restrict(sdomain_paving);
            targets.insert(make_pair(controller_cell,target_cells));
            for (auto const& tc : target_cells) {
                auto tref = backward_graph.find(tc);
                if (tref == backward_graph.end()) {
                    backward_graph.insert(make_pair(tc,Map<CCell,SPaving>()));
                    tref = backward_graph.find(tc);
                }
                auto sref = tref->second.find(controller_cell);
                if (sref == tref->second.end()) {
                    tref->second.insert(make_pair(controller_cell,SPaving(sgrid)));
                    sref = tref->second.find(controller_cell);
                }
                sref->second.adjoin(state_cell);
            }
        }
        forward_graph.insert(make_pair(state_cell,targets));
    }
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of constructing forward/backward graph: " << sw.elapsed_seconds() << " seconds (per state: " << sw.elapsed_seconds()/num_state_candidate_cells/num_controller_cells << ")")

    SizeType forward_transitions_size = 0;
    for (auto const& src : forward_graph) {
        for (auto const& ctrl : src.second) {
            forward_transitions_size += ctrl.second.size();
        }
    }
    CONCLOG_PRINTLN_AT(1,"Graph transitions: " << forward_transitions_size)

    SizeType backward_transitions_size = 0;
    for (auto const& tgt : backward_graph) {
        for (auto const& ctrl : tgt.second) {
            backward_transitions_size += ctrl.second.size();
        }
    }
    CONCLOG_PRINTLN_AT(1,"Backward starting cells: " << backward_graph.size())

    ARIADNE_ASSERT_EQUAL(forward_transitions_size,backward_transitions_size)

    print_directed_graph(forward_graph);
}