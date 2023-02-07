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

SizeType word_to_id(BinaryWord const& w) {
    SizeType result = 0;
    for (SizeType i=0; i<w.size(); ++i) result += static_cast<SizeType>((2<<i)*w.at(i));
    return result;
}

std::string to_identifier(GridCell const& cell, std::map<SizeType,SizeType> const& hashed) {
    SizeType expanded_id = word_to_id(cell.word());
    return to_string(hashed.at(expanded_id));
}

void print_directed_graph(DirectedGraph const& g, std::map<SizeType,SizeType> const& hashed_space, std::map<SizeType,SizeType> const& hashed_controller, bool show_boxes = false) {

    std::stringstream sstr;
    sstr << "{\n";

    for (auto const& src : g) {
        sstr << to_identifier(src.first,hashed_space);
        if (show_boxes) sstr << src.first.box();
        sstr << ":[";
        for (auto const& ctrl : src.second) {
            sstr << to_identifier(ctrl.first,hashed_controller) << "->(";
            for (auto const& tgt : ctrl.second) {
                sstr << to_identifier(tgt,hashed_space) << ",";
            }
            sstr.seekp(-1, std::ios_base::end);
            sstr << "),";
        }
        sstr.seekp(-1, std::ios_base::end);
        sstr << "]\n";
    }
    sstr << "}";

    CONCLOG_PRINTLN_AT(1,sstr.str())
}

class ReachAvoidGridding {
  public:
    ReachAvoidGridding(EffectiveVectorMultivariateFunction const& dynamics, SPaving const& state_paving, CPaving const& controller_paving) :
        _dynamics(dynamics), _state_paving(state_paving), _controller_paving(controller_paving), _unverified(state_paving),
        _obstacles(state_paving.grid()), _goals(state_paving.grid()) {
        for (auto const& c : state_paving)
            _state_ids.insert(make_pair(word_to_id(c.word()),_state_ids.size()));
        for (auto const& c : controller_paving)
            _controller_ids.insert(make_pair(word_to_id(c.word()),_controller_ids.size()));
    }

    ReachAvoidGridding& add_obstacle(ExactBoxType const& box) {
        SPaving obstacle_paving(_state_paving.grid());
        obstacle_paving.adjoin_outer_approximation(box,0);
        obstacle_paving.restrict(_state_paving);
        _unverified.remove(obstacle_paving);
        _obstacles.adjoin(obstacle_paving);
        return *this;
    }

    ReachAvoidGridding& add_goal(ExactBoxType const& box) {
        SPaving goal_paving(_state_paving.grid());
        goal_paving.adjoin_outer_approximation(box,0);
        goal_paving.restrict(_state_paving);
        _unverified.remove(goal_paving);
        _goals.adjoin(goal_paving);
        return *this;
    }

    SizeType state_size() const { return _state_paving.size(); }
    SizeType controller_size() const { return _controller_paving.size(); }
    SizeType obstacles_size() const { return _obstacles.size(); }
    SizeType goals_size() const { return _goals.size(); }
    SizeType unverified_size() const { return _unverified.size(); }
    SizeType forward_sources_size() const { return _forward_graph.size(); }

    void plot(ExactBoxType const& graphics_box, SizeType xaxis, SizeType yaxis) {
        Figure fig(graphics_box,xaxis,yaxis);

        SPaving safe = _state_paving;
        safe.remove(_obstacles);
        safe.remove(_goals);
        safe.remove(_unverified);
        CONCLOG_PRINTLN_AT(2,"Safe: " << safe.size())
        fig << fill_colour(white) << _state_paving << fill_colour(green) << _goals << fill_colour(blue) << _obstacles << fill_colour(yellow) << safe;
        fig.write("state_grid");
    }

    void print_goals() const {
        std::stringstream ss;
        for (auto const& g : _goals) ss << to_identifier(g,_state_ids) << " ";
        CONCLOG_PRINTLN_AT(2,"Goals: " << ss.str())
    }

    void print_obstacles() const {
        std::stringstream ss;
        for (auto const& o : _obstacles) ss << to_identifier(o,_state_ids) << " ";
        CONCLOG_PRINTLN_AT(2,"Obstacles: " << ss.str())
    }

    void print_forward_graph(bool show_boxes = false) const {
        print_directed_graph(_forward_graph,_state_ids,_controller_ids, show_boxes);
    }

    void print_backward_graph(bool show_boxes = false) const {
        print_directed_graph(_backward_graph,_state_ids,_controller_ids, show_boxes);
    }

    void reconstruct_reachability_graph() {
        CONCLOG_SCOPE_CREATE

        _forward_graph.clear();
        _backward_graph.clear();

        Stopwatch<Milliseconds> sw;

        ProgressIndicator indicator(_unverified.size());
        double indicator_value = 0;
        for (auto const& state_cell : _unverified) {
            CONCLOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "% ");
            Map<CCell,SPaving> targets;
            for (auto const& controller_cell : _controller_paving) {
                auto combined = product(state_cell.box(),controller_cell.box());
                SPaving target_cells(_state_paving.grid());
                target_cells.adjoin_outer_approximation(apply(_dynamics, combined),0);
                target_cells.mince(0);
                target_cells.restrict(_state_paving);
                if (not target_cells.is_empty())
                    targets.insert(make_pair(controller_cell,target_cells));
                for (auto const& tc : target_cells) {
                    auto tref = _backward_graph.find(tc);
                    if (tref == _backward_graph.end()) {
                        _backward_graph.insert(make_pair(tc,Map<CCell,SPaving>()));
                        tref = _backward_graph.find(tc);
                    }
                    auto sref = tref->second.find(controller_cell);
                    if (sref == tref->second.end()) {
                        tref->second.insert(make_pair(controller_cell,SPaving(_state_paving.grid())));
                        sref = tref->second.find(controller_cell);
                    }
                    sref->second.adjoin(state_cell);
                }
            }
            _forward_graph.insert(make_pair(state_cell,targets));
            indicator.update_current(++indicator_value);
        }
        sw.click();
        CONCLOG_PRINTLN_AT(1,"Time cost of constructing forward/backward graph: " << sw.elapsed_seconds() << " seconds (per state: " << sw.elapsed_seconds()/_unverified.size()/_controller_paving.size() << ")")
    }

    void compute_safe_forward_graph() {
        CONCLOG_SCOPE_CREATE

        Stopwatch<Milliseconds> sw;

        std::deque<SCell> unsafe;
        for (auto const& c : _obstacles) unsafe.push_back(c);

        double indicator_final_original = unsafe.size();
        double indicator_current_value = 0;
        ProgressIndicator indicator(indicator_final_original);
        while(not unsafe.empty()) {
            CONCLOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "% ");
            auto const& u = unsafe.front();
            auto const& bref = _backward_graph.find(u);
            if (bref != _backward_graph.end()) {
                for (auto const& trans : bref->second) {
                    for (auto const& src : trans.second) {
                        auto const& fref = _forward_graph.find(src);
                        if (fref != _forward_graph.end()) {
                            fref->second.erase(trans.first);
                            if (fref->second.empty()) {
                                unsafe.push_back(src);
                                indicator.update_final(++indicator_final_original);
                                _forward_graph.erase(src);
                            }
                        }
                    }
                }
            }
            unsafe.pop_front();
            indicator.update_current(++indicator_current_value);
        }
        sw.click();
        CONCLOG_PRINTLN_AT(1,"Time cost of reducing forward graph to safe one: " << sw.elapsed_seconds() << " seconds")

        for (auto const& entry : _forward_graph)
            _unverified.remove(entry.first);
    }

    SizeType forward_transitions() const {
        SizeType result = 0;
        for (auto const& src : _forward_graph) {
            for (auto const& ctrl : src.second) {
                result += ctrl.second.size();
            }
        }
        return result;
    }

    SizeType backward_transitions() const {
        SizeType result = 0;
        for (auto const& tgt : _backward_graph) {
            for (auto const& ctrl : tgt.second) {
                result += ctrl.second.size();
            }
        }
        return result;
    }

  private:

    EffectiveVectorMultivariateFunction const _dynamics;

    SPaving const _state_paving;
    CPaving const _controller_paving;

    SPaving _unverified;

    SPaving _obstacles;
    SPaving _goals;

    std::map<SizeType,SizeType> _state_ids, _controller_ids;

    DirectedGraph _forward_graph;
    DirectedGraph _backward_graph;
};

void ariadne_main()
{
    Real deltat=0.1_dec;
    Real v=3.3_dec;
    RealVariable x("x"), y("y"), theta("theta"), Kx("Kx"), Ky("Ky"), Kt("Kt"), b("b");
    IteratedMap heading({next(x)=x+deltat*v*cos(theta),next(y)=y+deltat*v*sin(theta),next(theta)= theta+deltat*(Kx*x+Ky*y+Kt*theta+b),
                         next(Kx)=Kx,next(Ky)=Ky,next(Kt)=Kt,next(b)=b});

    auto dynamics = heading.function().zeros(3,7);
    for (SizeType i=0; i<3; ++i)
        dynamics[i] = heading.function().get(i);
    CONCLOG_PRINTLN_VAR(dynamics);

    FloatDP eps(1e-10_x,DoublePrecision());

    Grid sgrid({0,0,eps.get_d()/2},{0.5,0.5,2*pi/8});
    ExactBoxType sdomain({{0+eps,5-eps},{0+eps,5-eps},{0+eps,2*pi-eps}});
    SPaving sdomain_paving(sgrid);
    sdomain_paving.adjoin_outer_approximation(sdomain,0);
    sdomain_paving.mince(0);

    Grid cgrid({0,0,0,0},{1.0_x,1.0_x,1.0_x,1.0_x});
    ExactBoxType cdomain({{-1+eps,1-eps},{-1+eps,1-eps},{-1+eps,1-eps},{-10+eps,10-eps}});
    CPaving cdomain_paving(cgrid);
    cdomain_paving.adjoin_outer_approximation(cdomain,0);
    cdomain_paving.mince(0);

    ReachAvoidGridding scs(dynamics,sdomain_paving, cdomain_paving);
    CONCLOG_PRINTLN_VAR_AT(1,scs.state_size())
    CONCLOG_PRINTLN_VAR_AT(1,scs.controller_size())

    scs.add_obstacle({{1+eps,3.5_x-eps},{4.5_x+eps,5-eps},{0+eps,2*pi-eps}});
    scs.add_obstacle({{0+eps,1-eps},{2+eps,3-eps},{0+eps,2*pi-eps}});
    scs.add_obstacle({{2.5_x+eps,5-eps},{2.0_x+eps,3-eps},{0+eps,2*pi-eps}});
    scs.add_obstacle({{0+eps,5-eps},{0+eps,0.5_x-eps},{0+eps,2*pi-eps}});
    scs.print_obstacles();
    CONCLOG_PRINTLN_VAR_AT(1,scs.obstacles_size())

    scs.add_goal({{4+eps,5-eps},{4.5_x+eps,5-eps},{0+eps,2*pi-eps}});
    scs.print_goals();
    CONCLOG_PRINTLN_VAR_AT(1,scs.goals_size())

    CONCLOG_PRINTLN_VAR_AT(1,scs.unverified_size())

    scs.reconstruct_reachability_graph();

    CONCLOG_PRINTLN_VAR_AT(1,scs.forward_transitions())
    CONCLOG_PRINTLN_VAR_AT(1,scs.backward_transitions())
/*
    CONCLOG_PRINTLN_AT(1,"Initial forward graph:")
    scs.print_forward_graph();
    CONCLOG_PRINTLN_AT(1,"Backward graph:")
    scs.print_backward_graph();
*/
    scs.compute_safe_forward_graph();

    CONCLOG_PRINTLN_AT(1,"Safe abstract states #: " << scs.forward_sources_size())

    CONCLOG_PRINTLN_VAR_AT(1,scs.unverified_size())

    scs.plot({{-1,6},{-1,6},{-1,8}},0,1);

    //CONCLOG_PRINTLN_AT(2,"Safe forward graph:")
    //CONCLOG_RUN_AT(1,scs.print_forward_graph(true))
}