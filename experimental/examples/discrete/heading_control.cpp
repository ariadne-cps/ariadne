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

template<class K, class V> using SizeTypeMap = Map<SizeType,Pair<K,V>>;
typedef Interval<double> BoundType;
typedef Vector<BoundType> BoundsBoxType;

ExactBoxType shrink(ExactBoxType const& bx, FloatDP const& eps) {
    ExactBoxType result(bx.dimension());
    for (SizeType i=0; i<bx.size(); ++i)
        result[i] = Interval<FloatDP>(bx[i].lower_bound()+eps,bx[i].upper_bound()-eps);
    return result;
}

ExactBoxType shrink(BoundsBoxType const& bx, FloatDP const& eps) {
    ExactBoxType result(bx.size());
    for (SizeType i=0; i<bx.size(); ++i)
        result[i] = Interval<FloatDP>(FloatDP(cast_exact(bx[i].lower_bound()),DoublePrecision())+eps,FloatDP(cast_exact(bx[i].upper_bound()),DoublePrecision())-eps);
    return result;
}

String word_to_id(BinaryWord const& w, SizeType size) {
    std::stringstream ss;
    auto size_offset = w.size()-size;
    for (SizeType i=0; i<size; ++i) ss << to_string(w.at(size_offset+i));
    return ss.str();
}

SizeType to_identifier(GridCell const& cell, SizeType default_extent, std::map<String,SizeType> const& hashed) {
    SizeType size_to_use = cell.word().size() - (cell.root_extent() - default_extent)*cell.dimension();
    String expanded_id = word_to_id(cell.word(),size_to_use);
    return hashed.at(expanded_id);
}

class DirectedHashedGraph {
public:
    typedef SizeTypeMap<SCell,SizeTypeMap<CCell,SPaving>>::iterator Iterator;
public:
    DirectedHashedGraph(EffectiveVectorMultivariateFunction const& dynamics, std::map<String,SizeType> const& hashed_space, std::map<String,SizeType> const& hashed_controller, SizeType default_space_extent, SizeType default_controller_extent) :
            _dynamics(dynamics), _hashed_space(hashed_space), _hashed_controller(hashed_controller),
            _default_space_extent(default_space_extent), _default_controller_extent(default_controller_extent) { }

    SizeType sources_size() const { return _graph.size(); }

    SizeType transitions_size() const {
        SizeType result = 0;
        for (auto const& src : _graph) {
            for (auto const& ctrl : src.second.second) {
                result += ctrl.second.second.size();
            }
        }
        return result;
    }

    //! \brief Insert a forward entry from \a source_cell using \a controller_cell with associated \a target_cells
    void insert_forward(SCell const& source_cell, CCell const& controller_cell, SPaving const& target_cells) {
        auto ctrl_id = to_identifier(controller_cell,_default_controller_extent,_hashed_controller);
        auto src_id = to_identifier(source_cell,_default_space_extent,_hashed_space);
        auto src_ref = _graph.find(src_id);
        if (src_ref == _graph.end()) {
            _graph.insert(make_pair(src_id,make_pair(source_cell,Map<SizeType,Pair<CCell,SPaving>>())));
            src_ref = _graph.find(src_id);
        }
        src_ref->second.second.insert(make_pair(ctrl_id,make_pair(controller_cell,target_cells)));
    }

    //! \brief Insert backward entries from each of \a target_cells to \a source_cell using \a controller_cell, hence
    //! hashing on the target and controller
    void insert_backward(SCell const& source_cell, CCell const& controller_cell, SPaving const& target_cells) {
        auto ctrl_id = to_identifier(controller_cell,_default_controller_extent,_hashed_controller);
        for (auto const& src : target_cells) {
            auto src_id = to_identifier(src,_default_space_extent,_hashed_space);
            auto src_ref = _graph.find(src_id);
            if (src_ref == _graph.end()) {
                _graph.insert(make_pair(src_id,make_pair(src,Map<SizeType,Pair<CCell,SPaving>>())));
                src_ref = _graph.find(src_id);
            }
            auto tgt_ref = src_ref->second.second.find(ctrl_id);
            if (tgt_ref == src_ref->second.second.end()) {
                src_ref->second.second.insert(make_pair(ctrl_id,make_pair(controller_cell,SPaving(source_cell.grid()))));
                tgt_ref = src_ref->second.second.find(ctrl_id);
            }
            tgt_ref->second.second.adjoin(source_cell);
        }
    }

    Iterator find(SCell const& source_cell) {
        return _graph.find(to_identifier(source_cell,_default_space_extent,_hashed_space));
    }

    bool contains(Iterator const& iterator) const {
        return iterator != _graph.end();
    }

    void erase(SCell const& source_cell) {
        _graph.erase(to_identifier(source_cell,_default_space_extent,_hashed_space));
    }

    //! \brief Subtract the source cells of the graph from \a paving
    void prune(SPaving& paving) const {
        for (auto const& entry : _graph) {
            paving.remove(entry.second.first);
        }
    }

    void clear() { _graph.clear(); }

    friend OutputStream& operator<<(OutputStream& os, DirectedHashedGraph const& g) {
        for (auto const& src : g._graph) {
            os << src.first << src.second.first.box() << ":[\n";
            for (auto const& ctrl : src.second.second) {
                os << ctrl.first << "->";
                auto combined = product(src.second.first.box(),ctrl.second.first.box());
                os << apply(g._dynamics,combined);
                os << "(";
                for (auto const& tgt : ctrl.second.second) {
                    os << to_identifier(tgt,g._default_space_extent,g._hashed_space) << ",";
                }
                os.seekp(-1, std::ios_base::end);
                os << ")\n";
            }
            os.seekp(-1, std::ios_base::end);
            os << "]\n";
        }
        os << "}";

        return os;
    }
private:
    EffectiveVectorMultivariateFunction const _dynamics;
    std::map<String,SizeType> const _hashed_space;
    std::map<String,SizeType> const _hashed_controller;
    SizeType const _default_space_extent;
    SizeType const _default_controller_extent;
    SizeTypeMap<SCell,SizeTypeMap<CCell,SPaving>>  _graph;
};

class ReachAvoidGridding {
public:
    ReachAvoidGridding(String const& name, EffectiveVectorMultivariateFunction const& dynamics, Grid const& state_grid, BoundsBoxType const& state_bounds, Grid const& control_grid, BoundsBoxType const& control_bounds, ExactDouble eps) :
            _name(name), _dynamics(dynamics), _eps({eps,DoublePrecision()}) {

        _state_paving = SPaving(state_grid);
        _state_paving.adjoin_outer_approximation(shrink(state_bounds,_eps),0);
        _state_paving.mince(0);

        _controller_paving = CPaving(control_grid);
        _controller_paving.adjoin_outer_approximation(shrink(control_bounds,_eps),0);
        _controller_paving.mince(0);

        _default_state_extent = _state_paving.begin()->root_extent();
        _default_controller_extent = _controller_paving.begin()->root_extent();

        _unverified = _state_paving;
        _obstacles = SPaving(_state_paving.grid());
        _goals = SPaving(_state_paving.grid());

        for (auto const& c : _state_paving)
            _state_ids.insert(make_pair(word_to_id(c.word(),_default_state_extent*_state_paving.dimension()),_state_ids.size()));
        for (auto const& c : _controller_paving) {
            _controller_ids.insert(make_pair(word_to_id(c.word(),_default_controller_extent*_controller_paving.dimension()),_controller_ids.size()));
        }

        _forward_graph.reset(new DirectedHashedGraph(_dynamics,_state_ids,_controller_ids,_default_state_extent,_default_controller_extent));
        _backward_graph.reset(new DirectedHashedGraph(_dynamics,_state_ids,_controller_ids,_default_state_extent,_default_controller_extent));
    }

    ReachAvoidGridding& add_obstacle(BoundsBoxType const& box) {
        SPaving obstacle_paving(_state_paving.grid());
        obstacle_paving.adjoin_outer_approximation(shrink(box,_eps),0);
        obstacle_paving.restrict(_state_paving);
        _unverified.remove(obstacle_paving);
        _obstacles.adjoin(obstacle_paving);
        return *this;
    }

    ReachAvoidGridding& add_goal(BoundsBoxType const& box) {
        SPaving goal_paving(_state_paving.grid());
        goal_paving.adjoin_outer_approximation(shrink(box,_eps),0);
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
    SizeType forward_sources_size() const { return _forward_graph->sources_size(); }

    //! \brief The percentage (in the 0-100 scale) of still unverified states
    double unverified_percentage() const {
        return static_cast<double>(_unverified.size())*100/(_state_paving.size()-_goals.size()-_obstacles.size());
    }

    void plot(ExactBoxType const& graphics_box, SizeType xaxis, SizeType yaxis) {
        Figure fig(graphics_box,xaxis,yaxis);

        SPaving safe = _state_paving;
        safe.remove(_obstacles);
        safe.remove(_goals);
        safe.remove(_unverified);
        fig << fill_colour(white) << _state_paving << fill_colour(green) << _goals << fill_colour(blue) << _obstacles << fill_colour(yellow) << safe;
        char num_char[64] = "";
        snprintf(num_char,64,"[%zu,%zu]",xaxis,yaxis);
        fig.write((_name+num_char).c_str());
    }

    void print_goals() const {
        std::stringstream ss;
        for (auto const& g : _goals) ss << to_identifier(g,_default_state_extent,_state_ids) << " ";
        CONCLOG_PRINTLN_AT(2,"Goals: " << ss.str())
    }

    void print_obstacles() const {
        std::stringstream ss;
        for (auto const& o : _obstacles) ss << to_identifier(o,_default_state_extent,_state_ids) << " ";
        CONCLOG_PRINTLN_AT(2,"Obstacles: " << ss.str())
    }

    void print_forward_graph() const {
        CONCLOG_PRINTLN_AT(1,"Forward graph:")
        CONCLOG_PRINTLN_AT(1,*_forward_graph)
    }

    void print_backward_graph() const {
        CONCLOG_PRINTLN_AT(1,"Backward graph:")
        CONCLOG_PRINTLN_AT(1,*_backward_graph)
    }

    void reconstruct_reachability_graph() {
        CONCLOG_SCOPE_CREATE

        _forward_graph->clear();
        _backward_graph->clear();

        Stopwatch<Milliseconds> sw;

        ProgressIndicator indicator(_unverified.size());
        double indicator_value = 0;
        for (auto const& state_cell : _unverified) {
            CONCLOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "% ");
            for (auto const& controller_cell : _controller_paving) {
                auto combined = product(state_cell.box(),controller_cell.box());
                SPaving target_cells(_state_paving.grid());
                target_cells.adjoin_outer_approximation(shrink(cast_exact_box(apply(_dynamics, combined).bounding_box()),_eps),0);
                target_cells.mince(0);
                target_cells.restrict(_state_paving);
                if (not target_cells.is_empty()) {
                    _forward_graph->insert_forward(state_cell,controller_cell,target_cells);
                    _backward_graph->insert_backward(state_cell,controller_cell,target_cells);
                }
            }
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
        while (not unsafe.empty()) {
            CONCLOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "% ");
            auto const& u = unsafe.front();
            auto const& bref = _backward_graph->find(u);
            if (_backward_graph->contains(bref)) {
                for (auto const& trans : bref->second.second) {
                    for (auto const& src : trans.second.second) {
                        auto const& fref = _forward_graph->find(src);
                        if (_forward_graph->contains(fref)) {
                            fref->second.second.erase(trans.first);
                            if (fref->second.second.empty()) {
                                unsafe.push_back(src);
                                indicator.update_final(++indicator_final_original);
                                _forward_graph->erase(src);
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

        _forward_graph->prune(_unverified);
    }

    SizeType forward_transitions() const { return _forward_graph->transitions_size(); }
    SizeType backward_transitions() const { return _backward_graph->transitions_size(); }

private:

    String const _name;
    EffectiveVectorMultivariateFunction const _dynamics;

    SPaving _state_paving;
    CPaving _controller_paving;

    SizeType _default_state_extent;
    SizeType _default_controller_extent;

    FloatDP const _eps;

    SPaving _unverified;

    SPaving _obstacles;
    SPaving _goals;

    std::map<String,SizeType> _state_ids, _controller_ids;

    Ariadne::SharedPointer<DirectedHashedGraph> _forward_graph;
    Ariadne::SharedPointer<DirectedHashedGraph> _backward_graph;
};

void ariadne_main()
{
    Real deltat=0.1_dec;
    Real v=3;
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
    Grid control_grid({pi/32});
    BoundsBoxType control_domain({{-pi_,pi_}});

    ReachAvoidGridding scs("heading",dynamics,state_grid,state_domain,control_grid,control_domain,1e-10_x);

    CONCLOG_PRINTLN_VAR_AT(1,scs.state_size())
    CONCLOG_PRINTLN_VAR_AT(1,scs.controller_size())

    scs.add_obstacle({{1,3.5},{4.5,5},theta_domain});
    scs.add_obstacle({{0,1},{2,3},theta_domain});
    scs.add_obstacle({{2.5,5},{2,3},theta_domain});
    scs.add_obstacle({{0,5},{0,0.5},theta_domain});
    scs.print_obstacles();
    CONCLOG_PRINTLN_VAR_AT(1,scs.obstacles_size())

    scs.add_goal({{4,5},{4.5,5},theta_domain});

    scs.print_goals();
    CONCLOG_PRINTLN_VAR_AT(1,scs.goals_size())

    CONCLOG_PRINTLN_VAR_AT(1,scs.unverified_size())

    scs.reconstruct_reachability_graph();

    CONCLOG_PRINTLN_VAR_AT(1,scs.forward_transitions())
    CONCLOG_PRINTLN_VAR_AT(1,scs.backward_transitions())

    //scs.print_forward_graph();
    //scs.print_backward_graph();

    scs.compute_safe_forward_graph();

    CONCLOG_PRINTLN_AT(1,"Safe abstract states: " << scs.forward_sources_size())

    CONCLOG_PRINTLN_AT(1,"Unverified abstract states: " << scs.unverified_size() << " (" << scs.unverified_percentage() << "% left)")

    scs.plot({{0,5},{0,5},{-6.28_x,6.28_x}},0,1);
    scs.plot({{0,5},{0,5},{-6.28_x,6.28_x}},0,2);
    scs.plot({{0,5},{0,5},{-6.28_x,6.28_x}},1,2);

    CONCLOG_RUN_AT(2,scs.print_forward_graph())
}