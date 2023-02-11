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

typedef GridCell NCell;
typedef GridCell ECell;
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
    typedef SizeTypeMap<NCell,SizeTypeMap<ECell,SPaving>>::iterator Iterator;
    typedef std::map<String,SizeType> HashTableType;
  public:
    DirectedHashedGraph(HashTableType const& hashed_vertices, HashTableType const& hashed_edges, SizeType default_vertex_extent, SizeType default_edge_extent) :
            _hashed_vertices(hashed_vertices), _hashed_edges(hashed_edges),
            _default_vertex_extent(default_vertex_extent), _default_edge_extent(default_edge_extent) { }

    SizeType num_sources() const { return _graph.size(); }

    SizeType num_transitions() const {
        SizeType result = 0;
        for (auto const& src : _graph) {
            for (auto const& ctrl : src.second.second) {
                result += ctrl.second.second.size();
            }
        }
        return result;
    }

    //! \brief Insert a forward entry from \a source_cell using \a transition_cell with associated \a target_cells
    void insert_forward(NCell const& source_cell, ECell const& transition_cell, SPaving const& target_cells) {
        auto ctrl_id = to_identifier(transition_cell, _default_edge_extent, _hashed_edges);
        auto src_id = to_identifier(source_cell, _default_vertex_extent, _hashed_vertices);
        auto src_ref = _graph.find(src_id);
        if (src_ref == _graph.end()) {
            _graph.insert(make_pair(src_id,make_pair(source_cell,Map<SizeType,Pair<ECell,SPaving>>())));
            src_ref = _graph.find(src_id);
        }
        src_ref->second.second.insert(make_pair(ctrl_id,make_pair(transition_cell, target_cells)));
    }

    //! \brief Insert backward entries from each of \a target_cells to \a source_cell using \a transition_cell, hence
    //! hashing on the target and transition
    void insert_backward(NCell const& source_cell, ECell const& transition_cell, SPaving const& target_cells) {
        auto ctrl_id = to_identifier(transition_cell, _default_edge_extent, _hashed_edges);
        for (auto const& src : target_cells) {
            auto src_id = to_identifier(src, _default_vertex_extent, _hashed_vertices);
            auto src_ref = _graph.find(src_id);
            if (src_ref == _graph.end()) {
                _graph.insert(make_pair(src_id,make_pair(src,Map<SizeType,Pair<ECell,SPaving>>())));
                src_ref = _graph.find(src_id);
            }
            auto tgt_ref = src_ref->second.second.find(ctrl_id);
            if (tgt_ref == src_ref->second.second.end()) {
                src_ref->second.second.insert(make_pair(ctrl_id,make_pair(transition_cell, SPaving(source_cell.grid()))));
                tgt_ref = src_ref->second.second.find(ctrl_id);
            }
            tgt_ref->second.second.adjoin(source_cell);
        }
    }

    Iterator find(NCell const& source_cell) {
        return _graph.find(to_identifier(source_cell, _default_vertex_extent, _hashed_vertices));
    }

    bool contains(Iterator const& iterator) const {
        return iterator != _graph.end();
    }

    void erase(NCell const& source_cell) {
        _graph.erase(to_identifier(source_cell, _default_vertex_extent, _hashed_vertices));
    }

    //! \brief Erase for a given \a source and \a transition its \a target
    void erase(NCell const& source, ECell const& transition, NCell const& target) {
        auto const& src_ref = _graph.find(to_identifier(source, _default_vertex_extent, _hashed_vertices));
        if (src_ref != _graph.end()) {
            auto const& trans_ref = src_ref->second.second.find(to_identifier(transition, _default_edge_extent, _hashed_edges));
            if (trans_ref != src_ref->second.second.end()) {
                trans_ref->second.second.remove(target);
            }
        }
    }

    //! \brief Subtract the source cells of the graph from \a paving
    void apply_to(SPaving& paving) const {
        for (auto const& entry : _graph) {
            paving.remove(entry.second.first);
        }
    }

    void clear() { _graph.clear(); }

    //! \brief Remove transitions for which the set of targets is empty
    //! and remove sources for which the set of transitions is empty
    void sweep() {
        auto src_it = _graph.begin();
        while (src_it != _graph.end()) {
            auto trans_it = src_it->second.second.begin();
            while (trans_it != src_it->second.second.end()) {
                if (trans_it->second.second.is_empty()) src_it->second.second.erase(trans_it++);
                else ++trans_it;
            }
            if (src_it->second.second.empty()) _graph.erase(src_it++);
            else  ++src_it;
        }
    }

    friend OutputStream& operator<<(OutputStream& os, DirectedHashedGraph const& g) {
        for (auto const& src : g._graph) {
            os << src.first << src.second.first.box() << ":[";
            for (auto const& ctrl : src.second.second) {
                os << ctrl.first << "->";
                os << "(";
                for (auto const& tgt : ctrl.second.second) {
                    os << to_identifier(tgt, g._default_vertex_extent, g._hashed_vertices) << ",";
                }
                os.seekp(-1, std::ios_base::end);
                os << ")";
            }
            os << "]\n";
        }
        os << "}";

        return os;
    }
private:
    std::map<String,SizeType> const _hashed_vertices;
    std::map<String,SizeType> const _hashed_edges;
    SizeType const _default_vertex_extent;
    SizeType const _default_edge_extent;
    SizeTypeMap<NCell,SizeTypeMap<ECell,SPaving>>  _graph;
};

class ReachabilityGraphInterface {
  public:
    virtual SizeType num_transitions() const = 0;
    virtual SizeType num_sources() const = 0;
    virtual SizeType num_targets() const = 0;

    virtual void insert(NCell const& source_cell, ECell const& transition_cell, SPaving const& target_cells) = 0;
    virtual void clear() = 0;
    virtual void refine_to_safety_graph(SPaving const& unsafe, SPaving& unverified) = 0;
    virtual ReachabilityGraphInterface* clone() const = 0;
    virtual void write(std::ostream& os) const = 0;
    virtual ~ReachabilityGraphInterface() = default;

    friend OutputStream& operator<<(OutputStream& os, ReachabilityGraphInterface const& g) {
        g.write(os);
        return os;
    }
};

class ForwardBackwardReachabilityGraph : public ReachabilityGraphInterface {
  public:
    typedef DirectedHashedGraph::HashTableType HashTableType;
  public:
    ForwardBackwardReachabilityGraph(HashTableType const& hashed_vertices, HashTableType const& hashed_edges, SizeType default_vertex_extent, SizeType default_edge_extent) :
        _forward_graph(hashed_vertices,hashed_edges,default_vertex_extent,default_edge_extent),
        _backward_graph(hashed_vertices,hashed_edges,default_vertex_extent,default_edge_extent) { }

    SizeType num_transitions() const override {
        return _forward_graph.num_transitions();
    }

    SizeType num_sources() const override {
        return _forward_graph.num_sources();
    }

    SizeType num_targets() const override {
        return _backward_graph.num_sources();
    }

    void insert(NCell const& source_cell, ECell const& transition_cell, SPaving const& target_cells) override {
        _forward_graph.insert_forward(source_cell,transition_cell,target_cells);
        _backward_graph.insert_backward(source_cell,transition_cell,target_cells);
    }

    void clear() override {
        _forward_graph.clear();
        _backward_graph.clear();
    }

    void refine_to_safety_graph(SPaving const& unsafe, SPaving& unverified) override {
        CONCLOG_SCOPE_CREATE

        std::deque<NCell> unsafe_cells_queue;
        for (auto const& c : unsafe) unsafe_cells_queue.push_back(c);

        double indicator_final_original = unsafe.size();
        double indicator_current_value = 0;
        ProgressIndicator indicator(indicator_final_original);
        while (not unsafe_cells_queue.empty()) {
            CONCLOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "% ");
            auto const& u = unsafe_cells_queue.front();
            auto const& bw_unsafe_ref = _backward_graph.find(u);
            if (_backward_graph.contains(bw_unsafe_ref)) {
                for (auto const& bw_unsafe_trans : bw_unsafe_ref->second.second) {
                    for (auto const& src_of_unsafe_trans : bw_unsafe_trans.second.second) {
                        auto const& fw_src_of_unsafe_trans = _forward_graph.find(src_of_unsafe_trans);
                        if (_forward_graph.contains(fw_src_of_unsafe_trans)) {
                            auto const& fw_trans = fw_src_of_unsafe_trans->second.second.at(bw_unsafe_trans.first);
                            for (auto const& tgt_to_prune : fw_trans.second) {
                                _backward_graph.erase(tgt_to_prune,bw_unsafe_trans.second.first,src_of_unsafe_trans);
                            }
                            fw_src_of_unsafe_trans->second.second.erase(bw_unsafe_trans.first);
                            if (fw_src_of_unsafe_trans->second.second.empty()) {
                                unsafe_cells_queue.push_back(src_of_unsafe_trans);
                                indicator.update_final(++indicator_final_original);
                                _forward_graph.erase(src_of_unsafe_trans);
                            }
                        }
                    }
                }
            }
            _backward_graph.erase(u);
            unsafe_cells_queue.pop_front();
            indicator.update_current(++indicator_current_value);
        }

        _backward_graph.sweep();

        _forward_graph.apply_to(unverified);
    }

    ReachabilityGraphInterface* clone() const override {
        return new ForwardBackwardReachabilityGraph(*this);
    }

    void write(std::ostream& os) const override {
        os << "Forward:\n" << _forward_graph << "\nBackward:\n" << _backward_graph;
    }

  private:
    DirectedHashedGraph _forward_graph;
    DirectedHashedGraph _backward_graph;
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

        _reachability_graph.reset(new ForwardBackwardReachabilityGraph(_state_ids, _controller_ids, _default_state_extent, _default_controller_extent));
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
    SizeType num_sources() const { return _reachability_graph->num_sources(); }
    SizeType num_targets() const { return _reachability_graph->num_targets(); }

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
        CONCLOG_PRINTLN("Goals: " << ss.str())
    }

    void print_obstacles() const {
        std::stringstream ss;
        for (auto const& o : _obstacles) ss << to_identifier(o,_default_state_extent,_state_ids) << " ";
        CONCLOG_PRINTLN("Obstacles: " << ss.str())
    }

    void print_graph() const {
        CONCLOG_PRINTLN(*_reachability_graph)
    }

    void compute_reachability_graph() {
        CONCLOG_SCOPE_CREATE

        _reachability_graph->clear();

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
                if (not target_cells.is_empty())
                    _reachability_graph->insert(state_cell,controller_cell,target_cells);
            }
            indicator.update_current(++indicator_value);
        }
    }

    void refine_to_safety_graph() {
        _reachability_graph->refine_to_safety_graph(_obstacles,_unverified);
    }

    SizeType num_transitions() const { return _reachability_graph->num_transitions(); }

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

    Ariadne::SharedPointer<ReachabilityGraphInterface> _reachability_graph;
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
    Grid control_grid({pi/4});
    BoundsBoxType control_domain({{-pi_,pi_}});

    ReachAvoidGridding scs("heading",dynamics,state_grid,state_domain,control_grid,control_domain,1e-10_x);

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
    CONCLOG_PRINTLN_AT(1,"Time cost of constructing forward/backward graph: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1,scs.num_transitions())

    sw.restart();
    scs.refine_to_safety_graph();
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Time cost of reducing forward graph to safe one: " << sw.elapsed_seconds() << " seconds")

    CONCLOG_PRINTLN_VAR_AT(1,scs.num_transitions())

    CONCLOG_PRINTLN_AT(1,"Safe abstract states: " << scs.num_sources())

    CONCLOG_PRINTLN_AT(1,"Unverified abstract states: " << scs.unverified_size() << " (" << scs.unverified_percentage() << "% left)")

    scs.plot({{0,5},{0,5},{-6.28_x,6.28_x}},0,1);
    scs.plot({{0,5},{0,5},{-6.28_x,6.28_x}},0,2);
    scs.plot({{0,5},{0,5},{-6.28_x,6.28_x}},1,2);

    CONCLOG_RUN_AT(3,scs.print_graph())
}