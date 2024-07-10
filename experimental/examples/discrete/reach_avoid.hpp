/***************************************************************************
 *            reach_avoid.hpp
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

#include "reachability_graph.hpp"

template<class K, class V> using SizeTypeMap = Map<SizeType,Pair<K,V>>;
typedef Interval<double> BoundType;
typedef Vector<BoundType> BoundsBoxType;

ExactBoxType shrink(ExactBoxType const& bx, FloatDP const& eps) {
    ExactBoxType result(bx.dimension());
    for (SizeType i=0; i<bx.size(); ++i) {
        result[i] = Interval<FloatDP>((eps+bx[i].lower_bound()).lower_raw(),(-eps+bx[i].upper_bound()).upper_raw());
    }
    return result;
}

ExactBoxType shrink(BoundsBoxType const& bx, FloatDP const& eps) {
    ExactBoxType result(bx.size());
    for (SizeType i=0; i<bx.size(); ++i)
        result[i] = Interval<FloatDP>((FloatDP(cast_exact(bx[i].lower_bound()),DoublePrecision())+eps).lower_raw(),(FloatDP(cast_exact(bx[i].upper_bound()),DoublePrecision())-eps).upper_raw());
    return result;
}

class ReachAvoid {
public:
    ReachAvoid(String const& name, EffectiveVectorMultivariateFunction const& dynamics, Grid const& state_grid, BoundsBoxType const& state_bounds, Grid const& control_grid, BoundsBoxType const& control_bounds, SizeType depth, ExactDouble eps, ProbabilityType probability_threshold) :
            _name(name), _dynamics(dynamics), _state_bounds(state_bounds), _control_bounds(control_bounds), _depth(depth), _eps({eps,DoublePrecision()}), _probability_threshold(probability_threshold) {

        _state_paving = SPaving(state_grid);
        _state_paving.adjoin_outer_approximation(shrink(state_bounds,_eps),depth);
        _state_paving.mince(depth);

        _control_paving = CPaving(control_grid);
        _control_paving.adjoin_outer_approximation(shrink(control_bounds,_eps),depth);
        _control_paving.mince(depth);

        SizeType default_vertex_extent = _state_paving.begin()->root_extent();
        SizeType default_edge_extent = _control_paving.begin()->root_extent();
        IdentifiedCellFactory::HashTableType vertex_ids;
        IdentifiedCellFactory::HashTableType edge_ids;

        _unverified = _state_paving;
        _obstacles = SPaving(_state_paving.grid());
        _goals = SPaving(_state_paving.grid());

        for (auto const& c : _state_paving)
            vertex_ids.insert(make_pair(word_to_id(c.word(),(default_vertex_extent+depth)*_state_paving.dimension()),vertex_ids.size()));
        for (auto const& c : _control_paving) {
            edge_ids.insert(make_pair(word_to_id(c.word(),(default_edge_extent+depth)*_control_paving.dimension()),edge_ids.size()));
        }

        _reachability_graph.reset(new ForwardBackwardReachabilityGraph(IdentifiedCellFactory(default_vertex_extent,vertex_ids),IdentifiedCellFactory(default_edge_extent,edge_ids)));
    }

    ReachAvoid& add_obstacle(BoundsBoxType const& box) {
        SPaving obstacle_paving(_state_paving.grid());
        obstacle_paving.adjoin_outer_approximation(shrink(box,_eps),_depth);
        obstacle_paving.restrict(_state_paving);
        _unverified.remove(obstacle_paving);
        _obstacles.adjoin(obstacle_paving);
        return *this;
    }

    ReachAvoid& add_goal(BoundsBoxType const& box) {
        SPaving goal_paving(_state_paving.grid());
        goal_paving.adjoin_outer_approximation(shrink(box,_eps),_depth);
        goal_paving.restrict(_state_paving);
        _unverified.remove(goal_paving);
        _goals.adjoin(goal_paving);
        return *this;
    }

    BoundsBoxType const& state_bounds() const { return _state_bounds; }
    BoundsBoxType const& control_bounds() const { return _control_bounds; }

    SizeType state_size() const { return _state_paving.size(); }
    SizeType control_size() const { return _control_paving.size(); }
    SizeType obstacles_size() const { return _obstacles.size(); }
    SizeType goals_size() const { return _goals.size(); }
    SizeType unverified_size() const { return _unverified.size(); }
    SizeType num_sources() const { return _reachability_graph->num_sources(); }
    SizeType num_destinations() const { return _reachability_graph->num_destinations(); }

    //! \brief The percentage (in the 0-100 scale) of still unverified states
    double unverified_percentage() const {
        return static_cast<double>(_unverified.size())*100/(_state_paving.size()-_goals.size()-_obstacles.size());
    }

    void plot(SizeType xaxis, SizeType yaxis) const {
        ExactBoxType graphics_box(_state_bounds.size());
        for (SizeType i=0; i<_state_bounds.size(); ++i)
            graphics_box[i] = FloatDPExactInterval({ExactDouble(_state_bounds[0].lower_bound()),DoublePrecision()},{ExactDouble(_state_bounds[0].upper_bound()),DoublePrecision()});

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

    void plot() const {
        for (SizeType i=0; i<_dynamics.result_size()-1; ++i)
            for (SizeType j=i+1; j<_dynamics.result_size(); ++j)
                plot(i,j);
    }

    void print_goals() const {
        std::stringstream ss;
        for (auto const& g : _goals) ss << _reachability_graph->vertex_id(g) << " ";
        CONCLOG_PRINTLN("Goals: " << ss.str())
    }

    void print_obstacles() const {
        std::stringstream ss;
        for (auto const& o : _obstacles) ss << _reachability_graph->vertex_id(o) << " ";
        CONCLOG_PRINTLN("Obstacles: " << ss.str())
    }

    void print_graph() const {
        CONCLOG_PRINTLN(*_reachability_graph)
    }

    void compute_reachability_graph() {
        CONCLOG_SCOPE_CREATE

        _reachability_graph->clear();

        ProgressIndicator indicator(_unverified.size());
        for (auto const& source_cell : _unverified) {
            CONCLOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "% ");
            for (auto const& controller_cell : _control_paving) {
                auto combined = product(source_cell.box(), controller_cell.box());
                SPaving destination_paving(_state_paving.grid());
                auto destination_box = shrink(cast_exact_box(apply(_dynamics, combined).bounding_box()), _eps);
                destination_paving.adjoin_outer_approximation(destination_box, _depth);
                destination_paving.mince(_depth);
                destination_paving.restrict(_state_paving);
                if (not destination_paving.is_empty()) {
                    List<Pair<NCell,ProbabilityType>> destination_cells;

                    double total_volume = 0;
                    for (auto const& cell : destination_paving) {
                        auto intersection_box = intersection(cell.box(),destination_box);

                        auto current_volume = intersection_box.volume().get_d();
                        destination_cells.append({cell,current_volume});
                        total_volume += current_volume;
                    }

                    double maximum_probability = 0;
                    for (auto& p : destination_cells) {
                        maximum_probability = std::max(maximum_probability,p.second);
                    }

                    List<Pair<NCell,ProbabilityType>> pruned_cells;

                    double total_probability = 0;
                    for (auto& p : destination_cells) {
                        p.second = p.second/total_volume;
                        if (p.second >= _probability_threshold*maximum_probability) {
                            pruned_cells.append(p);
                            total_probability += p.second;
                        }
                    }

                    CONCLOG_PRINTLN_AT(1,"Probabilities for state " << source_cell.box() << " and control " << controller_cell.box() << ":")
                    for (auto& p : pruned_cells) {
                        p.second = p.second/total_probability;
                        CONCLOG_PRINTLN_AT(1,p.first.box() << ": " << p.second)
                    }

                    _reachability_graph->insert(source_cell, controller_cell, pruned_cells);
                }
            }
            indicator.update_current(indicator.current_value()+1.0);
        }
    }

    void refine_to_safety_graph() {
        _reachability_graph->reduce_to_not_reaching(_obstacles);
    }

    void refine_to_goal_reachable_graph() {
        _reachability_graph->reduce_to_possibly_reaching(_goals);
    }

    void update_unverified() {
        _reachability_graph->apply_source_removal_to(_unverified);
    }

    SizeType num_transitions() const { return _reachability_graph->num_transitions(); }

private:

    String const _name;
    EffectiveVectorMultivariateFunction const _dynamics;

    BoundsBoxType const _state_bounds;
    BoundsBoxType const _control_bounds;

    SPaving _state_paving;
    CPaving _control_paving;

    SizeType const _depth;

    FloatDP const _eps;
    ProbabilityType const _probability_threshold;

    SPaving _unverified;

    SPaving _obstacles;
    SPaving _goals;

    Ariadne::SharedPointer<ReachabilityGraphInterface> _reachability_graph;
};