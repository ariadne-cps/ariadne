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

class ReachAvoid {
public:
    ReachAvoid(String const& name, EffectiveVectorMultivariateFunction const& dynamics, Grid const& state_grid, BoundsBoxType const& state_bounds, Grid const& control_grid, BoundsBoxType const& control_bounds, ExactDouble eps) :
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

    ReachAvoid& add_obstacle(BoundsBoxType const& box) {
        SPaving obstacle_paving(_state_paving.grid());
        obstacle_paving.adjoin_outer_approximation(shrink(box,_eps),0);
        obstacle_paving.restrict(_state_paving);
        _unverified.remove(obstacle_paving);
        _obstacles.adjoin(obstacle_paving);
        return *this;
    }

    ReachAvoid& add_goal(BoundsBoxType const& box) {
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