/***************************************************************************
 *            reach_avoid.cpp
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

#include "geometry/paving_interface.hpp"
#include "geometry/grid_cell.hpp"
#include "geometry/grid_paving.hpp"

#include "io/figure.hpp"

#include "conclog/logging.hpp"
#include "conclog/progress_indicator.hpp"

#include "reach_avoid.hpp"

using namespace ConcLog;

namespace Ariadne {

ExactBoxType shrink(ExactBoxType const& bx, FloatDP const& eps) {
    ExactBoxType result(bx.dimension());
    for (SizeType i=0; i<bx.size(); ++i) {
        result[i] = Interval<FloatDP>((eps+bx[i].lower_bound()).lower_raw(),(-eps+bx[i].upper_bound()).upper_raw());
    }
    return result;
}

ExactBoxType shrink(RealBox const& bx, FloatDP const& eps) {
    ExactBoxType result(bx.size());
    for (SizeType i=0; i<bx.size(); ++i)
        result[i] = Interval<FloatDP>((FloatDP(cast_exact(bx[i].lower_bound().get_d()),DoublePrecision())+eps).lower_raw(),(FloatDP(cast_exact(bx[i].upper_bound().get_d()),DoublePrecision())-eps).upper_raw());
    return result;
}

ReachAvoid::ReachAvoid(String const& name, EffectiveVectorMultivariateFunction const& dynamics, Grid const& state_grid, RealBox const& state_bounds, Grid const& control_grid, RealBox const& control_bounds, SizeType depth, ExactDouble eps) :
            _name(name), _dynamics(dynamics), _state_bounds(state_bounds), _control_bounds(control_bounds), _depth(depth), _eps({eps,DoublePrecision()}) {

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

    _vertex_factory.reset(new IdentifiedCellFactory(default_vertex_extent,vertex_ids));
    _edge_factory.reset(new IdentifiedCellFactory(default_edge_extent,edge_ids));
}

EffectiveVectorMultivariateFunction const& ReachAvoid::dynamics() const {
    return _dynamics;
}

Grid const& ReachAvoid::state_grid() const {
    return _state_paving.grid();
}

Grid const& ReachAvoid::control_grid() const {
    return _control_paving.grid();
}

SizeType ReachAvoid::grid_depth() const {
    return _depth;
}

IdentifiedCellFactory const& ReachAvoid::vertex_factory() const {
    return *_vertex_factory;
}

IdentifiedCellFactory const& ReachAvoid::edge_factory() const {
    return *_edge_factory;
}

SPaving const& ReachAvoid::state_paving() const {
    return _state_paving;
}

SPaving const& ReachAvoid::control_paving() const {
    return _control_paving;
}

SPaving const& ReachAvoid::goals() const {
    return _goals;
}

SPaving const& ReachAvoid::obstacles() const {
    return _obstacles;
}

SPaving const& ReachAvoid::unverified() const {
    return _unverified;
}

ReachAvoid& ReachAvoid::add_obstacle(RealBox const& box) {
    SPaving obstacle_paving(_state_paving.grid());
    obstacle_paving.adjoin_outer_approximation(shrink(box,_eps),_depth);
    obstacle_paving.restrict(_state_paving);
    _unverified.remove(obstacle_paving);
    _obstacles.adjoin(obstacle_paving);
    return *this;
}

ReachAvoid& ReachAvoid::add_goal(RealBox const& box) {
    SPaving goal_paving(_state_paving.grid());
    goal_paving.adjoin_outer_approximation(shrink(box,_eps),_depth);
    goal_paving.restrict(_state_paving);
    _unverified.remove(goal_paving);
    _goals.adjoin(goal_paving);
    return *this;
}

RealBox const& ReachAvoid::state_bounds() const {
    return _state_bounds;
}

RealBox const& ReachAvoid::control_bounds() const {
    return _control_bounds;
}

SizeType ReachAvoid::state_size() const {
    return _state_paving.size();
}

SizeType ReachAvoid::control_size() const {
    return _control_paving.size();
}

SizeType ReachAvoid::obstacles_size() const {
    return _obstacles.size();
}

SizeType ReachAvoid::goals_size() const {
    return _goals.size();
}

SizeType ReachAvoid::unverified_size() const {
    return _unverified.size();
}

SizeType ReachAvoid::num_sources() const {
    if (_bounded_domain_graph == nullptr)
        return state_size();
    else if (_avoiding_graph == nullptr)
        return _bounded_domain_graph->num_sources();
    else if (_possibly_reaching_graph == nullptr)
        return _avoiding_graph->num_sources();
    else
        return _possibly_reaching_graph->num_sources();
}

SizeType ReachAvoid::num_destinations() const {
    if (_bounded_domain_graph == nullptr)
        return state_size();
    return _bounded_domain_graph->num_destinations();
}

SizeType ReachAvoid::_vertex_id(NCell const& cell) const {
    return _vertex_factory->create(cell).id();
}

SizeType ReachAvoid::_edge_id(NCell const& cell) const {
    return _edge_factory->create(cell).id();
}

double ReachAvoid::unverified_percentage() const {
    return static_cast<double>(_unverified.size()) * 100 / (_state_paving.size() - _goals.size() - _obstacles.size());
}

void ReachAvoid::plot(SizeType xaxis, SizeType yaxis) const {
    ExactBoxType graphics_box(_state_bounds.size());
    for (SizeType i=0; i<_state_bounds.size(); ++i)
        graphics_box[i] = FloatDPExactInterval({ExactDouble(_state_bounds[i].lower_bound().get_d()),DoublePrecision()},{ExactDouble(_state_bounds[i].upper_bound().get_d()),DoublePrecision()});

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

void ReachAvoid::plot() const {
    for (SizeType i=0; i<_dynamics.result_size()-1; ++i)
        for (SizeType j=i+1; j<_dynamics.result_size(); ++j)
            plot(i,j);
}

void ReachAvoid::print_goals() const {
    std::stringstream ss;
    for (auto const& g : _goals) ss << _vertex_id(g) << " ";
    CONCLOG_PRINTLN("Goals: " << ss.str())
}

void ReachAvoid::print_obstacles() const {
    std::stringstream ss;
    for (auto const& o : _obstacles) ss << _vertex_id(o) << " ";
    CONCLOG_PRINTLN("Obstacles: " << ss.str())
}

PointType normal_direction(FloatDPPoint const& src, FloatDPPoint const& dst) {
    auto diff = dst-src;
    auto diff_norm = norm(diff).get_d();
    PointType result(diff.size());
    for (SizeType d=0; d < diff.size(); ++d)
        result.at(d) = diff.at(d).get_d()/diff_norm;
    return result;
}

double scalar_product(PointType const& v1, PointType const& v2) {
    double result = 0;
    for (SizeType d=0; d < v1.size(); ++d)
        result += v1.at(d)*v2.at(d);
    return result;
}

void ReachAvoid::compute_bounded_domain_graph() {
    CONCLOG_SCOPE_CREATE

    ForwardBackwardReachabilityGraph result(*_vertex_factory,*_edge_factory);

    SPaving excluded_border_cells(state_grid());

    for (auto const& source_cell : _state_paving) {
        bool has_feasible_controls = false;
        for (auto const& controller_cell : _control_paving) {
            auto combined = product(source_cell.box(), controller_cell.box());
            SPaving destination_paving(_state_paving.grid());
            auto image_box = shrink(cast_exact_box(apply(_dynamics, combined).bounding_box()), _eps);
            destination_paving.adjoin_outer_approximation(image_box, _depth);
            destination_paving.mince(_depth);
            SizeType original_size = destination_paving.size();
            destination_paving.restrict(_state_paving);
            SizeType restricted_size = destination_paving.size();
            if (restricted_size == original_size) {
                has_feasible_controls = true;
                break;
            }
        }
        if (not has_feasible_controls)
            excluded_border_cells.adjoin(source_cell);
    }

    SPaving feasibles = _state_paving;
    feasibles.remove(excluded_border_cells);

    ProgressIndicator indicator(static_cast<double>(feasibles.size()));
    for (auto const& source_cell : feasibles) {
        auto source_point = source_cell.box().midpoint();
        CONCLOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "% ")

        for (auto const& controller_cell : _control_paving) {
            auto combined = product(source_cell.box(), controller_cell.box());
            SPaving destination_paving(state_grid());
            auto image_box = shrink(cast_exact_box(apply(_dynamics, combined).bounding_box()), _eps);
            auto image_point = image_box.midpoint();

            auto image_dir = normal_direction(source_point,image_point);

            destination_paving.adjoin_outer_approximation(image_box, _depth);
            destination_paving.mince(_depth);
            destination_paving.restrict(feasibles);

            if (not destination_paving.is_empty()) {
                List<Pair<NCell,TargetScore>> destination_data;

                List<ProbabilityType> volumes;
                List<AlignmentType> alignments;
                ProbabilityType total_volume = 0.0;
                for (auto const& destination_cell : destination_paving) {
                    auto destination_dir = normal_direction(source_point,destination_cell.box().midpoint());
                    alignments.append(scalar_product(image_dir,destination_dir));
                    auto current_volume = intersection(image_box, destination_cell.box()).volume().get_d();
                    volumes.append(current_volume);
                    total_volume += current_volume;
                }

                SizeType i=0;
                for (auto const& destination_cell : destination_paving) {
                    destination_data.append({destination_cell, {volumes.at(i)/total_volume,alignments.at(i)}});
                    ++i;
                }

                CONCLOG_PRINTLN_AT(1,"Probabilities & alignments for state " << source_cell.box() << " and control " << controller_cell.box() << ":")
                for (auto& p : destination_data) {
                    CONCLOG_PRINTLN_AT(1,p.first.box() << ": " << p.second)
                }

                result.insert(source_cell, controller_cell, destination_data);
            }
        }

        indicator.update_current(indicator.current_value()+1.0);
    }

    _bounded_domain_graph.reset(new BoundedDomainRAG(result));
}

void ReachAvoid::compute_avoiding_graph() {
    SPaving remaining_obstacles = _obstacles;
    _bounded_domain_graph->apply_source_restriction_to(remaining_obstacles);
    _avoiding_graph.reset(_bounded_domain_graph->reduce_to_not_reaching(remaining_obstacles).clone());
}

SPaving ReachAvoid::safe_goals() const {
    SPaving result = _goals;
    if (not _avoiding_graph->is_empty())
        _avoiding_graph->apply_source_restriction_to(result);
    return result;
}

void ReachAvoid::compute_possibly_reaching_graph() {
    _possibly_reaching_graph.reset(_avoiding_graph->reduce_to_possibly_reaching(safe_goals()).clone());
    _possibly_reaching_graph->apply_source_removal_to(_unverified);
}

SizeType ReachAvoid::unconstrained_num_transitions() const {
    return _bounded_domain_graph->internal().num_transitions();
}

SizeType ReachAvoid::avoiding_num_transitions() const {
    return _avoiding_graph->internal().num_transitions();
}

SizeType ReachAvoid::possibly_reaching_num_transitions() const {
    return _possibly_reaching_graph->internal().num_transitions();
}

BoundedDomainRAG const& ReachAvoid::bounded_domain_graph() const {
    return *_bounded_domain_graph;
}

AvoidingRAG const& ReachAvoid::avoiding_graph() const {
    return *_avoiding_graph;
}

PossiblyReachingRAG const& ReachAvoid::possibly_reaching_graph() const {
    return *_possibly_reaching_graph;
}

}