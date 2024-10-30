/***************************************************************************
 *            reach_avoid_strategy.cpp
 *
 *  Copyright  2024  Luca Geretti
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

#include "reach_avoid_strategy.hpp"
#include "conclog/logging.hpp"

using namespace ConcLog;

namespace Ariadne {

AssignedControl::AssignedControl(IdentifiedCell const& control_cell, SPaving const& control_paving, IdentifiedCell const& target_cell) :
    _control_cell(control_cell), _control_paving(control_paving), _target_cell(target_cell) { }

IdentifiedCell const& AssignedControl::control_cell() const {
    return _control_cell;
}

SPaving const& AssignedControl::control_paving() const {
    return _control_paving;
}

IdentifiedCell const& AssignedControl::target_cell() const {
    return _target_cell;
}

ReachAvoidStrategyBuilder::ReachAvoidStrategyBuilder(EffectiveVectorMultivariateFunction const& dynamics, PossiblyReachingRAG const& rag) :
    _dynamics(dynamics), _rag(rag) {
}

double max_target_ratio(Box<FloatDPExactInterval> const& src_box, Box<FloatDPExactInterval> const& image_box) {
    double result = 0.0;

    auto src_midpoint = src_box.midpoint();
    auto image_midpoint = image_box.midpoint();

    for (SizeType i=0; i<image_box.size(); ++i) {
        double delta = image_box[i].radius().get_d() + (src_box[i].radius().get_d() - (image_midpoint.at(i).get_d()-src_midpoint.at(i).get_d()));
        double current_value = std::abs(1.0 + delta/(image_midpoint.at(i).get_d() - src_midpoint.at(i).get_d()));
        if (current_value > result) {
            result = current_value;
        }
    }
    return result;
}

ReachAvoidStrategy ReachAvoidStrategyBuilder::build() {
    Map<IdentifiedCell,AssignedControl> assignments;

    auto distances_to_goals = _rag.discrete_distances_to_goals();

    Map<IdentifiedCell,ScoreType> target_weights;
    for (auto const& d : distances_to_goals)
        target_weights.insert(d.first, 1.0/(1.0+static_cast<double>(d.second)));

    for (auto const& src : distances_to_goals.keys()) {
        auto const& trans = _rag.internal().forward_transitions(src);

        SPaving control_grid(trans.begin()->first.cell().grid());
        Map<IdentifiedCell,ScoreType> local_target_scores;
        Map<IdentifiedCell,ScoreType> local_control_scores;

        for (auto const& ctrl : trans) {
            local_control_scores.insert(ctrl.first, 0.0);
            control_grid.adjoin(ctrl.first.cell());
            for (auto const& tgt : ctrl.second) {
                ARIADNE_ASSERT_MSG(target_weights.contains(tgt.first),src.id() << " is not in the target weights")
                if (tgt.first.id() != src.id()) {
                    if (not local_target_scores.contains(tgt.first))
                        local_target_scores.insert(tgt.first, 0.0);
                    auto local_score = tgt.second.probability()*target_weights.at(tgt.first);
                    local_target_scores.at(tgt.first) = local_target_scores.at(tgt.first) + local_score;
                    local_control_scores.at(ctrl.first) = local_control_scores.at(ctrl.first) + local_score;
                }
            }
        }

        auto ctrl_it = local_control_scores.begin();
        auto best_control = ctrl_it;
        ++ctrl_it;
        for (; ctrl_it != local_control_scores.end(); ++ctrl_it) {
            if (ctrl_it->second > best_control->second)
                best_control = ctrl_it;
        }

        auto tgt_it = local_target_scores.begin();
        auto best_target = tgt_it;
        ++tgt_it;
        for (; tgt_it != local_target_scores.end(); ++tgt_it) {
            if (tgt_it->second > best_target->second)
                best_target = tgt_it;
        }

        assignments.insert(src,{best_control->first,control_grid,best_target->first});
    }

    return {assignments};
}

ReachAvoidStrategy::ReachAvoidStrategy(Map<IdentifiedCell,AssignedControl> const& assignments) :
    _assignments(assignments) { }

Map<IdentifiedCell,AssignedControl> const& ReachAvoidStrategy::assignments() const {
    return _assignments;
}

} // namespace Ariadne