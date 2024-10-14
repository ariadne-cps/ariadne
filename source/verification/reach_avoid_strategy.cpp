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

namespace Ariadne {

AssignedControl::AssignedControl(SPaving const& control_paving, IdentifiedCell const& target_cell) :
    _control_paving(control_paving), _target_cell(target_cell) { }

SPaving const& AssignedControl::control_paving() const {
    return _control_paving;
}

IdentifiedCell const& AssignedControl::target_cell() const {
    return _target_cell;
}

ReachAvoidStrategyBuilder::ReachAvoidStrategyBuilder(EffectiveVectorMultivariateFunction const& dynamics, PossiblyReachingRAG const& rag) :
    _dynamics(dynamics), _rag(rag) { }

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

    auto sets_equidistant_to_goal = _rag.sets_equidistant_to_goal();

    Map<IdentifiedCell,ScoreType> scores;
    for (auto const& s : sets_equidistant_to_goal.at(0))
        scores.insert(s,1.0);

    for (SizeType s_idx = 1; s_idx < sets_equidistant_to_goal.size(); ++s_idx) {
        auto const& s = sets_equidistant_to_goal.at(s_idx);
        for (auto const& src : s) {
            auto const& trans = _rag.internal().forward_transitions(src);

            SPaving control_grid(trans.begin()->first.cell().grid());
            Map<IdentifiedCell,ScoreType> target_scores;

            for (auto const& ctrl : trans) {
                control_grid.adjoin(ctrl.first.cell());
                for (auto const& tgt : ctrl.second) {
                    if (sets_equidistant_to_goal.at(s_idx-1).contains(tgt.first)) {
                        if (not target_scores.contains(tgt.first))
                            target_scores.insert(tgt.first,0.0);
                        target_scores.at(tgt.first) = target_scores.at(tgt.first) + tgt.second;
                    }
                }
            }

            auto best = target_scores.begin();
            for (auto it = target_scores.begin(); it != target_scores.end(); ++it) {
                if (it->second > best->second)
                    best = it;
            }

            assignments.insert(src,{control_grid,best->first});
        }
    }

    return {assignments};
}

ReachAvoidStrategy::ReachAvoidStrategy(Map<IdentifiedCell,AssignedControl> const& assignments) :
    _assignments(assignments) { }

Map<IdentifiedCell,AssignedControl> const& ReachAvoidStrategy::assignments() const {
    return _assignments;
}

} // namespace Ariadne