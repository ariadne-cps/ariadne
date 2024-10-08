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

AssignedControl::AssignedControl(IdentifiedCell const& control, PointType const& target_point) :
    _control(control), _target_point(target_point) { }

IdentifiedCell const& AssignedControl::control() const {
    return _control;
}

PointType const& AssignedControl::target_point() const {
    return _target_point;
}

ReachAvoidStrategyBuilder::ReachAvoidStrategyBuilder(EffectiveVectorMultivariateFunction const& dynamics, PossiblyReachingRAG const& rag) :
    _dynamics(dynamics), _rag(rag) { }

double min_target_ratio(Box<FloatDPExactInterval> const& src_box, Box<FloatDPExactInterval> const& image_box) {
    double result = std::numeric_limits<double>::infinity();

    auto src_midpoint = src_box.midpoint();
    auto image_midpoint = image_box.midpoint();

    for (SizeType i=0; i<image_box.size(); ++i) {
        double delta = image_box[i].radius().get_d() + (src_box[i].radius().get_d() - (image_midpoint.at(i).get_d()-src_midpoint.at(i).get_d()));
        double current_value = std::abs(1.0 + delta/(image_midpoint.at(i).get_d() - src_midpoint.at(i).get_d()));
        if (current_value < result) {
            result = current_value;
        }
    }
    return result;
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

    auto sets_equidistant_to_goal = _rag.sets_equidistant_to_goal();

    Map<IdentifiedCell,ScoreType> scores;
    for (auto const& s : sets_equidistant_to_goal.at(0))
        scores.insert(s,1.0);

    for (SizeType s_idx = 1; s_idx < sets_equidistant_to_goal.size(); ++s_idx) {
        auto const& s = sets_equidistant_to_goal.at(s_idx);
        for (auto const& src : s) {
            auto const& trans = _rag.internal().forward_transitions(src);
            IdentifiedCell best_control = trans.begin()->first;
            ScoreType best_score = -static_cast<ScoreType>(trans.size());
            for (auto const& ctrl : trans) {
                ScoreType current_score = 0.0;
                for (auto const& tgt : ctrl.second) {
                    if (sets_equidistant_to_goal.at(s_idx-1).contains(tgt.first))
                        current_score += scores.at(tgt.first) * tgt.second.probability();
                }
                if (current_score > best_score) {
                    best_control = ctrl.first;
                    best_score = current_score;
                }
            }
            scores.insert(src,best_score);

            auto combined_input = product(src.cell().box(), best_control.cell().box());
            auto image_box = cast_exact_box(apply(_dynamics, combined_input).bounding_box());
            auto image_midpoint = image_box.midpoint();
            PointType image_point(image_midpoint.dimension());
            for (SizeType i=0; i<image_point.size(); ++i) {
                image_point.at(i) = image_midpoint.at(i).get_d();
            }

            auto src_box = src.cell().box();
            auto src_midpoint = src_box.midpoint();

            auto ratio = max_target_ratio(src_box,image_box);

            PointType target_point(image_point.size());
            for (SizeType i=0; i<target_point.size(); ++i)
                target_point.at(i) = src_midpoint.at(i).get_d()+(image_point.at(i)-src_midpoint.at(i).get_d())*ratio;

            assignments.insert(src,{best_control,target_point});
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