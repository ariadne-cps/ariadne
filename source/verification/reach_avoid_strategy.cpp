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

AssignedControl::AssignedControl(IdentifiedCell const& source_, IdentifiedCell const& control_, DirectionType const& direction_) :
    source(source_), control(control_), direction(direction_) { }

ReachAvoidStrategyBuilder::ReachAvoidStrategyBuilder(PossiblyReachingRAG const& rag) :
    _rag(rag) { }

ReachAvoidStrategy ReachAvoidStrategyBuilder::build() {
    List<AssignedControl> assignments;

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
            PointType best_midpoint(src.cell().dimension());
            for (auto const& ctrl : trans) {
                ScoreType current_score = 0.0;
                for (auto const& tgt : ctrl.second) {
                    if (sets_equidistant_to_goal.at(s_idx-1).contains(tgt.first))
                        current_score += scores.at(tgt.first) * tgt.second.probability();
                }
                if (current_score > best_score) {
                    best_control = ctrl.first;
                    best_score = current_score;
                    best_midpoint = PointType(src.cell().dimension());
                    for (SizeType i=0; i<best_midpoint.size(); ++i) {
                        for (auto const& tgt : ctrl.second) {
                            best_midpoint.at(i) = best_midpoint.at(i) + tgt.second.probability()*tgt.second.point().at(i);
                        }
                    }
                }
            }
            scores.insert(src,best_score);

            DirectionType direction(src.cell().dimension());
            auto src_midpoint = src.cell().box().midpoint();
            double direction_norm = 0.0;
            for (SizeType i=0; i<direction.size(); ++i) {
                auto current_midpoint_dimension = src_midpoint.at(i).get_d();
                direction.at(i) = best_midpoint.at(i) - current_midpoint_dimension;
                direction_norm += current_midpoint_dimension*current_midpoint_dimension;
            }
            direction_norm = std::sqrt(direction_norm);
            for (SizeType i=0; i<direction.size(); ++i) {
                direction.at(i) = direction.at(i)/direction_norm;
            }

            assignments.append({src,best_control,direction});
        }
    }

    return {assignments};
}

ReachAvoidStrategy::ReachAvoidStrategy(List<AssignedControl> const& assignments) :
    _assignments(assignments) { }

List<AssignedControl> const& ReachAvoidStrategy::assignments() const {
    return _assignments;
}

} // namespace Ariadne