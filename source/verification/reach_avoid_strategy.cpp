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

AssignedControl::AssignedControl(IdentifiedCell const& source_, IdentifiedCell const& control_) :
    source(source_), control(control_) { }

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
            for (auto const& ctrl : trans) {
                ScoreType current_score = 0.0;
                for (auto const& tgt : ctrl.second) {
                    if (sets_equidistant_to_goal.at(s_idx-1).contains(tgt.first))
                        current_score += scores.at(tgt.first) * tgt.second;
                }
                if (current_score > best_score) {
                    best_control = ctrl.first;
                    best_score = current_score;
                }
            }
            scores.insert(src,best_score);
            assignments.append({src,best_control});
        }
    }

    return ReachAvoidStrategy(assignments);
}

ReachAvoidStrategy::ReachAvoidStrategy(List<AssignedControl> const& assignments) :
    _assignments(assignments) { }

List<AssignedControl> const& ReachAvoidStrategy::assignments() const {
    return _assignments;
}

} // namespace Ariadne