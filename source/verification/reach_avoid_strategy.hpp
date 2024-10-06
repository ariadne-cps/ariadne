/***************************************************************************
 *            reach_avoid_strategy.hpp
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

#ifndef ARIADNE_REACH_AVOID_STRATEGY_HPP
#define ARIADNE_REACH_AVOID_STRATEGY_HPP

#include "reachability_graph.hpp"

namespace Ariadne {

typedef PointType DirectionType;
typedef double ScoreType;

struct AssignedControl {
    AssignedControl(IdentifiedCell const& source, IdentifiedCell const& control, DirectionType const& direction);

    IdentifiedCell const& source;
    IdentifiedCell const& control;
    DirectionType const& direction;

    friend OutputStream& operator<<(OutputStream& os, AssignedControl const& ac) { os << "{" << ac.source.id() << ":" << ac.control.id() << "@" << ac.direction << "}"; return os; }
};

class ReachAvoidStrategy;

class ReachAvoidStrategyBuilder {
public:
    ReachAvoidStrategyBuilder(PossiblyReachingRAG const& rag);
    ReachAvoidStrategy build();
private:
    PossiblyReachingRAG const _rag;
};

class ReachAvoidStrategy {
   friend class ReachAvoidStrategyBuilder;
  protected:
    ReachAvoidStrategy(List<AssignedControl> const& assignments);
  public:
    List<AssignedControl> const& assignments() const;
  private:
    List<AssignedControl> const _assignments;
};

} // namespace Ariadne

#endif