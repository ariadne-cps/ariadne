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
#include "function/function.hpp"

namespace Ariadne {

typedef double ScoreType;

struct AssignedControl {
    AssignedControl(IdentifiedCell const& control, PointType const& target_point);

    IdentifiedCell const& control() const;
    PointType const& target_point() const;

    friend OutputStream& operator<<(OutputStream& os, AssignedControl const& ac) { os << "{" << ac.control().id() << "@" << ac.target_point() << "}"; return os; }

  private:

    IdentifiedCell const _control;
    PointType const _target_point;
};

class ReachAvoidStrategy;

class ReachAvoidStrategyBuilder {
public:
    ReachAvoidStrategyBuilder(EffectiveVectorMultivariateFunction const& dynamics, PossiblyReachingRAG const& rag);
    ReachAvoidStrategy build();
private:
    EffectiveVectorMultivariateFunction const _dynamics;
    PossiblyReachingRAG const _rag;
};

class ReachAvoidStrategy {
   friend class ReachAvoidStrategyBuilder;
  protected:
    ReachAvoidStrategy(Map<IdentifiedCell,AssignedControl> const& assignments);
  public:
    Map<IdentifiedCell,AssignedControl> const& assignments() const;
  private:
    Map<IdentifiedCell,AssignedControl> const _assignments;
};

} // namespace Ariadne

#endif