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

#include "reachability_graph.hpp"

namespace Ariadne {

class ReachAvoidStrategy;

class ReachAvoidStrategyBuilder {
  protected:
    friend class ReachAvoidStrategy;
  public:

    ReachAvoidStrategyBuilder(ReachabilityGraphInterface const& graph);

    ReachAvoidStrategy build();

  private:

    ReachabilityGraphInterface const& _graph;
};

class ReachAvoidStrategy {
  protected:

    ReachAvoidStrategy();
};

} // namespace Ariadne