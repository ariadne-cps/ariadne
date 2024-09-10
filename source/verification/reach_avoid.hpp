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

#include "function/function.hpp"
#include "utility/metaprogramming.hpp"

#include "reachability_graph.hpp"

namespace Ariadne {

template<class K, class V> using SizeTypeMap = Map<SizeType,Pair<K,V>>;
typedef Interval<double> BoundType;
typedef Vector<BoundType> BoundsBoxType;

ExactBoxType shrink(ExactBoxType const& bx, FloatDP const& eps);
ExactBoxType shrink(BoundsBoxType const& bx, FloatDP const& eps);

class ReachAvoid {
  public:
    ReachAvoid(String const& name, EffectiveVectorMultivariateFunction const& dynamics, Grid const& state_grid, BoundsBoxType const& state_bounds, Grid const& control_grid, BoundsBoxType const& control_bounds, SizeType depth, ExactDouble eps, ProbabilityType probability_threshold);

    ReachAvoid& add_obstacle(BoundsBoxType const& box);
    ReachAvoid& add_goal(BoundsBoxType const& box);

    BoundsBoxType const& state_bounds() const;
    BoundsBoxType const& control_bounds() const;

    SizeType state_size() const;
    SizeType control_size() const;
    SizeType obstacles_size() const;
    SizeType goals_size() const;
    SizeType unverified_size() const;

    SizeType num_sources() const;
    SizeType num_destinations() const;
    SizeType num_transitions() const;

    void plot(SizeType xaxis, SizeType yaxis) const;
    void plot() const;

    void print_goals() const;
    void print_obstacles() const;
    void print_graph() const;

    void compute_reachability_graph();

    void refine_to_safety_graph();

    void refine_to_goal_reachable_graph();

    void update_unverified();

    //! \brief The percentage (in the 0-100 scale) of still unverified states
    double unverified_percentage() const;

private:

    String const _name;
    EffectiveVectorMultivariateFunction const _dynamics;

    BoundsBoxType const _state_bounds;
    BoundsBoxType const _control_bounds;

    SPaving _state_paving;
    CPaving _control_paving;

    SizeType const _depth;

    FloatDP const _eps;
    ProbabilityType const _probability_threshold;

    SPaving _unverified;

    SPaving _obstacles;
    SPaving _goals;

    Ariadne::SharedPointer<ReachabilityGraphInterface> _reachability_graph;
};

} // namespace Ariadne