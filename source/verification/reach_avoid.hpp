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

#ifndef ARIADNE_REACH_AVOID_HPP
#define ARIADNE_REACH_AVOID_HPP

#include "function/function.hpp"
#include "utility/metaprogramming.hpp"

#include "reachability_graph.hpp"

namespace Ariadne {

template<class K, class V> using SizeTypeMap = Map<SizeType,Pair<K,V>>;

ExactBoxType shrink(ExactBoxType const& bx, FloatDP const& eps);
ExactBoxType shrink(RealBox const& bx, FloatDP const& eps);

class ReachAvoid {
  public:

    ReachAvoid(String const& name, EffectiveVectorMultivariateFunction const& dynamics, Grid const& state_grid, RealBox const& state_bounds, Grid const& control_grid, RealBox const& control_bounds, SizeType depth, ExactDouble eps);

    EffectiveVectorMultivariateFunction const& dynamics() const;

    Grid const& state_grid() const;
    Grid const& control_grid() const;
    SizeType grid_depth() const;

    IdentifiedCellFactory const& vertex_factory() const;
    IdentifiedCellFactory const& edge_factory() const;

    ReachAvoid& add_obstacle(RealBox const& box);
    ReachAvoid& add_goal(RealBox const& box);

    SPaving const& state_paving() const;
    CPaving const& control_paving() const;

    SPaving const& goals() const;
    SPaving const& obstacles() const;

    SPaving const& feasibles() const;
    SPaving const& unverified() const;

    RealBox const& state_bounds() const;
    RealBox const& control_bounds() const;

    SizeType state_size() const;
    SizeType control_size() const;
    SizeType obstacles_size() const;
    SizeType goals_size() const;
    SizeType unverified_size() const;

    SizeType num_sources() const;
    SizeType num_destinations() const;

    SizeType unconstrained_num_transitions() const;
    SizeType avoiding_num_transitions() const;
    SizeType possibly_reaching_num_transitions() const;

    void plot(SizeType xaxis, SizeType yaxis) const;
    void plot() const;

    void print_goals() const;
    void print_obstacles() const;

    void print_bounded_domain_graph() const;
    void print_avoiding_graph() const;
    void print_possibly_reaching_graph() const;

    //! \brief Return the goals that are still safe
    //! \return Returns the original goals if the avoid graph has not been computed yet
    SPaving safe_goals() const;

    void compute_bounded_domain_graph();

    void compute_avoiding_graph();

    void compute_possibly_reaching_graph();

    //! \brief The percentage (in the 0-100 scale) of still unverified states
    double unverified_percentage() const;

    PossiblyReachingRAG const& possibly_reaching_graph() const;

  private:

    SizeType _vertex_id(NCell const& cell) const;
    SizeType _edge_id(NCell const& cell) const;

  private:

    String const _name;
    EffectiveVectorMultivariateFunction const _dynamics;

    RealBox const _state_bounds;
    RealBox const _control_bounds;

    SPaving _state_paving;
    CPaving _control_paving;

    SizeType const _depth;

    FloatDP const _eps;

    SPaving _unverified;

    SPaving _obstacles;
    SPaving _goals;

    SharedPointer<IdentifiedCellFactory> _vertex_factory;
    SharedPointer<IdentifiedCellFactory> _edge_factory;

    BoundedDomainRAG _bounded_domain_graph;
    AvoidingRAG _avoiding_graph;
    PossiblyReachingRAG _possibly_reaching_graph;
};

} // namespace Ariadne

#endif