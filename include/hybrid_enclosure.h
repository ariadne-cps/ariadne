/***************************************************************************
 *            hybrid_enclosure.h
 *
 *  Copyright  2009-10  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file hybrid_enclosure.h
 *  \brief Enclosure sets for hybrid systems
 */

#ifndef ARIADNE_HYBRID_ENCLOSURE_H
#define ARIADNE_HYBRID_ENCLOSURE_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>
#include "discrete_location.h"
#include "discrete_event.h"
#include "taylor_set.h"
#include "graphics_interface.h"
#include "container.h"
#include "box.h"
#include "logging.h"

namespace Ariadne {

class Interval;
template<class X> class Vector;
template<class X> class LinearProgram;
class ScalarFunction;
class VectorFunction;
class ScalarTaylorFunction;
class VectorTaylorFunction;
class TaylorConstrainedImageSet;
class HybridEnclosure;
class Box;
class Grid;
class GridTreeSet;
class AffineSet;
class DiscreteEvent;
class Figure;
class CanvasInterface;

template<class ES> class ListSet;
template<class ES> class HybridListSet;
template<> class ListSet<HybridEnclosure>;

template<class BS> class HybridBasicSet;
typedef HybridBasicSet<Box> HybridBox;
class HybridGridTreeSet;

typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;
typedef ScalarTaylorFunction ScalarIntervalFunction;
typedef VectorTaylorFunction VectorIntervalFunction;

//! \brief A class representing an enclosure for a hybrid evolution.
//! Handles progress, activation and guard constraints correctly.
//! The set is represented as the image of a box \f$D\f$ under a function model \f$\hat{f}(s)\f$, under the constraints
//! \f$\hat{c}(s) \leq 0\f$ and \f$\hat{e}(s)=0\f$. Also keeps track of the current time \f$\hat{t}(s)\f$.
//! In other words, \f[ S=\{ \hat{f}(s);\  \hat{t}(s) \mid s\in D \mid \hat{c}(s) \leq 0 \ \wedge \hat{e}(s)=0 \} . \f]
//! In the following documentation, we sometimes write \f$\xi(s)\f$ for \f$\hat{f}(s)\f$ and \f$\tau(s)\f$ for \f$\hat{t}(s)\f$.
class HybridEnclosure
    : public DrawableInterface
{
    friend class SimpleHybridEvolver;
    friend class ConstraintHybridEvolver;
  public:
    typedef TaylorConstrainedImageSet ContinuousStateSetType;
  private:
    DiscreteLocation _location;
    List<DiscreteEvent> _events;
    IntervalVector _domain;
    VectorIntervalFunction _state;
    ScalarIntervalFunction _time;
    List<ScalarIntervalFunction> _negative_constraints;
    List<ScalarIntervalFunction> _zero_constraints;

    mutable Box _reduced_domain;
  public:
    //! \brief A box in location \a q
    HybridEnclosure(const DiscreteLocation& q, const Box& s);
    //! \brief Construct from a continuous state set.
    HybridEnclosure(const std::pair<DiscreteLocation,ContinuousStateSetType>&);
    ////! \brief A set in location \a q, constructed from a continuous state set.
    //HybridEnclosure(const DiscreteLocation&, const ContinuousStateSetType&);
    //! \brief Destructor.
    ~HybridEnclosure();
    //! \brief Create a dynamically-allocated copy.
    HybridEnclosure* clone() const;

    //! \brief The current location.
    const DiscreteLocation& location() const;
    //! \brief The current location.
    const List<DiscreteEvent>& previous_events() const;
    //! \brief The number of independent parameters.
    uint number_of_parameters() const;
    //! \brief The number of constraints.
    uint number_of_constraints() const;
    //! \brief The continuous state set.
    const IntervalVector& parameter_domain() const;
    //! \brief The continuous state set.
    const VectorIntervalFunction& space_function() const;
    //! \brief The continuous state set.
    const ScalarIntervalFunction& time_function() const;

    //! \brief The continuous state set.
    ContinuousStateSetType continuous_state_set() const;

    //! \brief Apply the reset map \a r corresponding to event \a e with target location \a q.
    //! Corresponds to replacing \f$\xi\f$ by \f$r\circ \xi\f$.
    void apply_reset(DiscreteEvent e, DiscreteLocation q, VectorFunction r);

    //! \brief Apply the flow \a phi over the time interval up to time \a h, i.e. over \f$[0,h]\f$.
    //! Corresponds to replacing \f$D\f$ with \f$D\times [0,h]\f$, \f$\xi\f$ with
    //! \f$(s,t)\mapsto\phi(\xi(s),t)\f$, and \f$\tau\f$ with \f$(s,t)\mapsto\tau(s)+t\f$.
    void apply_flow(VectorIntervalFunction phi, Float h);
    //! \brief Apply the flow \a phi over the time interval ending at \a eps, i.e. \f$[0,\varepsilon(x)]\f$.
    //! Assume that \f$\eps(x)\in[0,h]\f$ for all \f$x\in S\f$.
    //! Corresponds to replacing \f$D\f$ with \f$D\times T\f$, \f$f\f$ with
    //! \f$(s,t)\mapsto\phi(f(s),t)\f$, and introducing a new constraint \f$t\leq \epsilon(f(s))\f$.
    void apply_flow(VectorIntervalFunction phi, ScalarIntervalFunction eps);
    //! \brief Apply the flow \a phi over the time interval \f$[0,\omega(s)-\tau(s)]\f$ so that the final time is \a omega.
    //! Corresponds to replacing \f$D\f$ with \f$D\times T\f$, \f$\xi\f$ with
    //! \f$(s,t)\mapsto\phi(\xi(s),\omega(s)-\tau(s))\f$, and \f$\tau\f$ by \f$(s,t)\mapsto\tau(s)+t\f$.
    void apply_flow_to(VectorIntervalFunction phi, ScalarIntervalFunction omega);

    //! \brief Apply the flow \a phi at the time \f$t\f$.
    //! Corresponds to setting \f$\xi'(s,t) = \phi(\xi(s),t)\f$ for \f$t\leq \omega(s)-\tau(s)\f$ with final time \f$\omega(s)\f$.
    void apply_flow_step(VectorIntervalFunction phi, Float h);
    //! \brief Apply the flow \a phi for the time \a eps.
    //! Corresponds to replacing \f$f\f$ with
    //! \f$s\mapsto\phi(\xi(s),\varepsilon(\xi(s)))\f$.
    void apply_flow_step(VectorIntervalFunction, ScalarIntervalFunction eps);
    //! \brief Apply the flow \a phi until the total elapsed time equals \a tau, where tau is a function on the parameter domain.
    //! Corresponds to replacing \f$\xi\f$ with
    //! \f$s\mapsto\phi(\xi(s),\omega(s)-\tau(s))\f$.
    void apply_flow_step_to(VectorIntervalFunction, ScalarIntervalFunction omega);

    //! \brief Set the maximum time of evolution to \a tmax. \deprecated
    //! Corresponds to introducting the constraint \f$\tau(s)\leq t_{\max}\f$.
    void set_maximum_time(DiscreteEvent e, Float tmax);
    //! \brief Set the time of evolution to \a tmax.
    //! Corresponds to introducting the constraint \f$\tau(s) = t\f$.
    void set_time(DiscreteEvent e, Float t);
    //! \brief Set the time of evolution to \a tmax.
    //! Corresponds to introducting the constraint \f$\tau(s) = t(s)\f$.
    void set_time(DiscreteEvent e, ScalarFunction t);
    //! \brief Introduces the constraint \f$\tau(s)\leq t_{\max}(s)\f$.
    void new_time_bound(DiscreteEvent e, ScalarFunction tmax);
    //! \brief Introduces the constraint \f$\tau(s)\leq t_{\max}(s)\f$.
    void new_time_bound(DiscreteEvent e, Float tmax);

    //! \brief Set the current time-step to \f$h\f$. \deprecated
    void set_step_time(Float h);
    //! \brief \deprecated
    void new_time_step_bound(DiscreteEvent e, ScalarIntervalFunction tau);

    //! \brief Introduces a new constraint \f$C\f$.
    void new_constraint(DiscreteEvent e, NonlinearConstraint c);
    //! \brief Introduces the new invariant (progress predicate) \f$c(x)\leq0\f$.
    void new_invariant(DiscreteEvent e, ScalarFunction g);
    //! \brief Introduces the new activation condition \f$g(x)\geq0\f$ for the event \a e.
    void new_activation(DiscreteEvent e,ScalarFunction g);
    //! \brief Introduces the new guard condition \f$g(x)=0\f$ for the event \a e.
    //! More precisely, the continuous dynamics is restricted to \f$c(x)\leq0\f$, and the event happens when \f$c(x)\geq0\f$.
    void new_guard(DiscreteEvent e, ScalarFunction g);
    //! \brief Introduces the new guard condition \f$g(x)=0\f$ for the event \a e, with computed crossing time \f$\tau(s)\f$.
    void new_guard(DiscreteEvent e, ScalarFunction g, ScalarIntervalFunction ct);

    //! \brief The dimension of the set.
    uint dimension() const;
    //! \brief Tests whether the set is empty.
    tribool empty() const;
    //! \brief Tests whether the set satisfies the constraint \a c.
    tribool satisfies(NonlinearConstraint c) const;

    //! \brief Returns a bounding box for the set. Computed by a simple interval evaluation of \f$f(D)\f$.
    HybridBox bounding_box() const;
    //! \brief Tests whether the set is disjoint from the box \a hbx.
    tribool disjoint(const HybridBox& hbx) const;
    //! \brief Tests whether the set is a subset of the box \a hbx.
    tribool subset(const HybridBox& hbx) const;
    //! \brief Adjoins an outer approximation of the set to the grid-based set \a paving, with accuracy given by
    //! \a depth subdivisions in each component.
    void adjoin_outer_approximation_to(HybridGridTreeSet& paving, int depth) const;
    //! \brief Splits into two smaller subsets.
    Pair<HybridEnclosure,HybridEnclosure> split(uint dim) const;

    //! \brief Draws onto a canvas.
    virtual void draw(CanvasInterface&) const;
    //! \brief Write to an output stream.
    std::ostream& write(std::ostream&) const;
  private:
    // Compute the flow reach step xi'(s,t) = phi(xi(s),t) and tau'(s,t)=tau(s)+t for t in [0,h] and t <= eps(xi(s)) .
    void _apply_flow(VectorIntervalFunction phi, Float step, ScalarIntervalFunction elps);
    // Compute the flow evolve step \xi'(s) = phi(xi(s),eps(s)) and tau'(s)=tau(s)+eps(s)
    void _apply_flow_step(VectorIntervalFunction phi, ScalarIntervalFunction elps);
    // Compute constraints of the set
    List<NonlinearConstraint> _constraints() const;

};

inline std::ostream& operator<<(std::ostream& os, const HybridEnclosure& s) { return s.write(os); }

}

#include "hybrid_set.h"

namespace Ariadne {

template<>
class ListSet<HybridEnclosure>
    : public HybridListSet<HybridEnclosure::ContinuousStateSetType>
{
  public:
    ListSet() { }
    ListSet(const HybridEnclosure& hes) { this->adjoin(hes); }
    using HybridListSet<HybridEnclosure::ContinuousStateSetType>::adjoin;
    void adjoin(const HybridEnclosure& hes) { this->adjoin(hes.location(),hes.continuous_state_set()); }
};

inline std::ostream& operator<<(std::ostream& os, const ListSet<HybridEnclosure>& hls) {
    return os << static_cast<const HybridListSet<HybridEnclosure::ContinuousStateSetType>&>(hls);
}

} // namespace Ariadne

#endif // ARIADNE_HYBRID_ENCLOSURE_H
