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
//! Handles progress, activation and guard constraints internally.
//! The set is represented as the image of a box \f$D\f$ under a function model \f$\hat{f}(s)\f$, under the constraints
//! \f$\hat{c}(s) \leq 0\f$ and \f$\hat{e}(s)=0\f$. Also keeps track of the current time \f$\hat{t}(s)\f$.
//!
//! In other words,
//! \f[ S=\{ \hat{f}(s);\  \hat{t}(s) \mid s\in D \mid \hat{c}(s) \leq 0 \ \wedge \hat{e}(s)=0 \} . \f]
//! In the following documentation, we sometimes write \f$\xi(s)\f$ for \f$\hat{f}(s)\f$ and \f$\tau(s)\f$ for \f$\hat{t}(s)\f$.
//! yielding \f[ S=\{ \xi(s);\  \tau(s) \mid s\in D \mid \rho(s) \in C \} . \f]
//!
//! <b>Rationale:</b><\par>
//! <b>Parameterisation of the time specifier</b>
//! When computing a flow step, the evolution time can be given either as a
//! function \f$\varepsilon(x)\f$ of the state, or as an absolute final time
//! \f$\omega(s)\f$ of the parameters.
//! The rationale for this difference is that we often use the
//! \f$\varepsilon\f$ form when computing the flow up to a crossing with a
//! guard, in which case the time depends purely on the state, and it may
//! be possible to simplify the constraints by expressing them in terms of
//! the current states (in an implementation in which intermediate states are
//! considered explicitly).
//! The rationale for using a parameterised form with independent variable
//! \f$s\f$ is that this is needed to distinguished different states based on
//! crossings with the guard sets. The rationale for using a final time
//! \f$\omega(s)\f$ as a final time is that this allows us to remove
//! errors in the elapsed time \f$\tau(s)\f$.
//!
//! A possible extension would be to allow \f$\omega(s)\f$ to depend on
//! intermediate variables \f$x_i\f$ at previous stages. The reason for not
//! doing this is that it would enforce a policy for tracking intermediate
//! states, which we do not want to insist on during an implementation.
//!
//! <b>Replacing events</b>
//! If \f$g_e\f$ is increasing along trajectories
//! it is safe to remove a constraint \f$g_e(x_{i-1})\leq0\f$ and
//! replace it with the constraint \f$g_e(x_i)\leq 0\f$.
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
    //! \brief An empty enclosure.
    HybridEnclosure();
    //! \brief An enclosure corresponding to a box \a s in location \a q.
    HybridEnclosure(const DiscreteLocation& q, const Box& s);
    //! \brief An enclosure constructed from a continuous state set and a location.
    HybridEnclosure(const std::pair<DiscreteLocation,ContinuousStateSetType>&);
    ////! \brief A set in location \a q, constructed from a continuous state set.
    //HybridEnclosure(const DiscreteLocation&, const ContinuousStateSetType&);
    //! \brief Destructor.
    ~HybridEnclosure();
    //! \brief Create a dynamically-allocated copy.
    HybridEnclosure* clone() const;

    //! \brief The current location.
    const DiscreteLocation& location() const;
    //! \brief The list of previous events.
    const List<DiscreteEvent>& previous_events() const;
    //! \brief The number of independent parameters.
    uint number_of_parameters() const;
    //! \brief The number of constraints.
    uint number_of_constraints() const;
    //! \brief The continuous state set.
    const IntervalVector& parameter_domain() const;
    //! \brief The function related to space.
    const VectorIntervalFunction& space_function() const;
    //! \brief The function related to time.
    const ScalarIntervalFunction& time_function() const;

    //! \brief A bounding box for the space.
    IntervalVector space_bounding_box() const;
    //! \brief The range of times since the starting time that the set represents.
    Interval time_range() const;

    //! \brief The continuous state set.
    ContinuousStateSetType continuous_state_set() const;

    //! \brief Apply the reset map \a r corresponding to event \a e with target location \a q.
    //! Corresponds to replacing \f$\xi\f$ by \f$r\circ \xi\f$.
    void apply_reset(DiscreteEvent e, DiscreteLocation q, VectorFunction r);

    //! \brief Apply the flow \a phi over the time interval up to time \a h, i.e. over \f$[0,h]\f$.
    //! Corresponds to replacing \f$D\f$ with \f$D\times [0,h]\f$, \f$\xi\f$ with
    //! \f$(s,t)\mapsto\phi(\xi(s),t)\f$, and \f$\tau\f$ with \f$(s,t)\mapsto\tau(s)+t\f$.
    void apply_flow_for(VectorIntervalFunction phi, Float h);
    //! \brief Apply the flow \a phi over the time interval ending at \a eps, i.e. \f$[0,\varepsilon(x)]\f$.
    //! Assume that \f$\varepsilon(x)\in[0,h]\f$ for all \f$x\in S\f$.
    //! Corresponds to replacing \f$D\f$ with \f$D\times T\f$, \f$f\f$ with
    //! \f$(s,t)\mapsto\phi(f(s),t)\f$, and introducing a new constraint \f$t\leq \epsilon(f(s))\f$.
    //! The rationale for keeping the flow time in terms of the current state is to allow for an incremental
    //! approach to constructing the constraints.
    void apply_flow_for(VectorIntervalFunction phi, ScalarIntervalFunction eps);
    //! \brief Apply the flow \a phi over the time interval \f$[0,\omega(s)-\tau(s)]\f$ so that the final time is \a omega.
    //! Corresponds to replacing \f$D\f$ with \f$D\times T\f$, \f$\xi\f$ with
    //! \f$(s,t)\mapsto\phi(\xi(s),\omega(s)-\tau(s))\f$, and \f$\tau\f$ by \f$(s,t)\mapsto\tau(s)+t\f$.
    void apply_flow_to(VectorIntervalFunction phi, ScalarIntervalFunction omega);
    //! \brief Apply the flow \a phi over the time interval \f$[0,t_{\max}-\tau(s)]\f$ so that the final time bounded by a constant \a \f$t_{\max}\f$.
    //! Corresponds to replacing \f$D\f$ with \f$D\times T\f$, \f$\xi\f$ with
    //! \f$(s,t)\mapsto\phi(\xi(s),t)\f$, \f$\tau\f$ by \f$(s,t)\mapsto\tau(s)+t\f$ and constraint \f$\tau(s)+t\leq t_{\max}\f$.
    void apply_flow_to(VectorIntervalFunction phi, Float tmax);

    //! \brief Apply the flow \a phi at the time \f$h\f$.
    //! Corresponds to setting \f$\xi'(s) = \phi(\xi(s),h)\f$ and \f$\tau'(s) = \tau(s)+h\f$.
    void apply_flow_step_for(VectorIntervalFunction phi, Float h);
    //! \brief Apply the flow \a phi for the time \a eps.
    //! Corresponds to replacing \f$f\f$ with
    //! \f$s\mapsto\phi(\xi(s),\varepsilon(\xi(s)))\f$ and \f$\tau'(s) = \tau(s)+\varepsilon(\xi(s))\f$.
    //! The function \f$\varepsilon\f$ is required to satisfy \f$\varepsilon(x)\in[0,h]\f$ whenever \f$x\in S\f$.
    void apply_flow_step_for(VectorIntervalFunction phi, ScalarIntervalFunction eps);
    //! \brief Apply the flow \a phi until the total elapsed time equals \a omega, a function on the parameter domain.
    //! Corresponds to replacing \f$\xi\f$ with
    //! \f$s\mapsto\phi(\xi(s),\omega(s)-\tau(s))\f$ and setting \f$\tau'(s)=\omega(s)\f$.
    void apply_flow_step_to(VectorIntervalFunction phi, ScalarIntervalFunction omega);
    //! \brief Apply the flow \a phi until the total elapsed time equals a constant \a \f$t_{\max}\f$.
    //! Corresponds to replacing \f$\xi\f$ with
    //! \f$s\mapsto\phi(\xi(s),t_{\max}-\tau(s))\f$ and setting \f$\tau'(s)=t_{\max}\f$.
    void apply_flow_step_to(VectorIntervalFunction phi, Float tmax);

    //! \brief Set the time of evolution to \a \f$t_{\max}\f$.
    //! Corresponds to introducting the constraint \f$\tau(s) = t_{\max}\f$.
    void set_time(Float tmax);
    //! \brief Set the time of evolution to \a omega.
    //! Corresponds to introducting the constraint \f$\tau(s) = \omega(s)\f$.
    void set_time(ScalarFunction omega);
    //! \brief Introduces the constraint \f$\tau(s)\leq \omega(s)\f$.
    void bound_time(ScalarIntervalFunction omega);
    //! \brief Introduces the constraint \f$\tau(s)\leq \omega(s)\f$.
    void bound_time(ScalarFunction omega);
    //! \brief Introduces the constraint \f$\tau(s)\leq t_{\max}\f$.
    void bound_time(Float tmax);

    //! \brief Set the maximum time of evolution to \a \f$t_{\max}\f$. \deprecated
    //! Corresponds to introducting the constraint \f$\tau(s)\leq t_{\max}\f$.
    void set_maximum_time(DiscreteEvent event, Float tmax);
    //! \brief Set the current time-step to \f$h\f$. \deprecated
    void set_step_time(Float h);
    //! \brief \deprecated
    void new_time_step_bound(DiscreteEvent e, ScalarIntervalFunction tau);

    //! \brief Introduces a new state constraint \f$C\f$ on \f$x\f$. \deprecated
    void new_constraint(DiscreteEvent e, NonlinearConstraint c);
    //! \brief Introduces a new state constraint \f$C\f$ on \f$x\f$.
    void new_state_constraint(DiscreteEvent e, NonlinearConstraint c);
    //! \brief Introduces a new constraint \f$C\f$ on \f$s\f$.
    void new_parameter_constraint(DiscreteEvent e, NonlinearConstraint c);
    //! \brief Introduces the new invariant (progress predicate) \f$c(x)\leq0\f$.
    void new_invariant(DiscreteEvent e, ScalarFunction c);
    //! \brief Introduces the new invariant (progress predicate) \f$c(x)\leq0\f$.
    void new_invariant(DiscreteEvent e, ScalarIntervalFunction c);
    //! \brief Introduces the new activation condition \f$g(x)\geq0\f$ for the event \a e.
    void new_activation(DiscreteEvent e,ScalarFunction g);
    //! \brief Introduces the new guard condition \f$g(x)=0\f$ for the event \a e.
    //! More precisely, the continuous dynamics is restricted to \f$c(x)\leq0\f$, and the event happens when \f$c(x)\geq0\f$.
    void new_guard(DiscreteEvent e, ScalarFunction g);
    //! \brief Introduces the new guard condition \f$g(x)=0\f$ for the event \a e, with computed crossing time \f$\tau(s)\f$.
    void new_guard(DiscreteEvent e, ScalarFunction g, ScalarIntervalFunction ct);

    //! \brief Introduce a new independent variable with domain \a ivl.
    void new_variable(Interval ivl);

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
    // Compute the flow reach step xi'(s,t) = phi(xi(s),t) and tau'(s,t)=tau(s)+t for t in [0,h] .
    void _apply_flow(VectorIntervalFunction phi, Float step);
    // Compute the flow reach step xi'(s,t) = phi(xi(s),t) and tau'(s,t)=tau(s)+t for t in [0,h] and t <= eps(xi(s)) .
    void _apply_flow(VectorIntervalFunction phi, Float step, ScalarIntervalFunction elps);
    // Compute the flow evolve step \xi'(s) = phi(xi(s),eps(s)) and tau'(s)=tau(s)+eps(s)
    void _apply_flow_step(VectorIntervalFunction phi, ScalarIntervalFunction elps);
    // Compute constraints of the set
    List<NonlinearConstraint> constraints() const;

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
    ListSet(const List<HybridEnclosure>& hel) {
        for(List<HybridEnclosure>::const_iterator iter=hel.begin(); iter!=hel.end(); ++iter) { this->adjoin(*iter); } }
    using HybridListSet<HybridEnclosure::ContinuousStateSetType>::adjoin;
    void adjoin(const HybridEnclosure& hes) { this->adjoin(hes.location(),hes.continuous_state_set()); }
    void append(const HybridEnclosure& hes) { this->adjoin(hes.location(),hes.continuous_state_set()); }
    operator List<HybridEnclosure> () { return List<HybridEnclosure>(this->begin(),this->end()); }
};

inline std::ostream& operator<<(std::ostream& os, const ListSet<HybridEnclosure>& hls) {
    return os << static_cast<const HybridListSet<HybridEnclosure::ContinuousStateSetType>&>(hls);
}

inline tribool subset(const HybridEnclosure& e, const HybridBox& b) { return e.subset(b); }

} // namespace Ariadne

#endif // ARIADNE_HYBRID_ENCLOSURE_H
