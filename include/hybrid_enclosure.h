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

//! \brief A class representing an enclosure for a hybrid evolution.
//! Handles progress, activation and guard constraints correctly.
//! The set is represented as the image of a box \f$D\f$ under a function model \f$\hat{f}\f$, under the constraints
//! \f$c(s) \leq 0\f$ and \f$e(x)=0\f$. Also keeps track of the current time \f$t=\tau(s)\f$.
class HybridEnclosure
    : public DrawableInterface
{
    friend class ConstraintHybridEvolver;
  public:
    typedef TaylorConstrainedImageSet ContinuousStateSetType;
    private:
    DiscreteLocation _location;
    List<DiscreteEvent> _events;
    List< List<DiscreteEvent> > _constraint_events;
    ContinuousStateSetType _set;
    ScalarTaylorFunction _time;
  public:
    //! \brief A box in location \a q
    HybridEnclosure(const DiscreteLocation& q, const Box& s);
    HybridEnclosure(const std::pair<DiscreteLocation,ContinuousStateSetType>&);
    //! \brief A set in location \a q, constructed from a continuous state set.
    HybridEnclosure(const DiscreteLocation&, const ContinuousStateSetType&);
    //! \brief Destructor.
    ~HybridEnclosure();
    //! \brief Creaqte a dynamically-allocated copy.
    HybridEnclosure* clone() const;

    //! \brief The current location.
    const DiscreteLocation& location() const;
    //! \brief The continuous state set.
    const ContinuousStateSetType& continuous_state_set() const;
    //! \brief The continuous state set.
    const VectorTaylorFunction& space_function() const;
    //! \brief The continuous state set.
    const ScalarTaylorFunction& time_function() const;

    //! \brief Apply the reset map \a r corresponding to event \a e with target location \a q.
    //! Corresponds to replacing \f$f\f$ by \f$r\circ f\f$.
    void apply_reset(DiscreteEvent e, DiscreteLocation q, VectorFunction r);
    //! \brief Apply the flow \a phi over the time interval \f$T\f$.
    //! Corresponds to replacing \f$D\f$ with \f$D\times T\f$, and \f$f\f$ with
    //! \f$(s,t)\mapsto\phi(f(s),t)\f$.
    void apply_flow(VectorFunction phi, Interval T);
    //! \brief Apply the flow \a phi at the time \f$t\f$.
    //! Corresponds to replacing \f$f\f$ with \f$s\mapsto\phi(f(s),t)\f$.
    void apply_flow(VectorTaylorFunction, Float);
    void apply_flow(VectorTaylorFunction, Interval);

    //! \brief Set the maximum time of evolution to \a tmax.
    //! Corresponds to introducting the constraint \f$\tau(s)\leq t_{\max}\f$.
    void set_maximum_time(DiscreteEvent e, Float tmax);
    //! \brief Set the current time-step to \f$h\f$.
    //! Only valid if the time function \f$\tau(s)=t_0+s_n\f$, in which case the variable \f$s_n\f$ is replaced by \f$t_0+c\f$.
    void set_step_time(Float);
    void set_time(DiscreteEvent,Float);
    void set_dwell_time(DiscreteEvent,Float);
    void set_dwell_time(DiscreteEvent,ScalarFunction);
    void set_dwell_time(DiscreteEvent,ScalarTaylorFunction);

    //! \brief Introduces the new invariant (progress predicate) \f$c(x)\leq0\f$.
    void new_invariant(DiscreteEvent e, ScalarFunction g);
    //! \brief Introduces the new activation condition \f$g(x)\geq0\f$ for the event \a e.
    void new_activation(DiscreteEvent e,ScalarFunction g);
    //! \brief Introduces the new guard condition \f$g(x)=0\f$ for the event \a e.
    //! More precisely, the continuous dynamics is restricted to \f$c(x)\leq0\f$, and the event happens when \f$c(x)\geq0\f$.
    void new_guard(DiscreteEvent e, ScalarFunction g);
    //! \brief Introduces the new guard condition \f$g(x)=0\f$ for the event \a e, with computed crossing time \f$\tau(s)\f$.
    void new_guard(DiscreteEvent e, ScalarFunction g, ScalarTaylorFunction ct);
    void new_time_bound(DiscreteEvent e, ScalarFunction tau);

    //! \brief The dimension of the set.
    uint dimension() const;
    //! \brief Tests whether the set is empty.
    tribool empty() const;

    //! \brief Returns a bounding box for the set. Computed by a simple interval evaluation of \f$f(D)\f$.
    HybridBox bounding_box() const;
    //! \brief Tests whether the set is disjoint from the box \a hbx.
    tribool disjoint(const HybridBox& hbx) const;
    //! \brief Adjoins an outer approximation of the set to the grid-based set \a paving, with accuracy given by
    //! \a depth subdivisions in each component.
    void adjoin_outer_approximation_to(HybridGridTreeSet& paving, int depth) const;
    //! \brief Splits into two smaller subsets.
    Pair<HybridEnclosure,HybridEnclosure> split(uint dim) const;

    //! \brief Draws onto a canvas.
    virtual void draw(CanvasInterface&) const;
    //! \brief Write to an output stream.
    std::ostream& write(std::ostream&) const;
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
