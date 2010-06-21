/***************************************************************************
 *            hybrid_automaton_interface.h
 *
 *  Copyright  2010  Pieter Collins
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

/*! \file hybrid_automaton_interface.h
 *  \brief Interface for hybrid automata to be used by evolver classes.
 */

#ifndef ARIADNE_HYBRID_AUTOMATON_INTERFACE_H
#define ARIADNE_HYBRID_AUTOMATON_INTERFACE_H

#include "function.h"
#include "discrete_event.h"
#include "discrete_location.h"

namespace Ariadne {


class HybridTime;
class DiscreteEvent;
class DiscreteLocation;

class ScalarFunction;
class VectorFunction;
class Grid;

template<class T> class Set;

//! \brief Base interface for hybrid systems, to allow different types to be used in evolution routines.
class HybridAutomatonInterface {
  public:
    typedef HybridTime TimeType;
  public:
    //@{
    //! \name Data access and queries.

    //! \brief Virtual destructor.
    virtual ~HybridAutomatonInterface() { }

    //! \brief Test if the hybrid automaton has a discrete mode \a location.
    virtual bool has_mode(DiscreteLocation location) const = 0;

    //! \brief Test if the hybrid automaton has a discrete transition with \a event_id and \a source_id.
    virtual bool has_transition(DiscreteLocation source, DiscreteEvent event) const = 0;

    //! \brief The dimension of the state spacec in the given \a location.
    virtual uint dimension(DiscreteLocation location) const = 0;

    //! \brief The set of urgent events possible in the given \a location.
    virtual Set<DiscreteEvent> urgent_events(DiscreteLocation location) const = 0;

    //! \brief The set of permissive events possible in the given \a location.
    virtual Set<DiscreteEvent> permissive_events(DiscreteLocation location) const = 0;

    //! \brief The dynamic valid in the mode \a location.
    virtual VectorFunction dynamic_function(DiscreteLocation location) const = 0;

    //! \brief The invariants valid in the mode \a location.
    virtual ScalarFunction invariant_function(DiscreteLocation location, DiscreteEvent event) const = 0;

    //! \brief The guards active in the mode \a location.
    virtual ScalarFunction guard_function(DiscreteLocation location, DiscreteEvent event) const = 0;

    //! \brief The target location of \a event starting in the \a source location.
    virtual DiscreteLocation target(DiscreteLocation source, DiscreteEvent event) const = 0;

    //! \brief The dynamic valid in the mode \a location.
    virtual VectorFunction reset_function(DiscreteLocation location, DiscreteEvent event) const = 0;

    //! \brief The natural grid to use in the \a location.
    virtual Grid grid(DiscreteLocation location) const = 0;

    //@}

  public:
    virtual Set<DiscreteEvent> blocking_events(DiscreteLocation location) const = 0;
    virtual Set<DiscreteEvent> invariant_events(DiscreteLocation location) const = 0;
    virtual Set<DiscreteEvent> transition_events(DiscreteLocation location) const = 0;
    void reset(DiscreteEvent arg1);

};

} //namespace Ariadne

#endif // ARIADNE_HYBRID_AUTOMATON_INTERFACE_H
