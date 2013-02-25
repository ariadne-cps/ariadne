/***************************************************************************
 *            hybrid_evolver_interface.h
 *
 *  Copyright  2011  Luca Geretti
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

/*! \file hybrid_evolver_interface.h
 *  \brief Interface for evolver classes in the hybrid space.
 */

#ifndef ARIADNE_HYBRID_EVOLVER_INTERFACE_H
#define ARIADNE_HYBRID_EVOLVER_INTERFACE_H

#include "evolver_interface.h"

#include "hybrid_time.h"
#include "hybrid_set.h"
#include "hybrid_orbit.h"

#include "logging.h"

namespace Ariadne {

class HybridAutomatonInterface;
class HybridEnclosure;


//! \brief Interface for hybrid evolvers using HybridEnclosure as the enclosure type.
//! \details The class is loggable in order to allow verbosity tuning at the analyser layer.
class HybridEvolverInterface
    : public EvolverInterface<HybridAutomatonInterface,HybridEnclosure>
    , public Loggable
{
  public:
    //! \brief Make a dynamically-allocated copy.
    virtual HybridEvolverInterface* clone() const = 0;

    //@{
    //! \name Main evolution functions.

    //! \brief Compute an approximation to the orbit set using the given semantics, starting from an initial enclosure.
    //!   Useful for continuing a partially-computed orbit.
    virtual Orbit<EnclosureType> orbit(const EnclosureType& initial_enclosure,const TimeType& time,Semantics semantics) const = 0;
    //! \brief Compute an approximation to the orbit set using the given semantics, starting from a box.
    //!   Useful for computing the evolution starting from a cell of a grid.
    virtual Orbit<EnclosureType> orbit(const HybridBox& initial_box,const TimeType& time,Semantics semantics) const = 0;
    //! \brief Compute an approximation to the orbit set using the given semantics, starting from a set described by bounds and constraints.
    //!   Useful for computing the evolution starting from user-provided set.
    virtual Orbit<EnclosureType> orbit(const HybridSet& initial_set,const TimeType& time,Semantics semantics) const = 0;

    //@}

    //@{
    //! \name Auxiliary set conversion functionality

    //! \brief Set construct an enclosure from a box, such as one obtained from a grid.
    virtual EnclosureType enclosure(const HybridBox& initial_box) const = 0;
    //! \brief Set construct an enclosure from a user-provided set.
    virtual EnclosureType enclosure(const HybridSet& initial_set) const = 0;

    //@}

};


//! \brief Factory for hybrid evolver interface classes.
class HybridEvolverFactoryInterface
    : public EvolverFactoryInterface<HybridAutomatonInterface,HybridEnclosure>
{
  public:

    //! \brief Create a copy of the factory.
    virtual HybridEvolverFactoryInterface* clone() const = 0;
    //! \brief Create a hybrid evolver interface object around a \a system.
    virtual HybridEvolverInterface* create(const HybridAutomatonInterface& system) const = 0;
};


} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_INTERFACE_H
