/***************************************************************************
 *            evolver_base.h
 *
 *  Copyright  2008  Pieter Collins
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

/*! \file evolver_base.h
 *  \brief Interface for computing a single time step of the evolution of a system.
 */

#ifndef ARIADNE_EVOLVER_BASE_H
#define ARIADNE_EVOLVER_BASE_H

#include "evolver_interface.h"
#include "list_set.h"

namespace Ariadne {


//! \brief Base class for common evolver functionality.
template<class SYS, class ES, class TRM> class EvolverBase
    : public EvolverInterface<SYS,ES,TRM>
{
    typedef EvolverInterface<SYS,ES,TRM> Interface;

  public:
    typedef typename EvolverInterface<SYS,ES,TRM>::SystemType SystemType;
    typedef typename EvolverInterface<SYS,ES,TRM>::TimeType TimeType;
    typedef typename EvolverInterface<SYS,ES,TRM>::TerminationType TerminationType;
    typedef typename EvolverInterface<SYS,ES,TRM>::EnclosureType EnclosureType;
    typedef typename EvolverInterface<SYS,ES,TRM>::EnclosureListType EnclosureListType;

  public:

    virtual OutputStream& write(OutputStream& os) const {
        return os << "Evolver( ... )"; }

  public:
    //@{
    //! \name Main evolution functions.

    virtual Orbit<EnclosureType> orbit(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics) const = 0;

    EnclosureListType evolve(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,termination,semantics,false); return final; }

    EnclosureListType reach(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,termination,semantics,true); return reachable; }

    Pair<EnclosureListType,EnclosureListType> reach_evolve(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,termination,semantics,true); return std::make_pair(reachable,final); }

    //@}

  protected:
    //! \brief Main routine for computing the evolution.
    virtual void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate, const EnclosureType& initial,
                            const TerminationType& termination, Semantics semantics, bool reach) const = 0;
};


} // namespace Ariadne



#endif // ARIADNE_EVOLVER_BASE_H
