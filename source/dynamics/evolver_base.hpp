/***************************************************************************
 *            dynamics/evolver_base.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file dynamics/evolver_base.hpp
 *  \brief Interface for computing a single time step of the evolution of a system.
 */

#ifndef ARIADNE_EVOLVER_BASE_HPP
#define ARIADNE_EVOLVER_BASE_HPP

#include "../dynamics/evolver_interface.hpp"
#include "../geometry/list_set.hpp"

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

    virtual OutputStream& _write(OutputStream& os) const {
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
    virtual Void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate, const EnclosureType& initial,
                            const TerminationType& termination, Semantics semantics, Bool reach) const = 0;
};


} // namespace Ariadne



#endif // ARIADNE_EVOLVER_BASE_HPP
