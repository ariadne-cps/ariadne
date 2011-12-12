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


template<class SYS, class ES> class EvolverBase
    : public EvolverInterface<SYS,ES>
{
    typedef EvolverInterface<SYS,ES> Interface;
    typedef typename SYS::TimeType T;
    typedef ListSet<ES> ESL;
    typedef typename ESL::const_iterator ESLCI;
  public:

    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const {
        return os << "Evolver( ... )"; }

  public:
    //! \brief Compute an approximation to the evolution set under the given semantics.
    ESL evolve(const ES& initial_set, const T& time, Semantics semantics) const {
        ESL final; ESL reachable; ESL intermediate; this->_evolution(final,reachable,intermediate,initial_set,time,semantics,false); return final; }

    //! \brief Compute an approximation to the evolution set under the given semantics.
    ESL reach(const ES& initial_set, const T& time, Semantics semantics) const {
        ESL final; ESL reachable; ESL intermediate; this->_evolution(final,reachable,intermediate,initial_set,time,semantics,true); return reachable; }
    //! \brief Compute an approximation to the evolution set under the given semantics.
    std::pair<ESL,ESL> reach_evolve(const ES& initial_set, const T& time, Semantics semantics) const {
        ESL final; ESL reachable; ESL intermediate; this->_evolution(final,reachable,intermediate,initial_set,time,semantics,true); return std::make_pair(reachable,final); }

    //! \brief Compute an approximation to the evolution set under the given semantics.
    void evolution(ESL& final, const ES& initial, const T& time, Semantics semantics) const {
        ESL reachable; ESL intermediate; this->_evolution(final,reachable,intermediate,initial,time,semantics,false); }

    //! \brief Compute an approximation to the evolution set under the given semantics.
    void evolution(ESL& final, const ESL& initial, const T& time, Semantics semantics) const {
        ESL reachable; ESL intermediate; for(ESLCI iter=initial.begin(); iter!=initial.end(); ++iter) { this->_evolution(final,reachable,intermediate,ES(*iter),time,semantics,false); } }

    //! \brief Compute an approximation to the evolution set under the given semantics.
    void evolution(ESL& final, ESL& reachable, const ES& initial, const T& time, Semantics semantics) const {
        ESL intermediate; this->_evolution(final,reachable,intermediate,initial,time,semantics,true); }

    //! \brief Compute an approximation to the evolution set under the given semantics.
    void evolution(ESL& final, ESL& reachable, const ESL& initial, const T& time, Semantics semantics) const {
        ESL intermediate; for(ESLCI iter=initial.begin(); iter!=initial.end(); ++iter) { this->_evolution(final,reachable,intermediate,*iter,time,semantics,true); } }

    //! \brief Compute an approximation to the evolution set under the given semantics.
    void evolution(ESL& final, ESL& reachable, ESL& intermediate, const ESL& initial, const T& time, Semantics semantics) const {
        for(ESLCI iter=initial.begin(); iter!=initial.end(); ++iter) { this->_evolution(final,reachable,intermediate,*iter,time,semantics,true); } }
  protected:
    virtual void _evolution(ESL& final, ESL& reachable, ESL& intermediate, const ES& initial, const T& time, Semantics semantics, bool reach) const = 0;
};


} // namespace Ariadne



#endif // ARIADNE_EVOLVER_BASE_H
