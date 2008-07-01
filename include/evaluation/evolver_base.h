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

#include "evaluation/evolver_interface.h"

namespace Ariadne {
  
  
     template<class Sys, class ES> class EvolverBase
      : public EvolverInterface<Sys,ES>
    {
      typedef EvolverInterface<Sys,ES> Interface;
      typedef typename Sys::time_type T;
      typedef ListSet<ES> ESL;
      typedef typename ListSet<ES>::const_iterator ESLCI;
     public:
      /*! \brief Make a dynamically-allocated copy. */
      virtual EvolverBase<Sys,ES>* clone() const = 0;

      /*! \brief Compute an approximation to the evolution set under the given semantics. */
      ESL evolve(const Sys& system, const ES& initial_set, const T& time, Semantics semantics) const {
        ESL final; ESL intermediate; this->_evolution(final,intermediate,system,initial_set,time,semantics,false); return final; }
      /*! \brief Compute an approximation to the evolution set under the given semantics. */
      ESL reach(const Sys& system, const ES& initial_set, const T& time, Semantics semantics) const {
        ESL final; ESL intermediate; this->_evolution(final,intermediate,system,initial_set,time,semantics,true); return intermediate; }
      /*! \brief Compute an approximation to the evolution set under the given semantics. */
      std::pair<ESL,ESL> reach_evolve(const Sys& system, const ES& initial_set, const T& time, Semantics semantics) const {
        ESL final; ESL intermediate; this->_evolution(final,intermediate,system,initial_set,time,semantics,true); return std::make_pair(intermediate,final); }

      /*! \brief Compute an approximation to the evolution set under the given semantics. */
      void evolution(ESL& final, const Sys& system, const ES& initial, const T& time, Semantics semantics) const {
        ESL intermediate; this->_evolution(final,intermediate,system,initial,time,semantics,false); }
        
      /*! \brief Compute an approximation to the evolution set under the given semantics. */
      void evolution(ESL& final, const Sys& system, const ESL& initial, const T& time, Semantics semantics) const {
        ESL intermediate; for(ESLCI iter=initial.begin(); iter!=initial.end(); ++iter) { this->_evolution(final,intermediate,system,*iter,time,semantics,false); } }

      /*! \brief Compute an approximation to the evolution set under the given semantics. */
      void evolution(ESL& final, ESL& intermediate, const Sys& system, const ES& initial, const T& time, Semantics semantics) const {
        this->_evolution(final,intermediate,system,initial,time,semantics,true); }
        
      /*! \brief Compute an approximation to the evolution set under the given semantics. */
      void evolution(ESL& final, ESL& intermediate, const Sys& system, const ESL& initial, const T& time, Semantics semantics) const {
        for(ESLCI iter=initial.begin(); iter!=initial.end(); ++iter) { this->_evolution(final,intermediate,system,*iter,time,semantics,true); } }
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const {
        return os << "Evolver( ... )"; }
     protected:
      virtual void _evolution(ESL& final, ESL& intermediate, const Sys& system, const ES& initial, const T& time, Semantics semantics, bool reach) const = 0;
     };

  
} // namespace Ariadne



#endif /* ARIADNE_EVOLVER_BASE_H */
