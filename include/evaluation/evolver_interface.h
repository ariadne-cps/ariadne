/***************************************************************************
 *            evolver_interface.h
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
 
/*! \file evolver_interface.h
 *  \brief Interface for computing a single time step of the evolution of a system.
 */

#ifndef ARIADNE_EVOLVER_INTERFACE_H
#define ARIADNE_EVOLVER_INTERFACE_H

#include "evaluation/semantics.h"

namespace Ariadne {
  
    /*! \ingroup EvaluatorInterfaces \ingroup Evolvers
     *  \brief Interface for evolving a dynamic system.
     *  
     */
    template<class Sys, class ES> class EvolverInterface 
    {
      typedef typename Sys::time_type T;
      typedef ListSet<ES> ESL;
     public:
      /*! \brief Destructor. */
      virtual ~EvolverInterface<Sys,ES>() {};
      /*! \brief Make a dynamically-allocated copy. */
      virtual EvolverInterface<Sys,ES>* clone() const = 0;
     public:
      /*! \brief Compute an approximation to the evolved set under the given semantics. */
      virtual ESL evolve(const Sys& system, const ES& initial_set, const T& time, Semantics semantics) const = 0;
      /*! \brief Compute an approximation to the reachable set under the given semantics. */
      virtual ESL reach(const Sys& system, const ES& initial_set, const T& time, Semantics semantics) const = 0;
      /*! \brief Compute an approximation to the evolved and reachable sets under the given semantics. */
      virtual std::pair<ESL,ESL> reach_evolve(const Sys& system, const ES& initial_set, const T& time, Semantics semantics) const = 0;


      /*! \brief Compute an approximation to the evolved set under the given semantics. */
      virtual void evolution(ESL& final, const Sys& system, const ES& initial, const T& time, Semantics semantics) const = 0;
      /*! \brief Compute an approximation to the evolved and reachable sets under the given semantics. */
      virtual void evolution(ESL& final, ESL& intermediate, const Sys& system, const ES& initial, const T& time, Semantics semantics) const = 0;

      /*! \brief Compute an approximation to the evolved set under the given semantics, starting from a list of enclosure sets. */
      virtual void evolution(ESL& final, const Sys& system, const ESL& initial, const T& time, Semantics semantics) const = 0;
      /*! \brief Compute an approximation to the evolved and reachable sets under the given semantics starting from a list of enclosure sets. */
      virtual void evolution(ESL& final, ESL& intermediate, const Sys& system, const ESL& initial, const T& time, Semantics semantics) const = 0;
    };



} // namespace Ariadne



#endif /* ARIADNE_EVOLVER_INTERFACE_H */
