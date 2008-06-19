/***************************************************************************
 *            numerical_system.h
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
 
/*! \file numerical_system.h
 *  \brief A class under defined by orbits of systems.
 */

#ifndef ARIADNE_NUMERICAL_SYSTEM_H
#define ARIADNE_NUMERICAL_SYSTEM_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"

#include "geometry/euclidean_space.h"
#include "system/numerical_system_interface.h"
#include "evaluation/evolver_interface.h"

namespace Ariadne {
  

    /*! \ingroup NumericalSystem
     *  \brief A numerical evaluation of the evolution of a dynamic system on enclosure sets, exposing reach and evolve interfaces. 
     */
    template<class Sys, class ES>
    class NumericalSystem
      : public NumericalSystemInterface<typename Sys::time_type,ES>
    {
      typedef typename Sys::time_type T;
      typedef ListSet<ES> ESL;
     public:
      /*! \brief The type used to denote time. */
      typedef T time_type;
      /*! \brief The type of the system being discretised. */
      typedef Sys system_type;
      /*! \brief The type used to represent sets. */
      typedef ES enclosure_set_type;
     public:
      /*! \brief Construct from a system and an evolution algorithm. */
      NumericalSystem(const Sys& system, const EvolverInterface<Sys,ES>& evolver)
        : _system(system.clone()), _evolver(dynamic_cast< EvolverBase<Sys,ES>*>(evolver.clone())) { }
      /*! \brief Cloning operator. */
      virtual NumericalSystem<Sys,ES>* clone() const { 
        return new NumericalSystem<Sys,ES>(*this); }
      /*! \brief Return the discretised system. */
      const Sys& system() const { 
        return *this->_system; }
      /*! \brief The state space of the system. */
      virtual EuclideanSpace state_space() const { 
        return this->_system->state_space(); }

      /*! \brief Compute a lower-approximation to the evolved set under the system evolution. */
      virtual ESL evolve(const ES& initial, const T& time, Semantics semantics) const {
        return this->_evolver->evolve(this->system(),initial,time,semantics); }
      /*! \brief Compute a lower-approximation to the evolved set under the system evolution. */
      virtual ESL reach(const ES& initial, const T& time, Semantics semantics) const {
        return this->_evolver->reach(this->system(),initial,time,semantics); }
      /*! \brief Compute a lower-approximation to the reachable and evolved sets under the system evolution. */
      virtual std::pair<ESL,ESL> reach_evolve(const ES& initial, const T& time, Semantics semantics) const {
        return this->_evolver->reach_evolve(this->system(),initial,time,semantics); }
     private:
      boost::shared_ptr< Sys > _system;
      boost::shared_ptr< EvolverBase<Sys,ES> > _evolver;
    };

  
} // namespace Ariadne



#endif /* ARIADNE_NUMERICAL_SYSTEM_H */
