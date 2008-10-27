/***************************************************************************
 *            discretiser_interface.h
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
 
/*! \file discretiser_interface.h
 *  \brief Interface for computing a single time step of the evolution of a system.
 */

#ifndef ARIADNE_DISCRETISER_INTERFACE_H
#define ARIADNE_DISCRETISER_INTERFACE_H


namespace Ariadne {
  
/*! \ingroup EvaluatorInterfaces \ingroup Evolvers
 *  \brief Interface for evolving a dynamic system and discretising the result on a grid.
 */
template<class SYS>
class DiscretiserInterface;

typedef int DiscreteState;
class HybridGridCell;

class HybridTime;
class HybridGridTreeSet;
class HybridAutomaton;

template<class ES> class Orbit;

/*! \ingroup EvaluatorInterfaces \ingroup Evolvers
 *  \brief Interface for evolving a hybrid dynamic system and discretising the result on a grid.
 */
template<>
class DiscretiserInterface<HybridAutomaton>
{
 public:
  typedef HybridTime TimeType;
  typedef HybridAutomaton SystemType;
  typedef HybridGridCell BasicSetType;
  typedef HybridGridTreeSet ConcreteSetType;
  
  /*! \brief Destructor. */
  virtual ~DiscretiserInterface() { };
  
  /*! \brief Make a dynamically-allocated copy. */
  virtual DiscretiserInterface<SystemType>* clone() const = 0;
  
  /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
  virtual Orbit<BasicSetType> 
  lower_evolve(const SystemType& f, const BasicSetType& s, const TimeType& t, const int accuracy) const = 0;
  
  /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
  virtual Orbit<BasicSetType> 
  upper_evolve(const SystemType& f, const BasicSetType& s, const TimeType& t, const int accuracy) const = 0;
};



  
} // namespace Ariadne



#endif /* ARIADNE_DISCRETISER_INTERFACE_H */
