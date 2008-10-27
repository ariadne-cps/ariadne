/***************************************************************************
 *            discretiser.h
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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
 
/*! \file discretiser.h
 *  \brief Methods for computing the evolution of systems on grids/pavings.
 */

#ifndef ARIADNE_DISCRETISER_H
#define ARIADNE_DISCRETISER_H

#include <boost/smart_ptr.hpp>

#include "evolver_interface.h"
#include "discretiser_interface.h"

#include "logging.h"


namespace Ariadne {
  
class HybridGridCell;
class HybridGridCellListSet;

template<class SYS> class DiscretiserInterface;

typedef int DiscreteState;
class HybridAutomaton;

/*! \ingroup Evolvers
 *  \brief A class for computing the evolution of a discrete-time autonomous system.
 */
template<class ES>
class HybridDiscretiser
  : public DiscretiserInterface<HybridAutomaton>
  , public Loggable
{
  typedef HybridTime TimeType;
  typedef HybridAutomaton SystemType;
  typedef HybridGridCell BasicSetType;
  typedef HybridGridCellListSet DenotableSetType;
  typedef ES ContinuousEnclosureType;
  typedef std::pair<DiscreteState,ES> EnclosureType;
 private:
  boost::shared_ptr< EvolverInterface<SystemType,EnclosureType>  > _evolver;
 public:
  //@{
  //! \name Constructors and destructors
  
  //! \brief Construct from evolution parameters and a method for evolving basic sets, 
  //!  and a scheme for approximating sets.
  HybridDiscretiser(const EvolverInterface<SystemType,EnclosureType>& evolver)
    : _evolver(evolver.clone()) { }
      
  //! \brief Make a dynamically-allocated copy.
  HybridDiscretiser<ES>* clone() const { return new HybridDiscretiser<ES>(*this); }
  
  //@}
  
  //@{
  //! \name Evaluation on basic sets.
  
  //! \brief Compute lower approximations to the reachable and evolved sets 
  //! of \a system starting in \a initial_set over \a time. */
  virtual Orbit<BasicSetType> 
  lower_evolution(const SystemType& system, 
                  const BasicSetType& initial_set, 
                  const TimeType& time, const int accuracy) const;
  
  //! \brief Compute lower approximations to the reachable and evolved sets 
  //! of \a system starting in \a initial_set over \a time. */
  virtual Orbit<BasicSetType> 
  upper_evolution(const SystemType& system, 
                  const BasicSetType& initial_set, 
                  const TimeType& time, const int accuracy) const;

 private:
  EnclosureType _enclosure(const BasicSetType& bs) const;
  Orbit<BasicSetType> _discretise(const Orbit<EnclosureType>& orb,
                                  const BasicSetType& initial_set,
                                  const int accuracy) const;
  
};


} // namespace Ariadne


#endif /* ARIADNE_DISCRETISER_H */
