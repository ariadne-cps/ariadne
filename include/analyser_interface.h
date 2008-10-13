/***************************************************************************
 *            analyser_interface.h
 *
 *  Copyright  2006-7  Pieter Collins
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
 
/*! \file analyser_interface.h
 *  \brief Interface for computing the evolution of sets.
 */

#ifndef ARIADNE_ANALYSER_INTERFACE_H
#define ARIADNE_ANALYSER_INTERFACE_H

#include <boost/smart_ptr.hpp>

#include "hybrid_set_interface.h"

namespace Ariadne {

  

/*! \brief Interface for computing (chain) reachable sets of a dynamic system.
 *  \ingroup EvaluatorInterfaces \ingroup Analysers
 */
template<class SYS> class AnalyserInterface;

template<> 
class AnalyserInterface<HybridAutomaton> 
{
 public:  
  typedef HybridAutomaton SystemType;
  typedef HybridTime TimeType;
  typedef HybridGridTreeSet ConcreteSetType;
 public:
  /*! \brief Virtual destructor. */
  virtual ~AnalyserInterface() { }
  
  //@{
  //! \name Evaluation of maps on abstract sets
  
  //! \brief Compute an approximation to the set obtained by iterating \a steps times \a system starting in \a initial_set.
  virtual HybridGridTreeSet* 
  lower_evolve(const HybridAutomaton& system, 
               const HybridOvertSetInterface& initial_set, 
               const HybridTime& steps) const = 0;
  
  //! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a steps times.
  virtual HybridGridTreeSet* 
  lower_reach(const HybridAutomaton& system, 
              const HybridOvertSetInterface& initial_set, 
              const HybridTime& steps) const = 0;
  
  //! \brief Compute an approximation to the set obtained by iterating \a steps times \a system starting in \a initial_set.
  virtual HybridGridTreeSet*
  upper_evolve(const HybridAutomaton& system, 
               const HybridCompactSetInterface& initial_set, 
               const HybridTime& steps) const = 0;
  
  //! \brief Compute an approximation to the reachable set 
  //! of \a system starting in \a initial_set iterating at most \a steps times.
  virtual HybridGridTreeSet* 
  upper_reach(const HybridAutomaton& system, 
              const HybridCompactSetInterface& initial_set, 
              const HybridTime& steps) const = 0;
  
  //! \brief Compute an outer-approximation to the chain-reachable set 
  //! of \a system starting in \a initial_set.
  virtual HybridGridTreeSet* 
  chain_reach(const HybridAutomaton& system, 
              const HybridCompactSetInterface& initial_set,
              const HybridBoxes& bounding_set) const = 0;
  
  //! \brief Compute an outer-approximation to the viability kernel 
  //! of \a system within \a bounding_set.
  virtual HybridGridTreeSet* viable(const HybridAutomaton& system, 
                                    const HybridCompactSetInterface& bounding_set) const = 0;
  
  //! \brief Attempt to verify that the reachable set 
  //! of \a system starting in \a initial_set remains in \a safe_set.
  virtual tribool verify(const HybridAutomaton& system, const HybridLocatedSetInterface& initial_set, const HybridRegularSetInterface& safe_set) const = 0;
  //@}
  
};




} // namespace Ariadne




#endif // ARIADNE_ANALYSER_INTERFACE_H
