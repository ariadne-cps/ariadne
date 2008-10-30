/***************************************************************************
 *            reachability_analyser.h
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
 
/*! \file reachability_analyser.h
 *  \brief Methods for computing abstract reachable sets.
 */

#ifndef ARIADNE_REACHABILITY_ANALYSER_H
#define ARIADNE_REACHABILITY_ANALYSER_H

#include <boost/smart_ptr.hpp>


#include "hybrid_set_interface.h"

#include "evolver_interface.h"

#include "orbit.h"
#include "grid_set.h"
#include "hybrid_set.h"

#include "logging.h"


namespace Ariadne {
 
template<class ES> class Orbit;

typedef int DiscreteState;

class HybridGrid;
class HybridGridCell;
class HybridGridCellListSet;
class HybridGridTreeSet;

template<class ES> class HybridListSet;
class HybridEvolver;
template<class ES> class HybridDiscretiser;

typedef ApproximateTaylorModel DefaultEnclosureType;
typedef std::pair<DiscreteState,DefaultEnclosureType> DefaultHybridEnclosureType;



/*! \brief A class for performing reachability analysis on a hybrid system.
 *  \ingroup Analysers
 */
class HybridReachabilityAnalyser
  : public Loggable
{
 private:
  boost::shared_ptr< EvolutionParameters > _parameters;
  boost::shared_ptr< DiscreteEvolverInterface<HybridAutomaton,HybridGridCell> > _discretiser;
 public:
  typedef HybridAutomaton SystemType;
  typedef HybridTime TimeType;
  typedef HybridOpenSetInterface OpenSetType;
  typedef HybridOvertSetInterface OvertSetType;
  typedef HybridCompactSetInterface CompactSetType;
  typedef HybridRegularSetInterface RegularSetType;
  typedef HybridLocatedSetInterface LocatedSetType;
  typedef HybridGridTreeSet ConcreteSetType;
  typedef HybridBoxes BoundingSetType;
 public:
  //@{
  //! \name Constructors and destructors
  /*! \brief Virtual destructor */
  virtual ~HybridReachabilityAnalyser();

  /*! \brief Construct from evolution parameters and a method for evolving basic sets. */
  HybridReachabilityAnalyser(const EvolutionParameters& parameters, 
                             const EvolverInterface<HybridAutomaton,DefaultHybridEnclosureType>& evolver);

  /*! \brief Make a dynamically-allocated copy. */
  virtual HybridReachabilityAnalyser* clone() const { return new HybridReachabilityAnalyser(*this); }
  //@}
  
  //@{ 
  //! \name Methods to set and get the parameters controlling the accuracy.
  /*! \brief The parameters controlling the accuracy. */
  const EvolutionParameters& parameters() const { return *this->_parameters; }
  /*! \brief A reference to the parameters controlling the accuracy. */
  EvolutionParameters& parameters() { return *this->_parameters; }
  //@}
  
  
  //@{
  //! \name Evaluation of systems on abstract sets
  /*! \brief Compute a lower-approximation to the set obtained by evolving \a system for \a time starting in \a initial_set. */
  virtual HybridGridTreeSet* lower_evolve(const HybridAutomaton& system, 
                                          const HybridOvertSetInterface& initial_set, 
                                          const HybridTime& time) const;
  
  /*! \brief Compute a lower-approximation to the reachable set of \a system starting in \a initial_set up to \a time \a. */
  virtual HybridGridTreeSet*
  lower_reach(const HybridAutomaton& system, 
              const HybridOvertSetInterface& initial_set, 
              const HybridTime& time) const;
  
  /*! \brief Compute a lower-approximation to the reachable and evolved sets of \a system starting in \a initial_set up to \a time \a. */
  virtual std::pair<HybridGridTreeSet*,HybridGridTreeSet*>
  lower_reach_evolve(const HybridAutomaton& system, 
                     const HybridOvertSetInterface& initial_set, 
                     const HybridTime& time) const;
  
  /*! \brief Compute an approximation to the set obtained by iterating \a time times \a system starting in \a initial_set. */
  virtual HybridGridTreeSet* upper_evolve(const HybridAutomaton& system, 
                                         const HybridCompactSetInterface& initial_set, 
                                         const HybridTime& time) const;
  
  /*! \brief Compute an approximation to the reachable set of \a system starting in \a initial_set iterating at most \a time times. */
  virtual HybridGridTreeSet* upper_reach(const HybridAutomaton& system, 
                                        const HybridCompactSetInterface& initial_set, 
                                        const HybridTime& time) const;
  
  /*! \brief Compute an approximation to the reachable and evolved sets of \a system starting in \a initial_set iterating at most \a time times. */
  virtual std::pair<HybridGridTreeSet*,HybridGridTreeSet*>
  upper_reach_evolve(const HybridAutomaton& system, 
                     const HybridCompactSetInterface& initial_set, 
                     const HybridTime& time) const;
  
  /*! \brief Compute an outer-approximation to the chain-reachable set of \a system starting in \a initial_set. */
  virtual HybridGridTreeSet* chain_reach(const HybridAutomaton& system, 
                                         const HybridCompactSetInterface& initial_set, 
                                         const HybridBoxes& bounding_domain) const;
  
  /*! \brief Compute an outer-approximation to the viability kernel of \a system within \a bounding_set. */
  virtual HybridGridTreeSet* viable(const HybridAutomaton& system, 
                                   const HybridCompactSetInterface& bounding_set) const;
  
  /*! \brief Attempt to verify that the reachable set of \a system starting in \a initial_set remains in \a safe_set. */
  virtual tribool verify(const HybridAutomaton& system, 
                         const HybridLocatedSetInterface& initial_set, 
                         const HybridRegularSetInterface& safe_set) const;
  //@}
  
 public:
  typedef HybridTime T;
  typedef HybridAutomaton Sys;
  typedef HybridListSet<Box> BxLS;
  typedef HybridGrid Gr;
  typedef HybridGridCell GC;
  typedef HybridGridCellListSet GCLS;
  typedef HybridGridTreeSet GTS;
  typedef HybridOpenSetInterface OpSI;
  typedef HybridOvertSetInterface OvSI;
  typedef HybridCompactSetInterface CoSI;
 public:
  // Helper functions for operators on lists of sets.
  GTS _upper_reach(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
  GTS _upper_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
  std::pair<GTS,GTS> _upper_reach_evolve(const Sys& sys, const GTS& set, const T& time, const int accuracy) const;
 private:
  // Helper functions for approximating sets
  HybridGrid _grid(HybridSpace hspc) const { 
    return HybridGrid(hspc,this->_parameters->grid_length); }
};


} // namespace Ariadne


#endif // ARIADNE_REACHABILITY_ANALYSER_H
