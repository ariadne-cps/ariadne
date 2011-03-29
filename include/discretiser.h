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

#include "numeric.h"
#include "evolver_interface.h"
#include "discretiser_interface.h"
#include "hybrid_automaton.h"
#include "vector_field.h"
#include "evolution_parameters.h"

#include "logging.h"


namespace Ariadne {
  
class VectorField;
class Grid;
class GridCell;
class GridTreeSet;

class HybridAutomaton;
class DiscreteState;
template<class BS> class HybridBasicSet;

class HybridGrid;
class HybridGridCell;
class HybridGridTreeSet;


/*!  \brief A class for computing the evolution of a discrete-time autonomous system.
 */
template<class Sys, class ES>
class Discretiser
    : public Loggable
{
    typedef int AccuracyType;
    typedef Sys SystemType;
    typedef typename SystemType::TimeType TimeType;
    typedef GridCell BasicSetType;
    typedef GridTreeSet DenotableSetType;
    typedef ES EnclosureType;
  private:
    boost::shared_ptr< EvolverInterface<SystemType,EnclosureType>  > _evolver;
  public:
    //@{
    //! \name Constructors and destructors
  
    //! \brief Construct from evolution parameters and a method for evolving basic sets, 
    //!  and a scheme for approximating sets.
    Discretiser(const EvolverInterface<SystemType,EnclosureType>& evolver)
        : _evolver(evolver.clone()) { }

    /*! \brief Destructor. */
    virtual ~Discretiser() { }
      
    //! \brief Make a dynamically-allocated copy.
    Discretiser<Sys,ES>* clone() const { return new Discretiser<Sys,ES>(*this); }
  
    //@}
  
    //@{
    //! \name Evaluation on basic sets.
  
    //! \brief Compute approximations to the reachable and evolved sets 
    //! of \a system starting in \a initial_set over \a time. */
    virtual Orbit<BasicSetType> 
    evolution(const SystemType& system, 
              const BasicSetType& initial_set, 
              const TimeType& time, 
              const Grid& grid,
              const AccuracyType accuracy, 
              const Semantics semantics) const;
    
    /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
    virtual Orbit<BasicSetType> 
    lower_evolution(const SystemType& system, 
                    const BasicSetType& initial_set, 
                    const TimeType& time, 
                    const Grid& grid,
                    const AccuracyType accuracy) const;
  
    /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
    virtual Orbit<BasicSetType> 
    upper_evolution(const SystemType& system, 
                    const BasicSetType& initial_set, 
                    const TimeType& time, 
                    const Grid& grid,
                    const AccuracyType accuracy) const;  

  private:
    EnclosureType _enclosure(const BasicSetType& bs) const;

    Orbit<BasicSetType> _discretise(const Orbit<EnclosureType>& orb,
                                    const BasicSetType& initial_set,
                                    const Grid& grid,                                  
                                    const AccuracyType accuracy) const;
  
};


template<class ES>
HybridGridTreeSet 
outer_approximation(const ListSet<HybridBasicSet<ES> >& hls,
                    const HybridGrid& hgr,
                    const int accuracy);

template<class ES>
HybridGridTreeSet 
outer_approximation(const HybridBasicSet<ES>& hs,
                    const HybridGrid& hgr,
                    const int accuracy);

/*!  \brief A class for computing the evolution of a hybrid-time autonomous system.
 */
template<class ES>
class HybridDiscretiser
    : public Loggable
{
    typedef int AccuracyType;
    typedef HybridAutomaton::TimeType TimeType;
    typedef HybridAutomaton SystemType;
    typedef HybridGridCell BasicSetType;
    typedef HybridGridTreeSet DenotableSetType;
    typedef ES ContinuousEnclosureType;
    typedef HybridBasicSet<ES> EnclosureType;
	typedef ListSet<EnclosureType> EnclosureListType;
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
  
    //! \brief Return a shared pointer to the HybridEvolver.
    boost::shared_ptr< EvolverInterface<SystemType,EnclosureType> > evolver() const {
        return this->_evolver; }
  
    //@}

    //@{
    //! \name Gets and sets the continuous evolution parameters

	//! \brief Gets the evolution parameters from the evolver
	const ContinuousEvolutionParameters& parameters() const { return this->_evolver->parameters(); }

	//! \brief Gets a reference for setting the evolution parameters from the evolver
	ContinuousEvolutionParameters& parameters() { return this->_evolver->parameters(); }

	//@}

    //@{
    //! \name Evaluation on basic sets.
  
    //! \brief Compute approximations to the reachable and evolved sets 
    //! of \a system starting in \a initial_set over \a time. */
    virtual std::pair<DenotableSetType,DenotableSetType>
    evolution(const SystemType& system, 
              const EnclosureType& initial_set, 
              const TimeType& time,
              const AccuracyType accuracy,
              const Semantics semantics) const;

    //! \brief Compute approximations to the reachable and evolved sets 
    //! of \a system starting in \a initial_set over \a time, for upper 
	//! semantics with continuous evolution only. */
    virtual std::pair<DenotableSetType,DenotableSetType>
    upper_evolution_continuous(const SystemType& system, 
              			 	   const EnclosureType& initial_set, 
              			 	   const TimeType& time,
              			       const AccuracyType accuracy) const;

    //! \brief Compute approximations to the reachable set 
    //! of \a system starting in \a initial_set over \a time. */
    virtual DenotableSetType 
    reach(const SystemType& system, 
                const EnclosureType& initial_set, 
                const TimeType& time,
                const AccuracyType accuracy,
                const Semantics semantics) const;

    //! \brief Compute approximations to the evolved set 
    //! of \a system starting in \a initial_set over \a time. */
    virtual DenotableSetType 
    evolve(const SystemType& system, 
                 const EnclosureType& initial_set, 
                 const TimeType& time,
                 const AccuracyType accuracy,
                 const Semantics semantics) const;

    /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
    virtual std::pair<DenotableSetType,DenotableSetType>
    lower_evolution(const SystemType& system, 
                    const EnclosureType& initial_set, 
                    const TimeType& time,
                    const AccuracyType accuracy) const;
  
    /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
    virtual std::pair<DenotableSetType,DenotableSetType>
    upper_evolution(const SystemType& system, 
                    const EnclosureType& initial_set, 
                    const TimeType& time,
                    const AccuracyType accuracy) const; 

    /*! \brief Convert a cell of the grid into an enclosure set for computing evolution. */
    EnclosureType enclosure(const BasicSetType& bs) const;
    
  public:
    DenotableSetType _discretise(const ListSet<EnclosureType>& ls,
                                const HybridGrid& system_grid,
                                const AccuracyType accuracy) const; 
};

} // namespace Ariadne

#endif /* ARIADNE_DISCRETISER_H */
