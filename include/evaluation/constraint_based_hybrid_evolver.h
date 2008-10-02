/***************************************************************************
 *            constraint_based_hybrid_evolver.h
 *
 *  Copyright  2007-8  Pieter Collins
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
 
/*! \file constraint_based_hybrid_evolver.h
 *  \brief Evolver for hybrid systems defined using constraints.
 */

#ifndef ARIADNE_CONSTRAINT_BASED_HYBRID_EVOLVER_H
#define ARIADNE_CONSTRAINT_BASED_HYBRID_EVOLVER_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "base/tuple.h"

#include "system/hybrid_automaton.h"
#include "evaluation/evolver_interface.h"
#include "evaluation/evolver_base.h"

namespace Ariadne {  
  
    template<class Sys, class BS> class Evolver;

    template<class R> class ApproximateTaylorModel;
    template<class R> class ConstraintBasedHybridAutomaton;

    template<class R> class EvolutionParameters;
    template<class R> class DynamicalToolbox;
    
    class EvolutionProfiler;

    class HybridTime;

    template<class R> class Zonotope;
    template<class R> class TaylorSet;
    template<class BS> class HybridBasicSet;
    template<class BS> class TimeModelHybridBasicSet;
 
    enum CrossingKind { MISSING, TOUCHING, TRANSVERSE, CROSSING, GRAZING, UNKNOWN };

    /*! \ingroup Evolve 
     *  \brief A class for computing the evolution of a hybrid system. 
     *
     * The actual evolution steps are performed by the HybridEvolver class.
     */
    template< class R >
    class HybridEvolver
      : public EvolverBase< HybridAutomaton<R>, HybridBasicSet< TaylorSet<R> > >
    {
      typedef ApproximateTaylorModel<R> ModelType;
      typedef TimeModelHybridBasicSet< TaylorSet<R> > TimedSetType;
     public:
      typedef HybridAutomaton<R> SystemType;
      typedef HybridBasicSet< TaylorSet<R> > EnclosureType;
      typedef ListSet<EnclosureType> EnclosureListType;
      typedef HybridTime TimeType;
      typedef Rational ContinuousTimeType;
     public:

      /*! \brief The type used for real numbers. */
      typedef R RealType;

      /*! \brief Copy constructor. */
      HybridEvolver(const HybridEvolver<R>& evolver);

      /*! \brief Construct from parameters using a default integrator. */
      HybridEvolver(const EvolutionParameters<R>& parameters);

      /*! \brief Virtual destructor. */
      virtual ~HybridEvolver();


      //@{
      //! \name Parameters controlling the evolution.
      /*! \brief A reference to the parameters controlling the evolution. */
      EvolutionParameters<R>& parameters();
      const EvolutionParameters<R>& parameters() const;

      /*! \brief The maximum basic set radius before subdivision. */
      RealType maximum_enclosure_radius() const;
      /*! \brief The maximum step size for integration. */
      ContinuousTimeType maximum_step_size() const;
      /*! \brief The time before the sets are locked to the grid. */
      ContinuousTimeType lock_to_grid_time() const;

      //@}


      //@{
      //! \name Evolution using abstract sets.
      /*! \brief Compute an approximation to the evolution set using upper semantics. */
      EnclosureListType evolve(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate; 
        this->_evolution(final,reachable,intermediate,system,initial_set,time,upper_semantics,false); 
        return final; }
      /*! \brief Compute an approximation to the evolution set under upper semantics. */
      EnclosureListType reach(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate; 
        this->_evolution(final,reachable,intermediate,system,initial_set,time,upper_semantics,true); 
        return intermediate; }
     protected:
      virtual void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate, 
                              const SystemType& system, const EnclosureType& initial, const TimeType& time, 
                              Semantics semantics, bool reach) const;

      tuple< CrossingKind, ApproximateTaylorModel<R> > 
      _crossing(const ApproximateTaylorModel<R>& guard_model, const ApproximateTaylorModel<R>& flow_model);
     private:
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
      boost::shared_ptr< DynamicalToolbox<ModelType> > _toolbox;
      boost::shared_ptr< EvolutionProfiler >  _profiler;
   };


  
} // namespace Ariadne

#endif /* ARIADNE_CONSTRAINT_BASED_HYBRID_EVOLVER_H */
