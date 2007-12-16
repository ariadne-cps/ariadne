/***************************************************************************
 *            set_based_hybrid_evolver.h
 *
 *  Copyright  2006-7  Alberto Casagrande,  Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
/*! \file set_based_hybrid_evolver.h
 *  \brief Evolver for hybrid systems defined using invariant and activation sets.
 */

#ifndef ARIADNE_SET_BASED_HYBRID_EVOLVER_H
#define ARIADNE_SET_BASED_HYBRID_EVOLVER_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "base/types.h"
#include "base/tribool.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

#include "system/set_based_hybrid_automaton.h"

namespace Ariadne {  
  namespace Evaluation {
  

    /*! \ingroup Evolve
     *  \brief A class for computing the evolution of a hybrid system.
     */
    template< class R >
    class SetBasedHybridEvolver
    {
      typedef Numeric::Interval<R> I;
     public:
      //@{
      //! \name Constructors and destructors
      /*! \brief Construct from evolution parameters. */
      SetBasedHybridEvolver(const EvolutionParameters<R>& parameters);

      /*! \brief Construct from evolution parameters, an applicator and an integrator. */
      template<class BS>
      SetBasedHybridEvolver(const EvolutionParameters<R>& parameters, const ApplicatorInterface<BS>& applicator, const IntegratorInterface<BS>& integrator);

      /*! \brief Construct from a map evolver and a vector field evolver. */
      SetBasedHybridEvolver(const MapEvolver<R>& applicator, const VectorFieldEvolver<R>& integrator);

       /*! \brief Copy constructor. */
      SetBasedHybridEvolver(const SetBasedHybridEvolver<R>& evolver);

     /*! \brief Virtual destructor. */
      virtual ~SetBasedHybridEvolver();
      //@}

      //@{
      //! \name Evolution using abstract sets.
      /*! \brief Make a discrete step of the hybrid automaton, starting from initial set. */
      Geometry::HybridSet<R> discrete_step(const System::SetBasedHybridAutomaton<R>& automaton, 
                                           const Geometry::HybridSet<R>& initial_set);
     
      /*! \brief Evolve the hybrid automaton starting from the \a initial_set respecting invariants and without using discrete transitions. */
      Geometry::HybridSet<R> continuous_chainreach(const System::SetBasedHybridAutomaton<R>& automaton, 
                                                   const Geometry::HybridSet<R>& initial_set);
     
      /*! \brief Compute a lower approximation to the reachable set using lower semantics. (Not currently implemented) */
      Geometry::HybridSet<R> lower_reach(const System::SetBasedHybridAutomaton<R>& automaton, 
                                         const Geometry::HybridSet<R>& initial_set);

      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. */
      Geometry::HybridSet<R> chainreach(const System::SetBasedHybridAutomaton<R>& automaton, 
                                        const Geometry::HybridSet<R>& initial_set);


      //@}
      
      //@{
      //! \name Deprecated methods for evolution using abstract sets.
      /*! \brief Evolve the hybrid automaton starting from the \a initial_set respecting invariants and without using discrete transitions. (Deprecated; use continuous_chainreach(automaton,initial_set) instead, and set the bounding box in the accuracy parameters.) */
      Geometry::HybridSet<R> continuous_chainreach(const System::SetBasedHybridAutomaton<R>& automaton, 
                                                   const Geometry::HybridSet<R>& initial_set,
                                                   const Geometry::HybridSet<R>& bounding_set);
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. (Deprecated; use continuous_chainreach(automaton,initial_set) instead, and set the bounding box in the accuracy parameters.)  */
      Geometry::HybridSet<R> chainreach(const System::SetBasedHybridAutomaton<R>& automaton, 
                                        const Geometry::HybridSet<R>& initial_set, 
                                        const Geometry::HybridSet<R>& bounding_set);
      //}

    public:
      //@{
      //! \name Evolution using concrete sets. (Not available through Python interface.)
      /*! \brief Make a discrete step of the hybrid automaton, starting from initial set. */
      Geometry::HybridGridMaskSet<R> discrete_step(const System::SetBasedHybridAutomaton<R>& automaton, 
                                                       const Geometry::HybridGridMaskSet<R>& initial_set);
      /*! \brief Evolve the hybrid automaton within \a bounding_set starting from the \a initial_set respecting invariants and without using discrete transitions. */
      Geometry::HybridGridMaskSet<R> continuous_chainreach(const System::SetBasedHybridAutomaton<R>& automaton, 
                                                           const Geometry::HybridGridMaskSet<R>& initial_set,
                                                           const Geometry::HybridGridMaskSet<R>& bounding_set);

     
      /*! \brief Compute an over approximation to the chain-reachable set using lower semantics. */
      Geometry::HybridGridCellListSet<R> lower_reach(const System::SetBasedHybridAutomaton<R>& automaton, 
                                                     const Geometry::HybridGridCellListSet<R>& initial_set);

      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. */
      Geometry::HybridGridMaskSet<R> chainreach(const System::SetBasedHybridAutomaton<R>& automaton, 
                                                const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                const Geometry::HybridGridMaskSet<R>& bounding_set);

      //@}
     private:
      // Perform one step of the discrete time evolution using lower semantics. (No checking of arguments)
      Geometry::HybridListSet< Geometry::Rectangle<R> > _discrete_step(const System::SetBasedHybridAutomaton<R>& automaton, 
                                                                       const Geometry::HybridListSet< Geometry::Rectangle<R> >& initial_set);
       // Evolve the hybrid automaton starting from the initial_set without using discrete transitions. (No checking of arguments)
      Geometry::HybridListSet< Geometry::Rectangle<R> > _continuous_reach(const System::SetBasedHybridAutomaton<R>& automaton, 
                                                                          const Geometry::HybridListSet< Geometry::Rectangle<R> >& initial_set,
                                                                          const Numeric::Rational& maximum_time);

      // Perform one step of the discrete time evolution using upper semantics. (No checking of arguments)
      Geometry::HybridGridMaskSet<R> _discrete_step(const System::SetBasedHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridGridMaskSet<R>& initial_set,
                                                    const Geometry::HybridGridMaskSet<R>& domain_set);
      // Evolve the hybrid automaton within \a domains starting from the initial_set without using discrete transitions (no checking).
      Geometry::HybridGridMaskSet<R> _continuous_chainreach(const System::SetBasedHybridAutomaton<R>& automaton, 
                                                            const Geometry::HybridGridMaskSet<R>& initial_set,
                                                            const Geometry::HybridGridMaskSet<R>& domain_set);
     private:
      MapEvolver<R>* _applicator;
      VectorFieldEvolver<R>* _integrator;
    };


  }
}

#include "set_based_hybrid_evolver.inline.h"

#endif /* ARIADNE_SET_BASED_HYBRID_EVOLVER_H */
