/***************************************************************************
 *            constraint_based_hybrid_evolver.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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

#include "base/types.h"
#include "base/tribool.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

#include "geometry/discrete_state.h"
#include "system/discrete_event.h"
#include "evaluation/hybrid_time.h"

namespace Ariadne {  
  
  
    template<class R> class EvolutionParameters;
    template<class R> class ConstraintBasedHybridScheduler;


    /*! \ingroup Evolve 
     *  \brief A class for computing the evolution of a hybrid system. 
     *
     * The actual evolution steps are performed by the HybridEvolver class.
     */
    template< class R >
    class ConstraintBasedHybridEvolver
    {
      typedef Interval<R> I;
      typedef Zonotope<R,IntervalTag> BS;
     public:
      typedef ConstraintBasedDiscreteMode<R> mode_type;
      typedef ConstraintBasedDiscreteTransition<R> transition_type;
      typedef DiscreteState discrete_state_type;
      typedef Zonotope<R,IntervalTag> continuous_basic_set_type;
      typedef HybridBasicSet<continuous_basic_set_type> hybrid_basic_set_type;
      //typedef HybridTimedBasicSet<continuous_basic_set_type> timed_set_type;
      typedef TimeModelHybridBasicSet<continuous_basic_set_type> timed_set_type;
      typedef HybridListSet<continuous_basic_set_type> hybrid_list_set_type;
      typedef ConstraintInterface<R> constraint_type;
      typedef Box<R> bounding_box_type;
      typedef TimeModel<R> time_model_type;

      typedef boost::shared_ptr< const ConstraintInterface<R> > constraint_const_pointer;
     public:

      /*! \brief The type used for real numbers. */
      typedef R real_type;

      /*! \brief Copy constructor. */
      ConstraintBasedHybridEvolver(const ConstraintBasedHybridEvolver<R>& evolver);

      /*! \brief Construct from parameters using a default integrator. */
      ConstraintBasedHybridEvolver(const EvolutionParameters<R>& parameters);

      /*! \brief Construct from an applicator and an integrator. (Deprecated) */
      ConstraintBasedHybridEvolver(const EvolutionParameters<R>& parameters, const ApplicatorInterface<BS>& applicator, const IntegratorInterface<BS>& integrator);

      /*! \brief Construct from an applicator, an integrator and a detector. */
      ConstraintBasedHybridEvolver(const EvolutionParameters<R>& parameters, const ApplicatorInterface<BS>& applicator, const IntegratorInterface<BS>& integrator, const DetectorInterface<R>& detector);

      /*! \brief Virtual destructor. */
      virtual ~ConstraintBasedHybridEvolver();


      //@{
      //! \name Parameters controlling the evolution.
      /*! \brief A reference to the parameters controlling the evolution. */
      EvolutionParameters<R>& parameters();
      const EvolutionParameters<R>& parameters() const;

      /*! \brief The maximum step size for integration. */
      time_type maximum_step_size() const;
      /*! \brief The maximum basic set radius before subdivision. */
      real_type maximum_basic_set_radius() const;
      /*! \brief The time before the sets are locked to the grid. */
      time_type lock_to_grid_time() const;

      //@}


      //@{
      //! \name Evolution using abstract sets.
      /*! \brief Make a discrete step of the hybrid automaton, starting from initial set. (NOT CURRENTLY IMPLEMENTED) */
      HybridSet<R> discrete_step(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                           const HybridSet<R>& initial_set) const;
     
      /*! \brief Evolve the hybrid automaton within \a bounding_set starting from the \a initial_set respecting invariants and without using discrete transitions. (NOT CURRENTLY IMPLEMENTED) */
      HybridSet<R> continuous_chainreach(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                                   const HybridSet<R>& initial_set,
                                                   const HybridSet<R>& bounding_set) const;
     

      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      HybridSet<R> lower_evolve(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                          const HybridSet<R>&,
                                          time_type evolution_time,
                                          size_type maximum_number_of_events) const;
      
      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      HybridSet<R> upper_evolve(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                          const HybridSet<R>&,
                                          time_type evolution_time,
                                          size_type maximum_number_of_events) const;
      
      /*! \brief Compute a lower approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      HybridSet<R> lower_reach(const ConstraintBasedHybridAutomaton<R>&, 
                                         const HybridSet<R>&, 
                                         time_type initial_evolution_time, 
                                         time_type final_time, 
                                         size_type maximum_number_of_events) const;
      
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      HybridSet<R> upper_reach(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                         const HybridSet<R>& initial_set, 
                                         time_type initial_evolution_time, 
                                         time_type final_time, 
                                         size_type maximum_number_of_events) const;
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. */
      HybridSet<R> chainreach(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                        const HybridSet<R>& initial_set, 
                                        const HybridSet<R>& bounding_set) const;

      /*! \brief Compute the viability kernel of \a map within \a bounding_set. (NOT CURRENTLY IMPLEMENTED) */
      HybridSet<R> viable(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                    const HybridSet<R>& bounding_set) const;
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      tribool verify(const ConstraintBasedHybridAutomaton<R>& automaton, 
                     const HybridSet<R>& initial_set, 
                     const HybridSet<R>& safe_set) const;
      //@}

     public:
      //@{
      //! \name Evolution using concrete sets.
      /*! \brief Make a discrete step of the hybrid automaton, starting from initial set. */
      HybridGridMaskSet<R> discrete_step(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                                       const HybridGridMaskSet<R>& initial_set) const;
      /*! \brief Evolve the hybrid automaton within \a bounding_set starting from the \a initial_set respecting invariants and without using discrete transitions. */
      HybridGridMaskSet<R> continuous_chainreach(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                                           const HybridGridMaskSet<R>& initial_set,
                                                           const HybridGridMaskSet<R>& bounding_set) const;


      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using lower semantics. */
      HybridListSet<continuous_basic_set_type> lower_evolve(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                                                      const HybridListSet<continuous_basic_set_type>& initial_set, 
                                                                      time_type evolution_time,
                                                                      size_type maximum_number_of_events) const;
      
      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using upper semantics. */
      HybridListSet<continuous_basic_set_type> upper_evolve(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                                                      const HybridListSet<continuous_basic_set_type>& initial_set, 
                                                                      time_type evolution_time,
                                                                      size_type maximum_number_of_events) const;
      
      /*! \brief Compute the system evolution up to \a time with up to \a maximum_number_of_events using lower semantics. */
      HybridListSet<continuous_basic_set_type> lower_reach(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                                                      const HybridListSet<continuous_basic_set_type>& initial_set, 
                                                                      time_type evolution_time,
                                                                      size_type maximum_number_of_events) const;

      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics.*/
      HybridListSet<continuous_basic_set_type> 
      upper_reach(const ConstraintBasedHybridAutomaton<R>& automaton, 
                  const HybridListSet<continuous_basic_set_type>& initial_set, 
                  time_type evolution_time, 
                  size_type maximum_number_of_events) const;
      

     
     
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
      *  with up to \a maximum_number_of_events using upper semantics. */
      HybridGridCellListSet<R> 
      upper_evolve(const ConstraintBasedHybridAutomaton<R>& automaton, 
                   const HybridGridCell<R>& initial_set, 
                   const HybridGrid<R>& hybrid_grid, 
                   time_type evolution_time, 
                   size_type maximum_number_of_events) const;
      
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      HybridGridCellListSet<R> 
      upper_reach(const ConstraintBasedHybridAutomaton<R>& automaton, 
                  const HybridGridCell<R>& initial_set, 
                  const HybridGrid<R>& hybrid_grid, 
                  time_type evolution_time, 
                  size_type maximum_number_of_events) const;
      
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. */
      HybridGridMaskSet<R> upper_evolve(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                                  const HybridGridMaskSet<R>& initial_set, 
                                                  time_type evolution_time, 
                                                  size_type maximum_number_of_events) const;
     
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      HybridGridMaskSet<R> upper_reach(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                                 const HybridGridMaskSet<R>& initial_set, 
                                                 time_type evolution_time, 
                                                 size_type maximum_number_of_events) const;

     
      /*! \brief Compute an over-approximation to the set of points which remain in \a bounding_set under evolution of \a automaton. using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      HybridGridMaskSet<R> viable(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                            const HybridGridMaskSet<R>& bounding_set) const;
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. */
      HybridGridMaskSet<R> chainreach(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                                const HybridGridMaskSet<R>& initial_set, 
                                                const HybridGridMaskSet<R>& bounding_set) const;

      /*! \brief Attempt to verify that the reachable set of \a map starting in \a initial_set remains in \a safe_set. */
      tribool verify(const ConstraintBasedHybridAutomaton<R>& automaton, 
                     const HybridGridMaskSet<R>& initial_set, 
                     const HybridGridMaskSet<R>& safe_set) const;
      //@}

      //@{ 
      //! \name Traces and diagnostics
      /*! \brief A trace of the basic sets covered in the last evolution. */
      const std::vector<timed_set_type>& trace() const;
      //@}

     private:
      // Evolve the hybrid automaton within \a domains starting from the initial_set without using discrete transitions (no checking). */
      HybridGridMaskSet<R> _discrete_step(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                                    const HybridGridMaskSet<R>& initial_set,
                                                    const HybridGridMaskSet<R>& domain_set) const;
      // Evolve the hybrid automaton within \a domains starting from the initial_set without using discrete transitions (no checking). */
      HybridGridMaskSet<R> _continuous_chainreach(const ConstraintBasedHybridAutomaton<R>& automaton, 
                                                            const HybridGridMaskSet<R>& initial_set,
                                                            const HybridGridMaskSet<R>& domain_set) const;

     private:

      // Initialisation and finalisation routines
      typedef std::vector<timed_set_type> working_sets_type; 
      working_sets_type _compute_working_sets(const hybrid_list_set_type& set) const;
      hybrid_list_set_type _compute_list_set(const working_sets_type& working_sets, const HybridSpace& locations) const;

     private:
      EvolutionParameters<R>* _parameters;
      ConstraintBasedHybridScheduler<R>* _scheduler;
      mutable std::vector<timed_set_type> _trace;
   };


  
} // namespace Ariadne

#endif /* ARIADNE_CONSTRAINT_BASED_HYBRID_EVOLVER_H */
