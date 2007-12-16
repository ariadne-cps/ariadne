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
  namespace Evaluation {
  
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
      typedef Numeric::Interval<R> I;
      typedef Geometry::Zonotope<I,I> BS;
     public:
      typedef typename System::ConstraintBasedDiscreteMode<R> mode_type;
      typedef typename System::ConstraintBasedDiscreteTransition<R> transition_type;
      typedef Geometry::DiscreteState discrete_state_type;
      typedef typename Geometry::Zonotope<Numeric::Interval<R> > continuous_basic_set_type;
      typedef Geometry::HybridBasicSet<continuous_basic_set_type> hybrid_basic_set_type;
      //typedef Geometry::HybridTimedBasicSet<continuous_basic_set_type> timed_set_type;
      typedef TimeModelHybridBasicSet<continuous_basic_set_type> timed_set_type;
      typedef Geometry::HybridListSet<continuous_basic_set_type> hybrid_list_set_type;
      typedef Geometry::ConstraintInterface<R> constraint_type;
      typedef Geometry::Rectangle<R> bounding_box_type;
      typedef TimeModel<R> time_model_type;

      typedef boost::shared_ptr< const Geometry::ConstraintInterface<R> > constraint_const_pointer;
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
      Geometry::HybridSet<R> discrete_step(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                           const Geometry::HybridSet<R>& initial_set) const;
     
      /*! \brief Evolve the hybrid automaton within \a bounding_set starting from the \a initial_set respecting invariants and without using discrete transitions. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> continuous_chainreach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                   const Geometry::HybridSet<R>& initial_set,
                                                   const Geometry::HybridSet<R>& bounding_set) const;
     

      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> lower_evolve(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                          const Geometry::HybridSet<R>&,
                                          time_type evolution_time,
                                          size_type maximum_number_of_events) const;
      
      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> upper_evolve(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                          const Geometry::HybridSet<R>&,
                                          time_type evolution_time,
                                          size_type maximum_number_of_events) const;
      
      /*! \brief Compute a lower approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> lower_reach(const System::ConstraintBasedHybridAutomaton<R>&, 
                                         const Geometry::HybridSet<R>&, 
                                         time_type initial_evolution_time, 
                                         time_type final_time, 
                                         size_type maximum_number_of_events) const;
      
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> upper_reach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                         const Geometry::HybridSet<R>& initial_set, 
                                         time_type initial_evolution_time, 
                                         time_type final_time, 
                                         size_type maximum_number_of_events) const;
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. */
      Geometry::HybridSet<R> chainreach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                        const Geometry::HybridSet<R>& initial_set, 
                                        const Geometry::HybridSet<R>& bounding_set) const;

      /*! \brief Compute the viability kernel of \a map within \a bounding_set. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> viable(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                    const Geometry::HybridSet<R>& bounding_set) const;
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      tribool verify(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                     const Geometry::HybridSet<R>& initial_set, 
                     const Geometry::HybridSet<R>& safe_set) const;
      //@}

     public:
      //@{
      //! \name Evolution using concrete sets.
      /*! \brief Make a discrete step of the hybrid automaton, starting from initial set. */
      Geometry::HybridGridMaskSet<R> discrete_step(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                       const Geometry::HybridGridMaskSet<R>& initial_set) const;
      /*! \brief Evolve the hybrid automaton within \a bounding_set starting from the \a initial_set respecting invariants and without using discrete transitions. */
      Geometry::HybridGridMaskSet<R> continuous_chainreach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                           const Geometry::HybridGridMaskSet<R>& initial_set,
                                                           const Geometry::HybridGridMaskSet<R>& bounding_set) const;


      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using lower semantics. */
      Geometry::HybridListSet<continuous_basic_set_type> lower_evolve(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                                      const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                                      time_type evolution_time,
                                                                      size_type maximum_number_of_events) const;
      
      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using upper semantics. */
      Geometry::HybridListSet<continuous_basic_set_type> upper_evolve(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                                      const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                                      time_type evolution_time,
                                                                      size_type maximum_number_of_events) const;
      
      /*! \brief Compute the system evolution up to \a time with up to \a maximum_number_of_events using lower semantics. */
      Geometry::HybridListSet<continuous_basic_set_type> lower_reach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                                      const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                                                                      time_type evolution_time,
                                                                      size_type maximum_number_of_events) const;

      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics.*/
      Geometry::HybridListSet<continuous_basic_set_type> 
      upper_reach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                  const Geometry::HybridListSet<continuous_basic_set_type>& initial_set, 
                  time_type evolution_time, 
                  size_type maximum_number_of_events) const;
      

     
     
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
      *  with up to \a maximum_number_of_events using upper semantics. */
      Geometry::HybridGridCellListSet<R> 
      upper_evolve(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                   const Geometry::HybridGridCell<R>& initial_set, 
                   const Geometry::HybridGrid<R>& hybrid_grid, 
                   time_type evolution_time, 
                   size_type maximum_number_of_events) const;
      
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridGridCellListSet<R> 
      upper_reach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                  const Geometry::HybridGridCell<R>& initial_set, 
                  const Geometry::HybridGrid<R>& hybrid_grid, 
                  time_type evolution_time, 
                  size_type maximum_number_of_events) const;
      
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. */
      Geometry::HybridGridMaskSet<R> upper_evolve(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                  const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                  time_type evolution_time, 
                                                  size_type maximum_number_of_events) const;
     
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridGridMaskSet<R> upper_reach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                 const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                 time_type evolution_time, 
                                                 size_type maximum_number_of_events) const;

     
      /*! \brief Compute an over-approximation to the set of points which remain in \a bounding_set under evolution of \a automaton. using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridGridMaskSet<R> viable(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                            const Geometry::HybridGridMaskSet<R>& bounding_set) const;
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. */
      Geometry::HybridGridMaskSet<R> chainreach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                const Geometry::HybridGridMaskSet<R>& bounding_set) const;

      /*! \brief Attempt to verify that the reachable set of \a map starting in \a initial_set remains in \a safe_set. */
      tribool verify(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                     const Geometry::HybridGridMaskSet<R>& initial_set, 
                     const Geometry::HybridGridMaskSet<R>& safe_set) const;
      //@}

      //@{ 
      //! \name Traces and diagnostics
      /*! \brief A trace of the basic sets covered in the last evolution. */
      const std::vector<timed_set_type>& trace() const;
      //@}

     private:
      // Evolve the hybrid automaton within \a domains starting from the initial_set without using discrete transitions (no checking). */
      Geometry::HybridGridMaskSet<R> _discrete_step(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridGridMaskSet<R>& initial_set,
                                                    const Geometry::HybridGridMaskSet<R>& domain_set) const;
      // Evolve the hybrid automaton within \a domains starting from the initial_set without using discrete transitions (no checking). */
      Geometry::HybridGridMaskSet<R> _continuous_chainreach(const System::ConstraintBasedHybridAutomaton<R>& automaton, 
                                                            const Geometry::HybridGridMaskSet<R>& initial_set,
                                                            const Geometry::HybridGridMaskSet<R>& domain_set) const;

     private:

      // Initialisation and finalisation routines
      typedef std::vector<timed_set_type> working_sets_type; 
      working_sets_type _compute_working_sets(const hybrid_list_set_type& set) const;
      hybrid_list_set_type _compute_list_set(const working_sets_type& working_sets, const Geometry::HybridSpace& locations) const;

     private:
      EvolutionParameters<R>* _parameters;
      ConstraintBasedHybridScheduler<R>* _scheduler;
      mutable std::vector<timed_set_type> _trace;
   };


  }
}

#endif /* ARIADNE_CONSTRAINT_BASED_HYBRID_EVOLVER_H */
