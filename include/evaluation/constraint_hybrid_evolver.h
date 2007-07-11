/***************************************************************************
 *            constraint_hybrid_evolver.h
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
 
#ifndef ARIADNE_CONSTRAINT_HYBRID_EVOLVER_H
#define ARIADNE_CONSTRAINT_HYBRID_EVOLVER_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "../base/types.h"
#include "../base/tribool.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"
#include "../evaluation/declarations.h"

namespace Ariadne {  
  namespace Evaluation {
  
    template<class R> class LohnerIntegrator;


    /*! \ingroup Evolve
     *  \brief A class for computing the evolution of a hybrid system.
     */
    template< class R >
    class ConstraintHybridEvolver
    {
      typedef typename System::ConstraintDiscreteMode<R> mode_type;
      typedef typename System::ConstraintDiscreteTransition<R> transition_type;
      typedef typename Geometry::Zonotope<Numeric::Interval<R> > basic_set_type;
      typedef Geometry::HybridBasicSet<basic_set_type> hybrid_basic_set_type;
      typedef Geometry::TimedSet<hybrid_basic_set_type> timed_set_type;
      typedef boost::shared_ptr< const Geometry::ConstraintInterface<R> > constraint_pointer;
     public:

      /*! \brief Construct from an applicator and an integrator. */
      ConstraintHybridEvolver(Applicator<R>& applicator, Integrator<R>& integrator);

      /*! \brief Virtual destructor. */
      virtual ~ConstraintHybridEvolver();


      //@{
      //! \name Evolution steps.
      basic_set_type
      unforced_jump(const System::VectorFieldInterface<R>& dynamic1,
                    const System::VectorFieldInterface<R>& dynamic2,
                    const System::MapInterface<R>& reset,
                    const basic_set_type& basic_set,
                    const Geometry::ConstraintInterface<R>& activation,
                    const time_type& step_size);

      basic_set_type
      forced_jump(const System::VectorFieldInterface<R>& dynamic1,
                  const System::VectorFieldInterface<R>& dynamic2,
                  const System::MapInterface<R>& reset,
                  const basic_set_type& basic_set,
                  const Geometry::ConstraintInterface<R>& guard,
                  const time_type& step_size);

      std::vector<timed_set_type>
      lower_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                           const timed_set_type& initial_set,
                           time_type& step_size);

      std::vector<timed_set_type>
      upper_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                           const timed_set_type& initial_set,
                           time_type& step_size);
      //@}


      //@{
      //! \name Evolution using abstract sets.
      /*! \brief Make a discrete step of the hybrid automaton, starting from initial set. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> discrete_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                           const Geometry::HybridSet<R>& initial_set);
     
      /*! \brief Evolve the hybrid automaton within \a bounding_set starting from the \a initial_set respecting invariants and without using discrete transitions. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> continuous_chainreach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                   const Geometry::HybridSet<R>& initial_set,
                                                   const Geometry::HybridSet<R>& bounding_set);
     

      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      std::vector<hybrid_basic_set_type> lower_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                  const Geometry::HybridSet<R>&,
                                                  time_type evolution_time,
                                                  size_type maximum_number_of_events);
      
      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      std::vector<hybrid_basic_set_type> upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                  const Geometry::HybridSet<R>&,
                                                  time_type evolution_time,
                                                  size_type maximum_number_of_events);
      
      /*! \brief Compute a lower approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> lower_reach(const System::ConstraintHybridAutomaton<R>&, 
                                         const Geometry::HybridSet<R>&, 
                                         time_type initial_evolution_time, 
                                         time_type final_time, 
                                         size_type maximum_number_of_events);
      
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> upper_reach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                         const Geometry::HybridSet<R>& initial_set, 
                                         time_type initial_evolution_time, 
                                         time_type final_time, 
                                         size_type maximum_number_of_events);
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. */
      Geometry::HybridSet<R> chainreach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                        const Geometry::HybridSet<R>& initial_set, 
                                        const Geometry::HybridSet<R>& bounding_set);

      /*! \brief Compute the viability kernel of \a map within \a bounding_set. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridSet<R> viable(const System::ConstraintHybridAutomaton<R>& automaton, 
                                         const Geometry::HybridSet<R>& bounding_set);
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      tribool verify(const System::ConstraintHybridAutomaton<R>& automaton, 
                     const Geometry::HybridSet<R>& initial_set, 
                     const Geometry::HybridSet<R>& safe_set);
      //@}

     public:
      //@{
      //! \name Evolution using concrete sets.
      /*! \brief Make a discrete step of the hybrid automaton, starting from initial set. */
      Geometry::HybridGridMaskSet<R> discrete_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                       const Geometry::HybridGridMaskSet<R>& initial_set);
      /*! \brief Evolve the hybrid automaton within \a bounding_set starting from the \a initial_set respecting invariants and without using discrete transitions. */
      Geometry::HybridGridMaskSet<R> continuous_chainreach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                           const Geometry::HybridGridMaskSet<R>& initial_set,
                                                           const Geometry::HybridGridMaskSet<R>& bounding_set);

      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using lower semantics. */
      Geometry::HybridListSet<basic_set_type> lower_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                 const hybrid_basic_set_type& initial_set, 
                                                 time_type evolution_time,
                                                 size_type maximum_number_of_events);
      
      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using upper semantics. */
      Geometry::HybridListSet<basic_set_type> upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                 const hybrid_basic_set_type& initial_set, 
                                                 time_type evolution_time,
                                                 size_type maximum_number_of_events);
      
      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using lower semantics. */
      Geometry::HybridGridMaskSet<R> lower_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                  const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                  time_type evolution_time,
                                                  size_type maximum_number_of_events);
      
      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events using upper semantics. */
      Geometry::HybridGridMaskSet<R> upper_evolve(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                  const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                  time_type evolution_time,
                                                  size_type maximum_number_of_events);
      
      /*! \brief Compute a lower approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridGridMaskSet<R> lower_reach(const System::ConstraintHybridAutomaton<R>&, 
                                                 const Geometry::HybridGridMaskSet<R>&, 
                                                 time_type initial_evolution_time, 
                                                 time_type final_time, 
                                                 size_type maximum_number_of_events);
      
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridGridMaskSet<R> upper_reach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                 const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                 time_type initial_evolution_time, 
                                                 time_type final_time, 
                                                 size_type maximum_number_of_events);
     
      /*! \brief Compute an over-approximation to the set of points which remain in \a bounding_set under evolution of \a automaton. using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridGridMaskSet<R> viable(const System::ConstraintHybridAutomaton<R>& automaton, 
                                            const Geometry::HybridGridMaskSet<R>& bounding_set);
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. */
      Geometry::HybridGridMaskSet<R> chainreach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                const Geometry::HybridGridMaskSet<R>& bounding_set);

      /*! \brief Attempt to verify that the reachable set of \a map starting in \a initial_set remains in \a safe_set. */
      tribool verify(const System::ConstraintHybridAutomaton<R>& automaton, 
                     const Geometry::HybridGridMaskSet<R>& initial_set, 
                     const Geometry::HybridGridMaskSet<R>& safe_set);
      //@}
     private:
      // Evolve the hybrid automaton within \a domains starting from the initial_set without using discrete transitions (no checking). */
      Geometry::HybridGridMaskSet<R> _discrete_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                    const Geometry::HybridGridMaskSet<R>& initial_set,
                                                    const Geometry::HybridGridMaskSet<R>& domain_set);
      // Evolve the hybrid automaton within \a domains starting from the initial_set without using discrete transitions (no checking). */
      Geometry::HybridGridMaskSet<R> _continuous_chainreach(const System::ConstraintHybridAutomaton<R>& automaton, 
                                                            const Geometry::HybridGridMaskSet<R>& initial_set,
                                                            const Geometry::HybridGridMaskSet<R>& domain_set);
     private:
      Applicator<R>* _applicator;
      // FIXME: Allow arbitrary integrators
      LohnerIntegrator<R>* _integrator;
    };


  }
}

#endif /* ARIADNE_HYBRID_EVOLVER_H */
