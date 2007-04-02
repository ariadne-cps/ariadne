/***************************************************************************
 *            hybrid_evolver.h
 *
 *  Copyright  2006  Alberto Casagrande,  Pieter Collins
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
 
#ifndef ARIADNE_HYBRID_EVOLVER_H
#define ARIADNE_HYBRID_EVOLVER_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "../base/types.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"
#include "../evaluation/declarations.h"

namespace Ariadne {  
  namespace Evaluation {
  

    /*! \ingroup Evolve
     *  \brief A class for computing the evolution of a hybrid system.
     */
    template< class R >
    class HybridEvolver
    {
     public:
      /*! \brief Construct from an applicator and an integrator. */
     HybridEvolver(Applicator<R>& applicator, Integrator<R>& integrator);
      
      Geometry::HybridGridMaskSet<R> evolve(const System::HybridAutomaton<R>& automaton, 
                                            const Geometry::HybridGridMaskSet<R>& initial_set, 
                                            time_type evolution_time,
                                            size_type maximum_number_of_events);
      
      /*! \brief Compute the system evolution at \a time with up to \a maximum_number_of_events. */
      Geometry::HybridGridMaskSet<R> upper_evolve(const System::HybridAutomaton<R>& automaton, 
                                            const Geometry::HybridGridMaskSet<R>& initial_set, 
                                            time_type evolution_time,
                                            size_type maximum_number_of_events);
      
      /*! \brief Compute a lower approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using lower semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridGridMaskSet<R> lower_reach(const System::HybridAutomaton<R>&, 
                                                 const Geometry::HybridGridMaskSet<R>&, 
                                                 time_type initial_evolution_time, 
                                                 time_type final_time, 
                                                 size_type maximum_number_of_events);
      
      /*! \brief Compute an over approximation to the reachable set between \a initial_evolution_time and \a final_time
       *  with up to \a maximum_number_of_events using upper semantics. (NOT CURRENTLY IMPLEMENTED) */
      Geometry::HybridGridMaskSet<R> upper_reach(const System::HybridAutomaton<R>& automaton, 
                                                 const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                 time_type initial_evolution_time, 
                                                 time_type final_time, 
                                                 size_type maximum_number_of_events);
     
      /*! \brief Compute an over approximation to the chain-reachable set using upper semantics. */
      Geometry::HybridGridMaskSet<R> chainreach(const System::HybridAutomaton<R>& automaton, 
                                                const Geometry::HybridGridMaskSet<R>& initial_set, 
                                                const Geometry::HybridGridMaskSet<R>& bounding_set);
     private:
     public:
      /*! \brief Make a discrete step of the hybrid automaton, starting from initial set. */
      Geometry::HybridGridCellListSet<R> discrete_step(const System::HybridAutomaton<R>& automaton, 
                                                       const Geometry::HybridGridCellListSet<R>& initial_set);
      /*! \brief Evolve the hybrid automaton withing \a invariants starting from the \a initial_set without using discrete transitions. */
      Geometry::HybridGridMaskSet<R> continuous_chainreach(const System::HybridAutomaton<R>& automaton, 
                                                           const Geometry::HybridGridMaskSet<R>& initial_set,
                                                           const Geometry::HybridGridMaskSet<R>& invariants);
     private:
      Applicator<R>* _applicator;
      Integrator<R>* _integrator;
    };


  }
}

#endif /* ARIADNE_HYBRID_EVOLVER_H */
