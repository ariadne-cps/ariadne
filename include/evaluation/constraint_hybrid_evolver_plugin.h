/***************************************************************************
 *            constraint_hybrid_evolver_plugin.h
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
 
#ifndef ARIADNE_CONSTRAINT_HYBRID_EVOLVER_PLUGIN_H
#define ARIADNE_CONSTRAINT_HYBRID_EVOLVER_PLUGIN_H

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
#include "../evaluation/integrator.h"

#include "../evaluation/hybrid_time.h"

namespace Ariadne {  
  namespace Evaluation {
  
    template<class R> class ConstraintHybridEvolver;
    template<class R> class LohnerIntegrator;
  
    /*! \ingroup Evolve
     *  \brief A class for computing the evolution of a hybrid system.
     */
    template< class R >
    class ConstraintHybridEvolverPlugin
    {
      friend class ConstraintHybridEvolver<R>;
      typedef Numeric::Interval<R> I;
     public:
      typedef typename System::ConstraintDiscreteMode<R> mode_type;
      typedef typename System::ConstraintDiscreteTransition<R> transition_type;
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

      /*! \brief Construct from an applicator and an integrator. */
      ConstraintHybridEvolverPlugin(Applicator<R>& applicator, Integrator<R>& integrator);

      /*! \brief Copy constructor. */
      ConstraintHybridEvolverPlugin(const ConstraintHybridEvolverPlugin<R>& plugin);

      //@{
      //! \name Parameters controlling the evolution.
      /*! \brief The maximum step size for integration. */
      time_type maximum_step_size() const;
      /*! \brief The maximum step size for integration. */
      real_type maximum_basic_set_radius() const;

      //@}

      //! \name Evolution steps.

      /*! \brief Compute the possible states reached by a forced jump, given that the flow remains in \a bounding_box.
       */
      std::vector<timed_set_type>
      forced_jump_step(const transition_type& transition,
                       const timed_set_type& initial_set,
                       const bounding_box_type& bounding_box) const;


      /*! \brief Compute the possible states reached by a forced jump, given that the flow remains in \a bounding_box.
       */
      std::vector<timed_set_type>
      unforced_jump_step(const transition_type& transition,
                         const timed_set_type& initial_set,
                         const bounding_box_type& bounding_box,
                         const time_type& maximum_time) const;

      timed_set_type
      discrete_step(const transition_type& transition,
                    const timed_set_type& initial_set) const;

      timed_set_type
      integration_step(const mode_type& mode,
                       const timed_set_type& initial_set,
                       const bounding_box_type& bounding_box,
                       const time_model_type& time_step) const;

      timed_set_type
      integration_step(const mode_type& mode,
                       const timed_set_type& initial_set,
                       const bounding_box_type& bounding_box,
                       const time_type& final_time) const;

      std::vector<timed_set_type>
      lower_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                           const timed_set_type& initial_set,
                           time_type& maximum_time) const;

      std::vector<timed_set_type>
      upper_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                           const timed_set_type& initial_set,
                           time_type& maximum_time) const;

      std::vector<timed_set_type>
      lower_reachability_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                              const timed_set_type& initial_set,
                              time_type& maximum_time) const;

      std::vector<timed_set_type>
      upper_reachability_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                              const timed_set_type& initial_set,
                              time_type& maximum_time) const;


      /*! \brief Subdivide the set. */
      std::vector<timed_set_type> subdivide(const timed_set_type& ts) const;
      //@}
     private:
      Applicator<R>* _applicator;
      LohnerIntegrator<R>* _integrator;
    };


  }
}

#endif /* ARIADNE_CONSTRAINT_HYBRID_EVOLVER_PLUGIN_H */
