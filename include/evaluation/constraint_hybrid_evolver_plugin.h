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
#include <set>
#include <map>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "../numeric/float64.h"

#include "../base/types.h"
#include "../base/tribool.h"
#include "../base/reference_container.h"
#include "../geometry/declarations.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"
#include "../evaluation/declarations.h"
#include "../evaluation/integrator.h"

#include "../evaluation/hybrid_time.h"

namespace Ariadne {  
  namespace Evaluation {
  
    template<class R> class ConstraintHybridEvolver;
    template<class R> class LohnerIntegrator;
  
    /*! \brief The set appears to cross a constraint non-transversely. */
    class NonTransverseCrossingException : public std::exception { };
    /*! \brief The set appears to cross two or more constraints. */
    class CornerCollisionException : public std::exception { };
    /*! \brief A constraint is crossed during a time interval which does not contain the integration time window. */
    class PartiallyEnabledConstraintException : public std::exception { };

    /*! \brief The semantics used for the evolution trajectories (\c lower_semantics or \c upper_semantics). */
    enum EvolutionSemantics { lower_semantics, upper_semantics };
    /*! \brief The kind of set to approximate (\c compute_evolved_set or \c compute_reachable_set). */
    enum EvolutionKind { compute_evolved_set, compute_reachable_set };

  
    /*! \brief Data to hold information about a constraint crossing. */
    template<class R> 
    struct CrossingData { 
      Numeric::Interval<R> initial_constraint_value; 
      Numeric::Interval<R> final_constraint_value; 
      Numeric::Interval<R> normal_derivative;
      Numeric::Interval<R> crossing_time_bounds;
      TimeModel<R> crossing_time_model; 
      int direction;
    };
  

    /*! \brief Data to hold information about a constraint crossing. */
    template<class R> 
    struct IntegrationData { 
      typedef Numeric::Rational Time;
      Time final_time;
      Time maximum_time_step;
      TimeModel<R> spacial_integration_time_step_model;
      TimeModel<R> parameterised_integration_time_step_model;
      tribool terminating_event;
    };


    template<class R> 
    std::ostream& operator<<(std::ostream& os, const CrossingData<R>& cd);

    /*!\ingroup Evolve
     * \brief A class for computing the evolution of a hybrid system.
     *
     * The evolution of a hybrid system consists of continuous evolution interspersed
     * by discrete events. These events may be \em forced, which means they occur 
     * immediately that the \em guard constraint is satisfied, or \em unforced,
     * which means that they may occur anytime that a <em>activation</em> is to
     * be satisfied. Additionally, the continuous evolution may be \em blocked by an
     * \em invariant which fails to be satisfied.
     *
     * The semantics of activations is particularly challenging, since the evolved set
     * \f$\Psi(x,t)\f$ also depends on the event time \f$t_1\f$ for \f$t>t_1\f$.
     * In %Ariadne, we compute the set of points reached under <em>all</em> possible 
     * event times. 
     * (In a simulation tool, a single event time would probably be chosen,
     * possibly depending on some event probability distribution.)
     * This means that the number of generators of the evolved set is increased after 
     * such a transition.
     * 
     * The strategy used in %Ariadne for computing the evolution is as follows:
     * -# Estimate a time step size and a bounding box for the continuous evolution.
     * -# Compute the event which terminates the continuous evolution. 
     *    (This may be the "special event" that the time step size is reached.) \n
     *    Compute the time at which the evolution is terminated.
     * -# Compute the non-forced transitions which may occur during the continuous evolution. \n
     *    For each possible non-forced event, compute the time interval for which the event is possible.
     * -# For each event (including the terminating event), insert the the set of points reached immediately 
     *    after the event into the new working sets.
     * -# In the case of a reachability step, return the set of points reached by the continuous evolution.
     * .
     * We call the time at the end of the continuous evolution the \em terminal time.
     *
     * The above algorithm can be used directly for point-based approximate evolution routines. 
     * However, for set-based routines, there are a number of difficulties which need to be addresses:
     * -# The time of a forced event, or the time interval of an unforced event, depends on the initial point in the set.\n
     * -# The initial set may partially satisfy a guard condition.
     * -# The set at the end of an integration step may partially satisfy a guard condition.
     * -# The set at the maximum evolution time may partially satisfy a guard condition.
     * -# The set at the maximum evolution time may partially satisfy a guard condition.
     * -# Different forced transitions may terminate the evolution for different initial conditions.
     * -# An unforced transition may be partially activated at the initial/terminal time.
     * -# The flow may not be transverse to a constraint set, causing 
     *     -# Some points in the initial set to jump, and some points not to.
     *     -# A singularity in the jump map.
     *
     * There are two strategies for dealing with problem 1, namely the \em timed set strategy and the \em saltation map strategy.
     * - In the timed set strategy, the main data structure is records a set of space-time pairs \f$(y,t)\f$ which the flow can reach.
     *    This strategy has the advantage of storing the exact reachable set after a transision, 
     *    but the disadvantage of additional complexity in the switching logig.
     * - In the saltation map strategy, after every event, the flow is evolved forwards and backwards to the time of the centre point. 
     *    This has the advantage of keeping the time at each step fixed, but the disadvantage of requiring backwards integration. 
     *    As a result, a point may flow backwards into a guard set causing spurious jumps, 
     *    and also meaning that the set cannot be subdivided when considering forwards evolution, since points flowing backwards may be lost.
     *    
     *
     * The main difficult is event detection. As well as ordinary transitions,
     * we have the following special events 
     * -#  maximum_integration_time_reached
     * -#  integration_time_step_reached
     * 
     * The maximum integration time step \a maximum_time_step
     * is a constant determined by the time-stepping algorithm. 
     * Since for times higher than this value, the \a bounding_box
     * of the flow may not be reached, this time cannot be violated,
     * regardless of the semantics.
     *
     * The maximum evolution time \a maximum_time is an input variable
     * This time cannot be violated by lower semantics.
     * It may be violated with upper semantics by the action of a discrete event, 
     * but only if the set of initial point/time pairs for which the violation
     * occurs also give rise to a continuous evolution ending at \a maximum_time.
     *
     * The main evolution algorithm is performed by the following pseudocode:
     *
     * <c>\pseudocode
     *     Time maximum_time
     *     TimeModel maximum_time_step - maximum_time - initial_set.time()
     *
     *     // Compute all possible blocking events, including the special event
     *     map<Event,TimeModel> blocking_events = compute_blocking_events_and_time_steps()
     *
     *     if upper_semantics and size(blocking_events) > 1 then
     *         if radius>splitting_size then
     *             return subdivide(initial_set)
     *         else
     *     end
     *  
     *     // Some events may occur at times greater than the maximum integration time step
     *     //   
     *     foreach (event,time) in block do
     *     end 
     * \endpseudocode 
     */
    template< class R >
    class ConstraintHybridEvolverPlugin
    {
      friend class ConstraintHybridEvolver<R>;
      typedef Numeric::Interval<R> I;
     public:
      /*! \brief The type used for real numbers. */
      typedef R real_type;
      /*! \brief . */
      typedef Numeric::Rational time_type;
      /*! \brief . */
      typedef typename System::ConstraintDiscreteMode<R> mode_type;
      /*! \brief . */
      typedef typename System::ConstraintDiscreteTransition<R> transition_type;
      /*! \brief . */
      typedef typename Geometry::Zonotope<Numeric::Interval<R> > continuous_basic_set_type;
      /*! \brief . */
      typedef Geometry::HybridBasicSet<continuous_basic_set_type> hybrid_basic_set_type;
      /*! \brief . */
      typedef TimeModelHybridBasicSet<continuous_basic_set_type> timed_set_type;
      /*! \brief . */
      typedef Geometry::HybridListSet<continuous_basic_set_type> hybrid_list_set_type;
      /*! \brief . */
      typedef Geometry::ConstraintInterface<R> constraint_type;
      /*! \brief . */
      typedef Geometry::Rectangle<R> bounding_box_type;
      /*! \brief . */
      typedef System::MapInterface<R> map_type;
      /*! \brief . */
      typedef System::VectorFieldInterface<R> vector_field_type;
      /*! \brief . */
      typedef Evaluation::TimeModel<R> time_model_type;
      /*! \brief . */
      typedef std::pair<Evaluation::TimeModel<R>,Evaluation::TimeModel<R> > time_model_pair_type;
      /*! \brief . */
      typedef CrossingData<R> crossing_data_type;
  
      //! \name Constructors. */

      /*! \brief Construct from an applicator and an integrator. */
      ConstraintHybridEvolverPlugin(Applicator<R>& applicator, Integrator<R>& integrator);

      /*! \brief Copy constructor. */
      ConstraintHybridEvolverPlugin(const ConstraintHybridEvolverPlugin<R>& plugin);

      //! \name Parameters controlling the evolution.
      /*! \brief The maximum step size for integration. */
      time_type maximum_step_size() const;
      /*! \brief The maximum size of a basic set. */
      real_type maximum_basic_set_radius() const;
      /*! \brief The maximum size of a basic set which can be split across a transition boundary. */
      real_type maximum_splitting_set_radius() const;


      //! \name Evolution steps

      /*! \brief Compute the an approximation to the evolution of a timed basic set using given semantics. 
       * For reachability analysis, the last set in the list is the set reached by continuous evolution; 
       * the remaining sets should be continued as working sets.
       */
      std::vector<timed_set_type>
      evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                     const timed_set_type& initial_set,
                     const time_type& maximum_time,
                     EvolutionSemantics evolution_semantics,
                     EvolutionKind evolution_kind) const;

      /*! \brief Compute a lower-approximation to the evolution of a timed basic set using lower semantics. */
      std::vector<timed_set_type>
      lower_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                           const timed_set_type& initial_set,
                           const time_type& maximum_time) const;

      /*! \brief Compute the an over-approximation to the evolution of a timed basic set using upper semantics. */
      std::vector<timed_set_type>
      upper_evolution_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                           const timed_set_type& initial_set,
                           const time_type& maximum_time) const;

      /*! \brief Compute the possible states reached during an evolution step using lower semantics. */
      std::vector<timed_set_type>
      lower_reachability_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                              const timed_set_type& initial_set,
                              const time_type& maximum_time) const;

      /*! \brief Compute the possible states reached during an evolution step using upper semantics. */
      std::vector<timed_set_type>
      upper_reachability_step(const System::ConstraintHybridAutomaton<R>& automaton, 
                              const timed_set_type& initial_set,
                              const time_type& maximum_time) const;


      //! \name Single event evolution steps

      /*! \brief Compute the possible states reached by a forced jump occurring after time \a crossing_time, given that the flow remains in \a bounding_box. */
      timed_set_type
      forced_jump_step(const transition_type& transition,
                       const timed_set_type& initial_set,
                       const time_model_type& crossing_time_step,
                       const bounding_box_type& bounding_box) const;


      /*! \brief Compute the possible states reached by an unforced jump occurring between time \a lower_time_step and \a upper_time_step, given that the flow remains in \a bounding_box. */
      timed_set_type
      unforced_jump_step(const transition_type& transition,
                         const timed_set_type& initial_set,
                         const time_model_type& lower_time_step,
                         const time_model_type& upper_time_step,
                         const bounding_box_type& bounding_box) const;
      
      /*! \brief Compute the possible states reached by an instantaneous transition. */
      timed_set_type
      discrete_event_step(const transition_type& transition,
                          const timed_set_type& initial_set) const;


      /*! \brief Compute the possible states reached by an integration step with time given by time_step, given that the flow remains in \a bounding_box. */
      timed_set_type
      continuous_evolution_step(const mode_type& mode,
                                const timed_set_type& initial_set,
                                const time_model_type& time_step,
                                const bounding_box_type& bounding_box) const;

      /*! \brief Compute the possible states reached by a continuous evolution with times between \a lower_time_step and \a upper_time_step, given that the flow remains in \a bounding_box. */
      timed_set_type
      continuous_reachability_step(const mode_type& mode,
                                   const timed_set_type& initial_set,
                                   const time_model_type& lower_time_step,
                                   const time_model_type& upper_time_step,
                                   const bounding_box_type& bounding_box) const;


      /*! \brief Compute the possible states reached by continuous evolution at the time \a final_time, given that the flow remains in \a bounding_box. */
      timed_set_type
      final_continuous_evolution_step(const mode_type& mode,
                                      const timed_set_type& initial_set,
                                      const time_model_type& final_time,
                                      const bounding_box_type& bounding_box) const;

      /*! \brief Compute the possible states reached by continuous evolution up to \a final_time, given that the flow remains in \a bounding_box. */
      timed_set_type
      final_continuous_reachability_step(const mode_type& mode,
                                         const timed_set_type& initial_set,
                                         const time_model_type& final_time,
                                         const bounding_box_type& bounding_box) const;

      /*! \brief Compute an accurate approximation to the crossing time step. */
      time_model_type
      compute_crossing_time_step(const transition_type& transition,
                                 const timed_set_type& initial_set,
                                 const bounding_box_type& bounding_box) const;
                             

      //! \name Event sequenceing methods

      /*! Tests if one constraint forces another */
      tribool
      forces(const constraint_type& transition1,
             const constraint_type& transition2,
             const bounding_box_type& bounding_box) const;

      /*! \brief Test if an event is activated. */
      tribool
      enabled(const transition_type& transition,
              const bounding_box_type& bounding_box) const;
                             
      /*! \brief Test if an event is activated. */
      tribool
      enabled(const transition_type& transition,
              const timed_set_type& set) const;
                             
      /*! \brief Estimate the direction of a flow transverse to a constraint. */
      I
      estimate_normal_derivative(const transition_type& transition,
                                 const bounding_box_type& bounding_box) const;
                             
      /*! \brief Compute an estimation of the time step of a crossing with a constraint ifor a transition. */
      I
      estimate_crossing_time_step(const transition_type& transition,
                                  const timed_set_type& initial_set,
                                  const bounding_box_type& bounding_box) const;
                             
      /*! \brief Compute the crossings with a constraint for a transition. */
      crossing_data_type
      compute_crossing_data(const transition_type& transition,
                            const timed_set_type& initial_set,
                            const timed_set_type& final_set,
                            const bounding_box_type& bounding_box) const;

                             
      /*! \brief Compute the crossings with the guard set. */
      std::map<id_type, crossing_data_type>
      compute_crossing_data(const reference_set<const transition_type>& transitions,
                            const timed_set_type& initial_set,
                            const timed_set_type& final_set,
                            const bounding_box_type& bounding_box,
                            const time_type& maximum_time_step) const;


      /*! \brief Compute the crossings with the guard set. */
      std::map<id_type, time_model_type>
      compute_terminating_event_times(const mode_type& mode,
                                      const timed_set_type& initial_set,
                                      const timed_set_type& final_set,
                                      const time_type& maximum_time,
                                      const bounding_box_type& bounding_box) const;

      /*! \brief Compute the activations which are enabled. */
      std::map<id_type, time_model_pair_type>
      compute_enabled_activation_times(const mode_type& mode,
                                       const timed_set_type& initial_set,
                                       const time_model_type& maximum_time_step,
                                       const bounding_box_type& bounding_box) const;
 

      I compute_evolution_time_bounds(const std::map<id_type, time_model_type>& event_times) const;
      I compute_evolution_time_bounds(const std::map<id_type, crossing_data_type>& crossing_data) const;

      reference_set<transition_type>
      compute_possibly_enabled_transitions(const mode_type& mode,
                                           const bounding_box_type& bounding_box) const;
                           
      reference_set<transition_type>
      compute_terminating_time_step(const reference_key_map<transition_type, time_model_type>& terminating_event_times,
                                    const bounding_box_type& bounding_box) const;
                           
                           
                 
 

      //! \name Utility functions

      /*! \brief Compute the crossing times for a set of modes. */
      Geometry::Rectangle<R> 
      estimate_flow_bounds(const mode_type& mode, 
                           const timed_set_type& initial_set,
                           time_type& maximum_step_size) const;

       /*! \brief Compute the crossing times for a set of modes. */
      Geometry::Rectangle<R> 
      refine_flow_bounds(const mode_type& mode, 
                         const timed_set_type& initial_set,
                         const bounding_box_type& bounding_set,
                         const time_type& maximum_step_size) const;
     
      /*! \brief Subdivide the set. */
      std::vector<timed_set_type> 
      subdivide(const timed_set_type& ts) const;
       
      /*! \brief Regularise the set if the set is close to being singular. */
      timed_set_type 
      regularize(const timed_set_type& ts) const;
       
      //@}
     public:
      mutable std::vector<timed_set_type> trace;
     private:
      Applicator<R>* _applicator;
      LohnerIntegrator<R>* _integrator;
    };



    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const CrossingData<R>& cd)
    {
      return os << "{ initial_constraint_value=" << cd.initial_constraint_value
                << ", final_constraint_value=" << cd.final_constraint_value 
                << ", normal_derivative=" << cd.normal_derivative
                << ", crossing_time_bounds=" << cd.crossing_time_bounds
                << ", crossing_time_model=" << cd.crossing_time_model 
                << ", direction=" << cd.direction
                << " }";
    }


  }
}

#endif /* ARIADNE_CONSTRAINT_HYBRID_EVOLVER_PLUGIN_H */
