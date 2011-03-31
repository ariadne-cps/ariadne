/***************************************************************************
 *            hybrid_evolver-image.h
 *
 *  Copyright  2007-10  Alberto Casagrande, Pieter Collins, Luca Geretti
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

/*! \file hybrid_evolver-image.h
 *  \brief Evolver for hybrid systems using image sets.
 */

#ifndef ARIADNE_HYBRID_EVOLVER_IMAGE_H
#define ARIADNE_HYBRID_EVOLVER_IMAGE_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "tuple.h"

#include "hybrid_set.h"

#include "hybrid_automaton.h"
#include "evolver_interface.h"
#include "evolver_base.h"
#include "settings.h"

#include "logging.h"

namespace Ariadne {

template<class Sys, class BS> class Evolver;

class VectorTaylorFunction;
class TaylorSet;
class HybridAutomaton;
template<class ES> class Orbit;

class EvolutionSettings;
class TaylorModel;
template<class MDL> class CalculusInterface;

class HybridTime;

class DiscreteEvent;

/*! \brief A class for computing the evolution of a hybrid system.
 *
 * The actual evolution steps are performed by the HybridEvolver class.
 */
class ImageSetHybridEvolver
    : public EvolverBase<HybridAutomaton,HybridTaylorSet>
    , public Loggable
{
    typedef ScalarFunction ScalarFunctionType;
    typedef VectorFunction VectorFunctionType;
    typedef Vector<Interval> BoxType;
    typedef VectorTaylorFunction FunctionModelType;
    typedef VectorTaylorFunction MapModelType;
    typedef VectorTaylorFunction FlowModelType;
    typedef ScalarTaylorFunction ConstraintModelType;
    typedef TaylorModel TimeModelType;
    typedef TaylorSet FlowSetModelType;
    typedef TaylorSet SetModelType;
  public:
    typedef ContinuousEvolutionSettings EvolutionSettingsType;
    typedef HybridAutomaton::TimeType TimeType;
    typedef int IntegerType;
    typedef Float RealType;
    typedef std::vector<DiscreteEvent> EventListType;
    typedef HybridAutomaton SystemType;
    typedef TaylorSet ContinuousEnclosureType;
    typedef TaylorSet TimedSetModelType;
    typedef HybridBasicSet<TaylorSet> HybridEnclosureType;
    typedef HybridEnclosureType EnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<EnclosureType> EnclosureListType;
    typedef Float ContinuousTimeType;
    typedef tuple<DiscreteState, EventListType, SetModelType, TimeModelType> HybridTimedSetType;
    typedef std::map< DiscreteEvent,tuple<TaylorModel,TaylorModel> > ActivationTimesType;
  public:

    //! \brief Default constructor.
    ImageSetHybridEvolver();

    //! \brief Construct from parameters using a default integrator.
    ImageSetHybridEvolver(const EvolutionSettingsType& parameters);

    //! \brief Construct from integrator with default parameters.
    ImageSetHybridEvolver(const TaylorCalculus& tc);

    //! \brief Construct from parameters and integrator.
    ImageSetHybridEvolver(const EvolutionSettingsType& p,
   						  const TaylorCalculus& tc);

    /*! \brief Make a dynamically-allocated copy. */
    ImageSetHybridEvolver* clone() const { return new ImageSetHybridEvolver(*this); }

    const CalculusInterface<TaylorModel>& getCalculusInterface() const { return *this->_toolbox; }

    //@{
    //! \name Settings controlling the evolution.

    //! \brief A reference to the settings controlling the evolution.
    EvolutionSettingsType& settings() { return *this->_settings; }
	//! \brief A constant reference to the settings controlling the evolution.
    const EvolutionSettingsType& settings() const { return *this->_settings; }

    //@}

    //@{
    //! \name Evolution using abstract sets.

    //! \brief Compute an approximation to the orbit set using the given semantics.
    Orbit<EnclosureType> orbit(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const;

    //! \brief Compute an approximation to the orbit set for upper semantics, with continuous evolution only.
    Orbit<EnclosureType> upper_orbit_continuous(const SystemType& system, const EnclosureType& initial_set, const TimeType& time);

    //! \brief Compute an approximation to the evolution set using the given semantics.
    EnclosureListType evolve(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,system,initial_set,time,semantics);
        return final; }

    //! \brief Compute an approximation to the evolution set under the given semantics.
    EnclosureListType reach(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,system,initial_set,time,semantics);
        reachable.adjoin(final);
        return reachable; }

    //! \brief Compute an approximation to the evolution set under the given semantics, returning the reached and final sets.
    std::pair<EnclosureListType,EnclosureListType> reach_evolve(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=LOWER_SEMANTICS) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,system,initial_set,time,semantics);
        reachable.adjoin(final);
        return make_pair<EnclosureListType,EnclosureListType>(reachable,final); }

    //! \brief Compute an approximation to the evolution set under the given semantics, returning the reached and final sets, and the information
    //! on having disproved.
    tuple<EnclosureListType,EnclosureListType,DisproveData> lower_chain_reach_evolve_disprove(const SystemType& system, const EnclosureType& initial_set,
																					  const TimeType& time, const HybridBoxes& disprove_bounds,
																					  const bool& skip_if_disproved) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        DisproveData falsInfo = this->_lower_evolution_disprove(final,reachable,intermediate,system,initial_set,time,
														   	    disprove_bounds,skip_if_disproved);
        reachable.adjoin(final);
        return make_tuple<EnclosureListType,EnclosureListType,DisproveData>(reachable,final,falsInfo); }

  protected:
    virtual void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                            const SystemType& system, const EnclosureType& initial, const TimeType& time,
                            Semantics semantics) const;

    virtual DisproveData _lower_evolution_disprove(EnclosureListType& final, EnclosureListType& reachable,
										   EnclosureListType& intermediate, const SystemType& system,
										   const EnclosureType& initial, const TimeType& time,
										   const HybridBoxes& disprove_bounds, bool skip_if_disproved) const;

    virtual void _upper_evolution_continuous(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                            const SystemType& system, const EnclosureType& initial, const TimeType& time) const;

    virtual void _evolution_step(std::list< HybridTimedSetType >& working_sets,
                                  EnclosureListType& reachable, EnclosureListType& intermediate,
                                  const SystemType& system, const HybridTimedSetType& current_set, const TimeType& time,
                                  Semantics semantics) const;

    virtual void _upper_evolution_continuous_step(std::list< HybridTimedSetType >& working_sets,
                                  				  EnclosureListType& reachable, EnclosureListType& intermediate,
                                  				  const SystemType& system, const HybridTimedSetType& current_set, const TimeType& time) const;

    virtual DisproveData _lower_evolution_disprove_step(std::list< HybridTimedSetType >& working_sets,
												EnclosureListType& reachable, EnclosureListType& intermediate,
												const SystemType& system, const HybridTimedSetType& current_set, const TimeType& time,
												const HybridBoxes& disprove_bounds) const;

  protected:
    TimeModelType crossing_time(VectorFunction guard, const FlowSetModelType& flow_set) const;

    Interval normal_derivative(VectorFunction guard, const FlowSetModelType& flow_set, const TimeModelType& crossing_time) const;

    /*! \brief Computes the initially active events
     * \details Produces one entry for each possibly initially active event. Adds a convenience summary entry for the generic blocking_event,
     * which has activity FALSE iff no other entry exists, INDETERMINATE if other entries exist but no definitely active event exists, TRUE otherwise.
     */
    void compute_initially_active_events(std::map<DiscreteEvent,tribool>&,
                                         const std::map<DiscreteEvent,VectorFunction>&,
                                         const ContinuousEnclosureType&) const;

    bool has_nonnegative_crossing(const std::map<DiscreteEvent,VectorFunction>& blocking_guards,
								  const VectorFunction dynamic,
								  const Box& set_bounds) const;

    bool is_enclosure_to_be_discarded(const ContinuousEnclosureType& enclosure,
            					 	  const std::map<DiscreteEvent,VectorFunction>& blocking_guards,
            					 	  const VectorFunction& dynamic,
            					 	  Semantics semantics) const;

    void compute_flow_model(FlowSetModelType&, BoxType&, Float&, VectorFunction, 
                            const SetModelType&, const TimeModelType&, Float) const;

    void compute_crossing_time_and_direction(TimeModelType&, Interval&,
                                             VectorFunction guard, const FlowSetModelType& flow_set) const;

    void compute_eventBlockingTimes_and_nonTransverseEvents(std::map<DiscreteEvent,TimeModelType>&, std::set<DiscreteEvent>&,
                                 const std::map<DiscreteEvent,VectorFunction>& blocking_guards,
                                 const FlowSetModelType& flow_set_model) const;

    void compute_blockingTime_and_relatedEvents(std::set<DiscreteEvent>&, TimeModelType&,
                               const std::map<DiscreteEvent,TimeModelType>&) const;

    void compute_activationTimes(std::map<DiscreteEvent,tuple<TimeModelType,TimeModelType> >& activation_times,
                                  const std::map<DiscreteEvent,VectorFunction>& activations,
                                  const FlowSetModelType& flow_set_model,
                                  const TimeModelType& blocking_time_model,
                                  const Semantics semantics) const;

  private:

    bool _is_reachableSet_outside_disproveBounds(const uint numDivisions,
					   const TaylorSet& reachable_set,
					   const Box& disprove_bounds) const;

    bool _isEnclosureTooLarge(const SetModelType& initial_set_model) const;

    void _evolution_add_initialSet(std::list< HybridTimedSetType >& working_sets,
    							   const EnclosureType& initial_set) const;

    void _add_models_subdivisions_autoselect(std::list< HybridTimedSetType >& working_sets,
    		  	  	  	  	  	  	  		 const SetModelType& initial_set_model,
    		  	  	  	  	  	  	  		 const TimeModelType& initial_time_model,
    		  	  	  	  	  	  	  		 const DiscreteState& initial_location,
    		  	  	  	  	  	  	  		 const EventListType& initial_events) const;

    void _add_models_subdivisions_time(std::list< HybridTimedSetType >& working_sets,
    		  	  	  	  	  	  	   const SetModelType& initial_set_model,
    		  	  	  	  	  	  	   const TimeModelType& initial_time_model,
    		  	  	  	  	  	  	   const DiscreteState& initial_location,
    		  	  	  	  	  	  	   const EventListType& initial_events) const;

    void _add_subdivisions(std::list< HybridTimedSetType >& working_sets,
    					   const array< TimedSetModelType >& subdivisions,
    					   const DiscreteState& initial_location,
    					   const EventListType& initial_events,
    					   const uint dimension) const;

    void _logStepAtVerbosity1(const std::list<HybridTimedSetType>& working_sets,
    					 const EnclosureListType& reach_sets,
    					 const EventListType& initial_events,
    					 const TimeModelType& initial_time_model,
    					 const SetModelType& initial_set_model,
    					 const DiscreteState& initial_location) const;

    void _computeEvolutionForEvents(std::list< HybridTimedSetType >& working_sets,
			   	   	   	   	   	    EnclosureListType& intermediate_sets,
			   	   	   	   	   	    const SystemType& system,
			   	   	   	   	   	    const DiscreteState& location,
			   	   	   	   	   	    const std::set<DiscreteEvent>& blocking_events,
			   	   	   	   	   	    const EventListType& events,
			   	   	   	   	   	    const ActivationTimesType& activation_times,
			   	   	   	   	   	    const SetModelType& flow_set_model,
			   	   	   	   	   	    const TimeModelType& time_model,
			   	   	   	   	   	    const TimeModelType& blocking_time_model,
			   	   	   	   	   	    const Float& time_step) const;

    void _processInitiallyActiveBlockingEvents_continuous(EnclosureListType& reach_sets,
    												   	    EnclosureListType& intermediate_sets,
    												   	    const std::map<DiscreteEvent,VectorFunction>& invariants,
    												   	    const SetModelType& set_model,
    												   	    const DiscreteState& location) const;

    void _processInitiallyActiveBlockingEvents(std::list< HybridTimedSetType >& working_sets,
			  	  	  	  	  	  	  	  	  EnclosureListType& reach_sets,
        		   	   	   	   	   	   	   	  EnclosureListType& intermediate_sets,
        		   	   	   	   	   	   	   	  const DiscreteState& location,
        		   	   	   	   	   	   	   	  const EventListType& events,
        		   	   	   	   	   	   	   	  const SetModelType& set_model,
        		   	   	   	   	   	   	   	  const TimeModelType& time_model,
        		   	   	   	   	   	   	   	  const std::map<DiscreteEvent,VectorFunction>& invariants,
        		   	   	   	   	   	   	   	  const std::list<DiscreteTransition>& transitions) const;

    void _compute_blocking_info(std::set<DiscreteEvent>& non_transverse_events,
    				  	   std::set<DiscreteEvent>& blocking_events,
    				  	   TimeModelType& blocking_time_model,
    				  	   const TimeModelType& time_step_model,
    				  	   const SetModelType& flow_set_model,
    				  	   const std::map<DiscreteEvent,VectorFunction>& blocking_guards,
    				  	   double SMALL_RELATIVE_TIME) const;

    void _compute_activation_info(std::map<DiscreteEvent,VectorFunction>& activations,
    						 	  ActivationTimesType& activation_times,
    						 	  const std::set<DiscreteEvent>& non_transverse_events,
    						 	  const SetModelType& flow_set_model,
    						 	  const TimeModelType& blocking_time_model,
    						 	  const std::map<DiscreteEvent,VectorFunction>& blocking_guards,
    						 	  const Semantics semantics) const;

    void _compute_and_adjoin_reachableSet(EnclosureListType& reach_sets,
    									 SetModelType& reachable_set,
    									 const DiscreteState& location,
    									 const SetModelType& flow_set_model,
    									 const TimeModelType& zero_time_model,
    									 const TimeModelType& blocking_time_model) const;

    void _logEvolutionStepInitialState(const EventListType& events,
    							  	   const TimeModelType& time_model,
    							  	   const DiscreteState& location,
    							  	   const SetModelType& set_model,
    							  	   const VectorFunction& dynamic,
    							  	   const std::map<DiscreteEvent,VectorFunction>& invariants,
    							  	   const std::list<DiscreteTransition>& transitions,
    							  	   const std::map<DiscreteEvent,VectorFunction>& blocking_guards,
    							  	   const std::map<DiscreteEvent,VectorFunction>& permissive_guards) const;

  protected:
    // Special events
    static const DiscreteEvent starting_event;
    static const DiscreteEvent finishing_event;
    static const DiscreteEvent blocking_event;

 private:
    boost::shared_ptr< EvolutionSettingsType > _settings;
    boost::shared_ptr< CalculusInterface<TaylorModel> > _toolbox;
};

/*! \brief Whether a box \a set_bounds under a given \a dynamic positively crosses an \a activation. */
tribool positively_crossing(const Box& set_bounds,
							const VectorFunction& dynamic,
							const ScalarFunction& activation);


} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_IMAGE_H
