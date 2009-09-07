/***************************************************************************
 *            hybrid_evolver.h
 *
 *  Copyright  2007-8  Alberto Casagrande, Pieter Collins
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

/*! \file hybrid_evolver.h
 *  \brief Evolver for hybrid systems.
 */

#ifndef ARIADNE_HYBRID_EVOLVER_H
#define ARIADNE_HYBRID_EVOLVER_H

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
#include "evolution_parameters.h"

#include "logging.h"

namespace Ariadne {

template<class Sys, class BS> class Evolver;

class VectorTaylorFunction;
class TaylorSet;
typedef std::pair<DiscreteState,TaylorSet> HybridTaylorSet;
class HybridAutomaton;
template<class ES> class Orbit;

class EvolutionParameters;
class TaylorModel;
template<class MDL> class CalculusInterface;

class EvolutionProfiler;

class HybridTime;

typedef int DiscreteEvent;

/*! \brief A class for computing the evolution of a hybrid system.
 *
 * The actual evolution steps are performed by the HybridEvolver class.
 */
class HybridEvolver
    : public EvolverBase<HybridAutomaton,HybridTaylorSet>
    , public Loggable
{
    typedef boost::shared_ptr<const ScalarFunctionInterface> ScalarFunctionPtr;
    typedef boost::shared_ptr<const VectorFunctionInterface> VectorFunctionPtr;
    typedef ScalarFunctionInterface ScalarFunctionType;
    typedef VectorFunctionInterface VectorFunctionType;
    typedef Vector<Interval> BoxType;
    typedef VectorTaylorFunction FunctionModelType;
    typedef VectorTaylorFunction MapModelType;
    typedef VectorTaylorFunction FlowModelType;
    typedef ScalarTaylorFunction ConstraintModelType;
    typedef TaylorModel TimeModelType;
    typedef TaylorSet FlowSetModelType;
    typedef TaylorSet SetModelType;
    typedef TaylorSet TimedSetModelType;
  public:
    typedef ContinuousEvolutionParameters EvolutionParametersType;
    typedef HybridAutomaton::TimeType TimeType;
    typedef int IntegerType;
    typedef Float RealType;
    typedef std::vector<DiscreteEvent> EventListType;
    typedef HybridAutomaton SystemType;
    typedef TaylorSet ContinuousEnclosureType;
    typedef pair<DiscreteState,TaylorSet> HybridEnclosureType;
    typedef HybridEnclosureType EnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<EnclosureType> EnclosureListType;
    typedef Float ContinuousTimeType;
  public:

    //! \brief Default constructor.
    HybridEvolver();

    //! \brief Construct from parameters using a default integrator.
    HybridEvolver(const EvolutionParametersType& parameters);

    /*! \brief Make a dynamically-allocated copy. */
    HybridEvolver* clone() const { return new HybridEvolver(*this); }

    //@{
    //! \name Parameters controlling the evolution.
    //! \brief A reference to the parameters controlling the evolution.
    EvolutionParametersType& parameters() { return *this->_parameters; }
    const EvolutionParametersType& parameters() const { return *this->_parameters; }

    //@}


    //@{
    //! \name Evolution using abstract sets.
    //! \brief Compute an approximation to the orbit set using the given semantics.
    Orbit<EnclosureType> orbit(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const;


    //! \brief Compute an approximation to the evolution set using the given semantics.
    EnclosureListType evolve(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,system,initial_set,time,semantics,false);
        return final; }

    //! \brief Compute an approximation to the evolution set under the given semantics.
    EnclosureListType reach(const SystemType& system, const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,system,initial_set,time,semantics,true);
        return reachable; }

  protected:
    virtual void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                            const SystemType& system, const EnclosureType& initial, const TimeType& time,
                            Semantics semantics, bool reach) const;

    typedef tuple<DiscreteState, EventListType, SetModelType, TimeModelType> HybridTimedSetType;
    virtual void _evolution_step(std::vector< HybridTimedSetType >& working_sets,
                                  EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate,
                                  const SystemType& system, const HybridTimedSetType& current_set, const TimeType& time,
                                  Semantics semantics, bool reach) const;

  protected:
    TimeModelType crossing_time(VectorFunctionPtr guard, const FlowSetModelType& flow_set) const;

    Interval normal_derivative(VectorFunctionPtr guard, const FlowSetModelType& flow_set, const TimeModelType& crossing_time) const;

    void compute_initially_active_events(std::map<DiscreteEvent,tribool>&,
                                         const std::map<DiscreteEvent,VectorFunctionPtr>&,
                                         const ContinuousEnclosureType&) const;

    void compute_flow_model(FlowSetModelType&, BoxType&, Float&,
                            VectorFunctionPtr, const SetModelType&) const;
    void compute_flow_model(FunctionModelType&, BoxType&,
                            VectorFunctionPtr, const BoxType&) const;

    void compute_crossing_time_and_direction(TimeModelType&, Interval&,
                                             VectorFunctionPtr guard, const FlowSetModelType& flow_set) const;

    void compute_blocking_events(std::map<DiscreteEvent,TimeModelType>&, std::set<DiscreteEvent>&,
                                 const std::map<DiscreteEvent,VectorFunctionPtr>& guards,
                                 const FlowSetModelType& flow_set_model) const;

    void compute_blocking_time(std::set<DiscreteEvent>&, TimeModelType&,
                               const std::map<DiscreteEvent,TimeModelType>&) const;

    void compute_activation_events(std::map<DiscreteEvent,tuple<tribool,TimeModelType,tribool> >&,
                                  const std::map<DiscreteEvent,VectorFunctionPtr>& activations,
                                  const FlowSetModelType& flow_set_model) const;

    void compute_activation_times(std::map<DiscreteEvent,tuple<TimeModelType,TimeModelType> >&,
                                  const std::map<DiscreteEvent,VectorFunctionPtr>& activations,
                                  const FlowSetModelType& flow_set_model,
                                  const TimeModelType& blocking_time_model,
                                  const Semantics sematics) const;


  protected:
    // Special events
    static const DiscreteEvent starting_event;
    static const DiscreteEvent finishing_event;
    static const DiscreteEvent blocking_event;
    static const DiscreteEvent final_time_event;

 private:
    boost::shared_ptr< EvolutionParametersType > _parameters;
    boost::shared_ptr< CalculusInterface<TaylorModel> > _toolbox;
    //boost::shared_ptr< EvolutionProfiler >  _profiler;
};



} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_H
