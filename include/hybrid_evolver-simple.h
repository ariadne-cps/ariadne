/***************************************************************************
 *            hybrid_evolver-simple.h
 *
 *  Copyright  2009  Pieter Collins
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

/*! \file hybrid_evolver-simple.h
 *  \brief Evolver for hybrid systems with only transverse urgent guards.
 */

#ifndef ARIADNE_HYBRID_EVOLVER_SIMPLE_H
#define ARIADNE_HYBRID_EVOLVER_SIMPLE_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "tuple.h"

#include "hybrid_set.h"

#include "hybrid_enclosure.h"
#include "hybrid_automaton_interface.h"
#include "evolver_interface.h"
#include "evolver_base.h"
#include "evolution_parameters.h"

#include "logging.h"

namespace Ariadne {

struct TransitionData;

/*! \brief A class for computing the evolution of a hybrid system.
 *
 * The actual evolution steps are performed by the HybridEvolver class.
 */
class SimpleHybridEvolver
    : public EvolverBase<HybridAutomatonInterface,HybridEnclosure>
    , public Loggable
{
  public:
    typedef ContinuousEvolutionParameters EvolutionParametersType;
    typedef HybridAutomatonInterface::TimeType TimeType;
    typedef int IntegerType;
    typedef Float RealType;
    typedef std::vector<DiscreteEvent> EventListType;
    typedef HybridAutomatonInterface SystemType;
    typedef TaylorConstrainedImageSet ContinuousEnclosureType;
    typedef HybridEnclosure HybridEnclosureType;
    typedef HybridEnclosureType EnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<HybridEnclosure> EnclosureListType;
    typedef Float ContinuousTimeType;
  public:

    //! \brief Default constructor.
    SimpleHybridEvolver();

    //! \brief Construct from parameters using a default integrator.
    SimpleHybridEvolver(const EvolutionParametersType& parameters);

    /*! \brief Make a dynamically-allocated copy. */
    SimpleHybridEvolver* clone() const { return new SimpleHybridEvolver(*this); }

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

    virtual
    void
    _upper_evolution_step(List<HybridEnclosure>& working,
                          ListSet<HybridEnclosure>& final,
                          ListSet<HybridEnclosure>& reachable,
                          ListSet<HybridEnclosure>& intermediate,
                          SystemType const& system,
                          TimeType const& maximum_time) const;

    virtual
    VectorIntervalFunction
    _compute_flow(VectorFunction vector_field,
                  Box const& initial_set,
                  Float& step_size) const;

    virtual
    void
    _compute_active_events(Set<DiscreteEvent>& active_events,
                           VectorFunction const& dynamic,
                           Map<DiscreteEvent,ScalarFunction> const& guard_functions,
                           VectorIntervalFunction const& flow,
                           HybridEnclosure const& initial_set) const;

    virtual
    void
    _compute_crossing_times(Map<DiscreteEvent,ScalarIntervalFunction>& crossing_times,
                            Set<DiscreteEvent> const& active_events,
                            Map<DiscreteEvent,ScalarFunction> const& guard_functions,
                            VectorIntervalFunction const& flow,
                            IntervalVector const& domain) const;

    virtual
    HybridEnclosure
    _process_initial_events(List<HybridEnclosure>& working_sets,
                            HybridEnclosure const& starting_set,
                            Map<DiscreteEvent,TransitionData> const& transitions) const;

    virtual
    void
    _apply_time_step(HybridEnclosure& reach_set, HybridEnclosure& evolve_set,
                     ListSet<HybridEnclosure>& final_sets, List<HybridEnclosure>& jump_sets,
                     HybridEnclosure const& initial_set, Float const& final_time,
                     VectorIntervalFunction const& flow, ScalarIntervalFunction const& evolution_time,
                     Set<DiscreteEvent> const& active_events, Map<DiscreteEvent,ScalarIntervalFunction> const& crossing_times,
                     Map<DiscreteEvent,TransitionData> const& transitions) const;

    virtual
    ScalarIntervalFunction
    _evolution_time(ScalarIntervalFunction const& maximum_evolution_time,
                    Map<DiscreteEvent,ScalarIntervalFunction> const& crossing_times) const;

    ScalarIntervalFunction
    _crossing_time(const ScalarFunction& guard,
                   const VectorFunction& dynamic,
                   const Box& initial,
                   Float maximum_time);
  private:
    boost::shared_ptr< EvolutionParametersType > _parameters;
    //boost::shared_ptr< EvolutionProfiler >  _profiler;
};


class DeterministicHybridEvolver
    : public SimpleHybridEvolver
{
  protected:
    virtual void
    _upper_evolution_step(List<HybridEnclosure>& working,
                          ListSet<HybridEnclosure>& final,
                          ListSet<HybridEnclosure>& reachable,
                          ListSet<HybridEnclosure>& intermediate,
                          SystemType const& system,
                          TimeType const& maximum_time) const;

    virtual
    ScalarIntervalFunction
    _evolution_time(ScalarIntervalFunction const& maximum_evolution_time,
                    Map<DiscreteEvent,ScalarIntervalFunction> const& crossing_times) const;

};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_SIMPLE_H
