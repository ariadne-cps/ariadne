/***************************************************************************
 *            hybrid_evolver-working.h
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

/*! \file hybrid_evolver-working.h
 *  \brief Hybrid evolver classes under development
 */

#ifndef ARIADNE_HYBRID_EVOLVER_WORKING_H
#define ARIADNE_HYBRID_EVOLVER_WORKING_H

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

class FlowFunctionPatch
    : public VectorIntervalFunction
{
  public:
    FlowFunctionPatch(const VectorIntervalFunction& f) : VectorIntervalFunction(f) { }
    Float step_size() const { return this->time_domain().upper(); }
    Interval time_domain() const { return this->domain()[this->domain().size()-1]; }
    IntervalVector space_domain() const { return project(this->domain(),Ariadne::range(0,this->domain().size()-1)); }
    IntervalVector codomain() const { return this->VectorIntervalFunction::codomain(); }
};

//! \brief Interface for hybrid evolvers using HybridEnclosure as the enclosure type.
class HybridEvolverInterface
    : public EvolverBase<HybridAutomatonInterface,HybridEnclosure>
{
  public:
    virtual HybridEvolverInterface* clone() const = 0;
};

struct TransitionData;
struct CrossingData;
struct TimingData;
struct InitialData;
struct FinalData;
struct EvolutionData;


//! \brief Base routines for hybrid evolution.
//!
//! Includes routines for extracting system information in a mode,
//! applying initial events, computing crossing times,
//! computing the evolution time, and applying a time step.
class HybridEvolverBase
    : public HybridEvolverInterface
    , public Loggable
{
  public:
    typedef ContinuousEvolutionParameters EvolutionParametersType;
    typedef HybridAutomatonInterface::TimeType TimeType;
    typedef int IntegerType;
    typedef Float RealType;
    typedef std::vector<DiscreteEvent> EventListType;
    typedef HybridAutomatonInterface SystemType;
    typedef ConstrainedImageSet ContinuousEnclosureType;
    typedef HybridEnclosure HybridEnclosureType;
    typedef HybridEnclosureType EnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<HybridEnclosure> EnclosureListType;
    typedef Float ContinuousTimeType;
  public:

    //! \brief Default constructor.
    HybridEvolverBase();

    //! \brief Construct from parameters using a default integrator.
    HybridEvolverBase(const EvolutionParametersType& parameters);

    /*! \brief Make a dynamically-allocated copy. */
    HybridEvolverBase* clone() const = 0;

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
    virtual void _evolution(ListSet<HybridEnclosure>& final, ListSet<HybridEnclosure>& reachable, ListSet<HybridEnclosure>& intermediate,
                            const SystemType& system, const EnclosureType& initial, const TimeType& time,
                            Semantics semantics, bool reach) const;

    virtual
    void
    _upper_evolution_flow(EvolutionData& evolution_data,
                          SystemType const& system,
                          TimeType const& maximum_time) const;

    virtual
    void
    _upper_evolution_step(EvolutionData& evolution_data,
                          VectorFunction const& dynamic,
                          Map<DiscreteEvent,TransitionData> const& transitions,
                          Real const& final_time) const;

    //! \brief Output a one-line summary of the current evolution state to the logging stream.
    virtual
    void
    _log_summary(uint ws, uint rs, HybridEnclosure const& starting_set) const;

    virtual
    Map<DiscreteEvent,TransitionData>
    _extract_transitions(DiscreteLocation const& location,
                         HybridAutomatonInterface const& system) const;

    virtual
    VectorIntervalFunction
    _compute_flow(VectorFunction vector_field,
                  Box const& initial_set,
                  const Float& maximum_step_size) const;

    virtual
    Set<DiscreteEvent>
    _compute_active_events(VectorFunction const& dynamic,
                           Map<DiscreteEvent,ScalarFunction> const& guards,
                           VectorIntervalFunction const& flow,
                           HybridEnclosure const& initial_set) const;

    virtual
    Map<DiscreteEvent,CrossingData>
    _compute_crossings(Set<DiscreteEvent> const& active_events,
                       VectorFunction const& dynamic,
                       Map<DiscreteEvent,ScalarFunction> const& guards,
                       VectorIntervalFunction const& flow,
                       HybridEnclosure const& initial_set) const;

    virtual
    HybridEnclosure
    _compute_transverse_jump_set(HybridEnclosure const& starting_set,
                                DiscreteEvent const& event,
                                VectorIntervalFunction const& flow,
                                ScalarIntervalFunction const& evolution_time,
                                Float const& final_time,
                                Set<DiscreteEvent> const& active_events,
                                Map<DiscreteEvent,TransitionData> const& transitions,
                                Map<DiscreteEvent,CrossingData> const& crossing_data) const;

    virtual
    TimingData
    _compute_timing(Set<DiscreteEvent> const& active_events,
                    Real final_time,
                    VectorIntervalFunction const& flow,
                    Map<DiscreteEvent,CrossingData> const& crossings,
                    HybridEnclosure const& initial_set) const;


    virtual
    void
    _process_initial_events(EvolutionData& evolution_data,
                            HybridEnclosure const& initial_set,
                            Map<DiscreteEvent,TransitionData> const& transitions) const;

    virtual
    void
    _apply_time_step(EvolutionData& evolution_data,
                     HybridEnclosure const& starting_set,
                     VectorIntervalFunction const& flow,
                     TimingData const& timing_data,
                     Map<DiscreteEvent,CrossingData> const& crossing_data,
                     VectorFunction const& dynamic,
                     Map<DiscreteEvent,TransitionData> const& transitions) const = 0;

  private:
    boost::shared_ptr< EvolutionParametersType > _parameters;
    //boost::shared_ptr< EvolutionProfiler >  _profiler;
};


bool is_blocking(EventKind evk);
bool is_activating(EventKind evk);

//! \brief A data type used to store information about a transtion of a hybrid system.
//! \relates HybridEvolverInterface
struct TransitionData
{
    EventKind event_kind;
    ScalarFunction guard_function;
    DiscreteLocation target;
    VectorFunction reset_function;
    //TransitionData() { }
    //TransitionData(DiscreteLocation t, ScalarFunction g, VectorFunction r)
    //    : target(t), guard_function(g), reset_function(r) { }
};


//! \brief The way trajectories of the flow $\phi(x_0,t)$ cross the guard set $g(x)=0$.
//! \relates CrossingData
enum CrossingKind {
    DEGENERATE=0, //!< The crossing may be degenerate to second order.
    NEGATIVE, //!< The guard function is negative on the flow domain. No event occurs.
    POSITIVE, //!< The guard function is negative on the domain. The event occurs immediately (if urgent) or at all times (if permissive).
    INCREASING, //!< The guard function is strictly increasing along flow lines.
    DECREASING, //!< The guard function is strictly decreasing along flow lines.
    CONCAVE, //!< The guard function is positive over at most an interval. Implied by concavity along flow lines.
    CONVEX //!< The guard function is negative over at most an interval. Implied by convexity along flow lines.
};
std::ostream& operator<<(std::ostream& os, const CrossingKind& crk);

//! \brief A data type used to store information about the way flow lines cross a guard \f$g(x)=0\f$.
//! \relates HybridEvolverInterface
struct CrossingData
{
    CrossingData() : crossing_kind() { }
    CrossingData(CrossingKind crk) : crossing_kind(crk) { }
    CrossingData(CrossingKind crk, const ScalarIntervalFunction& crt)
        : crossing_kind(crk), crossing_time(crt) { }
    CrossingKind crossing_kind;
    ScalarIntervalFunction crossing_time;
    ScalarIntervalFunction critical_time;
};

//! \brief The kind of step taken in the evolution
//! \relates HybridEvolverInterface
enum StepKind {
    FULL_STEP, //!< The step is taken for a fixed time \a h. The actual step length depends only on the starting state.
    CREEP_STEP, //!< The step is taken for a time \f$\varepsilon(x)\f$ depending only on the starting state.
    UNWIND_STEP, //!< The step is taken up to a time \f$\omega(s)\f$ depending on the parameterisation of the starting set.
    FINAL_STEP, //!< The step is taken up to the specified evolution time t<sub>max</sub>. The actual step length depends on the parameterisation.
};
std::ostream& operator<<(std::ostream& os, const StepKind& crk);

//! \brief A data type used to store information about the kind of time step taken during hybrid evolution.
//! \relates HybridEvolverInterface
struct TimingData
{
    StepKind step_kind; //!< The kind of step taken in the evolution
    Float step_size; //!< The step size used in a \a FULL_STEP time step.
    ScalarIntervalFunction evolution_time; //!< The evolution time \f$\varepsilon(x)\f$ used in a \a CREEP_STEP time step.
    ScalarIntervalFunction finishing_time; //!< The time \f$\omega(s)\f$ reached after an \a UNWIND_STEP as a function of the parameters.
    Float final_time; //!< The time specified as the final time of the evolution trace, and used in a \a FINAL_STEP step.
};

//! \brief A data type used to store information about the kind of time step taken during hybrid evolution.
//! \relates HybridEvolverInterface
struct EvolutionData
{
    List<HybridEnclosure> working_sets;
    List<HybridEnclosure> starting_sets;

    List<HybridEnclosure> reach_sets;
    List<HybridEnclosure> evolve_sets;
    List<HybridEnclosure> intermediate_sets;
};

//! \brief A class for computing the evolution of a general hybrid system as
//! efficiently as possible. In particular, crossing times are computed
//! explicitly, and the number of additional constraints is minimised.
//!
//! This class is still under development and may not work properly.
class GeneralHybridEvolver
    : public HybridEvolverBase
{
  public:

    GeneralHybridEvolver() : HybridEvolverBase() { ARIADNE_WARN("GeneralHybridEvolver class still under development\n"); }
    GeneralHybridEvolver(const EvolutionParametersType& parameters) : HybridEvolverBase(parameters) { ARIADNE_WARN("GeneralHybridEvolver class still under development\n"); }
    GeneralHybridEvolver* clone() const { return new GeneralHybridEvolver(*this); }

  protected:
    virtual
    void
    _apply_time_step(EvolutionData& evolution_data,
                     HybridEnclosure const& starting_set,
                     VectorIntervalFunction const& flow,
                     TimingData const& timing_data,
                     Map<DiscreteEvent,CrossingData> const& crossing_data,
                     VectorFunction const& dynamic,
                     Map<DiscreteEvent,TransitionData> const& transitions) const;

};

//! \brief A class for computing the evolution of a hybrid system which never
//! computes crossing times explicitly. This means that there is no need to
//! solve algebraic equations. The restrictions on the enclosure sets are
//! specified purely in terms of constraints.
class PureConstraintHybridEvolver
    : public HybridEvolverBase
{
  public:

    PureConstraintHybridEvolver() : HybridEvolverBase() { }
    PureConstraintHybridEvolver(const EvolutionParametersType& parameters) : HybridEvolverBase(parameters) { };
    PureConstraintHybridEvolver* clone() const { return new PureConstraintHybridEvolver(*this); }

  protected:
    virtual
    void
    _apply_time_step(EvolutionData& evolution_data,
                     HybridEnclosure const& starting_set,
                     VectorIntervalFunction const& flow,
                     TimingData const& timing_data,
                     Map<DiscreteEvent,CrossingData> const& crossing_data,
                     VectorFunction const& dynamic,
                     Map<DiscreteEvent,TransitionData> const& transitions) const;

};




//! \brief A hybrid evolver which can only handle systems for which crossings
//! with the guard sets are transverse.
class TransverseHybridEvolver
    : public HybridEvolverBase
{
  public:
    TransverseHybridEvolver* clone() const { return new TransverseHybridEvolver(*this); }
  protected:
    virtual void
    _apply_time_step(EvolutionData& evolution_data,
                     HybridEnclosure const& starting_set,
                     VectorIntervalFunction const& flow,
                     TimingData const& timing_data,
                     Map<DiscreteEvent,CrossingData> const& crossing_data,
                     VectorFunction const& dynamic,
                     Map<DiscreteEvent,TransitionData> const& transitions) const;

};

//! \brief A hybrid evolver which can only handle deterministic systems
//! i.e. systems with only \a urgent or \a impact events.
class DeterministicHybridEvolver
    : public HybridEvolverBase
{
  public:
    DeterministicHybridEvolver* clone() const { return new DeterministicHybridEvolver(*this); }
  protected:
    virtual ScalarIntervalFunction
    _evolution_time(ScalarIntervalFunction const& maximum_evolution_time,
                    Map<DiscreteEvent,ScalarIntervalFunction> const& crossing_times) const;

    virtual void
    _apply_time_step(EvolutionData& evolution_data,
                     HybridEnclosure const& starting_set,
                     VectorIntervalFunction const& flow,
                     TimingData const& timing_data,
                     Map<DiscreteEvent,CrossingData> const& crossing_data,
                     VectorFunction const& dynamic,
                     Map<DiscreteEvent,TransitionData> const& transitions) const;

};


//! \brief A hybrid evolver which can only handle deterministic systems for which crossings
//! with the guard sets are transverse.
//! \sa DeterministicHybridEvolver TransverseHybridEvolver
class DeterministicTransverseHybridEvolver
    : public HybridEvolverBase
{
  public:
    DeterministicTransverseHybridEvolver* clone() const { return new DeterministicTransverseHybridEvolver(*this); }
  protected:
    virtual
    TimingData
    _compute_timing(Set<DiscreteEvent> const& active_events,
                    Real final_time,
                    VectorIntervalFunction const& flow,
                    Map<DiscreteEvent,CrossingData> const& crossings,
                    HybridEnclosure const& initial_set) const;

    virtual void
    _apply_time_step(EvolutionData& evolution_data,
                     HybridEnclosure const& starting_set,
                     VectorIntervalFunction const& flow,
                     TimingData const& timing_data,
                     Map<DiscreteEvent,CrossingData> const& crossing_data,
                     VectorFunction const& dynamic,
                     Map<DiscreteEvent,TransitionData> const& transitions) const;

};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_WORKING_H
