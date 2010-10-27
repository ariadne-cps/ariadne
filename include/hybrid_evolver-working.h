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

// A derived class of VectorIntervalFunction representing the flow $\phi(x,t)\f$ of a
// differential equations \f$\dot{x}=f(x)\f$.
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
    /*! \brief Make a dynamically-allocated copy. */
	virtual HybridEvolverInterface* clone() const = 0;
};

struct TransitionData;
struct CrossingData;
struct TimingData;
struct InitialData;
struct FinalData;
struct EvolutionData;


//! \brief Base routines for hybrid evolution.
//! Includes routines for extracting system information in a mode,
//! applying initial events, computing crossing times,
//! computing the evolution time, and applying a time step.
//!
//! \sa The \link hybrid_evolution_page Hybrid Evolution Methods \endlink page
//! for details of the basic algorithm.
//!
//! NOTE: This class currently contains both a skeleton for hybrid evolution,
//! and default implementations of many of the routines. It would be cleaner
//! to have an interface class without these default implementations.
//! However, which operations are needed may depend on the higher-level
//! operations involved, so these methods should also not go in the
//! interface.
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

    //@{
	//! \name Constructors and descructors

	//! \brief Default constructor.
    HybridEvolverBase();

    //! \brief Construct from parameters using a default integrator.
    HybridEvolverBase(const EvolutionParametersType& parameters);

    /*! \brief Make a dynamically-allocated copy. */
    HybridEvolverBase* clone() const = 0;

    //@}

	//@{
    //! \name Parameters controlling the evolution.

	//! \brief A reference to the parameters controlling the evolution.
    EvolutionParametersType& parameters() { return *this->_parameters; }
	//! \brief A constant reference to the parameters controlling the evolution.
	const EvolutionParametersType& parameters() const { return *this->_parameters; }

    //@}


    //@{
    //! \name Main evolution functions.

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
    //@}

  protected:
    //! \brief Internal wrapper routing for computing the evolution.
    //!
    //! \param final Output parameter containing the final (evolve) set.
    //! \param reachable Output parameter containing the reachable set.
    //! \param intermediate Output parameter containing the intermediate sets
    //!   attained after every time step.
    //! \param system The hybrid system under consideration.
	//! \param initial The inital set of the evolution.
	//! \param time The maximum time of evolution; either specifies the stopping time
    //!   or the maximum number of steps.
	//! \param semantics The semantics used for the solution trajectories.
    //!   Either \a #LOWER_SEMANTICS, in which case trajectories terminate at
    //!   discontinuities, or #UPPER_SEMANTICS, in which case all branches
	//!   are taken.
	//! \param reach A flag indicating whether the reachable sets should
	//!   be computed.
	virtual void _evolution(ListSet<HybridEnclosure>& final, ListSet<HybridEnclosure>& reachable, ListSet<HybridEnclosure>& intermediate,
                            const SystemType& system, const EnclosureType& initial, const TimeType& time,
                            Semantics semantics, bool reach) const;

	//! \brief Compute the evolution within a single discrete mode.
	//!
	//! Takes a set from evolution_data.working_sets and processes the initial
	//! events using _process_initial_events.
	//! The resulting set is then evolved using _upper_evolution_step until
	//! either the maximum continuous time is reached, or no furter continuous
	//! evolution is possible.
    virtual
    void
    _upper_evolution_flow(EvolutionData& evolution_data,
                          SystemType const& system,
                          TimeType const& maximum_time) const;

	//! \brief Performs an evolution step on one of the sets listed in \a
	//! evolution_date.starting_sets.
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

    //! \brief Extracts the transitions valid in \a location of \a system.
	//! \details For some systems, extracting the information about the
	//! guard and reset functions may be expensive. This routine collects all
	//! information in one place
	virtual
    Map<DiscreteEvent,TransitionData>
    _extract_transitions(DiscreteLocation const& location,
                         HybridAutomatonInterface const& system) const;

	//! \brief Compute the flow of the \a vector_field, starting in the
	//! given \a initial_set, and with a time step as large as
	//! \a maximum_time_set if possible. The result is a function patch
	//! defined on a domain \f$B\times [0,h]\f$, where \f$B\f$ is the bounding
	//! box for the set, and \f$h\f$ is the step size actually used.
    virtual
    VectorIntervalFunction
    _compute_flow(VectorFunction vector_field,
                  Box const& initial_set,
                  const Float& maximum_step_size) const;

    //! \brief Compute the active events for the \a flow \f$\phi\f$ with
    //! time step \f$h\f$ starting in the given \a starting_set \f$S\f$.
	//! A event \f$e\f$ is \em active if \f$g_e(\phi(x,t))\geq0\f$
	//! for some \f$x\in S\f$ and \f$t\in[0,h]\f$.
	//! The result is allowed to contain events which are not active
	//! during the time step, since the exact set may be too difficult to
	//! compute.
	virtual
    Set<DiscreteEvent>
    _compute_active_events(VectorFunction const& dynamic,
                           Map<DiscreteEvent,ScalarFunction> const& guards,
                           VectorIntervalFunction const& flow,
                           HybridEnclosure const& starting_set) const;

    //! \brief Compute data on how trajectories of the \a flow
	//! cross the given \a guards.
	//! For example, trajectories may have a #INCREASING_CROSSING crossing with
	//! a given crossing_time as a function of the initial state.
	//!
	//! \details The \a dynamic is given explicitly to facilitate computation
	//! of crossing times.
	//! A crossing_time computed must be valid over the \a initial_set, though
	//! will often be computed over the entire bounding box.
	//!
	//! If the crossing is #INCREASING_CROSSING, then the \a crossing_time \f$\gamma\f$
	//! should be computed, satisfying \f$g(\phi(x,\gamma(x)))=0\f$.
	//! If the crossing is #CONCAVE_CROSSING, the the \a critical_time (\a maximum_time) \f$\mu\f$
    //! should be computed, satisfying \f$L_{f}g(\phi(x,\mu(x)))=0\f$.
    virtual
    Map<DiscreteEvent,CrossingData>
    _compute_crossings(Set<DiscreteEvent> const& active_events,
                       VectorFunction const& dynamic,
                       Map<DiscreteEvent,ScalarFunction> const& guards,
                       VectorIntervalFunction const& flow,
                       HybridEnclosure const& initial_set) const;

    //! \brief Compute data related to the time of evolution related to the
	//! maximum time step allowed by the flow, the event times, and the
	//! final time which should be attained at the end of the evolution.
	//!
	//! \post
	//! If the step_kind is CREEP_STEP, then the evolution_time must evaluate
	//!   a value in [0,h] over the initial_set, where \f$h\f$ is
	//!   the \a step_size. In other words, \f$\varepsilon(\xi(s))\in[0,h]\f$
	//!   whenever \f$s\f$ is a valid parameter.
    //! If the step_kind is UNWIND_STEP, then \f$(\omega(s)-\tau(s))\in[0,h]\f$
	//!   whenever \f$s\f$ is a valid parameter.
	//! If the step_kind is FINAL_STEP, then \f$t_{\max}-\tau(s)\leq h\f$
	//!   whenever \f$s\f$ is a valid parameter.
	//!
	//! TODO: By reducing the step time, some events may become inactive.
	//! It should be possible to reflect this in the output.
	//!
	//! TODO: The FINAL_STEP choice might be inappropriate; maybe this should
	//! be replaced with an extra flag representing whether the final time
	//! could possibly be reached. However, it might be simpler to just check
	//! the final conditions latex, or
	virtual
    TimingData
    _compute_timing(Set<DiscreteEvent>& active_events,
                    Real final_time,
                    VectorIntervalFunction const& flow,
                    Map<DiscreteEvent,CrossingData> const& crossings,
                    HybridEnclosure const& initial_set) const;


    //! \brief Process the \a initial_set to find any
	//!  <em>initially active</em> events, and compute the relevant \em jump sets
	//!  and the <em>flowable</em> or <em>progress</em> set.
	//!
	//! First, and \em invariants \f$i_e(x)\leq0\f$ are processed to give
	//! an \em allowable or \em invariant
	//! \f$I=\{ x\in S \mid \forall e\in E_{\mathrm{inv}},\ i_e(x)\leq0\}\f$.
	//! Then for each urgent or permissive transition, the \em active set
	//! \f$A_e = \{ x\in I \mid g_e(x)\geq 0\}\f$ and jump set
	//! \f$ _e = r_e(A_e)\f$ are computed.
	//! Finally, any progress predicates \f$p_e(x)\leq0\f$ are processed
	//! (including those derived from urgent events),
	//! to give the flowable set
	//! \f$P=\{ x\in I \mid \forall e\in E_{\mathrm{tcp}},\ p_e(x)\leq0\}\f$.
	virtual
    void
    _process_initial_events(EvolutionData& evolution_data,
                            HybridEnclosure const& initial_set,
                            Map<DiscreteEvent,TransitionData> const& transitions) const;

    //! \brief Process the \a starting_set \f$S\f$
	//! based on the previously computed \a flow, \a timing_data
	//! and \a crossing_data to compute the successor sets in the evolution.
	//! Compute the resulting \a jump, \a finishing, \a final and \a reach
	//! sets after a single time step.
	//!
	//! \pre \f$p_e(x)\leq0\f$ for any invariant or progress predicate
	//! and any \f$x\in S\f$
	//! \pre The step time \f$\delta(s)\f$ in terms of the flow
	//! parameters satisfies \f$\delta(s)\in[0,h]\f$ whenever \f$\rho(s)\in C\f$.
	//!
	//! The evolution time (in terms of the flow parameters) is given by
	//! \f$\delta(s)=\varepsilon(\xi(s))\f$ for a \a CREEP_STEP
	//! or \f$\delta(s)=\omega(s)-\tau(s)\f$ for a \a FLOW_STEP.
	//!
	//! Let \f$E_\mathrm{blk}\f$ be set set of events \em blocking
	//! continuous evolution.
	//! - The \em reach set is given by
	//!   \f[ R=\{ \phi(\xi(s),t);\;\tau(s)+t \mid (x,t)\in D\times[0,h]
	//!       \mid \rho(s)\in C \ \wedge\ t\leq\delta(s)
	//!       \ \wedge\ \tau(s)+t\leq t_{\max}
	//!	      \ \wedge\ \forall e\in E_{\mathrm{blk}},
	//!         \,\max_{u\in[0,t]} g_e(\phi(x,u))\leq 0 \} \f]
	//! - The \em finishing set \f$E\f$ is given by the additional constraint
	//!   \f$t=\delta(s)\f$ on the reach set,
	//!   or equivalently, \f$\tau(s)+t=\omega(s)\f$.
	//! - The \em final set \f$F\f$ is given by the additional constraint
	//!   \f$\tau(s)+t=t_{\max}\f$ on the reach set.
	//! - The \em active sets \f$A_e\f$ are given by the additional constraint
	//!   \f$g(\phi(\xi(s),t))=0\f$, and the jump sets by \f$J_e=r_e(A_e)\f$.
	//! In the implementation, these sets are simplified to reduce the number
	//! of extra variables and the number of constraints.
	//!
	//! \note The \a jump, \a finishing, \a final and \a reach sets may need to
	//! be constructed as unions of Enclosure sets. The is because when dealing
	//! with #CONCAVE_CROSSING events, it may be the case that these sets are split into
	//! two.
	//!
	//! TODO: This method might be better split into even simpler events.
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
    //! \brief The kind of event (invariant, urgent, etc) determining what
	//! happens when the event is active.
	EventKind event_kind;
    //! \brief The guard function of the event, the event being active when
	//! \f$g(x)\geq0\f$.
	ScalarFunction guard_function;
    //! \brief The target location of the transition.
	DiscreteLocation target;
    //! \brief The reset function \f$x'=r(x)\f$ of the transition.
	VectorFunction reset_function;
    //TransitionData() { }
    //TransitionData(DiscreteLocation t, ScalarFunction g, VectorFunction r)
    //    : target(t), guard_function(g), reset_function(r) { }
};


//! \brief The way trajectories of the flow \f$\phi(x_0,t)\f$ of \f$\dot{x}=f(x)\f$ cross the guard set \f$g(x)=0\f$.
//! The rate of change of the guard function changes along flow lines is given by
//! \f$d g(x(t))/dt = L_{f}g(x(t))\f$ where the Lie derivative \f$L_{f}g\f$ is defined by \f$L_{f}g(x) = (\nabla g\cdot f)(x)\f$.
//!
//! Defaults to DEGENERATE_CROSSING whenever the crossing kind has not been resolved;
//! this need not necessarily mean that the crossing is degenerate, but may
//! also imply that the crossing information is too expensive or sensitive to
//! compute.
//! \relates HybridEvolverInterface \relates CrossingData
enum CrossingKind {
    DEGENERATE_CROSSING, //!< The crossing may be degenerate to second order.
    NEGATIVE_CROSSING, //!< The guard function is negative on the flow domain. No event occurs.
    POSITIVE_CROSSING, //!< The guard function is negative on the domain. The event occurs immediately (if urgent) or at all times (if permissive).
    INCREASING_CROSSING, //!< The guard function is strictly increasing along flow lines at a crossing point.
        //! i.e. \f$\frac{d}{dt}g(\phi(x_0,t))>0\f$ whenever \f$g(\phi(x_0,t))=0\f$.
        //! Implied by \f$L_{f}g \geq 0\f$.
        //! Given an initial point \f$x_0\f$, the crossing time \f$\gamma(x_0)\f$.
    DECREASING_CROSSING, //!< The guard function is strictly decreasing along flow lines.
    CONCAVE_CROSSING, //!< The guard function is positive over at most an interval.
        //! Implied by concavity along flow lines, which is equivalent to \f$L_{f}^{2} g < 0\f$ within the reached set.
        //! The time at which the guard function reaches a maximum for the evolution starting at \f$x_0\f$ is denoted by \f$\mu(x_0)\f$
    CONVEX_CROSSING //!< The guard function is negative over at most an interval. Implied by convexity along flow lines.
        //! Implied by convexity along flow lines, which is equivalent to \f$L_{f}^{2} g > 0\f$ within the reached set.
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
    //! \brief The way in which the guard function changes along trajectories
	//! during a crossing. e.g. INCREASING_CROSSING
	CrossingKind crossing_kind;
    //! \brief The time \f$\gamma(x)\f$ at which the crossing occurs,
    //! as a function of the initial point in space. Satisfies \f$g(\phi(x,\gamma(x)))=0\f$.
    ScalarIntervalFunction crossing_time;
    //! \brief The time \f$\mu(x)\f$ at which the guard function reaches a maximum or minimum
    //! i.e. \f$L_{f}g(\phi(x,\mu(x))) = 0\f$.
    ScalarIntervalFunction critical_time;
};
std::ostream& operator<<(std::ostream& os, const CrossingData& crk);

//! \brief The kind of step taken in the evolution. Determines how the evolution time is specified.
//! \details Assumes that the starting set and time is given as a subset of \f$\{ \xi(s); \tau(s) \mid s\in D\}\f$
//! where \f$s\in D\f$ is a parameter, \f$x=\xi(s)\f$ is the state corresponding to parameter \f$s\f$, and \f$t=\tau(s)\f$
//! is the time the point has so far been evolved for. Assumes that the flow is given by a function \f$x'=\phi(x,t)\f$,
//! typically only defined over a bounded set of space and time.
//! \relates TimingData
enum StepKind {
    FULL_STEP, //!< The step is taken for a fixed time \f$h\f$. The actual step length depends only on the starting state.
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),h)\f$ and \f$\tau'(s)=\tau(s)+h\f$.
    CREEP_STEP, //!< The step is taken for a time \f$\varepsilon(x)\f$ depending only on the starting state.
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),\varepsilon(\xi(s)))\f$ and \f$\tau'(s)=\tau(s)+\varepsilon(\xi(s))\f$.
    PARTIAL_STEP, //!< The step is taken for a time \f$\delta(s)\f$ depending on the parameterisation of the set.
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),\delta(s))\f$ and \f$\tau'(s)=\tau(s)+\delta(s)\f$.
    UNWIND_STEP, //!< The step is taken up to a time \f$\omega(s)\f$ depending on the parameterisation of the starting set.
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),\omega(s)-\tau(s))\f$ and \f$\tau'(s)=\omega(s)\f$.
    FINAL_STEP, //!< The step is taken up to the specified evolution time \f$t_{\max}\f$. The actual step length depends on the parameterisation.
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),t_{\max}-\tau(s))\f$ and \f$\tau'(s)=t_{\max}\f$.
};
std::ostream& operator<<(std::ostream& os, const StepKind& crk);

//! \brief A data type used to store information about the kind of time step taken during hybrid evolution.
//! \relates HybridEvolverInterface
struct TimingData
{
    StepKind step_kind; //!< The kind of step taken in the evolution
    Real final_time; //!< The time \f$t_{\max}\f$ specified as the final time of the evolution trace, and used in a \a FINAL_STEP step.
    Float step_size; //!< The step size \f$h\f$ used in a \a FULL_STEP time step.
    ScalarIntervalFunction spacial_evolution_time; //!< The evolution time \f$\varepsilon(x)\f$ used in a \a CREEP_STEP time step.
    ScalarIntervalFunction finishing_time; //!< The time \f$\omega(s)\f$ reached after an \a UNWIND_STEP as a function of the parameters.
    ScalarIntervalFunction evolution_time; //!< The time \f$\delta(s)\f$ used in an evolution step. Equal to \f$\varepsilon(\xi(s))\f$ for a \a CREEP_STEP and \f$\omega(s)-\varepsilon(s)\f$ for an \a UNWIND_STEP.
    ScalarIntervalFunction remaining_time; //!< The time \f$\rho(s)\f$ remaining in the evolution. Equal to \f$t_{\max}-\tau(s)\f$.
    Interval time_domain; //!< The time domain, equal to \f$[0,h]\f$.
    ScalarIntervalFunction time_coordinate; //!< The time coordinate function, equal to the identity on \f$[0,h]\f$.
};

//! \brief A data type used to store information about the kind of time step taken during hybrid evolution.
//! \relates HybridEvolverInterface
struct EvolutionData
{
	//! Sets which need to be processed, typically by finding initially active
	//! events, and moving the set into starting_sets.
	//! Typically, the initial set and any jump sets should be placed in this
	//! list.
    List<HybridEnclosure> working_sets;
    //! \brief Sets which have been processed for initially active events.
    //! Sets in this list can be assumed to satisfy all invariants and progress
    //! predicates.
    List<HybridEnclosure> starting_sets;

    //! \brief Sets which have been computed as reach sets for the current
	//! evolution.
	List<HybridEnclosure> reach_sets;
    //! \brief Sets which have been computed as final sets (i.e. satisfying
	//! either the final time or the maximum number of stesp).
	List<HybridEnclosure> evolve_sets;
    //! \brief Intermediate sets reached after each time step. Not relevant for
	//! the result, but useful for plotting, especially for debugging.
	List<HybridEnclosure> intermediate_sets;

    //! \brief The semantics used to compute the evolution. Defaults to UPPER_SEMANTICS.
    Semantics semantics;
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


    //! \brief Apply \a guard_function for \a event to each set in \a sets, using the computed \a crossing_data.
    //! \callgraph
    //! \details The sets are updated with the extra constraints in-place for efficiency.
    //! In the case of a concave tangency, it may be necessary to split the enclosure in two,
    //! one part containing points which eventually cross the guard, the other points which miss the guard.
    //! The extra sets are adjoined to the end of the list of sets to be constrained.
    virtual
    void
    _apply_guard(List<HybridEnclosure>& sets,
                 const VectorIntervalFunction& starting_state,
                 const VectorIntervalFunction& flow,
                 const ScalarIntervalFunction& elapsed_time,
                 const DiscreteEvent event,
                 const ScalarFunction& guard_function,
                 const CrossingData crossing_data,
                 const Semantics semantics) const;

};


} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_WORKING_H
