/***************************************************************************
 *            hybrid_evolver.h
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

/*! \file hybrid_evolver.h
 *  \brief Hybrid evolver classes.
 */

#ifndef ARIADNE_HYBRID_EVOLVER_H
#define ARIADNE_HYBRID_EVOLVER_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "tuple.h"

#include "hybrid_time.h"
#include "hybrid_set.h"

#include "configuration_interface.h"
#include "hybrid_enclosure.h"
#include "hybrid_orbit.h"
#include "hybrid_automaton_interface.h"
#include "hybrid_evolver_interface.h"

#include "logging.h"

namespace Ariadne {

template<class X> class FunctionModelFactoryInterface;
typedef FunctionModelFactoryInterface<Interval> IntervalFunctionModelFactoryInterface;

class IntegratorInterface;
class SolverInterface;
class HybridEvolverBaseConfiguration;
class GeneralHybridEvolverConfiguration;

// A derived class of IntervalVectorFunctionModel representing the flow $\phi(x,t)\f$ of a
// differential equations \f$\dot{x}=f(x)\f$.
class FlowFunctionPatch
    : public IntervalVectorFunctionModel
{
  public:
    FlowFunctionPatch(const IntervalVectorFunctionModel& f) : IntervalVectorFunctionModel(f) { }
    Float step_size() const { return this->time_domain().upper(); }
    Interval time_domain() const { return this->domain()[this->domain().size()-1]; }
    IntervalVector space_domain() const { return project(this->domain(),Ariadne::range(0,this->domain().size()-1)); }
    IntervalVector const codomain() const { return this->IntervalVectorFunctionModel::codomain(); }
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
{
    friend class HybridEvolverBaseConfiguration;
  public:
    typedef HybridEvolverBaseConfiguration ConfigurationType;
    typedef IntervalFunctionModelFactoryInterface FunctionFactoryType;
    typedef HybridAutomatonInterface SystemType;
    typedef SystemType::TimeType TimeType;
    typedef TimeType::ContinuousTimeType ContinuousTimeType;
    typedef List<DiscreteEvent> EventListType;
    typedef HybridEnclosure EnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<HybridEnclosure> EnclosureListType;

  public:

    //@{
    //! \name Constructors and destructors

    //! \brief Default constructor.
    HybridEvolverBase(const SystemType& system);

    //! \brief Construct from parameters using a default integrator.
    HybridEvolverBase(
    		const SystemType& system,
            const FunctionFactoryType& factory);

    /*! \brief Make a dynamically-allocated copy. */
    HybridEvolverBase* clone() const = 0;

    /*! \brief Get the internal system. */
    const SystemType& system() const;

    //@}

    //@{
    //! \name Settng and configuration for the class.

    ConfigurationType& configuration();
    const ConfigurationType& configuration() const;

    //! \brief Change the configuration from a \a domain and \a lengths (NOT IMPLEMENTED).
    virtual void reconfigure(const HybridBoxes& domain, const HybridFloatVector& lengths) { }

    /*! \brief The class which constructs functions for the enclosures. */
    const FunctionFactoryType& function_factory() const;
    /*! \brief Set the class which constructs functions for the enclosures. */
    void set_function_factory(const FunctionFactoryType& factory);

    /*! \brief Set the class which integrates the continuous dynamics. */
    void set_integrator(const IntegratorInterface& integrator);
    /*! \brief Set the class which integrates the continuous dynamics. */
    void set_solver(const SolverInterface& solver);

    virtual EnclosureType enclosure(const HybridBox& initial_box) const;
    virtual EnclosureType enclosure(const HybridExpressionSet& set) const;

    bool ALLOW_CREEP; //!< If true, a less-than-full evolution step may be taken to avoid splitting due to partially crossing a guard.
    bool ALLOW_UNWIND; //!< If true, a less-than-full evolution step may be taken to try to restore all time values over the parameter domain to the same value.
    //@}


    //@{
    //! \name Main evolution functions.

    //! \brief Compute an approximation to the orbit set using the given semantics.
    Orbit<EnclosureType> orbit(const EnclosureType& initial_enclosure, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const;
    Orbit<EnclosureType> orbit(const HybridBox& initial_box, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const;
    Orbit<EnclosureType> orbit(const HybridExpressionSet& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const;


    //! \brief Compute an approximation to the evolution set using the given semantics.
    EnclosureListType evolve(const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,time,semantics,false);
        return final; }

    //! \brief Compute an approximation to the evolution set under the given semantics.
    EnclosureListType reach(const EnclosureType& initial_set, const TimeType& time, Semantics semantics=UPPER_SEMANTICS) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate;
        this->_evolution(final,reachable,intermediate,initial_set,time,semantics,true);
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
                            const EnclosureType& initial, const TimeType& time,
                            Semantics semantics, bool reach) const;

    //! \brief Compute the evolution within a single discrete mode.
    //!
    //! Takes a set from evolution_data.initial_sets and processes the initial
    //! events using _process_initial_events.
    //! The resulting set is then evolved using _upper_evolution_step until
    //! either the maximum continuous time is reached, or no furter continuous
    //! evolution is possible.
    virtual
    void
    _evolution_in_mode(EvolutionData& evolution_data,
                       TimeType const& maximum_time) const;

    //! \brief Performs an evolution step on one of the sets listed in \a
    //! evolution_data.initial_sets.
    virtual
    void
    _evolution_step(EvolutionData& evolution_data,
                    RealVectorFunction const& dynamic,
                    Map<DiscreteEvent,TransitionData> const& transitions,
                    Real const& final_time) const;

    //! \brief Extracts the transitions valid in \a location of \a system.
    //! \details For some systems, extracting the information about the
    //! guard and reset functions may be expensive. This routine collects all
    //! information in one place
    virtual
    Map<DiscreteEvent,TransitionData>
    _extract_transitions(DiscreteLocation const& location) const;

    //! \brief Compute the flow of the \a vector_field, starting in the
    //! given \a initial_set, and with a time step as large as
    //! \a maximum_time_set if possible. The result is a function patch
    //! defined on a domain \f$B\times [0,h]\f$, where \f$B\f$ is the bounding
    //! box for the set, and \f$h\f$ is the step size actually used.
    virtual
    IntervalVectorFunctionModel
    _compute_flow(RealVectorFunction vector_field,
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
    _compute_active_events(RealVectorFunction const& dynamic,
                           Map<DiscreteEvent,RealScalarFunction> const& guards,
                           IntervalVectorFunctionModel const& flow,
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
                       RealVectorFunction const& dynamic,
                       Map<DiscreteEvent,RealScalarFunction> const& guards,
                       IntervalVectorFunctionModel const& flow,
                       HybridEnclosure const& initial_set) const;

    //! \brief Compute data related to the time of evolution related to the
    //! maximum time step allowed by the flow, the event times, and the
    //! final time which should be attained at the end of the evolution.
    //!
    //! \post
    //! If the step_kind is SPACETIME_DEPENDENT_EVOLUTION_TIME, then the evolution_time must evaluate
    //!   a value in [0,h] over the initial_set, where \f$h\f$ is
    //!   the \a step_size. In other words, \f$\varepsilon(\xi(s),\tau(s))\in[0,h]\f$
    //!   whenever \f$s\f$ is a valid parameter.
    //! If the step_kind is PARAMETER_DEPENDENT_FINISHING_TIME, then \f$(\omega(s)-\tau(s))\in[0,h]\f$
    //!   whenever \f$s\f$ is a valid parameter.
    //! If the step_kind is CONSTANT_FINISHING_TIME, then \f$t_{\max}-\tau(s)\leq h\f$
    //!   whenever \f$s\f$ is a valid parameter.
    //!
    //! TODO: By reducing the step time, some events may become inactive.
    //! It should be possible to reflect this in the output.
    //!
    //! TODO: The CONSTANT_FINISHING_TIME choice might be inappropriate; maybe this should
    //! be replaced with an extra flag representing whether the final time
    //! could possibly be reached. However, it might be simpler to just check
    //! the final conditions latex, or
    virtual
    TimingData
    _estimate_timing(Set<DiscreteEvent>& active_events,
                     Real final_time,
                     IntervalVectorFunctionModel const& flow,
                     Map<DiscreteEvent,CrossingData>& crossings,
                     Map<DiscreteEvent,TransitionData> const& transitions,
                     HybridEnclosure const& initial_set) const;


    //! \brief Simplify the description of a set or allow for reconditioning
    //! of numerical error.
    virtual
    void
    _recondition(HybridEnclosure& set) const;

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

    //! \brief Apply the \a flow to the \a set up to the time specified by \a timing_data
    //! to obtain the single-step unconstrained reachable set.
    virtual
    void
    _apply_reach_step(HybridEnclosure& set,
                      IntervalVectorFunctionModel const& flow,
                      TimingData const& timing_data) const;

    //! \brief Apply the \a flow to the \a set for the time specified by \a timing_data
    //! to obtain the single-step unconstrained evolved set.
    virtual
    void
    _apply_evolve_step(HybridEnclosure& set,
                       IntervalVectorFunctionModel const& flow,
                       TimingData const& timing_data) const;

    //! \brief Apply the \a flow to the \a set for to reach the
    //! guard specified by \a transition_data using guidelines provided by \a crossing_data.
    //!
    //! \details
    //! NOTE: The \a dynamic is used to compute the Lie derivative of
    //! the flow to ensure positive crossing in the degenerate case. It should not really be necessary.
    virtual
    void
    _apply_guard_step(HybridEnclosure& set,
                      RealVectorFunction const& dynamic,
                      IntervalVectorFunctionModel const& flow,
                      TimingData const& timing_data,
                      TransitionData const& transition_data,
                      CrossingData const& crossing_data,
                      const Semantics semantics) const;

    //! \brief Apply the invariants in \a transition_data to the set \a set.
    //! \callgraph
    virtual
    void
    _apply_invariants(HybridEnclosure& set,
                      const Map<DiscreteEvent,TransitionData>& transition_data) const;

    //! \brief Apply \a guard_function for \a event to each set in \a sets, using the computed \a crossing_data.
    //! \callgraph
    //! \details The sets are updated with the extra constraints in-place for efficiency.
    //! In the case of a concave tangency, it may be necessary to split the enclosure in two,
    //! one part containing points which eventually cross the guard, the other points which miss the guard.
    //! The extra sets are adjoined to the end of the list of sets to be constrained.
    //! NOTE: The \a evolve_time should probably be in crossing_data, and should not really be necessary.
    virtual
    void
    _apply_guard(List<HybridEnclosure>& sets,
                 const HybridEnclosure& starting_set,
                 const IntervalVectorFunctionModel& flow,
                 const IntervalScalarFunctionModel& evolve_time,
                 const TransitionData& transition_data,
                 const CrossingData crossing_data,
                 const Semantics semantics) const;

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
    //! \f$\delta(s)=\varepsilon(\xi(s))\f$ for a \a SPACE_DEPENDENT_EVOLUTION_TIME
    //! or \f$\delta(s)=\omega(s)-\tau(s)\f$ for an \a PARAMETER_DEPENDENT_FINISHING_TIME.
    //!
    //! Let \f$E_\mathrm{blk}\f$ be set set of events \em blocking
    //! continuous evolution.
    //! - The \em reach set is given by
    //!   \f[ R=\{ \phi(\xi(s),t);\;\tau(s)+t \mid (x,t)\in D\times[0,h]
    //!       \mid \rho(s)\in C \ \wedge\ t\leq\delta(s)
    //!       \ \wedge\ \tau(s)+t\leq t_{\max}
    //!          \ \wedge\ \forall e\in E_{\mathrm{blk}},
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
    _apply_evolution_step(EvolutionData& evolution_data,
                          HybridEnclosure const& starting_set,
                          IntervalVectorFunctionModel const& flow,
                          TimingData const& timing_data,
                          Map<DiscreteEvent,CrossingData> const& crossing_data,
                          RealVectorFunction const& dynamic,
                          Map<DiscreteEvent,TransitionData> const& transitions) const;

    //! \brief Output a one-line summary of the current evolution state to the logging stream.
    virtual
    void
    _log_summary(EvolutionData const& evolution_data, HybridEnclosure const& starting_set) const;

  protected:
    void _create(const SystemType& system, FunctionFactoryType* factory);
    boost::shared_ptr< IntegratorInterface > _integrator_ptr;
  private:
    boost::shared_ptr< FunctionFactoryType > _function_factory_ptr;

    boost::shared_ptr< SolverInterface > _solver_ptr;
    boost::shared_ptr< SystemType > _sys_ptr;
    //boost::shared_ptr< EvolutionProfiler >  _profiler;
  protected:
    boost::shared_ptr< ConfigurationType > _configuration_ptr;
};


bool is_blocking(EventKind evk);
bool is_activating(EventKind evk);

//! \brief A data type used to store information about a transition of a hybrid system.
//! \relates HybridEvolverInterface
struct TransitionData
{
    //! \brief The event label.
    DiscreteEvent event;
    //! \brief The kind of event (invariant, urgent, etc) determining what
    //! happens when the event is active.
    EventKind event_kind;
    //! \brief The guard function of the event, the event being active when
    //! \f$g(x)\geq0\f$.
    RealScalarFunction guard_function;
    //! \brief The Lie derivative of the guard function along the flow.
    RealScalarFunction guard_flow_derivative_function;
    //! \brief The target location of the transition.
    DiscreteLocation target;
    //! \brief The reset function \f$x'=r(x)\f$ of the transition.
    RealVectorFunction reset_function;
    //! \brief The state space in the target location.
    RealSpace target_space;
    //TransitionData() { }
    //TransitionData(DiscreteLocation t, IntervalScalarFunction g, IntervalVectorFunction r)
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
    DECREASING_CROSSING, //!< The guard function is strictly decreasing along flow lines.
    CONCAVE_CROSSING, //!< The guard function is positive over at most an interval.
        //! Implied by concavity along flow lines, which is equivalent to \f$L_{f}^{2} g < 0\f$ within the reached set.
    CONVEX_CROSSING, //!< The guard function is negative over at most an interval. Implied by convexity along flow lines.
        //! Implied by convexity along flow lines, which is equivalent to \f$L_{f}^{2} g > 0\f$ within the reached set.
    TRANSVERSE_CROSSING, //!< The guard function is strictly increasing along flow lines, and
        //! the crossing time can be computed as a smooth function of the initial state.
        //! Given an initial point \f$x_0\f$, the crossing time is \f$\gamma(x_0)\f$.
    GRAZING_CROSSING //!< The time at which the guard function reaches a maximum along flow lines is a smooth function \f$\mu(x_0)\f$ of the initial state \f$x_0\f$.
        //! Implies by concavity along flow lines, which is equivalent to \f$L_{f}^{2} g < 0\f$ within the reached set.
};
std::ostream& operator<<(std::ostream& os, const CrossingKind& crk);

//! \brief A data type used to store information about the way flow lines cross a guard \f$g(x)=0\f$.
//! \relates HybridEvolverInterface
struct CrossingData
{
    CrossingData() : crossing_kind() { }
    CrossingData(CrossingKind crk) : crossing_kind(crk) { }
    CrossingData(CrossingKind crk, const IntervalScalarFunctionModel& crt)
        : crossing_kind(crk), crossing_time(crt) { }
    //! \brief The way in which the guard function changes along trajectories
    //! during a crossing. e.g. INCREASING_CROSSING
    CrossingKind crossing_kind;
    //! \brief The time \f$\gamma(x)\f$ at which the crossing occurs,
    //! as a function of the initial point in space. Satisfies \f$g(\phi(x,\gamma(x)))=0\f$.
    IntervalScalarFunctionModel crossing_time;
    //! \brief The time \f$\mu(x)\f$ at which the guard function reaches a maximum or minimum
    //! i.e. \f$L_{f}g(\phi(x,\mu(x))) = 0\f$.
    IntervalScalarFunctionModel critical_time;
};
std::ostream& operator<<(std::ostream& os, const CrossingData& crk);

//! \brief The kind of step taken in the evolution. Determines how the evolution time is specified.
//! \details Assumes that the starting set and time is given as a subset of \f$\{ \xi(s); \tau(s) \mid s\in D\}\f$
//! where \f$s\in D\f$ is a parameter, \f$x=\xi(s)\f$ is the state corresponding to parameter \f$s\f$, and \f$t=\tau(s)\f$
//! is the time the point has so far been evolved for. Assumes that the flow is given by a function \f$x'=\phi(x,t)\f$,
//! typically only defined over a bounded set of space and time.
//! \relates TimingData
enum StepKind {
    CONSTANT_EVOLUTION_TIME, //!< The step is taken for a fixed time \f$h\f$. The actual step length depends only on the starting state.
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),h)\f$ and \f$\tau'(s)=\tau(s)+h\f$.
    SPACE_DEPENDENT_EVOLUTION_TIME, //!< The step is taken for a time \f$\varepsilon(x)\f$ depending only on the starting state.
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),\varepsilon(\xi(s)))\f$ and \f$\tau'(s)=\tau(s)+\varepsilon(\xi(s))\f$.
    TIME_DEPENDENT_EVOLUTION_TIME, //!< The step is taken for a time \f$\varepsilon(t)\f$ depending only on the starting time.
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),\varepsilon(\tau(s)))\f$ and \f$\tau'(s)=\tau(s)+\varepsilon(\tau(s))\f$.
    SPACETIME_DEPENDENT_EVOLUTION_TIME, //!< The step is taken for a time \f$\varepsilon(x,t)\f$ depending only on the current state (space and time).
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),\varepsilon(\xi(s),\tau(s)))\f$ and \f$\tau'(s)=\tau(s)+\varepsilon(\xi(s),\tau(s))\f$.
    PARAMETER_DEPENDENT_EVOLUTION_TIME, //!< The step is taken for a time \f$\delta(s)\f$ depending on the parameterisation of the set.
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),\delta(s))\f$ and \f$\tau'(s)=\tau(s)+\delta(s)\f$.
    PARAMETER_DEPENDENT_FINISHING_TIME, //!< The step is taken up to a time \f$\omega(s)\f$ depending on the parameterisation of the starting set.
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),\omega(s)-\tau(s))\f$ and \f$\tau'(s)=\omega(s)\f$.
    SPACETIME_DEPENDENT_FINISHING_TIME, //!< The step is taken up to a time \f$\omega(x,t)\f$ depending only on the current state (space and time).
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),\omega(\xi(s),\tau(s)-\tau(s)))\f$ and \f$\tau'(s)=\omega(\xi(s),\tau(s))\f$.
    CONSTANT_FINISHING_TIME, //!< The step is taken up to the specified evolution time \f$t_{\max}\f$. The actual step length depends on the parameterisation.
      //! After the step, we have \f$\xi'(s) = \phi(\xi(s),t_{\max}-\tau(s))\f$ and \f$\tau'(s)=t_{\max}\f$.
};
std::ostream& operator<<(std::ostream& os, const StepKind& crk);

//! \brief The way in which the step interacts with the final evolution time.
//! \details Used to control whether final time constraints need to be added,
//! whether further evolution is possible, and whether and how sets should be put in the list of final sets.
//! Needed since the final time may be an arbitrary real number, at it may not be possible to determine
//! whether a given enclosure is exactly at the final time or not
//! \relates TimingData
enum FinishingKind {
    BEFORE_FINAL_TIME, //!< At the end of the step, the final time has definitely not been reached by any point.
    AT_FINAL_TIME, //!< At the end of the step, the final time is reached exactly. No more evolution is possible.
    AFTER_FINAL_TIME, //!< At the end of the step, the final time has definitely been passed by every point. No more evolution is possible.
    STRADDLE_FINAL_TIME, //!< At the end of the step, some points may have passed the final time, and some may not reached it.
        //!< Usable as a default value.
};
std::ostream& operator<<(std::ostream& os, const FinishingKind& crk);

//! \brief A data type used to store information about the kind of time step taken during hybrid evolution.
//! \relates HybridEvolverInterface
struct TimingData
{
    StepKind step_kind; //!< The kind of step taken in the evolution
    FinishingKind finishing_kind; //!< The relationship between the finishing time of the step, and the final time of the evolution trace.
    Real final_time; //!< The time \f$t_{\max}\f$ specified as the final time of the evolution trace.
    Float step_size; //!< The maximum step size \f$h\f$ allowed by the computed flow function.
    IntervalScalarFunctionModel spacetime_dependent_evolution_time;
        //!< The evolution time \f$\varepsilon(x,t)\f$ used in a \a SPACETIME_DEPENDENT_EVOLUTION_TIME step.
    IntervalScalarFunctionModel spacetime_dependent_finishing_time;
        //!< The final time \f$\omega(x,t)\f$ used in a \a SPACETIME_DEPENDENT_FINISHING_TIME step.
    IntervalScalarFunctionModel parameter_dependent_finishing_time;
        //!< The time \f$\omega(s)\f$ reached after an \a PARAMETER_DEPENDENT_FINISHING_TIME as a function of the parameters.
    IntervalScalarFunctionModel parameter_dependent_evolution_time;
        //!< The time \f$\delta(s)\f$ used in a \a PARAMETER_DEPENDENT_EVOLUTION_TIME step.
        //! Set equal to \f$\varepsilon(\xi(s))\f$ for a \a SPACE_DEPENDENT_EVOLUTION_TIME
        //! and \f$\omega(s)-\varepsilon(s)\f$ for an \a PARAMETER_DEPENDENT_FINISHING_TIME.
    Interval evolution_time_domain; //!< The time domain of the flow function, equal to \f$[0,h]\f$.
    IntervalScalarFunctionModel evolution_time_coordinate; //!< The time coordinate of the flow function, equal to the identity on \f$[0,h]\f$.
};

//! \brief A data type used to store information about the kind of time step taken during hybrid evolution.
//! \relates HybridEvolverInterface
struct EvolutionData
{
    //! \brief Sets for which the evolution starts in a new mode, initially or after a jump.
    //! The set needs to be processed by finding initially active events,
    //! and moving the set into working_sets.
    List<HybridEnclosure> initial_sets;
    //! \brief Sets which have been processed for initially active events.
    //! Sets in this list can be assumed to satisfy all invariants and progress
    //! predicates.
    List<HybridEnclosure> working_sets;

    //! \brief Sets which have been computed as reach sets for the current
    //! evolution.
    List<HybridEnclosure> reach_sets;
    //! \brief Sets which have been computed as final sets (i.e. satisfying
    //! either the final time or the maximum number of stesp).
    List<HybridEnclosure> final_sets;
    //! \brief Intermediate sets reached after each time step. Not relevant for
    //! the result, but useful for plotting, especially for debugging.
    List<HybridEnclosure> intermediate_sets;

    //! \brief The semantics used to compute the evolution. Defaults to UPPER_SEMANTICS.
    Semantics semantics;
};


//! \brief Configuration for HybridEvolverBase, for controlling the accuracy of continuous evolution methods.
class HybridEvolverBaseConfiguration : public ConfigurationInterface
{
  public:
    typedef uint UnsignedIntType;
    typedef double RealType;

  protected:

    HybridEvolverBaseConfiguration(HybridEvolverBase& evolver);

  private:

    HybridEvolverBase& _evolver;

  private:

    //! \brief The accuracy required for integration.
    //! Decreasing this value increases the accuracy of the computation.
    RealType _flow_accuracy;

    //! \brief The maximum allowable step size for integration.
    //! Decreasing this value increases the accuracy of the computation.
    RealType _maximum_step_size;

    //! \brief The maximum allowable radius of a basic set during integration.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    //! Decreasing this value reduces the actual evolution time of a lower-approximation.
    RealType _maximum_enclosure_radius;

    //! \brief The maximum allowable approximation error in the parameter-to-space mapping of an enclosure set.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType _maximum_spacial_error;

    //! \brief Enable reconditioning of basic sets (false by default).
    bool _enable_reconditioning;

    //! \brief Enable subdivision of basic sets (false by default), only for upper semantics.
    bool _enable_subdivisions;

  public:

    const RealType& flow_accuracy() const { return _flow_accuracy; }
    //! \brief Construct the _integrator of the evolver, then set the _flow_accuracy.
    void set_flow_accuracy(const RealType value);

    const RealType& maximum_step_size() const { return _maximum_step_size; }
    void set_maximum_step_size(const RealType value) { _maximum_step_size = value; }

    const RealType& maximum_enclosure_radius() const { return _maximum_enclosure_radius; }
    void set_maximum_enclosure_radius(const RealType value) { _maximum_enclosure_radius = value; }

    const RealType& maximum_spacial_error() const { return _maximum_spacial_error; }
    void set_maximum_spacial_error(const RealType value) { _maximum_spacial_error = value; }

    const bool& enable_reconditioning() const { return _enable_reconditioning; }
    void set_enable_reconditioning(const bool value) { _enable_reconditioning = value; }

    const bool& enable_subdivisions() const { return _enable_subdivisions; }
    void set_enable_subdivisions(const bool value) { _enable_subdivisions = value; }

  public:

    virtual std::ostream& write(std::ostream& os) const;
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

    GeneralHybridEvolver(const SystemType& system);
    GeneralHybridEvolver(const SystemType& system,
                         const IntervalFunctionModelFactoryInterface& factory);
    GeneralHybridEvolver* clone() const { return new GeneralHybridEvolver(*this); }

  protected:
    virtual
    TimingData
    _estimate_timing(Set<DiscreteEvent>& active_events,
                     Real final_time,
                     IntervalVectorFunctionModel const& flow,
                     Map<DiscreteEvent,CrossingData>& crossings,
                     Map<DiscreteEvent,TransitionData> const& transitions,
                     HybridEnclosure const& initial_set) const;

};


//! \brief Configuration for GeneralHybridEvolver, for controlling the accuracy of continuous evolution methods.
class GeneralHybridEvolverConfiguration : public HybridEvolverBaseConfiguration
{
  public:

    GeneralHybridEvolverConfiguration(GeneralHybridEvolver& evolver);

    virtual ~GeneralHybridEvolverConfiguration() { }
};

//! \brief Factory for GeneralHybridEvolver objects.
class GeneralHybridEvolverFactory
    : public HybridEvolverFactoryInterface
{
  private:

    boost::shared_ptr<IntervalFunctionModelFactoryInterface> _function_factory;

  public:

    GeneralHybridEvolverFactory();

    GeneralHybridEvolverFactory(const IntervalFunctionModelFactoryInterface& factory);

    virtual GeneralHybridEvolverFactory* clone() const { return new GeneralHybridEvolverFactory(*this); }

    //! \brief Create an hybrid evolver around a \a system.
    virtual GeneralHybridEvolver* create(const HybridAutomatonInterface& system) const;
};


} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_H
