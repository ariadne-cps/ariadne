/***************************************************************************
 *            hybrid_evolver.hpp
 *
 *  Copyright  2009  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file hybrid_evolver.hpp
 *  \brief Hybrid evolver classes.
 */

#ifndef ARIADNE_HYBRID_EVOLVER_HPP
#define ARIADNE_HYBRID_EVOLVER_HPP

#include <string>
#include <vector>
#include <list>
#include <iostream>


#include "../utility/tuple.hpp"

#include "../hybrid/hybrid_time.hpp"
#include "../hybrid/hybrid_set.hpp"

#include "../solvers/configuration_interface.hpp"
#include "../hybrid/hybrid_enclosure.hpp"
#include "../hybrid/hybrid_orbit.hpp"
#include "../hybrid/hybrid_automaton_interface.hpp"
#include "../hybrid/hybrid_evolver_interface.hpp"

#include "../output/logging.hpp"

namespace Ariadne {

typedef Map< DiscreteLocation, Vector<FloatDPValue> > HybridExactFloatVector;
typedef Dyadic StepSizeType;

class IntegratorInterface;
class SolverInterface;
class HybridEvolverBaseConfiguration;
class GeneralHybridEvolverConfiguration;

//! \ingroup FunctionModule
//! \brief A class representing the flow \f$\phi(x,t)\f$ of a differential equation \f$\frac{dx}{dt}=f(x)\f$.
class FlowFunctionModel
    : public ValidatedVectorFunctionModelDP
{
  public:
    FlowFunctionModel(const ValidatedVectorFunctionModelDP& f) : ValidatedVectorFunctionModelDP(f) { }
    StepSizeType step_size() const { return static_cast<StepSizeType>(this->time_domain().upper()); }
    ExactIntervalType time_domain() const { return this->domain()[this->domain().size()-1]; }
    ExactBoxType space_domain() const { return ExactBoxType(project(this->domain(),Ariadne::range(0,this->domain().size()-1))); }
    ExactBoxType const codomain() const { return this->ValidatedVectorFunctionModelDP::codomain(); }
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
    typedef ValidatedFunctionModelDPFactoryInterface FunctionFactoryType;
    typedef HybridAutomatonInterface SystemType;
    typedef SystemType::TimeType TimeType;
    typedef TimeType::ContinuousTimeType ContinuousTimeType;
    typedef List<DiscreteEvent> EventListType;
    typedef HybridEnclosure EnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<HybridEnclosure> EnclosureListType;

  public:

    //! \brief Default constructor.
    HybridEvolverBase(const SystemType& system);

    //! \brief Construct from parameters using a default integrator.
    HybridEvolverBase(
    		const SystemType& system,
            const FunctionFactoryType& factory);

    //! \brief Make a dynamically-allocated copy.
    HybridEvolverBase* clone() const = 0;

    //! \brief Get the internal system.
    const SystemType& system() const;


    //@{
    //! \name Settng and configuration for the class.

    //! \brief A reference to the evolver's configuration parameters.
    ConfigurationType& configuration();
    const ConfigurationType& configuration() const;

    //! \brief Change the configuration from a \a domain and \a lengths (NOT IMPLEMENTED).
    virtual Void reconfigure(const HybridExactBoxes& domain, const HybridExactFloatVector& lengths) { }

    //! \brief The class which constructs functions for the enclosures.
    const FunctionFactoryType& function_factory() const;
    //! \brief Set the class which constructs functions for the enclosures.
    Void set_function_factory(const FunctionFactoryType& factory);

    //! \brief Set the class which integrates the continuous dynamics.
    Void set_integrator(const IntegratorInterface& integrator);
    //! \brief Set the class which integrates the continuous dynamics.
    Void set_solver(const SolverInterface& solver);

    Bool ALLOW_CREEP; //!< If true, a less-than-full evolution step may be taken to avoid splitting due to partially crossing a guard.
    Bool ALLOW_UNWIND; //!< If true, a less-than-full evolution step may be taken to try to restore all time values over the parameter domain to the same value.
    //@}

    //@{
    //! \name Main evolution functions.

    Orbit<EnclosureType> orbit(const HybridExactBoxType& initial_box, const TerminationType& termination, Semantics semantics=Semantics::UPPER) const;
    Orbit<EnclosureType> orbit(const HybridBoxSet& initial_box, const TerminationType& termination, Semantics semantics=Semantics::UPPER) const;
    Orbit<EnclosureType> orbit(const HybridBoundedConstraintSet& initial_set, const TerminationType& termination, Semantics semantics=Semantics::UPPER) const;

    //! \brief Compute an approximation to the orbit set using the given semantics, starting from an initial enclosure.
    Orbit<EnclosureType> orbit(const EnclosureType& initial_enclosure, const TerminationType& termination, Semantics semantics=Semantics::UPPER) const;

    //! \brief Compute an approximation to the evolution set using the given semantics.
    EnclosureListType evolve(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics=Semantics::UPPER) const;

    //! \brief Compute an approximation to the evolution set under the given semantics.
    EnclosureListType reach(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics=Semantics::UPPER) const;

    //! \brief Compute an approximation to the evolution set under the given semantics.
    Pair<EnclosureListType,EnclosureListType> reach_evolve(const EnclosureType& initial_set, const TerminationType& termination, Semantics semantics=Semantics::UPPER) const;
    //@}

    //@{
    //! \name Auxiliary set conversion functionality

    //! \brief Set construct an enclosure from a box, such as one obtained from a grid.
    virtual EnclosureType enclosure(const HybridExactBox& initial_box) const;
    //! \brief Set construct an enclosure from a user-provided set.
    virtual EnclosureType enclosure(const HybridBoundedConstraintSet& initial_set) const;

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
    //!   Either \a #Semantics::LOWER, in which case trajectories terminate at
    //!   discontinuities, or #Semantics::UPPER, in which case all branches
    //!   are taken.
    //! \param reach A flag indicating whether the reachable sets should
    //!   be computed.
    virtual Void _evolution(ListSet<HybridEnclosure>& final, ListSet<HybridEnclosure>& reachable, ListSet<HybridEnclosure>& intermediate,
                            const EnclosureType& initial, const TerminationType& termination,
                            Semantics semantics, Bool reach) const;

    //! \brief Compute the evolution within a single discrete mode.
    //!
    //! Takes a set from evolution_data.starting_sets and processes the active
    //! events using _process_starting_sets.
    //! The resulting set is then evolved using _upper_evolution_step until
    //! either the maximum continuous time is reached, or no furter continuous
    //! evolution is possible.
    virtual
    Void
    _evolution_in_mode(EvolutionData& evolution_data,
                       TerminationType const& maximum_time) const;

    //! \brief Performs an evolution step on one of the sets listed in \a
    //! evolution_data.initial_sets.
    virtual
    Void
    _evolution_step(EvolutionData& evolution_data,
                    EffectiveVectorFunction const& dynamic,
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
    ValidatedVectorFunctionModelDP
    _compute_flow(EffectiveVectorFunction vector_field,
                  ExactBoxType const& initial_set,
                  const StepSizeType& maximum_step_size) const;

    //! \brief Compute the active events for the \a flow \f$\phi\f$ with
    //! time step \f$h\f$ starting in the given \a starting_set \f$S\f$.
    //! A event \f$e\f$ is \em active if \f$g_e(\phi(x,t))\geq0\f$
    //! for some \f$x\in S\f$ and \f$t\in[0,h]\f$.
    //! The result is allowed to contain events which are not active
    //! during the time step, since the exact set may be too difficult to
    //! compute.
    virtual
    Set<DiscreteEvent>
    _compute_active_events(EffectiveVectorFunction const& dynamic,
                           Map<DiscreteEvent,EffectiveScalarFunction> const& guards,
                           ValidatedVectorFunctionModelDP const& flow,
                           HybridEnclosure const& starting_set) const;

    //! \brief Compute data on how trajectories of the \a flow
    //! cross the given \a guards.
    //! For example, trajectories may have a #INCREASING crossing with
    //! a given crossing_time as a function of the initial state.
    //!
    //! \details The \a dynamic is given explicitly to facilitate computation
    //! of crossing times.
    //! A crossing_time computed must be valid over the \a initial_set, though
    //! will often be computed over the entire bounding box.
    //!
    //! If the crossing is #INCREASING, then the \a crossing_time \f$\gamma\f$
    //! should be computed, satisfying \f$g(\phi(x,\gamma(x)))=0\f$.
    //! If the crossing is #CONCAVE, the the \a critical_time (\a maximum_time) \f$\mu\f$
    //! should be computed, satisfying \f$L_{f}g(\phi(x,\mu(x)))=0\f$.
    virtual
    Map<DiscreteEvent,CrossingData>
    _compute_crossings(Set<DiscreteEvent> const& active_events,
                       EffectiveVectorFunction const& dynamic,
                       Map<DiscreteEvent,EffectiveScalarFunction> const& guards,
                       FlowFunctionModel const& flow,
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
                     FlowFunctionModel const& flow,
                     Map<DiscreteEvent,CrossingData>& crossings,
                     Map<DiscreteEvent,TransitionData> const& transitions,
                     HybridEnclosure const& initial_set) const;


    //! \brief Simplify the description of a set or allow for reconditioning
    //! of numerical error.
    virtual
    Void
    _recondition(HybridEnclosure& set) const;

    //! \brief Process the \a starting_set to find any
    //!  <em>immediately active</em> events, and compute the relevant \em jump sets
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
    Void
    _process_starting_events(EvolutionData& evolution_data,
                            HybridEnclosure const& starting_set,
                            Map<DiscreteEvent,TransitionData> const& transitions) const;

    //! \brief Apply the \a flow to the \a set up to the time specified by \a timing_data
    //! to obtain the single-step unconstrained reachable set.
    virtual
    Void
    _apply_reach_step(HybridEnclosure& set,
                      ValidatedVectorFunctionModelDP const& flow,
                      TimingData const& timing_data) const;

    //! \brief Apply the \a flow to the \a set for the time specified by \a timing_data
    //! to obtain the single-step unconstrained evolved set.
    virtual
    Void
    _apply_evolve_step(HybridEnclosure& set,
                       ValidatedVectorFunctionModelDP const& flow,
                       TimingData const& timing_data) const;

    //! \brief Apply the \a flow to the \a set for to reach the
    //! guard specified by \a transition_data using guidelines provided by \a crossing_data.
    //!
    //! \details
    //! NOTE: The \a dynamic is used to compute the Lie derivative of
    //! the flow to ensure positive crossing in the degenerate case. It should not really be necessary.
    virtual
    Void
    _apply_guard_step(HybridEnclosure& set,
                      EffectiveVectorFunction const& dynamic,
                      ValidatedVectorFunctionModelDP const& flow,
                      TimingData const& timing_data,
                      TransitionData const& transition_data,
                      CrossingData const& crossing_data,
                      const Semantics semantics) const;

    //! \brief Apply the invariants in \a transition_data to the set \a set.
    //! \callgraph
    virtual
    Void
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
    Void
    _apply_guard(List<HybridEnclosure>& sets,
                 const ValidatedScalarFunctionModelDP& sets_elapsed_time,
                 const HybridEnclosure& starting_set,
                 const ValidatedVectorFunctionModelDP& flow,
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
    //! with #CONCAVE events, it may be the case that these sets are split into
    //! two.
    //!
    //! TODO: This method might be better split into even simpler events.
    virtual
    Void
    _apply_evolution_step(EvolutionData& evolution_data,
                          HybridEnclosure const& starting_set,
                          ValidatedVectorFunctionModelDP const& flow,
                          TimingData const& timing_data,
                          Map<DiscreteEvent,CrossingData> const& crossing_data,
                          EffectiveVectorFunction const& dynamic,
                          Map<DiscreteEvent,TransitionData> const& transitions) const;

    //! \brief Output a one-line summary of the current evolution state to the logging stream.
    virtual
    Void
    _log_summary(EvolutionData const& evolution_data, HybridEnclosure const& starting_set) const;

  protected:
    Void _create(const SystemType& system, FunctionFactoryType* factory);
    std::shared_ptr< IntegratorInterface > _integrator_ptr;
  private:
    std::shared_ptr< FunctionFactoryType > _function_factory_ptr;
  protected:
    std::shared_ptr< SolverInterface > _solver_ptr;
    std::shared_ptr< SystemType > _sys_ptr;
    //std::shared_ptr< EvolutionProfiler >  _profiler;
  protected:
    std::shared_ptr< ConfigurationType > _configuration_ptr;
};


Bool is_blocking(EventKind evk);
Bool is_activating(EventKind evk);

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
    EffectiveScalarFunction guard_function;
    //! \brief The Lie derivative of the guard function along the flow.
    EffectiveScalarFunction guard_flow_derivative_function;
    //! \brief The target location of the transition.
    DiscreteLocation target;
    //! \brief The reset function \f$x'=r(x)\f$ of the transition.
    EffectiveVectorFunction reset_function;
    //! \brief The state space in the target location.
    RealSpace target_space;
    //TransitionData() { }
    //TransitionData(DiscreteLocation t, ValidatedScalarFunction g, ValidatedVectorFunction r)
    //    : target(t), guard_function(g), reset_function(r) { }
};
OutputStream& operator<<(OutputStream& os, const TransitionData& transition);

//! \relates HybridEvolverBase
//! \brief Information on how a flow tube crosses a hypersurface.
enum class DirectionKind : std::uint8_t {
    POSITIVE, //!< The guard function is strictly positive on the flow range.
        //! The event occurs immediately (if urgent) or at any time (if permissive).
    NEGATIVE, //!< The guard function is strictly negative on the flow range. No event occurs.
    INCREASING, //!< The guard function is strictly increasing along flow lines.
        //! i.e. \f$\frac{d}{dt}g(\phi(x_0,t))>0\f$. Implied by \f$L_{f}g > 0\f$.
    DECREASING, //!< The guard function is strictly decreasing along flow lines.
        //! i.e. \f$\frac{d}{dt}g(\phi(x_0,t))<0\f$. Implied by \f$L_{f}g < 0\f$.
    CONCAVE, //!< The guard function varies concavely along flow lines, which is equivalent to \f$L_{f}^{2} g < 0\f$.
    CONVEX, //!< The guard function varies convexly along flow lines, which is equivalent to \f$L_{f}^{2} g > 0\f$.
    INDETERMINATE //!< Neither the guard function, nor its first or second Lie derivatives, has a definite sign.
};
OutputStream& operator<<(OutputStream& os, const DirectionKind& dir);

//! \relates HybridEvolverBase
//! \brief The way trajectories of the flow \f$\phi(x_0,t)\f$ of \f$\frac{dx}{dt}=f(x)\f$ cross the guard set \f$g(x)=0\f$.
//! The rate of change of the guard function changes along flow lines is given by
//! \f$d g(x(t))/dt = L_{f}g(x(t))\f$ where the Lie derivative \f$L_{f}g\f$ is defined by \f$L_{f}g(x) = (\nabla g\cdot f)(x)\f$.
//!
//! Defaults to DEGENERATE whenever the crossing kind has not been resolved;
//! this need not necessarily mean that the crossing is degenerate, but may
//! also imply that the crossing information is too expensive or sensitive to
//! compute.
//! \relates HybridEvolverInterface \relates CrossingData
enum class CrossingKind : std::uint8_t {
    DEGENERATE, //!< The crossing may be degenerate to second order.
    NEGATIVE, //!< The guard function is negative on the flow domain. No event occurs.
    POSITIVE, //!< The guard function is negative on the domain. The event occurs immediately (if urgent) or at all times (if permissive).
    MONOTONE, //!< The guard function is strictly increasing or decreasing along flow lines at a crossing point.
        //! i.e. \f$\frac{d}{dt}g(\phi(x_0,t))>0\f$ whenever \f$g(\phi(x_0,t))=0\f$.
        //! Implied by \f$L_{f}g \geq 0\f$.
    INCREASING, //!< The guard function is strictly increasing along flow lines at a crossing point.
        //! i.e. \f$\frac{d}{dt}g(\phi(x_0,t))>0\f$ whenever \f$g(\phi(x_0,t))=0\f$.
        //! Implied by \f$L_{f}g \geq 0\f$.
    DECREASING, //!< The guard function is strictly decreasing along flow lines.
    CONCAVE, //!< The guard function is positive over at most an interval.
        //! Implied by concavity along flow lines, which is equivalent to \f$L_{f}^{2} g < 0\f$ within the reached set.
    CONVEX, //!< The guard function is negative over at most an interval. Implied by convexity along flow lines.
        //! Implied by convexity along flow lines, which is equivalent to \f$L_{f}^{2} g > 0\f$ within the reached set.
    TRANSVERSE, //!< The guard function is strictly increasing along flow lines, and
        //! the crossing time can be computed as a smooth function of the initial state.
        //! Given an initial point \f$x_0\f$, the crossing time is \f$\gamma(x_0)\f$.
    GRAZING //!< The time at which the guard function reaches a maximum along flow lines is a smooth function \f$\mu(x_0)\f$ of the initial state \f$x_0\f$.
        //! Implies by concavity along flow lines, which is equivalent to \f$L_{f}^{2} g < 0\f$ within the reached set.
};
OutputStream& operator<<(OutputStream& os, const CrossingKind& crk);

//! \brief A data type used to store information about the way flow lines cross a guard \f$g(x)=0\f$.
//! \relates HybridEvolverInterface
struct CrossingData
{
    CrossingData() : crossing_kind() { }
    CrossingData(CrossingKind crk) : crossing_kind(crk) { }
    CrossingData(CrossingKind crk, const ValidatedScalarFunctionModelDP& crt)
        : crossing_kind(crk), crossing_time(crt) { }
    //! \brief The way in which the guard function changes along trajectories
    DirectionKind direction_kind;
    //! \brief The way in which the guard function changes along trajectories
    //! during a crossing. e.g. INCREASING
    CrossingKind crossing_kind;
    //! \brief The range of times at which the crossing may occur.
    ExactIntervalType crossing_time_range;
    //! \brief The time \f$\gamma(x)\f$ at which the crossing occurs,
    //! as a function of the initial point in space. Satisfies \f$g(\phi(x,\gamma(x)))=0\f$.
    ValidatedScalarFunctionModelDP crossing_time;
    //! \brief The time \f$\mu(x)\f$ at which the guard function reaches a maximum or minimum
    //! i.e. \f$L_{f}g(\phi(x,\mu(x))) = 0\f$.
    ValidatedScalarFunctionModelDP critical_time;
};
OutputStream& operator<<(OutputStream& os, const CrossingData& crk);

//! \brief The kind of step taken in the evolution. Determines how the evolution time is specified.
//! \details Assumes that the starting set and time is given as a subset of \f$\{ \xi(s); \tau(s) \mid s\in D\}\f$
//! where \f$s\in D\f$ is a parameter, \f$x=\xi(s)\f$ is the state corresponding to parameter \f$s\f$, and \f$t=\tau(s)\f$
//! is the time the point has so far been evolved for. Assumes that the flow is given by a function \f$x'=\phi(x,t)\f$,
//! typically only defined over a singleton set of space and time.
//! \relates TimingData
enum class StepKind : std::uint8_t {
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
OutputStream& operator<<(OutputStream& os, const StepKind& crk);

//! \brief The way in which the step interacts with the final evolution time.
//! \details Used to control whether final time constraints need to be added,
//! whether further evolution is possible, and whether and how sets should be put in the list of final sets.
//! Needed since the final time may be an arbitrary real number, at it may not be possible to determine
//! whether a given enclosure is exactly at the final time or not
//! \relates TimingData
enum class FinishingKind : std::uint8_t {
    BEFORE_FINAL_TIME, //!< At the end of the step, the final time has definitely not been reached by any point.
    AT_FINAL_TIME, //!< At the end of the step, the final time is reached exactly. No more evolution is possible.
    AFTER_FINAL_TIME, //!< At the end of the step, the final time has definitely been passed by every point. No more evolution is possible.
    STRADDLE_FINAL_TIME, //!< At the end of the step, some points may have passed the final time, and some may not reached it.
        //!< Usable as a default value.
};
OutputStream& operator<<(OutputStream& os, const FinishingKind& crk);

//! \brief A data type used to store information about the kind of time step taken during hybrid evolution.
//! \relates HybridEvolverInterface
struct TimingData
{
    StepKind step_kind; //!< The kind of step taken in the evolution
    FinishingKind finishing_kind; //!< The relationship between the finishing time of the step, and the final time of the evolution trace.
    Real final_time; //!< The time \f$t_{\max}\f$ specified as the final time of the evolution trace.
    FloatDPValue step_size; //!< The maximum step size \f$h\f$ allowed by the computed flow function.
    ValidatedScalarFunctionModelDP spacetime_dependent_evolution_time;
        //!< The evolution time \f$\varepsilon(x,t)\f$ used in a \a SPACETIME_DEPENDENT_EVOLUTION_TIME step.
    ValidatedScalarFunctionModelDP spacetime_dependent_finishing_time;
        //!< The final time \f$\omega(x,t)\f$ used in a \a SPACETIME_DEPENDENT_FINISHING_TIME step.
    ValidatedScalarFunctionModelDP parameter_dependent_finishing_time;
        //!< The time \f$\omega(s)\f$ reached after an \a PARAMETER_DEPENDENT_FINISHING_TIME as a function of the parameters.
    ValidatedScalarFunctionModelDP parameter_dependent_evolution_time;
        //!< The time \f$\delta(s)\f$ used in a \a PARAMETER_DEPENDENT_EVOLUTION_TIME step.
        //! Set equal to \f$\varepsilon(\xi(s))\f$ for a \a SPACE_DEPENDENT_EVOLUTION_TIME
        //! and \f$\omega(s)-\varepsilon(s)\f$ for an \a PARAMETER_DEPENDENT_FINISHING_TIME.
    ExactIntervalType evolution_time_domain; //!< The time domain of the flow function, equal to \f$[0,h]\f$.
    ValidatedScalarFunctionModelDP evolution_time_coordinate; //!< The time coordinate of the flow function, equal to the identity on \f$[0,h]\f$.
};
OutputStream& operator<<(OutputStream& os, const TimingData& timing);

//! \brief A data type used to store information about the kind of time step taken during hybrid evolution.
//! \relates HybridEvolverInterface
//! \relates HybridEvolverBase
struct EvolutionData
{
    //! \brief Sets for which the evolution starts in a new mode, initially or after a jump.
    //! The set needs to be processed by finding initially active events,
    //! and moving the set into working_sets.
    List<HybridEnclosure> starting_sets;
    //! \brief Sets which have been processed for initially active events.
    //! Sets in this list can be assumed to satisfy all invariants and progress
    //! predicates.
    List<HybridEnclosure> working_sets;

    //! \brief The initial sets for the evolution.
    List<HybridEnclosure> initial_sets;
    //! \brief Sets which have been computed as reach sets for the current
    //! evolution.
    List<HybridEnclosure> reach_sets;
    //! \brief Sets which have been computed as final sets (i.e. satisfying
    //! either the final time or the maximum number of stesp).
    List<HybridEnclosure> final_sets;
    //! \brief Intermediate sets reached after each time step. Not relevant for
    //! the result, but useful for plotting, especially for debugging.
    List<HybridEnclosure> intermediate_sets;

    //! \brief The semantics used to compute the evolution. Defaults to Semantics::UPPER.
    Semantics semantics;
};


//! \relates HybridEvolverBase
//! \brief Configuration for HybridEvolverBase, for controlling the accuracy of continuous evolution methods.
class HybridEvolverBaseConfiguration : public ConfigurationInterface
{
  public:
    typedef Nat UnsignedIntType;
    typedef FloatDPValue RealType;
    typedef double RawRealType;

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
    StepSizeType _maximum_step_size;

    //! \brief The maximum allowable radius of a basic set during integration.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    //! Decreasing this value reduces the actual evolution time of a lower-approximation.
    RealType _maximum_enclosure_radius;

    //! \brief The maximum allowable approximation error in the parameter-to-space mapping of an enclosure set.
    //! Decreasing this value increases the accuracy of the computation of an over-approximation.
    RealType _maximum_spacial_error;

    //! \brief Enable reconditioning of basic sets (false by default).
    Bool _enable_reconditioning;

    //! \brief Enable subdivision of basic sets (false by default), only for upper semantics.
    Bool _enable_subdivisions;

  public:

    const RealType& flow_accuracy() const { return _flow_accuracy; }
    //! \brief Construct the _integrator of the evolver, then set the _flow_accuracy.
    Void set_flow_accuracy(const RawRealType value);

    const StepSizeType& maximum_step_size() const { return _maximum_step_size; }
    Void set_maximum_step_size(const StepSizeType value) { _maximum_step_size = value; }
    Void set_maximum_step_size(const double value) { _maximum_step_size = static_cast<StepSizeType>(value); }

    const RealType& maximum_enclosure_radius() const { return _maximum_enclosure_radius; }
    Void set_maximum_enclosure_radius(const RawRealType value) { _maximum_enclosure_radius = RealType(value); }

    const RealType& maximum_spacial_error() const { return _maximum_spacial_error; }
    Void set_maximum_spacial_error(const RawRealType value) { _maximum_spacial_error = RealType(value); }

    const Bool& enable_reconditioning() const { return _enable_reconditioning; }
    Void set_enable_reconditioning(const Bool value) { _enable_reconditioning = value; }

    const Bool& enable_subdivisions() const { return _enable_subdivisions; }
    Void set_enable_subdivisions(const Bool value) { _enable_subdivisions = value; }

  public:

    virtual OutputStream& write(OutputStream& os) const;
};


//! \ingroup AnalysisModule
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
                         const ValidatedFunctionModelDPFactoryInterface& factory);
    virtual GeneralHybridEvolver* clone() const { return new GeneralHybridEvolver(*this); }
    virtual OutputStream& write(OutputStream& os) const { return os << "GeneralHybridEvolver( " << this->configuration() << ")"; }

  protected:
    virtual
    TimingData
    _estimate_timing(Set<DiscreteEvent>& active_events,
                     Real final_time,
                     FlowFunctionModel const& flow,
                     Map<DiscreteEvent,CrossingData>& crossings,
                     Map<DiscreteEvent,TransitionData> const& transitions,
                     HybridEnclosure const& initial_set) const;

/*
    virtual
    Void
    _compute_finishing_time(TimingData& result,
                            Real final_time,
                            FlowFunctionModel const& flow,
                            HybridEnclosure const& initial_set) const;

    virtual
    Void
    _compute_unwind_time(TimingData& result,
                         Real final_time,
                         FlowFunctionModel const& flow,
                         HybridEnclosure const& initial_set) const;
    virtual
    Void
    _compute_creep_time(TimingData& result,
                        Real final_time,
                        FlowFunctionModel const& flow,
                        Map<DiscreteEvent,CrossingData>& crossings,
                        Map<DiscreteEvent,TransitionData> const& transitions,
                        HybridEnclosure const& initial_set) const;
*/
};


//! \brief Configuration for GeneralHybridEvolver, for controlling the accuracy of continuous evolution methods.
//! \relates GeneralHybridEvolver
class GeneralHybridEvolverConfiguration : public HybridEvolverBaseConfiguration
{
  public:

    GeneralHybridEvolverConfiguration(GeneralHybridEvolver& evolver);

    virtual ~GeneralHybridEvolverConfiguration() = default;
};

//! \brief Factory for GeneralHybridEvolver objects.
//! \deprecated Use cloning instead
class GeneralHybridEvolverFactory
    : public HybridEvolverFactoryInterface
{
  private:

    std::shared_ptr<ValidatedFunctionModelDPFactoryInterface> _function_factory;

  public:

    GeneralHybridEvolverFactory();

    GeneralHybridEvolverFactory(const ValidatedFunctionModelDPFactoryInterface& factory);

    virtual GeneralHybridEvolverFactory* clone() const { return new GeneralHybridEvolverFactory(*this); }

    //! \brief Create an hybrid evolver around a \a system.
    virtual GeneralHybridEvolver* create(const HybridAutomatonInterface& system) const;
};


} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_HPP
