/***************************************************************************
 *            hybrid_enclosure.hpp
 *
 *  Copyright  2009-10  Pieter Collins
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

/*! \file hybrid_enclosure.hpp
 *  \brief Enclosure sets for hybrid systems
 */

#ifndef ARIADNE_HYBRID_ENCLOSURE_HPP
#define ARIADNE_HYBRID_ENCLOSURE_HPP

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "../output/logging.hpp"
#include "../utility/declarations.hpp"
#include "../utility/pointer.hpp"
#include "../utility/container.hpp"

#include "../hybrid/hybrid_set.decl.hpp"
#include "../hybrid/discrete_location.hpp"
#include "../hybrid/discrete_event.hpp"
#include "../hybrid/hybrid_graphics_interface.hpp"

#include "../geometry/box.hpp"
#include "../geometry/enclosure.hpp"

namespace Ariadne {

template<class X> class Vector;
template<class X> struct LinearProgram;

class Enclosure;
class Grid;
class GridTreePaving;
class AffineSet;
class DiscreteEvent;
class Figure;
class CanvasInterface;

template<class L, class R> class Assignment;
typedef Assignment<RealVariable,RealExpression> RealAssignment;
template<class X> class Space;
typedef Space<Real> RealSpace;
template<class ES> class ListSet;
template<class ES> class HybridListSet;

class HybridEnclosure;
template<> class ListSet<HybridEnclosure>;

enum class EnclosureVariableType : std::uint8_t { INITIAL, TEMPORAL, PARAMETER, INPUT, NOISE, ERROR, UNKNOWN };

//! \brief A class representing an enclosure for a hybrid evolution.
//! Handles progress, activation and guard constraints internally.
//! The set is represented as the image of a box \f$D\f$ under a function model \f$\hat{f}(s)\f$, under the constraints
//! \f$\hat{c}(s) \leq 0\f$ and \f$\hat{e}(s)=0\f$. Also keeps track of the current time \f$\hat{t}(s)\f$.
//!
//! In other words,
//! \f[ S=\{ \hat{f}(s);\  \hat{t}(s) \mid s\in D \mid \hat{c}(s) \leq 0 \ \wedge \hat{e}(s)=0 \} . \f]
//! In the following documentation, we sometimes write \f$\xi(s)\f$ for \f$\hat{f}(s)\f$ and \f$\tau(s)\f$ for \f$\hat{t}(s)\f$.
//! yielding \f[ S=\{ \xi(s);\  \tau(s) \mid s\in D \mid \rho(s) \in C \} . \f]
//!
//! <b>Rationale:</b><\par>
//! <b>Parameterisation of the time specifier</b>
//! When computing a flow step, the evolution time can be given either as a
//! function \f$\varepsilon(x)\f$ of the state, or as an absolute final time
//! \f$\omega(s)\f$ of the parameters.
//! The rationale for this difference is that we often use the
//! \f$\varepsilon\f$ form when computing the flow up to a crossing with a
//! guard, in which case the time depends purely on the state, and it may
//! be possible to simplify the constraints by expressing them in terms of
//! the current states (in an implementation in which intermediate states are
//! considered explicitly).
//! The rationale for using a parameterised form with independent variable
//! \f$s\f$ is that this is needed to distinguished different states based on
//! crossings with the guard sets. The rationale for using a final time
//! \f$\omega(s)\f$ as a final time is that this allows us to remove
//! errors in the elapsed time \f$\tau(s)\f$.
//!
//! A possible extension would be to allow \f$\omega(s)\f$ to depend on
//! intermediate variables \f$x_i\f$ at previous stages. The reason for not
//! doing this is that it would enforce a policy for tracking intermediate
//! states, which we do not want to insist on during an implementation.
//!
//! <b>Replacing events</b>
//! If \f$g_e\f$ is increasing along trajectories
//! it is safe to remove a constraint \f$g_e(x_{i-1})\leq0\f$ and
//! replace it with the constraint \f$g_e(x_i)\leq 0\f$.
class HybridEnclosure
    : public HybridDrawableInterface
    , public Loggable
{
    friend class SimpleHybridEvolver;
    friend class ConstraintHybridEvolver;
  public:
    typedef Enclosure ContinuousStateSetType;
  private:
    DiscreteLocation _location;
    List<DiscreteEvent> _events;
    List<RealVariable> _state_space;
    List<RealVariable> _auxiliary_space;
    Enclosure _set;
    List<EnclosureVariableType> _variables;
  public:
    //! \brief An empty enclosure.
    HybridEnclosure();
    //! \brief An enclosure corresponding to a hybrid box \a hbx with variables canonically ordered by \a spc.
    HybridEnclosure(const HybridBoxSet& hbx, const RealSpace& spc, const ValidatedFunctionModelDPFactoryInterface& fac);
    //! \brief An enclosure corresponding to the hybrid set \a set using \a space to order the continuous variables.
    HybridEnclosure(const HybridBoundedConstraintSet& set, const RealSpace& space, const ValidatedFunctionModelDPFactoryInterface& factory);

    //! \brief An enclosure corresponding to a Euclidean box \a bx in location \a q with variables ordered by \a spc.
    HybridEnclosure(const DiscreteLocation& q, const RealSpace& spc, const RealBox& bx, const ValidatedFunctionModelDPFactoryInterface& fac);
    //! \brief An enclosure corresponding to a hybrid box \a hbx.
    explicit HybridEnclosure(const HybridRealBox& hbx, const ValidatedFunctionModelDPFactoryInterface& fac);

    //! \brief An enclosure corresponding to a hybrid box \a hbx.
    explicit HybridEnclosure(const HybridExactBoxType& hbx, const ValidatedFunctionModelDPFactoryInterface& fac);
    //! \brief An enclosure corresponding to a hybrid box \a hbx.
    explicit HybridEnclosure(const HybridExactBoxType& hbx, List<RealAssignment> aux, const ValidatedFunctionModelDPFactoryInterface& fac);
    //! \brief An enclosure constructed from a location \a q, a real space \a spc, and a (timed) enclosure \a es.
    explicit HybridEnclosure(const DiscreteLocation& q, const RealSpace& spc, const Enclosure& es);

    //! \brief Destructor.
    ~HybridEnclosure();
    //! \brief Create a dynamically-allocated copy.
    HybridEnclosure* clone() const;

    //! \brief The algorithms used to compute with the set.
    const EnclosureConfiguration& configuration() const;
    //! \brief The current location.
    const DiscreteLocation& location() const;
    //! \brief The Euclidean space of the location, including state, time and auxiliary functions.
    const RealSpace state_time_auxiliary_space() const;
    //! \brief The Euclidean state space of the location.
    const RealSpace state_space() const;
    //! \brief The global evolution time variable.
    const RealVariable time_variable() const;
    //! \brief The Euclidean state space of the location.
    const RealSpace auxiliary_space() const;
    //! \brief The factory used to create functions.
    const ValidatedFunctionModelDPFactoryInterface& function_factory() const;
    //! \brief The list of previous events.
    const List<DiscreteEvent>& previous_events() const;
    //! \brief The number of independent parameters.
    SizeType number_of_parameters() const;
    //! \brief The number of constraints.
    SizeType number_of_constraints() const;
    //! \brief The continuous state set.
    const ExactBoxType parameter_domain() const;
    //! \brief The function related to the state space.
    const ValidatedVectorFunctionModelDP& state_function() const;
    //! \brief The function related to time.
    const ValidatedScalarFunctionModelDP& time_function() const;
    //! \brief The function giving the time since the last event.
    const ValidatedScalarFunctionModelDP& dwell_time_function() const;
    //! \brief The function related to the auxiliary space.
    const ValidatedVectorFunctionModelDP auxiliary_function() const;
    //! \brief The function related to the auxiliary space.
    const ValidatedVectorFunctionModelDP state_time_auxiliary_function() const;
    //! \brief The function related to the variable \a var.
    const ValidatedScalarFunctionModelDP function(RealVariable var) const;

    //! \brief Set the evolution time function to \a omega.
    Void set_time_function(const ValidatedScalarFunctionModelDP& omega);

    //! \brief A bounding box for the space.
    UpperBoxType state_bounding_box() const;
    //! \brief The range of times since the starting time that the set represents.
    UpperIntervalType time_range() const;
    //! \brief The range of times since the last event.
    UpperIntervalType dwell_time_range() const;

    //! \brief The continuous state set.
    friend HybridBasicSet<Enclosure> project(HybridEnclosure const&, RealSpace const& spc);
    //! \brief The continuous state set.
    HybridBasicSet<Enclosure> state_set() const;
    //! \brief The continuous state set including time.
    HybridBasicSet<Enclosure> state_time_set() const;
    //! \brief The continuous state set including auxiliary variables.
    HybridBasicSet<Enclosure> state_auxiliary_set() const;
    //! \brief The continuous enclosure.
    const ContinuousStateSetType& continuous_set() const;

    //! \brief Set the time function to zero.
    Void clear_time();
    //! \brief Clears the list of previous events.
    Void clear_events();

    //! \brief Apply the reset map \a r corresponding to event \a e with target location \a q.
    //! Corresponds to replacing \f$\xi\f$ by \f$r\circ \xi\f$.
    Void apply_reset(DiscreteEvent e, DiscreteLocation q, RealSpace s, const ValidatedVectorFunction& r);
    //! \brief Apply the evolve step \xi'(s) = phi(xi(s),eps) and tau'(s)=tau(s)+eps
    Void apply_fixed_evolve_step(const ValidatedVectorFunctionModelDP& phi, const FloatDPValue& eps);
    //! \brief Apply the evolve step \xi'(s) = phi(xi(s),eps(xi(s))) and tau'(s)=tau(s)+eps(xi(s))
    Void apply_space_evolve_step(const ValidatedVectorFunctionModelDP& phi, const ValidatedScalarFunctionModelDP& eps);
    //! \brief Apply the evolve step \xi'(s) = phi(xi(s),eps(xi(s),tau(s))) and tau'(s)=tau(s)+eps(xi(s),tau(s))
    Void apply_spacetime_evolve_step(const ValidatedVectorFunctionModelDP& phi, const ValidatedScalarFunctionModelDP& eps);
    //! \brief Apply the reach step \xi'(s) = phi(xi(s),t-tau(s)) and tau'(s)=tau(s)+t for 0<=t<=eps(xi(s),tau(s))
    Void apply_spacetime_reach_step(const ValidatedVectorFunctionModelDP& phi, const ValidatedScalarFunctionModelDP& eps);
    // Compute the evolve step \xi'(s) = phi(xi(s),eps(s)) and tau'(s)=tau(s)+eps(s)
    Void apply_parameter_evolve_step(const ValidatedVectorFunctionModelDP& phi, const ValidatedScalarFunctionModelDP& eps);
    // Compute the evolve step \xi'(s) = phi(xi(s),\omega(s)-tau(s)) and tau'(s)=omega(s)
    Void apply_finishing_parameter_evolve_step(const ValidatedVectorFunctionModelDP& phi, const ValidatedScalarFunctionModelDP& omega);
    //! \brief Compute the reach step xi'(s,t) = phi(xi(s),t) and tau'(s,t)=tau(s)+t for t in [0,h] and t <= eps(s) , assuming eps(s)<= h throughout.
    Void apply_parameter_reach_step(const ValidatedVectorFunctionModelDP& phi, const ValidatedScalarFunctionModelDP& eps);
    //! \brief Compute the reach step xi'(s,t) = phi(xi(s),t) and tau'(s,t)=tau(s)+t for t in [0,h].
    Void apply_full_reach_step(const ValidatedVectorFunctionModelDP& phi);


    //! \brief Set the time of evolution to \a \f$t_{\max}\f$.
    //! Corresponds to introducting the constraint \f$\tau(s) = t_{\max}\f$.
    Void set_time(Real tmax);
    //! \brief Set the time of evolution to \a omega.
    //! Corresponds to introducting the constraint \f$\tau(s) = \omega(s)\f$.
    Void set_time(ValidatedScalarFunction omega);
    //! \brief Introduces the constraint \f$\tau(s)\leq \omega(s)\f$.
    Void bound_time(ValidatedScalarFunction omega);
    //! \brief Introduces the constraint \f$\tau(s)\leq t_{\max}\f$.
    Void bound_time(Real tmax);

    //! \brief Set the maximum time of evolution to \a \f$t_{\max}\f$. \deprecated
    //! Corresponds to introducting the constraint \f$\tau(s)\leq t_{\max}\f$.
    Void set_maximum_time(DiscreteEvent event, RawFloatDP tmax);
    //! \brief Set the current time-step to \f$h\f$. \deprecated
    Void set_step_time(FloatDPValue h);
    //! \brief \deprecated
    Void new_time_step_bound(DiscreteEvent e, ValidatedScalarFunction tau);

    //! \brief Sets the auxiliary variables and functions.
    Void set_auxiliary(List<RealVariable> vars, EffectiveVectorFunction aux);

    //! \brief Introduces a new parameter with domain \a ivl.
    Void new_parameter(ExactIntervalType ivl, EnclosureVariableType);
    //! \brief Introduce a new independent variable with domain \a ivl.
    Void new_variable(ExactIntervalType ivl, EnclosureVariableType);
    //! \brief Introduces a new state constraint \f$C\f$ on \f$x\f$. \deprecated
    Void new_constraint(DiscreteEvent e, ValidatedConstraint c);
    //! \brief Introduces a new state constraint \f$C\f$ on \f$x\f$.
    Void new_state_constraint(DiscreteEvent e, ValidatedConstraint c);
    //! \brief Introduces a new state constraint \f$C\f$ on \f$(x,t)\f$.
    Void new_state_time_constraint(DiscreteEvent e, ValidatedConstraint c);
    //! \brief Introduces a new constraint \f$C\f$ on \f$s\f$.
    Void new_parameter_constraint(DiscreteEvent e, ValidatedConstraint c);
    //! \brief Introduces the new constraint \f$t\leq\gamma(x)\f$.
    Void new_state_time_bound(DiscreteEvent e, ValidatedScalarFunction gamma);
    //! \brief Introduces the new invariant (progress predicate) \f$c(x)\leq0\f$.
    Void new_invariant(DiscreteEvent e, ValidatedScalarFunction c);
    //! \brief Introduces the new activation condition \f$g(x)\geq0\f$ for the event \a e.
    Void new_activation(DiscreteEvent e,ValidatedScalarFunction g);
    //! \brief Introduces the new guard condition \f$g(x)=0\f$ for the event \a e.
    //! More precisely, the continuous dynamics is restricted to \f$c(x)\leq0\f$, and the event happens when \f$c(x)\geq0\f$.
    Void new_guard(DiscreteEvent e, ValidatedScalarFunction g);
    //! \brief Introduces the new guard condition \f$g(x)=0\f$ for the event \a e, with computed crossing time \f$\tau(s)\f$.
    Void new_guard(DiscreteEvent e, ValidatedScalarFunction g, ValidatedScalarFunction ct);


    //! \brief The dimension of the set. Returns the state-space dimension.
    DimensionType dimension() const;
    //! \brief The dimension of the state space of the set.
    DimensionType state_dimension() const;
    //! \brief Tests whether the set is empty.
    ValidatedLowerKleenean is_empty() const;
    //! \brief Tests whether the set satisfies the constraint \a c.
    ValidatedKleenean satisfies(EffectiveConstraint c) const;

    //! \brief Returns a bounding box for the set. Computed by a simple interval evaluation of \f$f(D)\f$.
    HybridUpperBoxType bounding_box() const;
    //! \brief Returns an over-approximation to the range of \a g over the set.
    UpperIntervalType range_of(EffectiveScalarFunction const& g) const;
    //! \brief Tests whether the set is disjoint from the box \a hbx.
    ValidatedLowerKleenean separated(const HybridExactBox& hbx) const;
    //! \brief Tests whether the set is a subset of the interior of the box \a hbx.
    ValidatedLowerKleenean inside(const HybridExactBox& hbx) const;
    //! \brief Restricts to a subdomain of the \em parameter domain.
    Void restrict(const ExactBoxType& subdomain);
    //! \brief Adjoins an outer approximation of the set to the grid-based set \a paving, with accuracy given by
    //! \a depth subdivisions in each component.
    Void adjoin_outer_approximation_to(HybridGridTreePaving& paving, Nat depth) const;

    //! \brief Splits into two smaller subsets along parameter direction \a dim.
    Pair<HybridEnclosure,HybridEnclosure> split(Nat dim) const;
    //! \brief Splits into smaller subsets.
    List<HybridEnclosure> split() const;

    //! \brief Reduce the size of the domain by constraint propagation, if possible.
    Void reduce();
    //! \brief Simplifies the representation.
    Void recondition();
    //! \brief Simplifies the representation by changing all uniform errors into independent variables.
    Void uniform_error_recondition();
    //! \brief Simplifies the representation by choosing most significant independent variables to keep, and merging the rest into a single error for each component.
    Void kuhn_recondition();

    //! \brief Draws onto a canvas.
    virtual Void draw(CanvasInterface&, const Set<DiscreteLocation>&, const Variables2d&) const;
    //! \brief Write to an output stream.
    OutputStream& write(OutputStream&) const;
    //! \brief Write an abbreviated representation to an output stream.
    OutputStream& print(OutputStream&) const;
    //! \brief Write a full representation to an output stream which can be used in a constructor.
    OutputStream& repr(OutputStream&) const;
  private:
  public:
    // Compute the flow reach step xi'(s,t) = phi(xi(s),t) and tau'(s,t)=tau(s)+t for t in [0,h] .
    Void _apply_flow(ValidatedVectorFunction phi, FloatDPValue step);
    // Compute the flow reach step xi'(s,t) = phi(xi(s),t) and tau'(s,t)=tau(s)+t for t in [0,h] and t <= eps(xi(s)) .
    Void _apply_flow(ValidatedVectorFunction phi, FloatDPValue step, ValidatedScalarFunction elps);
    // Compute the flow evolve step \xi'(s) = phi(xi(s),eps(s)) and tau'(s)=tau(s)+eps(s)
    Void _apply_flow_step(ValidatedVectorFunction phi, ValidatedScalarFunction elps);
    Void _check() const; // Check that set is well-formed.
    // Compute constraints of the set
    List<ValidatedConstraint> constraints() const;

};

ValidatedLowerKleenean inside(const HybridEnclosure& he, const HybridRealBox& hbx);
inline ValidatedLowerKleenean inside(const HybridEnclosure& he, const HybridExactBox& hbx) { return he.inside(hbx); }
inline ValidatedLowerKleenean separated(const HybridEnclosure& he, const HybridExactBox& hbx) { return he.separated(hbx); }

inline OutputStream& operator<<(OutputStream& os, const HybridEnclosure& s) { return s.write(os); }
inline OutputStream& operator<<(OutputStream& os, const Representation<HybridEnclosure>& s) { return s.pointer->repr(os); }


class HybridGrid;
class HybridGridTreePaving;

template<>
class ListSet<HybridEnclosure>
    : public HybridDrawableInterface
{
  public:
    typedef List<HybridEnclosure>::Iterator Iterator;
    typedef List<HybridEnclosure>::ConstIterator ConstIterator;
  public:
    operator const List<HybridEnclosure>& () const { return _list; }
    operator List<HybridEnclosure>& () { return _list; }
    ListSet() { }
    ListSet(const HybridEnclosure& hes) { this->adjoin(hes); }
    ListSet(const List<HybridEnclosure>& hel) : _list(hel) { }
    Void adjoin(const HybridEnclosure& hes) { this->_list.append(hes); }
    Void adjoin(const ListSet<HybridEnclosure>& hels) {
        for(List<HybridEnclosure>::ConstIterator iter=hels.begin(); iter!=hels.end(); ++iter) {
            this->adjoin(*iter); } }
    Void append(const HybridEnclosure& hes) { this->_list.append(hes); }
    SizeType size() const { return _list.size(); }
    const HybridEnclosure& operator[](Nat i) const { return _list[i]; }
    ListSet<HybridEnclosure::ContinuousStateSetType> operator[](const DiscreteLocation& loc) const;
    Void reduce() { for(auto& set : _list ) { set.reduce(); } }
    Void clear() { this->_list.clear(); }
    HybridEnclosure& front() { return this->_list.front(); }
    const HybridEnclosure& front() const { return this->_list.front(); }
    HybridEnclosure& back() { return this->_list.back(); }
    const HybridEnclosure& back() const { return this->_list.back(); }
    Iterator begin() { return _list.begin(); }
    Iterator end() { return _list.end(); }
    ConstIterator begin() const { return _list.begin(); }
    ConstIterator end() const { return _list.end(); }
    Void draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d& v) const {
        for(Nat i=0; i!=_list.size(); ++i) { _list[i].draw(c,l,v); } }

    friend ValidatedLowerKleenean inside(ListSet<HybridEnclosure> const& set, HybridExactBox const& bx) {
        ValidatedLowerKleenean result=true; for(auto iter=set.begin(); iter!=set.end(); ++iter) { result = result && inside(*iter,bx); } return result; }
    friend OutputStream& operator<<(OutputStream& os, const ListSet<HybridEnclosure>& hls) {
        return os << hls._list; };
    friend HybridGridTreePaving outer_approximation(const ListSet<HybridEnclosure>& hls, const HybridGrid& g, Int d);
  private:
    List<HybridEnclosure> _list;
};

HybridGridTreePaving outer_approximation(const ListSet<HybridEnclosure>& hls, const HybridGrid& g, Nat d);


} // namespace Ariadne

#endif // ARIADNE_HYBRID_ENCLOSURE_HPP
