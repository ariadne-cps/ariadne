/***************************************************************************
 *            hybrid_evolver-constrained.h
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
 *  \brief Evolver for hybrid systems.
 */

#ifndef ARIADNE_HYBRID_EVOLVER_CONSTRAINED_H
#define ARIADNE_HYBRID_EVOLVER_CONSTRAINED_H

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

template<class Sys, class ES> class Evolver;

class ConstrainedImageSet;
class HybridAutomaton;

class EvolutionParameters;

class EvolutionProfiler;

class HybridTime;

class DiscreteEvent;

template<class ES> class ListSet;

class ConstrainedImageSet
    : public DrawableInterface
{
    class Data;
    typedef Vector<Interval> IntervalVector;
  public:
    ConstrainedImageSet();
    ConstrainedImageSet(IntervalVector);
    ConstrainedImageSet(IntervalVector, VectorFunction);

    ConstrainedImageSet(const ConstrainedImageSet&);
    ConstrainedImageSet& operator=(const ConstrainedImageSet&);
    ConstrainedImageSet* clone() const;

    void apply_map(VectorFunction);
    void apply_flow(VectorFunction, Interval);
    void apply_flow(VectorTaylorFunction);

    void new_invariant(DiscreteEvent,ScalarFunction,ScalarFunction);
    void new_activation(DiscreteEvent,ScalarFunction,ScalarFunction);
    void new_guard(DiscreteEvent,ScalarFunction,ScalarFunction);
    void new_equation(DiscreteEvent,ScalarFunction);

    uint dimension() const;
    uint number_of_parameters() const;
    Box bounding_box() const;
    tribool disjoint(Box bx);
    GridTreeSet outer_approximation(Grid g, int depth) const;

    void draw(CanvasInterface&) const;
    std::ostream& write(std::ostream&) const;
  private:
    shared_ptr<Data> _data; 
};

inline std::ostream& operator<<(std::ostream& os, const ConstrainedImageSet& s) { return s.write(os); }

typedef HybridBasicSet<ConstrainedImageSet> HybridConstrainedImageSet;

template<>
struct HybridBasicSet<ConstrainedImageSet>
{
    typedef DiscreteState first_type;
    typedef ConstrainedImageSet ContinuousStateSetType;

    HybridBasicSet(const DiscreteState& q_, const Box& s_)
        : _events(), _location(q_), _set(s_) { }
    HybridBasicSet(const std::pair<DiscreteState,ConstrainedImageSet>& hs_)
        : _events(), _location(hs_.first), _set(hs_.second) { }
    HybridBasicSet(const DiscreteState& q_, const ConstrainedImageSet& s_)
        : _events(), _location(q_), _set(s_) { }
    HybridBasicSet(const List<DiscreteEvent>& e_, const DiscreteState& q_, const ConstrainedImageSet& s_)
        : _events(e_), _location(q_), _set(s_) { }
    uint dimension() const { return _set.dimension(); }
    DiscreteState location() const { return _location; }
    ConstrainedImageSet const& continuous_state_set() const { return _set; }

    List<DiscreteEvent> _events;
    DiscreteState _location;
    ConstrainedImageSet _set;
};
 
inline std::ostream& operator<<(std::ostream& os, const HybridConstrainedImageSet& hs) {
    return os << hs._events << ":" << hs._location <<"," << hs._set;
}


struct TimedHybridConstrainedImageSet : public HybridConstrainedImageSet {
    TimedHybridConstrainedImageSet(const HybridConstrainedImageSet& hs_)
        : HybridConstrainedImageSet(hs_), _time(ScalarFunction::constant(hs_.dimension(),0.0)) { }
    TimedHybridConstrainedImageSet(DiscreteState q_, const ConstrainedImageSet& s_)
        : HybridConstrainedImageSet(q_,s_), _time(ScalarFunction::constant(s_.dimension(),0.0)) { }
    TimedHybridConstrainedImageSet(List<DiscreteEvent> e_, DiscreteState q_, const ConstrainedImageSet& s_, const ScalarFunction& t_)
        : HybridConstrainedImageSet(e_,q_,s_), _time(t_) { }
    ScalarFunction _time;
};


/*! \brief A class for computing the evolution of a hybrid system.
 *
 * The actual evolution steps are performed by the HybridEvolver class.
 */
class ConstrainedImageSetHybridEvolver
    : public EvolverBase<HybridAutomaton,HybridConstrainedImageSet>
//    , public Loggable
{
  public:
    typedef ContinuousEvolutionParameters EvolutionParametersType;
    typedef HybridAutomaton::TimeType TimeType;
    typedef int IntegerType;
    typedef Float RealType;
    typedef std::vector<DiscreteEvent> EventListType;
    typedef HybridAutomaton SystemType;
    typedef TaylorSet ContinuousEnclosureType;
    typedef HybridConstrainedImageSet HybridEnclosureType;
    typedef HybridEnclosureType EnclosureType;
    typedef Orbit<EnclosureType> OrbitType;
    typedef ListSet<HybridConstrainedImageSet> EnclosureListType;
    typedef Float ContinuousTimeType;
  public:

    //! \brief Default constructor.
    ConstrainedImageSetHybridEvolver() : _parameters() { }

    //! \brief Construct from parameters using a default integrator.
    ConstrainedImageSetHybridEvolver(const EvolutionParametersType& parameters) :
        _parameters(new EvolutionParametersType(parameters)) { }

    /*! \brief Make a dynamically-allocated copy. */
    ConstrainedImageSetHybridEvolver* clone() const { return new ConstrainedImageSetHybridEvolver(*this); }

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

    void
    _upper_evolution_step(List<TimedHybridConstrainedImageSet>&,
                          ListSet<HybridConstrainedImageSet>&,
                          ListSet<HybridConstrainedImageSet>&,
                          ListSet<HybridConstrainedImageSet>&,
                          HybridAutomaton const& system,
                          HybridTime const& maximum_time) const;

 private:
    boost::shared_ptr< EvolutionParametersType > _parameters;
    boost::shared_ptr< EvolutionProfiler >  _profiler;
};



} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_CONSTRAINED_H
