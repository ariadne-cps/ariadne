/***************************************************************************
 *            hybrid/hybrid_simulator.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "../function/functional.hpp"

#include "../algebra/algebra.hpp"

#include "../config.hpp"

#include "../hybrid/hybrid_simulator.hpp"

#include "../utility/array.hpp"
#include "../utility/container.hpp"
#include "../utility/tuple.hpp"
#include "../utility/stlio.hpp"
#include "../symbolic/valuation.hpp"
#include "../symbolic/assignment.hpp"
#include "../symbolic/space.hpp"

#include "../function/function.hpp"
#include "../function/formula.hpp"
#include "../function/taylor_model.hpp"

#include "../solvers/runge_kutta_integrator.hpp"

#include "../output/logging.hpp"

#include "../hybrid/hybrid_set.hpp"
#include "../hybrid/hybrid_orbit.hpp"
#include "../hybrid/hybrid_time.hpp"
#include "../hybrid/hybrid_automaton_interface.hpp"


namespace Ariadne {

template<class T> class Orbit;

class DegenerateCrossingException { };

HybridSimulator::HybridSimulator()
    : _step_size(0.125)
{
    this->charcode="s";
}

Void HybridSimulator::set_step_size(double h)
{
    this->_step_size=h;
}

namespace {

Map<DiscreteEvent,EffectiveScalarMultivariateFunction> guard_functions(const HybridAutomatonInterface& system, const DiscreteLocation& location) {
    Set<DiscreteEvent> events=system.events(location);
    Map<DiscreteEvent,EffectiveScalarMultivariateFunction> guards;
    for(Set<DiscreteEvent>::ConstIterator iter=events.begin(); iter!=events.end(); ++iter) {
        EventKind kind = system.event_kind(location,*iter);
        if (kind != EventKind::INVARIANT && kind != EventKind::PROGRESS)
            guards.insert(*iter,system.guard_function(location,*iter));
    }
    return guards;
}

Bool satisfies_invariants(const HybridAutomatonInterface& system, const DiscreteLocation& location, const Point<FloatDPApproximation>& point) {
    Bool result=true;
    Set<DiscreteEvent> events=system.events(location);
    for(Set<DiscreteEvent>::ConstIterator iter=events.begin(); iter!=events.end(); ++iter) {
        EventKind kind = system.event_kind(location,*iter);
        EffectiveScalarMultivariateFunction guard = system.guard_function(location,*iter);
        if ((kind == EventKind::INVARIANT || kind == EventKind::PROGRESS) && probably(evaluate(guard,point)<0)) {
            result=false;
            break;
        }
    }
    return result;
}

template<class X> Point<X> make_point(const HybridPoint<X>& hpt, const RealSpace& sspc) {
    if(hpt.space()==sspc) { return hpt.point(); }
    Map<RealVariable,X> values=hpt.values();
    Point<X> pt(sspc.dimension());
    for(Nat i=0; i!=pt.size(); ++i) {
        pt[i]=values[sspc.variable(i)];
    }
    return pt;
}

template<class X> HybridPoint<X> make_hybrid_state_auxiliary_point(const DiscreteLocation& location, const Point<X>& spt,
        const RealSpace& sspc, const RealSpace& aspc, const RealSpace& saspc, const EffectiveVectorMultivariateFunction& auxiliary_function) {
    Point<X> sapt(saspc.dimension());
    Point<X> apt = evaluate(auxiliary_function,spt);
    for(Nat i=0; i!=sapt.size(); ++i) {
        RealVariable var = saspc.variable(i);
        sapt[i]= sspc.contains(var) ? spt[sspc[var]] : apt[aspc[var]];
    }
    return HybridPoint<X>(location,saspc,sapt);
}

}

inline FloatDPApproximation evaluate(const EffectiveScalarMultivariateFunction& f, const Vector<FloatDPApproximation>& x) { return f(x); }
inline Vector<FloatDPApproximation> evaluate(const EffectiveVectorMultivariateFunction& f, const Vector<FloatDPApproximation>& x) { return f(x); }

auto HybridSimulator::orbit(const HybridAutomatonInterface& system, const HybridRealPoint& init_pt, const TerminationType& termination) const
    -> Orbit<HybridApproximatePointType>
{
    return orbit(system,HybridApproximatePoint(init_pt,dp),termination);
}

auto HybridSimulator::orbit(const HybridAutomatonInterface& system,
                            const HybridApproximatePointType& init_pt,
                            const TerminationType& termination) const
    -> Orbit<HybridApproximatePointType>
{
    DoublePrecision pr;
    HybridTime t(0.0,0);
    Dyadic h(ExactDouble(this->_step_size.get_d()));
    HybridTime tmax(termination.maximum_time(),termination.maximum_steps());

    DiscreteLocation location=init_pt.location();
    RealSpace continuous_state_space=system.continuous_state_space(location);
    RealSpace continuous_auxiliary_space = system.continuous_auxiliary_space(location);
    RealSpace continuous_state_auxiliary_space = system.state_auxiliary_space()[location];
    EffectiveVectorMultivariateFunction auxiliary_function = system.auxiliary_function(location);
    ApproximatePointType point=make_point(init_pt,continuous_state_space);
    ApproximatePointType next_point;
    List<DiscreteEvent> event_trace;

    Orbit<HybridApproximatePoint> orbit(make_hybrid_state_auxiliary_point(location,point,continuous_state_space,
            continuous_auxiliary_space,continuous_state_auxiliary_space,system.auxiliary_function(location)));

    EffectiveVectorMultivariateFunction dynamic=system.dynamic_function(location);
    Map<DiscreteEvent,EffectiveScalarMultivariateFunction> guards=guard_functions(system,location);

    RungeKutta4Integrator integrator(this->_step_size.get_d());

    while(possibly(t<tmax) && (event_trace.empty() || !termination.terminating_events().contains(event_trace.back()))) {
        Int old_precision = std::clog.precision();
        ARIADNE_LOG(1, (verbosity == 1 ? "\r" : "")
                << "t=" << std::setw(4) << std::left << t.continuous_time().lower().get(pr)
                << " #e=" << std::left << t.discrete_time()
                << " p=" << point
                << " l=" << std::left << location
                << " e=" << std::left << event_trace
                << " \n" << std::setprecision(old_precision));

        if (not satisfies_invariants(system, location, point)) {
            ARIADNE_LOG(2,"invariant/progress condition not satisfied, stopping evolution.\n");
            break;
        }

        Bool enabled=false;
        DiscreteEvent event;
        for(Map<DiscreteEvent,EffectiveScalarMultivariateFunction>::ConstIterator guard_iter=guards.begin(); guard_iter!=guards.end(); ++guard_iter) {
            if(probably(evaluate(guard_iter->second,point)>0)) {
                enabled=true;
                event=guard_iter->first;
                break;
            }
        }

        if(enabled) {
            DiscreteLocation target=system.target(location,event);
            EffectiveVectorMultivariateFunction reset=system.reset_function(location,event);
            location=target;
            continuous_state_space=system.continuous_state_space(location);
            continuous_auxiliary_space=system.continuous_auxiliary_space(location);
            continuous_state_auxiliary_space=system.state_auxiliary_space()[location];
            auxiliary_function=system.auxiliary_function(location);
            next_point=reset(point);
            event_trace.push_back(event);

            ARIADNE_LOG(2,"event " << event << " enabled: next point " << next_point << ", on location " << target << "\n");

            dynamic=system.dynamic_function(location);
            guards=guard_functions(system,location);

            t._discrete_time += 1;
        } else {
            next_point = ApproximatePointType(integrator.step(dynamic,point,this->_step_size));
            t._continuous_time += h;
        }
        point=next_point;
        orbit.insert(t,make_hybrid_state_auxiliary_point(location,point,
                continuous_state_space,continuous_auxiliary_space,continuous_state_auxiliary_space,auxiliary_function));
    }

    return orbit;
}



}  // namespace Ariadne

