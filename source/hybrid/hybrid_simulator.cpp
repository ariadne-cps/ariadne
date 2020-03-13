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

template<class X> Point<X> make_point(const HybridPoint<X>& hpt, const RealSpace& spc) {
    if(hpt.space()==spc) { return hpt.point(); }
    Map<RealVariable,X> values=hpt.values();
    Point<X> pt(spc.dimension());
    for(Nat i=0; i!=pt.size(); ++i) {
        pt[i]=values[spc.variable(i)];
    }
    return pt;
}

}

inline FloatDPApproximation evaluate(const EffectiveScalarMultivariateFunction& f, const Vector<FloatDPApproximation>& x) { return f(x); }
inline Vector<FloatDPApproximation> evaluate(const EffectiveVectorMultivariateFunction& f, const Vector<FloatDPApproximation>& x) { return f(x); }

auto HybridSimulator::orbit(const HybridAutomatonInterface& system, const HybridRealPoint& init_pt, const HybridTime& tmax) const
    -> Orbit<HybridApproximatePointType>
{
    return orbit(system,HybridApproximatePoint(init_pt,dp),tmax);
}

auto HybridSimulator::orbit(const HybridAutomatonInterface& system,
                            const HybridApproximatePointType& init_pt,
                            const HybridTime& tmax) const
    -> Orbit<HybridApproximatePointType>
{
    DoublePrecision pr;
    HybridTime t(0.0,0);
    Dyadic h(ExactDouble(this->_step_size.get_d()));

    DiscreteLocation location=init_pt.location();
    RealSpace space=system.continuous_state_space(location);
    ApproximatePointType point=make_point(init_pt,space);
    ApproximatePointType next_point;
    List<DiscreteEvent> event_trace;

    Orbit<HybridApproximatePoint> orbit(HybridApproximatePoint(location,space,cast_exact(point)));

    EffectiveVectorMultivariateFunction dynamic=system.dynamic_function(location);
    Map<DiscreteEvent,EffectiveScalarMultivariateFunction> guards=guard_functions(system,location);

    while(possibly(check(t<tmax,Effort::get_default()))) {
        Int old_precision = std::clog.precision();
        ARIADNE_LOG(1,(verbosity==1?"\r":"")
                <<"t="<<std::setw(4)<<std::left<<t.continuous_time().lower().get(pr)
                <<" #e="<<std::setw(4)<<std::left<<t.discrete_time()
                <<" p="<<std::setw(5)<<std::left<<point
                <<" l="<<std::left<<location
                <<" e="<<std::left<<event_trace
                <<" \n"<<std::setprecision(old_precision));

        Bool enabled=false;
        DiscreteEvent event;
        for(Map<DiscreteEvent,EffectiveScalarMultivariateFunction>::ConstIterator guard_iter=guards.begin(); guard_iter!=guards.end(); ++guard_iter) {
            if(probably(evaluate(guard_iter->second,point)>0)) {
                enabled=true;
                event=guard_iter->first;
                ARIADNE_LOG(2,"event " << event << " enabled.\n");
                break;
            }
        }

        if(enabled) {
            DiscreteLocation target=system.target(location,event);
            EffectiveVectorMultivariateFunction reset=system.reset_function(location,event);
            location=target;
            space=system.continuous_state_space(location);
            next_point=reset(point);
            event_trace.push_back(event);

            ARIADNE_LOG(2,"next point: " << next_point << ", on location " << target << "\n");

            dynamic=system.dynamic_function(location);
            guards=guard_functions(system,location);
            t._discrete_time+=1;
        } else {
            FloatDPApproximationVector k1,k2,k3,k4;
            ApproximatePointType pt1,pt2,pt3,pt4;

            ApproximatePointType const& pt=point;
            k1=evaluate(dynamic,pt);
            pt1=pt+h*k1;

            k2=evaluate(dynamic,pt1);
            pt2=pt1+(h/2)*k2;

            k3=evaluate(dynamic,pt2);
            pt3=pt1+(h/2)*k3;

            k4=evaluate(dynamic,pt3);

            next_point=pt+(h/6)*(k1+FloatDPApproximation(2.0)*(k2+k3)+k4);
            t._continuous_time += h;
        }
        point=next_point;
        orbit.insert(t,HybridApproximatePoint(location,space,point));
    }

    return orbit;

}



}  // namespace Ariadne

