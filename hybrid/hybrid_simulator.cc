/***************************************************************************
 *            hybrid_simulator.cc
 *
 *  Copyright  2008-11  Pieter Collins
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

#include "function/functional.h"
#include "config.h"

#include "hybrid/hybrid_simulator.h"

#include "utility/array.h"
#include "utility/container.h"
#include "utility/tuple.h"
#include "utility/stlio.h"
#include "expression/valuation.h"
#include "expression/assignment.h"
#include "expression/space.h"

#include "function/function.h"

#include "utility/logging.h"

#include "hybrid/hybrid_set.h"
#include "hybrid/hybrid_orbit.h"
#include "hybrid/hybrid_time.h"
#include "hybrid/hybrid_automaton_interface.h"


namespace Ariadne {

class HybridPoint;
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


ExactPoint make_point(const HybridPoint& hpt, const RealSpace& spc) {
    if(hpt.space()==spc) { return hpt.point(); }
    Map<RealVariable,Float64Value> values=hpt.values();
    ExactPoint pt(spc.dimension());
    for(Nat i=0; i!=pt.size(); ++i) {
        pt[i]=values[spc.variable(i)];
    }
    return pt;
}

inline Float64Approximation evaluate(const EffectiveScalarFunction& f, const Vector<Float64Approximation>& x) { return f(x); }
inline Vector<Float64Approximation> evaluate(const EffectiveVectorFunction& f, const Vector<Float64Approximation>& x) { return f(x); }

Map<DiscreteEvent,EffectiveScalarFunction> guard_functions(const HybridAutomatonInterface& system, const DiscreteLocation& location) {
    Set<DiscreteEvent> events=system.events(location);
    Map<DiscreteEvent,EffectiveScalarFunction> guards;
    for(Set<DiscreteEvent>::ConstIterator iter=events.begin(); iter!=events.end(); ++iter) {
        guards.insert(*iter,system.guard_function(location,*iter));
    }
    return guards;
}

Orbit<HybridPoint>
HybridSimulator::orbit(const HybridAutomatonInterface& system, const HybridPoint& init_pt, const HybridTime& tmax) const
{
    HybridTime t(0.0,0);
    Float64Approximation h=this->_step_size;

    DiscreteLocation location=init_pt.location();
    RealSpace space=system.continuous_state_space(location);
    ApproximatePoint point=make_point(init_pt,space);
    ApproximatePoint next_point;

    Orbit<HybridPoint> orbit(HybridPoint(location,space,cast_exact(point)));

    EffectiveVectorFunction dynamic=system.dynamic_function(location);
    Map<DiscreteEvent,EffectiveScalarFunction> guards=guard_functions(system,location);

    while(possibly(t<tmax)) {

        Bool enabled=false;
        DiscreteEvent event;
        for(Map<DiscreteEvent,EffectiveScalarFunction>::ConstIterator guard_iter=guards.begin(); guard_iter!=guards.end(); ++guard_iter) {
            if(probably(evaluate(guard_iter->second,point)>0)) {
                enabled=true;
                event=guard_iter->first;
                break;
            }
        }

        if(enabled) {
            DiscreteLocation target=system.target(location,event);
            EffectiveVectorFunction reset=system.reset_function(location,event);
            location=target;
            space=system.continuous_state_space(location);
            next_point=reset(point);

            dynamic=system.dynamic_function(location);
            guards=guard_functions(system,location);
            t._discrete_time+=1;
        } else {
            FloatApproximationVector k1,k2,k3,k4;
            ApproximatePoint pt1,pt2,pt3,pt4;

            ApproximatePoint const& pt=point;
            k1=evaluate(dynamic,pt);
            pt1=pt+h*k1;

            k2=evaluate(dynamic,pt1);
            pt2=pt1+(h/2)*k2;

            k3=evaluate(dynamic,pt2);
            pt3=pt1+(h/2)*k3;

            k4=evaluate(dynamic,pt3);

            next_point=pt+(h/6)*(k1+Float64Approximation(2.0)*(k2+k3)+k4);
            t._continuous_time += Real(Float64Value(h.raw()));
        }
        point=next_point;
        orbit.insert(t,HybridPoint(location,space,cast_exact(point)));
    }

    return orbit;

}



}  // namespace Ariadne

