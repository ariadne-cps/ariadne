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
#include "config.h"

#include "hybrid_simulator.h"

#include "array.h"
#include "container.h"
#include "tuple.h"
#include "stlio.h"
#include "valuation.h"
#include "assignment.h"
#include "space.h"

#include "function.h"

#include "logging.h"

#include "hybrid_set.h"
#include "hybrid_orbit.h"
#include "hybrid_time.h"
#include "hybrid_automaton_interface.h"


namespace Ariadne {

class HybridPoint;
template<class T> class Orbit;

class DegenerateCrossingException { };


HybridSimulator::HybridSimulator()
    : _step_size(0.125)
{
}

void HybridSimulator::set_step_size(double h)
{
    this->_step_size=h;
}


Point make_point(const HybridPoint& hpt, const RealSpace& spc) {
    if(hpt.space()==spc) { return hpt.point(); }
    Map<RealVariable,Float> values=hpt.values();
    Point pt(spc.dimension());
    for(uint i=0; i!=pt.size(); ++i) {
        pt[i]=values[spc.variable(i)];
    }
    return pt;
}

inline Float evaluate(const ScalarFunction<Real>& f, const Vector<Float>& x) { return f(x); }
inline Vector<Float> evaluate(const VectorFunction<Real>& f, const Vector<Float>& x) { return f(x); }

Map<DiscreteEvent,RealScalarFunction> guard_functions(const HybridAutomatonInterface& system, const DiscreteLocation& location) {
    Set<DiscreteEvent> events=system.events(location);
    Map<DiscreteEvent,RealScalarFunction> guards;
    for(Set<DiscreteEvent>::const_iterator iter=events.begin(); iter!=events.end(); ++iter) {
        guards.insert(*iter,system.guard_function(location,*iter));
    }
    return guards;
}

Orbit<HybridPoint>
HybridSimulator::orbit(const HybridAutomatonInterface& system, const HybridPoint& init_pt, const HybridTime& tmax) const
{
    HybridTime t(0.0,0);
    Float h=this->_step_size;

    DiscreteLocation location=init_pt.location();
    RealSpace space=system.continuous_state_space(location);
    Point point=make_point(init_pt,space);
    Point next_point;

    Orbit<HybridPoint> orbit(HybridPoint(location,space,point));

    RealVectorFunction dynamic=system.dynamic_function(location);
    Map<DiscreteEvent,RealScalarFunction> guards=guard_functions(system,location);

    while(t<tmax) {

        bool enabled=false;
        DiscreteEvent event;
        for(Map<DiscreteEvent,RealScalarFunction>::const_iterator guard_iter=guards.begin(); guard_iter!=guards.end(); ++guard_iter) {
            if(evaluate(guard_iter->second,point)>0) {
                enabled=true;
                event=guard_iter->first;
                break;
            }
        }

        if(enabled) {
            DiscreteLocation target=system.target(location,event);
            RealVectorFunction reset=system.reset_function(location,event);
            location=target;
            space=system.continuous_state_space(location);
            next_point=reset(point);

            dynamic=system.dynamic_function(location);
            guards=guard_functions(system,location);
            t._discrete_time+=1;
        } else {
            FloatVector k1,k2,k3,k4;
            Point pt1,pt2,pt3,pt4;
            Point const& pt=point;
            k1=evaluate(dynamic,pt);
            pt1=pt+h*k1;

            k2=evaluate(dynamic,pt1);
            pt2=pt1+(h/2)*k2;

            k3=evaluate(dynamic,pt2);
            pt3=pt1+(h/2)*k3;

            k4=evaluate(dynamic,pt3);

            next_point=pt+(h/6)*(k1+2.0*(k2+k3)+k4);
            t._continuous_time += Real(ExactFloat(h));
        }
        point=next_point;
        orbit.insert(t,HybridPoint(location,space,point));
    }

    return orbit;

}



}  // namespace Ariadne

