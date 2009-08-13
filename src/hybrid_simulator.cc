/***************************************************************************
 *            hybrid_simulator.cc
 *
 *  Copyright  2008  Alberto Casagrande, Pieter Collins
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

#include "macros.h"
#include "array.h"
#include "tuple.h"
#include "stlio.h"
#include "valuation.h"

#include "function_interface.h"

#include "logging.h"

#include "hybrid_time.h"
#include "hybrid_system.h"
#include "hybrid_simulator.h"

class B1 { public: void foo(int) { } };
class B2 { public: void foo(double) { } };
class D : public B1, public B2 {
  public:
    //void foo(int n) { B1::foo(n); }
    //void foo(double x) { B2::foo(x); }
    using B1::foo;
    using B2::foo;
};

void f1(D d) { d.foo(3); }
void f2(D d) { d.foo(1.4); }

namespace Ariadne {

template<class T> class Orbit;

template<>
class Orbit<HybridPoint>
    : public std::map<HybridTime,HybridPoint>
{
};

EnumeratedValue evaluate(const Expression<EnumeratedValue>& e, const DiscreteValuation& loc);
Float evaluate(const Expression<Real>& e, const ContinuousValuation<Float>& pt);
tribool evaluate(const Expression<tribool>& e, const ContinuousValuation<Float>& pt);

class DegenerateCrossingException { };


Simulator<HybridSystem>::Simulator()
{
}


Simulator<HybridSystem>*
Simulator<HybridSystem>::clone() const
{
    return new Simulator<HybridSystem>(*this);
}



class HybridPoint : public Valuation<Float> { };

class HybridVector : public ContinuousValuation<Float> { };

HybridPoint& operator+=(HybridPoint& hpt, const HybridVector& hv) {
    Map<Identifier,Float>& ptmp=hpt.real_values();
    const Map<Identifier,Float>& vmp=hv.real_values();
    for(Map<Identifier,Float>::const_iterator iter=vmp.begin(); iter!=vmp.end(); ++iter) {
        ptmp[iter->first]+=iter->second;
    }
    return hpt;
}

HybridVector& operator+=(HybridVector& hv1, const HybridVector& hv2) {
    Map<Identifier,Float>& v1mp=hv1.real_values();
    for(Map<Identifier,Float>::const_iterator iter=hv2.real_values().begin(); iter!=hv2.real_values().end(); ++iter) {
        v1mp[iter->first]+=iter->second;
    }
    return hv1;
}

HybridVector& operator*=(HybridVector& hv, const Float& s) {
    Map<Identifier,Float>& vmp=hv.real_values();
    for(Map<Identifier,Float>::iterator iter=vmp.begin(); iter!=vmp.end(); ++iter) {
        iter->second*=s;
    }
    return hv;
}

HybridPoint operator+(const HybridPoint& hpt, const HybridVector& hv) {
    HybridPoint r(hpt); r+=hv; return r;
}

HybridVector operator+(const HybridVector& hv1, const HybridVector& hv2) {
    HybridVector r(hv1); r+=hv2; return r;
}

HybridVector operator*(const HybridVector& hv, const Float& s) {
    HybridVector r(hv); r*=s; return r;
}

HybridVector operator*(const Float& s, const HybridVector& hv) {
    HybridVector r(hv); r*=s; return r;
}




HybridVector evaluate(const std::vector<RealDynamic>& dyn, const HybridPoint& pt) {
    HybridVector r;
    for(std::vector<RealDynamic>::const_iterator dyn_iter=dyn.begin(); dyn_iter!=dyn.end(); ++dyn_iter) {
        r[dyn_iter->lhs.base]=evaluate(dyn_iter->rhs,pt);
    }
    return r;
}

void evaluate(const std::vector<RealAssignment>& alg, HybridPoint& pt) {
    for(std::vector<RealAssignment>::const_iterator alg_iter=alg.begin(); alg_iter!=alg.end(); ++alg_iter) {
        pt[alg_iter->lhs]=evaluate(alg_iter->rhs,pt);
    }
    return;
}



Orbit<HybridPoint>
Simulator<HybridSystem>::orbit(const HybridSystem& sys, const HybridPoint& init_pt, const HybridTime& tmax) const
{
    HybridTime t(0.0,0);
    Orbit<HybridPoint> orbit;
    Float h;

    HybridPoint pt(init_pt);
    HybridPoint next_pt;

    std::map<Event,ContinuousPredicate> guards=sys.guards(pt);
    std::vector<RealAssignment> algebraic_assignments=sys.equations(pt);
    std::vector<RealDynamic> differential_assignments=sys.dynamic(pt);

    while(t<tmax) {

        bool enabled=false;
        Event event;
        for(std::map<Event,ContinuousPredicate>::const_iterator guard_iter=guards.begin(); guard_iter!=guards.end(); ++guard_iter) {
            if(evaluate(guard_iter->second,pt)) {
                enabled=true;
                event=guard_iter->first;
                break;
            }
        }

        if(enabled) {
            std::vector<EnumeratedUpdate> switchings=sys.switching(event,pt);
            for(std::vector<EnumeratedUpdate>::const_iterator iter=switchings.begin(); iter!=switchings.end(); ++iter) {
                next_pt.set(iter->lhs.base,evaluate(iter->rhs,pt));
            }
            std::vector<RealUpdate> real_updates=sys.reset(event,pt);
            for(std::vector<RealUpdate>::const_iterator iter=real_updates.begin(); iter!=real_updates.end(); ++iter) {
                next_pt.set(iter->lhs.base,evaluate(iter->rhs,pt));
            }
            algebraic_assignments=sys.equations(pt);
            for(std::vector<RealAssignment>::const_iterator iter=algebraic_assignments.begin(); iter!=algebraic_assignments.end(); ++iter) {
                next_pt.set(iter->lhs,evaluate(iter->rhs,next_pt));
            }
            differential_assignments=sys.dynamic(next_pt);
            guards=sys.guards(next_pt);
            t.discrete_time+=1;
        } else {
            HybridVector k1,k2,k3,k4;
            HybridPoint pt1,pt2,pt3,pt4;
            k1=evaluate(differential_assignments,pt);
            pt1=pt+h*k1;
            evaluate(algebraic_assignments,pt1);

            k2=evaluate(differential_assignments,pt1);
            pt2=pt1+(h/2)*k2;
            evaluate(algebraic_assignments,pt2);

            k3=evaluate(differential_assignments,pt2);
            pt3=pt1+(h/2)*k3;
            evaluate(algebraic_assignments,pt3);

            k4=evaluate(differential_assignments,pt3);

            next_pt=pt+(h/6)*(k1+2*(k2+k3)+k4);
            t.continuous_time+=h;
        }
        pt=next_pt;
        orbit[t]=pt;
    }

    return orbit;

}



}  // namespace Ariadne

