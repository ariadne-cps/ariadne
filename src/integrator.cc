/***************************************************************************
 *            integrator.cc
 *
 *  Copyright  2006-10  Pieter Collins
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

#include <iomanip>

#include "integrator.h"

#include "logging.h"
#include "container.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "differential.h"
#include "function.h"
#include "taylor_function.h"

#include "polynomial.h"

namespace Ariadne {

template<class X> Vector<X> join(const Vector<X>& v1, const Vector<X>& v2, const X& s3) {
    Vector<X> r(v1.size()+v2.size()+1u);
    for(uint i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(uint i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    r[v1.size()+v2.size()]=s3;
    return r;
}


Pair<Float,IntervalVector>
IntegratorBase::flow_bounds(const RealVectorFunction& vf, const IntervalVector& dx, const Float& hmax) const
{

    ARIADNE_ASSERT_MSG(vf.result_size()==dx.size(),"vector_field="<<vf<<", states="<<dx);
    ARIADNE_ASSERT_MSG(vf.argument_size()==dx.size(),"vector_field="<<vf<<", states="<<dx);
    ARIADNE_ASSERT(hmax>0);

    // Set up constants of the method.
    // TODO: Better estimates of constants
    const double INITIAL_MULTIPLIER=2.0;
    const double MULTIPLIER=1.125;
    const double BOX_RADIUS_MULTIPLIER=1.25;
    const uint EXPANSION_STEPS=8;
    const uint REDUCTION_STEPS=8;
    const uint REFINEMENT_STEPS=4;

    IntervalVector delta=(dx-midpoint(dx))*(BOX_RADIUS_MULTIPLIER-1);

    // Compute the Lipschitz constant over the initial box
    Float lip = norm(vf.jacobian(dx)).upper();
    Float hlip = this->_lipschitz_tolerance/lip;

    Float hmin=hmax/(1<<REDUCTION_STEPS);
    Float h=max(hmin,min(hmax,hlip));
    ARIADNE_LOG(4,"L="<<lip<<", hL="<<hlip<<", hmax="<<hmax<<"\n");

    bool success=false;
    IntervalVector bx,nbx,df;
    Interval ih(0,h);

    while(!success) {
        ARIADNE_ASSERT_MSG(h>=hmin," h="<<h<<", hmin="<<hmin);
        bx=dx+INITIAL_MULTIPLIER*ih*vf.evaluate(dx)+delta;
        for(uint i=0; i!=EXPANSION_STEPS; ++i) {
            df=vf.evaluate(bx);
            nbx=dx+delta+ih*df;
            if(subset(nbx,bx)) {
                success=true;
                break;
            } else {
                bx=dx+delta+MULTIPLIER*ih*df;
            }
        }
        if(!success) {
            h/=2;
            ih=Interval(0,h);
        }
    }

    ARIADNE_ASSERT(subset(nbx,bx));

    IntervalVector vfbx;
    vfbx=vf.evaluate(bx);

    for(uint i=0; i!=REFINEMENT_STEPS; ++i) {
        bx=nbx;
        vfbx=vf.evaluate(bx);
        nbx=dx+delta+ih*vfbx;
        ARIADNE_ASSERT_MSG(subset(nbx,bx),std::setprecision(20)<<"refinement "<<i<<": "<<nbx<<" is not a inside of "<<bx);
    }


    // Check result of operation
    // We use subset rather than inner subset here since the bound may touch
    ARIADNE_ASSERT(subset(nbx,bx));

    bx=nbx;


    ARIADNE_ASSERT(subset(dx,bx));

    ARIADNE_ASSERT_MSG(subset(dx+h*vf.evaluate(bx),bx),
        "d="<<dx<<"\nh="<<h<<"\nf(b)="<<vf.evaluate(bx)<<"\nd+hf(b)="<<IntervalVector(dx+h*vf.evaluate(bx))<<"\nb="<<bx<<"\n");

    return std::make_pair(h,bx);
}




VectorTaylorFunction
IntegratorBase::flow(const RealVectorFunction& vf, const IntervalVector& dx0, const Real& tmax) const
{
    const uint n=dx0.size(); // Dimension of the state space
    VectorTaylorFunction flow=VectorTaylorFunction::identity(dx0);
    Float t=0.0;
    while(t<tmax) {
        IntervalVector dx=flow.range();
        Float h=tmax-t;
        IntervalVector bx;
        make_lpair(h,bx) = this->flow_bounds(vf,dx,h);
        VectorTaylorFunction step=this->flow_step(vf,dx,h,bx);
        step=partial_evaluate(step,n,h);
        flow=compose(step,flow);
        t=t+h;
    }
    return flow;
}


VectorTaylorFunction
IntegratorBase::flow(const RealVectorFunction& vf, const IntervalVector& dx0, const Interval& dt) const
{
    Real dtl=Real(dt.lower());
    VectorTaylorFunction evolve=this->flow(vf,dx0,dtl);
    Float dtw=dt.upper()-dt.lower();
    IntervalVector dx=evolve.range();
    Float h;
    IntervalVector bx;
    make_lpair(h,bx) = this->flow_bounds(vf,dx,h);
    ARIADNE_ASSERT_MSG(dtw==h,"Width of time interval "<<dt<<" cannot be covered in a single flow step; maximum flow step "<<h<<" over domain "<<dx);
    VectorTaylorFunction step=this->flow_step(vf,dx,h,bx);
    ScalarTaylorFunction time=ScalarTaylorFunction::identity(dt)-dt.lower();
    VectorTaylorFunction flow=compose(step,combine(evolve,time));
    return flow;
}



VectorTaylorFunction
IntegratorBase::flow_step(const RealVectorFunction& vf, const IntervalVector& dx, const Float& hmax) const
{
    Float h;
    IntervalVector bx;
    make_lpair(h,bx)=this->flow_bounds(vf,dx,hmax);
    return this->flow_step(vf,dx,h,bx);
}

VectorTaylorFunction
TaylorIntegrator::flow_step(const RealVectorFunction& f, const IntervalVector& dx, const Float& h, const IntervalVector& bx) const
{
    ARIADNE_LOG(2,"f="<<f<<" dx="<<dx<<" h="<<h<<" bx="<<bx<<"\n");
    const uint nx=dx.size();
    const double sw=this->_sweep_threshold;

    IntervalVector dom=join(dx,Interval(-h,h));
    ARIADNE_LOG(3,"dom="<<dom<<"\n");

/*
    VectorTaylorFunction phi0(VectorTaylorFunction(nx,ScalarTaylorFunction(dom)));
    for(uint i=0; i!=nx; ++i) { phi0[i]=ScalarTaylorFunction::coordinate(dom,i); }
    ARIADNE_LOG(0,"phi0="<<phi0<<"\n");
*/
    VectorTaylorFunction phi0(nx,ScalarTaylorFunction(dom));
    for(uint i=0; i!=nx; ++i) { phi0[i]=ScalarTaylorFunction::coordinate(dom,i); phi0[i].set_sweep_threshold(sw); }
    ARIADNE_LOG(3,"phi0="<<phi0<<"\n");

    VectorTaylorFunction phi(nx,ScalarTaylorFunction(dom));
    for(uint i=0; i!=nx; ++i) { phi[i]=ScalarTaylorFunction::constant(dom,bx[i]); phi0[i].set_sweep_threshold(sw); }

    /*
    for(uint k=0; k!=this->_temporal_order; ++k) {
        for(uint i=0; i!=nx; ++i) {
            phi[np+i]=antiderivative(compose(vf[i],phi),np+nx)+phi0[i];
        }
    }
    */

    ARIADNE_LOG(4,"phi="<<phi<<"\n");
    for(uint k=0; k!=this->temporal_order(); ++k) {
        bool last_step=(phi.error()<this->maximum_error());
        VectorTaylorFunction fphi=compose(f,phi);
        ARIADNE_LOG(4,"fphi="<<fphi<<"\n");
        for(uint i=0; i!=nx; ++i) {
            phi[i]=antiderivative(fphi[i],nx)+phi0[i];
        }
        ARIADNE_LOG(3,"phi="<<phi<<"\n");
        if(last_step) { break; }
    }

    if(phi.error()>this->maximum_error()) {
        ARIADNE_WARN("Integration of "<<f<<" starting in "<<dx<<" for time "<<h<<" has error "<<phi.error()<<" after "<<this->temporal_order()<<" iterations, which exceeds maximum error "<<this->maximum_error()<<"\n");
    }

    VectorTaylorFunction res(nx,ScalarTaylorFunction(dom));
    for(uint i=0; i!=nx; ++i) { res[i]=phi[i]; res[i].set_sweep_threshold(sw); res[i].sweep(); }
    ARIADNE_LOG(3,"res="<<res<<"\n");
    return res;

}


} // namespace Ariadne
