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

#include "config.h"

#include "integrator.h"

#include "logging.h"
#include "container.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "differential.h"
#include "function.h"
#include "function_model.h"

#include "polynomial.h"
#include <include/function.h>

namespace Ariadne {

FunctionModelFactoryInterface<Interval>* make_taylor_function_factory();
FunctionModelFactoryInterface<Interval>* make_taylor_function_factory(double s);

IntegratorBase::IntegratorBase(MaximumError e, LipschitzConstant l)
    :  _maximum_error(e), _lipschitz_tolerance(l), _function_factory_ptr(make_taylor_function_factory())
{
    ARIADNE_PRECONDITION(e>0.0);
    ARIADNE_PRECONDITION(l>0.0)
}

IntegratorBase::IntegratorBase(MaximumError e, SweepThreshold s, LipschitzConstant l)
    :  _maximum_error(e), _lipschitz_tolerance(l), _function_factory_ptr(make_taylor_function_factory(s))
{
    ARIADNE_PRECONDITION(e>0.0);
    ARIADNE_PRECONDITION(l>0.0);
}

void
IntegratorBase::set_function_factory(const IntervalFunctionModelFactoryInterface& factory)
{
    this->_function_factory_ptr=FunctionFactoryPointer(factory.clone());
}

const IntervalFunctionModelFactoryInterface&
IntegratorBase::function_factory() const
{
    return *this->_function_factory_ptr;
}


Pair<Float,IntervalVector>
IntegratorBase::flow_bounds(const IntervalVectorFunction& vf, const IntervalVector& dx, const Float& hmax) const
{
    ARIADNE_LOG(3,"IntegratorBase::flow_bounds(IntervalVectorFunction vf, IntervalVector dx, Float hmax)\n");
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





IntervalVectorFunctionModel
IntegratorBase::flow(const IntervalVectorFunction& vf, const IntervalVector& dx0, const Real& tmax) const
{
    ARIADNE_LOG(1,"IntegratorBase::flow(IntervalVectorFunction vf, IntervalVector dx0, Real tmax)\n");
    ARIADNE_LOG(2,"vf="<<vf<<"\n");
    ARIADNE_LOG(2,"dom(x0)="<<dx0<<" tmax="<<tmax<<"\n");
    const uint n=dx0.size(); // Dimension of the state space
    IntervalVectorFunctionModel flow_function=this->function_factory().create_identity(dx0);
    Float t=0.0;
    IntervalVectorFunctionModel step_function;
    while(t<tmax) {
        IntervalVector dx=flow_function.range();
        Float h=tmax-t;
        IntervalVector bx;
        make_lpair(h,bx) = this->flow_bounds(vf,dx,h);
        bool flow_successfully_computed=false;
        while(!flow_successfully_computed) {
            try {
                step_function=this->flow_step(vf,dx,h,bx);
                flow_successfully_computed=true;
            } catch(FlowTimeStepException e) {
                h/=2;
            }
        }
        step_function=partial_evaluate(step_function,n,numeric_cast<Interval>(h));
        flow_function=compose(step_function,flow_function);
        t=t+h;
    }
    return flow_function;
}


IntervalVectorFunctionModel
IntegratorBase::flow(const IntervalVectorFunction& vf, const IntervalVector& dx0, const Interval& dt) const
{
    ARIADNE_LOG(1,"IntegratorBase::flow(IntervalVectorFunction vf, IntervalVector dx0, Interval dt)\n");
    Real dtl=Real(ExactFloat(dt.lower()));
    IntervalVectorFunctionModel evolve=this->flow(vf,dx0,dtl);
    Float dtw=dt.upper()-dt.lower();
    IntervalVector dx=evolve.range();
    Float h;
    IntervalVector bx;
    make_lpair(h,bx) = this->flow_bounds(vf,dx,dtw);
    if(dtw!=h) {
        ARIADNE_THROW(FlowTimeStepException,"IntegratorBase::flow","Width of time interval "<<dt<<" cannot be covered in a single flow step; maximum flow step "<<h<<" over domain "<<dx);
    }
    IntervalVectorFunctionModel step=this->flow_step(vf,dx,h,bx);
    IntervalScalarFunctionModel time=this->function_factory().create_identity(dt)-ExactFloat(dt.lower());
    IntervalVectorFunctionModel flow=compose(step,combine(evolve,time));
    return flow;
}



IntervalVectorFunctionModel
IntegratorBase::flow_step(const IntervalVectorFunction& vf, const IntervalVector& dx, Float& h) const
{
    ARIADNE_LOG(3,"IntegratorBase::flow_step(IntervalVectorFunction vf, IntervalVector dx, Float hmax)\n");
    Float hmax=h;
    IntervalVector bx;
    make_lpair(h,bx)=this->flow_bounds(vf,dx,hmax);
    while(true) {
        try {
            return this->flow_step(vf,dx,h,bx);
        } catch(FlowTimeStepException e) {
            h/=2;
        }
    }
}

IntervalVectorFunctionModel
TaylorPicardIntegrator::flow_step(const IntervalVectorFunction& f, const IntervalVector& dx, const Float& h, const IntervalVector& bx) const
{
    ARIADNE_LOG(3,"TaylorPicardIntegrator::flow_step(IntervalVectorFunction vf, IntervalVector dx, Float h, IntervalVector bx)\n");
    ARIADNE_LOG(3," dx="<<dx<<" h="<<h<<" bx="<<bx<<"\n");
    const uint nx=dx.size();
    Sweeper sweeper(new ThresholdSweeper(this->_step_sweep_threshold));

    IntervalVector dom=join(dx,Interval(-h,h));
    ARIADNE_LOG(7,"dom="<<dom<<"\n");

    IntervalVectorFunctionModel phi0=this->function_factory().create_zeros(nx,dom);
    for(uint i=0; i!=nx; ++i) { phi0[i]=this->function_factory().create_coordinate(dom,i); }
    ARIADNE_LOG(5,"phi0="<<phi0<<"\n");

    IntervalVectorFunctionModel phi=this->function_factory().create_zeros(nx,dom);
    for(uint i=0; i!=nx; ++i) { phi[i]=this->function_factory().create_constant(dom,bx[i]); }

    ARIADNE_LOG(5,"phi="<<phi<<"\n");
    for(uint k=0; k!=this->_maximum_temporal_order; ++k) {
        bool last_step=(phi.error()<this->maximum_error());
        IntervalVectorFunctionModel fphi=compose(f,phi);
        ARIADNE_LOG(5,"fphi="<<fphi<<"\n");
        for(uint i=0; i!=nx; ++i) {
            phi[i]=antiderivative(fphi[i],nx)+phi0[i];
        }
        ARIADNE_LOG(4,"phi="<<phi<<"\n");
        if(last_step) { break; }
    }

    if(phi.error()>this->step_maximum_error()) {
        ARIADNE_THROW(FlowTimeStepException,"TaylorPicardIntegrator::flow_step","Integration of "<<f<<" starting in "<<dx<<" for time "<<h<<" has error "<<phi.error()<<" after "<<this->_maximum_temporal_order<<" iterations, which exceeds maximum error "<<this->maximum_error()<<"\n");
    }

    IntervalVectorFunctionModel res=this->function_factory().create_zeros(nx,dom);
    ARIADNE_LOG(4,"res_init="<<res<<"\n");
    for(uint i=0; i!=nx; ++i) { res[i]=phi[i]; }
    //res.sweep();
    ARIADNE_LOG(4,"res="<<res<<"\n");
    return res;

}

void TaylorPicardIntegrator::write(std::ostream& os) const {
    os << "TaylorPicardIntegrator"
       << "(maximum_error = " << this->maximum_error()
       << ", function_factory = " << this->function_factory()
       << ", lipschitz_tolerance = " << this->lipschitz_tolerance()
       << ", step_maximum_error = " << this->step_maximum_error()
       << ", step_sweep_threshold = " << this->step_sweep_threshold()
       << ", maximum_temporal_order = " << this->maximum_temporal_order()
       << " )";
}

} // namespace Ariadne
#include "graded.h"
#include "procedure.h"
#include <include/taylor_function.h>
namespace Ariadne {


template<class F> Graded< Differential<Interval> > flow(const F& f, const Interval& c, Nat M, Nat N) {
    Differential<Interval> x=make_differential_variable(1u,M,c,0u);
    Graded< Differential<Interval> > y=make_graded(x);
    Graded< Differential<Interval> > t=create_graded(x);

    for(Nat n=0; n!=N; ++n) {
        t=f(y);
        y=antidifferential(t);
    }

    return y;
}

Vector< Graded< Differential<Interval> > > flow(const Vector< Procedure<Interval> >& f, const Vector<Interval>& c, Nat M, Nat N) {
    Graded< Differential<Interval> > null;
    Vector< Graded< Differential<Interval> > > y(f.result_size(),null);
    Vector< Graded< Differential<Interval> > > fy(f.result_size(),null);
    List< Graded< Differential<Interval> > > t(f.temporaries_size(),null);
    for(Nat i=0; i!=y.size(); ++i) {
        y[i]=Graded< Differential<Interval> >(Differential<Interval>::variable(y.size(),M,c[i],i));
    }

    for(Nat n=0; n!=N; ++n) {
        Ariadne::compute(f,fy,t,y);
        for(Nat i=0; i!=y.size(); ++i) {
            y[i]=antidifferential(fy[i]);
        }
    }

    return y;
}

class FormulaFunction;
typedef Procedure<Interval> IntervalProcedure;
typedef Differential<Interval> IntervalDifferential;
typedef Polynomial<Interval> IntervalPolynomial;
typedef Graded<IntervalDifferential> GradedIntervalDifferential;
bool operator<(const MultiIndex& a1, const MultiIndex& a2);

template<class X> inline void append_join(Expansion<X>& e, const MultiIndex& a1, const uint a2, const X& c) {
    MultiIndex a(a1.size()+1);
    for(uint i=0; i!=a1.size(); ++i) { a[i]=a1[i]; }
    a[a1.size()]=a2;
    e.append(a,c);
}

inline Vector<GradedIntervalDifferential> graded_variables(int so, const IntervalVector& x) {
    Vector<GradedIntervalDifferential> r(x.size(),GradedIntervalDifferential());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=GradedIntervalDifferential(IntervalDifferential::variable(x.size(),so,x[i],i));
    }
    return r;
}


Vector<IntervalFormula> formula(const IntervalVectorFunction& f) {
    return f.evaluate(IntervalFormula::identity(f.argument_size()));
}

void flow_init(const Vector<IntervalProcedure>& p,
               Vector<GradedIntervalDifferential>& fy, List<GradedIntervalDifferential>& t, Vector<GradedIntervalDifferential>& y,
               const Vector<Interval>& x, const Vector<Interval>& r,  uint so)
{
    Graded< Differential<Interval> > null;
    y=Vector< Graded< Differential<Interval> > >(p.result_size(),null);
    fy=Vector< Graded< Differential<Interval> > >(p.result_size(),null);
    t=List< Graded< Differential<Interval> > >(p.temporaries_size(),null);
    for(Nat i=0; i!=y.size(); ++i) {
        y[i]=Graded< Differential<Interval> >(Differential<Interval>::variable(y.size(),so,Interval(0.0),i)*r[i]+x[i]);
    }
}

void flow_iterate(const Vector<IntervalProcedure>& p, Float h,
                  Vector<GradedIntervalDifferential>& fy, List<GradedIntervalDifferential>& t, Vector<GradedIntervalDifferential>& y)
{
    Ariadne::compute(p,fy,t,y);
    for(Nat i=0; i!=y.size(); ++i) {
        y[i]=antidifferential(fy[i]);
        y[i]*=h;
    }
}


Vector<IntervalDifferential> flow_differential(Vector<GradedIntervalDifferential> const& dphia, Vector<GradedIntervalDifferential> const& dphib,
                                               Vector<GradedIntervalDifferential> const& dphic, Vector<GradedIntervalDifferential> const& dphid,
                                               uint so, uint to, int verbosity=0)
{
    uint nx=dphia.size();
    Vector<GradedIntervalDifferential> gdphi(nx,GradedIntervalDifferential(List<IntervalDifferential>(to+1u,IntervalDifferential(nx,so))));
    for(uint i=0; i!=nx; ++i) {
        for(uint j=0; j!=to; ++j) {
            for(Differential<Interval>::const_iterator iter=dphic[i][j].begin(); iter!=dphic[i][j].end(); ++iter) {
                if(iter->key().degree()<so) { gdphi[i][j].expansion().append(iter->key(),iter->data()); }
            }
            for(Differential<Interval>::const_iterator iter=dphid[i][j].begin(); iter!=dphid[i][j].end(); ++iter) {
                if(iter->key().degree()==so) { gdphi[i][j].expansion().append(iter->key(),iter->data()); }
            }
        }
        uint j=to;
        for(Differential<Interval>::const_iterator iter=dphia[i][j].begin(); iter!=dphia[i][j].end(); ++iter) {
            if(iter->key().degree()<so) { gdphi[i][j].expansion().append(iter->key(),iter->data()); }
        }
        for(Differential<Interval>::const_iterator iter=dphib[i][j].begin(); iter!=dphib[i][j].end(); ++iter) {
            if(iter->key().degree()==so) { gdphi[i][j].expansion().append(iter->key(),iter->data()); }
        }
    }
    ARIADNE_LOG(4,"gdphi="<<gdphi<<"\n");

    Vector<IntervalDifferential> dphi(nx,IntervalDifferential(nx+1u,so+to));
    for(uint i=0; i!=nx; ++i) {
        Expansion<Interval>& component=dphi[i].expansion();
        for(uint j=0; j<=to; ++j) {
            const Expansion<Interval>& expansion=gdphi[i][j].expansion();
            for(Expansion<Interval>::const_iterator iter=expansion.begin(); iter!=expansion.end(); ++iter) {
                append_join(component,iter->key(),j,iter->data());
            }
        }
    }

    ARIADNE_LOG(4,"dphi="<<dphi<<"\n");
    return dphi;
}

VectorTaylorFunction flow_function(const Vector<IntervalDifferential>& dphi, const IntervalVector& dx, const Float& h, double swpt, int verbosity=0) {
    const uint n=dphi.size();
    Sweeper sweeper(new ThresholdSweeper(swpt));
    VectorTaylorFunction tphi(n,join(dx,Interval(-h,+h)),sweeper);

    for(uint i=0; i!=n; ++i) {
        IntervalTaylorModel& model=tphi.model(i);
        Expansion<Float>& expansion=model.expansion();
        Float& error=model.error();
        error=0.0;
        expansion.reserve(dphi[i].expansion().number_of_nonzeros());

        Differential<Interval>::const_iterator iter=dphi[i].begin();
        while(iter!=dphi[i].end()) {
            MultiIndex const a=iter->key();
            Interval ivl=iter->data();
            Float x=ivl.midpoint();
            set_rounding_upward();
            Float e=max(sub_rnd(ivl.upper(),x),sub_rnd(x,ivl.lower()));
            error=add_rnd(error,e);
            set_rounding_to_nearest();
            expansion.append(a,x);
            ++iter;
        }
        model.unique_sort();
        model.sweep();
    }
    return tphi;
}

IntervalVectorFunctionModel
differential_flow_step(const IntervalVectorFunction& f, const IntervalVector& dx, const Float& h, const IntervalVector& bx,
                       double swpt, uint so, uint to, uint verbosity=0)
{
    uint n=f.result_size();
    Vector<IntervalDifferential> idc(n,n+1,so+to);
    Vector<IntervalDifferential> idb(n,n+1,so+to);
    Vector<IntervalDifferential> dphic(n,n+1,so+to);
    Vector<IntervalDifferential> dphib(n,n+1,so+to);
    for(uint i=0; i!=n; ++i) {
        idc[i]=IntervalDifferential::variable(n+1,so+to,0.0,i)*Interval(dx[i].radius())+Interval(midpoint(dx[i]));
        idb[i]=IntervalDifferential::variable(n+1,so+to,0.0,i)*Interval(dx[i].radius())+Interval(bx[i]);
        dphic[i]=idc[i];
        dphib[i]=idb[i];
    }
    for(uint i=0; i!=so+to; ++i) {
        dphic=antiderivative(f(dphic),n)*h+idc;
        dphib=antiderivative(f(dphib),n)*h+idb;
    }

    VectorTaylorFunction tphi(n,join(dx,Interval(-h,+h)),ThresholdSweeper(swpt));
    for(uint i=0; i!=n; ++i) {
        IntervalTaylorModel& model=tphi.model(i);
        Expansion<Float>& expansion=model.expansion();
        Float& error=model.error();
        error=0.0;
        expansion.reserve(dphic[i].expansion().number_of_nonzeros());

        Differential<Interval>::const_iterator citer=dphic[i].begin();
        Differential<Interval>::const_iterator biter=dphib[i].begin();
        while(citer!=dphic[i].end() && biter!=dphib[i].end()) {
            assert(citer->key()==biter->key());
            MultiIndex const a=citer->key();
            Interval ivl;
            if (a.degree()==so+to) {
                ivl=biter->data();
            } else {
                ivl=citer->data();
            }
            Float x=ivl.midpoint();
            set_rounding_upward();
            Float e=max(sub_rnd(ivl.upper(),x),sub_rnd(x,ivl.lower()));
            error=add_rnd(error,e);
            set_rounding_to_nearest();
            expansion.append(a,x);
            ++citer;
            ++biter;
        }
        model.unique_sort();
        model.sweep();
    }
    return tphi;
}

IntervalVectorFunctionModel
differential_space_time_flow_step(const IntervalVectorFunction& f, const IntervalVector& dx, const Float& h, const IntervalVector& bx,
                                  double swpt, uint so, uint to, uint verbosity=0)
{
    uint n=f.result_size();
    Vector<IntervalDifferential> idc(n,n+1,so+to);
    Vector<IntervalDifferential> idb(n,n+1,so+to);
    Vector<IntervalDifferential> dphic(n,n+1,so+to);
    Vector<IntervalDifferential> dphib(n,n+1,so+to);
    for(uint i=0; i!=n; ++i) {
        idc[i]=IntervalDifferential::variable(n+1,so+to,0.0,i)*Interval(dx[i].radius())+Interval(midpoint(dx[i]));
        idb[i]=IntervalDifferential::variable(n+1,so+to,0.0,i)*Interval(dx[i].radius())+Interval(bx[i]);
        dphic[i]=idc[i];
        dphib[i]=idb[i];
    }
    for(uint i=0; i!=so+to; ++i) {
        dphic=antiderivative(f(dphic),n)*h+idc;
        dphib=antiderivative(f(dphib),n)*h+idb;
    }

    VectorTaylorFunction tphi(n,join(dx,Interval(-h,+h)),ThresholdSweeper(swpt));
    for(uint i=0; i!=n; ++i) {
        IntervalTaylorModel& model=tphi.model(i);
        Expansion<Float>& expansion=model.expansion();
        Float& error=model.error();
        error=0.0;
        expansion.reserve(dphic[i].expansion().number_of_nonzeros());

        Differential<Interval>::const_iterator citer=dphic[i].begin();
        Differential<Interval>::const_iterator biter=dphib[i].begin();
        while(citer!=dphic[i].end() && biter!=dphib[i].end()) {
            assert(citer->key()==biter->key());
            MultiIndex const a=citer->key();
            Interval ivl;
            if (a[n]<=to && a.degree()<=so+a[n]) {
                if(a[n]<to && a.degree()<so+a[n]) {
                    ivl=citer->data();
                } else {
                    ivl=biter->data();
                }
                Float x=ivl.midpoint();
                set_rounding_upward();
                Float e=max(sub_rnd(ivl.upper(),x),sub_rnd(x,ivl.lower()));
                error=add_rnd(error,e);
                set_rounding_to_nearest();
                expansion.append(a,x);
            }
            ++citer;
            ++biter;
        }
        model.unique_sort();
        model.sweep();
    }
    return tphi;
}

IntervalVectorFunctionModel
series_flow_step(const IntervalVectorFunction& f, const IntervalVector& dx, const Float& h, const IntervalVector& bx,
                 double max_err, double swpt, uint init_so, uint init_to, uint max_so, uint max_to, uint verbosity)
{
    static const double TRY_SPACIAL_ORDER_INCREASE_FACTOR=4;

    Vector<IntervalFormula> ff = formula(f);
    Vector< Procedure<Interval> > p(ff);
    ARIADNE_LOG(4,"p="<<p<<"\n");

    IntervalVector cx=midpoint(dx);
    IntervalVector ax=cx+Interval(0,h)*evaluate(p,bx);
    ax=cx+Interval(0,h)*evaluate(p,ax);

    uint so=init_so;
    uint to=init_to;

    uint nso=0;
    uint nto=0;

    uint n=dx.size();
    Vector<Interval> rdx(n);
    for(uint i=0; i!=n; ++i) { rdx[i]=rad_ivl(dx[i].lower(),dx[i].upper()); }

    Vector<GradedIntervalDifferential> dphia,fdphia,dphib,fdphib,dphic,fdphic,dphid,fdphid;
    List<GradedIntervalDifferential> tdphia,tdphib,tdphic,tdphid;
    Vector<GradedIntervalDifferential> ndphia,nfdphia,ndphib,nfdphib,ndphic,nfdphic,ndphid,nfdphid;
    List<GradedIntervalDifferential> ntdphia,ntdphib,ntdphic,ntdphid;

    flow_init(p,fdphia,tdphia,dphia,ax,rdx,so);
    flow_init(p,fdphib,tdphib,dphib,bx,rdx,so);
    flow_init(p,fdphic,tdphic,dphic,cx,rdx,so);
    flow_init(p,fdphid,tdphid,dphid,dx,rdx,so);

    for(uint i=0; i!=to; ++i) {
        Ariadne::flow_iterate(p,h,fdphia,tdphia,dphia);
        Ariadne::flow_iterate(p,h,fdphib,tdphib,dphib);
        Ariadne::flow_iterate(p,h,fdphic,tdphic,dphic);
        Ariadne::flow_iterate(p,h,fdphid,tdphid,dphid);
    }

    Vector<IntervalDifferential> dphi=flow_differential(dphia,dphib,dphic,dphid,so,to,verbosity);
    ARIADNE_LOG(5,"dphi="<<dphi<<"\n");

    VectorTaylorFunction tphi=flow_function(dphi,dx,h,swpt,verbosity);
    ARIADNE_LOG(5,"phi="<<tphi<<"\n");

    Float old_error=tphi.error()*TRY_SPACIAL_ORDER_INCREASE_FACTOR*2;

    while(tphi.error()>max_err && (so<max_so || to<max_to) ) {
        uint nnz=0; for(uint i=0; i!=tphi.size(); ++i) { nnz+=tphi.model(i).number_of_nonzeros(); }
        ARIADNE_LOG(3,"so="<<so<<" to="<<to<<" nnz="<<nnz<<" err="<<tphi.error()<<"\n");

        if( (so<max_so) && (tphi.error()*TRY_SPACIAL_ORDER_INCREASE_FACTOR > old_error) ) {
            // try increasing spacial degree
            if(nto==0) {
                // Initialise higher spacial order
                nso=so+1;
                nto=to-1;

                flow_init(p,nfdphia,ntdphia,ndphia,ax,rdx,nso);
                flow_init(p,nfdphib,ntdphib,ndphib,bx,rdx,nso);
                flow_init(p,nfdphic,ntdphic,ndphic,cx,rdx,nso);
                flow_init(p,nfdphid,ntdphid,ndphid,dx,rdx,nso);

                for(uint i=0; i!=nto; ++i) {
                    Ariadne::flow_iterate(p,h,nfdphia,ntdphia,ndphia);
                    Ariadne::flow_iterate(p,h,nfdphib,ntdphib,ndphib);
                    Ariadne::flow_iterate(p,h,nfdphic,ntdphic,ndphic);
                    Ariadne::flow_iterate(p,h,nfdphid,ntdphid,ndphid);
                }
            }
            while(nto+1<to) {
                ++nto;
                Ariadne::flow_iterate(p,h,nfdphia,ntdphia,ndphia);
                Ariadne::flow_iterate(p,h,nfdphib,ntdphib,ndphib);
                Ariadne::flow_iterate(p,h,nfdphic,ntdphic,ndphic);
                Ariadne::flow_iterate(p,h,nfdphid,ntdphid,ndphid);
            }
            Vector<IntervalDifferential> ndphi=flow_differential(ndphia,ndphib,ndphic,ndphid,nso,nto,verbosity);
            VectorTaylorFunction ntphi=flow_function(ndphi,dx,h,swpt,verbosity);

            uint nnnz=0; for(uint i=0; i!=tphi.size(); ++i) { nnnz+=tphi.model(i).number_of_nonzeros(); }
            ARIADNE_LOG(3,"nso="<<nso<<" nto="<<nto<<" nnnz="<<nnnz<<" nerr="<<ntphi.error()<<"\n");

            if( to==max_to || ntphi.error()<tphi.error()) {
                dphia=ndphia; dphib=ndphib; dphic=ndphic; dphid=ndphid;
                fdphia=nfdphia; fdphib=nfdphib; fdphic=nfdphic; fdphid=nfdphid;
                tdphia=ntdphia; tdphib=ntdphib; tdphic=ntdphic; tdphid=ntdphid;
                dphi=ndphi;
                tphi=ntphi;
                so=nso;
                to=nto;
                nso=0;
                nto=0;
            }
        }

        old_error=tphi.error();

        ++to;
        Ariadne::flow_iterate(p,h,fdphia,tdphia,dphia);
        Ariadne::flow_iterate(p,h,fdphib,tdphib,dphib);
        Ariadne::flow_iterate(p,h,fdphic,tdphic,dphic);
        Ariadne::flow_iterate(p,h,fdphid,tdphid,dphid);
        dphi=flow_differential(dphia,dphib,dphic,dphid,so,to,verbosity);
        tphi=flow_function(dphi,dx,h,swpt,verbosity);
    }
    uint nnz=0; for(uint i=0; i!=tphi.size(); ++i) { nnz+=tphi.model(i).number_of_nonzeros(); }
    ARIADNE_LOG(2,"so="<<so<<" to="<<to<<" nnz="<<nnz<<" err="<<tphi.error()<<"\n");
    ARIADNE_LOG(4,"phi="<<tphi<<"\n");
    return tphi;
}

IntervalVectorFunctionModel
TaylorSeriesIntegrator::flow_step(const IntervalVectorFunction& f, const IntervalVector& dx, const Float& h, const IntervalVector& bx) const
{
    IntervalVectorFunctionModel tphi=Ariadne::series_flow_step(f,dx,h,bx,
        this->step_maximum_error(),this->step_sweep_threshold(),
        this->minimum_spacial_order(),this->minimum_temporal_order(),
        this->maximum_spacial_order(),this->maximum_temporal_order(),this->verbosity);

/*
    IntervalVectorFunctionModel tphi=Ariadne::differential_space_time_flow_step(f,dx,h,bx,
        this->step_sweep_threshold(),6,4,this->verbosity);
*/

    if(tphi.error()>this->step_maximum_error()) {
        ARIADNE_THROW(FlowTimeStepException,"TaylorSeriesIntegrator::flow_step",
                      "Integration of "<<f<<" over "<<dx<<" for time "<<h<<" has error "<<tphi.errors()<<
                      " using spacial order "<<this->maximum_spacial_order()<<" and temporal order "<<this->maximum_temporal_order()<<
                      ", which exceeds maximum single-step error "<<this->step_maximum_error()<<"\n");
    }

    return tphi;
}

Pair<Float,IntervalVector>
TaylorSeriesIntegrator::flow_bounds(const IntervalVectorFunction& vf, const IntervalVector& dx, const Float& hmax) const
{
    ARIADNE_LOG(3,"TaylorSeriesIntegrator::flow_bounds(IntervalVectorFunction vf, IntervalVector dx, Float hmax)\n");
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

    Float hmin=hmax/(1<<REDUCTION_STEPS);
    Float h=hmax;
    ARIADNE_LOG(4,"hmax="<<hmax<<"\n");

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

void TaylorSeriesIntegrator::write(std::ostream& os) const {
    os << "TaylorSeriesIntegrator"
       << "( function_factory = " << this->function_factory()
       << ", maximum_error = " << this->maximum_error()
       << ", lipschitz_tolerance = " << this->lipschitz_tolerance()
       << ", step_maximum_error = " << this->step_maximum_error()
       << ", step_sweep_threshold = " << this->step_sweep_threshold()
       << ", minimum_spacial_order = " << this->minimum_spacial_order()
       << ", minimum_temporal_order = " << this->minimum_temporal_order()
       << ", maximum_temporal_order = " << this->maximum_temporal_order()
       << ", maximum_spacial_order = " << this->maximum_spacial_order()
       << " )";
}


template<class X> void truncate(Differential<X>& x, uint spacial_order, uint temporal_order) {
    uint n=x.argument_size()-1;
    typename Differential<X>::iterator write_iter=x.begin();
    typename Differential<X>::const_iterator read_iter=x.begin();
    while(read_iter!=x.end()) {
        const MultiIndex& index = read_iter->key();
        if(index[n]>temporal_order || index[n]+spacial_order<index.degree()) {
        } else {
            *write_iter=*read_iter;
            ++write_iter;
        }
        ++read_iter;
    }
    x.expansion().resize(write_iter-x.begin());
}

template<class X> void truncate(Vector< Differential<X> >& x, uint spacial_order, uint temporal_order) {
    for(uint i=0; i!=x.size(); ++i) { truncate(x[i],spacial_order,temporal_order); }
}


IntervalDifferentialVector
AffineIntegrator::flow_derivative(const IntervalVectorFunction& f, const IntervalVector& dom) const
{
    IntervalDifferentialVector dx=IntervalDifferential::variables(this->_spacial_order+this->_temporal_order,join(dom,Interval(0.0)));
    dx[dom.size()]=0.0;
    IntervalDifferentialVector dphi = dx;
    for(uint i=0; i!=_temporal_order; ++i) {
        dphi = antiderivative(f.evaluate(dphi),dom.size())+dx;
    }
    truncate(dphi,this->_spacial_order,this->_temporal_order);
    return dphi;
}

IntervalVectorFunctionModel
AffineIntegrator::flow_step(const IntervalVectorFunction& f, const IntervalVector& dom, const Float& h, const IntervalVector& bbox) const
{
    IntervalVector mid = IntervalVector(midpoint(dom));

    IntervalDifferentialVector mdphi = this->flow_derivative(f,mid);
    IntervalDifferentialVector bdphi = this->flow_derivative(f,bbox);

    ARIADNE_WARN("AffineIntegrator may compute overly optimistic error bounds.");

    const uint n=dom.size();
    Vector<Float> err(n);

    set_rounding_upward();
    Vector<Float> rad(n+1);
    for(uint i=0; i!=n; ++i) {
        rad[i] = max(dom[i].upper()-mid[i].lower(),mid[i].upper()-dom[i].lower());
    }
    rad[n] = h;

    for(uint i=0; i!=n; ++i) {
        for(Expansion<Interval>::const_iterator iter=bdphi[i].begin(); iter!=bdphi[i].end(); ++iter) {
            const MultiIndex& a=iter->key();
            if(a[n]==this->_temporal_order && a[n]+this->_spacial_order==a.degree()) {
                const Interval& rng = iter->data();
                const Interval& mid = mdphi[i][a];
                ARIADNE_ASSERT(rng.lower()<=mid.lower() && mid.upper()<=rng.upper());
                Float mag = max(rng.upper()-mid.upper(),mid.lower()-rng.lower());
                for(uint j=0; j!=n+1; ++j) { mag *= pow(rad[j],a[j]); }
                err[i] += mag;
            }
        }
    }
    set_rounding_to_nearest();

    IntervalVector flow_domain = join(dom,Interval(0,h));

    IntervalVectorFunctionModel id = this->function_factory().create_identity(flow_domain);
    IntervalVectorFunctionModel res = this->function_factory().create_zeros(n,flow_domain);
    for(uint i=0; i!=n; ++i) {
        IntervalScalarFunctionModel res_model = res[i] + mdphi[i].expansion()[MultiIndex::zero(n+1)];
        for(uint j=0; j!=mdphi[i].argument_size()-1; ++j) { res_model+=mdphi[i].expansion()[MultiIndex::unit(n+1,j)]*(id[j]-Interval(midpoint(flow_domain[j]))); }
        uint j=mdphi[i].argument_size()-1; { res_model+=mdphi[i].expansion()[MultiIndex::unit(n+1,j)]*id[j]; }
        res_model += Interval(-err[i],+err[i]);
        res[i]=res_model;
    }
    return res;
}

void AffineIntegrator::write(std::ostream& os) const {
    os << "AffineIntegrator"
       << "( function_factory = " << this->function_factory()
       << ", maximum_error = " << this->maximum_error()
       << ", lipschitz_tolerance = " << this->lipschitz_tolerance()
       << ", spacial_order = " << this->spacial_order()
       << ", temporal_order = " << this->temporal_order()
       << " )";
}



} // namespace Ariadne
