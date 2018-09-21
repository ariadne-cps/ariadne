/***************************************************************************
 *            integrator.cpp
 *
 *  Copyright  2006-10  Pieter Collins
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
#include "../config.hpp"

#include <iomanip>

#include "../solvers/integrator.hpp"
#include "../solvers/bounder.hpp"

#include "../output/logging.hpp"
#include "../utility/container.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/differential.hpp"
#include "../algebra/sweeper.hpp"
#include "../algebra/algebra.hpp"
#include "../function/function.hpp"
#include "../function/function_model.hpp"
#include "../function/formula.hpp"
#include "../function/taylor_function.hpp"
#include "../function/taylor_model.hpp"

#include "../function/polynomial.hpp"
#include "../geometry/interval.hpp"

#include "../algebra/expansion.inl.hpp"

namespace Ariadne {

static const FloatDPValue zero=FloatDPValue(0,dp);

inline UpperBoxType operator+(Vector<ExactIntervalType> bx, Vector<UpperIntervalType> const& ex) {
    return Vector<UpperIntervalType>(bx) + ex;
}

inline UpperBoxType operator+(Vector<UpperIntervalType> bx, Vector<FloatDPBounds> const& v) {
    return bx + Vector<UpperIntervalType>(v);
}

inline UpperBoxType operator+(Vector<ExactIntervalType> bx, Vector<FloatDPBounds> const& v) {
    return Vector<UpperIntervalType>(bx) + Vector<UpperIntervalType>(v);
}


IntegratorBase::IntegratorBase(MaximumError e, LipschitzConstant l)
    :  _maximum_error(e), _lipschitz_tolerance(l), _maximum_step_size(16), _function_factory_ptr(make_taylor_function_factory())
{
    ARIADNE_PRECONDITION(e>0.0);
    ARIADNE_PRECONDITION(l>0.0)
}

IntegratorBase::IntegratorBase(MaximumError e, SweepThreshold s, LipschitzConstant l)
    :  _maximum_error(e), _lipschitz_tolerance(l), _maximum_step_size(16), _function_factory_ptr(make_taylor_function_factory(s))
{
    ARIADNE_PRECONDITION(e>0.0);
    ARIADNE_PRECONDITION(l>0.0);
}

Void
IntegratorBase::set_function_factory(const ValidatedFunctionModelDPFactoryInterface& factory)
{
    this->_function_factory_ptr=FunctionFactoryPointer(factory.clone());
}

const ValidatedFunctionModelDPFactoryInterface&
IntegratorBase::function_factory() const
{
    return *this->_function_factory_ptr;
}

Pair<FloatDPValue,UpperBoxType>
IntegratorBase::flow_bounds(const ValidatedVectorFunction& vf, const ExactBoxType& domx, const RawFloatDP& hsug) const {
    return EulerBounder().flow_bounds(vf,domx,PositiveFloatDPApproximation(hsug));
}





ValidatedVectorFunctionModelDP
IntegratorBase::flow_to(const ValidatedVectorFunction& vf, const ExactBoxType& dx0, const Real& tmax) const
{
    ARIADNE_LOG(1,"IntegratorBase::flow_to(ValidatedVectorFunction vf, ExactBoxType dx0, Real tmax)\n");
    ARIADNE_LOG(2,"vf="<<vf<<"\n");
    ARIADNE_LOG(2,"dom(x0)="<<dx0<<" tmax="<<tmax<<"\n");
    const Nat n=dx0.size(); // Dimension of the state space
    ValidatedVectorFunctionModelDP flow_function=this->function_factory().create_identity(dx0);
    Rational t=0;
    ValidatedVectorFunctionModelDP step_function;
    while(possibly(t<tmax)) {
        ExactBoxType dx=flow_function.codomain();
        FloatDPBounds h_max=FloatDPBounds(tmax,dp)-t;
        FloatDPValue h;
        UpperBoxType bx;
        make_lpair(h,bx) = this->flow_bounds(vf,dx,h_max.upper().raw());
        Bool flow_successfully_computed=false;
        while(!flow_successfully_computed) {
            try {
                step_function=this->flow_step(vf,dx,h,bx);
                flow_successfully_computed=true;
            } catch(const FlowTimeStepException& e) {
                h=hlf(h);
            }
        }
        step_function=partial_evaluate(step_function,n,h);
        flow_function=compose(step_function,flow_function);
        t=t+Rational(h);
    }
    return flow_function;
}


List<ValidatedVectorFunctionModelDP>
IntegratorBase::flow(const ValidatedVectorFunction& vf, const ExactBoxType& dx0, const Real& tmin, const Real& tmax) const
{
    ARIADNE_LOG(1,"IntegratorBase::flow(ValidatedVectorFunction vf, ExactBoxType dx0, Real tmin, Real tmax)\n");
    DoublePrecision precision;
    FloatDPLowerBound tminl=FloatDPBounds(tmin,precision).lower();
    FloatDPUpperBound tmaxu=FloatDPBounds(tmax,precision).upper();
    ValidatedVectorFunctionModelDP evolve_function=this->flow_to(vf,dx0,tmin);
    FloatDPValue t=cast_exact(tminl);
    List<ValidatedVectorFunctionModelDP> result;

    while(possibly(t<tmax)) {
        ExactBoxType dx=evolve_function.codomain();
        FloatDPValue h=cast_exact(tmaxu-t);
        UpperBoxType bx;
        make_lpair(h,bx) = this->flow_bounds(vf,dx,h.raw());
        ValidatedVectorFunctionModelDP flow_step_function=this->flow_step(vf,dx,h,bx);
        FloatDPValue new_t=cast_exact((t+h).lower());
        ExactIntervalType dt(t,new_t);
        ValidatedScalarFunctionModelDP step_time_function=this->function_factory().create_identity(dt)-FloatDPValue(t);
        ValidatedVectorFunctionModelDP flow_function=compose(flow_step_function,combine(evolve_function,step_time_function));
        ARIADNE_ASSERT(flow_function.domain()[dx0.size()].upper()==new_t);
        result.append(flow_function);
        evolve_function=partial_evaluate(flow_function,dx0.size(),FloatDPValue(new_t));
        t=new_t;
    }
    return result;
}

List<ValidatedVectorFunctionModelDP>
IntegratorBase::flow(const ValidatedVectorFunction& vf, const ExactBoxType& dx0, const Real& tmax) const
{
    return flow(vf,dx0,Real(0),tmax);
}



ValidatedVectorFunctionModelDP
IntegratorBase::flow_step(const ValidatedVectorFunction& vf, const ExactBoxType& dx, FloatDP& hmax) const
{
    ARIADNE_LOG(3,"IntegratorBase::flow_step(ValidatedVectorFunction vf, ExactBoxType dx, FloatDP hmax)\n");
    FloatDPValue& h=reinterpret_cast<FloatDPValue&>(hmax);
    UpperBoxType bx;
    make_lpair(h,bx)=this->flow_bounds(vf,dx,hmax);
    while(true) {
        try {
            return this->flow_step(vf,dx,h,bx);
        } catch(const FlowTimeStepException& e) {
            h=hlf(h);
        }
    }
}

ValidatedVectorFunctionModelDP
TaylorPicardIntegrator::flow_step(const ValidatedVectorFunction& f, const ExactBoxType& dx, const FloatDPValue& h, const UpperBoxType& bx) const
{
    ARIADNE_LOG(3,"TaylorPicardIntegrator::flow_step(ValidatedVectorFunction vf, ExactBoxType dx, FloatDPValue h, UpperBoxType bx)\n");
    ARIADNE_LOG(3," dx="<<dx<<" h="<<h<<" bx="<<bx<<"\n");
    const Nat nx=dx.size();
    Sweeper<FloatDP> sweeper(new ThresholdSweeper<FloatDP>(dp,this->_step_sweep_threshold));

    ExactBoxType dom=join(dx,ExactIntervalType(-h,h));
    ARIADNE_LOG(7,"dom="<<dom<<"\n");

    ValidatedVectorFunctionModelDP phi0=this->function_factory().create_zeros(nx,dom);
    for(Nat i=0; i!=nx; ++i) { phi0[i]=this->function_factory().create_coordinate(dom,i); }
    ARIADNE_LOG(5,"phi0="<<phi0<<"\n");

    ValidatedVectorFunctionModelDP phi=this->function_factory().create_zeros(nx,dom);
    for(Nat i=0; i!=nx; ++i) { phi[i]=this->function_factory().create_constant(dom,cast_singleton(bx[i])); }

    ARIADNE_LOG(5,"phi="<<phi<<"\n");
    for(Nat k=0; k!=this->_maximum_temporal_order; ++k) {
        Bool last_step=(phi.error().raw()<this->maximum_error());
        ValidatedVectorFunctionModelDP fphi=compose(f,phi);
        ARIADNE_LOG(5,"fphi="<<fphi<<"\n");
        for(Nat i=0; i!=nx; ++i) {
            phi[i]=antiderivative(fphi[i],nx)+phi0[i];
        }
        ARIADNE_LOG(4,"phi="<<phi<<"\n");
        if(last_step) { break; }
    }

    if(phi.error().raw()>this->step_maximum_error()) {
        ARIADNE_THROW(FlowTimeStepException,"TaylorPicardIntegrator::flow_step","Integration of "<<f<<" starting in "<<dx<<" for time "<<h<<" has error "<<phi.error()<<" after "<<this->_maximum_temporal_order<<" iterations, which exceeds maximum error "<<this->maximum_error()<<"\n");
    }

    ValidatedVectorFunctionModelDP res=this->function_factory().create_zeros(nx,dom);
    ARIADNE_LOG(4,"res_init="<<res<<"\n");
    for(Nat i=0; i!=nx; ++i) { res[i]=phi[i]; }
    //res.sweep();
    ARIADNE_LOG(4,"res="<<res<<"\n");
    return res;

}

Void TaylorPicardIntegrator::write(OutputStream& os) const {
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


#include "../algebra/graded.hpp"
#include "../function/procedure.hpp"
namespace Ariadne {

class FormulaFunction;
typedef Procedure<ValidatedNumber> ValidatedProcedure;
typedef Differential<FloatDPBounds> ValidatedDifferential;
typedef Polynomial<FloatDPBounds> ValidatedPolynomial;
typedef Graded<ValidatedDifferential> GradedValidatedDifferential;
Bool operator<(const MultiIndex& a1, const MultiIndex& a2);


TaylorSeriesIntegrator::TaylorSeriesIntegrator(MaximumError err)
    : TaylorSeriesIntegrator(err,SweepThreshold(err/1024),LipschitzConstant(0.5))
{ }

TaylorSeriesIntegrator::TaylorSeriesIntegrator(MaximumError err, SweepThreshold swp, LipschitzConstant lip)
    : TaylorSeriesIntegrator(err,swp,lip,StepMaximumError(err/128),StepSweepThreshold(swp/1024),MaximumTemporalOrder(12))
{ }

TaylorSeriesIntegrator::TaylorSeriesIntegrator(MaximumError err, SweepThreshold swp, LipschitzConstant lip,
                        StepMaximumError stperr, StepSweepThreshold stpswp, MaximumTemporalOrder maxto)
    : TaylorSeriesIntegrator(err,swp,lip,stperr,stpswp,MinimumSpacialOrder(1),MinimumTemporalOrder(4),MaximumSpacialOrder(4),maxto)
{
}

TaylorSeriesIntegrator::TaylorSeriesIntegrator(
        MaximumError err, SweepThreshold swp, LipschitzConstant lip,
        StepMaximumError stperr, StepSweepThreshold stpswp,
        MinimumSpacialOrder minso, MinimumTemporalOrder minto,
        MaximumSpacialOrder maxso, MaximumTemporalOrder maxto)
    : IntegratorBase(err,swp,lip), _step_maximum_error(stperr), _step_sweep_threshold(stpswp)
    , _minimum_spacial_order(minso), _minimum_temporal_order(minto), _maximum_spacial_order(maxso), _maximum_temporal_order(maxto)
{ }


namespace {


template<class F> GradedValidatedDifferential flow(const F& f, const ExactIntervalType& c, Nat M, Nat N) {
    ValidatedDifferential x=make_differential_variable(1u,M,cast_singleton(c),0u);
    GradedValidatedDifferential y=make_graded(x);
    GradedValidatedDifferential t=create_graded(x);

    for(Nat n=0; n!=N; ++n) {
        t=f(y);
        y=antidifferential(t);
    }

    return y;
}

template<class X> Void append_join(Expansion<MultiIndex,X>& e, const MultiIndex& a1, const Nat a2, const X& c) {
    MultiIndex a(a1.size()+1);
    for(Nat i=0; i!=a1.size(); ++i) { a[i]=a1[i]; }
    a[a1.size()]=a2;
    e.append(a,c);
}


inline Vector<ValidatedFormula> formula(const ValidatedVectorFunction& f) {
    return f.evaluate(ValidatedFormula::identity(f.argument_size()));
}

Void flow_init(const Vector<ValidatedProcedure>& p,
               Vector<GradedValidatedDifferential>& fy, List<GradedValidatedDifferential>& t, Vector<GradedValidatedDifferential>& y,
               const Vector<ValidatedNumericType>& x, const Vector<ValidatedNumericType>& r,  Nat so)
{
    GradedValidatedDifferential null;
    y=Vector< GradedValidatedDifferential >(p.result_size(),null);
    fy=Vector< GradedValidatedDifferential >(p.result_size(),null);
    t=List< GradedValidatedDifferential >(p.temporaries_size(),null);
    for(Nat i=0; i!=y.size(); ++i) {
        y[i]=GradedValidatedDifferential(Differential<ValidatedNumericType>::variable(y.size(),so,ValidatedNumericType(0,0),i)*r[i]+x[i]);
    }
}

Void flow_iterate(const Vector<ValidatedProcedure>& p, FloatDPValue h,
                  Vector<GradedValidatedDifferential>& fy, List<GradedValidatedDifferential>& t, Vector<GradedValidatedDifferential>& y)
{
    Ariadne::compute(p,fy,t,y);
    for(Nat i=0; i!=y.size(); ++i) {
        y[i]=antidifferential(fy[i]);
        y[i]*=h;
    }
}


Vector<ValidatedDifferential> flow_differential(Vector<GradedValidatedDifferential> const& dphia, Vector<GradedValidatedDifferential> const& dphib,
                                               Vector<GradedValidatedDifferential> const& dphic, Vector<GradedValidatedDifferential> const& dphid,
                                               Nat so, Nat to, Nat verbosity=0)
{
    Nat nx=dphia.size();
    Vector<GradedValidatedDifferential> gdphi(nx,GradedValidatedDifferential(List<ValidatedDifferential>(to+1u,ValidatedDifferential(nx,so))));
    for(Nat i=0; i!=nx; ++i) {
        for(Nat j=0; j!=to; ++j) {
            for(ValidatedDifferential::ConstIterator iter=dphic[i][j].begin(); iter!=dphic[i][j].end(); ++iter) {
                if(iter->index().degree()<so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
            }
            for(ValidatedDifferential::ConstIterator iter=dphid[i][j].begin(); iter!=dphid[i][j].end(); ++iter) {
                if(iter->index().degree()==so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
            }
        }
        Nat j=to;
        for(ValidatedDifferential::ConstIterator iter=dphia[i][j].begin(); iter!=dphia[i][j].end(); ++iter) {
            if(iter->index().degree()<so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
        }
        for(ValidatedDifferential::ConstIterator iter=dphib[i][j].begin(); iter!=dphib[i][j].end(); ++iter) {
            if(iter->index().degree()==so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
        }
    }
    ARIADNE_LOG(4,"gdphi="<<gdphi<<"\n");

    Vector<ValidatedDifferential> dphi(nx,ValidatedDifferential(nx+1u,so+to));
    for(Nat i=0; i!=nx; ++i) {
        Expansion<MultiIndex,FloatDPBounds>& component=dphi[i].expansion();
        for(Nat j=0; j<=to; ++j) {
            const Expansion<MultiIndex,FloatDPBounds>& expansion=gdphi[i][j].expansion();
            for(Expansion<MultiIndex,FloatDPBounds>::ConstIterator iter=expansion.begin(); iter!=expansion.end(); ++iter) {
                append_join(component,iter->index(),j,iter->coefficient());
            }
        }
    }

    ARIADNE_LOG(4,"dphi="<<dphi<<"\n");
    return dphi;
}

ValidatedVectorTaylorFunctionModelDP flow_function(const Vector<ValidatedDifferential>& dphi, const ExactBoxType& dx, const FloatDPValue& h, double swpt, Nat verbosity=0) {
    const Nat n=dphi.size();
    Sweeper<FloatDP> sweeper(new ThresholdSweeper<FloatDP>(dp,swpt));
    ValidatedVectorTaylorFunctionModelDP tphi(n,join(dx,ExactIntervalType(-h,+h)),sweeper);

    for(Nat i=0; i!=n; ++i) {
        ValidatedTaylorModelDP& model=tphi.model(i);
        Expansion<MultiIndex,FloatDPValue>& expansion=model.expansion();
        FloatDPError& error=model.error();
        error=0u;
        expansion.reserve(dphi[i].expansion().number_of_nonzeros());

        Differential<FloatDPBounds>::ConstIterator iter=dphi[i].begin();
        while(iter!=dphi[i].end()) {
            MultiIndex const a=iter->index();
            FloatDPBounds coef=iter->coefficient();
            FloatDPValue x=coef.value();
            error+=coef.error();
            expansion.append(a,x);
            ++iter;
        }
        model.cleanup();
    }
    return tphi;
}

ValidatedVectorFunctionModelDP
series_flow_step(const ValidatedVectorFunction& f, const ExactBoxType& bdx, const FloatDPValue& h, const UpperBoxType& bbx,
                 double max_err, double swpt, Nat init_so, Nat init_to, Nat max_so, Nat max_to, Nat verbosity)
{
    static const double TRY_SPACIAL_ORDER_INCREASE_FACTOR=4;

    Vector<ValidatedFormula> ff = formula(f);
    Vector<ValidatedProcedure> p(ff);
    ARIADNE_LOG(4,"p="<<p<<"\n");

    Vector<ValidatedNumericType> dx=cast_singleton(bdx);
    Vector<ValidatedNumericType> bx=cast_singleton(bbx);
    Vector<ValidatedNumericType> cx=midpoint(bdx);
    Vector<ValidatedNumericType> ax=cx+ValidatedNumericType(0,h)*evaluate(p,bx);
    ax=cx+ValidatedNumericType(0,h)*evaluate(p,ax);

    Nat so=init_so;
    Nat to=init_to;

    Nat nso=0;
    Nat nto=0;

    Nat n=dx.size();
    Vector<ValidatedNumericType> rdx(n);
    for(Nat i=0; i!=n; ++i) { rdx[i]=bdx[i].radius(); }

    Vector<GradedValidatedDifferential> dphia,fdphia,dphib,fdphib,dphic,fdphic,dphid,fdphid;
    List<GradedValidatedDifferential> tdphia,tdphib,tdphic,tdphid;
    Vector<GradedValidatedDifferential> ndphia,nfdphia,ndphib,nfdphib,ndphic,nfdphic,ndphid,nfdphid;
    List<GradedValidatedDifferential> ntdphia,ntdphib,ntdphic,ntdphid;

    flow_init(p,fdphia,tdphia,dphia,ax,rdx,so);
    flow_init(p,fdphib,tdphib,dphib,bx,rdx,so);
    flow_init(p,fdphic,tdphic,dphic,cx,rdx,so);
    flow_init(p,fdphid,tdphid,dphid,dx,rdx,so);

    for(Nat i=0; i!=to; ++i) {
        Ariadne::flow_iterate(p,h,fdphia,tdphia,dphia);
        Ariadne::flow_iterate(p,h,fdphib,tdphib,dphib);
        Ariadne::flow_iterate(p,h,fdphic,tdphic,dphic);
        Ariadne::flow_iterate(p,h,fdphid,tdphid,dphid);
    }
    ARIADNE_LOG(7,"dphia="<<dphia<<"\n");
    ARIADNE_LOG(7,"dphib="<<dphib<<"\n");
    ARIADNE_LOG(7,"dphic="<<dphic<<"\n");
    ARIADNE_LOG(7,"dphid="<<dphid<<"\n");

    Vector<ValidatedDifferential> dphi=flow_differential(dphia,dphib,dphic,dphid,so,to,verbosity);
    ARIADNE_LOG(5,"dphi="<<dphi<<"\n");

    ValidatedVectorTaylorFunctionModelDP tphi=flow_function(dphi,bdx,h,swpt,verbosity);
    ARIADNE_LOG(5,"phi="<<tphi<<"\n");

    FloatDPError old_error=tphi.error()*FloatDPError(TRY_SPACIAL_ORDER_INCREASE_FACTOR*2);

    while(tphi.error().raw()>max_err && (so<max_so || to<max_to) ) {
        Nat nnz=0; for(Nat i=0; i!=tphi.size(); ++i) { nnz+=tphi.model(i).number_of_nonzeros(); }
        ARIADNE_LOG(3,"so="<<so<<" to="<<to<<" nnz="<<nnz<<" err="<<tphi.error()<<"\n");

        if( (so<max_so) && ((tphi.error()*FloatDPError(TRY_SPACIAL_ORDER_INCREASE_FACTOR)).raw() > old_error.raw()) ) {
            // try increasing spacial degree
            if(nto==0) {
                // Initialise higher spacial order
                nso=so+1;
                nto=to-1;

                flow_init(p,nfdphia,ntdphia,ndphia,ax,rdx,nso);
                flow_init(p,nfdphib,ntdphib,ndphib,bx,rdx,nso);
                flow_init(p,nfdphic,ntdphic,ndphic,cx,rdx,nso);
                flow_init(p,nfdphid,ntdphid,ndphid,dx,rdx,nso);

                for(Nat i=0; i!=nto; ++i) {
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
            Vector<ValidatedDifferential> ndphi=flow_differential(ndphia,ndphib,ndphic,ndphid,nso,nto,verbosity);
            ValidatedVectorTaylorFunctionModelDP ntphi=flow_function(ndphi,bdx,h,swpt,verbosity);

            Nat nnnz=0; for(Nat i=0; i!=tphi.size(); ++i) { nnnz+=tphi.model(i).number_of_nonzeros(); }
            ARIADNE_LOG(3,"nso="<<nso<<" nto="<<nto<<" nnnz="<<nnnz<<" nerr="<<ntphi.error()<<"\n");

            if( to==max_to || ntphi.error().raw()<tphi.error().raw()) {
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
        tphi=flow_function(dphi,bdx,h,swpt,verbosity);
    }
    Nat nnz=0; for(Nat i=0; i!=tphi.size(); ++i) { nnz+=tphi.model(i).number_of_nonzeros(); }
    ARIADNE_LOG(2,"so="<<so<<" to="<<to<<" nnz="<<nnz<<" err="<<tphi.error()<<"\n");
    ARIADNE_LOG(4,"phi="<<tphi<<"\n");
    return tphi;
}

}

ValidatedVectorFunctionModelDP
TaylorSeriesIntegrator::flow_step(const ValidatedVectorFunction& f, const ExactBoxType& dx, const FloatDPValue& h, const UpperBoxType& bx) const
{
    ARIADNE_LOG(3,"TaylorSeriesIntegrator::flow_step(ValidatedVectorFunction f, ExactBoxType dx, FloatDPValue h, const UpperBoxType& bx)\n");
    ValidatedVectorFunctionModelDP tphi=Ariadne::series_flow_step(f,dx,h,bx,
        this->step_maximum_error(),this->step_sweep_threshold(),
        this->minimum_spacial_order(),this->minimum_temporal_order(),
        this->maximum_spacial_order(),this->maximum_temporal_order(),this->verbosity);

//    ValidatedVectorFunctionModelDP tphi=Ariadne::differential_space_time_flow_step(f,dx,h,bx,
//        this->step_sweep_threshold(),6,4,this->verbosity);

    if(tphi.error().raw()>this->step_maximum_error()) {
        ARIADNE_THROW(FlowTimeStepException,"TaylorSeriesIntegrator::flow_step",
                      "Integration of "<<f<<" over "<<dx<<" for time "<<h<<" has error "<<tphi.errors()<<
                      " using spacial order "<<this->maximum_spacial_order()<<" and temporal order "<<this->maximum_temporal_order()<<
                      ", which exceeds maximum single-step error "<<this->step_maximum_error()<<"\n");
    }

    return tphi;
}

Pair<FloatDPValue,UpperBoxType>
TaylorSeriesIntegrator::flow_bounds(const ValidatedVectorFunction& vf, const ExactBoxType& dx, const RawFloatDP& hmax) const
{
    return this->IntegratorBase::flow_bounds(vf,dx,hmax);
}

Void TaylorSeriesIntegrator::write(OutputStream& os) const {
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



template<class X> Void truncate(Differential<X>& x, Nat spacial_order_, Nat temporal_order_) {
    Nat n=x.argument_size()-1;
    typename Differential<X>::Iterator write_iter=x.begin();
    typename Differential<X>::ConstIterator read_iter=x.begin();
    while(read_iter!=x.end()) {
        ConstReferenceType<MultiIndex> index = read_iter->index();
        if(index[n]>temporal_order_ || index[n]+spacial_order_<index.degree()) {
        } else {
            *write_iter=*read_iter;
            ++write_iter;
        }
        ++read_iter;
    }
    x.expansion().resize(static_cast<SizeType>(write_iter-x.begin()));
}

template<class X> Void truncate(Vector< Differential<X> >& x, Nat spacial_order_, Nat temporal_order_) {
    for(Nat i=0; i!=x.size(); ++i) { truncate(x[i],spacial_order_,temporal_order_); }
}

AffineIntegrator::AffineIntegrator(MaximumError maximum_error_, TemporalOrder temporal_order_)
    : IntegratorBase(maximum_error_,lipschitz_constant=0.5), _spacial_order(1u), _temporal_order(temporal_order_) { }

AffineIntegrator::AffineIntegrator(MaximumError maximum_error_, SpacialOrder spacial_order_, TemporalOrder temporal_order_)
    : IntegratorBase(maximum_error_,lipschitz_constant=0.5), _spacial_order(spacial_order_), _temporal_order(temporal_order_) { }

Vector<ValidatedDifferential>
AffineIntegrator::flow_derivative(const ValidatedVectorFunction& f, const Vector<ValidatedNumericType>& dom) const
{
    ARIADNE_LOG(5,"AffineIntegrator::flow_derivative(ValidatedVectorFunction f, ValidatedBoxType dom)\n");
    Vector<ValidatedDifferential> dx=
        ValidatedDifferential::variables(this->_spacial_order+this->_temporal_order,
                                         join(dom,zero));
    dx[dom.size()]=zero;
    Vector<ValidatedDifferential> dphi = dx;
    for(Nat i=0; i!=_temporal_order; ++i) {
        dphi = antiderivative(f.evaluate(dphi),dom.size())+dx;
    }
    truncate(dphi,this->_spacial_order,this->_temporal_order);
    return dphi;
}

ValidatedVectorFunctionModelDP
AffineIntegrator::flow_step(const ValidatedVectorFunction& f, const ExactBoxType& dom, const FloatDPValue& h, const UpperBoxType& bbox) const
{
    ARIADNE_LOG(3,"AffineIntegrator::flow_step(ValidatedVectorFunction f, ExactBoxType dom, FloatDP h, UpperBoxType bbox)\n");
    Vector<ValidatedNumericType> dmid = Vector<ValidatedNumericType>(midpoint(dom));

    Vector<ValidatedDifferential> mdphi = this->flow_derivative(f,dmid);
    Vector<ValidatedDifferential> bdphi = this->flow_derivative(f,cast_singleton(bbox));

    ARIADNE_WARN("AffineIntegrator may compute overly optimistic error bounds.");

    const Nat n=dom.size();
    DoublePrecision prec;
    FloatDPError zero_err(prec);

    Vector<FloatDPError> err(n,zero_err);

    Vector<FloatDPError> rad(n+1,zero_err);
    for(Nat i=0; i!=n; ++i) {
        rad[i] = cast_positive(max(dom[i].upper()-dmid[i].lower(),dmid[i].upper()-dom[i].lower()));
    }
    rad[n] = mag(h);

    for(Nat i=0; i!=n; ++i) {
        for(Expansion<MultiIndex,ValidatedNumericType>::ConstIterator iter=bdphi[i].begin(); iter!=bdphi[i].end(); ++iter) {
            ConstReferenceType<MultiIndex> a=iter->index();
            if(a[n]==this->_temporal_order && a[n]+this->_spacial_order==a.degree()) {
                ConstReferenceType<ValidatedNumericType> rng = iter->coefficient();
                ConstReferenceType<ValidatedNumericType> mid = mdphi[i][a];
                ARIADNE_ASSERT(rng.lower().raw()<=mid.lower().raw() && mid.upper().raw()<=rng.upper().raw());
                FloatDPError mag = FloatDPError(max(rng.upper()-mid.lower(),mid.upper()-rng.lower()));
                for(Nat j=0; j!=n+1; ++j) { mag *= pow(rad[j],Nat(a[j])); }
                err[i] += mag;
            }
        }
    }

    ExactBoxType flow_domain = join(dom,ExactIntervalType(0,h));

    ValidatedVectorFunctionModelDP id = this->function_factory().create_identity(flow_domain);
    ValidatedVectorFunctionModelDP res = this->function_factory().create_zeros(n,flow_domain);
    for(Nat i=0; i!=n; ++i) {
        ValidatedScalarFunctionModelDP res_model = res[i] + mdphi[i].expansion()[MultiIndex::zero(n+1)];
        for(Nat j=0; j!=mdphi[i].argument_size()-1; ++j) { res_model+=mdphi[i].expansion()[MultiIndex::unit(n+1,j)]*(id[j]-ValidatedNumericType(midpoint(flow_domain[j]))); }
        Nat j=mdphi[i].argument_size()-1; { res_model+=mdphi[i].expansion()[MultiIndex::unit(n+1,j)]*id[j]; }
        res_model += FloatDPBounds(-err[i],+err[i]);
        res[i]=res_model;
    }
    return res;
}

Void AffineIntegrator::write(OutputStream& os) const {
    os << "AffineIntegrator"
       << "( function_factory = " << this->function_factory()
       << ", maximum_error = " << this->maximum_error()
       << ", lipschitz_tolerance = " << this->lipschitz_tolerance()
       << ", spacial_order = " << this->spacial_order()
       << ", temporal_order = " << this->temporal_order()
       << " )";
}



} // namespace Ariadne
