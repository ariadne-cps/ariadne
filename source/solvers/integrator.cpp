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
#include "../function/scaling.hpp"
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
    :  _maximum_error(e), _lipschitz_tolerance(l), _maximum_step_size(16), _function_factory_ptr(make_taylor_function_factory()), _bounder_ptr(new EulerBounder())
{
    ARIADNE_PRECONDITION(e>0.0);
    ARIADNE_PRECONDITION(l>0.0)
}

IntegratorBase::IntegratorBase(MaximumError e, SweepThreshold s, LipschitzConstant l)
    :  _maximum_error(e), _lipschitz_tolerance(l), _maximum_step_size(16), _function_factory_ptr(make_taylor_function_factory(s)), _bounder_ptr(new EulerBounder())
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

Void
IntegratorBase::set_bounder(const BounderInterface& bounder)
{
    this->_bounder_ptr=BounderPointer(bounder.clone());
}

const BounderInterface&
IntegratorBase::bounder() const
{
    return *this->_bounder_ptr;
}

Pair<StepSizeType,UpperBoxType>
IntegratorBase::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, const StepSizeType& hsug) const {
    return this->_bounder_ptr->compute(vf,D,hsug);
}

Pair<StepSizeType,UpperBoxType>
IntegratorBase::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, const ExactBoxType& A, const StepSizeType& hsug) const {
    return this->_bounder_ptr->compute(vf,D,A,hsug);
}

Pair<StepSizeType,UpperBoxType>
IntegratorBase::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, StepSizeType const& t, const ExactBoxType& A, const StepSizeType& hsug) const {
    return this->_bounder_ptr->compute(vf,D,t,A,hsug);
}


ValidatedVectorMultivariateFunctionModelDP
IntegratorBase::flow_to(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx0, const Real& tmax) const
{
    ARIADNE_LOG(1,"IntegratorBase::flow_to(ValidatedVectorMultivariateFunction vf, ExactBoxType dx0, Real tmax)\n");
    ARIADNE_LOG(2,"vf="<<vf<<"\n");
    ARIADNE_LOG(2,"dom(x0)="<<dx0<<" tmax="<<tmax<<"\n");
    const Nat n=dx0.size(); // Dimension of the state space
    ValidatedVectorMultivariateFunctionModelDP flow_function=this->function_factory().create_identity(dx0);
    StepSizeType t=0;
    ValidatedVectorMultivariateFunctionModelDP step_function;
    StepSizeType tmaxu = tmax.compute(Effort(0)).get().upper_raw();
    while(t<tmaxu) {
        ExactBoxType dx=flow_function.codomain();
        StepSizeType h_max=tmaxu-t;
        StepSizeType h;
        UpperBoxType bx;
        make_lpair(h,bx) = this->flow_bounds(vf,dx,h_max);
        Bool flow_successfully_computed=false;
        while(!flow_successfully_computed) {
            try {
                step_function=this->flow_step(vf,dx,h,bx);
                flow_successfully_computed=true;
            } catch(const FlowTimeStepException& e) {
                h=hlf(h);
            }
        }
        // FIXME: Should be able to partial evaluate FunctionModel on generic number
        // and to extract precision from a FunctionModel
        step_function=partial_evaluate(step_function,n,h);
        flow_function=compose(step_function,flow_function);
        t=t+h;
    }
    return flow_function;
}


List<ValidatedVectorMultivariateFunctionModelDP>
IntegratorBase::flow(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx0, const Real& tmin, const Real& tmax) const
{
    ARIADNE_LOG(1,"IntegratorBase::flow(ValidatedVectorMultivariateFunction vf, ExactBoxType dx0, Real tmin, Real tmax)\n");
    StepSizeType tminl = tmin.compute(Effort(0)).get().lower_raw();
    StepSizeType tmaxu = tmax.compute(Effort(0)).get().upper_raw();
    ValidatedVectorMultivariateFunctionModelDP evolve_function=this->flow_to(vf,dx0,tmin);
    StepSizeType t=tminl;
    List<ValidatedVectorMultivariateFunctionModelDP> result;

    while(possibly(t<tmax)) {
        ExactBoxType dx=evolve_function.codomain();
        StepSizeType h=tmaxu-t;
        UpperBoxType bx;
        make_lpair(h,bx) = this->flow_bounds(vf,dx,h);
        ValidatedVectorMultivariateFunctionModelDP flow_step_function=this->flow_step(vf,dx,h,bx);
        StepSizeType new_t=t+h;
        ExactIntervalType dt(t,new_t);
        ValidatedScalarMultivariateFunctionModelDP step_time_function=this->function_factory().create_identity(dt)-t;
        ValidatedVectorMultivariateFunctionModelDP flow_function=compose(flow_step_function,combine(evolve_function,step_time_function));
        ARIADNE_ASSERT(flow_function.domain()[dx0.size()].upper()==new_t);
        result.append(flow_function);
        // FIXME:
        evolve_function=partial_evaluate(flow_function,dx0.size(),new_t);
        t=new_t;
    }
    return result;
}

List<ValidatedVectorMultivariateFunctionModelDP>
IntegratorBase::flow(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx0, const Real& tmax) const
{
    return flow(vf,dx0,Real(0),tmax);
}



ValidatedVectorMultivariateFunctionModelDP
IntegratorBase::flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, StepSizeType& hmax) const
{
    ARIADNE_LOG(3,"IntegratorBase::flow_step(ValidatedVectorMultivariateFunction vf, ExactBoxType dx, StepSizeType hmax)\n");
    StepSizeType& h=hmax;
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

ValidatedVectorMultivariateFunctionModelDP
TaylorPicardIntegrator::flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const StepSizeType& h, const UpperBoxType& bx) const
{
    ARIADNE_PRECONDITION(vf.result_size()==dx.dimension());
    ARIADNE_PRECONDITION(vf.argument_size()==dx.dimension());
    ARIADNE_PRECONDITION(bx.dimension()==dx.dimension());
    return this->_flow_step(vf,dx,IntervalDomainType(0,h),BoxDomainType(0u),bx);
}

ValidatedVectorMultivariateFunctionModelDP
TaylorPicardIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& D, const Interval<StepSizeType>& T, const ExactBoxType& A, const UpperBoxType& B) const
{
    ARIADNE_PRECONDITION(f.result_size()==D.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==D.dimension()+T.dimension()+A.dimension());
    ARIADNE_PRECONDITION(B.dimension()==D.dimension());
    return this->_flow_step(f,D,IntervalDomainType(T),A,B);
}

ValidatedVectorMultivariateFunctionModelDP
TaylorPicardIntegrator::_flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& D, const ExactIntervalType& T, const ExactBoxType& A, const UpperBoxType& B) const
{
    ARIADNE_LOG(3,"TaylorPicardIntegrator::flow_step(ValidatedVectorMultivariateFunction f, ExactBoxType D, ExactIntervalType T, ExactBoxType A, UpperBoxType B)\n");
    ARIADNE_LOG(3," f="<<f);
    ARIADNE_LOG(3," D="<<D<<" T="<<T<<", A="<<A<<", B="<<B<<"\n");

    const bool is_autonomous = (f.argument_size()==D.dimension()+A.dimension());

    const Nat nx=D.size();
    const Nat na=A.size();

    Range tarng = is_autonomous ? Range(nx+1u,nx+1u+na) : Range(nx,nx+1u+na);

    StepSizeType t=static_cast<StepSizeType>(T.lower());
    StepSizeType h=static_cast<StepSizeType>(T.upper())-t;

    // Time interval centred on initial time, which will make the antiderivative more efficient
    ExactIntervalType wT(t-h,t+h);
    ARIADNE_ASSERT(t==med(wT));

    ExactBoxType dom=join(D,T,A);
    ExactBoxType wdom=join(D,wT,A);
    UpperBoxType const& bx=B;
    ARIADNE_LOG(7,"dom="<<dom<<", wdom="<<wdom<<"\n");

    ValidatedVectorMultivariateFunctionModelDP phi0=this->function_factory().create_projection(wdom,range(0,nx));
    ARIADNE_LOG(5,"phi0="<<phi0<<"\n");
    ValidatedVectorMultivariateFunctionModelDP phi=this->function_factory().create_constants(wdom,cast_singleton(bx));
    ValidatedVectorMultivariateFunctionModelDP ta=this->function_factory().create_projection(wdom,tarng);

    ARIADNE_LOG(5,"phi="<<phi<<"\n");
    for(Nat k=0; k!=this->_maximum_temporal_order; ++k) {
        Bool last_step=(phi.error().raw()<this->step_maximum_error());
        ValidatedVectorMultivariateFunctionModelDP fphi=compose(f,join(std::move(phi),ta));
        ARIADNE_LOG(5,"fphi="<<fphi<<"\n");
        // NOTE: In principle safer to use antiderivative(fphi,nx,t) here,
        // but since t is the midpoint of wdom, the (standard) antiderivative works
        // TODO: Change based antiderivative to be efficient when t is midpoint of domain
        phi=antiderivative(fphi,nx)+phi0;
        ARIADNE_LOG(4,"phi="<<phi<<"\n");
        if(last_step) { break; }
    }

    if(phi.error().raw()>this->step_maximum_error()) {
        ARIADNE_THROW(FlowTimeStepException,"TaylorPicardIntegrator::flow_step","Integration of "<<f<<" starting in "<<D<<" over time interval "<<T<<" of length "<<h<<" has error "<<phi.error()<<" after "<<this->_maximum_temporal_order<<" iterations, which exceeds step maximum error "<<this->step_maximum_error()<<"\n");
    }

    ValidatedVectorMultivariateFunctionModelDP res=restrict(phi,dom);

    //for(Nat i=0; i!=nx; ++i) { res[i]=restrict(phi[i],dom); }
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

ExactIntervalType forwards_backwards_time_domain(ExactIntervalType domt) {
    Dyadic t0(domt.lower());
    Dyadic tf(domt.upper());
    return ExactIntervalType(t0-(tf-t0),tf);
}

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
    for(SizeType i=0; i!=a1.size(); ++i) { a[i]=a1[i]; }
    a[a1.size()]=a2;
    e.append(a,c);
}


Void autonomous_flow_init(const Vector<ValidatedProcedure>& p,
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

Void autonomous_flow_iterate(const Vector<ValidatedProcedure>& p, StepSizeType h,
                  Vector<GradedValidatedDifferential>& fy, List<GradedValidatedDifferential>& t, Vector<GradedValidatedDifferential>& y)
{
    Ariadne::compute(p,fy,t,y);
    for(Nat i=0; i!=y.size(); ++i) {
        y[i]=antidifferential(fy[i]);
        y[i]*=h;
    }
}


Vector<ValidatedDifferential>
autonomous_flow_differential(Vector<GradedValidatedDifferential> const& dphia, Vector<GradedValidatedDifferential> const& dphib,
                             Vector<GradedValidatedDifferential> const& dphic, Vector<GradedValidatedDifferential> const& dphid,
                             DegreeType so, DegreeType to, Nat verbosity=0)
{
    SizeType nx=dphia.size();
    Vector<GradedValidatedDifferential> gdphi(nx,GradedValidatedDifferential(List<ValidatedDifferential>(to+1u,ValidatedDifferential(nx,so))));
    for(SizeType i=0; i!=nx; ++i) {
        for(DegreeType j=0; j!=to; ++j) {
            for(ValidatedDifferential::ConstIterator iter=dphic[i][j].begin(); iter!=dphic[i][j].end(); ++iter) {
                if(iter->index().degree()<so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
            }
            for(ValidatedDifferential::ConstIterator iter=dphid[i][j].begin(); iter!=dphid[i][j].end(); ++iter) {
                if(iter->index().degree()==so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
            }
        }
        DegreeType j=to;
        for(ValidatedDifferential::ConstIterator iter=dphia[i][j].begin(); iter!=dphia[i][j].end(); ++iter) {
            if(iter->index().degree()<so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
        }
        for(ValidatedDifferential::ConstIterator iter=dphib[i][j].begin(); iter!=dphib[i][j].end(); ++iter) {
            if(iter->index().degree()==so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
        }
    }
    ARIADNE_LOG(4,"gdphi="<<gdphi<<"\n");

    Vector<ValidatedDifferential> dphi(nx,ValidatedDifferential(nx+1u,so+to));
    for(SizeType i=0; i!=nx; ++i) {
        Expansion<MultiIndex,FloatDPBounds>& component=dphi[i].expansion();
        for(DegreeType j=0; j<=to; ++j) {
            const Expansion<MultiIndex,FloatDPBounds>& expansion=gdphi[i][j].expansion();
            for(Expansion<MultiIndex,FloatDPBounds>::ConstIterator iter=expansion.begin(); iter!=expansion.end(); ++iter) {
                append_join(component,iter->index(),j,iter->coefficient());
            }
        }
    }

    ARIADNE_LOG(4,"dphi="<<dphi<<"\n");
    return dphi;
}

ValidatedVectorMultivariateTaylorFunctionModelDP
autonomous_flow_function(const Vector<Differential<FloatBounds<DP>>>& dphi, const ExactBoxType& xh, Sweeper<FloatDP> swp) {
    const SizeType n=dphi.size();
    ValidatedVectorMultivariateTaylorFunctionModelDP tphi(n,xh,swp);

    for(SizeType i=0; i!=n; ++i) {
        ValidatedTaylorModelDP& model=tphi.model(i);
        Expansion<MultiIndex,FloatDPValue>& expansion=model.expansion();
        FloatDPError& error=model.error();
        error=0u;
        expansion.reserve(dphi[i].expansion().number_of_nonzeros());

        typename Differential<FloatDPBounds>::ConstIterator iter=dphi[i].begin();
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

ValidatedVectorMultivariateTaylorFunctionModelDP
autonomous_flow_function(const Vector<Differential<FloatBounds<DP>>>& dphi, const ExactBoxType& dx, const StepSizeType& h, Sweeper<FloatDP> swp) {
    return autonomous_flow_function(dphi,join(dx,ExactIntervalType(-h,+h)),swp);
}

ValidatedVectorMultivariateTaylorFunctionModelDP
autonomous_flow_function(const Vector<Differential<FloatBounds<DP>>>& dphi, const ExactBoxType& dx, const StepSizeType& h, double swpt) {
    ThresholdSweeper<FloatDP> swp(DP(),swpt);
    return autonomous_flow_function(dphi,dx,h,swp);
}



Void flow_init(const Vector<ValidatedProcedure>& f,
               Vector<GradedValidatedDifferential>& fy, List<GradedValidatedDifferential>& tmp,
               Vector<GradedValidatedDifferential>& yta,
               const Vector<ValidatedNumericType>& x, const ValidatedNumericType& t0, const Vector<ValidatedNumericType>& a,
               DegreeType so, DegreeType to, Nat verbosity=0)
{
    ARIADNE_LOG(2,"flow_init(f,...)\n");
    const SizeType xs=x.size();
    const SizeType as=a.size();
    const SizeType ress=f.result_size();
    const SizeType tmps=f.temporaries_size();
    const SizeType args=f.argument_size();

    assert(args==xs+1u+as or args==xs+as);
    const bool is_autonomous = (args==xs+as);

    GradedValidatedDifferential null;
    fy=Vector< GradedValidatedDifferential >(ress,null);
    tmp=List< GradedValidatedDifferential >(tmps,null);
    yta=Vector< GradedValidatedDifferential >(args,null);
    for(SizeType i=0; i!=xs; ++i) {
        yta[i]=GradedValidatedDifferential(ValidatedDifferential::variable(xs+as,so,x[i],i));
    }
    if (is_autonomous) {
        for(SizeType i=0; i!=as; ++i) {
            yta[xs+i]=GradedValidatedDifferential(ValidatedDifferential::variable(xs+as,so,a[i],xs+i));
        }
    } else {
        yta[xs]=GradedValidatedDifferential(ValidatedDifferential::constant(xs+as,so,t0));
        for(SizeType i=0; i!=as; ++i) {
            yta[xs+1u+i]=GradedValidatedDifferential(ValidatedDifferential::variable(xs+as,so,a[i],xs+i));
        }
    }
    ARIADNE_LOG(3,"fy="<<fy<<", tmp="<<tmp<<", yta="<<yta<<"\n");
}


Void flow_iterate(const Vector<ValidatedProcedure>& p,
                  Vector<GradedValidatedDifferential>& fy, List<GradedValidatedDifferential>& tmp, Vector<GradedValidatedDifferential>& yta,
                  Nat verbosity=0)
{
    ARIADNE_LOG(4,"flow_iterate(f,...)\n");
    ARIADNE_LOG(5,"degree="<<yta[0].degree()<<"\n");
    const bool is_autonomous = (p.argument_size()==yta[0][0].argument_size());
    const SizeType n=p.result_size();

    ValidatedDifferential z(yta[0][0].argument_size(),yta[0][0].degree());

    Ariadne::compute(p,fy,tmp,yta);
    for(SizeType i=0; i!=n; ++i) {
        yta[i]=antidifferential(fy[i]);
    }

    if (is_autonomous) {
        for(SizeType i=n; i!=yta.size(); ++i) { yta[i].append(z); }
    } else {
        GradedValidatedDifferential& t=yta[n];
        if(t.degree()==0 && not is_autonomous) { t.append(z+1); } else { t.append(z); }
        for(SizeType i=n+1u; i!=yta.size(); ++i) { yta[i].append(z); }
    }
}




Vector<GradedValidatedDifferential>
graded_flow_differential(Vector<GradedValidatedDifferential> const& dphic, Vector<GradedValidatedDifferential> const& dphib,
                         DegreeType so, DegreeType to, Nat verbosity=0)
{
    const SizeType nx=dphic.size();
    const SizeType m=dphic[0][0].argument_size();
    ARIADNE_PRECONDITION(dphic.size()==dphib.size());
    for(SizeType i=0; i!=dphic.size(); ++i) {
        ARIADNE_PRECONDITION(dphic[i].degree()+1u>=to);
        for(DegreeType j=0; j<to; ++j) {
            ARIADNE_ASSERT(dphic[i][j].argument_size()==m);
            ARIADNE_ASSERT(dphic[i][j].degree()>=so);
        }
        ARIADNE_PRECONDITION(dphib[i].degree()>=to);
        for(DegreeType j=0; j<=to; ++j) {
            ARIADNE_ASSERT(dphib[i][j].argument_size()==m);
            ARIADNE_ASSERT(dphib[i][j].degree()>=so);
        }
    }

    Vector<GradedValidatedDifferential> gdphi(nx,GradedValidatedDifferential(List<ValidatedDifferential>(to+1u,ValidatedDifferential(m,so))));
    for(SizeType i=0; i!=nx; ++i) {
        for(DegreeType j=0; j!=to; ++j) {
            for(ValidatedDifferential::ConstIterator iter=dphic[i][j].begin(); iter!=dphic[i][j].end(); ++iter) {
                if(iter->index().degree()<so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
            }
            for(ValidatedDifferential::ConstIterator iter=dphib[i][j].begin(); iter!=dphib[i][j].end(); ++iter) {
                if(iter->index().degree()==so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
            }
        }
        DegreeType j=to;
        for(ValidatedDifferential::ConstIterator iter=dphib[i][j].begin(); iter!=dphib[i][j].end(); ++iter) {
            gdphi[i][j].expansion().append(iter->index(),iter->coefficient());
        }
    }
    ARIADNE_LOG(4,"gdphi="<<gdphi<<"\n");

    return gdphi;
}

Vector<ValidatedDifferential>
differential(Vector<GradedValidatedDifferential> const& gdphi, SizeType gind,
             DegreeType so, DegreeType to, Nat verbosity=0)
{
    SizeType nx=gdphi.size();
    SizeType na=gdphi[0][0].argument_size();

    Vector<ValidatedDifferential> dphi(nx,na+1u,so+to);
    MultiIndex a(na+1u);
    for(SizeType i=0; i!=nx; ++i) {
        Expansion<MultiIndex,FloatDPBounds>& component=dphi[i].expansion();
        for(DegreeType j=0; j<=to; ++j) {
            a[gind]=j;
            const Expansion<MultiIndex,FloatDPBounds>& expansion=gdphi[i][j].expansion();
            for(auto term : expansion) {
                for(SizeType k=0; k!=gind; ++k) { a[k]=term.index()[k]; }
                for(SizeType k=gind; k!=na; ++k) { a[k+1u]=term.index()[k]; }
                component.append(a,term.coefficient());
            }
        }
    }
    ARIADNE_LOG(4,"dphi="<<dphi<<"\n");
    return dphi;
}

Vector<ValidatedDifferential>
flow_differential(Vector<GradedValidatedDifferential> const& dphic, Vector<GradedValidatedDifferential> const& dphib,
                  DegreeType so, DegreeType to, Nat verbosity=0)
{
    Vector<GradedValidatedDifferential> gdphi=graded_flow_differential(dphic,dphib,so,to ,verbosity);
    return differential(gdphi, dphic.size(),so,to, verbosity);
}


ValidatedVectorMultivariateTaylorFunctionModelDP make_taylor_function_model(const Vector<Differential<FloatBounds<DP>>>& df, const ExactBoxType& dom, Sweeper<FloatDP> swp, Nat verbosity=0u) {
    ARIADNE_ASSERT(df.argument_size()==dom.dimension());
    const SizeType n=df.size();
    const SizeType m=dom.dimension();
    const DegreeType deg = df.degree();
    ValidatedVectorMultivariateTaylorFunctionModelDP tf(n,dom,swp);

    Vector<Differential<FloatBounds<DP>>> ds=scale(Differential<FloatBounds<DP>>::variables(deg,Vector<FloatBounds<DP>>(m,dp)),dom);
    ARIADNE_LOG(5,"ds="<<ds<<"\n");
    Vector<Differential<FloatBounds<DP>>> dfs = compose(df,ds);

    for(SizeType i=0; i!=n; ++i) {
        ValidatedTaylorModelDP& model=tf.model(i);
        Expansion<MultiIndex,FloatDPValue>& expansion=model.expansion();
        FloatDPError& error=model.error();
        error=0u;
        expansion.reserve(dfs[i].expansion().number_of_nonzeros());

        typename Differential<FloatDPBounds>::ConstIterator iter=dfs[i].begin();
        while(iter!=dfs[i].end()) {
            MultiIndex const a=iter->index();
            FloatDPBounds coef=iter->coefficient();
            FloatDPValue x=coef.value();
            error+=coef.error();
            expansion.append(a,x);
            ++iter;
        }
        model.cleanup();
    }
    return tf;
}

ValidatedVectorMultivariateTaylorFunctionModelDP flow_function(const Vector<Differential<FloatBounds<DP>>>& dphi, const ExactBoxType& dx, const ExactIntervalType& dt, const ExactBoxType& da, Sweeper<FloatDP> swp) {
    StepSizeType t=static_cast<StepSizeType>(dt.lower());
    StepSizeType h=static_cast<StepSizeType>(dt.upper())-t;
    ExactIntervalType wdt(t-h,t+h);

    return restriction(make_taylor_function_model(dphi,join(dx,wdt,da),swp),join(dx,dt,da));
}

ValidatedVectorMultivariateTaylorFunctionModelDP flow_function(const Vector<Differential<FloatBounds<DP>>>& dphi, const ExactBoxType& dx, const ExactIntervalType& dt, const ExactBoxType& da, double swpt) {
    ThresholdSweeper<FloatDP> swp(DP(),swpt);
    return flow_function(dphi,dx,dt,da,swp);
}

} // namespace



// DEPRECATED
ValidatedVectorMultivariateFunctionModelDP
autonomous_series_flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& bdx, const StepSizeType& h, const UpperBoxType& bbx,
                            double max_err, double swpt, Nat init_so, Nat init_to, Nat max_so, Nat max_to, Nat verbosity)
{
    ARIADNE_LOG(2,"autonomous_series_flow_step(f,bdx,h,max_erro,swpt,init_so,init_to,max_so,max_to)\n");
    static const double TRY_SPACIAL_ORDER_INCREASE_FACTOR=4;

    Vector<ValidatedProcedure> p(f);
    ARIADNE_LOG(4,"p="<<p<<"\n");

    Vector<ValidatedNumericType> dx=cast_singleton(bdx);
    Vector<ValidatedNumericType> bx=cast_singleton(bbx);
    Vector<ValidatedNumericType> cx=midpoint(bdx);
    Vector<ValidatedNumericType> ax=cx+ValidatedNumericType(0,h,DoublePrecision())*evaluate(p,bx);
    ax=cx+ValidatedNumericType(0,h,DoublePrecision())*evaluate(p,ax);

    DegreeType so=init_so;
    DegreeType to=init_to;

    DegreeType nso=0;
    DegreeType nto=0;

    SizeType n=dx.size();
    Vector<ValidatedNumericType> rdx(n);
    for(SizeType i=0; i!=n; ++i) { rdx[i]=bdx[i].radius(); }

    Vector<GradedValidatedDifferential> dphia,fdphia,dphib,fdphib,dphic,fdphic,dphid,fdphid;
    List<GradedValidatedDifferential> tdphia,tdphib,tdphic,tdphid;
    Vector<GradedValidatedDifferential> ndphia,nfdphia,ndphib,nfdphib,ndphic,nfdphic,ndphid,nfdphid;
    List<GradedValidatedDifferential> ntdphia,ntdphib,ntdphic,ntdphid;

    autonomous_flow_init(p,fdphia,tdphia,dphia,ax,rdx,so);
    autonomous_flow_init(p,fdphib,tdphib,dphib,bx,rdx,so);
    autonomous_flow_init(p,fdphic,tdphic,dphic,cx,rdx,so);
    autonomous_flow_init(p,fdphid,tdphid,dphid,dx,rdx,so);

    for(DegreeType i=0; i!=to; ++i) {
        Ariadne::autonomous_flow_iterate(p,h,fdphia,tdphia,dphia);
        Ariadne::autonomous_flow_iterate(p,h,fdphib,tdphib,dphib);
        Ariadne::autonomous_flow_iterate(p,h,fdphic,tdphic,dphic);
        Ariadne::autonomous_flow_iterate(p,h,fdphid,tdphid,dphid);
    }
    ARIADNE_LOG(7,"dphia="<<dphia<<"\n");
    ARIADNE_LOG(7,"dphib="<<dphib<<"\n");
    ARIADNE_LOG(7,"dphic="<<dphic<<"\n");
    ARIADNE_LOG(7,"dphid="<<dphid<<"\n");

    Vector<ValidatedDifferential> dphi=autonomous_flow_differential(dphia,dphib,dphic,dphid,so,to);
    ARIADNE_LOG(5,"dphi="<<dphi<<"\n");

    ValidatedVectorMultivariateTaylorFunctionModelDP tphi=autonomous_flow_function(dphi,bdx,h,swpt);
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

                autonomous_flow_init(p,nfdphia,ntdphia,ndphia,ax,rdx,nso);
                autonomous_flow_init(p,nfdphib,ntdphib,ndphib,bx,rdx,nso);
                autonomous_flow_init(p,nfdphic,ntdphic,ndphic,cx,rdx,nso);
                autonomous_flow_init(p,nfdphid,ntdphid,ndphid,dx,rdx,nso);

                for(Nat i=0; i!=nto; ++i) {
                    Ariadne::autonomous_flow_iterate(p,h,nfdphia,ntdphia,ndphia);
                    Ariadne::autonomous_flow_iterate(p,h,nfdphib,ntdphib,ndphib);
                    Ariadne::autonomous_flow_iterate(p,h,nfdphic,ntdphic,ndphic);
                    Ariadne::autonomous_flow_iterate(p,h,nfdphid,ntdphid,ndphid);
                }
            }
            while(nto+1<to) {
                ++nto;
                Ariadne::autonomous_flow_iterate(p,h,nfdphia,ntdphia,ndphia);
                Ariadne::autonomous_flow_iterate(p,h,nfdphib,ntdphib,ndphib);
                Ariadne::autonomous_flow_iterate(p,h,nfdphic,ntdphic,ndphic);
                Ariadne::autonomous_flow_iterate(p,h,nfdphid,ntdphid,ndphid);
            }
            Vector<ValidatedDifferential> ndphi=autonomous_flow_differential(ndphia,ndphib,ndphic,ndphid,nso,nto);
            ValidatedVectorMultivariateTaylorFunctionModelDP ntphi=autonomous_flow_function(ndphi,bdx,h,swpt);

            SizeType nnnz=0; for(SizeType i=0; i!=tphi.size(); ++i) { nnnz+=tphi.model(i).number_of_nonzeros(); }
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
        Ariadne::autonomous_flow_iterate(p,h,fdphia,tdphia,dphia);
        Ariadne::autonomous_flow_iterate(p,h,fdphib,tdphib,dphib);
        Ariadne::autonomous_flow_iterate(p,h,fdphic,tdphic,dphic);
        Ariadne::autonomous_flow_iterate(p,h,fdphid,tdphid,dphid);
        dphi=autonomous_flow_differential(dphia,dphib,dphic,dphid,so,to);
        tphi=autonomous_flow_function(dphi,bdx,h,swpt);
    }
    SizeType nnz=0; for(SizeType i=0; i!=tphi.size(); ++i) { nnz+=tphi.model(i).number_of_nonzeros(); }
    ARIADNE_LOG(2,"so="<<so<<" to="<<to<<" nnz="<<nnz<<" err="<<tphi.error()<<"\n");
    ARIADNE_LOG(4,"phi="<<tphi<<"\n");
    return tphi;
}


// Flow step using graded differential with fixed degree
ValidatedVectorMultivariateTaylorFunctionModelDP
graded_series_flow_step(const Vector<ValidatedProcedure>& f,
                        const ExactBoxType& bdx, const ExactIntervalType& idt, const ExactBoxType& bda, const UpperBoxType& bbx,
                        double swpt, DegreeType so, DegreeType to, Nat verbosity=0u)
{
    ARIADNE_LOG(2,"graded_series_flow_step(f,bdx,idt,bda,bbx,swpt,so,to)\n");
    ARIADNE_LOG(3,"f="<<f<<"\n");
    ARIADNE_LOG(3,"bdx="<<bdx<<", idt="<<idt<<", bda="<<bda<<", bbx="<<bbx<<"\n");
    ARIADNE_LOG(3,"swpt="<<swpt<<", so="<<so<<", to="<<to<<"\n");

    ARIADNE_PRECONDITION(f.result_size()==bdx.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==bdx.dimension()+bda.dimension() || f.argument_size()==bdx.dimension()+1u+bda.dimension());

    SizeType n=bdx.dimension();
    SizeType m=bdx.dimension()+idt.dimension()+bda.dimension();

    StepSizeType t=static_cast<StepSizeType>(idt.lower());
    StepSizeType h=static_cast<StepSizeType>(idt.upper())-t;

    ExactIntervalType widt(t-h,t+h);

    Vector<ValidatedNumericType> dx=cast_singleton(bdx);
    Scalar<ValidatedNumericType> dt=cast_singleton(idt);
    Scalar<ValidatedNumericType> wdt=cast_singleton(widt);
    Vector<ValidatedNumericType> da=cast_singleton(bda);
    Vector<ValidatedNumericType> bx=cast_singleton(bbx);

    Vector<ValidatedNumericType> mdx=midpoint(dx);
    Scalar<ValidatedNumericType> mdt=midpoint(widt);
    Vector<ValidatedNumericType> mda=midpoint(da);

    Vector<ValidatedNumericType> domc=midpoint(join(dx,wdt,da));
    Vector<ValidatedNumericType> domb=join(bx,dt,da);
    ExactBoxType bdomc=join(bdx,idt,bda);

    Vector<ValidatedNumericType> rdx(m);
    for(SizeType i=0; i!=m; ++i) { rdx[i]=bdomc[i].radius(); }

    ARIADNE_LOG(4,"dx="<<dx<<", dt="<<dt<<", da="<<da<<", wdt="<<wdt<<", bx="<<bx<<"\n");

    Vector<GradedValidatedDifferential> dphic,fdphic,dphib,fdphib;
    List<GradedValidatedDifferential> tdphic,tdphib;

    Ariadne::flow_init(f,fdphic,tdphic,dphic,mdx,mdt,mda,so,to, verbosity);
    Ariadne::flow_init(f,fdphib,tdphib,dphib,bx,dt,da,so,to, verbosity);

    for(DegreeType i=0; i!=to; ++i) {
        Ariadne::flow_iterate(f,fdphic,tdphic,dphic, verbosity);
        Ariadne::flow_iterate(f,fdphib,tdphib,dphib, verbosity);
    }
    ARIADNE_LOG(7,"dphic="<<dphic<<"\n");
    ARIADNE_LOG(7,"dphib="<<dphib<<"\n");

    dphic=project(dphic,range(0,n));
    dphib=project(dphib,range(0,n));
    ARIADNE_LOG(7,"dphic="<<dphic<<"\n");
    ARIADNE_LOG(7,"dphib="<<dphib<<"\n");

    Vector<ValidatedDifferential> dphi=Ariadne::flow_differential(dphic,dphib,so,to,verbosity);
    ARIADNE_LOG(5,"dphi="<<dphi<<"\n");

    ValidatedVectorMultivariateTaylorFunctionModelDP tphi=Ariadne::flow_function(dphi,bdx,idt,bda,swpt);
    ARIADNE_LOG(5,"phi="<<tphi<<"\n");

    return tphi;
}

// Flow step using graded differential with varying degree and specified maximum error
ValidatedVectorMultivariateFunctionModelDP
graded_series_flow_step(const Vector<ValidatedProcedure>& f,
                        const ExactBoxType& domx, const ExactIntervalType& domt, const ExactBoxType& doma, const UpperBoxType& bndx,
                        double max_err, double swpt, DegreeType init_so, DegreeType init_to, DegreeType max_so, DegreeType max_to, Nat verbosity)
{
    ARIADNE_LOG(2,"graded_series_flow_step(f,domx,domt,doma,bndx,max_err,swpt,init_so,init_to,max_so,max_to)\n");
    ARIADNE_LOG(3,"f="<<f<<"\n");
    ARIADNE_LOG(3,"domx="<<domx<<", domt="<<domt<<", doma="<<doma<<", bndx="<<bndx<<"\n");
    ARIADNE_LOG(3,"max_err="<<max_err<<", swpt="<<swpt<<", "<<
                  "init_so="<<init_so<<", init_to="<<init_to<<", max_so="<<max_so<<", max_to="<<max_to<<"\n");

    ARIADNE_PRECONDITION(f.result_size()==domx.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==domx.dimension()+doma.dimension() || f.argument_size()==domx.dimension()+1u+doma.dimension());

    static const double TRY_SPACIAL_ORDER_INCREASE_FACTOR=4;

    ARIADNE_LOG(4,"f="<<f<<"\n");

    DegreeType so=init_so;
    DegreeType to=init_to;

    ValidatedVectorMultivariateTaylorFunctionModelDP phi=graded_series_flow_step(f,domx,domt,doma,bndx, swpt,so,to, verbosity);

    ARIADNE_LOG(5,"phi="<<phi<<"\n");
    SizeType nnz=0; for(SizeType i=0; i!=phi.size(); ++i) { nnz+=phi.model(i).number_of_nonzeros(); }
    ARIADNE_LOG(3,"so="<<so<<" to="<<to<<" nnz="<<nnz<<" err="<<phi.error()<<"\n");

    FloatDPError old_error=phi.error()*FloatDPError(TRY_SPACIAL_ORDER_INCREASE_FACTOR*2);

    while(phi.error().raw()>max_err && (so<max_so || to<max_to) ) {

        old_error = phi.error();

        if( (so<max_so) && ((phi.error()*FloatDPError(TRY_SPACIAL_ORDER_INCREASE_FACTOR)).raw() > old_error.raw()) ) {
            // try increasing spacial degree
            ++so;
            to=init_to;
        } else {
            ++to;
        }

        phi=graded_series_flow_step(f,domx,domt,doma,bndx, swpt,so,to, verbosity);

        SizeType nnnz=0; for(SizeType i=0; i!=phi.size(); ++i) { nnnz+=phi.model(i).number_of_nonzeros(); }
        ARIADNE_LOG(3,"so="<<so<<" to="<<to<<" nnz="<<nnz<<" err="<<phi.error()<<"\n");

    }
    ARIADNE_LOG(4,"phi="<<phi<<"\n");
    return phi;
}


ValidatedVectorMultivariateFunctionModelDP
graded_series_flow_step(const ValidatedVectorMultivariateFunction& f,
                        const ExactBoxType& domx, const Interval<StepSizeType>& domt, const ExactBoxType& doma, const UpperBoxType& bndx,
                        double max_err, double swpt, DegreeType init_so, DegreeType init_to, DegreeType max_so, DegreeType max_to, Nat verbosity)
{
    Vector<ValidatedProcedure> p(f);
    ExactIntervalType idt(domt);

    return graded_series_flow_step(p,domx,idt,doma,bndx,max_err,swpt,init_so,init_to, max_so, max_to, verbosity);
}


// FIXME: Should not be necessary, as should be able to construct FloatBounds<DP> from (FloatValue<DP>,DP)
FloatBounds<DoublePrecision> cast_singleton(ExactIntervalType const& ivl, DoublePrecision pr) {
    return FloatBounds<DoublePrecision>(ivl.lower(),ivl.upper()); }


namespace {
// Compute the midpoint of x, and add the error to e
template<class F, class FE> Value<F> med(Bounds<F> const& x, Error<FE>& e) {
    e+=x.error(); return x.value(); }
} // namespace


template<class FLT> ValidatedVectorMultivariateTaylorFunctionModel<FLT>
make_taylor_function_model(ExactBoxType domain, Vector<Differential<Bounds<FLT>>> centre_derivatives, Vector<Differential<Bounds<FLT>>> derivative_ranges, Sweeper<FLT> swp) {
    ARIADNE_ASSERT(centre_derivatives.result_size()==derivative_ranges.result_size());
    ARIADNE_ASSERT(centre_derivatives.argument_size()==domain.dimension());
    ARIADNE_ASSERT(derivative_ranges.argument_size()==domain.dimension());
    ARIADNE_ASSERT(centre_derivatives.degree()==derivative_ranges.degree() or centre_derivatives.degree()+1u==derivative_ranges.degree());

    using PR = PrecisionType<FLT>;
    using X = FloatBounds<PR>;

    const SizeType m=derivative_ranges.result_size();
    const SizeType n=derivative_ranges.argument_size();
    const DegreeType deg = derivative_ranges.degree();
    FloatBounds<PR> z=derivative_ranges.zero_element().zero_coefficient();

    auto scalings = Vector<Differential<FloatBounds<PR>>>(n,[&](SizeType i){return Differential<X>::variable(n,deg,z,i)*rad(domain[i]);});
    auto scaled_centre_derivatives = compose(centre_derivatives,scalings);
    auto scaled_derivative_ranges = compose(derivative_ranges,scalings);


    // Make the models
    ValidatedVectorMultivariateTaylorFunctionModel<FLT> tf(m,domain,swp);

    for(SizeType i=0; i!=m; ++i) {
        Differential<FloatBounds<PR>> const& dc = scaled_centre_derivatives[i];
        Differential<FloatBounds<PR>> const& dr = scaled_derivative_ranges[i];
        ValidatedTaylorModel<FLT>& model=tf.model(i);

        Expansion<MultiIndex,FloatDPValue>& expansion=model.expansion();
        FloatError<PR>& error=model.error();
        error=0u;
        expansion.reserve(centre_derivatives[i].expansion().number_of_nonzeros());

        auto riter=dr.begin();
        FloatBounds<PR> coef;

        // Since coefficients are stored in increasing total degree, can first do centre and then ranges
        for(auto centre_iter=dc.begin(); centre_iter!=dc.end() && centre_iter->index().degree()<deg; ++centre_iter) {
            expansion.append(centre_iter->index(),med(centre_iter->coefficient(),error));
        }

        if(not (dr.expansion().empty() or dr.expansion().back().index().degree()<deg)) {
            auto range_iter=dr.begin(); while (range_iter->index().degree()<deg) { ++range_iter; }
            for ( ; range_iter!=dr.end(); ++range_iter) {
                expansion.append(range_iter->index(),med(range_iter->coefficient(),error));
            }
        }

        model.cleanup();
    }
    return tf;
}


// Solve \f$\dt{\phi}(x,t,a)=f(\phi(x,t),t,a)\f$ for x in dx, t in dt, and a in da, assuming x remains in bx.
ValidatedVectorMultivariateFunctionModelDP
series_flow_step(const ValidatedVectorMultivariateFunction& f,
                 const ExactBoxType& domx,
                 const ExactIntervalType& domt,
                 const ExactBoxType& doma,
                 const UpperBoxType& bndbx,
                 DegreeType deg,
                 Sweeper<FloatDP> swp,
                 Nat verbosity=0)
{
    ARIADNE_PRECONDITION(f.result_size()==domx.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==domx.dimension()+doma.dimension()
                            || f.argument_size()==domx.dimension()+domt.dimension()+doma.dimension());
    const bool is_autonomous = f.argument_size()==domx.dimension()+doma.dimension();

    typedef DoublePrecision PR;
    typedef FloatBounds<PR> X;
    PR pr;

    // Extend time domain from [t:t+h] to [t-h:t+h]
    auto wide_domt = forwards_backwards_time_domain(domt);
    ExactBoxType domxta=product(domx,wide_domt,doma);

    Vector<X> cx(midpoint(domx),pr);
    X t0(domt.lower(),pr);
    Vector<X> ca(midpoint(doma),pr);
    Vector<Differential<X>> cdf = is_autonomous ? f.differential(join(cx,ca),deg) : f.differential(join(cx,t0,ca),deg);
    Vector<Differential<X>> centre_flow_derivatives = is_autonomous ? flow(cdf, cx,ca) : flow(cdf, cx,t0,ca);

    Vector<X> bndx=cast_singleton(bndbx);
    Vector<X> rngx=cast_singleton(domx,pr);
    X rngt=cast_singleton(domt,pr);
    Vector<X> rnga=cast_singleton(doma,pr);
    Vector<Differential<X>> rngdf = is_autonomous ? f.differential(join(bndx,rnga),deg) : f.differential(join(bndx,rngt,rnga),deg);
    Vector<Differential<X>> range_flow_derivatives = is_autonomous ? flow(rngdf, bndx,rnga) : flow(rngdf, bndx,rngt,rnga);

    ValidatedVectorMultivariateFunctionModelDP forwards_backwards_taylor_function_model = make_taylor_function_model(domxta, centre_flow_derivatives, range_flow_derivatives, swp);
    domxta[domx.size()]=ExactIntervalType(domt);
    ValidatedVectorMultivariateFunctionModelDP forwards_taylor_function_model=restriction(forwards_backwards_taylor_function_model,domxta);
    return forwards_taylor_function_model;
}



ValidatedVectorMultivariateFunctionModelDP
TaylorSeriesIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& domx, const StepSizeType& h, const UpperBoxType& bndx) const
{
    ARIADNE_LOG(3,"TaylorSeriesIntegrator::flow_step(ValidatedVectorMultivariateFunction f, ExactBoxType domx, StepSizeType h, const UpperBoxType& bndx)\n");

    Interval<StepSizeType> domt(0,h);
    ExactBoxType doma;

    return this->flow_step(f, domx,domt,doma, bndx);
}

ValidatedVectorMultivariateFunctionModelDP
TaylorSeriesIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& domx, const Interval<StepSizeType>& rngt, const ExactBoxType& doma, const UpperBoxType& bndx) const
{
    ExactIntervalType domt(rngt);
    double max_err=this->step_maximum_error();
    double swpt=this->step_sweep_threshold();

    DegreeType init_so=this->minimum_spacial_order();
    DegreeType init_to=this->minimum_temporal_order();
    DegreeType max_so=this->maximum_spacial_order();
    DegreeType max_to=this->maximum_temporal_order();

    Vector<ValidatedProcedure> p(f);
    ValidatedVectorMultivariateFunctionModelDP tphi=Ariadne::graded_series_flow_step(p,domx,domt,doma,bndx,
        max_err,swpt, init_so,init_to,max_so,max_to, this->verbosity);

/*
    ThresholdSweeper<FloatDP> swp(DP(),this->step_sweep_threshold());
    DegreeType deg = this->maximum_temporal_order();
    ValidatedVectorMultivariateFunctionModelDP tphi=Ariadne::series_flow_step(f,domx,domt,doma,bndx,deg,swp,this->verbosity);
*/

    if(tphi.error().raw()>this->step_maximum_error()) {
        ARIADNE_THROW(FlowTimeStepException,"TaylorSeriesIntegrator::flow_step",
                      "Integration of "<<f<<" over "<<domx<<" for time interval "<<domt<<" has error "<<tphi.errors()<<
                      " using spacial order "<<max_so<<" and temporal order "<<max_to<<
                      ", which exceeds maximum single-step error "<<max_err<<"\n");
    }

    return tphi;
}

Pair<StepSizeType,UpperBoxType>
TaylorSeriesIntegrator::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const StepSizeType& hmax) const
{
    return this->bounder().compute(vf,dx,hmax);
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
        UniformConstReference<MultiIndex> index = read_iter->index();
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
AffineIntegrator::flow_derivative(const ValidatedVectorMultivariateFunction& f, const Vector<ValidatedNumericType>& dom) const
{
    ARIADNE_LOG(5,"AffineIntegrator::flow_derivative(ValidatedVectorMultivariateFunction f, ValidatedBoxType dom)\n");
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

ValidatedVectorMultivariateFunctionModelDP
AffineIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& dom, const StepSizeType& h, const UpperBoxType& bbox) const
{
    ARIADNE_LOG(3,"AffineIntegrator::flow_step(ValidatedVectorMultivariateFunction f, ExactBoxType dom, StepSizeType h, UpperBoxType bbox)\n");
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
    rad[n] = abs(h);

    for(Nat i=0; i!=n; ++i) {
        for(Expansion<MultiIndex,ValidatedNumericType>::ConstIterator iter=bdphi[i].begin(); iter!=bdphi[i].end(); ++iter) {
            UniformConstReference<MultiIndex> a=iter->index();
            if(a[n]==this->_temporal_order && a[n]+this->_spacial_order==a.degree()) {
                UniformConstReference<ValidatedNumericType> rng = iter->coefficient();
                UniformConstReference<ValidatedNumericType> mid = mdphi[i][a];
                ARIADNE_ASSERT(rng.lower().raw()<=mid.lower().raw() && mid.upper().raw()<=rng.upper().raw());
                FloatDPError mag = FloatDPError(max(rng.upper()-mid.lower(),mid.upper()-rng.lower()));
                for(Nat j=0; j!=n+1; ++j) { mag *= pow(rad[j],Nat(a[j])); }
                err[i] += mag;
            }
        }
    }

    ExactBoxType flow_domain = join(dom,ExactIntervalType(0,h));

    ValidatedVectorMultivariateFunctionModelDP id = this->function_factory().create_identity(flow_domain);
    ValidatedVectorMultivariateFunctionModelDP res = this->function_factory().create_zeros(n,flow_domain);
    for(Nat i=0; i!=n; ++i) {
        ValidatedScalarMultivariateFunctionModelDP res_model = res[i] + mdphi[i].expansion()[MultiIndex::zero(n+1)];
        for(Nat j=0; j!=mdphi[i].argument_size()-1; ++j) { res_model+=mdphi[i].expansion()[MultiIndex::unit(n+1,j)]*(id[j]-ValidatedNumericType(midpoint(flow_domain[j]))); }
        Nat j=mdphi[i].argument_size()-1; { res_model+=mdphi[i].expansion()[MultiIndex::unit(n+1,j)]*id[j]; }
        res_model += FloatDPBounds(-err[i],+err[i]);
        res[i]=res_model;
    }
    return res;
}

ValidatedVectorMultivariateFunctionModelDP
AffineIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& D, const Interval<StepSizeType>& T, const ExactBoxType& A, const UpperBoxType& B) const
{
    ARIADNE_NOT_IMPLEMENTED;
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
