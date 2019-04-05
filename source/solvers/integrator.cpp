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
    const SizeType n=dx0.size(); // Dimension of the state space
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

    const SizeType nx=D.size();
    const SizeType na=A.size();

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
    for(DegreeType k=0; k!=this->_maximum_temporal_order; ++k) {
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

    //for(SizeType i=0; i!=nx; ++i) { res[i]=restrict(phi[i],dom); }
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
typedef MultivariatePolynomial<FloatDPBounds> ValidatedMultivariatelynomial;
typedef Graded<ValidatedDifferential> GradedValidatedDifferential;
Bool operator<(const MultiIndex& a1, const MultiIndex& a2);


TaylorSeriesIntegrator::TaylorSeriesIntegrator(
        MaximumError err, SweepThreshold swp, LipschitzConstant lip,
        StepSweepThreshold stpswp, Order ord)
    : IntegratorBase(err,swp,lip), _step_sweep_threshold(stpswp), _order(ord)
{ }

TaylorSeriesIntegrator::TaylorSeriesIntegrator(
        MaximumError err, Order ord)
    : TaylorSeriesIntegrator(err,SweepThreshold(err/1024),LipschitzConstant(0.5),StepSweepThreshold(err/1048576),ord)
{ }


GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator(MaximumError err)
    : GradedTaylorSeriesIntegrator(err,SweepThreshold(err/1024),LipschitzConstant(0.5))
{ }

GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator(MaximumError err, SweepThreshold swp, LipschitzConstant lip)
    : GradedTaylorSeriesIntegrator(err,swp,lip,StepMaximumError(err/128),StepSweepThreshold(swp/1024),MaximumTemporalOrder(12))
{ }

GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator(MaximumError err, SweepThreshold swp, LipschitzConstant lip,
                        StepMaximumError stperr, StepSweepThreshold stpswp, MaximumTemporalOrder maxto)
    : GradedTaylorSeriesIntegrator(err,swp,lip,stperr,stpswp,MinimumSpacialOrder(1),MinimumTemporalOrder(4),MaximumSpacialOrder(4),maxto)
{
}

GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator(
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

template<class F> GradedValidatedDifferential flow(const F& f, const ExactIntervalType& c, DegreeType so, DegreeType to) {
    ValidatedDifferential x=make_differential_variable(1u,so,cast_singleton(c),0u);
    GradedValidatedDifferential y=make_graded(x);
    GradedValidatedDifferential t=create_graded(x);

    for(DegreeType n=0; n!=to; ++n) {
        t=f(y);
        y=antidifferential(t);
    }

    return y;
}

template<class X> Void append_join(Expansion<MultiIndex,X>& e, const MultiIndex& a1, const DegreeType a2, const X& c) {
    MultiIndex a(a1.size()+1);
    for(SizeType i=0; i!=a1.size(); ++i) { a[i]=a1[i]; }
    a[a1.size()]=a2;
    e.append(a,c);
}



Void graded_flow_init(const Vector<ValidatedProcedure>& f,
               Vector<GradedValidatedDifferential>& fy, List<GradedValidatedDifferential>& tmp,
               Vector<GradedValidatedDifferential>& yta,
               const Vector<ValidatedNumericType>& x, const ValidatedNumericType& t0, const Vector<ValidatedNumericType>& a,
               DegreeType so, DegreeType to, Nat verbosity=0)
{
    ARIADNE_LOG(2,"graded_flow_init(f,...)\n");
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


Void graded_flow_iterate(const Vector<ValidatedProcedure>& p,
                         Vector<GradedValidatedDifferential>& fy, List<GradedValidatedDifferential>& tmp, Vector<GradedValidatedDifferential>& yta,
                         Nat verbosity=0)
{
    ARIADNE_LOG(4,"graded_flow_iterate(f,...)\n");
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
    const SizeType rs=dphic.size();
    const SizeType as=dphic[0][0].argument_size();

    Vector<GradedValidatedDifferential> gdphi(rs,GradedValidatedDifferential(List<ValidatedDifferential>(to+1u,ValidatedDifferential(as,so))));
    for(SizeType i=0; i!=rs; ++i) {
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
    SizeType rs=gdphi.size();
    SizeType as=gdphi[0][0].argument_size();

    Vector<ValidatedDifferential> dphi(rs,as+1u,so+to);
    MultiIndex a(as+1u);
    for(SizeType i=0; i!=rs; ++i) {
        Expansion<MultiIndex,FloatDPBounds>& component=dphi[i].expansion();
        for(DegreeType j=0; j<=to; ++j) {
            a[gind]=j;
            const Expansion<MultiIndex,FloatDPBounds>& expansion=gdphi[i][j].expansion();
            for(auto term : expansion) {
                for(SizeType k=0; k!=gind; ++k) { a[k]=term.index()[k]; }
                for(SizeType k=gind; k!=as; ++k) { a[k+1u]=term.index()[k]; }
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
    const SizeType rs=df.size();
    const SizeType as=dom.dimension();
    const DegreeType deg = df.degree();
    ValidatedVectorMultivariateTaylorFunctionModelDP tf(rs,dom,swp);

    Vector<Differential<FloatBounds<DP>>> ds=scale(Differential<FloatBounds<DP>>::variables(deg,Vector<FloatBounds<DP>>(as,dp)),dom);
    ARIADNE_LOG(5,"ds="<<ds<<"\rs");
    Vector<Differential<FloatBounds<DP>>> dfs = compose(df,ds);

    for(SizeType i=0; i!=rs; ++i) {
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

ValidatedVectorMultivariateTaylorFunctionModelDP flow_function(const Vector<Differential<FloatBounds<DP>>>& dphi, const ExactBoxType& domx, const ExactIntervalType& domt, const ExactBoxType& doma, Sweeper<FloatDP> swp) {
    StepSizeType t=static_cast<StepSizeType>(domt.lower());
    StepSizeType h=static_cast<StepSizeType>(domt.upper())-t;
    ExactIntervalType wdt(t-h,t+h);

    return restriction(make_taylor_function_model(dphi,join(domx,wdt,doma),swp),join(domx,domt,doma));
}

ValidatedVectorMultivariateTaylorFunctionModelDP flow_function(const Vector<Differential<FloatBounds<DP>>>& dphi, const ExactBoxType& domx, const ExactIntervalType& domt, const ExactBoxType& doma, double swpt) {
    ThresholdSweeper<FloatDP> swp(DP(),swpt);
    return flow_function(dphi,domx,domt,doma,swp);
}

} // namespace


// Flow step using graded differential with fixed degree
ValidatedVectorMultivariateTaylorFunctionModelDP
graded_series_flow_step(const Vector<ValidatedProcedure>& f,
                        const ExactBoxType& domx, const ExactIntervalType& domt, const ExactBoxType& doma, const UpperBoxType& bndx,
                        double swpt, DegreeType so, DegreeType to, Nat verbosity=0u)
{
    ARIADNE_LOG(2,"graded_series_flow_step(f,domx,domt,doma,bndx,swpt,so,to)\n");
    ARIADNE_LOG(3,"f="<<f<<"\n");
    ARIADNE_LOG(3,"domx="<<domx<<", domt="<<domt<<", doma="<<doma<<", bndx="<<bndx<<"\n");
    ARIADNE_LOG(3,"swpt="<<swpt<<", so="<<so<<", to="<<to<<"\n");

    ARIADNE_PRECONDITION(f.result_size()==domx.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==domx.dimension()+doma.dimension() || f.argument_size()==domx.dimension()+1u+doma.dimension());

    SizeType nx=domx.dimension();

    StepSizeType t=static_cast<StepSizeType>(domt.lower());
    StepSizeType h=static_cast<StepSizeType>(domt.upper())-t;

    ExactIntervalType widt(t-h,t+h);

    Vector<ValidatedNumericType> dx=cast_singleton(domx);
    Scalar<ValidatedNumericType> dt=cast_singleton(domt);
    Scalar<ValidatedNumericType> wdt=cast_singleton(widt);
    Vector<ValidatedNumericType> da=cast_singleton(doma);
    Vector<ValidatedNumericType> bx=cast_singleton(bndx);

    Vector<ValidatedNumericType> mdx=midpoint(dx);
    Scalar<ValidatedNumericType> mdt=midpoint(widt);
    Vector<ValidatedNumericType> mda=midpoint(da);

    Vector<ValidatedNumericType> dc=midpoint(join(dx,wdt,da));
    Vector<ValidatedNumericType> db=join(bx,dt,da);
    ExactBoxType domc=join(domx,domt,doma);

    ARIADNE_LOG(4,"dx="<<dx<<", dt="<<dt<<", da="<<da<<", wdt="<<wdt<<", bx="<<bx<<"\n");

    Vector<GradedValidatedDifferential> dphic,fdphic,dphib,fdphib;
    List<GradedValidatedDifferential> tmpdphic,tmpdphib;

    Ariadne::graded_flow_init(f,fdphic,tmpdphic,dphic,mdx,mdt,mda,so,to, verbosity);
    Ariadne::graded_flow_init(f,fdphib,tmpdphib,dphib,bx,dt,da,so,to, verbosity);

    for(DegreeType i=0; i!=to; ++i) {
        Ariadne::graded_flow_iterate(f,fdphic,tmpdphic,dphic, verbosity);
        Ariadne::graded_flow_iterate(f,fdphib,tmpdphib,dphib, verbosity);
    }
    ARIADNE_LOG(7,"dphic="<<dphic<<"\n");
    ARIADNE_LOG(7,"dphib="<<dphib<<"\n");

    dphic=project(dphic,range(0,nx));
    dphib=project(dphib,range(0,nx));
    ARIADNE_LOG(7,"dphic="<<dphic<<"\n");
    ARIADNE_LOG(7,"dphib="<<dphib<<"\n");

    Vector<ValidatedDifferential> dphi=Ariadne::flow_differential(dphic,dphib,so,to,verbosity);
    ARIADNE_LOG(5,"dphi="<<dphi<<"\n");

    ValidatedVectorMultivariateTaylorFunctionModelDP tphi=Ariadne::flow_function(dphi,domx,domt,doma,swpt);
    ARIADNE_LOG(5,"phi="<<tphi<<"\n");

    return tphi;
}

// Flow step using graded differential with varying degree and specified maximum error
ValidatedVectorMultivariateFunctionModelDP
graded_series_flow_step(const Vector<ValidatedProcedure>& f,
                        const ExactBoxType& domx, const ExactIntervalType& domt, const ExactBoxType& doma, const UpperBoxType& bndx,
                        double max_err, double swpt, DegreeType init_so, DegreeType init_to, DegreeType max_so, DegreeType max_to, Nat verbosity=0)
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
    ExactIntervalType idomt(domt);

    return graded_series_flow_step(p,domx,idomt,doma,bndx,max_err,swpt,init_so,init_to, max_so, max_to, verbosity);
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


// Solve \f$\dt{\phi}(x,t,a)=f(\phi(x,t),t,a)\f$ for x in domx, t in domt, and a in doma, assuming x remains in bndx.
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
    ThresholdSweeper<FloatDP> swp(DP(),this->step_sweep_threshold());
    DegreeType deg = this->order();
    ValidatedVectorMultivariateFunctionModelDP tphi=Ariadne::series_flow_step(f,domx,domt,doma,bndx,deg,swp,this->verbosity);
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
       << ", step_sweep_threshold = " << this->step_sweep_threshold()
       << ", order = " << this->order()
       << " )";
}



ValidatedVectorMultivariateFunctionModelDP
GradedTaylorSeriesIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& domx, const StepSizeType& h, const UpperBoxType& bndx) const
{
    ARIADNE_LOG(3,"GradedTaylorSeriesIntegrator::flow_step(ValidatedVectorMultivariateFunction f, ExactBoxType domx, StepSizeType h, const UpperBoxType& bndx)\n");

    Interval<StepSizeType> domt(0,h);
    ExactBoxType doma;

    return this->flow_step(f, domx,domt,doma, bndx);
}

ValidatedVectorMultivariateFunctionModelDP
GradedTaylorSeriesIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& domx, const Interval<StepSizeType>& rngt, const ExactBoxType& doma, const UpperBoxType& bndx) const
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
        ARIADNE_THROW(FlowTimeStepException,"GradedTaylorSeriesIntegrator::flow_step",
                      "Integration of "<<f<<" over "<<domx<<" for time interval "<<domt<<" has error "<<tphi.errors()<<
                      " using spacial order "<<max_so<<" and temporal order "<<max_to<<
                      ", which exceeds maximum single-step error "<<max_err<<"\n");
    }

    return tphi;
}

Pair<StepSizeType,UpperBoxType>
GradedTaylorSeriesIntegrator::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const StepSizeType& hmax) const
{
    return this->bounder().compute(vf,dx,hmax);
}

Void GradedTaylorSeriesIntegrator::write(OutputStream& os) const {
    os << "GradedTaylorSeriesIntegrator"
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



template<class X> Void truncate(Differential<X>& x, DegreeType spacial_order_, DegreeType temporal_order_) {
    SizeType n=x.argument_size()-1;
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

template<class X> Void truncate(Vector< Differential<X> >& x, DegreeType spacial_order_, DegreeType temporal_order_) {
    for(DegreeType i=0; i!=x.size(); ++i) { truncate(x[i],spacial_order_,temporal_order_); }
}

template<class X> inline Vector<Differential<X>> flow_differential(Vector<Differential<X>> const& df, Vector<X> const& x0) {
    return flow(df,x0); }


TaylorPolynomialIntegrator::TaylorPolynomialIntegrator(MaximumError e, LipschitzConstant l, Order o)
    : IntegratorBase(e,l)
    , _order(o)
{
}

ValidatedVectorMultivariateFunctionModelDP
TaylorPolynomialIntegrator::flow_step(const ValidatedVectorMultivariateFunction& vector_field,
                                      const ExactBoxType& state_domain,
                                      StepSizeType& suggested_time_step) const
{
    std::cerr<<"vector_field="<<"<rs="<<vector_field.result_size()<<",as="<<vector_field.argument_size()<<">"<< vector_field<<"\n";

    StepSizeType const& time_step = suggested_time_step;
    IntervalDomainType time_domain(0,time_step);
    IntervalDomainType symmetric_time_domain(-time_step,+time_step);

    SizeType time_index = state_domain.dimension();

    Sweeper<FloatDP> sweeper(ThresholdSweeper<FloatDP>(double_precision, 1e-10));

    Vector<FloatDPBounds> x0 = state_domain.centre();
    FloatDPBounds t0 = symmetric_time_domain.centre();
    Vector<FloatDPBounds> x=cast_singleton(state_domain);
    Vector<FloatDPBounds> wx=2*x-x0;
    Vector<FloatDPBounds> xt0 = join(x0,t0);
    Vector<FloatDPBounds> xt = cast_singleton(join(state_domain,time_domain));
    DegreeType deg=this->_order;

    deg=5;

    Vector<Differential<FloatDPBounds>> df0 = vector_field.differential(x0,deg);
    Vector<Differential<FloatDPBounds>> df  = vector_field.differential(x,deg);

    FloatDPValue widen(2.0_dy,double_precision);
    Vector<Differential<FloatDPBounds>> centre_flow_differential = flow_differential(df0, x0);
    Vector<Differential<FloatDPBounds>> bounds_flow_differential = flow_differential(df, wx);

    Vector<Differential<FloatDPBounds>> widened_bounds_flow_differential = (bounds_flow_differential-centre_flow_differential)*widen+centre_flow_differential;

    std::cerr<<"centre_flow_differential="<<centre_flow_differential<<"\n";
    ExactBoxType flow_domain = product(state_domain,time_domain);
    ExactBoxType symmetric_flow_domain = product(state_domain,symmetric_time_domain);
    ValidatedVectorMultivariateTaylorFunctionModelDP tfm
        = make_taylor_function_model(symmetric_flow_domain,centre_flow_differential,widened_bounds_flow_differential,sweeper);
    std::cerr<<"flow_model="<<tfm<<"\n";

    tfm.restrict(flow_domain);

    Vector<FloatDPBounds> flow_range_estimate = cast_singleton(tfm.range());
    widened_bounds_flow_differential = flow_differential(df, flow_range_estimate);
    tfm = make_taylor_function_model(symmetric_flow_domain,centre_flow_differential,widened_bounds_flow_differential,sweeper);


    const SizeType nx=time_index;
#warning Simplify construction
    ValidatedVectorMultivariateTaylorFunctionModelDP phi0(nx,flow_domain,sweeper);
    for(SizeType i=0; i!=nx; ++i) { phi0[i] = ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(flow_domain,i,sweeper); }

    //ValidatedVectorMultivariateFunctionModelDP phi0=factory(tfm).create_projection(range(0,time_index));

    //std::cerr<<"compose(vector_field,tfm)="<<compose(tfflow_model="<<tfm<<"\n";

    ValidatedVectorMultivariateTaylorFunctionModelDP pic_tfm = phi0 + antiderivative( compose(vector_field, tfm) , time_index, 0);

    tfm=pic_tfm; pic_tfm = phi0 + antiderivative( compose(vector_field, tfm) , time_index, 0);

    auto diff_tfm = pic_tfm-tfm; diff_tfm.clobber();
    std::cerr<<"\ndiff="<<diff_tfm << "\n\nnorm="<<norm(diff_tfm) << ", e1="<<pic_tfm.error()<<", e0="<<tfm.error()<<"\n";
    std::cerr<<"refines="<<std::flush; std::cerr<<refines(pic_tfm,tfm)<<std::endl;
    std::cerr<<"pic_tfm="<<pic_tfm<<"\n";
    if (refines(pic_tfm,tfm)) {
        return pic_tfm;
    } else {
        auto diff_tfm = pic_tfm-tfm; diff_tfm.clobber();
        std::cerr<<"\ndiff="<<diff_tfm << "\n\nnorm="<<norm(pic_tfm-tfm) << ", e1="<<pic_tfm.error()<<", e0="<<tfm.error()<<"\n\n";
        ARIADNE_THROW(std::runtime_error, "TaylorPolynomialIntegrator::flow_step(...)","Test flow function model "<<tfm<<" not verified by Picard iterate "<<pic_tfm<<"\n");
    }
}


ValidatedVectorMultivariateFunctionModelDP
TaylorPolynomialIntegrator::flow_step(const ValidatedVectorMultivariateFunction& vector_field,
            const ExactBoxType& state_domain,
            const StepSizeType& time_step,
            const UpperBoxType& bounding_box) const
{
    StepSizeType suggested_time_step = time_step;
    return this->flow_step(vector_field,state_domain,suggested_time_step);
#warning Check time step has not changed
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
    for(DegreeType i=0; i!=_temporal_order; ++i) {
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

    const SizeType n=dom.size();
    DoublePrecision prec;
    FloatDPError zero_err(prec);

    Vector<FloatDPError> err(n,zero_err);

    Vector<FloatDPError> rad(n+1,zero_err);
    for(SizeType i=0; i!=n; ++i) {
        rad[i] = cast_positive(max(dom[i].upper()-dmid[i].lower(),dmid[i].upper()-dom[i].lower()));
    }
    rad[n] = abs(h);

    for(SizeType i=0; i!=n; ++i) {
        for(Expansion<MultiIndex,ValidatedNumericType>::ConstIterator iter=bdphi[i].begin(); iter!=bdphi[i].end(); ++iter) {
            UniformConstReference<MultiIndex> a=iter->index();
            if(a[n]==this->_temporal_order && a[n]+this->_spacial_order==a.degree()) {
                UniformConstReference<ValidatedNumericType> rng = iter->coefficient();
                UniformConstReference<ValidatedNumericType> mid = mdphi[i][a];
                ARIADNE_ASSERT(rng.lower().raw()<=mid.lower().raw() && mid.upper().raw()<=rng.upper().raw());
                FloatDPError mag = FloatDPError(max(rng.upper()-mid.lower(),mid.upper()-rng.lower()));
                for(SizeType j=0; j!=n+1u; ++j) { mag *= pow(rad[j],static_cast<DegreeType>(a[j])); }
                err[i] += mag;
            }
        }
    }

    ExactBoxType flow_domain = join(dom,ExactIntervalType(0,h));

    ValidatedVectorMultivariateFunctionModelDP id = this->function_factory().create_identity(flow_domain);
    ValidatedVectorMultivariateFunctionModelDP res = this->function_factory().create_zeros(n,flow_domain);
    for(SizeType i=0; i!=n; ++i) {
        ValidatedScalarMultivariateFunctionModelDP res_model = res[i] + mdphi[i].expansion()[MultiIndex::zero(n+1)];
        for(SizeType j=0; j!=mdphi[i].argument_size()-1; ++j) { res_model+=mdphi[i].expansion()[MultiIndex::unit(n+1,j)]*(id[j]-ValidatedNumericType(midpoint(flow_domain[j]))); }
        SizeType j=mdphi[i].argument_size()-1u; { res_model+=mdphi[i].expansion()[MultiIndex::unit(n+1,j)]*id[j]; }
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
