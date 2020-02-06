/***************************************************************************
 *            inclusion_integrator.cpp
 *
 *  Copyright  2008-18  Luca Geretti, Pieter Collins, Sanja Zivanovic
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

#include "../function/taylor_function.hpp"
#include "../solvers/integrator.hpp"
#include "../solvers/bounder.hpp"
#include "../algebra/expansion.inl.hpp"
#include "inclusion_integrator.hpp"

namespace Ariadne {

template<> ErrorType wstar_multiplier<ZeroApproximation>() { return ErrorType(0u); }
template<> ErrorType wstar_multiplier<ConstantApproximation>() { return ErrorType(1u); }
template<> ErrorType wstar_multiplier<AffineApproximation>() { return ErrorType(5.0/3u); }
template<> ErrorType wstar_multiplier<SinusoidalApproximation>() { return ErrorType(1.3645_upper); }
template<> ErrorType wstar_multiplier<PiecewiseApproximation>() { return ErrorType(5.0/4u); }

ErrorType vstar(BoxDomainType const& inputs) {
    ErrorType result(0u);
    for (auto i : range(0,inputs.size())) {
        result = max(result,cast_positive(cast_exact(abs(inputs[i]).upper())));
    }
    return result;
}

ErrorType zeroparam_worstcase_error(ErrorConstants const& n, ErrorType const& r, PositiveFloatDPValue const& h);
ErrorType oneparam_worstcase_error(ErrorConstants const& n, ErrorType const& r, PositiveFloatDPValue const& h);
template<class R> ErrorType twoparam_worstcase_error(ErrorConstants const& n, ErrorType const& r, PositiveFloatDPValue const& h);
ErrorType zeroparam_component_error(ErrorConstants const& n, ErrorType const& v, ErrorType const& w, PositiveFloatDPValue const& h, SizeType j);
ErrorType oneparam_component_error(ErrorConstants const& n, ErrorType const& v, ErrorType const& w, PositiveFloatDPValue const& h, SizeType j);
template<class R> ErrorType twoparam_component_error(ErrorConstants const& n, ErrorType const& v, ErrorType const& w, PositiveFloatDPValue const& h, SizeType j);

template<class A, class R, EnableIf<IsSame<A,ZeroApproximation>> = dummy> ErrorType worstcase_error(ErrorConstants const& n, PositiveFloatDPValue const& h) { return zeroparam_worstcase_error(n, wstar_multiplier<A>(), h); }
template<class A, class R, EnableIf<IsSame<A,ConstantApproximation>> = dummy> ErrorType worstcase_error(ErrorConstants const& n, PositiveFloatDPValue const& h) { return oneparam_worstcase_error(n, wstar_multiplier<A>(), h); }
template<class A, class R, EnableIf<Not<Or<IsSame<A,ZeroApproximation>,IsSame<A,ConstantApproximation>>>> = dummy> ErrorType worstcase_error(ErrorConstants const& n, PositiveFloatDPValue const& h) { return twoparam_worstcase_error<R>(n, wstar_multiplier<A>(), h); }
template<class A, class R, EnableIf<IsSame<A,ZeroApproximation>> = dummy> ErrorType component_error(ErrorConstants const& n, BoxDomainType const& inputs, PositiveFloatDPValue const& h, SizeType j) { return zeroparam_component_error(n, vstar(inputs), wstar<A>(inputs), h, j); }
template<class A, class R, EnableIf<IsSame<A,ConstantApproximation>> = dummy> ErrorType component_error(ErrorConstants const& n, BoxDomainType const& inputs, PositiveFloatDPValue const& h, SizeType j) { return oneparam_component_error(n, vstar(inputs), wstar<A>(inputs), h, j); }
template<class A, class R, EnableIf<Not<Or<IsSame<A,ZeroApproximation>,IsSame<A,ConstantApproximation>>>> = dummy> ErrorType component_error(ErrorConstants const& n, BoxDomainType const& inputs, PositiveFloatDPValue const& h, SizeType j) { return twoparam_component_error<R>(n, vstar(inputs), wstar<A>(inputs), h, j); }

ErrorConstants::ErrorConstants(ErrorType const& K_, Vector<ErrorType> const& Kj_, ErrorType const& pK_, ErrorType const& pKv_, ErrorType const& pKw_, Vector<ErrorType> const& pKj_, Vector<ErrorType> const& pKjv_, Vector<ErrorType> const& pKjw_,
                               ErrorType const& L_, Vector<ErrorType> const& Lj_, ErrorType const& pL_, ErrorType const& pLv_, ErrorType const& pLw_, Vector<ErrorType> const& pLj_, Vector<ErrorType> const& pLjv_, Vector<ErrorType> const& pLjw_,
                               ErrorType const& H_, Vector<ErrorType> const& Hj_, ErrorType const& pH_, ErrorType const& pHv_, ErrorType const& pHw_, Vector<ErrorType> const& pHj_, Vector<ErrorType> const& pHjv_, Vector<ErrorType> const& pHjw_,
                               FloatDPUpperBound const& Lambda_, ErrorType const& expLambda_, ErrorType const& expL_)
 : K(K_), Kj(Kj_), pK(pK_), pKv(pKv_), pKw(pKw_), pKj(pKj_), pKjv(pKjv_), pKjw(pKjw_),
   L(L_), Lj(Lj_), pL(pL_), pLv(pLv_), pLw(pLw_), pLj(pLj_), pLjv(pLjv_), pLjw(pLjw_),
   H(H_), Hj(Hj_), pH(pH_), pHv(pHv_), pHw(pHw_), pHj(pHj_), pHjv(pHjv_), pHjw(pHjw_),
   Lambda(Lambda_), expLambda(expLambda_), expL(expL_) {
    _dimension = Kj.size();
    assert(Kj.size() == _dimension and pKj.size() == _dimension and pKjv.size() == _dimension and pKjw.size() == _dimension and
           Lj.size() == _dimension and pLj.size() == _dimension and pLjv.size() == _dimension and pLjw.size() == _dimension and
           Hj.size() == _dimension and pHj.size() == _dimension and pHjv.size() == _dimension and pHjw.size() == _dimension);
}

auto ErrorConstants::values() const
    -> Tuple<ErrorType,Vector<ErrorType>,ErrorType,ErrorType,ErrorType,Vector<ErrorType>,Vector<ErrorType>,Vector<ErrorType>,
      ErrorType,Vector<ErrorType>,ErrorType,ErrorType,ErrorType,Vector<ErrorType>,Vector<ErrorType>,Vector<ErrorType>,
      ErrorType,Vector<ErrorType>,ErrorType,ErrorType,ErrorType,Vector<ErrorType>,Vector<ErrorType>,Vector<ErrorType>,
      FloatDPUpperBound,ErrorType,ErrorType>
{
    return std::tie(this->K,this->Kj,this->pK,this->pKv,this->pKw,this->pKj,this->pKjv,this->pKjw,this->L,this->Lj,this->pL,this->pLv,this->pLw,this->pLj,this->pLjv,this->pLjw,this->H,this->Hj,this->pH,this->pHv,this->pHw,this->pHj,this->pHjv,this->pHjw,this->Lambda,this->expLambda,this->expL);
}

ErrorConstants
compute_constants(EffectiveVectorMultivariateFunction const& noise_independent_component, Vector<EffectiveVectorMultivariateFunction> const& input_derivatives, BoxDomainType const& inputs, PositiveFloatDPValue const& h, UpperBoxType const& B) {

    auto n = noise_independent_component.result_size();
    auto m = input_derivatives.size();
    DoublePrecision pr;
    ErrorType ze(pr);
    ErrorType K=ze, pK=ze, pKv=ze, pKw=ze, L=ze, pL=ze, pLv=ze, pLw=ze, H=ze, pH=ze, pHv=ze, pHw=ze;
    Vector<ErrorType> Kj(n), pKj(n), pKjv(n), pKjw(n), Lj(n), pLj(n), pLjv(n), pLjw(n), Hj(n), pHj(n), pHjv(n), pHjw(n);
    FloatDPUpperBound Lambda=ze;

    auto Df=noise_independent_component.differential(cast_singleton(B),2);
    for (auto j : range(n)) {
        auto Df_j=Df[j].expansion();
        ErrorType K_j=ze, L_j=ze, H_j=ze; FloatDPUpperBound Lambda_j=ze;
        for (auto ac : Df_j) {
            UniformReference<MultiIndex> a=ac.index();
            UniformReference<FloatDPBounds> c=ac.coefficient();
            if (a.degree()==0) {
                K_j = mag(c);
            } else if (a.degree()==1) {
                L_j += mag(c);
                if (a[j]==1) { Lambda_j += c.upper(); }
                else { Lambda_j += mag(c); }
            } else {
                assert(a.degree()==2);
                H_j += mag(c);
            }
        }
        K=max(K,K_j); L=max(L,L_j); H=max(H,H_j); Lambda=max(Lambda,Lambda_j);
        Kj[j] = K_j;
        Lj[j] = L_j;
        Hj[j] = H_j;
    }

    Matrix<ErrorType> pK_ij(m,n), pL_ij(m,n), pH_ij(m,n);

    for (auto i : range(m)) {
        auto Diff_g_i=input_derivatives[i].differential(cast_singleton(B),2);
        ErrorType pK_i=ze, pL_i=ze, pH_i=ze;
        for (auto j : range(n)) {
            auto Diff_g_ij=Diff_g_i[j].expansion();
            ErrorType g_ij=ze, Dg_ij=ze, D2g_ij=ze;
            for (auto ac : Diff_g_ij) {
                UniformReference<MultiIndex> a=ac.index();
                FloatDPBounds const& c=ac.coefficient();
                auto val = mag(c);
                if (a.degree()==0) {
                    g_ij = val;
                    pK_ij[i][j] = val;
                } else if (a.degree()==1) {
                    Dg_ij += val;
                    pL_ij[i][j] = max(pL_ij[i][j],ErrorType(val));
                } else {
                    assert(a.degree()==2);
                    D2g_ij += val;
                    pH_ij[i][j] = max(pH_ij[i][j],ErrorType(val));
                }
            }
            pK_i=max(pK_i,g_ij); pL_i=max(pL_i,Dg_ij); pH_i=max(pH_i,D2g_ij);
        }
        pK+=pK_i; pL+=pL_i; pH+=pH_i;
    }

    for (auto j : range(n)) {
        for (auto i : range(m)) {
            ErrorType Vi(abs(inputs[i]).upper());
            pKjv[j] += Vi*pK_ij[i][j]; pLjv[j] += Vi*pL_ij[i][j]; pHjv[j] += Vi*pH_ij[i][j];
        }
        pKv = max(pKv,pKjv[j]); pLv = max(pLv,pLjv[j]); pHv = max(pHv,pHjv[j]);
    }

    ErrorType expLambda = (possibly(Lambda>0)) ? ErrorType(dexp(Lambda*h)) : ErrorType(1u,pr);
    ErrorType expL = cast_positive(exp(L*h));

    return ErrorConstants(K, Kj, pK, pKv, pKw, pKj, pKjv, pKjw,
                          L, Lj, pL, pLv, pLw, pLj, pLjv, pLjw,
                          H, Hj, pH, pHv, pHw, pHj, pHjv, pHjw,
                          Lambda, expLambda, expL);
}

ErrorType zeroparam_worstcase_error(ErrorConstants const& n, ErrorType const& r, PositiveFloatDPValue const& h) {
    return min(n.pK*n.expLambda*h, (n.K*2u+n.pK)*h); }
ErrorType zeroparam_component_error(ErrorConstants const& n, ErrorType const& v, ErrorType const& w, PositiveFloatDPValue const& h, SizeType j) {
    return min(n.pK*n.expL*h*v, min( (n.Kj[j]*2u+n.pKj[j]*v)*h, n.pK*v*n.expLambda*h));
}

ErrorType oneparam_worstcase_error(ErrorConstants const& n, ErrorType const& r, PositiveFloatDPValue const& h) {
    return pow(h,2u)*((n.K+n.pK)*n.pL/3u + n.pK*2u*(n.L+n.pL)*n.expLambda); }
ErrorType oneparam_component_error(ErrorConstants const& n, ErrorType const& v, ErrorType const& w, PositiveFloatDPValue const& h, SizeType j) {
    return pow(h,2u)*(n.pLj[j]*v*(n.K+v*n.pK)/2u + (n.Lj[j]+n.pLj[j])*(v+w)*n.pK*cast_positive(cast_exact(n.expLambda-1u))/(cast_positive(cast_exact(n.Lambda))*h) ); }

template<> ErrorType twoparam_worstcase_error<AffineInputs>(ErrorConstants const& n, ErrorType const& r, PositiveFloatDPValue const& h) {
    return ((r*r+1u)*n.pL*n.pK +
            (r+1u)*h*n.pK*((n.pH*2u*r + n.H)*(n.K+r*n.pK)+pow(n.L,2u)+(n.L*3u*r+n.pL*r*r*2u)*n.pL)*n.expLambda +
            (r+1u)/6u*h*(n.K+n.pK)*((n.H*n.pK+n.L*n.pL)*3u+(n.pH*n.K+n.L*n.pL)*4u)
           )/cast_positive(+1u-h*n.L/2u-h*n.pL*r)*pow(h,2u)/4u; }
template<> ErrorType twoparam_worstcase_error<AdditiveInputs>(ErrorConstants const& n, ErrorType const& r, PositiveFloatDPValue const& h) {
    return (n.H*(n.K+n.pK)/2u +
            (pow(n.L,2u)+n.H*(n.K+r*n.pK))*n.expLambda
           )/cast_positive(+1u-h*n.L/2u)*(r+1u)*n.pK*pow(h,3u)/4u; }
template<> ErrorType twoparam_worstcase_error<SingularInput>(ErrorConstants const& n, ErrorType const& r, PositiveFloatDPValue const& h) {
    return ((r+1u)*n.pK*((n.pH*2u*r+n.H)*(n.K+r*n.pK)+pow(n.L,2u)+(n.L*3u*r+pow(r,2u)*2u*n.pL)*n.pL)*n.expLambda +
            (n.K+n.pK)/6u*((r+1u)*((n.H*n.pK+n.L*n.pL)*3u +(n.pH*n.K+n.L*n.pL)*4u) +
                           (n.pH*n.pK+pow(n.pL,2u))*8u*(r*r+1u))
           )*pow(h,3u)/4u/cast_positive(+1u-h*n.L/2u-h*n.pL*r); }
template<> ErrorType twoparam_worstcase_error<DualInputs>(ErrorConstants const& n, ErrorType const& r, PositiveFloatDPValue const& h) {
    return twoparam_worstcase_error<AffineInputs>(n, r, h); }

template<> ErrorType twoparam_component_error<AffineInputs>(ErrorConstants const& n, ErrorType const& v, ErrorType const& w, PositiveFloatDPValue const& h, SizeType j) {
    return (
        pow(h,2u) * (sqr(v)+sqr(w)) * n.pK*n.pLj[j]/2u +
        pow(h,3u) * (v+w) * (
            (n.K+v*n.pK) * ( (n.Hj[j]*n.pK+n.Lj[j]*n.pL)/8u + (n.pHj[j]*n.K+n.L*n.pLj[j])/6u ) +
            n.pK * ( n.Lj[j]*n.L + w*n.pL*n.Lj[j] + n.Hj[j]*(n.K+w*n.pK) ) * psi0(n.Lambda*h) +
            n.pK * ( n.pLj[j]*n.L + w*n.pL*n.pLj[j] +n.pHj[j]*(n.K+w*n.pK) ) * psi1(n.Lambda*h)  )
      ) / cast_positive(1u - h * (n.Lj[j]/2u+w*n.pLj[j]) );
}
template<> ErrorType twoparam_component_error<AdditiveInputs>(ErrorConstants const& n, ErrorType const& v, ErrorType const& w, PositiveFloatDPValue const& h, SizeType j) {
    return (n.Hj[j]*(n.K+v)/8u +
            (n.Lj[j]*n.L+n.Hj[j]*(n.K+w))*cast_positive(cast_exact((n.L*h*(1+exp(n.Lambda*h))/2u -exp(n.Lambda*h)+1u))/pow(cast_exact(n.Lambda*h),3u))
           )*pow(h,3u)*(v+w)/cast_positive(+1u-h*n.Lj[j]/2u); }
template<> ErrorType twoparam_component_error<SingularInput>(ErrorConstants const& n, ErrorType const& v, ErrorType const& w, PositiveFloatDPValue const& h, SizeType j) {
    return (pow(h,3u)*((n.K+v*n.pK)*(v+w)*(n.Hj[j]*n.pK+n.Lj[j]*n.pL)*3u+(n.pHj[j]*n.K*(n.K+v*n.pK)+n.L*n.pLj[j]*(n.K+w*n.pK))*4u +
            (pow(v,2)+pow(w,2))*8u*((n.K+w*n.pK)*n.pHj[j]*n.pK + (n.K+v*n.pK)*n.pL*n.pLj[j]))/24u +
            n.pK*(v+w)*((n.Lj[j]*n.L+w*n.pL*n.Lj[j]+n.Hj[j]*(n.K+w*n.pK))*cast_positive(cast_exact((n.L*h*(1+exp(n.Lambda*h))/2u -exp(n.Lambda*h)+1u))/cast_exact(pow(cast_exact(n.Lambda),3u)))) +
            n.pK*w*(v+w)*(n.pLj[j]*n.L+w*n.pL*n.pLj[j]+n.pHj[j]*(n.K+w*n.pK))*cast_positive(cast_exact((exp(n.Lambda*h)*(n.Lambda*h-1u)-pow(cast_exact(n.Lambda),2u)*pow(h,2u)/2u+1u)/pow(cast_exact(n.Lambda),3u)))
           )/cast_positive(1u-h*n.Lj[j]/2u-h*w*n.pLj[j]); }
template<> ErrorType twoparam_component_error<DualInputs>(ErrorConstants const& n, ErrorType const& v, ErrorType const& w, PositiveFloatDPValue const& h, SizeType j) {
    ErrorType acc(0u);
    for (auto l : range(0,n.Lj.size())) {
        acc += n.pLj[l]*(n.K+w*n.pK)*n.pLj[j]+n.pKj[l]*n.pHj[j]*(n.K+v*n.pK);
    }
    return (pow(h,3u)*(pow(v,2)+pow(w,2))*4u/6u*acc +
            pow(h,3u)*(n.K+v*n.pK)*(v+w)*((n.Hj[j]*n.pK+n.Lj[j]*n.pL)/8u+(n.pHj[j]*n.K+n.L*n.pLj[j])/6u) +
            n.pK*(v+w)*((n.Lj[j]*n.L+w*n.pL*n.Lj[j]+n.Hj[j]*(n.K+w*n.pK))*cast_positive(cast_exact((n.L*h*(1+exp(n.Lambda*h))/2u -exp(n.Lambda*h)+1u))/cast_exact(pow(cast_exact(n.Lambda),3u)))) +
            n.pK*w*(v+w)*(n.pLj[j]*n.L+w*n.pL*n.pLj[j]+n.pHj[j]*(n.K+w*n.pK))*cast_positive(cast_exact((exp(n.Lambda*h)*(n.Lambda*h-1u)-pow(cast_exact(n.Lambda),2u)*pow(h,2u)/2u+1u)/pow(cast_exact(n.Lambda),3u)))
           )/cast_positive(1u-h*n.Lj[j]/2u-h*w*n.pLj[j]); }

template<class A, class R> Vector<ErrorType> ApproximationErrorProcessor<A,R>::process(ErrorConstants const& n, PositiveFloatDPValue const& h) const {

    Vector<ErrorType> result_worstcase(n.dimension(),worstcase_error<A,R>(n,h));

    Vector<ErrorType> result_componentwise(n.dimension());
    for (auto j: range(n.dimension()))
        result_componentwise[j] = component_error<A,R>(n,_inputs,h,j);

    Vector<ErrorType> result(n.dimension());
    for (auto j: range(n.dimension())) {
        result[j] = min(result_worstcase[j],result_componentwise[j]);
    }

    return result;
}

template<class A, class R> Vector<ErrorType> ApproximationErrorProcessor<A,R>::process(PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    ErrorConstants norms = compute_constants(Ariadne::noise_independent_component(_f, _inputs.size()),
                                             Ariadne::input_derivatives(_f, _inputs.size()), _inputs, h, B);
    ARIADNE_LOG(7,"norms: " << norms << "\n");
    Set<Nat> input_idx;
    for (Nat i : range(_f.result_size(),_f.result_size()+_inputs.size())) { input_idx.insert(i); }
    if (is_additive_in(_f,input_idx))
        norms.pK=mag(norm(_inputs));
    return process(norms,h);
}

InclusionIntegratorHandle
InclusionIntegratorFactory::create(EffectiveVectorMultivariateFunction const& f, BoxDomainType const& inputs, InputApproximation const& approximation) const {
    if (approximation.handles(ZeroApproximation())) return InclusionIntegratorHandle(SharedPointer<InclusionIntegratorInterface>(new InclusionIntegrator<ZeroApproximation>(f,inputs,_integrator)));
    if (approximation.handles(ConstantApproximation())) return InclusionIntegratorHandle(SharedPointer<InclusionIntegratorInterface>(new InclusionIntegrator<ConstantApproximation>(f,inputs,_integrator)));
    if (approximation.handles(AffineApproximation())) return InclusionIntegratorHandle(SharedPointer<InclusionIntegratorInterface>(new InclusionIntegrator<AffineApproximation>(f,inputs,_integrator)));
    if (approximation.handles(SinusoidalApproximation())) return InclusionIntegratorHandle(SharedPointer<InclusionIntegratorInterface>(new InclusionIntegrator<SinusoidalApproximation>(f,inputs,_integrator)));
    if (approximation.handles(PiecewiseApproximation())) return InclusionIntegratorHandle(SharedPointer<InclusionIntegratorInterface>(new InclusionIntegrator<PiecewiseApproximation>(f,inputs,_integrator)));
    ARIADNE_FAIL_MSG("Unhandled input approximation " << approximation << "\n");
}


template<class A> Bool
InclusionIntegrator<A>::operator==(const InclusionIntegratorInterface& rhs) const {
    return instance_of<InclusionIntegrator<A>>(&rhs);
}

template<class A> Bool
InclusionIntegrator<A>::operator<(const InclusionIntegratorInterface& rhs) const {
    return this->index() < rhs.index();
}

template<class A> List<ValidatedVectorMultivariateFunctionModelType>
InclusionIntegrator<A>::reach(BoxDomainType const& domx, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const {

    TimeStepType new_t = lower_bound(t+h);

    Interval<TimeStepType> domt(t,new_t);

    auto e=this->compute_errors(h,B);
    ARIADNE_LOG(6,"approximation errors:"<<e<<"\n");
    auto doma = this->build_parameter_domain(_inputs);

    auto w = build_w_functions<A>(domt,doma,_f.result_size(),_inputs.size());
    ARIADNE_LOG(6,"w:"<<w<<"\n");
    auto Fw = substitute_v_with_w(_f, w);
    ARIADNE_LOG(6,"Fw:"<<Fw<<"\n");
    auto phi = this->_integrator->flow_step(Fw,domx,domt,doma,B);
    add_errors(phi,e);

    List<ValidatedVectorMultivariateFunctionModelType> result;
    result.append(this->build_reach_function(evolve_function, phi, t, new_t));

    return result;
}

template<class A> Vector<EffectiveScalarMultivariateFunction> InclusionIntegrator<A>::build_secondhalf_piecewise_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) const {
    auto zero = EffectiveScalarMultivariateFunction::zero(n+1+2*m);
    auto one = EffectiveScalarMultivariateFunction::constant(n+1+2*m,1_z);

    auto result = Vector<EffectiveScalarMultivariateFunction>(m);
    for (auto i : range(m)) {
        auto Vi = ExactNumber(doma[i].upper());
        auto p0 = EffectiveScalarMultivariateFunction::coordinate(n+1+2*m,n+1+i);
        auto p1 = EffectiveScalarMultivariateFunction::coordinate(n+1+2*m,n+1+m+i);
        result[i] = (definitely (doma[i].upper() == 0.0_exact) ? zero : p0+(one-p0*p0/Vi/Vi)*p1);
    }
    return result;
}

template<class A> ValidatedVectorMultivariateFunctionModelDP InclusionIntegrator<A>::build_secondhalf_piecewise_reach_function(
        ValidatedVectorMultivariateFunctionModelDP const& evolve_function, ValidatedVectorMultivariateFunctionModelDP const& Phi, TimeStepType const& t,
        TimeStepType const& new_t) const {

    // Evolve function is e(x,a,b) at s; Flow is phi(x,h,b)
    // Want (x,t,a,b):->phi(e(x,a,b),t-s,b))

    auto swp = dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(Phi.reference()).properties();

    SizeType n=evolve_function.result_size();
    SizeType b=Phi.argument_size()-(n+1);

    SizeType a=evolve_function.argument_size()-n-b;

    BoxDomainType X=evolve_function.domain()[range(0,n)];
    BoxDomainType PA=evolve_function.domain()[range(n,n+a)];
    BoxDomainType PB=Phi.domain()[range(n+1,n+1+b)];

    auto Tau=IntervalDomainType(t,new_t);
    BoxDomainType XTP = join(X,Tau,PA,PB);
    ValidatedVectorMultivariateTaylorFunctionModelDP xf=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(0,n),swp);
    ValidatedScalarMultivariateTaylorFunctionModelDP tf=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(XTP,n,swp);
    ValidatedVectorMultivariateTaylorFunctionModelDP af=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(n+1,n+1+a),swp);
    ValidatedVectorMultivariateTaylorFunctionModelDP bf=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(n+1+a,n+1+a+b),swp);

    ValidatedVectorMultivariateTaylorFunctionModelDP ef=compose(evolve_function,join(xf,af,bf));

    return compose(Phi,join(ef,tf,bf));
}

template<class A> ValidatedVectorMultivariateFunctionModelDP InclusionIntegrator<A>::build_reach_function(
        ValidatedVectorMultivariateFunctionModelDP const& evolve_function, ValidatedVectorMultivariateFunctionModelDP const& Phi, TimeStepType const& t,
        TimeStepType const& new_t) const {

    // Evolve function is e(x,a) at s; flow is phi(x,h,b)
    // Want (x,t,a,b):->phi(e(x,a),t-s,b))

    auto swp = dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(Phi.reference()).properties();

    SizeType n=evolve_function.result_size();

    SizeType a=evolve_function.argument_size()-n;
    SizeType b=Phi.argument_size()-(n+1);

    BoxDomainType X=evolve_function.domain()[range(0,n)];
    BoxDomainType PA=evolve_function.domain()[range(n,n+a)];
    BoxDomainType PB=Phi.domain()[range(n+1,n+1+b)];

    auto Tau=IntervalDomainType(t,new_t);
    BoxDomainType XTP = join(X,Tau,PA,PB);
    ValidatedVectorMultivariateTaylorFunctionModelDP xf=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(0,n),swp);
    ValidatedScalarMultivariateTaylorFunctionModelDP tf=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(XTP,n,swp);
    ValidatedVectorMultivariateTaylorFunctionModelDP af=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(n+1,n+1+a),swp);
    ValidatedVectorMultivariateTaylorFunctionModelDP bf=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(n+1+a,n+1+a+b),swp);

    ValidatedVectorMultivariateTaylorFunctionModelDP ef=compose(evolve_function,join(xf,af));

    return compose(Phi,join(ef,tf,bf));
}

template<class A> ValidatedVectorMultivariateFunctionModelDP InclusionIntegrator<A>::evolve(ValidatedVectorMultivariateFunctionModelDP const& reach_function, TimeStepType const& t) const {
    return partial_evaluate(reach_function,reach_function.result_size(),t);
}

Void add_errors(ValidatedVectorMultivariateFunctionModelDP& phi, Vector<ErrorType> const& e) {
    assert(phi.result_size()==e.size());
    ValidatedVectorMultivariateTaylorFunctionModelDP& tphi = dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP&>(phi.reference());
    for (auto i : range(e.size())) {
        tphi[i].add_error(e[i]);
    }
}

EffectiveVectorMultivariateFunction substitute_v_with_w(EffectiveVectorMultivariateFunction const& F, Vector<EffectiveScalarMultivariateFunction> const& w) {

    auto n = F.result_size();
    auto m = w.size();
    auto p = w[0].argument_size();

    auto coordinates = EffectiveVectorMultivariateFunction::coordinates(p);

    auto substitution = EffectiveVectorMultivariateFunction::zeros(n+m,p);
    for (auto i : range(n)) {
        substitution.set(i,coordinates[i]);
    }
    for (auto i : range(m)) {
        substitution.set(n+i,w[i]);
    }

    return compose(F,substitution);
}


template<class A> Pair<StepSizeType,UpperBoxType> InclusionIntegrator<A>::flow_bounds(BoxDomainType const& domx, BoxDomainType const& doma, StepSizeType const& hsug) const {
    return EulerBounder().compute(_f,domx,doma,hsug);
}

template<class A> BoxDomainType InclusionIntegrator<A>::build_parameter_domain(BoxDomainType const& V) const {
    BoxDomainType result(0u);
    for (Nat i=0; i<this->_num_params_per_input; ++i)
        result = product(result,V);
    return result;
}

template<> Vector<EffectiveScalarMultivariateFunction> build_w_functions<ZeroApproximation>(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) {
    auto result = Vector<EffectiveScalarMultivariateFunction>(m);
    for (auto i : range(0,m))
        result[i] = EffectiveScalarMultivariateFunction::zero(n+1);
    return result;
}


template<> Vector<EffectiveScalarMultivariateFunction> build_w_functions<ConstantApproximation>(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) {
    auto result = Vector<EffectiveScalarMultivariateFunction>(m);
    for (auto i : range(0,m))
        result[i] = EffectiveScalarMultivariateFunction::coordinate(n+1+m,n+1+i);
    return result;
}


template<> Vector<EffectiveScalarMultivariateFunction> build_w_functions<AffineApproximation>(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) {
    auto zero = EffectiveScalarMultivariateFunction::zero(n+1+2*m);
    auto one = EffectiveScalarMultivariateFunction::constant(n+1+2*m,1_z);
    auto three = EffectiveScalarMultivariateFunction::constant(n+1+2*m,3_z);
    auto t = EffectiveScalarMultivariateFunction::coordinate(n+1+2*m,n);
    auto tk = EffectiveScalarMultivariateFunction::constant(n+1+2*m,domt.lower());
    auto hc = EffectiveScalarMultivariateFunction::constant(n+1+2*m,domt.width());

    auto result = Vector<EffectiveScalarMultivariateFunction>(m);
    for (auto i : range(m)) {
        auto Vi = ExactNumber(doma[i].upper());
        auto p0 = EffectiveScalarMultivariateFunction::coordinate(n+1+2*m,n+1+i);
        auto p1 = EffectiveScalarMultivariateFunction::coordinate(n+1+2*m,n+1+m+i);
        result[i] = (definitely (doma[i].upper() == 0.0_exact) ? zero : p0+three*(one-p0*p0/Vi/Vi)*p1*(t-tk-hc/2)/hc);
    }
    return result;
}


template<> Vector<EffectiveScalarMultivariateFunction> build_w_functions<SinusoidalApproximation>(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) {
    auto zero = EffectiveScalarMultivariateFunction::zero(n+1+2*m);
    auto one = EffectiveScalarMultivariateFunction::constant(n+1+2*m,1_z);
    auto pgamma = EffectiveScalarMultivariateFunction::constant(n+1+2*m,1.1464_dec);
    auto gamma = EffectiveScalarMultivariateFunction::constant(n+1+2*m,4.162586_dec);
    auto t = EffectiveScalarMultivariateFunction::coordinate(n+1+2*m,n);
    auto tk = EffectiveScalarMultivariateFunction::constant(n+1+2*m,domt.lower());
    auto hc = EffectiveScalarMultivariateFunction::constant(n+1+2*m,domt.width());

    auto result = Vector<EffectiveScalarMultivariateFunction>(m);
    for (auto i : range(m)) {
        auto Vi = ExactNumber(doma[i].upper());
        auto p0 = EffectiveScalarMultivariateFunction::coordinate(n+1+2*m,n+1+i);
        auto p1 = EffectiveScalarMultivariateFunction::coordinate(n+1+2*m,n+1+m+i);
        result[i] = (definitely (doma[i].upper() == 0.0_exact) ? zero : p0+(one-p0*p0/Vi/Vi)*pgamma*p1*sin((t-tk-hc/2)*gamma/hc));
    }
    return result;
}


template<> Vector<EffectiveScalarMultivariateFunction> build_w_functions<PiecewiseApproximation>(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) {
    auto zero = EffectiveScalarMultivariateFunction::zero(n+1+2*m);
    auto one = EffectiveScalarMultivariateFunction::constant(n+1+2*m,1_z);

    auto result = Vector<EffectiveScalarMultivariateFunction>(m);
    for (auto i : range(m)) {
        auto Vi = ExactNumber(doma[i].upper());
        auto p0 = EffectiveScalarMultivariateFunction::coordinate(n+1+2*m,n+1+i);
        auto p1 = EffectiveScalarMultivariateFunction::coordinate(n+1+2*m,n+1+m+i);
        result[i] = (definitely (doma[i].upper() == 0.0_exact) ? zero : p0-(one-p0*p0/Vi/Vi)*p1);
    }
    return result;
}

template<> List<ValidatedVectorMultivariateFunctionModelType>
InclusionIntegrator<PiecewiseApproximation>::reach(BoxDomainType const& domx, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const {

    List<ValidatedVectorMultivariateFunctionModelType> result;

    auto n = _f.result_size();
    auto m = _inputs.size();

    TimeStepType intermediate_t = lower_bound(t+hlf(h));
    TimeStepType new_t = lower_bound(t+h);

    Interval<TimeStepType> domt_first(t,intermediate_t);
    Interval<TimeStepType> domt_second(intermediate_t,new_t);
    auto doma = this->build_parameter_domain(_inputs);

    auto e=this->compute_errors(h,B);
    ARIADNE_LOG(6,"approximation errors:"<<e<<"\n");

    auto w_hlf = build_w_functions<PiecewiseApproximation>(domt_first,doma,n,m);
    ARIADNE_LOG(6,"w_hlf:"<<w_hlf<<"\n");
    auto Fw_hlf = substitute_v_with_w(_f, w_hlf);
    ARIADNE_LOG(6,"Fw_hlf:" << Fw_hlf << "\n");

    auto phi_hlf = this->_integrator->flow_step(Fw_hlf,domx,domt_first,doma,B);
    auto intermediate_reach=this->build_reach_function(evolve_function, phi_hlf, t, intermediate_t);

    result.append(intermediate_reach);

    auto intermediate_evolve=this->evolve(intermediate_reach,intermediate_t);

    auto domx_second = cast_exact_box(intermediate_evolve.range());

    auto w = this->build_secondhalf_piecewise_w_functions(domt_second,doma,n,m);
    ARIADNE_LOG(6,"w:"<<w<<"\n");
    auto Fw = substitute_v_with_w(_f, w);
    ARIADNE_LOG(6,"Fw:"<<Fw<<"\n");
    auto phi = this->_integrator->flow_step(Fw,domx_second,domt_second,doma,B);
    add_errors(phi,e);

    result.append(this->build_secondhalf_piecewise_reach_function(intermediate_evolve, phi, intermediate_t, new_t));

    return result;
}

} // namespace Ariadne;
