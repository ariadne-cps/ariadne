/***************************************************************************
 *            differential_inclusion.cpp
 *
 *  Copyright  2008-17  Pieter Collins, Sanja Zivanovic
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

#include "differential_inclusion.hpp"
#include "function/taylor_function.hpp"
#include "solvers/integrator.hpp"

namespace Ariadne {

#define ARIADNE_LOG_PRINT(level, expr) { ARIADNE_LOG(level,#expr << "=" << (expr) << "\n"); }

Box<UpperIntervalType> apply(VectorFunction<ValidatedTag>const& f, const Box<ExactIntervalType>& bx) {
    return apply(f,Box<UpperIntervalType>(bx));
}

Boolean is_identity_matrix(Vector<ValidatedVectorFunction> const& g, UpperBoxType const& B) {

    for (SizeType m : range(g.size())) {
        for (SizeType n: range(g[m].result_size())) {
            if (m == n) {
                if (definitely(g[m][n].evaluate(cast_singleton(B)) != 1.0_exact))
                    return false;
            } else {
                if (definitely(g[m][n].evaluate(cast_singleton(B)) != 0.0_exact))
                    return false;
            }
        }
    }
    return true;
}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
compute_norms_LC(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) {

    DoublePrecision pr;
    FloatDPError ze(pr);
    FloatDPError K=ze, Kp=ze, L=ze, Lp=ze, H=ze, Hp=ze; FloatDPUpperBound Lambda=ze;

    auto Df=f.differential(cast_singleton(B),1);
    for (auto n : range(f.result_size())) {
        auto Df_n=Df[n].expansion();
        FloatDPError K_n=ze, L_n=ze; FloatDPUpperBound Lambda_n=ze;
        for (auto ac : Df_n) {
            MultiIndex const& a=ac.index();
            FloatDPBounds const& c=ac.coefficient();
            if (a.degree()==0) {
                K_n += mag(c);
            } else {
                assert(a.degree()==1);
                if (a[n]==1) { Lambda_n += c.upper(); }
                else { Lambda_n += mag(c); }
            }
        }
        K=max(K,K_n); L=max(L,L_n); Lambda=max(Lambda,Lambda_n);
    }

    for (auto m : range(g.size())) {
        auto g_m=g[m];
        auto Dg_m=g_m.differential(cast_singleton(B),1);
        FloatDPError Vm(abs(V[m]).upper());
        FloatDPError Kp_m=ze;
        for (auto n : range(g_m.result_size())) {
            auto Dg_mn=Dg_m[n].expansion();
            FloatDPError Kp_mn=ze;
            for (auto ac : Dg_mn) {
                MultiIndex const& a=ac.index();
                FloatDPBounds const& c=ac.coefficient();
                if (a.degree()==0) {
                    Kp_mn += mag(c);
                }
            }
            Kp_m=max(Kp_m,Kp_mn);
        }

        Kp+=Vm*Kp_m;
    }

    FloatDPError expLambda = (possibly(Lambda>0)) ? FloatDPError(dexp(Lambda*h)) : FloatDPError(1u,pr);

    return std::tie(K,Kp,L,Lp,H,Hp,expLambda);
}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
compute_norms_LC_additive(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) {
    //! For additive noise, K'=|V| while Lp = Hp = 0
    FloatDPError Kp=mag(norm(V));
    FloatDPError Lp, Hp;

    auto Df=f.differential(cast_singleton(B),1);
    DoublePrecision pr;
    FloatDPError ze(pr);
    FloatDPError K=ze, L=ze, H=ze; FloatDPUpperBound Lambda=ze;
    for (auto n : range(f.result_size())) {
        auto Df_n=Df[n].expansion();
        FloatDPError K_n=ze; FloatDPUpperBound Lambda_n=ze;
        for (auto ac : Df_n) {
            MultiIndex const& a=ac.index();
            FloatDPBounds const& c=ac.coefficient();
            if (a.degree()==0) {
                K_n += mag(c);
            } else {
                assert(a.degree()==1);
                if (a[n]==1) { Lambda_n += c.upper(); }
                else { Lambda_n += mag(c); }
            }
        }
        K=max(K,K_n); Lambda=max(Lambda,Lambda_n);
    }

    FloatDPError expLambda = (possibly(Lambda>0)) ? FloatDPError(dexp(Lambda*h)) : FloatDPError(1u,pr);

    return std::tie(K,Kp,L,Lp,H,Hp,expLambda);
}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
compute_norms_C1_additive(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) {
    //! For additive noise, K'=|V| while Lp = Hp = 0
    FloatDPError Kp=mag(norm(V));
    FloatDPError Lp, Hp;

    auto Df=f.differential(cast_singleton(B),1);
    DoublePrecision pr;
    FloatDPError ze(pr);
    FloatDPError K=ze, L=ze, H=ze; FloatDPUpperBound Lambda=ze;
    for (auto n : range(f.result_size())) {
        auto Df_n=Df[n].expansion();
        FloatDPError K_n=ze, L_n=ze; FloatDPUpperBound Lambda_n=ze;
        for (auto ac : Df_n) {
            MultiIndex const& a=ac.index();
            FloatDPBounds const& c=ac.coefficient();
            if (a.degree()==0) {
                K_n += mag(c);
            } else {
                assert(a.degree()==1);
                L_n += mag(c);
                if (a[n]==1) { Lambda_n += c.upper(); }
                else { Lambda_n += mag(c); }
            }
        }
        K=max(K,K_n); L=max(L,L_n); Lambda=max(Lambda,Lambda_n);
    }

    FloatDPError expLambda = (possibly(Lambda>0)) ? FloatDPError(dexp(Lambda*h)) : FloatDPError(1u,pr);

    return std::tie(K,Kp,L,Lp,H,Hp,expLambda);
}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
compute_norms_C1(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) {

    DoublePrecision pr;
    FloatDPError ze(pr);
    FloatDPError K=ze, Kp=ze, L=ze, Lp=ze, H=ze, Hp=ze; FloatDPUpperBound Lambda=ze;

    auto Df=f.differential(cast_singleton(B),1);
    for (auto n : range(f.result_size())) {
        auto Df_n=Df[n].expansion();
        FloatDPError K_n=ze, L_n=ze; FloatDPUpperBound Lambda_n=ze;
        for (auto ac : Df_n) {
            MultiIndex const& a=ac.index();
            FloatDPBounds const& c=ac.coefficient();
            if (a.degree()==0) {
                K_n += mag(c);
            } else {
                assert(a.degree()==1);
                L_n += mag(c);
                if (a[n]==1) { Lambda_n += c.upper(); }
                else { Lambda_n += mag(c); }
            }
        }
        K=max(K,K_n); L=max(L,L_n); Lambda=max(Lambda,Lambda_n);
    }

    for (auto m : range(g.size())) {
        auto g_m=g[m];
        auto Dg_m=g_m.differential(cast_singleton(B),1);
        FloatDPError Vm(abs(V[m]).upper());
        FloatDPError Kp_m=ze, Lp_m=ze, Hp_m=ze;
        for (auto n : range(g_m.result_size())) {
            auto Dg_mn=Dg_m[n].expansion();
            FloatDPError Kp_mn=ze, Lp_mn=ze, Hp_mn=ze;
            for (auto ac : Dg_mn) {
                MultiIndex const& a=ac.index();
                FloatDPBounds const& c=ac.coefficient();
                if (a.degree()==0) {
                    Kp_mn += mag(c);
                } else {
                    assert(a.degree()==1);
                    Lp_mn += mag(c);
                }
            }
            Kp_m=max(Kp_m,Kp_mn); Lp_m=max(Lp_m,Lp_mn); Hp_m=max(Hp,Hp_mn);
        }

        Kp+=Vm*Kp_m; Lp+=Vm*Lp_m; Hp+=Vm*Hp_m;
    }

    FloatDPError expLambda = (possibly(Lambda>0)) ? FloatDPError(dexp(Lambda*h)) : FloatDPError(1u,pr);

    return std::tie(K,Kp,L,Lp,H,Hp,expLambda);
}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
compute_norms_C2(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) {

    DoublePrecision pr;
    FloatDPError ze(pr);
    FloatDPError K=ze, Kp=ze, L=ze, Lp=ze, H=ze, Hp=ze; FloatDPUpperBound Lambda=ze;

    auto Df=f.differential(cast_singleton(B),2);
    for (auto n : range(f.result_size())) {
        auto Df_n=Df[n].expansion();
        FloatDPError K_n=ze, L_n=ze, H_n=ze; FloatDPUpperBound Lambda_n=ze;
        for (auto ac : Df_n) {
            MultiIndex const& a=ac.index();
            FloatDPBounds const& c=ac.coefficient();
            if (a.degree()==0) {
                K_n += mag(c);
            } else if (a.degree()==1) {
                L_n += mag(c);
                if (a[n]==1) { Lambda_n += c.upper(); }
                else { Lambda_n += mag(c); }
            } else {
                assert(a.degree()==2);
                H_n += mag(c);
            }
        }
        K=max(K,K_n); L=max(L,L_n); H=max(H,H_n); Lambda=max(Lambda,Lambda_n);
    }

    for (auto m : range(g.size())) {
        auto g_m=g[m];
        auto Dg_m=g_m.differential(cast_singleton(B),2);
        FloatDPError Vm(abs(V[m]).upper());
        FloatDPError Kp_m=ze, Lp_m=ze, Hp_m=ze;
        for (auto n : range(g_m.result_size())) {
            auto Dg_mn=Dg_m[n].expansion();
            FloatDPError Kp_mn=ze, Lp_mn=ze, Hp_mn=ze;
            for (auto ac : Dg_mn) {
                MultiIndex const& a=ac.index();
                FloatDPBounds const& c=ac.coefficient();
                if (a.degree()==0) {
                    Kp_mn += mag(c);
                } else if (a.degree()==1) {
                    Lp_mn += mag(c);
                } else {
                    assert(a.degree()==2);
                    Hp_mn += mag(c);
                }
            }
            Kp_m=max(Kp_m,Kp_mn); Lp_m=max(Lp_m,Lp_mn); Hp_m=max(Hp,Hp_mn);
        }

        Kp+=Vm*Kp_m; Lp+=Vm*Lp_m; Hp+=Vm*Hp_m;
    }

    FloatDPError expLambda = (possibly(Lambda>0)) ? FloatDPError(dexp(Lambda*h)) : FloatDPError(1u,pr);

    return std::tie(K,Kp,L,Lp,H,Hp,expLambda);
}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
compute_norms_C2_additive(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) {
    //! For additive noise, K'=|V| while Lp = Hp = 0
    FloatDPError Kp=mag(norm(V));
    FloatDPError Lp, Hp;

    auto Df=f.differential(cast_singleton(B),2);
    DoublePrecision pr;
    FloatDPError ze(pr);
    FloatDPError K=ze, L=ze, H=ze; FloatDPUpperBound Lambda=ze;
    for (auto n : range(f.result_size())) {
        auto Df_n=Df[n].expansion();
        FloatDPError K_n=ze, L_n=ze, H_n=ze; FloatDPUpperBound Lambda_n=ze;
        for (auto ac : Df_n) {
            MultiIndex const& a=ac.index();
            FloatDPBounds const& c=ac.coefficient();
            if (a.degree()==0) {
                K_n += mag(c);
            } else if (a.degree()==1) {
                L_n += mag(c);
                if (a[n]==1) { Lambda_n += c.upper(); }
                else { Lambda_n += mag(c); }
            } else {
                assert(a.degree()==2);
                H_n += mag(c);
            }
        }
        K=max(K,K_n); L=max(L,L_n); H=max(H,H_n); Lambda=max(Lambda,Lambda_n);
    }

    FloatDPError expLambda = (possibly(Lambda>0)) ? FloatDPError(dexp(Lambda*h)) : FloatDPError(1u,pr);

    return std::tie(K,Kp,L,Lp,H,Hp,expLambda);
}

InclusionErrorProcessor::InclusionErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
    : _f(f), _g(g), _V(V), _h(h), _B(B) { }

ErrorType InclusionErrorProcessor::process() const {

    FloatDPError K, Kp, L, Lp, H, Hp, expLambda;

    std::tie(K,Kp,L,Lp,H,Hp,expLambda) = compute_norms(_f,_g,_V,_h,_B);

    return compute_error(K,Kp,L,Lp,H,Hp,expLambda,_h);
}

ZeroErrorProcessor::ZeroErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
        : InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
ZeroErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    return compute_norms_LC(f,g,V,h,B);
}

ErrorType ZeroErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {
    FloatDPError result1 = Kp*expLambda*h;
    FloatDPError result2 = (K*2u+Kp)*h;
    return min(result1,result2);
}

AdditiveZeroErrorProcessor::AdditiveZeroErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
: InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
AdditiveZeroErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
return compute_norms_LC_additive(f,g,V,h,B);
}

ErrorType AdditiveZeroErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {
    FloatDPError result1 = Kp*expLambda*h;
    FloatDPError result2 = (K*2u+Kp)*h;
    return min(result1,result2);
}

ConstantErrorProcessor::ConstantErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
        : InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
ConstantErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    return compute_norms_C2(f,g,V,h,B);
}

ErrorType ConstantErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {
    FloatDPError result = (pow(h,2u)*(Kp*Lp*expLambda*2u + Lp*(K+Kp)/3u + Kp*L/2u)+ pow(h,3u)*Kp*(L*Lp + L*L + H*(K+Kp))/2u*expLambda)/cast_positive(1u-(h*L/2u));
    return result;
}

/*

ErrorType ConstantErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {
    FloatDPError result = pow(h,2u)*((K+Kp)*Lp/3u + Kp*(L+Lp)*expLambda*2u);
    return result;
}
*/


AdditiveConstantErrorProcessor::AdditiveConstantErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
    : InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
AdditiveConstantErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    return compute_norms_C1_additive(f,g,V,h,B);
}

ErrorType AdditiveConstantErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {
    FloatDPError result = pow(h,2u)*(Kp*L*expLambda);
    return result;
}


PiecewiseErrorProcessor::PiecewiseErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
        : InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
PiecewiseErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    return compute_norms_C2(f,g,V,h,B);
}

ErrorType PiecewiseErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {

    FloatDPError r(5.0/4u);
    FloatDPError result = ((r*r+1u)*Lp*Kp + (r+1u)*h*Kp*((Hp*2u*r + H)*(K+r*Kp)+L*L+(L*3u*r+Lp*r*r*2u)*Lp)*expLambda + (r+1u)/6u*h*(K+Kp)*((H*Kp+L*Lp)*3u+(Hp*K+L*Lp)*4u))/cast_positive(1u-h*L/2u-h*Lp*r)*pow(h,2u)/4u;

    return result;
}

SingleInputPiecewiseErrorProcessor::SingleInputPiecewiseErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
        : InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
SingleInputPiecewiseErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    return compute_norms_C2(f,g,V,h,B);
}

ErrorType SingleInputPiecewiseErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {

    FloatDPError r(5.0/4u);
    FloatDPError result = ((r+1u)*Kp*((Hp*2u*r+H)*(K+r*Kp)+L*L+(L*3u*r+Lp*r*r*2u)*Lp)*expLambda + (r+1u)/6u*(K+Kp)*((r+1u)*((H*Kp+L*Lp)*3u +(Hp*K+L*Lp)*4u) + (Hp*Kp+Lp*Lp)*8u*(r*r+1u)))/cast_positive(1u-h*L/2u-h*Lp*r)*pow(h,3u)/4u;

    return result;
}


AdditivePiecewiseErrorProcessor::AdditivePiecewiseErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
        : InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
AdditivePiecewiseErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    return compute_norms_C2_additive(f,g,V,h,B);
}

ErrorType AdditivePiecewiseErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {

    FloatDPError r(5.0/4u);
    FloatDPError result = (Kp*(H*(K+r*Kp)+L*L)*expLambda + (K+Kp)*H*Kp/2u)/cast_positive(1u-h*L/2u)*(r+1u)*pow(h,3u)/4u;

    return result;
}


AffineErrorProcessor::AffineErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
    : InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
AffineErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    return compute_norms_C2(f,g,V,h,B);
}

ErrorType AffineErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {

    FloatDPError r(5.0/3u);
    FloatDPError result = ((r*r+1u)*Lp*Kp + (r+1u)*h*Kp*((Hp*2u*r + H)*(K+r*Kp)+L*L+(L*3u*r+Lp*r*r*2u)*Lp)*expLambda + (r+1u)/6u*h*(K+Kp)*((H*Kp+L*Lp)*3u+(Hp*K+L*Lp)*4u))/cast_positive(1u-h*L/2u-h*Lp*r)*pow(h,2u)/4u;

    return result;
}

SingleInputAffineErrorProcessor::SingleInputAffineErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
    : InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
SingleInputAffineErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    return compute_norms_C2(f,g,V,h,B);
}

ErrorType SingleInputAffineErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {

    FloatDPError r(5.0/3u);
    FloatDPError result = ((r+1u)*Kp*((Hp*2u*r+H)*(K+r*Kp)+L*L+(L*3u*r+Lp*r*r*2u)*Lp)*expLambda + (r+1u)/6u*(K+Kp)*((r+1u)*((H*Kp+L*Lp)*3u +(Hp*K+L*Lp)*4u) + (Hp*Kp+Lp*Lp)*8u*(r*r+1u)))/cast_positive(1u-h*L/2u-h*Lp*r)*pow(h,3u)/4u;

    return result;
}


AdditiveAffineErrorProcessor::AdditiveAffineErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
        : InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
AdditiveAffineErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    return compute_norms_C2_additive(f,g,V,h,B);
}

ErrorType AdditiveAffineErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {

    FloatDPError r(5.0/3u);
    FloatDPError result = (Kp*(H*(K+r*Kp)+L*L)*expLambda + (K+Kp)*H*Kp/2u)/cast_positive(1u-h*L/2u)*(r+1u)*pow(h,3u)/4u;

    return result;
}


SinusoidalErrorProcessor::SinusoidalErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
        : InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
SinusoidalErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    return compute_norms_C2(f,g,V,h,B);
}

ErrorType SinusoidalErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {

    FloatDPError r(1.0093);
    FloatDPError result = ((r*r+1u)*Lp*Kp + (r+1u)*h*Kp*((Hp*2u*r + H)*(K+r*Kp)+L*L+(L*3u*r+Lp*r*r*2u)*Lp)*expLambda + (r+1u)/6u*h*(K+Kp)*((H*Kp+L*Lp)*3u+(Hp*K+L*Lp)*4u))/cast_positive(1u-h*L/2u-h*Lp*r)*pow(h,2u)/4u;

    return result;
}

SingleInputSinusoidalErrorProcessor::SingleInputSinusoidalErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
        : InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
SingleInputSinusoidalErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    return compute_norms_C2(f,g,V,h,B);
}

ErrorType SingleInputSinusoidalErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {

    FloatDPError r(1.0093);
    FloatDPError result = ((r+1u)*Kp*((Hp*2u*r+H)*(K+r*Kp)+L*L+(L*3u*r+Lp*r*r*2u)*Lp)*expLambda + (r+1u)/6u*(K+Kp)*((r+1u)*((H*Kp+L*Lp)*3u +(Hp*K+L*Lp)*4u) + (Hp*Kp+Lp*Lp)*8u*(r*r+1u)))/cast_positive(1u-h*L/2u-h*Lp*r)*pow(h,3u)/4u;

    return result;
}


AdditiveSinusoidalErrorProcessor::AdditiveSinusoidalErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B)
        : InclusionErrorProcessor(f,g,V,h,B) {}

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError>
AdditiveSinusoidalErrorProcessor::compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    return compute_norms_C2_additive(f,g,V,h,B);
}

ErrorType AdditiveSinusoidalErrorProcessor::compute_error(FloatDPError const& K,FloatDPError const& Kp,FloatDPError const& L,FloatDPError const& Lp,FloatDPError const& H,FloatDPError const& Hp,FloatDPError const& expLambda,PositiveFloatDPValue const& h) const {

    FloatDPError r(1.0093);
    FloatDPError result = (Kp*(H*(K+r*Kp)+L*L)*expLambda + (K+Kp)*H*Kp/2u)/cast_positive(1u-h*L/2u)*(r+1u)*pow(h,3u)/4u;

    return result;
}


ValidatedVectorTaylorFunctionModelDP get_time_derivative(ValidatedVectorTaylorFunctionModelDP phi, ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, ValidatedVectorTaylorFunctionModelDP wf) {

    auto n=f.result_size();
    auto m=g.size();

    ValidatedVectorTaylorFunctionModelDP result(n);

    for (auto i : range(n)) {
        result[i] = compose(f[i], phi);
        for (auto j : range(m)) {
            result[i] = result[i] + compose(g[j][i],phi) * wf[j];
        }
    }

    return result;
}


InclusionIntegrator::InclusionIntegrator(List<SharedPointer<InclusionIntegratorApproximation>> approximations, SweeperDP sweeper, StepSize step_size)
    : _approximations(approximations)
    , _sweeper(sweeper)
    , _step_size(step_size)
    , _number_of_steps_between_simplifications(8)
    , _number_of_variables_to_keep(4)
{
    assert(approximations.size()>0);
}

List<ValidatedVectorFunctionModelDP> InclusionIntegrator::flow(ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, BoxDomainType V, BoxDomainType X0, Real tmax) {
    ARIADNE_LOG(1,"\nf:"<<f<<"\ng:"<<g<<"\nV:"<<V<<"\nX0:"<<X0<<"\ntmax:"<<tmax<<"\n");

    // Ensure all arguments have the correct size;
    auto n=X0.size();
    DoublePrecision pr;
    assert(f.result_size()==n);
    assert(f.argument_size()==n);
    assert(V.size()==g.size());

    PositiveFloatDPValue hsug(this->_step_size);

    ValidatedVectorFunctionModelDP evolve_function = ValidatedVectorTaylorFunctionModelDP::identity(X0,this->_sweeper);
    auto t=PositiveFloatDPValue(0.0);

    List<ValidatedVectorFunctionModelDP> result;

    auto step = 0;

    FloatDPError total_error(0);
    while (possibly(t<FloatDPBounds(tmax,pr))) {
        ARIADNE_LOG(2,"step#:"<<step<<", t:"<<t<<", hsug:"<<hsug << "\n");
        if(possibly(t+hsug>FloatDPBounds(tmax,pr))) {  //FIXME: Check types for timing;
            hsug=cast_positive(cast_exact((tmax-t).upper()));
        }

        ARIADNE_LOG(3,"n. of parameters="<<evolve_function.argument_size()<<"\n");

        auto D = cast_exact_box(evolve_function.range());
        UpperBoxType B;
        PositiveFloatDPValue h;
        std::tie(h,B)=this->flow_bounds(f,g,V,D,hsug);
        ARIADNE_LOG(2,"flow bounds = "<<B<<" (using h = " << h << ")\n");

        PositiveFloatDPValue new_t=cast_positive(cast_exact((t+h).lower()));

        ValidatedVectorFunctionModelDP reach_function;
        ValidatedVectorFunctionModelDP best_reach_function, best_evolve_function;
        DIApproximationKind best;
        FloatDP best_total_diameter(0);

        SizeType i = 0;
        for (auto i : range(_approximations.size())) {
            this->_approximation = _approximations.at(i);
            ARIADNE_LOG(4,"checking approximation "<<this->_approximation->getKind()<<"\n");

            auto e=this->_approximation->compute_error(f,g,V,h,B);
            total_error += e;

            ValidatedVectorFunctionModelDP current_reach_function;
            ValidatedVectorFunctionModelDP current_evolve_function;

            if (this->_approximation->getKind() != DIApproximationKind::PIECEWISE) {
                auto Phi = this->compute_flow_function(f,g,V,D,h,B);
                ARIADNE_LOG(5,"Phi="<<Phi<<"\n");
                assert(Phi.domain()[n].upper()==h);

                current_reach_function=build_reach_function(evolve_function, Phi, t, new_t);
                ARIADNE_LOG(5,"current_reach_function="<<current_reach_function<<"\n");

                current_evolve_function=partial_evaluate(current_reach_function,n,new_t);
                ARIADNE_LOG(5,"current_evolve_function="<<current_evolve_function<<"\n");
            } else {
                InclusionIntegratorPiecewiseApproximation& approx = dynamic_cast<InclusionIntegratorPiecewiseApproximation&>(*this->_approximation);
                auto m=V.size();
                auto e=approx.compute_error(f,g,V,h,B);
                ARIADNE_LOG(6,"approximation error:"<<e<<"\n");
                auto swp=this->_sweeper;
                auto FD=approx.build_flow_domain(D, hlf(h), V);
                ARIADNE_LOG(6,"FD:"<<FD<<"\n");
                auto x0f=ValidatedVectorTaylorFunctionModelDP::projection(FD,range(0,n),swp);

                auto w1 =approx.build_firsthalf_approximating_function(FD, n, m);
                ValidatedVectorTaylorFunctionModelDP& w1f = dynamic_cast<ValidatedVectorTaylorFunctionModelDP&>(w1.reference());
                ARIADNE_LOG(6,"w1f:"<<w1f<<"\n");

                auto phi1=ValidatedVectorTaylorFunctionModelDP(n,FD,swp);
                for (auto i : range(n)) { phi1[i]=phi1[i]+cast_singleton(B[i]); }
                ARIADNE_LOG(6,"phi01:"<<phi1<<"\n");

                for (auto i : range(6)) {
                    phi1=antiderivative(get_time_derivative(phi1,f,g,w1f),n)+x0f;
                }
                ARIADNE_LOG(7,(derivative(phi1,n) - get_time_derivative(phi1,f,g,w1f)).range()<<"\n");
                PositiveFloatDPValue intermediate_t=cast_positive(cast_exact((t+hlf(h)).lower()));

                current_reach_function=build_reach_function(evolve_function, phi1, t, intermediate_t);
                ARIADNE_LOG(5,"current_reach_function="<<current_reach_function<<"\n");

                current_evolve_function=partial_evaluate(current_reach_function,n,intermediate_t);
                ARIADNE_LOG(5,"current_evolve_function="<<current_evolve_function<<"\n");

                auto D2 = cast_exact_box(current_evolve_function.range());
                auto FD2=approx.build_flow_domain(D2, hlf(h), V);
                auto w2 =approx.build_secondhalf_approximating_function(FD2, n, m);
                ValidatedVectorTaylorFunctionModelDP& w2f = dynamic_cast<ValidatedVectorTaylorFunctionModelDP&>(w2.reference());
                ARIADNE_LOG(6,"w2f:"<<w2f<<"\n");

                auto x0f2=ValidatedVectorTaylorFunctionModelDP::projection(FD2,range(0,n),swp);

                auto phi2=ValidatedVectorTaylorFunctionModelDP(n,FD2,swp);
                for (auto i : range(n)) { phi2[i]=phi2[i]+cast_singleton(B[i]); }
                ARIADNE_LOG(6,"phi02:"<<phi2<<"\n");

                for (auto i : range(6)) {
                    phi2=antiderivative(get_time_derivative(phi2,f,g,w2f),n)+x0f2;
                }

                for (auto i : range(n)) {
                    phi2[i].set_error(phi2[i].error()+e);
                }
                ARIADNE_LOG(7,(derivative(phi2,n) - get_time_derivative(phi2,f,g,w2f)).range()<<"\n");

                current_reach_function=build_reach_function(current_evolve_function, phi2, intermediate_t, new_t);
                ARIADNE_LOG(5,"current_reach_function="<<current_reach_function<<"\n");

                current_evolve_function=partial_evaluate(current_reach_function,n,new_t);
                ARIADNE_LOG(5,"current_evolve_function="<<current_evolve_function<<"\n");
            }

            if (i == 0) {
                best_reach_function = current_reach_function;
                best_evolve_function = current_evolve_function;
                best = this->_approximation->getKind();
                for (auto i: range(n)) {
                    best_total_diameter += best_evolve_function.range()[i].width().raw();
                }
            } else {
                FloatDP current_total_diameter(0);
                for (auto i: range(n)) {
                    current_total_diameter += current_evolve_function.range()[i].width().raw();
                }
                if (current_total_diameter < best_total_diameter) {
                    best = this->_approximation->getKind();
                    ARIADNE_LOG(5,"best approximation: " << this->_approximation->getKind() << "\n");
                    best_reach_function = current_reach_function;
                    best_evolve_function = current_evolve_function;
                    best_total_diameter = current_total_diameter;
                }
            }
        }

        ARIADNE_LOG(2,"chosen approximation: " << best << "\n");

        reach_function = best_reach_function;
        evolve_function = best_evolve_function;

        if (step%this->_number_of_steps_between_simplifications==0) {
            this->_reconditioner->simplify(evolve_function);
            ARIADNE_LOG(5,"simplified_evolve_function="<<evolve_function<<"\n");
        }

        evolve_function = this->_reconditioner->expand_errors(evolve_function);

        step+=1;

        t=new_t;
        result.append(reach_function);
    }

    std::cout << "    " << evolve_function.range()[0].upper() << " " << evolve_function.range()[0].upper() - 2*V[0].upper() << " " << total_error << " " << total_error.raw()/step << std::endl;

    return result;
}

ValidatedVectorFunctionModelDP InclusionIntegrator::build_reach_function(
        ValidatedVectorFunctionModelDP evolve_function, ValidatedVectorFunctionModelDP Phi, PositiveFloatDPValue t,
        PositiveFloatDPValue new_t) const {

    // Evolve function is xi(x,a) at s; Flow is phi(x,h,b)
    // Want (x,t,a,b):->phi(xi(x,a),t-s,b))

    SizeType n=evolve_function.result_size();

    SizeType npxE=evolve_function.argument_size()-n;
    SizeType npxP=Phi.argument_size()-(n+1);

    BoxDomainType D=evolve_function.domain()[range(0,n)];
    BoxDomainType PXE=evolve_function.domain()[range(n,n+npxE)];
    BoxDomainType PXP=Phi.domain()[range(n+1,n+1+npxP)];

    auto swp=this->_sweeper;
    auto Tau=IntervalDomainType(t,new_t);
    BoxDomainType DTP = join(D,Tau,PXE,PXP);
    ValidatedVectorTaylorFunctionModelDP xf=ValidatedVectorTaylorFunctionModelDP::projection(DTP,range(0,n),swp);
    ValidatedScalarTaylorFunctionModelDP tf=ValidatedScalarTaylorFunctionModelDP::coordinate(DTP,n,swp);
    ValidatedScalarTaylorFunctionModelDP hf=tf-t;
    ValidatedVectorTaylorFunctionModelDP a0f=ValidatedVectorTaylorFunctionModelDP::projection(DTP,range(n+1,n+1+npxE),swp);
    ValidatedVectorTaylorFunctionModelDP a1f=ValidatedVectorTaylorFunctionModelDP::projection(DTP,range(n+1+npxE,n+1+npxE+npxP),swp);

    ValidatedVectorTaylorFunctionModelDP ef=compose(evolve_function,join(xf,a0f));

    return compose(Phi,join(ef,hf,a1f));
}

//! Computes q(D), where q = f + g * V
UpperBoxType apply(ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, UpperBoxType V, UpperBoxType D) {

    UpperBoxType result = apply(f,D);
    for (auto i : range(g.size())) {
        result = result + apply(g[i],D) * V[i];
    }
    return result;
}

Pair<PositiveFloatDPValue,UpperBoxType> InclusionIntegrator::flow_bounds(ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, UpperBoxType V, BoxDomainType D, PositiveFloatDPApproximation hsug) const {

    //! Compute a bound B for the differential inclusion dot(x) in f(x) + G(x) * V, for x(0) in D for step size h;
    ARIADNE_LOG(5,"D:"<<D);

    PositiveFloatDPValue h=cast_exact(hsug);
    UpperBoxType wD = D + (D-D.midpoint());
    UpperBoxType B = wD + 2*IntervalDomainType(0,h)*apply(f,g,V,D);

    while(not refines(D+h*apply(f,g,V,B),B)) {
        h=hlf(h);
    }
    for(auto i : range(4)) {
        B=D+IntervalDomainType(0,h)*apply(f,g,V,B);
    }

    return std::make_pair(h,B);
}


ValidatedVectorFunctionModelDP InclusionIntegrator::
compute_flow_function(ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, BoxDomainType V, BoxDomainType D,
                      PositiveFloatDPValue h, UpperBoxType B) const {
    auto n=D.size();
    auto m=V.size();
    auto e=_approximation->compute_error(f,g,V,h,B);
    ARIADNE_LOG(6,"approximation error:"<<e<<"\n");
    auto swp=this->_sweeper;
    auto FD=_approximation->build_flow_domain(D, h, V);
    ARIADNE_LOG(6,"FD:"<<FD<<"\n");
    auto x0f=ValidatedVectorTaylorFunctionModelDP::projection(FD,range(0,n),swp);

    auto w =_approximation->build_approximating_function(FD, n, m);
    ValidatedVectorTaylorFunctionModelDP& wf = dynamic_cast<ValidatedVectorTaylorFunctionModelDP&>(w.reference());
    ARIADNE_LOG(6,"wf:"<<wf<<"\n");

    auto phi=ValidatedVectorTaylorFunctionModelDP(n,FD,swp);
    for (auto i : range(n)) { phi[i]=phi[i]+cast_singleton(B[i]); }
    ARIADNE_LOG(6,"phi0:"<<phi<<"\n");

    for (auto i : range(6)) {
        phi=antiderivative(get_time_derivative(phi,f,g,wf),n)+x0f;
    }
    ARIADNE_LOG(7,(derivative(phi,n) - get_time_derivative(phi,f,g,wf)).range()<<"\n");

    for (auto i : range(n)) {
        phi[i].set_error(phi[i].error()+e);
    }

    return phi;
}


ErrorType InclusionIntegratorZeroApproximation::compute_error(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType V, PositiveFloatDPValue h, UpperBoxType const& B) const {
    if (is_identity_matrix(g,B))
        return AdditiveZeroErrorProcessor(f,g,V,h,B).process();
    else
        return ZeroErrorProcessor(f,g,V,h,B).process();
}


ErrorType InclusionIntegratorConstantApproximation::compute_error(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType V, PositiveFloatDPValue h, UpperBoxType const& B) const {
    if (is_identity_matrix(g,B))
        return AdditiveConstantErrorProcessor(f,g,V,h,B).process();
    else
        return ConstantErrorProcessor(f,g,V,h,B).process();
}

ErrorType InclusionIntegratorPiecewiseApproximation::compute_error(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType V, PositiveFloatDPValue h, UpperBoxType const& B) const {
    if (is_identity_matrix(g,B))
        return AdditivePiecewiseErrorProcessor(f,g,V,h,B).process();
    else if (g.size() == 1)
        return SingleInputPiecewiseErrorProcessor(f,g,V,h,B).process();
    else
        return PiecewiseErrorProcessor(f,g,V,h,B).process();
}

ErrorType InclusionIntegratorAffineApproximation::compute_error(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType V, PositiveFloatDPValue h, UpperBoxType const& B) const {
    if (is_identity_matrix(g,B))
        return AdditiveAffineErrorProcessor(f,g,V,h,B).process();
    else if (g.size() == 1)
        return SingleInputAffineErrorProcessor(f,g,V,h,B).process();
    else
        return AffineErrorProcessor(f,g,V,h,B).process();
}


ErrorType InclusionIntegratorSinusoidalApproximation::compute_error(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType V, PositiveFloatDPValue h, UpperBoxType const& B) const {
    if (is_identity_matrix(g,B))
        return AdditiveSinusoidalErrorProcessor(f,g,V,h,B).process();
    else if (g.size() == 1)
        return SingleInputSinusoidalErrorProcessor(f,g,V,h,B).process();
    else
        return SinusoidalErrorProcessor(f,g,V,h,B).process();
}


BoxDomainType InclusionIntegratorZeroApproximation::build_flow_domain(BoxDomainType D, PositiveFloatDPValue h, BoxDomainType V) const {

    auto Ht=IntervalDomainType(-h,+h);

    return join(D,Ht);
}


BoxDomainType InclusionIntegratorConstantApproximation::build_flow_domain(BoxDomainType D, PositiveFloatDPValue h, BoxDomainType V) const {

    auto Ht=IntervalDomainType(-h,+h);
    auto P0=V;

    return join(D,Ht,V);
}


BoxDomainType InclusionIntegratorPiecewiseApproximation::build_flow_domain(BoxDomainType D, PositiveFloatDPValue h, BoxDomainType V) const {

    auto Ht=IntervalDomainType(-h,+h);
    auto P0=V;
    auto P1=V;

    return join(D,Ht,P0,P1);
}


BoxDomainType InclusionIntegratorAffineApproximation::build_flow_domain(BoxDomainType D, PositiveFloatDPValue h, BoxDomainType V) const {

    auto Ht=IntervalDomainType(-h,+h);
    auto P0=V;
    auto P1=V;

    return join(D,Ht,P0,P1);
}


BoxDomainType InclusionIntegratorSinusoidalApproximation::build_flow_domain(BoxDomainType D, PositiveFloatDPValue h, BoxDomainType V) const {

    auto Ht=IntervalDomainType(-h,+h);
    auto P0=V;
    auto P1=V;

    return join(D,Ht,P0,P1);
}

ValidatedVectorFunctionModelType InclusionIntegratorZeroApproximation::build_approximating_function(BoxDomainType FD, SizeType n, SizeType m) const {

    auto swp=this->_sweeper;

    auto zero=ValidatedScalarTaylorFunctionModelDP::zero(FD,swp);

    auto result=ValidatedVectorTaylorFunctionModelDP(m,FD,swp);
    for (auto i : range(m)) { result[i]=zero; }

    return result;
}


ValidatedVectorFunctionModelType InclusionIntegratorConstantApproximation::build_approximating_function(BoxDomainType FD, SizeType n, SizeType m) const {

    auto swp=this->_sweeper;

    auto p0f=ValidatedVectorTaylorFunctionModelDP::projection(FD,range(n+1,n+1+m),swp);

    auto result=ValidatedVectorTaylorFunctionModelDP(m,FD,swp);
    for (auto i : range(m)) { result[i]=p0f[i]; }

    return result;
}


ValidatedVectorFunctionModelType InclusionIntegratorPiecewiseApproximation::build_approximating_function(BoxDomainType FD, SizeType n, SizeType m) const {

    return build_firsthalf_approximating_function(FD,n,m);
}



ValidatedVectorFunctionModelType InclusionIntegratorPiecewiseApproximation::build_firsthalf_approximating_function(BoxDomainType FD, SizeType n, SizeType m) const {

    auto swp = this->_sweeper;

    auto p0f = ValidatedVectorTaylorFunctionModelDP::projection(FD, range(n+1, n+1+m), swp);
    auto p1f = ValidatedVectorTaylorFunctionModelDP::projection(FD, range(n+1+m, n+1+2*m), swp);

    auto one = ValidatedScalarTaylorFunctionModelDP::constant(FD, 1.0_exact, swp);


    auto result = ValidatedVectorTaylorFunctionModelDP(m, FD, swp);
    for (auto i : range(m)) {
        auto Vi=FD[n+1+i].upper();
        result[i] = p0f[i] - (one - p0f[i]*p0f[i]/Vi/Vi)*p1f[i];
    }

    return result;
}


ValidatedVectorFunctionModelType InclusionIntegratorPiecewiseApproximation::build_secondhalf_approximating_function(BoxDomainType FD, SizeType n, SizeType m) const {

    auto swp = this->_sweeper;

    auto p0f = ValidatedVectorTaylorFunctionModelDP::projection(FD, range(n+1, n+1+m), swp);
    auto p1f = ValidatedVectorTaylorFunctionModelDP::projection(FD, range(n+1+m, n+1+2*m), swp);

    auto one = ValidatedScalarTaylorFunctionModelDP::constant(FD, 1.0_exact, swp);

    auto result = ValidatedVectorTaylorFunctionModelDP(m, FD, swp);
    for (auto i : range(m)) {
        auto Vi = FD[n+1+i].upper();
        result[i] = p0f[i] + (one - p0f[i]*p0f[i]/Vi/Vi)*p1f[i];
    }

    return result;
}

ValidatedVectorFunctionModelType InclusionIntegratorAffineApproximation::build_approximating_function(BoxDomainType FD, SizeType n, SizeType m) const {

    auto swp=this->_sweeper;

    auto tf=ValidatedScalarTaylorFunctionModelDP::coordinate(FD,n,swp);
    auto p0f=ValidatedVectorTaylorFunctionModelDP::projection(FD,range(n+1,n+1+m),swp);
    auto p1f=ValidatedVectorTaylorFunctionModelDP::projection(FD,range(n+1+m,n+1+2*m),swp);

    auto one=ValidatedScalarTaylorFunctionModelDP::constant(FD,1.0_exact,swp);
    auto three=ValidatedScalarTaylorFunctionModelDP::constant(FD,3.0_exact,swp);

    auto h = FD[n].upper();

    auto result=ValidatedVectorTaylorFunctionModelDP(m,FD,swp);
    for (auto i : range(m)) {
        auto Vi = FD[n+1+i].upper();
        result[i]=p0f[i]+three*(one-p0f[i]*p0f[i]/Vi/Vi)*p1f[i]*(tf-h/2)/h;
    }

    return result;
}


ValidatedVectorFunctionModelType InclusionIntegratorSinusoidalApproximation::build_approximating_function(BoxDomainType FD, SizeType n, SizeType m) const {

    auto swp=this->_sweeper;

    auto tf=ValidatedScalarTaylorFunctionModelDP::coordinate(FD,n,swp);
    auto p0f=ValidatedVectorTaylorFunctionModelDP::projection(FD,range(n+1,n+1+m),swp);
    auto p1f=ValidatedVectorTaylorFunctionModelDP::projection(FD,range(n+1+m,n+1+2*m),swp);

    auto h = FD[n].upper();
    auto one=ValidatedScalarTaylorFunctionModelDP::constant(FD,1.0_exact,swp);
    auto pgamma=1.1464_dec;
    auto gamma=4.162586_dec;

    auto result=ValidatedVectorTaylorFunctionModelDP(m,FD,swp);
    for (auto i : range(m)) {
        auto Vi = FD[n+1+i].upper();
        result[i]=p0f[i]+(one-p0f[i]*p0f[i]/Vi/Vi)*pgamma*p1f[i]*sin((tf-h/2)*gamma/h);
    }

    return result;
}


LohnerReconditioner::LohnerReconditioner(SweeperDP sweeper, Nat number_of_variables_to_keep)
    : _sweeper(sweeper), _number_of_variables_to_keep(number_of_variables_to_keep) {
    this->verbosity = 0;
}

ValidatedVectorFunctionModelDP LohnerReconditioner::expand_errors(ValidatedVectorFunctionModelDP Phi) const {
    BoxDomainType domain=Phi.domain();
    BoxDomainType errors=cast_exact(cast_exact(Phi.errors())*FloatDPUpperInterval(-1,+1)); // FIXME: Avoid cast;

    ARIADNE_LOG(6,"Uniform errors:"<<errors<<"\n");
    for(SizeType i=0; i!=Phi.result_size(); ++i) { Phi[i].set_error(0); }
    ValidatedVectorFunctionModelDP error_function=ValidatedVectorTaylorFunctionModelDP::identity(errors,this->_sweeper);
    return embed(Phi,errors)+embed(domain,error_function);
}

struct IndexedFloatDPError
{
    SizeType index;
    FloatDPError value;

    IndexedFloatDPError() : index(0), value(FloatDPError()) {}
};

OutputStream& operator<<(OutputStream& os, IndexedFloatDPError const& ifl) {
    return os << "(" << ifl.index << ":" << std::scientific << ifl.value.raw() << std::fixed << ")"; }

struct IndexedFloatDPErrorComparator
{
    inline bool operator() (const IndexedFloatDPError& ifl1, const IndexedFloatDPError& ifl2)
    {
        return (ifl1.value.raw() < ifl2.value.raw());
    }
};

Void LohnerReconditioner::simplify(ValidatedVectorFunctionModelDP& phi) const {
    ARIADNE_LOG(6,"simplifying\n");
    ARIADNE_LOG(6,"phi="<<phi<<"\n");

    auto m=phi.argument_size();
    auto n=phi.result_size();

    ARIADNE_LOG(6,"num.parameters="<<m<<", to keep="<< this->_number_of_variables_to_keep <<"\n");

    ValidatedVectorTaylorFunctionModelDP& tphi = dynamic_cast<ValidatedVectorTaylorFunctionModelDP&>(phi.reference());

    // Compute effect of error terms, but not of original variables;
    Matrix<FloatDPError> C(m,n);
    for (auto i : range(n)) {
        auto p=tphi[i].model().expansion();

        for (auto ac : p) {
            MultiIndex const& a=ac.index();
            FloatDPValue& c=ac.coefficient();
            for (auto j : range(m)) {
                if (a[j]!=0) {
                    C[j][i] += mag(c);
                }
            }
        }
    }

    ARIADNE_LOG(6,"C"<<C<<"\n");

    Array<IndexedFloatDPError> Ce(m);
    for (auto j : range(m)) {
        Ce[j].index = j;
        for (auto i : range(n)) {
            Ce[j].value += C[j][i];
        }
    }
    ARIADNE_LOG(6,"Ce:"<<Ce<<"\n");
    auto SCe=Ce;
    std::sort(SCe.begin(),SCe.end(),IndexedFloatDPErrorComparator());
    ARIADNE_LOG(6,"SortedCe:"<<SCe<<"\n");
    List<SizeType> keep_indices;
    List<SizeType> remove_indices;
    int number_of_variables_to_remove = m - this->_number_of_variables_to_keep;
    ARIADNE_LOG(6, "Number of variables to remove:" << number_of_variables_to_remove<<"\n");

    FloatDPError total_sum_SCe(0);
    for (int j : range(m))
        total_sum_SCe += SCe[j].value;

    FloatDP coeff(1.0/50.0);

    bool skip = false;
    FloatDPError current_sum_SCe(0);
    for (int j : range(m)) {
        current_sum_SCe += SCe[j].value;
        if (!skip && current_sum_SCe.raw() < total_sum_SCe.raw() * coeff) {
            remove_indices.append(SCe[j].index);
        } else {
            keep_indices.append(SCe[j].index);
            skip = true;
        }
    }

    /*
    for (int j : range(m)) {
        if (j < number_of_variables_to_remove) {
            remove_indices.append(SCe[j].index);
        } else {
            keep_indices.append(SCe[j].index);
        }
    }
     */

    ARIADNE_LOG(2,"number of kept parameters: " << keep_indices.size() << "/" << m << "\n");

    ARIADNE_LOG(6,"keep_indices:"<<keep_indices<<"\n");
    ARIADNE_LOG(6,"remove_indices:"<<remove_indices<<"\n");

    for (int i : range(n)) {
        ErrorType error = tphi[i].error();
        for(SizeType k=0; k!=remove_indices.size(); ++k) {
            error += mag(C[remove_indices[k]][i]);
        }
        tphi[i].set_error(error);
    }

    auto old_domain=phi.domain();
    auto new_domain=BoxDomainType(Vector<IntervalDomainType>(keep_indices.size(),[&old_domain,&keep_indices](SizeType j){return old_domain[keep_indices[j]];}));
    auto projection=ValidatedVectorTaylorFunctionModelDP(m,new_domain,this->_sweeper);
    for (auto i : range(new_domain.size())) { projection[keep_indices[i]]=ValidatedScalarTaylorFunctionModelDP::coordinate(new_domain,i,this->_sweeper); }
    for (auto i : range(remove_indices.size())) {
        auto j=remove_indices[i]; auto cj=old_domain[j].midpoint();
        projection[j]=ValidatedScalarTaylorFunctionModelDP::constant(new_domain,cj,this->_sweeper); }
    phi=compose(phi,projection);
}

} // namespace Ariadne;


/*

#include "geometry/zonotope.hpp"

namespace Ariadne {

ValidatedVectorTaylorFunctionModelDP lohner_approximation(ValidatedVectorTaylorFunctionModelDP f) {
    auto n=f.result_size();
    auto models=f.models();
    DoublePrecision pr;
    PositiveFloatDPValue zero(pr);
    Vector<FloatDPValue> b=Vector<FloatDPValue>(n,zero);
    Vector<FloatDPError> e=Vector<FloatDPError>(n,zero);
    Matrix<FloatDPValue> A=Matrix<FloatDPValue>(n,models[0].argument_size(),zero);
    for (auto i : range(n)) {
        b[i]=models[i].value();
        for (auto j : range(models[0].argument_size())) {
            A[i][j]=models[i].gradient_value(j);
        }
        e[i]=models[i].error();
    }
    auto z=Zonotope(b,A,e);
    // print z.error();
    z=orthogonal_approximation(z);

    b=reinterpret_cast<Vector<FloatDPValue>const&>(z.centre());
    A=reinterpret_cast<Matrix<FloatDPValue>const&>(z.generators());
    e=reinterpret_cast<Vector<FloatDPError>const&>(z.error());
    auto p=z.number_of_generators();
    Vector<ValidatedTaylorModelDP> r(n,ValidatedTaylorModelDP(p,f.properties()));
    for (auto i : range(n)) {
        r[i].set_value(b[i]);
        for (auto j : range(p)) {
            r[i].set_gradient(j,A[i][j]);
        }
        r[i].set_error(e[i]);
    }

    return ValidatedVectorTaylorFunctionModelDP(BoxDomainType(n,IntervalDomainType(-1,+1)),r);
}



} // namespace Ariadne;

*/
