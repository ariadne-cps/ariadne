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

struct ScheduledApproximation
{
    SizeType step;
    SharedPointer<InputApproximation> approximation;

    ScheduledApproximation(SizeType step, SharedPointer<InputApproximation> approximation) : step(step), approximation(approximation) {}
};

OutputStream& operator<<(OutputStream& os, ScheduledApproximation const& sa) {
    return os << "(" << sa.step << ":" << sa.approximation->getKind() << ")"; }

struct ScheduledApproximationComparator
{
    inline bool operator() (const ScheduledApproximation& sa1, const ScheduledApproximation& sa2)
    {
        return (sa1.step > sa2.step);
    }
};

FloatDPError get_r(InputApproximationKind approx_kind) {
    switch (approx_kind) {
    case InputApproximationKind::AFFINE:
        return FloatDPError(5.0/3u);
    case InputApproximationKind::SINUSOIDAL:
        return FloatDPError(5.0/4u);
    case InputApproximationKind::PIECEWISE:
        return FloatDPError(1.3645_upper);
    default:
        ARIADNE_FAIL_MSG("A value of 'r' does not exist for kind " << approx_kind << "\n");
    }
}

Box<UpperIntervalType> apply(VectorFunction<ValidatedTag>const& f, const Box<ExactIntervalType>& bx) {
    return apply(f,Box<UpperIntervalType>(bx));
}

FloatDP volume(Vector<ApproximateIntervalType> const& box) {
    FloatDP result = 1.0;
    for (auto i: range(box.size())) {
        result *= box[i].width().raw();
    }
    return result;
}


Boolean inputs_are_additive(Vector<ValidatedVectorFunction> const &g, UpperBoxType const &B) {

    SizeType m = g.size();
    SizeType n = g[0].result_size();

    if (m > n)
        return false;

    for (SizeType j: range(n)) {
        bool foundOne = false;
        for (SizeType i : range(m)) {
            auto eval = g[i][j].evaluate(cast_singleton(ExactBoxType(n,ExactIntervalType(-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()))));
            if (!foundOne) {
                if (definitely(eval == 1.0_exact))
                    foundOne = true;
                else if (definitely(eval != 0.0_exact))
                    return false;
            } else {
                if (possibly(eval != 0.0_exact))
                    return false;
            }
        }
    }
    return true;
}

ValidatedVectorFunction construct_function_affine_in_input(ValidatedVectorFunction const &f,
                                                           Vector<ValidatedVectorFunction> const &g,
                                                           Vector<ValidatedScalarFunction> const &u);

Norms::Norms(FloatDPError const& K,Vector<FloatDPError> const& Kj,FloatDPError const& pK,Vector<FloatDPError> const& pKj,
             FloatDPError const& L,Vector<FloatDPError> const& Lj,FloatDPError const& pL,Vector<FloatDPError> const& pLj,
             FloatDPError const& H,Vector<FloatDPError> const& Hj,FloatDPError const& pH,Vector<FloatDPError> const& pHj,
             FloatDPError const& expLambda,FloatDPError const& expL)
 : K(K), Kj(Kj), pK(pK), pKj(pKj), L(L), Lj(Lj), pL(pL), pLj(pLj), H(H), Hj(Hj), pH(pH), pHj(pHj), expLambda(expLambda), expL(expL) {
    _dimension = Kj.size();
    assert(Kj.size() == _dimension and pKj.size() == _dimension and Lj.size() == _dimension and pKj.size() == _dimension and Hj.size() == _dimension && pHj.size() == _dimension);
}

Tuple<FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,FloatDPError>
Norms::values() const {
    return std::tie(this->K,this->Kj,this->pK,this->pKj,this->L,this->Lj,this->pL,this->pLj,this->H,this->Hj,this->pH,this->pHj,this->expLambda,this->expL);
}


Norms
compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B) {

    auto n = f.result_size();
    auto m = g.size();
    DoublePrecision pr;
    FloatDPError ze(pr);
    FloatDPError K=ze, pK=ze, L=ze, pL=ze, H=ze, pH=ze;
    Vector<FloatDPError> Kj(n), pKj(n), Lj(n), pLj(n), Hj(n), pHj(n);
    FloatDPUpperBound Lambda=ze;

    auto Df=f.differential(cast_singleton(B),2);
    for (auto j : range(n)) {
        auto Df_j=Df[j].expansion();
        FloatDPError K_j=ze, L_j=ze, H_j=ze; FloatDPUpperBound Lambda_j=ze;
        for (auto ac : Df_j) {
            MultiIndex const& a=ac.index();
            FloatDPBounds const& c=ac.coefficient();
            if (a.degree()==0) {
                K_j += mag(c);
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

    Matrix<FloatDPError> pK_matrix(m,n), pL_matrix(m,n), pH_matrix(m,n);

    for (auto i : range(m)) {
        auto Dg_i=g[i].differential(cast_singleton(B),2);
        FloatDPError Vi(abs(V[i]).upper());
        FloatDPError pK_i=ze, pL_i=ze, pH_i=ze;
        for (auto j : range(n)) {
            auto Dg_ij=Dg_i[j].expansion();
            FloatDPError pK_ij=ze, pL_ij=ze, pH_ij=ze;
            for (auto ac : Dg_ij) {
                MultiIndex const& a=ac.index();
                FloatDPBounds const& c=ac.coefficient();
                if (a.degree()==0) {
                    pK_ij += mag(c);
                } else if (a.degree()==1) {
                    pL_ij += mag(c);
                } else {
                    assert(a.degree()==2);
                    pH_ij += mag(c);
                }
            }
            pK_i=max(pK_i,pK_ij); pL_i=max(pL_i,pL_ij); pH_i=max(pH_i,pH_ij);
            pK_matrix[i][j] += pK_ij; pL_matrix[i][j] += pL_ij; pH_matrix[i][j] += pH_ij;
        }

        pK+=Vi*pK_i; pL+=Vi*pL_i; pH+=Vi*pH_i;
    }

    for (auto j : range(n)) {
        pKj[j] = ze; pLj[j] = ze; pHj[j] = ze;
        for (auto i : range(m)) {
            FloatDPError Vi(abs(V[i]).upper());
            pKj[j] += Vi*pK_matrix[i][j]; pLj[j] += Vi*pL_matrix[i][j]; pHj[j] += Vi*pH_matrix[i][j];
        }
    }

    FloatDPError expLambda = (possibly(Lambda>0)) ? FloatDPError(dexp(Lambda*h)) : FloatDPError(1u,pr);
    FloatDPError expL = cast_positive(exp(L*h));

    return Norms(K,Kj,pK,pKj,L,Lj,pL,pLj,H,Hj,pH,pHj,expLambda,expL);
}


Vector<ErrorType> ApproximationErrorProcessor::process(PositiveFloatDPValue const& h, UpperBoxType const& B) const {

    Norms norms = compute_norms(_f,_g,_V,h,B);

    if (inputs_are_additive(_g,B))
        norms.pK=mag(norm(_V));

    return compute_errors(norms,h);
}

ErrorType noparam_error_option1(Norms const& n, PositiveFloatDPValue const& h) {
    return n.pK*n.expLambda*h;
}

ErrorType noparam_error_option2(Norms const& n, PositiveFloatDPValue const& h) {
    return (n.K*2u+n.pK)*h;
}

ErrorType oneparam_additive_error(Norms const& n, PositiveFloatDPValue const& h) {
    return pow(h,2u)*(n.pK*n.L*n.expLambda);
}

ErrorType oneparam_generic_error(Norms const& n, PositiveFloatDPValue const& h) {
    return (pow(h,2u)*(n.pK*n.pL*n.expLambda*2u + n.pL*(n.K+n.pK)/3u + n.pK*n.L/2u)+ pow(h,3u)*n.pK*(n.L*n.pL + n.L*n.L + n.H*(n.K+n.pK))/2u*n.expLambda)/cast_positive(1u-(h*n.L/2u));
}

ErrorType twoparam_additive_error(Norms const& n, PositiveFloatDPValue const& h, FloatDPError const& r) {
    return (n.pK*(n.H*(n.K+r*n.pK)+n.L*n.L)*n.expLambda + (n.K+n.pK)*n.H*n.pK/2u)/cast_positive(1u-h*n.L/2u)*(r+1u)*pow(h,3u)/4u;
}

ErrorType twoparam_singleinput_error(Norms const& n, PositiveFloatDPValue const& h, FloatDPError const& r) {
    return ((r+1u)*n.pK*((n.pH*2u*r+n.H)*(n.K+r*n.pK)+n.L*n.L+(n.L*3u*r+n.pL*r*r*2u)*n.pL)*n.expLambda + (r+1u)/6u*(n.K+n.pK)*((r+1u)*((n.H*n.pK+n.L*n.pL)*3u +(n.pH*n.K+n.L*n.pL)*4u) + (n.pH*n.pK+n.pL*n.pL)*8u*(r*r+1u)))/cast_positive(1u-h*n.L/2u-h*n.pL*r)*pow(h,3u)/4u;
}

ErrorType twoparam_generic_error(Norms const& n, PositiveFloatDPValue const& h, FloatDPError const& r) {
    return ((r*r+1u)*n.pL*n.pK + (r+1u)*h*n.pK*((n.pH*2u*r + n.H)*(n.K+r*n.pK)+n.L*n.L+(n.L*3u*r+n.pL*r*r*2u)*n.pL)*n.expLambda + (r+1u)/6u*h*(n.K+n.pK)*((n.H*n.pK+n.L*n.pL)*3u+(n.pH*n.K+n.L*n.pL)*4u))/cast_positive(1u-h*n.L/2u-h*n.pL*r)*pow(h,2u)/4u;
}

Vector<ErrorType> ZeroErrorProcessor::compute_errors(Norms const& n, PositiveFloatDPValue const& h) const {
    FloatDPError result1 = noparam_error_option1(n,h);
    FloatDPError result2 = noparam_error_option2(n,h);
    return Vector<ErrorType>(n.dimension(),min(result1,result2));
}

Vector<ErrorType> AdditiveConstantErrorProcessor::compute_errors(Norms const& n, PositiveFloatDPValue const& h) const {
    return Vector<ErrorType>(n.dimension(),oneparam_additive_error(n,h));
}

Vector<ErrorType> ConstantErrorProcessor::compute_errors(Norms const& n, PositiveFloatDPValue const& h) const {
    return Vector<ErrorType>(n.dimension(),oneparam_generic_error(n,h));
}

Vector<ErrorType> AdditiveAffineErrorProcessor::compute_errors(Norms const& n, PositiveFloatDPValue const& h) const {
    return Vector<ErrorType>(n.dimension(),twoparam_additive_error(n,h,get_r(_kind)));
}

Vector<ErrorType> SingleInputAffineErrorProcessor::compute_errors(Norms const& n, PositiveFloatDPValue const& h) const {
    return Vector<ErrorType>(n.dimension(),twoparam_singleinput_error(n,h,get_r(_kind)));
}

Vector<ErrorType> AffineErrorProcessor::compute_errors(Norms const& n, PositiveFloatDPValue const& h) const {
    return Vector<ErrorType>(n.dimension(),twoparam_generic_error(n,h,get_r(_kind)));
}

Vector<ErrorType> AdditiveSinusoidalErrorProcessor::compute_errors(Norms const& n, PositiveFloatDPValue const& h) const {
    return Vector<ErrorType>(n.dimension(),twoparam_additive_error(n,h,get_r(_kind)));
}

Vector<ErrorType> SingleInputSinusoidalErrorProcessor::compute_errors(Norms const& n, PositiveFloatDPValue const& h) const {
    return Vector<ErrorType>(n.dimension(),twoparam_singleinput_error(n,h,get_r(_kind)));
}

Vector<ErrorType> SinusoidalErrorProcessor::compute_errors(Norms const& n, PositiveFloatDPValue const& h) const {
    return Vector<ErrorType>(n.dimension(),twoparam_generic_error(n,h,get_r(_kind)));
}

Vector<ErrorType> AdditivePiecewiseErrorProcessor::compute_errors(Norms const& n, PositiveFloatDPValue const& h) const {
    return Vector<ErrorType>(n.dimension(),twoparam_additive_error(n,h,get_r(_kind)));
}

Vector<ErrorType> SingleInputPiecewiseErrorProcessor::compute_errors(Norms const& n, PositiveFloatDPValue const& h) const {
    return Vector<ErrorType>(n.dimension(),twoparam_singleinput_error(n,h,get_r(_kind)));
}

Vector<ErrorType> PiecewiseErrorProcessor::compute_errors(Norms const& n, PositiveFloatDPValue const& h) const {
    return Vector<ErrorType>(n.dimension(),twoparam_generic_error(n,h,get_r(_kind)));
}


ValidatedVectorTaylorFunctionModelDP build_f_plus_Gw(ValidatedVectorTaylorFunctionModelDP phi,
                                                     ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g,
                                                     ValidatedVectorTaylorFunctionModelDP wf) {
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


InclusionIntegrator::InclusionIntegrator(List<SharedPointer<InputApproximation>> approximations, SweeperDP sweeper, StepSize step_size)
    : _approximations(approximations)
    , _sweeper(sweeper)
    , _step_size(step_size)
    , _number_of_steps_between_simplifications(8)
    , _number_of_variables_to_keep(4)
{
    assert(approximations.size()>0);
}

static const SizeType NUMBER_OF_PICARD_ITERATES=6;

List<ValidatedVectorFunctionModelDP> InclusionIntegrator::flow(const List<DottedRealAssignment>& dynamics, const RealVariablesBox& inputs, const RealVariablesBox& initial, ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, BoxDomainType V, BoxDomainType X0, Real tmax) {
    ARIADNE_LOG(1,"\nf:"<<f<<"\ng:"<<g<<"\nV:"<<V<<"\nX0:"<<X0<<"\ntmax:"<<tmax<<"\n");

    // Ensure all arguments have the correct size;
    auto n=X0.size();
    auto m=V.size();
    auto freq=this->_number_of_steps_between_simplifications;
    DoublePrecision pr;
    assert(f.result_size()==n);
    assert(f.argument_size()==n);
    assert(g.size()==m);

    auto number_of_states = n;
    auto number_of_inputs = m;
    auto state_variables = range(0,n);

    PositiveFloatDPValue hsug(this->_step_size);

    ValidatedVectorFunctionModelDP evolve_function = ValidatedVectorTaylorFunctionModelDP::identity(X0,this->_sweeper);
    auto t=PositiveFloatDPValue(0.0);

    Map<InputApproximationKind,SizeType> approximation_global_frequencies, approximation_local_frequencies;
    for (auto appro: _approximations) {
        approximation_global_frequencies[appro->getKind()] = 0;
        approximation_local_frequencies[appro->getKind()] = 0;
    }

    List<ValidatedVectorFunctionModelDP> result;

    auto step = 0;

    List<ScheduledApproximation> schedule;
    Map<SharedPointer<InputApproximation>,Nat> delays;
    for (auto appro: _approximations) {
        schedule.push_back(ScheduledApproximation(SizeType(step),appro));
        delays[appro] = 0;
    }

    while (possibly(t<FloatDPBounds(tmax,pr))) {
        ARIADNE_LOG(2,"step#:"<<step<<", t:"<<t<<", hsug:"<<hsug << "\n");

        List<SharedPointer<InputApproximation>> approximations_to_use;
        while (!schedule.empty()) {
            auto entry = schedule.back();
            if (entry.step == step) {
                approximations_to_use.push_back(entry.approximation);
                schedule.pop_back();
            } else if (entry.step > step) {
                break;
            }
        }

        if(possibly(t+hsug>FloatDPBounds(tmax,pr))) {  //FIXME: Check types for timing;
            hsug=cast_positive(cast_exact((tmax-t).upper()));
        }

        ARIADNE_LOG(3,"n. of parameters="<<evolve_function.argument_size()<<"\n");

        auto D = cast_exact_box(evolve_function.range());
        UpperBoxType B;
        PositiveFloatDPValue h;

        auto n = f.result_size();
        auto m = g.size();
        Vector<ValidatedScalarFunction> v(m);
        for (auto i : range(0,m))
            v[i] = ValidatedScalarFunction::coordinate(n+m,n+i);

        ValidatedVectorFunction fgv = construct_function_affine_in_input(f, g, v);

        std::tie(h,B)=this->flow_bounds(fgv,V,D,hsug);
        ARIADNE_LOG(2,"flow bounds = "<<B<<" (using h = " << h << ")\n");

        PositiveFloatDPValue new_t=cast_positive(cast_exact((t+h).lower()));

        ValidatedVectorFunctionModelDP reach_function;
        ValidatedVectorFunctionModelDP best_reach_function, best_evolve_function;
        SharedPointer<InputApproximation> best;
        FloatDP best_volume(0);

        ARIADNE_LOG(3,"n. of approximations to use="<<approximations_to_use.size()<<"\n");

        SizeType i = 0;
        for (auto i : range(approximations_to_use.size())) {
            this->_approximation = approximations_to_use.at(i);
            ARIADNE_LOG(4,"checking approximation "<<this->_approximation->getKind()<<"\n");

            ValidatedVectorFunctionModelDP current_reach_function;
            ValidatedVectorFunctionModelDP current_evolve_function;

            if (this->_approximation->getKind() != InputApproximationKind::PIECEWISE) {

                auto Phi = this->compute_flow_function(dynamics,inputs,initial,f,g,V,D,h,B);

                ARIADNE_LOG(5,"Phi="<<Phi<<"\n");
                assert(Phi.domain()[Phi.argument_size()-1].upper()==h);

                current_reach_function=build_reach_function(evolve_function, Phi, t, new_t);
                ARIADNE_LOG(5,"current_reach_function="<<current_reach_function<<"\n");

                current_evolve_function=partial_evaluate(current_reach_function,current_reach_function.argument_size()-1,new_t);
                ARIADNE_LOG(5,"current_evolve_function="<<current_evolve_function<<"\n");
            } else {
                PiecewiseInputApproximation& approx = dynamic_cast<PiecewiseInputApproximation&>(*this->_approximation);

                auto n=D.size();
                auto m=V.size();
                auto number_of_states = n;
                auto number_of_inputs = m;
                auto state_variables = range(0,n);
                auto e=approx.compute_errors(h,B);
                ARIADNE_LOG(6,"approximation errors:"<<e<<"\n");
                auto swp=this->_sweeper;
                auto FD1 = approx.build_flow_domain(D,V,hlf(h));
                auto w1 = approx.build_firsthalf_approximating_function(FD1, number_of_states, number_of_inputs);
                ARIADNE_LOG(6,"FD1:"<<FD1<<"\n");
                ARIADNE_LOG(6,"w1:"<<w1<<"\n");

                auto fgw1 = construct_function_affine_in_input(f, g, w1);

                ARIADNE_LOG(2,"fgw:" << fgw1 << "\n");

                auto x0f1=ValidatedVectorTaylorFunctionModelDP::projection(FD1,state_variables,swp);
                auto af1=ValidatedVectorTaylorFunctionModelDP::projection(FD1,range(n,fgw1.argument_size()),swp);

                auto phi1=ValidatedVectorTaylorFunctionModelDP(number_of_states,FD1,swp);
                phi1=phi1+cast_singleton(B);

                for (auto i : range(NUMBER_OF_PICARD_ITERATES)) {
                    auto f_of_phi1 = compose(fgw1,join(phi1,af1));
                    phi1=antiderivative(f_of_phi1,f_of_phi1.argument_size()-1)+x0f1;
                }
                PositiveFloatDPValue intermediate_t=cast_positive(cast_exact((t+hlf(h)).lower()));

                current_reach_function=build_reach_function(evolve_function, phi1, t, intermediate_t);
                ARIADNE_LOG(5,"current_reach_function="<<current_reach_function<<"\n");

                current_evolve_function=partial_evaluate(current_reach_function,current_reach_function.argument_size()-1,intermediate_t);
                ARIADNE_LOG(5,"current_evolve_function="<<current_evolve_function<<"\n");

                auto D2 = cast_exact_box(current_evolve_function.range());

                auto FD2=approx.build_flow_domain(D2,V,hlf(h));
                auto w2 =approx.build_secondhalf_approximating_function(FD2, number_of_states, number_of_inputs);
                ARIADNE_LOG(6,"w2:"<<w2<<"\n");

                auto x0f2=ValidatedVectorTaylorFunctionModelDP::projection(FD2,state_variables,swp);

                auto fgw2 = construct_function_affine_in_input(f, g, w2);

                ARIADNE_LOG(2,"fgw:" << fgw1 << "\n");

                auto x0f=ValidatedVectorTaylorFunctionModelDP::projection(FD2,state_variables,swp);
                auto af2=ValidatedVectorTaylorFunctionModelDP::projection(FD2,range(n,fgw2.argument_size()),swp);

                auto phi2=ValidatedVectorTaylorFunctionModelDP(number_of_states,FD2,swp);
                phi2=phi2+cast_singleton(B);

                for (auto i : range(NUMBER_OF_PICARD_ITERATES)) {
                    auto f_of_phi2 = compose(fgw2,join(phi2,af2));
                    phi2=antiderivative(f_of_phi2,f_of_phi2.argument_size()-1)+x0f2;
                }

                for (auto i : state_variables) {
                    phi2[i].add_error(e[i]);
                }

                current_reach_function=build_secondhalf_piecewise_reach_function(current_evolve_function, phi2, m, intermediate_t, new_t);
                ARIADNE_LOG(5,"current_reach_function="<<current_reach_function<<"\n");

                current_evolve_function=partial_evaluate(current_reach_function,current_reach_function.argument_size()-1,new_t);
                ARIADNE_LOG(5,"current_evolve_function="<<current_evolve_function<<"\n");
            }

            if (i == 0) {
                best_reach_function = current_reach_function;
                best_evolve_function = current_evolve_function;
                best = this->_approximation;
                best_volume = volume(best_evolve_function.range());
            } else {
                FloatDP current_volume = volume(current_evolve_function.range());
                if (current_volume < best_volume) {
                    best = this->_approximation;
                    ARIADNE_LOG(5,"best approximation: " << best->getKind() << "\n");
                    best_reach_function = current_reach_function;
                    best_evolve_function = current_evolve_function;
                    best_volume = current_volume;
                }
            }
        }

        if (approximations_to_use.size() > 1)
            ARIADNE_LOG(3,"chosen approximation: " << best->getKind() << "\n");

        for (auto appro : approximations_to_use) {
            if (best->getKind() == appro->getKind())
                delays[appro] = 0;
            else
                delays[appro]++;

            Nat offset = 1<<delays[appro];
            schedule.push_back(ScheduledApproximation(step+offset,appro));
        }
        std::sort(schedule.begin(),schedule.end(),ScheduledApproximationComparator());

        ARIADNE_LOG(3,"updated schedule: " << schedule << "\n");

        approximation_global_frequencies[best->getKind()] += 1;
        approximation_local_frequencies[best->getKind()] += 1;

        reach_function = best_reach_function;
        evolve_function = best_evolve_function;

        if (step%freq==freq-1) {

            double base = 0;
            double rho = 6.0;
            for (auto appro: approximation_local_frequencies) {
                SizeType ppi;
                switch (appro.first) {
                    case InputApproximationKind::ZERO:
                        ppi = 0;
                        break;
                    case InputApproximationKind::CONSTANT:
                        ppi = 1;
                        break;
                    default:
                        ppi = 2;
                }
                double partial = n + rho*(n+2*m) + (freq-1)*m*(2 - ppi);
                base += partial*appro.second/freq;
            }
            LohnerReconditioner& lreconditioner = dynamic_cast<LohnerReconditioner&>(*this->_reconditioner);

            Nat num_variables_to_keep(base);
            ARIADNE_LOG(4,"simplifying to "<<num_variables_to_keep<<" variables\n");
            lreconditioner.set_number_of_variables_to_keep(num_variables_to_keep);
            this->_reconditioner->simplify(evolve_function);
            ARIADNE_LOG(5,"simplified_evolve_function="<<evolve_function<<"\n");
            for (auto appro: _approximations) {
                approximation_local_frequencies[appro->getKind()] = 0;
            }
        }

        evolve_function = this->_reconditioner->expand_errors(evolve_function);

        ARIADNE_LOG(2,"evolve bounds="<<evolve_function.range()<<"\n");

        step+=1;

        t=new_t;
        result.append(reach_function);

    }

    ARIADNE_LOG(1,"frequencies="<<approximation_global_frequencies<<"\n");

    return result;
}

ValidatedVectorFunctionModelDP InclusionIntegrator::build_secondhalf_piecewise_reach_function(
        ValidatedVectorFunctionModelDP evolve_function, ValidatedVectorFunctionModelDP Phi, SizeType m, PositiveFloatDPValue t,
        PositiveFloatDPValue new_t) const {

    // Evolve function is e(x,a,2*m) at s; Flow is phi(x,h,b,2*m)
    // Want (x,a,b,2*m,t):->phi(e(x,a,2*m),b,2*m,t-s))

    SizeType n=evolve_function.result_size();

    SizeType a=evolve_function.argument_size()-n-2*m;
    SizeType b=Phi.argument_size()-(n+1)-2*m;

    BoxDomainType X=evolve_function.domain()[range(0,n)];
    BoxDomainType PA=evolve_function.domain()[range(n,n+a)];
    BoxDomainType PB=Phi.domain()[range(n,n+b)];
    BoxDomainType PM=Phi.domain()[range(n+b,n+b+2*m)];

    auto swp=this->_sweeper;
    auto Tau=IntervalDomainType(t,new_t);
    BoxDomainType XPT = join(X,PA,PB,PM,Tau);
    ValidatedVectorTaylorFunctionModelDP xf=ValidatedVectorTaylorFunctionModelDP::projection(XPT,range(0,n),swp);
    ValidatedVectorTaylorFunctionModelDP af=ValidatedVectorTaylorFunctionModelDP::projection(XPT,range(n,n+a),swp);
    ValidatedVectorTaylorFunctionModelDP bf=ValidatedVectorTaylorFunctionModelDP::projection(XPT,range(n+a,n+a+b),swp);
    ValidatedVectorTaylorFunctionModelDP mf=ValidatedVectorTaylorFunctionModelDP::projection(XPT,range(n+a+b,n+a+b+2*m),swp);
    ValidatedScalarTaylorFunctionModelDP tf=ValidatedScalarTaylorFunctionModelDP::coordinate(XPT,n+a+b+2*m,swp);
    ValidatedScalarTaylorFunctionModelDP hf=tf-t;

    ValidatedVectorTaylorFunctionModelDP ef=compose(evolve_function,join(xf,af,mf));

    return compose(Phi,join(ef,bf,mf,hf));
}

ValidatedVectorFunctionModelDP InclusionIntegrator::build_reach_function(
        ValidatedVectorFunctionModelDP evolve_function, ValidatedVectorFunctionModelDP Phi, PositiveFloatDPValue t,
        PositiveFloatDPValue new_t) const {

    // Evolve function is e(x,a) at s; flow is phi(x,b,h)
    // Want (x,a,b,t):->phi(e(x,a),b,t-s))

    SizeType n=evolve_function.result_size();

    SizeType a=evolve_function.argument_size()-n;
    SizeType b=Phi.argument_size()-(n+1);

    BoxDomainType X=evolve_function.domain()[range(0,n)];
    BoxDomainType PA=evolve_function.domain()[range(n,n+a)];
    BoxDomainType PB=Phi.domain()[range(n,n+b)];

    auto swp=this->_sweeper;
    auto Tau=IntervalDomainType(t,new_t);
    BoxDomainType XPT = join(X,PA,PB,Tau);
    ValidatedVectorTaylorFunctionModelDP xf=ValidatedVectorTaylorFunctionModelDP::projection(XPT,range(0,n),swp);
    ValidatedVectorTaylorFunctionModelDP af=ValidatedVectorTaylorFunctionModelDP::projection(XPT,range(n,n+a),swp);
    ValidatedVectorTaylorFunctionModelDP bf=ValidatedVectorTaylorFunctionModelDP::projection(XPT,range(n+a,n+a+b),swp);
    ValidatedScalarTaylorFunctionModelDP tf=ValidatedScalarTaylorFunctionModelDP::coordinate(XPT,n+a+b,swp);
    ValidatedScalarTaylorFunctionModelDP hf=tf-t;

    ValidatedVectorTaylorFunctionModelDP ef=compose(evolve_function,join(xf,af));

    return compose(Phi,join(ef,bf,hf));
}


ValidatedVectorFunction construct_function_affine_in_input(ValidatedVectorFunction const &f,
                                                           Vector<ValidatedVectorFunction> const &g,
                                                           Vector<ValidatedScalarFunction> const &u) {

    auto n = f.result_size();
    auto m = g.size();
    auto p = u[0].argument_size();
    auto infinity_box = cast_singleton(ExactBoxType(p,ExactIntervalType(-std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity())));

    auto coordinates = ValidatedVectorFunction::coordinates(p);

    auto extension = ValidatedVectorFunction::zeros(n,p);
    for (auto i : range(0,n)) {
        extension.set(i,coordinates[i]);
    }

    auto fext = compose(f,extension);
    Vector<ValidatedVectorFunction> gext(m);
    for (Nat j : range(0,m)) {
        gext[j] = compose(g[j],extension);
    }

    ValidatedVectorFunction result = ValidatedVectorFunction::zeros(n,p);
    for (Nat i : range(0,n)) {
        result[i] = fext[i];
        for (Nat j : range(0,m)) {
            auto eval = gext[j][i].evaluate(infinity_box);
            if (definitely(eval == 1.0_exact)) {
                result[i] = result[i] + u[j];
            } else if (definitely(eval == 0.0_exact)) {
            } else {
                result[i] = result[i] + gext[j][i]*u[j];
            }
        }
    }

    return result;
}

ValidatedVectorFunction construct_f_plus_gw_squared(ValidatedVectorFunction const &f, Vector<ValidatedVectorFunction> const &g, Vector<ValidatedScalarFunction> const& w) {

    auto n = f.result_size();
    auto m = g.size();
    auto p = w[0].argument_size();

    auto coordinates = ValidatedVectorFunction::coordinates(p);

    auto extension = ValidatedVectorFunction::zeros(n,p);
    for (auto i : range(0,n)) {
        extension.set(i,coordinates[i]);
    }

    auto fext = compose(f,extension);
    Vector<ValidatedVectorFunction> gext(m);
    for (Nat j : range(0,m)) {
        gext[j] = compose(g[j],extension);
    }

    ValidatedVectorFunction result = ValidatedVectorFunction::zeros(p,p);
    for (Nat i : range(0,n)) {
        result[i] = fext[i];
        for (Nat j : range(0,m)) {
            result[i] = result[i] + gext[j][i]*w[j];
        }
    }
    result[p-1] = ValidatedScalarFunction::constant(p,1_z);
    return result;
}

Pair<PositiveFloatDPValue,UpperBoxType> InclusionIntegrator::flow_bounds(ValidatedVectorFunction f, BoxDomainType V, BoxDomainType D, PositiveFloatDPApproximation hsug) const {

    //! Compute a bound B for the differential inclusion dot(x) in f(x) + G(x) * V, for x(0) in D for step size h;
    ARIADNE_LOG(5,"D:"<<D);

    PositiveFloatDPValue h=cast_exact(hsug);
    UpperBoxType wD = D + (D-D.midpoint());

    ExactBoxType DV = join(D,V);

    UpperBoxType B = wD + 2*IntervalDomainType(0,h)*apply(f,DV);

    UpperBoxType BV = join(B,UpperBoxType(V));

    while(not refines(D+IntervalDomainType(0,h)*apply(f,BV),B)) {
        h=hlf(h);
    }

    for(auto i : range(4)) {
        B=D+IntervalDomainType(0,h)*apply(f,BV);
        BV = join(B,UpperBoxType(V));
    }

    return std::make_pair(h,B);
}


ValidatedVectorFunctionModelDP InclusionIntegrator::
compute_flow_function(const List<DottedRealAssignment>& dynamics, const RealVariablesBox& inputs, const RealVariablesBox& initial, ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, BoxDomainType V, BoxDomainType D,
                      PositiveFloatDPValue h, UpperBoxType B) const {
    auto n=D.size();
    auto m=V.size();
    auto number_of_states = n;
    auto number_of_inputs = m;
    auto state_variables = range(0,n);
    auto e=_approximation->compute_errors(h,B);
    ARIADNE_LOG(6,"approximation errors:"<<e<<"\n");
    auto swp=this->_sweeper;
    auto DVh =_approximation->build_flow_domain(D,V,h);
    auto w = _approximation->build_w_functions(DVh, number_of_states, number_of_inputs);
    ARIADNE_LOG(6,"DVh:"<<DVh<<"\n");
    ARIADNE_LOG(6,"w:"<<w<<"\n");

    auto fgw = construct_function_affine_in_input(f, g, w);

    ARIADNE_LOG(2,"fgw:" << fgw << "\n");

    auto x0f=ValidatedVectorTaylorFunctionModelDP::projection(DVh,state_variables,swp);
    auto af=ValidatedVectorTaylorFunctionModelDP::projection(DVh,range(n,fgw.argument_size()),swp);

    auto picardPhi=ValidatedVectorTaylorFunctionModelDP(number_of_states,DVh,swp);
    picardPhi=picardPhi+cast_singleton(B);

    for (auto i : range(NUMBER_OF_PICARD_ITERATES)) {
        auto f_of_phi = compose(fgw,join(picardPhi,af));
        picardPhi=antiderivative(f_of_phi,f_of_phi.argument_size()-1)+x0f;
    }

    for (auto i : state_variables) {
        picardPhi[i].add_error(e[i]);
    }


/*    TaylorSeriesIntegrator integrator(MaximumError(1e-4),SweepThreshold(1e-8),LipschitzConstant(0.5));

    auto fgws = construct_f_plus_gw_squared(f,g,w);
    auto BVh =_approximation->build_flow_domain(cast_exact_box(B), V, h);

    auto squaredSeriesPhi = integrator.flow_step(fgws,DVh,h,BVh);

    ValidatedVectorTaylorFunctionModelDP& tsquaredSeriesPhi = dynamic_cast<ValidatedVectorTaylorFunctionModelDP&>(squaredSeriesPhi.reference());
    auto seriesPhi=ValidatedVectorTaylorFunctionModelDP(n,squaredSeriesPhi.domain(),swp);
    for (auto i : state_variables) {
        seriesPhi[i] = tsquaredSeriesPhi[i];
        seriesPhi[i].add_error(e);
    }

    if (volume(picardPhi.range()) < volume(seriesPhi.range())) {
        ARIADNE_LOG(2,"Picard flow function chosen\n");
        return picardPhi;

    } else {
        ARIADNE_LOG(2,"Series flow function chosen\n");
        return seriesPhi;
    }*/

    return picardPhi;
}

Vector<ErrorType> InputApproximation::compute_errors(PositiveFloatDPValue h, UpperBoxType const& B) const {
    if (inputs_are_additive(_g, B))
        return _additive_processor->process(h,B);
    else if (_g.size() == 1)
        return _single_input_processor->process(h,B);
    else
        return _generic_processor->process(h,B);
}


BoxDomainType InputApproximation::build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h) const {
    auto result = D;

    for (Nat i : range(this->_num_params_per_input))
        result = product(result,V);

    return product(result,IntervalDomainType(-h,+h));
}

Vector<ValidatedScalarFunction> ZeroInputApproximation::build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const {
    auto result = Vector<ValidatedScalarFunction>(m);

    for (auto i : range(0,m))
        result[i] = ValidatedScalarFunction::zero(n+1);

    return result;
}


Vector<ValidatedScalarFunction> ConstantInputApproximation::build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const {
    auto result = Vector<ValidatedScalarFunction>(m);

    for (auto i : range(0,m))
        result[i] = ValidatedScalarFunction::coordinate(n+m+1,n+i);

    return result;
}


Vector<ValidatedScalarFunction> AffineInputApproximation::build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const {
    auto result = Vector<ValidatedScalarFunction>(m);

    auto zero = ValidatedScalarFunction::zero(n+2*m+1);
    auto one = ValidatedScalarFunction::constant(n+2*m+1,1_z);
    auto three = ValidatedScalarFunction::constant(n+2*m+1,3_z);
    auto t = ValidatedScalarFunction::coordinate(n+2*m+1,n+2*m);

    auto h = ValidatedScalarFunction::constant(n+2*m+1,ExactNumber(DVh[n+2*m].upper()));

    for (auto i : range(m)) {
        auto Vi = ExactNumber(DVh[n+i].upper());
        auto p0 = ValidatedScalarFunction::coordinate(n+2*m+1,n+i);
        auto p1 = ValidatedScalarFunction::coordinate(n+2*m+1,n+m+i);
        result[i] = (definitely (DVh[n+i].upper() == 0.0_exact) ? zero : p0+three*(one-p0*p0/Vi/Vi)*p1*(t-h/2)/h);
    }

    return result;
}


Vector<ValidatedScalarFunction> SinusoidalInputApproximation::build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const {

    auto result = Vector<ValidatedScalarFunction>(m);

    auto zero = ValidatedScalarFunction::zero(n+2*m+1);
    auto one = ValidatedScalarFunction::constant(n+2*m+1,1_z);
    auto three = ValidatedScalarFunction::constant(n+2*m+1,3_z);
    auto t = ValidatedScalarFunction::coordinate(n+2*m+1,n+2*m);

    auto h = ValidatedScalarFunction::constant(n+2*m+1,ExactNumber(DVh[n+2*m].upper()));
    auto pgamma = ValidatedScalarFunction::constant(n+2*m+1,1.1464_dec);
    auto gamma = ValidatedScalarFunction::constant(n+2*m+1,4.162586_dec);

    for (auto i : range(m)) {
        auto Vi = ExactNumber(DVh[n+i].upper());
        auto p0 = ValidatedScalarFunction::coordinate(n+2*m+1,n+i);
        auto p1 = ValidatedScalarFunction::coordinate(n+2*m+1,n+m+i);
        result[i] = (definitely (DVh[n+i].upper() == 0.0_exact) ? zero : p0+(one-p0*p0/Vi/Vi)*pgamma*p1*sin((t-h/2)*gamma/h));
    }

    return result;
}


Vector<ValidatedScalarFunction> PiecewiseInputApproximation::build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const {

    return build_firsthalf_approximating_function(DVh,n,m);
}


Vector<ValidatedScalarFunction> PiecewiseInputApproximation::build_firsthalf_approximating_function(BoxDomainType DVh, SizeType n, SizeType m) const {
    auto result = Vector<ValidatedScalarFunction>(m);

    auto zero = ValidatedScalarFunction::zero(n+2*m+1);
    auto one = ValidatedScalarFunction::constant(n+2*m+1,1_z);

    for (auto i : range(m)) {
        auto Vi = ExactNumber(DVh[n+i].upper());
        auto p0 = ValidatedScalarFunction::coordinate(n+2*m+1,n+i);
        auto p1 = ValidatedScalarFunction::coordinate(n+2*m+1,n+m+i);
        result[i] = (definitely (DVh[n+i].upper() == 0.0_exact) ? zero : p0-(one-p0*p0/Vi/Vi)*p1);
    }
    
    return result;
}


Vector<ValidatedScalarFunction> PiecewiseInputApproximation::build_secondhalf_approximating_function(BoxDomainType DVh, SizeType n, SizeType m) const {
    auto result = Vector<ValidatedScalarFunction>(m);

    auto zero = ValidatedScalarFunction::zero(n+2*m+1);
    auto one = ValidatedScalarFunction::constant(n+2*m+1,1_z);

    for (auto i : range(m)) {
        auto Vi = ExactNumber(DVh[n+i].upper());
        auto p0 = ValidatedScalarFunction::coordinate(n+2*m+1,n+i);
        auto p1 = ValidatedScalarFunction::coordinate(n+2*m+1,n+m+i);
        result[i] = (definitely (DVh[n+i].upper() == 0.0_exact) ? zero : p0+(one-p0*p0/Vi/Vi)*p1);
    }
    
    return result;
}

LohnerReconditioner::LohnerReconditioner(SweeperDP sweeper, Nat number_of_variables_to_keep)
    : _sweeper(sweeper), _number_of_variables_to_keep(number_of_variables_to_keep) {
    this->verbosity = 0;
}

ValidatedVectorFunctionModelDP LohnerReconditioner::expand_errors(ValidatedVectorFunctionModelDP f) const {
    BoxDomainType domain=f.domain();
    BoxDomainType errors=cast_exact(cast_exact(f.errors())*FloatDPUpperInterval(-1,+1)); // FIXME: Avoid cast;

    ARIADNE_LOG(6,"Uniform errors:"<<errors<<"\n");
    for(SizeType i=0; i!=f.result_size(); ++i) { f[i].set_error(0); }
    ValidatedVectorFunctionModelDP error_function=ValidatedVectorTaylorFunctionModelDP::identity(errors,this->_sweeper);
    return embed(f,errors)+embed(domain,error_function);
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

Void LohnerReconditioner::simplify(ValidatedVectorFunctionModelDP& f) const {
    ARIADNE_LOG(6,"simplifying\n");
    ARIADNE_LOG(6,"f="<<f<<"\n");

    auto m=f.argument_size();
    auto n=f.result_size();

    ARIADNE_LOG(6,"num.parameters="<<m<<", to keep="<< this->_number_of_variables_to_keep <<"\n");

    ValidatedVectorTaylorFunctionModelDP& tf = dynamic_cast<ValidatedVectorTaylorFunctionModelDP&>(f.reference());

    // Compute effect of error terms, but not of original variables;
    Matrix<FloatDPError> C(m,n);
    for (auto i : range(n)) {
        auto p=tf[i].model().expansion();

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

    if (number_of_variables_to_remove <= 0)
        return;

    /*
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
    */

    for (int j : range(number_of_variables_to_remove)) {
        remove_indices.append(SCe[j].index);
    }

    for (int j : range(number_of_variables_to_remove,m)) {
        keep_indices.append(SCe[j].index);
    }

    ARIADNE_LOG(2,"number of kept parameters: " << keep_indices.size() << "/" << m << "\n");

    ARIADNE_LOG(6,"keep_indices:"<<keep_indices<<"\n");
    ARIADNE_LOG(6,"remove_indices:"<<remove_indices<<"\n");

    for (int i : range(n)) {
        ErrorType error = tf[i].error();
        for(SizeType k=0; k!=remove_indices.size(); ++k) {
            error += mag(C[remove_indices[k]][i]);
        }
        tf[i].set_error(error);
    }

    auto old_domain=f.domain();
    auto new_domain=BoxDomainType(Vector<IntervalDomainType>(keep_indices.size(),[&old_domain,&keep_indices](SizeType j){return old_domain[keep_indices[j]];}));
    auto projection=ValidatedVectorTaylorFunctionModelDP(m,new_domain,this->_sweeper);
    for (auto i : range(new_domain.size())) { projection[keep_indices[i]]=ValidatedScalarTaylorFunctionModelDP::coordinate(new_domain,i,this->_sweeper); }
    for (auto i : range(remove_indices.size())) {
        auto j=remove_indices[i]; auto cj=old_domain[j].midpoint();
        projection[j]=ValidatedScalarTaylorFunctionModelDP::constant(new_domain,cj,this->_sweeper); }
    f=compose(f,projection);
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
