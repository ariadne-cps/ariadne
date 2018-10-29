/***************************************************************************
 *            differential_inclusion.cpp
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

#include "differential_inclusion.hpp"
#include "../function/taylor_function.hpp"
#include "../solvers/integrator.hpp"
#include "../solvers/bounder.hpp"
#include "../algebra/expansion.inl.hpp"

namespace Ariadne {

BoxDomainType initial_ranges_to_box(RealVariablesBox const& var_ranges) {
    auto vars = var_ranges.variables();
    List<IntervalDomainType> result;
    for (auto v : vars) {
        result.push_back(cast_exact(widen(IntervalDomainType(var_ranges[v].lower().get_d(),var_ranges[v].upper().get_d()))));
    }
    return Vector<IntervalDomainType>(result);
}

#define ARIADNE_LOG_PRINT(level, expr) { ARIADNE_LOG(level,#expr << "=" << (expr) << "\n"); }

struct ScheduledApproximator
{
    Nat step;
    InputApproximator approximator;

    ScheduledApproximator(Nat s, InputApproximator a) : step(s), approximator(a) {}
};

inline OutputStream& operator<<(OutputStream& os, ScheduledApproximator const& sa) {
    return os << "(" << sa.step << ":" << sa.approximator.kind() << ")"; }

struct ScheduledApproximatorComparator {
    inline bool operator() (const ScheduledApproximator& sa1, const ScheduledApproximator& sa2) {
        return (sa1.step > sa2.step); }
};

inline char activity_symbol(SizeType step) {
    switch (step % 4) {
    case 0: return '\\';
    case 1: return '|';
    case 2: return '/';
    default: return '-';
    }
}

inline Box<UpperIntervalType> apply(VectorMultivariateFunction<ValidatedTag>const& f, const Box<ExactIntervalType>& bx) {
    return apply(f,Box<UpperIntervalType>(bx));
}

inline Map<InputApproximationKind,FloatDP> convert_to_percentages(const Map<InputApproximationKind,SizeType>& approximation_global_frequencies) {

    SizeType total_steps(0);
    for (auto entry: approximation_global_frequencies) {
        total_steps += entry.second;
    }

    Map<InputApproximationKind,FloatDP> result;
    for (auto entry: approximation_global_frequencies) {
        result[entry.first] = 1.0/total_steps*entry.second;
    }

    return result;
}

FloatDP volume(Vector<ApproximateIntervalType> const& box) {
    FloatDP result = 1.0;
    for (auto i: range(box.size())) {
        result *= box[i].width().raw();
    }
    return result;
}


C1Norms::C1Norms(FloatDPError const& K_,Vector<FloatDPError> const& Kj_,FloatDPError const& pK_,Vector<FloatDPError> const& pKj_,
             FloatDPError const& L_,Vector<FloatDPError> const& Lj_,FloatDPError const& pL_,Vector<FloatDPError> const& pLj_,
             FloatDPError const& H_,Vector<FloatDPError> const& Hj_,FloatDPError const& pH_,Vector<FloatDPError> const& pHj_,
             FloatDPError const& expLambda_,FloatDPError const& expL_)
 : K(K_), Kj(Kj_), pK(pK_), pKj(pKj_), L(L_), Lj(Lj_), pL(pL_), pLj(pLj_), H(H_), Hj(Hj_), pH(pH_), pHj(pHj_), expLambda(expLambda_), expL(expL_) {
    _dimension = Kj.size();
    assert(Kj.size() == _dimension and pKj.size() == _dimension and Lj.size() == _dimension and pKj.size() == _dimension and Hj.size() == _dimension && pHj.size() == _dimension);
}

Tuple<FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,FloatDPError>
C1Norms::values() const {
    return std::tie(this->K,this->Kj,this->pK,this->pKj,this->L,this->Lj,this->pL,this->pLj,this->H,this->Hj,this->pH,this->pHj,this->expLambda,this->expL);
}


C1Norms
compute_norms(InclusionVectorField const& ivf, PositiveFloatDPValue const& h, UpperBoxType const& B) {

    auto n = ivf.dimension();
    auto m = ivf.number_of_inputs();
    auto V = ivf.inputs();
    DoublePrecision pr;
    FloatDPError ze(pr);
    FloatDPError K=ze, pK=ze, L=ze, pL=ze, H=ze, pH=ze;
    Vector<FloatDPError> Kj(n), pKj(n), Lj(n), pLj(n), Hj(n), pHj(n);
    FloatDPUpperBound Lambda=ze;

    auto Df=ivf.noise_independent_component().differential(cast_singleton(B),2);
    for (auto j : range(n)) {
        auto Df_j=Df[j].expansion();
        FloatDPError K_j=ze, L_j=ze, H_j=ze; FloatDPUpperBound Lambda_j=ze;
        for (auto ac : Df_j) {
            UniformReference<MultiIndex> a=ac.index();
            UniformReference<FloatDPBounds> c=ac.coefficient();
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

    auto input_derivatives = ivf.input_derivatives();

    for (auto i : range(m)) {
        auto Dg_i=input_derivatives[i].differential(cast_singleton(B),2);
        FloatDPError Vi(abs(V[i]).upper());
        FloatDPError pK_i=ze, pL_i=ze, pH_i=ze;
        for (auto j : range(n)) {
            auto Dg_ij=Dg_i[j].expansion();
            FloatDPError pK_ij=ze, pL_ij=ze, pH_ij=ze;
            for (auto ac : Dg_ij) {
                UniformReference<MultiIndex> a=ac.index();
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

    return C1Norms(K,Kj,pK,pKj,L,Lj,pL,pLj,H,Hj,pH,pHj,expLambda,expL);
}

ErrorType zeroparam_worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h) {
    return min(n.pK*n.expLambda*h,
               (n.K*2u+n.pK)*h); }
ErrorType zeroparam_component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j) {
    return min(n.pK*n.expL*h,
               (n.Kj[j]*2u+n.pKj[j])*h); }

ErrorType oneparam_worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h) {
    return pow(h,2u)*((n.K+n.pK)*n.pL/3u +
                      n.pK*2u*(n.L+n.pL)*n.expLambda); }
ErrorType oneparam_component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j) {
    return n.pLj[j]*(n.K+n.pK)*pow(h,2u)/3u +
           ((n.Lj[j]+n.pLj[j])*2u*n.pK)*cast_positive(cast_exact((n.L*n.expL*h+1u-n.expL)/pow(n.L,2u))); }

template<> ErrorType twoparam_worstcase_error<AffineInputs>(C1Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r) {
    return ((r*r+1u)*n.pL*n.pK +
            (r+1u)*h*n.pK*((n.pH*2u*r + n.H)*(n.K+r*n.pK)+pow(n.L,2u)+(n.L*3u*r+n.pL*r*r*2u)*n.pL)*n.expLambda +
            (r+1u)/6u*h*(n.K+n.pK)*((n.H*n.pK+n.L*n.pL)*3u+(n.pH*n.K+n.L*n.pL)*4u)
           )/cast_positive(+1u-h*n.L/2u-h*n.pL*r)*pow(h,2u)/4u; }
template<> ErrorType twoparam_worstcase_error<AdditiveInputs>(C1Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r) {
    return (n.H*(n.K+n.pK)/2u +
            (pow(n.L,2u)+n.H*(n.K+r*n.pK))*n.expLambda
           )/cast_positive(+1u-h*n.L/2u)*(r+1u)*n.pK*pow(h,3u)/4u; }
template<> ErrorType twoparam_worstcase_error<SingularInput>(C1Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r) {
    return ((r+1u)*n.pK*((n.pH*2u*r+n.H)*(n.K+r*n.pK)+pow(n.L,2u)+(n.L*3u*r+pow(r,2u)*2u*n.pL)*n.pL)*n.expLambda +
            (n.K+n.pK)/6u*((r+1u)*((n.H*n.pK+n.L*n.pL)*3u +(n.pH*n.K+n.L*n.pL)*4u) +
                           (n.pH*n.pK+pow(n.pL,2u))*8u*(r*r+1u))
           )*pow(h,3u)/4u/cast_positive(+1u-h*n.L/2u-h*n.pL*r); }

template<> ErrorType twoparam_component_error<AffineInputs>(C1Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r, SizeType j) {
    return (pow(h,2u)*(pow(r,2u)+1u)*n.pK*n.pLj[j]/2u +
            pow(h,3u)*(n.K+n.pK)*(r+1u)*((n.Hj[j]*n.pK+n.Lj[j]*n.pL)/8u+(n.pHj[j]*n.K+n.L*n.pLj[j])/6u) +
            n.pK*(r+1u)*((n.Lj[j]*n.L+r*n.pL*n.Lj[j]+n.Hj[j]*(n.K+r*n.pK))/2u*cast_positive(cast_exact((n.expL*(pow(h*n.L,2u)*3u+4u-h*n.L*5u)+h*n.L-4u)/pow(n.L,3u))) +
                                                 (n.pLj[j]*n.L+r*n.pL*n.pLj[j]+n.pHj[j]*(n.K+r*n.pK))*r*cast_positive(cast_exact((n.expL*(pow(h*n.L,2u)+2u-h*n.L*2u)-2u)/pow(n.L,3u))))
           )/cast_positive(1u-h*n.Lj[j]/2u-h*r*n.pLj[j]); }
template<> ErrorType twoparam_component_error<AdditiveInputs>(C1Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r, SizeType j) {
    return (n.Hj[j]*(n.K+n.pK)/2u +
            (n.Lj[j]*n.L+n.Hj[j]*(n.K+r*n.pK))*cast_positive(cast_exact((n.expL*(pow(h*n.L,2u)*3u+4u-h*n.L*5u)+h*n.L-4u)/pow(n.L*h,3u)))/2u
           )*n.pK*pow(h,3u)/4u*(r+1u)/cast_positive(+1u-h*n.Lj[j]/2u); }
template<> ErrorType twoparam_component_error<SingularInput>(C1Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r, SizeType j) {
    return (pow(h,3u)*(n.K+n.pK)/24u*((r+1u)*((n.Hj[j]*n.pK+n.Lj[j]*n.pL)*3u+(n.pHj[j]*n.K+n.L*n.pLj[j])*4u) +
                                      (n.pHj[j]*n.pK+n.pL+n.pLj[j])*(pow(r,2u)+1u)*8u) +
            n.pK*(r+1u)*((n.Lj[j]*n.L+r*n.pL*n.Lj[j]+n.Hj[j]*(n.K+r*n.pK))/2u*cast_positive(cast_exact((n.expL*(pow(h*n.L,2u)*3u+4u-h*n.L*5u)+h*n.L-4u)/pow(n.L,3u))) +
                         (n.pLj[j]*n.L+r*n.pL*n.pLj[j]+n.pHj[j]*(n.K+r*n.pK))*r*cast_positive(cast_exact((n.expL*(pow(h*n.L,2u)+2u-h*n.L*2u)-2u)/pow(n.L,3u))))
           )/cast_positive(1u-h*n.Lj[j]/2u-h*r*n.pLj[j]); }


template<class A, class R> Vector<ErrorType> ApproximationErrorProcessor<A,R>::process(C1Norms const& n, PositiveFloatDPValue const& h) const {

    Vector<ErrorType> result(n.dimension(),worstcase_error<A,R>(n,h));
    
    if (_enable_componentwise_error) {
        for (auto j: range(n.dimension()))
            result[j] = min(result[j],component_error<A,R>(n,h,j));
    }
    return result;
}

template<class A, class R> Vector<ErrorType> ApproximationErrorProcessor<A,R>::process(PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    C1Norms norms = compute_norms(_ivf,h,B);
    ARIADNE_LOG(7,"norms: " << norms << "\n");
    if (_ivf.is_input_additive())
        norms.pK=mag(norm(_ivf.inputs()));
    return process(norms,h);
}

InputApproximator
InputApproximatorFactory::create(InclusionVectorField const& ivf, InputApproximationKind kind, SweeperDP sweeper) const {
    switch(kind) {
    case InputApproximationKind::ZERO : return InputApproximator(SharedPointer<InputApproximatorInterface>(new InputApproximatorBase<ZeroApproximation>(ivf,sweeper)));
    case InputApproximationKind::CONSTANT : return InputApproximator(SharedPointer<InputApproximatorInterface>(new InputApproximatorBase<ConstantApproximation>(ivf,sweeper)));
    case InputApproximationKind::AFFINE : return InputApproximator(SharedPointer<InputApproximatorInterface>(new InputApproximatorBase<AffineApproximation>(ivf,sweeper)));
    case InputApproximationKind::SINUSOIDAL: return InputApproximator(SharedPointer<InputApproximatorInterface>(new InputApproximatorBase<SinusoidalApproximation>(ivf,sweeper)));
    case InputApproximationKind::PIECEWISE : return InputApproximator(SharedPointer<InputApproximatorInterface>(new InputApproximatorBase<PiecewiseApproximation>(ivf,sweeper)));
    default:
        ARIADNE_FAIL_MSG("Unexpected input approximation kind "<<kind<<"\n");
    }
}


InclusionIntegrator::InclusionIntegrator(List<InputApproximationKind> approximations, SweeperDP sweeper, StepSize step_size_)
    : _approximations(approximations)
    , _sweeper(sweeper)
    , _step_size(step_size_)
    , _number_of_steps_between_simplifications(8)
    , _number_of_variables_to_keep(4)
{
    assert(approximations.size()>0);
}

static const SizeType NUMBER_OF_PICARD_ITERATES=6;

Bool InclusionIntegrator::must_recondition(Nat step) const {
    return (step%this->_number_of_steps_between_simplifications == this->_number_of_steps_between_simplifications-1);
}

List<ValidatedVectorMultivariateFunctionModelDP> InclusionIntegrator::flow(InclusionVectorField const& ivf, BoxDomainType const& initial, Real tmax) {

    ARIADNE_LOG(2,"\nDynamics: "<<ivf<<"\n");
    ARIADNE_LOG(2,"\nInitial: "<<initial<<"\n");

    const ValidatedVectorMultivariateFunction& F = ivf.function();
    const BoxDomainType& V = ivf.inputs();
    const BoxDomainType& X0 = initial;

    auto n=ivf.dimension();
    auto m=ivf.number_of_inputs();
    auto freq=this->_number_of_steps_between_simplifications;
    DoublePrecision pr;

    PositiveFloatDPValue hsug(this->_step_size);

    ValidatedVectorMultivariateFunctionModelDP evolve_function = ValidatedVectorMultivariateTaylorFunctionModelDP::identity(X0,this->_sweeper);
    auto t=PositiveFloatDPValue(0.0);

    Map<InputApproximationKind,SizeType> approximation_global_frequencies, approximation_local_frequencies;
    Map<InputApproximationKind,Nat> delays;
    for (auto appro: _approximations) {
        approximation_global_frequencies[appro] = 0;
        approximation_local_frequencies[appro] = 0;
        delays[appro] = 0;
    }

    List<ScheduledApproximator> schedule;

    List<InputApproximator> approximations;
    InputApproximatorFactory factory;
    for (auto appro : _approximations)
        approximations.append(factory.create(ivf,appro,_sweeper));

    for (auto appro: approximations) {
        schedule.push_back(ScheduledApproximator(0u,appro));
    }

    List<ValidatedVectorMultivariateFunctionModelDP> result;
    Nat step = 0u;
    while (possibly(t<FloatDPBounds(tmax,pr))) {

        if (verbosity == 1)
            std::cout << "\r[" << activity_symbol(step) << "] " << static_cast<int>(std::round(100*t.get_d()/tmax.get_d())) << "% " << std::flush;

        ARIADNE_LOG(3,"step#:"<<step<<", t:"<<t<<", hsug:"<<hsug << "\n");

        List<InputApproximator> approximators_to_use;
        while (!schedule.empty()) {
            auto entry = schedule.back();
            if (entry.step == step) {
                approximators_to_use.push_back(entry.approximator);
                schedule.pop_back();
            } else if (entry.step > step) {
                break;
            }
        }

        if(possibly(t+hsug>FloatDPBounds(tmax,pr))) {  //FIXME: Check types for timing;
            hsug=cast_positive(cast_exact((tmax-t).upper()));
        }

        ARIADNE_LOG(4,"n. of parameters="<<evolve_function.argument_size()<<"\n");

        auto D = cast_exact_box(evolve_function.range());
        UpperBoxType B;
        StepSizeType h;

        BoxDomainType dom = product(D,V);

        std::tie(h,B)=this->flow_bounds(F,dom,static_cast<StepSizeType>(cast_exact(hsug)));
        ARIADNE_LOG(3,"flow bounds = "<<B<<" (using h = " << h << ")\n");

        PositiveFloatDPValue new_t=cast_positive(cast_exact((t+h).lower()));

        ValidatedVectorMultivariateFunctionModelDP reach_function;
        ValidatedVectorMultivariateFunctionModelDP best_reach_function, best_evolve_function;
        SharedPointer<InputApproximator> best;
        FloatDP best_volume(0);

        ARIADNE_LOG(4,"n. of approximations to use="<<approximators_to_use.size()<<"\n");

        for (auto i : range(approximators_to_use.size())) {
            this->_approximator = SharedPointer<InputApproximator>(new InputApproximator(approximators_to_use.at(i)));
            ARIADNE_LOG(5,"checking "<<this->_approximator->kind()<<" approximation\n");

            auto current_reach=reach(ivf,D,evolve_function,B,t,PositiveFloatDPValue(h,DoublePrecision()));
            auto current_evolve=evaluate_evolve_function(current_reach,new_t);

            if (i == 0) {
                best_reach_function = current_reach;
                best_evolve_function = current_evolve;
                best = this->_approximator;
                best_volume = volume(best_evolve_function.range());
            } else {
                FloatDP current_volume = volume(current_evolve.range());
                if (current_volume < best_volume) {
                    best = this->_approximator;
                    ARIADNE_LOG(6,"best approximation: " << best->kind() << "\n");
                    best_reach_function = current_reach;
                    best_evolve_function = current_evolve;
                    best_volume = current_volume;
                }
            }
        }

        if (approximators_to_use.size() > 1)
            ARIADNE_LOG(4,"chosen approximation: " << best->kind() << "\n");

        for (auto appro : approximators_to_use) {
            if (best->kind() == appro.kind())
                delays[appro.kind()] = 0;
            else
                delays[appro.kind()]++;

            Nat offset = 1u<<delays[appro.kind()];
            schedule.push_back(ScheduledApproximator(step+offset,appro));
        }
        std::sort(schedule.begin(),schedule.end(),ScheduledApproximatorComparator());

        ARIADNE_LOG(4,"updated schedule: " << schedule << "\n");

        approximation_global_frequencies[best->kind()] += 1;
        approximation_local_frequencies[best->kind()] += 1;

        reach_function = best_reach_function;
        evolve_function = best_evolve_function;

        if (must_recondition(step)) {

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

            Nat num_variables_to_keep(base);
            ARIADNE_LOG(5,"simplifying to "<<num_variables_to_keep<<" variables\n");
            LohnerReconditioner& lreconditioner = dynamic_cast<LohnerReconditioner&>(*this->_reconditioner);
            lreconditioner.set_number_of_variables_to_keep(num_variables_to_keep);
            lreconditioner.simplify(evolve_function);

            for (auto appro: _approximations) {
                approximation_local_frequencies[appro] = 0;
            }
        }

        evolve_function = this->_reconditioner->expand_errors(evolve_function);

        ARIADNE_LOG(3,"evolve bounds="<<evolve_function.range()<<"\n");

        step+=1;

        t=new_t;
        result.append(reach_function);

    }

    ARIADNE_LOG(2,"\napproximation % ="<<convert_to_percentages(approximation_global_frequencies)<<"\n");

    return result;
}

ValidatedVectorMultivariateFunctionModelType
InclusionIntegrator::reach(InclusionVectorField const& ivf, BoxDomainType D, ValidatedVectorMultivariateFunctionModelType evolve_function, UpperBoxType B, PositiveFloatDPValue t, PositiveFloatDPValue h) const {

    auto n = ivf.dimension();
    auto m = ivf.number_of_inputs();
    auto F = ivf.function();
    auto V = ivf.inputs();
    PositiveFloatDPValue new_t=cast_positive(cast_exact((t+h).lower()));

    ValidatedVectorMultivariateFunctionModelType result;

    if (this->_approximator->kind() != InputApproximationKind::PIECEWISE) {
        auto e=this->_approximator->compute_errors(h,B);
        ARIADNE_LOG(6,"approximation errors:"<<e<<"\n");
        auto DHV = this->_approximator->build_flow_domain(D,V,h);
        ARIADNE_LOG(6,"DHV:"<<DHV<<"\n");
        auto w = this->_approximator->build_w_functions(DHV, n, m);
        ARIADNE_LOG(6,"w:"<<w<<"\n");
        auto Fw = build_Fw(F,w);
        ARIADNE_LOG(6,"Fw:"<<Fw<<"\n");
        auto phi = this->compute_flow_function(Fw,DHV,B);
        phi = add_errors(phi,e);
        result=build_reach_function(evolve_function, phi, t, new_t);
    } else {
        auto e=this->_approximator->compute_errors(h,B);
        ARIADNE_LOG(6,"approximation errors:"<<e<<"\n");
        auto DHV_hlf = this->_approximator->build_flow_domain(D,V,hlf(h));
        ARIADNE_LOG(6,"DHV_hlf:"<<DHV_hlf<<"\n");
        auto w_hlf = this->_approximator->build_w_functions(DHV_hlf, n, m);
        ARIADNE_LOG(6,"w_hlf:"<<w_hlf<<"\n");
        auto Fw_hlf = build_Fw(F,w_hlf);
        ARIADNE_LOG(6,"Fw_hlf:" << Fw_hlf << "\n");
        auto phi_hlf = this->compute_flow_function(Fw_hlf,DHV_hlf,B);
        PositiveFloatDPValue intermediate_t=cast_positive(cast_exact((t+hlf(h)).lower()));
        auto intermediate_reach=build_reach_function(evolve_function, phi_hlf, t, intermediate_t);
        auto intermediate_evolve=evaluate_evolve_function(intermediate_reach,intermediate_t);

        auto D_int = cast_exact_box(intermediate_evolve.range());

        auto DHV=this->_approximator->build_flow_domain(D_int,V,hlf(h));
        auto w = build_secondhalf_piecewise_w_functions(DHV, n, m);
        ARIADNE_LOG(6,"w:"<<w<<"\n");
        auto Fw = build_Fw(F, w);
        ARIADNE_LOG(6,"Fw:"<<Fw<<"\n");
        auto phi = this->compute_flow_function(Fw,DHV,B);
        phi = add_errors(phi,e);
        result = build_secondhalf_piecewise_reach_function(intermediate_evolve, phi, intermediate_t, new_t);
    }
    return result;
}

Vector<ValidatedScalarMultivariateFunction> InclusionIntegrator::build_secondhalf_piecewise_w_functions(BoxDomainType DHV, SizeType n, SizeType m) const {
    auto zero = ValidatedScalarMultivariateFunction::zero(n+1+2*m);
    auto one = ValidatedScalarMultivariateFunction::constant(n+1+2*m,1_z);

    auto result = Vector<ValidatedScalarMultivariateFunction>(m);
    for (auto i : range(m)) {
        auto Vi = ExactNumber(DHV[n+1+i].upper());
        auto p0 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+i);
        auto p1 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+m+i);
        result[i] = (definitely (DHV[n+1+i].upper() == 0.0_exact) ? zero : p0+(one-p0*p0/Vi/Vi)*p1);
    }
    return result;
}

ValidatedVectorMultivariateFunctionModelDP InclusionIntegrator::build_secondhalf_piecewise_reach_function(
        ValidatedVectorMultivariateFunctionModelDP evolve_function, ValidatedVectorMultivariateFunctionModelDP Phi, PositiveFloatDPValue t,
        PositiveFloatDPValue new_t) const {

    // Evolve function is e(x,a,b) at s; Flow is phi(x,h,b)
    // Want (x,t,a,b):->phi(e(x,a,b),t-s,b))

    SizeType n=evolve_function.result_size();
    SizeType b=Phi.argument_size()-(n+1);

    SizeType a=evolve_function.argument_size()-n-b;

    BoxDomainType X=evolve_function.domain()[range(0,n)];
    BoxDomainType PA=evolve_function.domain()[range(n,n+a)];
    BoxDomainType PB=Phi.domain()[range(n+1,n+1+b)];

    auto swp=this->_sweeper;
    auto Tau=IntervalDomainType(t,new_t);
    BoxDomainType XTP = join(X,Tau,PA,PB);
    ValidatedVectorMultivariateTaylorFunctionModelDP xf=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(0,n),swp);
    ValidatedScalarMultivariateTaylorFunctionModelDP tf=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(XTP,n,swp);
    ValidatedVectorMultivariateTaylorFunctionModelDP af=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(n+1,n+1+a),swp);
    ValidatedVectorMultivariateTaylorFunctionModelDP bf=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(n+1+a,n+1+a+b),swp);

    ValidatedScalarMultivariateTaylorFunctionModelDP hf=tf-t;

    ValidatedVectorMultivariateTaylorFunctionModelDP ef=compose(evolve_function,join(xf,af,bf));

    return compose(Phi,join(ef,hf,bf));
}

ValidatedVectorMultivariateFunctionModelDP InclusionIntegrator::build_reach_function(
        ValidatedVectorMultivariateFunctionModelDP evolve_function, ValidatedVectorMultivariateFunctionModelDP Phi, PositiveFloatDPValue t,
        PositiveFloatDPValue new_t) const {

    // Evolve function is e(x,a) at s; flow is phi(x,h,b)
    // Want (x,t,a,b):->phi(e(x,a),t-s,b))

    SizeType n=evolve_function.result_size();

    SizeType a=evolve_function.argument_size()-n;
    SizeType b=Phi.argument_size()-(n+1);

    BoxDomainType X=evolve_function.domain()[range(0,n)];
    BoxDomainType PA=evolve_function.domain()[range(n,n+a)];
    BoxDomainType PB=Phi.domain()[range(n+1,n+1+b)];

    auto swp=this->_sweeper;
    auto Tau=IntervalDomainType(t,new_t);
    BoxDomainType XTP = join(X,Tau,PA,PB);
    ValidatedVectorMultivariateTaylorFunctionModelDP xf=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(0,n),swp);
    ValidatedScalarMultivariateTaylorFunctionModelDP tf=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(XTP,n,swp);
    ValidatedVectorMultivariateTaylorFunctionModelDP af=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(n+1,n+1+a),swp);
    ValidatedVectorMultivariateTaylorFunctionModelDP bf=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(n+1+a,n+1+a+b),swp);

    ValidatedScalarMultivariateTaylorFunctionModelDP hf=tf-t;

    ValidatedVectorMultivariateTaylorFunctionModelDP ef=compose(evolve_function,join(xf,af));

    return compose(Phi,join(ef,hf,bf));
}

ValidatedVectorMultivariateFunctionModelDP InclusionIntegrator::evaluate_evolve_function(ValidatedVectorMultivariateFunctionModelDP reach_function, PositiveFloatDPValue t) const {
    return partial_evaluate(reach_function,reach_function.result_size(),t);
}

ValidatedVectorMultivariateFunctionModelDP add_errors(ValidatedVectorMultivariateFunctionModelDP phi, Vector<ErrorType> const& e) {
    assert(phi.result_size()==e.size());
    ValidatedVectorMultivariateTaylorFunctionModelDP& tphi = dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP&>(phi.reference());
    for (auto i : range(e.size())) {
        tphi[i].add_error(e[i]);
    }
    return phi;
}

ValidatedVectorMultivariateFunction build_Fw(ValidatedVectorMultivariateFunction const& F, Vector<ValidatedScalarMultivariateFunction> const& w) {

    auto n = F.result_size();
    auto m = w.size();
    auto p = w[0].argument_size();

    auto coordinates = ValidatedVectorMultivariateFunction::coordinates(p);

    auto substitution = ValidatedVectorMultivariateFunction::zeros(n+m,p);
    for (auto i : range(n)) {
        substitution.set(i,coordinates[i]);
    }
    for (auto i : range(m)) {
        substitution.set(n+i,w[i]);
    }

    return compose(F,substitution);
}


Pair<StepSizeType,UpperBoxType> InclusionIntegrator::flow_bounds(ValidatedVectorMultivariateFunction f, BoxDomainType dom, StepSizeType hsug) const {
    BoxDomainType D=project(dom,range(0,f.result_size()));
    BoxDomainType A=project(dom,range(f.result_size(),f.argument_size()));

    return EulerBounder().compute(f,D,A,hsug);
}


ValidatedVectorMultivariateFunctionModelDP InclusionIntegrator::
compute_flow_function(ValidatedVectorMultivariateFunction const& dyn, BoxDomainType const& domain, UpperBoxType const& B) const {
    auto n=dyn.result_size();
    auto swp=this->_sweeper;

    auto x0f=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(domain,range(n),swp);
    auto af=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(domain,range(n,dyn.argument_size()),swp);

    auto picardPhi=ValidatedVectorMultivariateTaylorFunctionModelDP(n,domain,swp);
    picardPhi=picardPhi+cast_singleton(B);

    for(Nat i=0; i<NUMBER_OF_PICARD_ITERATES; ++i) {
        auto dyn_of_phi = compose(dyn,join(picardPhi,af));
        picardPhi=antiderivative(dyn_of_phi,dyn_of_phi.result_size())+x0f;
    }

    return picardPhi;
}

template<class A> BoxDomainType InputApproximatorBase<A>::build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h) const {
    auto result = D;
    result = product(result,IntervalDomainType(-h,+h));
    for (Nat i=0; i<this->_num_params_per_input; ++i)
        result = product(result,V);
    return result;
}

template<> Vector<ValidatedScalarMultivariateFunction> InputApproximatorBase<ZeroApproximation>::build_w_functions(BoxDomainType DHV, SizeType n, SizeType m) const {
    auto result = Vector<ValidatedScalarMultivariateFunction>(m);
    for (auto i : range(0,m))
        result[i] = ValidatedScalarMultivariateFunction::zero(n+1);
    return result;
}


template<> Vector<ValidatedScalarMultivariateFunction> InputApproximatorBase<ConstantApproximation>::build_w_functions(BoxDomainType DHV, SizeType n, SizeType m) const {
    auto result = Vector<ValidatedScalarMultivariateFunction>(m);
    for (auto i : range(0,m))
        result[i] = ValidatedScalarMultivariateFunction::coordinate(n+1+m,n+1+i);
    return result;
}


template<> Vector<ValidatedScalarMultivariateFunction> InputApproximatorBase<AffineApproximation>::build_w_functions(BoxDomainType DHV, SizeType n, SizeType m) const {
    auto zero = ValidatedScalarMultivariateFunction::zero(n+1+2*m);
    auto one = ValidatedScalarMultivariateFunction::constant(n+1+2*m,1_z);
    auto three = ValidatedScalarMultivariateFunction::constant(n+1+2*m,3_z);
    auto t = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n);
    auto h = ValidatedScalarMultivariateFunction::constant(n+1+2*m,ExactNumber(DHV[n].upper()));

    auto result = Vector<ValidatedScalarMultivariateFunction>(m);
    for (auto i : range(m)) {
        auto Vi = ExactNumber(DHV[n+1+i].upper());
        auto p0 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+i);
        auto p1 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+m+i);
        result[i] = (definitely (DHV[n+1+i].upper() == 0.0_exact) ? zero : p0+three*(one-p0*p0/Vi/Vi)*p1*(t-h/2)/h);
    }
    return result;
}


template<> Vector<ValidatedScalarMultivariateFunction> InputApproximatorBase<SinusoidalApproximation>::build_w_functions(BoxDomainType DHV, SizeType n, SizeType m) const {
    auto zero = ValidatedScalarMultivariateFunction::zero(n+1+2*m);
    auto one = ValidatedScalarMultivariateFunction::constant(n+1+2*m,1_z);
    auto pgamma = ValidatedScalarMultivariateFunction::constant(n+1+2*m,1.1464_dec);
    auto gamma = ValidatedScalarMultivariateFunction::constant(n+1+2*m,4.162586_dec);
    auto t = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n);
    auto h = ValidatedScalarMultivariateFunction::constant(n+1+2*m,ExactNumber(DHV[n].upper()));

    auto result = Vector<ValidatedScalarMultivariateFunction>(m);
    for (auto i : range(m)) {
        auto Vi = ExactNumber(DHV[n+1+i].upper());
        auto p0 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+i);
        auto p1 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+m+i);
        result[i] = (definitely (DHV[n+1+i].upper() == 0.0_exact) ? zero : p0+(one-p0*p0/Vi/Vi)*pgamma*p1*sin((t-h/2)*gamma/h));
    }
    return result;
}


template<> Vector<ValidatedScalarMultivariateFunction> InputApproximatorBase<PiecewiseApproximation>::build_w_functions(BoxDomainType DHV, SizeType n, SizeType m) const {
    auto zero = ValidatedScalarMultivariateFunction::zero(n+1+2*m);
    auto one = ValidatedScalarMultivariateFunction::constant(n+1+2*m,1_z);

    auto result = Vector<ValidatedScalarMultivariateFunction>(m);
    for (auto i : range(m)) {
        auto Vi = ExactNumber(DHV[n+1+i].upper());
        auto p0 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+i);
        auto p1 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+m+i);
        result[i] = (definitely (DHV[n+1+i].upper() == 0.0_exact) ? zero : p0-(one-p0*p0/Vi/Vi)*p1);
    }
    return result;
}


LohnerReconditioner::LohnerReconditioner(SweeperDP sweeper, Nat number_of_variables_to_keep_)
    : _sweeper(sweeper), _number_of_variables_to_keep(number_of_variables_to_keep_) {
    this->verbosity = 0;
}

ValidatedVectorMultivariateFunctionModelDP LohnerReconditioner::expand_errors(ValidatedVectorMultivariateFunctionModelDP f) const {
    BoxDomainType domain=f.domain();
    BoxDomainType errors=cast_exact(cast_exact(f.errors())*FloatDPUpperInterval(-1,+1)); // FIXME: Avoid cast;

    ARIADNE_LOG(6,"Uniform errors:"<<errors<<"\n");
    for(SizeType i=0; i!=f.result_size(); ++i) { f[i].set_error(0); }
    ValidatedVectorMultivariateFunctionModelDP error_function=ValidatedVectorMultivariateTaylorFunctionModelDP::identity(errors,this->_sweeper);
    return embed(f,errors)+embed(domain,error_function);
}

struct IndexedFloatDPError
{
    SizeType index;
    FloatDPError value;

    IndexedFloatDPError() : index(0), value(FloatDPError()) {}
};

inline OutputStream& operator<<(OutputStream& os, IndexedFloatDPError const& ifl) {
    return os << "(" << ifl.index << ":" << std::scientific << ifl.value.raw() << std::fixed << ")"; }

struct IndexedFloatDPErrorComparator
{
    inline bool operator() (const IndexedFloatDPError& ifl1, const IndexedFloatDPError& ifl2)
    {
        return (ifl1.value.raw() < ifl2.value.raw());
    }
};

Void LohnerReconditioner::simplify(ValidatedVectorMultivariateFunctionModelDP& f) const {
    ARIADNE_LOG(6,"simplifying\n");
    ARIADNE_LOG(6,"f="<<f<<"\n");

    auto m=f.argument_size();
    auto n=f.result_size();

    ARIADNE_LOG(6,"num.parameters="<<m<<", to keep="<< this->_number_of_variables_to_keep <<"\n");

    ValidatedVectorMultivariateTaylorFunctionModelDP& tf = dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP&>(f.reference());

    // Compute effect of error terms, but not of original variables;
    Matrix<FloatDPError> C(m,n);
    for (auto i : range(n)) {
        auto p=tf[i].model().expansion();

        for (auto ac : p) {
            UniformConstReference<MultiIndex> a=ac.index();
            UniformReference<FloatDPValue> c=ac.coefficient();
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

    if (m <= this->_number_of_variables_to_keep) {
        ARIADNE_LOG(6, "Insufficient number of variables, not simplifying\n");
        return;
    }

    Nat number_of_variables_to_remove = m - this->_number_of_variables_to_keep;
    ARIADNE_LOG(6, "Number of variables to remove:" << number_of_variables_to_remove<<"\n");

    for (auto j : range(number_of_variables_to_remove)) {
        remove_indices.append(SCe[j].index);
    }

    for (auto j : range(number_of_variables_to_remove,m)) {
        keep_indices.append(SCe[j].index);
    }

    ARIADNE_LOG(2,"number of kept parameters: " << keep_indices.size() << "/" << m << "\n");

    ARIADNE_LOG(6,"keep_indices:"<<keep_indices<<"\n");
    ARIADNE_LOG(6,"remove_indices:"<<remove_indices<<"\n");

    for (auto i : range(n)) {
        ErrorType error = tf[i].error();
        for(SizeType k=0; k!=remove_indices.size(); ++k) {
            error += mag(C[remove_indices[k]][i]);
        }
        tf[i].set_error(error);
    }

    auto old_domain=f.domain();
    auto new_domain=BoxDomainType(Vector<IntervalDomainType>(keep_indices.size(),[&old_domain,&keep_indices](SizeType j){return old_domain[keep_indices[j]];}));
    auto projection=ValidatedVectorMultivariateTaylorFunctionModelDP(m,new_domain,this->_sweeper);
    for (auto i : range(new_domain.size())) { projection[keep_indices[i]]=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(new_domain,i,this->_sweeper); }
    for (auto i : range(remove_indices.size())) {
        auto j=remove_indices[i]; auto cj=old_domain[j].midpoint();
        projection[j]=ValidatedScalarMultivariateTaylorFunctionModelDP::constant(new_domain,cj,this->_sweeper); }
    f=compose(f,projection);
}

} // namespace Ariadne;
