/***************************************************************************
 *            dynamics/differential_inclusion.cpp
 *
 *  Copyright  2008-20  Luca Geretti, Pieter Collins, Sanja Zivanovic
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
#include "inclusion_evolver.hpp"

namespace Ariadne {

TimeStepType lower_bound(TimeStepType const& t) {
    return TimeStepType(FloatDPLowerBound(t,DoublePrecision()).raw());
}

TimeStepType lower_bound(Real const& t) {
    return TimeStepType(FloatDPLowerBound(t,DoublePrecision()).raw());
}

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
    InclusionApproximatorHandle approximator;

    ScheduledApproximator(Nat s, InclusionApproximatorHandle a) : step(s), approximator(a) {}
};

inline OutputStream& operator<<(OutputStream& os, ScheduledApproximator const& sa) {
    return os << "(" << sa.step << ":" << sa.approximator << ")"; }

struct ScheduledApproximatorComparator {
    inline bool operator() (ScheduledApproximator const& sa1, ScheduledApproximator const& sa2) {
        return (sa1.step > sa2.step); }
};


class InclusionEvolverState {
  private:
    Nat _step;
    List<ScheduledApproximator> _schedule;
    Map<InclusionApproximatorHandle,Nat> _approximator_check_delay;
    Map<InclusionApproximatorHandle,Nat> _approximator_global_optima_count;
    Map<InclusionApproximatorHandle,Nat> _approximator_local_optima_count;
  public:
    InclusionEvolverState(InclusionVectorField const& ivf, List<InputApproximation> const& approximations, IntegratorInterface const& integrator)
        : _step(0u)
    {
        InclusionApproximatorFactory factory(integrator);
        for (auto appro : approximations) {
            InclusionApproximatorHandle approximator = factory.create(ivf,appro);
            _schedule.push_back(ScheduledApproximator(0u,approximator));
            _approximator_global_optima_count[approximator] = 0;
            _approximator_local_optima_count[approximator] = 0;
            _approximator_check_delay[approximator] = 0;
        }
    }

    Nat step() const { return _step; }
    List<ScheduledApproximator> const& schedule() const { return _schedule; }
    Map<InclusionApproximatorHandle,Nat> const& check_delay() const { return _approximator_check_delay; }
    Map<InclusionApproximatorHandle,Nat> const& global_optima_count() const { return _approximator_global_optima_count; }
    Map<InclusionApproximatorHandle,Nat> const& local_optima_count() const { return _approximator_local_optima_count; }

    Void reset_local_optima_count() {
        for (auto entry : _approximator_local_optima_count)
            _approximator_local_optima_count[entry.first] = 0;
    }

    Void update_with_best(InclusionApproximatorHandle const& best) {
        for (auto appro : approximators_to_use()) {
            if (best == appro) _reset_check_delay(appro);
            else _increase_check_delay(appro);

            Nat offset = 1u<<_approximator_check_delay[appro];
            _schedule.push_back(ScheduledApproximator(_step+offset,appro));
        }
        std::sort(_schedule.begin(),_schedule.end(),ScheduledApproximatorComparator());

        _increase_optima_count(best);
    }

    Void append_to_schedule(Nat const& step, InclusionApproximatorHandle const& approximator) {
        _schedule.push_back(ScheduledApproximator(step,approximator));
    }

    List<InclusionApproximatorHandle> approximators_to_use() const {
        List<InclusionApproximatorHandle> result;
        for (auto i = _schedule.size(); i > 0; --i) {
            auto entry = _schedule[i-1];
            if (entry.step == _step) {
                result.push_back(entry.approximator);
            } else if (entry.step > _step) {
                break;
            }
        }
        return result;
    }

    Void next_step() {
        while (!_schedule.empty()) {
            auto entry = _schedule.back();
            if (entry.step == _step) {
                _schedule.pop_back();
            } else if (entry.step > _step) {
                break;
            }
        }
        _step++;
    }

  private:
    Void _reset_check_delay(InclusionApproximatorHandle const& approximator) { _approximator_check_delay[approximator] = 0; }
    Void _increase_check_delay(InclusionApproximatorHandle const& approximator) { _approximator_check_delay[approximator]++; }
    Void _increase_optima_count(InclusionApproximatorHandle const& approximator) {
        _approximator_global_optima_count[approximator]++;
        _approximator_local_optima_count[approximator]++;
    }
};


inline char activity_symbol(SizeType step) {
    switch (step % 4) {
    case 0: return '\\';
    case 1: return '|';
    case 2: return '/';
    default: return '-';
    }
}


inline Map<InclusionApproximatorHandle,FloatDP> convert_to_percentages(Map<InclusionApproximatorHandle,Nat> const& approximation_global_frequencies) {

    Nat total_steps(0);
    for (auto entry: approximation_global_frequencies) {
        total_steps += entry.second;
    }

    Map<InclusionApproximatorHandle,FloatDP> result;
    for (auto entry: approximation_global_frequencies) {
        result[entry.first] = 1.0/total_steps*entry.second;
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

FloatDPError zeroparam_worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h) {
    return min(n.pK*n.expLambda*h,
               (n.K*2u+n.pK)*h); }
FloatDPError zeroparam_component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j) {
    return min(n.pK*n.expL*h,
               (n.Kj[j]*2u+n.pKj[j])*h); }

FloatDPError oneparam_worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h) {
    return pow(h,2u)*((n.K+n.pK)*n.pL/3u +
                      n.pK*2u*(n.L+n.pL)*n.expLambda); }
FloatDPError oneparam_component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j) {
    return n.pLj[j]*(n.K+n.pK)*pow(h,2u)/3u +
           ((n.Lj[j]+n.pLj[j])*2u*n.pK)*cast_positive(cast_exact((n.L*n.expL*h+1u-n.expL)/pow(n.L,2u))); }

template<> FloatDPError twoparam_worstcase_error<AffineInputs>(C1Norms const& n, PositiveFloatDPValue const& h, FloatDPError const& r) {
    return ((r*r+1u)*n.pL*n.pK +
            (r+1u)*h*n.pK*((n.pH*2u*r + n.H)*(n.K+r*n.pK)+pow(n.L,2u)+(n.L*3u*r+n.pL*r*r*2u)*n.pL)*n.expLambda +
            (r+1u)/6u*h*(n.K+n.pK)*((n.H*n.pK+n.L*n.pL)*3u+(n.pH*n.K+n.L*n.pL)*4u)
           )/cast_positive(+1u-h*n.L/2u-h*n.pL*r)*pow(h,2u)/4u; }
template<> FloatDPError twoparam_worstcase_error<AdditiveInputs>(C1Norms const& n, PositiveFloatDPValue const& h, FloatDPError const& r) {
    return (n.H*(n.K+n.pK)/2u +
            (pow(n.L,2u)+n.H*(n.K+r*n.pK))*n.expLambda
           )/cast_positive(+1u-h*n.L/2u)*(r+1u)*n.pK*pow(h,3u)/4u; }
template<> FloatDPError twoparam_worstcase_error<SingularInput>(C1Norms const& n, PositiveFloatDPValue const& h, FloatDPError const& r) {
    return ((r+1u)*n.pK*((n.pH*2u*r+n.H)*(n.K+r*n.pK)+pow(n.L,2u)+(n.L*3u*r+pow(r,2u)*2u*n.pL)*n.pL)*n.expLambda +
            (n.K+n.pK)/6u*((r+1u)*((n.H*n.pK+n.L*n.pL)*3u +(n.pH*n.K+n.L*n.pL)*4u) +
                           (n.pH*n.pK+pow(n.pL,2u))*8u*(r*r+1u))
           )*pow(h,3u)/4u/cast_positive(+1u-h*n.L/2u-h*n.pL*r); }

template<> FloatDPError twoparam_component_error<AffineInputs>(C1Norms const& n, PositiveFloatDPValue const& h, FloatDPError const& r, SizeType j) {
    return (pow(h,2u)*(pow(r,2u)+1u)*n.pK*n.pLj[j]/2u +
            pow(h,3u)*(n.K+n.pK)*(r+1u)*((n.Hj[j]*n.pK+n.Lj[j]*n.pL)/8u+(n.pHj[j]*n.K+n.L*n.pLj[j])/6u) +
            n.pK*(r+1u)*((n.Lj[j]*n.L+r*n.pL*n.Lj[j]+n.Hj[j]*(n.K+r*n.pK))/2u*cast_positive(cast_exact((n.expL*(pow(h*n.L,2u)*3u+4u-h*n.L*5u)+h*n.L-4u)/pow(n.L,3u))) +
                                                 (n.pLj[j]*n.L+r*n.pL*n.pLj[j]+n.pHj[j]*(n.K+r*n.pK))*r*cast_positive(cast_exact((n.expL*(pow(h*n.L,2u)+2u-h*n.L*2u)-2u)/pow(n.L,3u))))
           )/cast_positive(1u-h*n.Lj[j]/2u-h*r*n.pLj[j]); }
template<> FloatDPError twoparam_component_error<AdditiveInputs>(C1Norms const& n, PositiveFloatDPValue const& h, FloatDPError const& r, SizeType j) {
    return (n.Hj[j]*(n.K+n.pK)/2u +
            (n.Lj[j]*n.L+n.Hj[j]*(n.K+r*n.pK))*cast_positive(cast_exact((n.expL*(pow(h*n.L,2u)*3u+4u-h*n.L*5u)+h*n.L-4u)/pow(n.L*h,3u)))/2u
           )*n.pK*pow(h,3u)/4u*(r+1u)/cast_positive(+1u-h*n.Lj[j]/2u); }
template<> FloatDPError twoparam_component_error<SingularInput>(C1Norms const& n, PositiveFloatDPValue const& h, FloatDPError const& r, SizeType j) {
    return (pow(h,3u)*(n.K+n.pK)/24u*((r+1u)*((n.Hj[j]*n.pK+n.Lj[j]*n.pL)*3u+(n.pHj[j]*n.K+n.L*n.pLj[j])*4u) +
                                      (n.pHj[j]*n.pK+n.pL+n.pLj[j])*(pow(r,2u)+1u)*8u) +
            n.pK*(r+1u)*((n.Lj[j]*n.L+r*n.pL*n.Lj[j]+n.Hj[j]*(n.K+r*n.pK))/2u*cast_positive(cast_exact((n.expL*(pow(h*n.L,2u)*3u+4u-h*n.L*5u)+h*n.L-4u)/pow(n.L,3u))) +
                         (n.pLj[j]*n.L+r*n.pL*n.pLj[j]+n.pHj[j]*(n.K+r*n.pK))*r*cast_positive(cast_exact((n.expL*(pow(h*n.L,2u)+2u-h*n.L*2u)-2u)/pow(n.L,3u))))
           )/cast_positive(1u-h*n.Lj[j]/2u-h*r*n.pLj[j]); }


template<class A, class R> Vector<FloatDPError> ApproximationErrorProcessor<A,R>::process(C1Norms const& n, PositiveFloatDPValue const& h) const {

    Vector<FloatDPError> result(n.dimension(),worstcase_error<A,R>(n,h));

    if (_enable_componentwise_error) {
        for (auto j: range(n.dimension()))
            result[j] = min(result[j],component_error<A,R>(n,h,j));
    }
    return result;
}

template<class A, class R> Vector<FloatDPError> ApproximationErrorProcessor<A,R>::process(PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    C1Norms norms = compute_norms(_ivf,h,B);
    ARIADNE_LOG(7,"norms: " << norms << "\n");
    if (_ivf.is_input_additive())
        norms.pK=mag(norm(_ivf.inputs()));
    return process(norms,h);
}

InclusionApproximatorHandle
InclusionApproximatorFactory::create(InclusionVectorField const& ivf, InputApproximation const& approximation) const {
    if (approximation.handles(ZeroApproximation())) return InclusionApproximatorHandle(SharedPointer<InclusionApproximatorInterface>(new InclusionApproximatorBase<ZeroApproximation>(ivf,_integrator)));
    if (approximation.handles(ConstantApproximation())) return InclusionApproximatorHandle(SharedPointer<InclusionApproximatorInterface>(new InclusionApproximatorBase<ConstantApproximation>(ivf,_integrator)));
    if (approximation.handles(AffineApproximation())) return InclusionApproximatorHandle(SharedPointer<InclusionApproximatorInterface>(new InclusionApproximatorBase<AffineApproximation>(ivf,_integrator)));
    if (approximation.handles(SinusoidalApproximation())) return InclusionApproximatorHandle(SharedPointer<InclusionApproximatorInterface>(new InclusionApproximatorBase<SinusoidalApproximation>(ivf,_integrator)));
    if (approximation.handles(PiecewiseApproximation())) return InclusionApproximatorHandle(SharedPointer<InclusionApproximatorInterface>(new InclusionApproximatorBase<PiecewiseApproximation>(ivf,_integrator)));
    ARIADNE_FAIL_MSG("Unhandled input approximation " << approximation << "\n");
}


InclusionEvolver::InclusionEvolver(SystemType const& system, SweeperDP const& sweeper, IntegratorInterface const& integrator, Reconditioner const& reconditioner)
    : _system(system)
    , _sweeper(sweeper)
    , _integrator(integrator.clone())
    , _reconditioner(reconditioner)
    , _configuration(new ConfigurationType())
{
    assert(system.inputs().size >= 0);
}

Void InclusionEvolver::_recondition_and_update(ValidatedVectorMultivariateFunctionModelType& function, InclusionEvolverState& state) {
    if (_reconditioner.must_reduce_parameters(state)) {
        _reconditioner.update_from(state);
        _reconditioner.reduce_parameters(function);
        state.reset_local_optima_count();
    }

    if (_reconditioner.must_incorporate_errors(state)) {
        function = _reconditioner.incorporate_errors(function);
    }
}

List<ValidatedVectorMultivariateFunctionModelDP> InclusionEvolver::reach(BoxDomainType const& initial, Real const& tmax) {

    ARIADNE_LOG(2,"System: "<<_system<<"\n");
    ARIADNE_LOG(2,"Initial: "<<initial<<"\n");

    const ValidatedVectorMultivariateFunction& F = _system.function();
    const BoxDomainType& V = _system.inputs();
    const BoxDomainType& X0 = initial;

    StepSizeType hsug(_configuration->maximum_step_size());

    ValidatedVectorMultivariateFunctionModelDP evolve_function = ValidatedVectorMultivariateTaylorFunctionModelDP::identity(X0,this->_sweeper);

    TimeStepType t;

    InclusionEvolverState state(_system,_configuration->approximations(),*_integrator);

    List<ValidatedVectorMultivariateFunctionModelDP> result;

    while (possibly(t<lower_bound(tmax))) {

        if (verbosity == 1)
            std::cout << "\r[" << activity_symbol(state.step()) << "] " << static_cast<int>(std::round(100*t.get_d()/tmax.get_d())) << "% " << std::flush;

        ARIADNE_LOG(3,"step#:"<<state.step()<<", t:"<<t<<", hsug:"<<hsug << "\n");

        ARIADNE_LOG(4,"n. of parameters="<<evolve_function.argument_size()<<"\n");

        auto approximators_to_use = state.approximators_to_use();

        ARIADNE_LOG(4,"approximators to use="<<approximators_to_use<<"\n");

        auto domx = cast_exact_box(evolve_function.range());

        UpperBoxType B;
        StepSizeType h;
        std::tie(h,B)=approximators_to_use.at(0).flow_bounds(F,domx,V,hsug);
        ARIADNE_LOG(3,"flow bounds = "<<B<<" (using h = " << h << ")\n");

        TimeStepType new_t = lower_bound(t+h);

        ValidatedVectorMultivariateFunctionModelDP reach_function;
        ValidatedVectorMultivariateFunctionModelDP best_reach_function, best_evolve_function;
        InclusionApproximatorHandle best = approximators_to_use.at(0);
        FloatDP best_volume = FloatDP::inf(DoublePrecision());

        for (auto approximator : approximators_to_use) {
            ARIADNE_LOG(5,"checking "<<approximator<<" approximator\n");

            auto current_reach=approximator.reach(_system,domx,evolve_function,B,t,h);
            auto current_evolve=approximator.evolve(current_reach,new_t);

            FloatDP current_volume = volume(current_evolve.range());
            if (current_volume < best_volume) {
                best = approximator;
                ARIADNE_LOG(6,"best approximator: " << best << "\n");
                best_reach_function = current_reach;
                best_evolve_function = current_evolve;
                best_volume = current_volume;
            }
        }

        if (approximators_to_use.size() > 1)
            ARIADNE_LOG(4,"chosen approximator: " << best << "\n");

        state.update_with_best(best);
        ARIADNE_LOG(4,"updated schedule: " << state.schedule() << "\n");

        reach_function = best_reach_function;
        evolve_function = best_evolve_function;

        ARIADNE_LOG(3,"evolve bounds="<<evolve_function.range()<<"\n");

        result.append(reach_function);

        this->_recondition_and_update(evolve_function,state);

        state.next_step();
        t=new_t;
    }

    ARIADNE_LOG(2,"approximation % ="<<convert_to_percentages(state.global_optima_count())<<"\n");

    return result;
}

template<class A> Bool
InclusionApproximatorBase<A>::operator==(const InclusionApproximatorInterface& rhs) const {
    return instance_of<InclusionApproximatorBase<A>>(&rhs);
}

template<class A> Bool
InclusionApproximatorBase<A>::operator<(const InclusionApproximatorInterface& rhs) const {
    return this->index() < rhs.index();
}

template<class A> ValidatedVectorMultivariateFunctionModelType
InclusionApproximatorBase<A>::reach(InclusionVectorField const& ivf, BoxDomainType const& domx, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const {

    TimeStepType new_t = lower_bound(t+h);

    Interval<TimeStepType> domt(t,new_t);

    auto e=this->compute_errors(h,B);
    ARIADNE_LOG(6,"approximation errors:"<<e<<"\n");
    auto doma = this->build_parameter_domain(ivf.inputs());

    auto w = this->build_w_functions(domt,doma,ivf.dimension(),ivf.number_of_inputs());
    ARIADNE_LOG(6,"w:"<<w<<"\n");
    auto Fw = build_Fw(ivf.function(),w);
    ARIADNE_LOG(6,"Fw:"<<Fw<<"\n");
    auto phi = this->_integrator->flow_step(Fw,domx,domt,doma,B);
    add_errors(phi,e);

    return this->build_reach_function(evolve_function, phi, t, new_t);
}

template<class A> Vector<ValidatedScalarMultivariateFunction> InclusionApproximatorBase<A>::build_secondhalf_piecewise_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) const {
    auto zero = ValidatedScalarMultivariateFunction::zero(n+1+2*m);
    auto one = ValidatedScalarMultivariateFunction::constant(n+1+2*m,1_z);

    auto result = Vector<ValidatedScalarMultivariateFunction>(m);
    for (auto i : range(m)) {
        auto Vi = ExactNumber(doma[i].upper());
        auto p0 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+i);
        auto p1 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+m+i);
        result[i] = (definitely (doma[i].upper() == 0.0_exact) ? zero : p0+(one-p0*p0/Vi/Vi)*p1);
    }
    return result;
}

template<class A> ValidatedVectorMultivariateFunctionModelDP InclusionApproximatorBase<A>::build_secondhalf_piecewise_reach_function(
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
    BoxDomainType XTP = product(X,Tau,PA,PB);
    ValidatedVectorMultivariateTaylorFunctionModelDP xf=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(0,n),swp);
    ValidatedScalarMultivariateTaylorFunctionModelDP tf=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(XTP,n,swp);
    ValidatedVectorMultivariateTaylorFunctionModelDP af=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(n+1,n+1+a),swp);
    ValidatedVectorMultivariateTaylorFunctionModelDP bf=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(n+1+a,n+1+a+b),swp);

    ValidatedVectorMultivariateTaylorFunctionModelDP ef=compose(evolve_function,join(xf,af,bf));

    return compose(Phi,join(ef,tf,bf));
}

template<class A> ValidatedVectorMultivariateFunctionModelDP InclusionApproximatorBase<A>::build_reach_function(
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
    BoxDomainType XTP = product(X,Tau,PA,PB);
    ValidatedVectorMultivariateTaylorFunctionModelDP xf=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(0,n),swp);
    ValidatedScalarMultivariateTaylorFunctionModelDP tf=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(XTP,n,swp);
    ValidatedVectorMultivariateTaylorFunctionModelDP af=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(n+1,n+1+a),swp);
    ValidatedVectorMultivariateTaylorFunctionModelDP bf=ValidatedVectorMultivariateTaylorFunctionModelDP::projection(XTP,range(n+1+a,n+1+a+b),swp);

    ValidatedVectorMultivariateTaylorFunctionModelDP ef=compose(evolve_function,join(xf,af));

    return compose(Phi,join(ef,tf,bf));
}

template<class A> ValidatedVectorMultivariateFunctionModelDP InclusionApproximatorBase<A>::evolve(ValidatedVectorMultivariateFunctionModelDP const& reach_function, TimeStepType const& t) const {
    return partial_evaluate(reach_function,reach_function.result_size(),t);
}

Void add_errors(ValidatedVectorMultivariateFunctionModelDP& phi, Vector<ErrorType> const& e) {
    assert(phi.result_size()==e.size());
    ValidatedVectorMultivariateTaylorFunctionModelDP& tphi = dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP&>(phi.reference());
    for (auto i : range(e.size())) {
        tphi[i].add_error(e[i]);
    }
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


template<class A> Pair<StepSizeType,UpperBoxType> InclusionApproximatorBase<A>::flow_bounds(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& domx, BoxDomainType const& doma, StepSizeType const& hsug) const {
    return EulerBounder().compute(f,domx,doma,hsug);
}

template<class A> BoxDomainType InclusionApproximatorBase<A>::build_parameter_domain(BoxDomainType const& V) const {
    BoxDomainType result(0u);
    for (Nat i=0; i<this->_num_params_per_input; ++i)
        result = product(result,V);
    return result;
}

template<> Vector<ValidatedScalarMultivariateFunction> InclusionApproximatorBase<ZeroApproximation>::build_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) const {
    auto result = Vector<ValidatedScalarMultivariateFunction>(m);
    for (auto i : range(0,m))
        result[i] = ValidatedScalarMultivariateFunction::zero(n+1);
    return result;
}


template<> Vector<ValidatedScalarMultivariateFunction> InclusionApproximatorBase<ConstantApproximation>::build_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) const {
    auto result = Vector<ValidatedScalarMultivariateFunction>(m);
    for (auto i : range(0,m))
        result[i] = ValidatedScalarMultivariateFunction::coordinate(n+1+m,n+1+i);
    return result;
}


template<> Vector<ValidatedScalarMultivariateFunction> InclusionApproximatorBase<AffineApproximation>::build_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) const {
    auto zero = ValidatedScalarMultivariateFunction::zero(n+1+2*m);
    auto one = ValidatedScalarMultivariateFunction::constant(n+1+2*m,1_z);
    auto three = ValidatedScalarMultivariateFunction::constant(n+1+2*m,3_z);
    auto t = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n);
    auto tk = ValidatedScalarMultivariateFunction::constant(n+1+2*m,domt.lower());
    auto hc = ValidatedScalarMultivariateFunction::constant(n+1+2*m,domt.width());

    auto result = Vector<ValidatedScalarMultivariateFunction>(m);
    for (auto i : range(m)) {
        auto Vi = ExactNumber(doma[i].upper());
        auto p0 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+i);
        auto p1 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+m+i);
        result[i] = (definitely (doma[i].upper() == 0.0_exact) ? zero : p0+three*(one-p0*p0/Vi/Vi)*p1*(t-tk-hc/2)/hc);
    }
    return result;
}


template<> Vector<ValidatedScalarMultivariateFunction> InclusionApproximatorBase<SinusoidalApproximation>::build_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) const {
    auto zero = ValidatedScalarMultivariateFunction::zero(n+1+2*m);
    auto one = ValidatedScalarMultivariateFunction::constant(n+1+2*m,1_z);
    auto pgamma = ValidatedScalarMultivariateFunction::constant(n+1+2*m,1.1464_dec);
    auto gamma = ValidatedScalarMultivariateFunction::constant(n+1+2*m,4.162586_dec);
    auto t = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n);
    auto tk = ValidatedScalarMultivariateFunction::constant(n+1+2*m,domt.lower());
    auto hc = ValidatedScalarMultivariateFunction::constant(n+1+2*m,domt.width());

    auto result = Vector<ValidatedScalarMultivariateFunction>(m);
    for (auto i : range(m)) {
        auto Vi = ExactNumber(doma[i].upper());
        auto p0 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+i);
        auto p1 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+m+i);
        result[i] = (definitely (doma[i].upper() == 0.0_exact) ? zero : p0+(one-p0*p0/Vi/Vi)*pgamma*p1*sin((t-tk-hc/2)*gamma/hc));
    }
    return result;
}


template<> Vector<ValidatedScalarMultivariateFunction> InclusionApproximatorBase<PiecewiseApproximation>::build_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) const {
    auto zero = ValidatedScalarMultivariateFunction::zero(n+1+2*m);
    auto one = ValidatedScalarMultivariateFunction::constant(n+1+2*m,1_z);

    auto result = Vector<ValidatedScalarMultivariateFunction>(m);
    for (auto i : range(m)) {
        auto Vi = ExactNumber(doma[i].upper());
        auto p0 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+i);
        auto p1 = ValidatedScalarMultivariateFunction::coordinate(n+1+2*m,n+1+m+i);
        result[i] = (definitely (doma[i].upper() == 0.0_exact) ? zero : p0-(one-p0*p0/Vi/Vi)*p1);
    }
    return result;
}

template<> ValidatedVectorMultivariateFunctionModelType
InclusionApproximatorBase<PiecewiseApproximation>::reach(InclusionVectorField const& ivf, BoxDomainType const& domx, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const {

    auto n = ivf.dimension();
    auto m = ivf.number_of_inputs();
    auto F = ivf.function();
    auto V = ivf.inputs();

    TimeStepType intermediate_t = lower_bound(t+hlf(h));
    TimeStepType new_t = lower_bound(t+h);

    Interval<TimeStepType> domt_first(t,intermediate_t);
    Interval<TimeStepType> domt_second(intermediate_t,new_t);
    auto doma = this->build_parameter_domain(V);

    auto e=this->compute_errors(h,B);
    ARIADNE_LOG(6,"approximation errors:"<<e<<"\n");

    auto w_hlf = this->build_w_functions(domt_first,doma,n,m);
    ARIADNE_LOG(6,"w_hlf:"<<w_hlf<<"\n");
    auto Fw_hlf = build_Fw(F,w_hlf);
    ARIADNE_LOG(6,"Fw_hlf:" << Fw_hlf << "\n");

    auto phi_hlf = this->_integrator->flow_step(Fw_hlf,domx,domt_first,doma,B);
    auto intermediate_reach=this->build_reach_function(evolve_function, phi_hlf, t, intermediate_t);
    auto intermediate_evolve=this->evolve(intermediate_reach,intermediate_t);

    auto domx_second = cast_exact_box(intermediate_evolve.range());

    auto w = this->build_secondhalf_piecewise_w_functions(domt_second,doma,n,m);
    ARIADNE_LOG(6,"w:"<<w<<"\n");
    auto Fw = build_Fw(F, w);
    ARIADNE_LOG(6,"Fw:"<<Fw<<"\n");
    auto phi = this->_integrator->flow_step(Fw,domx_second,domt_second,doma,B);
    add_errors(phi,e);

    return this->build_secondhalf_piecewise_reach_function(intermediate_evolve, phi, intermediate_t, new_t);
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

Bool LohnerReconditioner::must_reduce_parameters(InclusionEvolverState const& state) const {
    return (state.step()%_number_of_steps_between_simplifications == _number_of_steps_between_simplifications-1);
}

Bool LohnerReconditioner::must_incorporate_errors(InclusionEvolverState const& state) const {
    return true;
}

Void LohnerReconditioner::update_from(InclusionEvolverState const& state) {
    auto freq = _number_of_steps_between_simplifications;
    auto n = _number_of_variables;
    auto m = _number_of_inputs;

    FloatDP npk = 0;
    FloatDP rho = _ratio_of_parameters_to_keep;
    for (auto entry: state.local_optima_count()) {
        SizeType ppi = entry.first.num_params_per_input();
        FloatDP partial = n + rho*(n+2*m) + (freq-1)*m*(2 - ppi);
        npk += partial*entry.second/freq;
    }

    _number_of_parameters_to_keep = static_cast<Nat>(round(npk).get_d());
}

ValidatedVectorMultivariateFunctionModelDP LohnerReconditioner::incorporate_errors(ValidatedVectorMultivariateFunctionModelDP const& f) const {

    ValidatedVectorMultivariateTaylorFunctionModelDP const& tf = dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(f.reference());

    BoxDomainType domain=f.domain();
    BoxDomainType errors=cast_exact(cast_exact(f.errors())*FloatDPUpperInterval(-1,+1)); // FIXME: Avoid cast;

    ARIADNE_LOG(6,"Uniform errors:"<<errors<<"\n");

    ValidatedVectorMultivariateFunctionModelDP error_function=ValidatedVectorMultivariateTaylorFunctionModelDP::identity(errors,tf.properties());
    ValidatedVectorMultivariateFunctionModelDP result = embed(f,errors)+embed(domain,error_function);
    for(SizeType i=0; i!=result.result_size(); ++i) { result[i].set_error(0); }
    return result;
}

Void LohnerReconditioner::reduce_parameters(ValidatedVectorMultivariateFunctionModelDP& f) const {
    ARIADNE_LOG(6,"simplifying\n");
    ARIADNE_LOG(6,"f="<<f<<"\n");

    auto m=f.argument_size();
    auto n=f.result_size();

    ARIADNE_LOG(6,"num.parameters="<<m<<", to keep="<< this->_number_of_parameters_to_keep <<"\n");

    ValidatedVectorMultivariateTaylorFunctionModelDP& tf = dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP&>(f.reference());

    auto sweeper = tf.properties();

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

    if (m <= this->_number_of_parameters_to_keep) {
        ARIADNE_LOG(6, "Insufficient number of variables, not simplifying\n");
        return;
    }

    Nat number_of_variables_to_remove = m - this->_number_of_parameters_to_keep;
    ARIADNE_LOG(6, "Number of parameters to remove:" << _number_of_parameters_to_keep<<"\n");

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
        FloatDPError error = tf[i].error();
        for(SizeType k=0; k!=remove_indices.size(); ++k) {
            error += mag(C[remove_indices[k]][i]);
        }
        tf[i].set_error(error);
    }

    auto old_domain=f.domain();
    auto new_domain=BoxDomainType(Vector<IntervalDomainType>(keep_indices.size(),[&old_domain,&keep_indices](SizeType j){return old_domain[keep_indices[j]];}));
    auto projection=ValidatedVectorMultivariateTaylorFunctionModelDP(m,new_domain,sweeper);
    for (auto i : range(new_domain.size())) { projection[keep_indices[i]]=ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate(new_domain,i,sweeper); }
    for (auto i : range(remove_indices.size())) {
        auto j=remove_indices[i]; auto cj=old_domain[j].midpoint();
        projection[j]=ValidatedScalarMultivariateTaylorFunctionModelDP::constant(new_domain,cj,sweeper); }
    f=compose(f,projection);
}


InclusionEvolverConfiguration::InclusionEvolverConfiguration()
{
    maximum_step_size(1);
    maximum_enclosure_radius(100.0);
    enable_parameter_reduction(true);
    approximations({ZeroApproximation()});
}


OutputStream&
InclusionEvolverConfiguration::write(OutputStream& os) const
{
    os << "InclusionEvolverConfiguration"
       << ",\n  maximum_step_size=" << maximum_step_size()
       << ",\n  maximum_enclosure_radius=" << maximum_enclosure_radius()
       << ",\n  enable_parameter_reduction=" << enable_parameter_reduction()
       << ",\n  approximations=" << approximations()
       << "\n)\n";
    return os;
}

} // namespace Ariadne;
