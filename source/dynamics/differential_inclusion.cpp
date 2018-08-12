/***************************************************************************
 *            differential_inclusion.cpp
 *
 *  Copyright  2008-18  Luca Geretti, Pieter Collins, Sanja Zivanovic
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

Tuple<ValidatedVectorFunction,Vector<ValidatedVectorFunction>,BoxDomainType> expression_to_function(DottedRealAssignments const& dynamics, const RealVariablesBox& inputs) {

    auto transformations = centered_variables_transformation(inputs);

    DottedRealAssignments substituted_dynamics;
    for (auto dyn : dynamics) {
        substituted_dynamics.push_back(DottedRealAssignment(dyn.left_hand_side(),substitute(dyn.right_hand_side(),transformations.first)));
    }

    BoxDomainType V = bounds_to_domain(transformations.second);

    Map<RealVariable,Map<RealVariable,RealExpression>> gs;
    for (auto in : inputs.variables()) {
        Map<RealVariable,RealExpression> g;
        for (auto dyn : substituted_dynamics) {
            g[dyn.left_hand_side().base()] = simplify(derivative(dyn.right_hand_side(),in));
        }
        gs[in] = g;
    }

    Map<RealVariable,RealExpression> f_expr;
    RealAssignments subs;
    for (auto in : inputs.variables()) {
        subs.push_back(RealAssignment(in,RealExpression::constant(0)));
    }
    for (auto dyn : substituted_dynamics) {
        f_expr[dyn.left_hand_side().base()] = simplify(substitute(dyn.right_hand_side(),subs));
    }

    RealSpace var_spc(left_hand_sides(substituted_dynamics));

    Vector<RealExpression> f_dyn(var_spc.dimension());
    for (auto var : var_spc.indices()) {
        f_dyn[var.second] = f_expr[var.first];
    }
    ValidatedVectorFunction f = make_function(var_spc,f_dyn);

    Vector<ValidatedVectorFunction> g(gs.size());

    SizeType i = 0;
    for (auto in : inputs.variables()) {
        Vector<RealExpression> g_dyn(var_spc.dimension());
        for (auto var : var_spc.indices()) {
            g_dyn[var.second] = gs[in][var.first];
        }
        g[i++] = ValidatedVectorFunction(make_function(var_spc,g_dyn));
    }

    return make_tuple(f,g,V);
}

struct ScheduledApproximation
{
    SizeType step;
    InputApproximator approximation;

    ScheduledApproximation(SizeType step, InputApproximator approximation) : step(step), approximation(approximation) {}
};

OutputStream& operator<<(OutputStream& os, ScheduledApproximation const& sa) {
    return os << "(" << sa.step << ":" << sa.approximation.kind() << ")"; }

struct ScheduledApproximationComparator
{
    inline bool operator() (const ScheduledApproximation& sa1, const ScheduledApproximation& sa2)
    {
        return (sa1.step > sa2.step);
    }
};

char activity_symbol(SizeType step) {
    switch (step % 4) {
    case 0: return '\\';
    case 1: return '|';
    case 2: return '/';
    default: return '-';
    }
}

BoxDomainType bounds_to_domain(RealVariablesBox const& var_bounds) {
    auto vars = var_bounds.variables();
    List<IntervalDomainType> result;
    for (auto v : vars) {
        result.push_back(cast_exact(widen(IntervalDomainType(var_bounds[v].lower().get_d(),var_bounds[v].upper().get_d()))));
    }
    return Vector<IntervalDomainType>(result);
}

Pair<RealAssignment,RealInterval> centered_variable_transformation(RealVariable const& v, RealInterval const& bounds) {
    if (same(bounds.lower(),-bounds.upper())) return Pair<RealAssignment,RealInterval>(RealAssignment(v,v),bounds);
    else return Pair<RealAssignment,RealInterval>(RealAssignment(v,v+bounds.midpoint()),RealInterval(bounds.lower()-bounds.midpoint(),bounds.upper()-bounds.midpoint()));
}

Pair<RealAssignments,RealVariablesBox> centered_variables_transformation(RealVariablesBox const& inputs) {
    RealAssignments assignments;
    List<RealVariableInterval> new_bounds;
    for (auto entry : inputs.bounds()) {
        auto tr = centered_variable_transformation(entry.first,entry.second);
        assignments.push_back(tr.first);
        new_bounds.push_back(RealVariableInterval(entry.first,tr.second));
    }
    return Pair<RealAssignments,RealVariablesBox>(assignments,new_bounds);
}

FloatDPError get_r(InputApproximation approx_kind) {
    switch (approx_kind) {
    case InputApproximation::AFFINE:
        return FloatDPError(5.0/3u);
    case InputApproximation::SINUSOIDAL:
        return FloatDPError(5.0/4u);
    case InputApproximation::PIECEWISE:
        return FloatDPError(1.3645_upper);
    default:
        ARIADNE_FAIL_MSG("A value of 'r' does not exist for kind " << approx_kind << "\n");
    }
}

Box<UpperIntervalType> apply(VectorFunction<ValidatedTag>const& f, const Box<ExactIntervalType>& bx) {
    return apply(f,Box<UpperIntervalType>(bx));
}

Map<InputApproximation,FloatDP> convert_to_percentages(const Map<InputApproximation,SizeType>& approximation_global_frequencies) {

    SizeType total_steps(0);
    for (auto entry: approximation_global_frequencies) {
        total_steps += entry.second;
    }

    Map<InputApproximation,FloatDP> result;
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


Boolean inputs_are_additive(Vector<ValidatedVectorFunction> const &g) {

    SizeType m = g.size();
    SizeType n = g[0].result_size();

    if (m > n)
        return false;

    for (SizeType j: range(n)) {
        bool foundOne = false;
        for (SizeType i : range(m)) {
            auto eval = g[i][j].evaluate(cast_singleton(ExactBoxType(n,ExactIntervalType(std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()))));
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
compute_norms(DifferentialInclusion const& di, PositiveFloatDPValue const& h, UpperBoxType const& B) {

    auto n = di.f.result_size();
    auto m = di.g.size();
    DoublePrecision pr;
    FloatDPError ze(pr);
    FloatDPError K=ze, pK=ze, L=ze, pL=ze, H=ze, pH=ze;
    Vector<FloatDPError> Kj(n), pKj(n), Lj(n), pLj(n), Hj(n), pHj(n);
    FloatDPUpperBound Lambda=ze;

    auto Df=di.f.differential(cast_singleton(B),2);
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
        auto Dg_i=di.g[i].differential(cast_singleton(B),2);
        FloatDPError Vi(abs(di.V[i]).upper());
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
            FloatDPError Vi(abs(di.V[i]).upper());
            pKj[j] += Vi*pK_matrix[i][j]; pLj[j] += Vi*pL_matrix[i][j]; pHj[j] += Vi*pH_matrix[i][j];
        }
    }

    FloatDPError expLambda = (possibly(Lambda>0)) ? FloatDPError(dexp(Lambda*h)) : FloatDPError(1u,pr);
    FloatDPError expL = cast_positive(exp(L*h));

    return Norms(K,Kj,pK,pKj,L,Lj,pL,pLj,H,Hj,pH,pHj,expLambda,expL);
}

InputApproximator
InputApproximatorFactory::create(DifferentialInclusion const& di, InputApproximation kind, SweeperDP sweeper) const {

    switch(kind) {
    case InputApproximation::ZERO : return InputApproximator(SharedPointer<InputApproximatorInterface>(new ZeroInputApproximator(di,sweeper)));
    case InputApproximation::CONSTANT : return InputApproximator(SharedPointer<InputApproximatorInterface>(new ConstantInputApproximator(di,sweeper)));
    case InputApproximation::AFFINE : return InputApproximator(SharedPointer<InputApproximatorInterface>(new AffineInputApproximator(di,sweeper)));
    case InputApproximation::SINUSOIDAL: return InputApproximator(SharedPointer<InputApproximatorInterface>(new SinusoidalInputApproximator(di,sweeper)));
    case InputApproximation::PIECEWISE : return InputApproximator(SharedPointer<InputApproximatorInterface>(new PiecewiseInputApproximator(di,sweeper)));
    }
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


InclusionIntegrator::InclusionIntegrator(List<InputApproximation> approximations, SweeperDP sweeper, StepSize step_size)
    : _approximations(approximations)
    , _sweeper(sweeper)
    , _step_size(step_size)
    , _number_of_steps_between_simplifications(8)
    , _number_of_variables_to_keep(4)
{
    assert(approximations.size()>0);
}

static const SizeType NUMBER_OF_PICARD_ITERATES=6;

List<ValidatedVectorFunctionModelDP> InclusionIntegrator::flow(DifferentialInclusionIVP const& ivp, Real tmax) {
    ARIADNE_LOG(1,"\n"<<ivp<<"\n");

    const DifferentialInclusion& di = ivp.di;
    const ValidatedVectorFunction& f = di.f;
    const Vector<ValidatedVectorFunction>& g = di.g;
    const BoxDomainType& V = di.V;
    const BoxDomainType& X0 = ivp.X0;

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

    Map<InputApproximation,SizeType> approximation_global_frequencies, approximation_local_frequencies;
    for (auto appro: _approximations) {
        approximation_global_frequencies[appro] = 0;
        approximation_local_frequencies[appro] = 0;
    }

    List<ValidatedVectorFunctionModelDP> result;

    auto step = 0;

    List<ScheduledApproximation> schedule;
    Map<InputApproximation,Nat> delays;

    List<InputApproximator> approximations;
    InputApproximatorFactory factory;
    for (auto appro : _approximations)
        approximations.append(factory.create(di,appro,_sweeper));

    for (auto appro: approximations) {
        schedule.push_back(ScheduledApproximation(SizeType(step),appro));
        delays[appro.kind()] = 0;
    }

    while (possibly(t<FloatDPBounds(tmax,pr))) {

        if (verbosity == 1)
            std::cout << "\r[" << activity_symbol(step) << "] " << static_cast<int>(round(100*t.get_d()/tmax.get_d())) << "%" << std::flush;

        ARIADNE_LOG(2,"step#:"<<step<<", t:"<<t<<", hsug:"<<hsug << "\n");

        List<InputApproximator> approximations_to_use;
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
        SharedPointer<InputApproximator> best;
        FloatDP best_volume(0);

        ARIADNE_LOG(3,"n. of approximations to use="<<approximations_to_use.size()<<"\n");

        SizeType i = 0;
        for (auto i : range(approximations_to_use.size())) {
            this->_approximator = SharedPointer<InputApproximator>(new InputApproximator(approximations_to_use.at(i)));
            ARIADNE_LOG(4,"checking approximation "<<this->_approximator->kind()<<"\n");

            auto current_reach=reach(di,D,evolve_function,B,t,h);
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
                    ARIADNE_LOG(5,"best approximation: " << best->kind() << "\n");
                    best_reach_function = current_reach;
                    best_evolve_function = current_evolve;
                    best_volume = current_volume;
                }
            }
        }

        if (approximations_to_use.size() > 1)
            ARIADNE_LOG(3,"chosen approximation: " << best->kind() << "\n");

        for (auto appro : approximations_to_use) {
            if (best->kind() == appro.kind())
                delays[appro.kind()] = 0;
            else
                delays[appro.kind()]++;

            Nat offset = 1<<delays[appro.kind()];
            schedule.push_back(ScheduledApproximation(step+offset,appro));
        }
        std::sort(schedule.begin(),schedule.end(),ScheduledApproximationComparator());

        ARIADNE_LOG(3,"updated schedule: " << schedule << "\n");

        approximation_global_frequencies[best->kind()] += 1;
        approximation_local_frequencies[best->kind()] += 1;

        reach_function = best_reach_function;
        evolve_function = best_evolve_function;

        if (step%freq==freq-1) {

            double base = 0;
            double rho = 6.0;
            for (auto appro: approximation_local_frequencies) {
                SizeType ppi;
                switch (appro.first) {
                    case InputApproximation::ZERO:
                        ppi = 0;
                        break;
                    case InputApproximation::CONSTANT:
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
                approximation_local_frequencies[appro] = 0;
            }
        }

        evolve_function = this->_reconditioner->expand_errors(evolve_function);

        ARIADNE_LOG(2,"evolve bounds="<<evolve_function.range()<<"\n");

        step+=1;

        t=new_t;
        result.append(reach_function);

    }

    ARIADNE_LOG(1,"\napproximation % ="<<convert_to_percentages(approximation_global_frequencies)<<"\n");

    return result;
}

ValidatedVectorFunctionModelType
InclusionIntegrator::reach(DifferentialInclusion const& di, BoxDomainType D, ValidatedVectorFunctionModelType evolve_function, UpperBoxType B, PositiveFloatDPValue t, PositiveFloatDPValue h) const {

    auto n = di.f.result_size();
    auto m = di.g.size();
    PositiveFloatDPValue new_t=cast_positive(cast_exact((t+h).lower()));

    ValidatedVectorFunctionModelType result;

    if (this->_approximator->kind() != InputApproximation::PIECEWISE) {
        auto e=this->_approximator->compute_errors(h,B);
        ARIADNE_LOG(6,"approximation errors:"<<e<<"\n");
        auto DVh = this->_approximator->build_flow_domain(D,di.V,h);
        ARIADNE_LOG(6,"DVh:"<<DVh<<"\n");
        auto w = this->_approximator->build_w_functions(DVh, n, m);
        ARIADNE_LOG(6,"w:"<<w<<"\n");
        auto dyn_vf = construct_function_affine_in_input(di.f, di.g, w);
        ARIADNE_LOG(6,"dyn VF:"<<dyn_vf<<"\n");
        auto phi = this->compute_flow_function(dyn_vf,DVh,B);
        phi = add_errors(phi,e);
        result=build_reach_function(evolve_function, phi, t, new_t);
    } else {
        auto e=this->_approximator->compute_errors(h,B);
        ARIADNE_LOG(6,"approximation errors:"<<e<<"\n");
        auto DVh_hlf = this->_approximator->build_flow_domain(D,di.V,hlf(h));
        ARIADNE_LOG(6,"DVh_hlf:"<<DVh_hlf<<"\n");
        auto w_hlf = this->_approximator->build_w_functions(DVh_hlf, n, m);
        ARIADNE_LOG(6,"w_hlf:"<<w_hlf<<"\n");
        auto dyn_vf_hlf = construct_function_affine_in_input(di.f, di.g, w_hlf);
        ARIADNE_LOG(6,"dyn_vf_hlf:" << dyn_vf_hlf << "\n");
        auto phi_hlf = this->compute_flow_function(dyn_vf_hlf,DVh_hlf,B);
        PositiveFloatDPValue intermediate_t=cast_positive(cast_exact((t+hlf(h)).lower()));
        auto intermediate_reach=build_reach_function(evolve_function, phi_hlf, t, intermediate_t);
        auto intermediate_evolve=evaluate_evolve_function(intermediate_reach,intermediate_t);

        auto D_int = cast_exact_box(intermediate_evolve.range());

        auto DVh=this->_approximator->build_flow_domain(D_int,di.V,hlf(h));
        auto w =this->_approximator->build_secondhalf_w_functions(DVh, n, m);
        ARIADNE_LOG(6,"w:"<<w<<"\n");
        auto dyn_vf = construct_function_affine_in_input(di.f, di.g, w);
        ARIADNE_LOG(6,"dyn VF:"<<dyn_vf<<"\n");
        auto phi = this->compute_flow_function(dyn_vf,DVh,B);
        phi = add_errors(phi,e);
        result=build_secondhalf_piecewise_reach_function(intermediate_evolve, phi, m, intermediate_t, new_t);
    }
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

ValidatedVectorFunctionModelDP InclusionIntegrator::evaluate_evolve_function(ValidatedVectorFunctionModelDP reach_function, PositiveFloatDPValue t) const {
    return partial_evaluate(reach_function,reach_function.argument_size()-1,t);
}

ValidatedVectorFunctionModelDP add_errors(ValidatedVectorFunctionModelDP phi, Vector<ErrorType> const& e) {
    assert(phi.result_size()==e.size());
    ValidatedVectorTaylorFunctionModelDP& tphi = dynamic_cast<ValidatedVectorTaylorFunctionModelDP&>(phi.reference());
    for (auto i : range(e.size())) {
        tphi[i].add_error(e[i]);
    }
    return phi;
}

ValidatedVectorFunction construct_function_affine_in_input(ValidatedVectorFunction const &f,
                                                           Vector<ValidatedVectorFunction> const &g,
                                                           Vector<ValidatedScalarFunction> const &u) {

    auto n = f.result_size();
    auto m = g.size();
    auto p = u[0].argument_size();
    auto infinity_box = cast_singleton(ExactBoxType(p,ExactIntervalType(std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity())));

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
compute_flow_function(ValidatedVectorFunction const& dyn, BoxDomainType const& domain, UpperBoxType const& B) const {
    auto n=dyn.result_size();
    auto swp=this->_sweeper;

    auto x0f=ValidatedVectorTaylorFunctionModelDP::projection(domain,range(n),swp);
    auto af=ValidatedVectorTaylorFunctionModelDP::projection(domain,range(n,dyn.argument_size()),swp);

    auto picardPhi=ValidatedVectorTaylorFunctionModelDP(n,domain,swp);
    picardPhi=picardPhi+cast_singleton(B);

    for (auto i : range(NUMBER_OF_PICARD_ITERATES)) {
        auto dyn_of_phi = compose(dyn,join(picardPhi,af));
        picardPhi=antiderivative(dyn_of_phi,dyn_of_phi.argument_size()-1)+x0f;
    }

    /*
    TaylorSeriesIntegrator integrator(MaximumError(1e-4),SweepThreshold(1e-8),LipschitzConstant(0.5));

    auto fgws = construct_f_plus_gw_squared(f,g,w);
    auto BVh =_approximator->build_flow_domain(cast_exact_box(B), V, h);

    auto squaredSeriesPhi = integrator.flow_step(fgws,DVh,h,BVh);

    ValidatedVectorTaylorFunctionModelDP& tsquaredSeriesPhi = dynamic_cast<ValidatedVectorTaylorFunctionModelDP&>(squaredSeriesPhi.reference());
    auto seriesPhi=ValidatedVectorTaylorFunctionModelDP(n,squaredSeriesPhi.domain(),swp);
    for (auto i : state_variables) {
        seriesPhi[i] = tsquaredSeriesPhi[i];
    }

    if (volume(picardPhi.range()) < volume(seriesPhi.range())) {
        ARIADNE_LOG(2,"Picard flow function chosen\n");
        return picardPhi;

    } else {
        ARIADNE_LOG(2,"Series flow function chosen\n");
        return seriesPhi;
    }
    */


    return picardPhi;
}


BoxDomainType InputApproximatorBase::build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h) const {
    auto result = D;

    for (Nat i : range(this->_num_params_per_input))
        result = product(result,V);

    return product(result,IntervalDomainType(-h,+h));
}

Vector<ValidatedScalarFunction> InputApproximator::build_secondhalf_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const {
    assert(kind() != InputApproximation::PIECEWISE);
    PiecewiseInputApproximator& approx = dynamic_cast<PiecewiseInputApproximator&>(*this->_impl);
    return approx.build_secondhalf_w_functions(DVh,n,m);
}

Vector<ValidatedScalarFunction> ZeroInputApproximator::build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const {
    auto result = Vector<ValidatedScalarFunction>(m);

    for (auto i : range(0,m))
        result[i] = ValidatedScalarFunction::zero(n+1);

    return result;
}


Vector<ValidatedScalarFunction> ConstantInputApproximator::build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const {
    auto result = Vector<ValidatedScalarFunction>(m);

    for (auto i : range(0,m))
        result[i] = ValidatedScalarFunction::coordinate(n+m+1,n+i);

    return result;
}


Vector<ValidatedScalarFunction> AffineInputApproximator::build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const {
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


Vector<ValidatedScalarFunction> SinusoidalInputApproximator::build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const {

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


Vector<ValidatedScalarFunction> PiecewiseInputApproximator::build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const {
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


Vector<ValidatedScalarFunction> PiecewiseInputApproximator::build_secondhalf_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const {
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
