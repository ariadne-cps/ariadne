/***************************************************************************
 *            verification.cpp
 *
 *  Copyright  2023  Luca Geretti
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

#include "solvers/integrator.hpp"
#include "dynamics/orbit.hpp"
#include "dynamics/vector_field_evolver.hpp"
#include "dynamics/inner_approximation.hpp"
#include "helper/stopwatch.hpp"
#include "helper/container.hpp"
#include "pexplore/task_runner.tpl.hpp"
#include "verification.hpp"

namespace Ariadne {

using std::ostream;

using Helper::List;
using pExplore::ConstrainingState;
using pExplore::ConstraintBuilder;
using pExplore::ConstraintObjectiveImpact;
using pExplore::ConstraintSuccessAction;
using pExplore::ConstraintFailureKind;
using pExplore::TaskInput;
using pExplore::TaskOutput;
using ConstraintIndexType = size_t;

EvaluationSequence::EvaluationSequence(List<TimedVolume> const& tv, Vector<HardConstraintPrescription> const& usages) : _sequence(tv), _prescriptions(usages) { }

size_t EvaluationSequence::number_of_constraints() const { return _prescriptions.size(); }

TimedVolume const& EvaluationSequence::at(size_t idx) const { return _sequence.at(idx); }

double const& EvaluationSequence::near(double time) const {
    size_t lower = 0;
    size_t upper = size()-1;

    if (_sequence.at(lower).time >= time) return _sequence.at(lower).volume;
    if (_sequence.at(upper).time <= time) return _sequence.at(upper).volume;

    while (upper - lower > 1) {
        auto mid = (lower+upper)/2;
        if (_sequence.at(mid).time < time) lower = mid;
        else upper = mid;
    }
    if (time - _sequence.at(lower).time > _sequence.at(upper).time - time) return _sequence.at(upper).volume;
    else return _sequence.at(lower).volume;
}

HardConstraintPrescription const& EvaluationSequence::usage(ConstraintIndexType idx) const { return _prescriptions.at(idx); }

size_t EvaluationSequence::size() const { return _sequence.size(); }

ostream& operator<<(ostream& os, EvaluationSequence const& es) {
    os << "timed_volumes:{";
    for (size_t i=0; i<es.size()-1; ++i) os << es.at(i) << ",";
    return os << es.at(es.size()-1) << "}, usages:" << es._prescriptions;
}

EvaluationSequenceBuilder::EvaluationSequenceBuilder(size_t N, EffectiveVectorMultivariateFunction const& h) : _N(N), _M(h.result_size()), _h(h), _prescriptions(h.result_size(),SatisfactionPrescription::TRUE),
                                          _max_robustness_false(h.result_size(),{{SatisfactionPrescription::FALSE_FOR_ALL,0.0},{SatisfactionPrescription::FALSE_FOR_SOME,0.0}}),
                                          _t_false(h.result_size(),{{SatisfactionPrescription::FALSE_FOR_ALL,0.0},{SatisfactionPrescription::FALSE_FOR_SOME,0.0}}) { }

void EvaluationSequenceBuilder::add_from(LabelledEnclosure const& e) {
    auto eval = evaluate_from_function(_h,e);

    auto bb = e.euclidean_set().bounding_box();
    add({e.time_function().range().midpoint().get_d(),reinterpret_cast<Vector<FloatDPBounds>const&>(bb),eval});
}

void EvaluationSequenceBuilder::add(TimedBoxEvaluation const& tbe) {
    HELPER_PRECONDITION(_timed_box_evaluations.empty() or _timed_box_evaluations.at(_timed_box_evaluations.size()-1).time < tbe.time)
    _timed_box_evaluations.push_back(tbe);
    for (ConstraintIndexType m=0; m<_M; ++m) {
        auto const& evaluation = tbe.evaluation[m];
        auto upper = evaluation.upper().get_d();
        auto lower = evaluation.lower().get_d();
        if (_prescriptions[m] != SatisfactionPrescription::FALSE_FOR_ALL) {
            if (upper < 0) _prescriptions[m] = SatisfactionPrescription::FALSE_FOR_ALL;
            else if (lower < 0) _prescriptions[m] = SatisfactionPrescription::FALSE_FOR_SOME;
        }
        SatisfactionPrescription prescription_to_check = SatisfactionPrescription::TRUE;
        if (upper < 0) prescription_to_check = SatisfactionPrescription::FALSE_FOR_ALL;
        else if (lower < 0) prescription_to_check = SatisfactionPrescription::FALSE_FOR_SOME;

        if (prescription_to_check != SatisfactionPrescription::TRUE) {
            auto chi = get_chi(tbe.box,_h[m],prescription_to_check);
            auto robustness = get_rho(chi,_N,volume(tbe.box),prescription_to_check);
            if (_max_robustness_false[m].get(prescription_to_check) < robustness) {
                _max_robustness_false[m].at(prescription_to_check) = robustness;
                _t_false[m].at(prescription_to_check) = tbe.time;
            }
        }
    }
}

EvaluationSequence EvaluationSequenceBuilder::build() const {
    HELPER_PRECONDITION(not _timed_box_evaluations.empty())

    Vector<double> T_star(_M, 0.0);
    for (ConstraintIndexType m=0; m<_M; ++m) {
        switch (_prescriptions[m]) {
            case SatisfactionPrescription::TRUE :           T_star[m] = _timed_box_evaluations.at(_timed_box_evaluations.size()-1).time; break;
            case SatisfactionPrescription::FALSE_FOR_ALL :  T_star[m] = _t_false[m].get(SatisfactionPrescription::FALSE_FOR_ALL); break;
            case SatisfactionPrescription::FALSE_FOR_SOME : T_star[m] = _t_false[m].get(SatisfactionPrescription::FALSE_FOR_SOME); break;
            default: HELPER_FAIL_MSG("Unhandled SatisfactionPrescription value")
        }
    }

    List<TimedVolume> timed_volumes;

    double v0 = volume(_timed_box_evaluations.at(0).box);
    timed_volumes.push_back({_timed_box_evaluations.at(0).time,v0});

    Vector<double> alpha(_M, std::numeric_limits<double>::max());
    for (size_t i=1; i < _timed_box_evaluations.size(); ++i) {
        auto const& tbe = _timed_box_evaluations.at(i);
        double v = volume(tbe.box);
        for (ConstraintIndexType m=0; m<_M; ++m) {
            if (tbe.time <= T_star[m]) {
                if (_prescriptions[m] == SatisfactionPrescription::TRUE or tbe.time == T_star[m]) {
                    double chi = get_chi(tbe.box,_h[m],_prescriptions[m]);
                    double rho = get_rho(_N,chi,v,_prescriptions[m]);
                    auto this_alpha = rho/(v-v0);
                    CONCLOG_PRINTLN_AT(1,"i="<< i <<",m="<< m << ",chi="<< chi << ",rho="<<rho<<",v-v0=" << v-v0<<",alpha="<<this_alpha)
                    if (this_alpha>0) alpha[m] = std::min(alpha[m],this_alpha);
                }
            }
        }
        timed_volumes.push_back({tbe.time,v});
    }

    Vector<HardConstraintPrescription> constants(_M, {SatisfactionPrescription::TRUE, 0.0, 0.0});
    for (ConstraintIndexType m=0; m<_M; ++m)
        constants[m] = {_prescriptions[m],T_star[m],alpha[m]};

    return {timed_volumes,constants};
}

size_t EvaluationSequenceBuilder::_list_index(double time) const {
    auto initial_time = _timed_box_evaluations.at(0).time;
    auto final_time = _timed_box_evaluations.at(_timed_box_evaluations.size()-1).time;
    HELPER_PRECONDITION(initial_time <= time)
    HELPER_PRECONDITION(final_time >= time)

    size_t lower = 0;
    size_t upper = _timed_box_evaluations.size()-1;

    if (time == initial_time) return 0;
    if (time == final_time) return _timed_box_evaluations.size()-1;

    while (upper - lower > 0) {
        auto mid = (lower+upper)/2;
        auto midtime = _timed_box_evaluations.at(mid).time;
        if (midtime == time) return mid;
        else if (midtime < time) lower = mid;
        else upper = mid;
    }
    HELPER_ASSERT_MSG(_timed_box_evaluations.at(lower).time == time, "The time " << time << " could not be found in the list")
    return lower;
}


Vector<FloatDPBounds> resize(Vector<FloatDPBounds> const& bx, double chi, bool widen) {
    double factor = (widen ? chi - 1 : 1 - chi);
    Vector<FloatDPBounds> result(bx.size(),DoublePrecision());
    for(size_t i=0; i!=bx.size(); ++i) {
        auto width = bx[i].upper()-bx[i].lower();
        auto lower = cast_exact((bx[i].lower()-width*factor).get_d());
        auto upper = cast_exact((bx[i].upper()+width*factor).get_d());
        result[i] = FloatDPBounds(lower,upper,DoublePrecision());
    }
    return result;
}

double volume(Vector<FloatDPBounds> const& bnds) {
    double result = 1;
    for (size_t i=0; i<bnds.size(); ++i) result *= (bnds[i].upper()-bnds[i].lower()).get_d();
    return result;
}

Vector<FloatDPBounds> widen(Vector<FloatDPBounds> const& bx, double chi) { return resize(bx,chi,true); }
Vector<FloatDPBounds> shrink(Vector<FloatDPBounds> const& bx, double chi) { return resize(bx,chi,false); }

void search_chi_for_true(double& lower, double& upper, Vector<FloatDPBounds>const& bnds, EffectiveScalarMultivariateFunction const& constraint, double factor) {
    while (true) {
        if (constraint(widen(bnds, upper)).lower().get_d() <= 0) { lower = upper/2; break; }
        else upper*=2;
    }
    double margin = (upper-lower)/factor;
    while (upper-lower > margin) {
        auto midpoint = (upper+lower)/2;
        if (constraint(widen(bnds, midpoint)).lower().get_d() <= 0) upper = midpoint;
        else lower = midpoint;
    }
}

void search_chi_for_false_all(double& lower, double& upper, Vector<FloatDPBounds>const& bnds, EffectiveScalarMultivariateFunction const& constraint, double factor) {
    while (true) {
        if (constraint(widen(bnds, upper)).upper().get_d() >= 0) { lower = upper/2; break; }
        else upper*=2;
    }
    double margin = (upper-lower)/factor;
    while (upper-lower > margin) {
        auto midpoint = (upper+lower)/2;
        if (constraint(widen(bnds, midpoint)).upper().get_d() >= 0) upper = midpoint;
        else lower = midpoint;
    }
}

void search_chi_for_false_some(double& lower, double& upper, Vector<FloatDPBounds>const& bnds, EffectiveScalarMultivariateFunction const& constraint, double factor) {
    while (true) {
        auto eval = constraint(shrink(bnds, upper));
        if (eval.upper().get_d() < 0) { lower = upper = std::numeric_limits<double>::max(); break; }
        if (eval.lower().get_d() >= 0) { lower = upper/2; break; }
        else upper*=2;
    }
    double margin = (upper-lower)/factor;
    while (upper-lower > margin) {
        auto midpoint = (upper+lower)/2;
        if (constraint(shrink(bnds, midpoint)).lower().get_d() >= 0) upper = midpoint;
        else lower = midpoint;
    }
}

double get_chi(Vector<FloatDPBounds>const& bnds, EffectiveScalarMultivariateFunction const& constraint, SatisfactionPrescription prescription) {
    static const double FACTOR = 100;
    double lower = 1.0;
    double upper = 1.0;
    switch (prescription) {
        case SatisfactionPrescription::TRUE : HELPER_PRECONDITION(constraint(bnds).lower().get_d() > 0) search_chi_for_true(lower, upper, bnds, constraint, FACTOR); break;
        case SatisfactionPrescription::FALSE_FOR_ALL : HELPER_PRECONDITION(constraint(bnds).upper().get_d() < 0) search_chi_for_false_all(lower, upper, bnds, constraint, FACTOR); break;
        case SatisfactionPrescription::FALSE_FOR_SOME : HELPER_PRECONDITION(constraint(bnds).lower().get_d() < 0) search_chi_for_false_some(lower, upper, bnds, constraint, FACTOR); break;
        default: HELPER_FAIL_MSG("Unhandled SatisfactionPrescription value")
    }
    return lower;
}

double get_rho(double chi, size_t N, double volume, SatisfactionPrescription prescription) {
    if (prescription == SatisfactionPrescription::FALSE_FOR_SOME) {
        auto pw = pow(chi,N);
        return (pw == std::numeric_limits<double>::infinity() ? volume : volume*(1.0-1.0/pw));

    } else return (pow(chi,N)-1)*volume;
}

EffectiveVectorMultivariateFunction to_function(List<RealExpression> const& constraints, RealSpace const& spc) {
    EffectiveVectorMultivariateFunction result(constraints.size(),spc.size());
    for (ConstraintIndexType m=0; m<constraints.size(); ++m)
        result[m] = make_function(spc,constraints.at(m));
    return result;
}

FloatDPBounds evaluate_from_function(EffectiveScalarMultivariateFunction const& function, LabelledEnclosure const& enclosure) {
    auto bb = enclosure.euclidean_set().bounding_box();
    return function(reinterpret_cast<Vector<FloatDPBounds>const&>(bb));
}

Vector<FloatDPBounds> evaluate_from_function(EffectiveVectorMultivariateFunction const& function, LabelledEnclosure const& enclosure) {
    auto bb = enclosure.euclidean_set().bounding_box();
    return function(reinterpret_cast<Vector<FloatDPBounds>const&>(bb));
}

EvaluationSequence evaluate_approximate_orbit(Orbit<LabelledEnclosure> const& orbit, EffectiveVectorMultivariateFunction const& h) {
    CONCLOG_SCOPE_CREATE

    EvaluationSequenceBuilder sb(orbit.initial().dimension(),h);
    sb.add_from(orbit.initial());
    for (auto const& e : orbit.intermediate()) { sb.add_from(e); }
    return sb.build();
}

List<pExplore::Constraint<VectorFieldEvolver>> build_task_constraints(EvaluationSequence const& evaluation, EffectiveVectorMultivariateFunction const& h) {
    using A = VectorFieldEvolver;
    using I = TaskInput<A>;
    using O = TaskOutput<A>;

    List<pExplore::Constraint<A>> result;

    for (ConstraintIndexType m=0; m<evaluation.number_of_constraints(); ++m) {
        result.push_back(ConstraintBuilder<A>([evaluation,m](I const&, O const& o) {
            auto v0 = evaluation.at(0).volume;
            auto v_approx = evaluation.near(o.time.get_d());
            auto alpha = evaluation.usage(m).alpha;
            return (alpha+1)*v_approx - alpha*v0 - o.evolve.euclidean_set().bounding_box().volume().get_d();
        }).set_name("objective&soft#"+to_string(m)).set_group_id(m).set_objective_impact(ConstraintObjectiveImpact::UNSIGNED).set_failure_kind(ConstraintFailureKind::SOFT).build()
        );
        result.push_back(ConstraintBuilder<A>([evaluation,h,m](I const&, O const& o) {
            if (evaluation.usage(m).sigma == SatisfactionPrescription::TRUE)
                return evaluate_from_function(h[m],o.reach).lower().get_d();
            else return evaluation.usage(m).t_star - o.time.get_d();
        }).set_name("hard#"+to_string(m)).set_group_id(m).set_failure_kind(ConstraintFailureKind::HARD).build()
        );
        if (evaluation.usage(m).sigma != SatisfactionPrescription::TRUE) {
            result.push_back(ConstraintBuilder<A>([evaluation,h,m](I const&, O const& o) {
                if (evaluation.usage(m).sigma == SatisfactionPrescription::FALSE_FOR_ALL) {
                    return -evaluate_from_function(h[m], o.evolve).upper().get_d();
                } else {
                    if (evaluate_from_function(h[m],o.evolve).lower().get_d() < 0) {
                        try {
                            auto approximator = NonlinearCandidateValidationInnerApproximator(ParallelLinearisationContractor(GLPKSimplex(),2,1));
                            auto inner_evolve = approximator.compute_from(o.evolve);
                            return -evaluate_from_function(h[m],inner_evolve).lower().get_d();
                        } catch (std::exception&) { }
                    }
                    return -1.0;
                }
            }).set_name("falsify/success#"+to_string(m)).set_group_id(m).set_success_action(ConstraintSuccessAction::DEACTIVATE).build()
            );
        }
    }

    return result;
}

Vector<Kleenean> synthesise_outcomes(EvaluationSequence const& preanalysis, ConstrainingState<VectorFieldEvolver> const& constraining) {
    auto M = preanalysis.number_of_constraints();
    Vector<Kleenean> result(M,indeterminate);
    for (auto const& state : constraining.states()) {
        auto const& sigma = preanalysis.usage(state.constraint().group_id()).sigma;
        if (sigma == SatisfactionPrescription::TRUE and state.constraint().failure_kind() == ConstraintFailureKind::HARD and not state.has_failed())
            result.at(state.constraint().group_id()) = true;
        else if (sigma != SatisfactionPrescription::TRUE and state.constraint().success_action() == ConstraintSuccessAction::DEACTIVATE and state.has_succeeded()) {
            result.at(state.constraint().group_id()) = false;
        }
    }
    return result;
}

Tuple<Orbit<LabelledEnclosure>,Orbit<LabelledEnclosure>,Vector<Kleenean>> constrained_evolution(VectorField const& dynamics, RealExpressionBoundedConstraintSet const& initial_set, Real const& evolution_time,
                                                                      List<RealExpression> const& constraints, Configuration<VectorFieldEvolver> const& configuration) {
    CONCLOG_SCOPE_CREATE

    VectorFieldEvolver approximate_evolver(dynamics,Configuration<VectorFieldEvolver>().set_integrator(AffineIntegrator(Configuration<AffineIntegrator>())));

    Helper::Stopwatch<std::chrono::microseconds> sw;
    CONCLOG_PRINTLN("Computing approximate evolution...")
    auto approximate_orbit = approximate_evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    sw.click();
    CONCLOG_PRINTLN_AT(0,"Done in " << sw.elapsed_seconds() << " seconds.")

    sw.restart();
    CONCLOG_PRINTLN("Processing approximate evolution for constraints... ")

    auto h = to_function(constraints, dynamics.state_space());
    CONCLOG_PRINTLN_VAR(h)

    auto analysis = evaluate_approximate_orbit(approximate_orbit, h);

    auto task_constraints = build_task_constraints(analysis, h);

    sw.click();
    CONCLOG_PRINTLN_AT(0,"Done in " << sw.elapsed_seconds() << " seconds.")

    CONCLOG_PRINTLN_VAR_AT(1,analysis)

    auto integr = dynamic_cast<TaylorPicardIntegrator const&>(configuration.integrator());
    CONCLOG_PRINTLN_VAR_AT(1,integr.configuration())
    VectorFieldEvolver rigorous_evolver(dynamics,configuration);
    sw.restart();

    CONCLOG_PRINTLN("Computing rigorous evolution... ")
    rigorous_evolver.set_constraints(task_constraints);
    auto rigorous_orbit = rigorous_evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    sw.click();
    CONCLOG_PRINTLN_AT(0,"Done in " << sw.elapsed_seconds() << " seconds.")

    auto constraining_state = rigorous_evolver.constraining_state();

    for (auto state : constraining_state.states()) {
        CONCLOG_PRINTLN_AT(1,state)
    }

    auto outcomes = synthesise_outcomes(analysis,constraining_state);

    return {approximate_orbit,rigorous_orbit,outcomes};
}

} // namespace Ariadne
