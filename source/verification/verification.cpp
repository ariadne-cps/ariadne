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
#include "pexplore/robustness_controller.hpp"
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
using pExplore::TimeProgressLinearRobustnessController;
using ConstraintIndexType = size_t;

double nthroot(double value, size_t n) {
    static const size_t NUM_ITERATIONS = 10;
    auto number_of_sqrt = static_cast<size_t>(std::round(sqrt(n)));
    double result = value;
    for (size_t i=0; i<number_of_sqrt; ++i)
        result = sqrt(result);

    for (size_t i=0; i<NUM_ITERATIONS; ++i) {
        result = result - (pow(result,n) - value)/(n*pow(result,n-1));
    }
    return result;
}

EvaluationSequence::EvaluationSequence(List<TimedMeasurement> const& tv, Vector<ControlSpecification> const& usages) : _sequence(tv), _prescriptions(usages) { }

size_t EvaluationSequence::number_of_constraints() const { return _prescriptions.size(); }

TimedMeasurement const& EvaluationSequence::at(size_t idx) const { return _sequence.at(idx); }

TimedMeasurement const& EvaluationSequence::near(double time) const {
    size_t lower = 0;
    size_t upper = size()-1;

    if (_sequence.at(lower).time >= time) return _sequence.at(lower);
    if (_sequence.at(upper).time <= time) return _sequence.at(upper);

    while (upper - lower > 1) {
        auto mid = (lower+upper)/2;
        if (_sequence.at(mid).time < time) lower = mid;
        else upper = mid;
    }
    if (time - _sequence.at(lower).time > _sequence.at(upper).time - time) return _sequence.at(upper);
    else return _sequence.at(lower);
}

ControlSpecification const& EvaluationSequence::usage(ConstraintIndexType idx) const { return _prescriptions.at(idx); }

size_t EvaluationSequence::size() const { return _sequence.size(); }

ostream& operator<<(ostream& os, EvaluationSequence const& es) {
    os << "timed_beta_B:{";
    for (size_t i=0; i<es.size()-1; ++i) os << es.at(i) << ", ";
    return os << es.at(es.size()-1) << "}, usages:" << es._prescriptions;
}

EvaluationSequenceBuilder::EvaluationSequenceBuilder(size_t N, EffectiveVectorMultivariateFunction const& h) : _N(N), _M(h.result_size()), _h(h), _prescriptions(h.result_size(), SatisfactionPrescriptionKind::TRUE),
                                          _max_robustness_false(h.result_size(),{{SatisfactionPrescriptionKind::FALSE_FOR_ALL, 0.0}, {SatisfactionPrescriptionKind::FALSE_FOR_SOME, 0.0}}),
                                          _t_false(h.result_size(),{{SatisfactionPrescriptionKind::FALSE_FOR_ALL, 0.0}, {SatisfactionPrescriptionKind::FALSE_FOR_SOME, 0.0}}) { }

void EvaluationSequenceBuilder::add_from(LabelledEnclosure const& approximate, LabelledEnclosure const& rigorous) {
    auto approximate_evaluation = evaluate_from_function(_h,approximate);

    auto approximate_box = approximate.euclidean_set().bounding_box();
    auto rigorous_box = rigorous.euclidean_set().bounding_box();
    add({approximate.time_function().range().midpoint().get_d(),reinterpret_cast<Vector<FloatDPBounds>const&>(approximate_box),reinterpret_cast<Vector<FloatDPBounds>const&>(rigorous_box),approximate_evaluation});
}

void EvaluationSequenceBuilder::add(TimedBoxEvaluation const& tbe) {
    HELPER_PRECONDITION(_timed_box_evaluations.empty() or _timed_box_evaluations.at(_timed_box_evaluations.size()-1).time < tbe.time)
    _timed_box_evaluations.push_back(tbe);
    for (ConstraintIndexType m=0; m<_M; ++m) {
        auto const& approximate_evaluation = tbe.approximate_evaluation[m];
        auto upper = approximate_evaluation.upper().get_d();
        auto lower = approximate_evaluation.lower().get_d();
        if (_prescriptions[m] != SatisfactionPrescriptionKind::FALSE_FOR_ALL) {
            if (upper < 0) _prescriptions[m] = SatisfactionPrescriptionKind::FALSE_FOR_ALL;
            else if (lower < 0) _prescriptions[m] = SatisfactionPrescriptionKind::FALSE_FOR_SOME;
        }
        SatisfactionPrescriptionKind prescription_to_check = SatisfactionPrescriptionKind::TRUE;
        if (upper < 0) prescription_to_check = SatisfactionPrescriptionKind::FALSE_FOR_ALL;
        else if (lower < 0) prescription_to_check = SatisfactionPrescriptionKind::FALSE_FOR_SOME;

        if (prescription_to_check != SatisfactionPrescriptionKind::TRUE) {
            auto chi = get_chi(tbe.approximate_box,_h[m],prescription_to_check);
            auto robustness = get_rho(chi,get_beta(tbe.approximate_box,_N),prescription_to_check);
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
            case SatisfactionPrescriptionKind::TRUE : T_star[m] = _timed_box_evaluations.at(_timed_box_evaluations.size() - 1).time; break;
            case SatisfactionPrescriptionKind::FALSE_FOR_ALL : T_star[m] = _t_false[m].get(SatisfactionPrescriptionKind::FALSE_FOR_ALL); break;
            case SatisfactionPrescriptionKind::FALSE_FOR_SOME : T_star[m] = _t_false[m].get(SatisfactionPrescriptionKind::FALSE_FOR_SOME); break;
            default: HELPER_FAIL_MSG("Unhandled SatisfactionPrescriptionKind value")
        }
    }

    List<TimedMeasurement> timed_measurements;

    double approximate_beta0 = get_beta(_timed_box_evaluations.at(0).approximate_box,_N);
    double rigorous_beta0 = get_beta(_timed_box_evaluations.at(0).rigorous_box,_N);
    timed_measurements.push_back({_timed_box_evaluations.at(0).time,approximate_beta0,rigorous_beta0});

    Vector<double> alpha(_M, std::numeric_limits<double>::max());
    for (size_t i=1; i < _timed_box_evaluations.size(); ++i) {
        auto const& tbe = _timed_box_evaluations.at(i);
        double approximate_beta = get_beta(tbe.approximate_box,_N);
        double rigorous_beta = get_beta(tbe.rigorous_box,_N);
        for (ConstraintIndexType m=0; m<_M; ++m) {
            if (tbe.time <= T_star[m]) {
                if (_prescriptions[m] == SatisfactionPrescriptionKind::TRUE or tbe.time == T_star[m]) {
                    double chi = get_chi(tbe.approximate_box,_h[m],_prescriptions[m]);
                    double rho = get_rho(chi,approximate_beta,_prescriptions[m]);
                    auto this_alpha = rho/(rigorous_beta-approximate_beta);
                    CONCLOG_PRINTLN_AT(1,"i="<< i <<",m="<< m << ",chi="<< chi << ",rho="<<rho<<",e="<< rigorous_beta-approximate_beta << ",alpha="<<this_alpha)
                    if (this_alpha>0) alpha[m] = std::min(alpha[m],this_alpha);
                }
            }
        }
        timed_measurements.push_back({tbe.time,approximate_beta,rigorous_beta});
    }

    Vector<ControlSpecification> constants(_M, {SatisfactionPrescriptionKind::TRUE, 0.0, 0.0});
    for (ConstraintIndexType m=0; m<_M; ++m)
        constants[m] = {_prescriptions[m],T_star[m],alpha[m]};

    return {timed_measurements,constants};
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

ConstraintSatisfaction::ConstraintSatisfaction(List<RealExpression> const& cs, RealSpace const& spc) :
    _cs(cs), _space(spc), _prescriptions(Vector<SatisfactionPrescriptionKind>(cs.size(),SatisfactionPrescriptionKind::UNSPECIFIED)), _outcomes(Vector<LogicalValue>(cs.size(),LogicalValue::INDETERMINATE)) { }

size_t ConstraintSatisfaction::dimension() const {
    return _cs.size();
}

EffectiveVectorMultivariateFunction ConstraintSatisfaction::indeterminate_constraints_function() const {
    List<EffectiveScalarMultivariateFunction> result;
    for (ConstraintIndexType m=0; m<_cs.size(); ++m)
        if (Detail::is_indeterminate(_outcomes[m]))
            result.push_back(make_function(_space,_cs.at(m)));
    HELPER_ASSERT(not result.empty())
    return result;
}

void ConstraintSatisfaction::merge_from_uncontrolled(ConstrainingState<VectorFieldEvolver> const& state) {
    auto indeterminates = indeterminate_indexes();
    for (auto const& s : state.states()) {
        auto const& m = indeterminates.at(s.constraint().group_id());
        if (s.constraint().failure_kind() == ConstraintFailureKind::HARD and not s.has_failed())
            set_outcome(m,true);
        else if (s.constraint().success_action() == ConstraintSuccessAction::DEACTIVATE and s.has_succeeded())
            set_outcome(m,false);
    }
}

void ConstraintSatisfaction::merge_from_controlled(ConstrainingState<VectorFieldEvolver> const& state, EvaluationSequence const& preanalysis) {
    auto indeterminates = indeterminate_indexes();
    for (auto const& s : state.states()) {
        auto const& sigma = preanalysis.usage(s.constraint().group_id()).sigma;
        auto const& m = indeterminates.at(s.constraint().group_id());
        if (sigma == SatisfactionPrescriptionKind::TRUE and s.constraint().failure_kind() == ConstraintFailureKind::HARD and not s.has_failed())
            set_outcome(m,true);
        else if (sigma != SatisfactionPrescriptionKind::TRUE and s.constraint().success_action() == ConstraintSuccessAction::DEACTIVATE and s.has_succeeded())
            set_outcome(m,false);
        reset_prescription(m,sigma);
    }
}

List<size_t> ConstraintSatisfaction::indeterminate_indexes() const {
    List<size_t> result;
    for (size_t m=0; m<dimension(); ++m)
        if (is_indeterminate(_outcomes[m]))
            result.push_back(m);
    return result;
}

bool ConstraintSatisfaction::completed() const {
    return indeterminate_indexes().size() == 0;
}

RealExpression const& ConstraintSatisfaction::expression(size_t m) const {
    HELPER_PRECONDITION(m<dimension())
    return _cs[m];
}

SatisfactionPrescriptionKind const& ConstraintSatisfaction::prescription(size_t m) const {
    HELPER_PRECONDITION(m<dimension())
    return _prescriptions[m];
}

LogicalValue const& ConstraintSatisfaction::outcome(size_t m) const {
    HELPER_PRECONDITION(m<dimension())
    return _outcomes[m];
}

ostream& operator<<(ostream& os, ConstraintSatisfaction const& cs) {
    os << "{";
    for (size_t m=0; m<cs.dimension()-1; ++m) {
        os << cs.expression(m) << " >= 0 : " << cs.outcome(m);
        if (cs.prescription(m) != SatisfactionPrescriptionKind::TRUE or Detail::is_indeterminate(cs.outcome(m))) {
            os << " (" << cs.prescription(m) << ")";
        }
        os << ", ";
    }
    if (cs.dimension() > 0) {
        os << cs.expression(cs.dimension()-1) << " >= 0 : " << cs.outcome(cs.dimension()-1);
        if (cs.prescription(cs.dimension()-1) != SatisfactionPrescriptionKind::TRUE or Detail::is_indeterminate(cs.outcome(cs.dimension()-1))) {
            os << " (" << cs.prescription(cs.dimension()-1) << ")";
        }
    }
    return os << "}";
}

void ConstraintSatisfaction::reset_prescription(size_t m, SatisfactionPrescriptionKind prescription) {
    HELPER_PRECONDITION(m<dimension())
    _prescriptions[m] = prescription;
}

void ConstraintSatisfaction::set_outcome(size_t m, bool outcome) {
    HELPER_PRECONDITION(m<dimension())
    HELPER_PRECONDITION(is_indeterminate(_outcomes[m]))
    _outcomes[m] = Detail::make_logical_value(outcome);
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

double get_beta(Vector<FloatDPBounds> const& bnds, size_t n) {
    double result = 1;
    for (size_t i=0; i<bnds.size(); ++i) result *= (bnds[i].upper()-bnds[i].lower()).get_d();
    return nthroot(result,n);
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

double get_chi(Vector<FloatDPBounds>const& bnds, EffectiveScalarMultivariateFunction const& constraint, SatisfactionPrescriptionKind prescription) {
    static const double FACTOR = 100;
    double lower = 1.0;
    double upper = 1.0;
    switch (prescription) {
        case SatisfactionPrescriptionKind::TRUE : HELPER_PRECONDITION(constraint(bnds).lower().get_d() > 0) search_chi_for_true(lower, upper, bnds, constraint, FACTOR); break;
        case SatisfactionPrescriptionKind::FALSE_FOR_ALL : HELPER_PRECONDITION(constraint(bnds).upper().get_d() < 0) search_chi_for_false_all(lower, upper, bnds, constraint, FACTOR); break;
        case SatisfactionPrescriptionKind::FALSE_FOR_SOME : HELPER_PRECONDITION(constraint(bnds).lower().get_d() < 0) search_chi_for_false_some(lower, upper, bnds, constraint, FACTOR); break;
        default: HELPER_FAIL_MSG("Unhandled SatisfactionPrescriptionKind value")
    }
    return lower;
}

double get_rho(double chi, double beta, SatisfactionPrescriptionKind prescription) {
    if (prescription == SatisfactionPrescriptionKind::FALSE_FOR_SOME) {
        return (chi == std::numeric_limits<double>::infinity() ? beta : beta*(1.0-1.0/chi));
    } else return (chi-1.0)*beta;
}

FloatDPBounds evaluate_from_function(EffectiveScalarMultivariateFunction const& function, LabelledEnclosure const& enclosure) {
    auto bb = enclosure.euclidean_set().bounding_box();
    return function(reinterpret_cast<Vector<FloatDPBounds>const&>(bb));
}

Vector<FloatDPBounds> evaluate_from_function(EffectiveVectorMultivariateFunction const& function, LabelledEnclosure const& enclosure) {
    auto bb = enclosure.euclidean_set().bounding_box();
    return function(reinterpret_cast<Vector<FloatDPBounds>const&>(bb));
}

EvaluationSequence evaluate_singleton_orbits(Orbit<LabelledEnclosure> const& approximate, Orbit<LabelledEnclosure> const& rigorous, EffectiveVectorMultivariateFunction const& h) {
    CONCLOG_SCOPE_CREATE

    EvaluationSequenceBuilder sb(approximate.initial().dimension(),h);
    sb.add_from(approximate.initial(),rigorous.initial());
    size_t r_idx = 0;
    for (auto const& a : approximate.intermediate()) {
        auto a_t = a.time_function().range().midpoint().get_d();
        LabelledEnclosure r = rigorous.intermediate()[r_idx];
        while (r_idx < rigorous.intermediate().size()-1) {
            auto r_t = r.time_function().range().midpoint().get_d();
            LabelledEnclosure r_next = rigorous.intermediate()[r_idx+1];
            auto r_t_next = r_next.time_function().range().midpoint().get_d();
            if (abs(r_t_next - a_t) < abs(r_t - a_t)) {
                ++r_idx;
                r = r_next;
            } else break;
        }
        sb.add_from(a,r);
    }
    return sb.build();
}

List<pExplore::Constraint<VectorFieldEvolver>> build_uncontrolled_task_constraints(EffectiveVectorMultivariateFunction const& h) {
    using A = VectorFieldEvolver;
    using I = TaskInput<A>;
    using O = TaskOutput<A>;

    List<pExplore::Constraint<A>> result;

    for (ConstraintIndexType m=0; m<h.result_size(); ++m) {
        result.push_back(ConstraintBuilder<A>([h,m](I const&, O const& o) {
            return evaluate_from_function(h[m],o.reach).lower().get_d();
        }).set_name("true#"+to_string(m)).set_group_id(m).set_failure_kind(ConstraintFailureKind::HARD).build()
        );
        result.push_back(ConstraintBuilder<A>([h,m](I const&, O const& o) {
            return -evaluate_from_function(h[m], o.evolve).upper().get_d();
        }).set_name("false#"+to_string(m)).set_group_id(m).set_success_action(ConstraintSuccessAction::DEACTIVATE).build()
        );
    }

    return result;
}

List<pExplore::Constraint<VectorFieldEvolver>> build_controlled_task_constraints(EffectiveVectorMultivariateFunction const& h, EvaluationSequence const& evaluation) {
    using A = VectorFieldEvolver;
    using I = TaskInput<A>;
    using O = TaskOutput<A>;

    List<pExplore::Constraint<A>> result;

    for (ConstraintIndexType m=0; m<evaluation.number_of_constraints(); ++m) {
        result.push_back(ConstraintBuilder<A>([evaluation,m](I const&, O const& o) {
            auto t = o.time.get_d();
            auto near_evaluation = evaluation.near(t);
            auto const& approximate_beta = near_evaluation.approximate_beta;
            auto const& rigorous_beta = near_evaluation.rigorous_beta;
            auto alpha = evaluation.usage(m).alpha;
            return (alpha+1.0)*approximate_beta - alpha*rigorous_beta - nthroot(o.evolve.euclidean_set().bounding_box().volume().get_d(),o.evolve.dimension());
        }).set_name("objective&soft#"+to_string(m)).set_group_id(m).set_objective_impact(ConstraintObjectiveImpact::UNSIGNED).set_failure_kind(ConstraintFailureKind::SOFT)
          .set_controller(TimeProgressLinearRobustnessController<VectorFieldEvolver>([](I const&, O const& o){ return o.time.get_d(); },evaluation.usage(m).t_star))                                                                                                                                                                                                                                                                                                                                                          .set_controller(TimeProgressLinearRobustnessController<VectorFieldEvolver>([](I const&, O const& o){ return o.time.get_d(); },evaluation.usage(m).t_star))
          .build()
        );
        result.push_back(ConstraintBuilder<A>([evaluation,h,m](I const&, O const& o) {
            if (evaluation.usage(m).sigma == SatisfactionPrescriptionKind::TRUE)
                return evaluate_from_function(h[m],o.reach).lower().get_d();
            else return evaluation.usage(m).t_star - o.time.get_d();
        }).set_name("hard#"+to_string(m)).set_group_id(m).set_failure_kind(ConstraintFailureKind::HARD).build()
        );
        if (evaluation.usage(m).sigma != SatisfactionPrescriptionKind::TRUE) {
            result.push_back(ConstraintBuilder<A>([evaluation,h,m](I const&, O const& o) {
                if (evaluation.usage(m).sigma == SatisfactionPrescriptionKind::FALSE_FOR_ALL) {
                    return -evaluate_from_function(h[m], o.evolve).upper().get_d();
                } else {
                    if (evaluate_from_function(h[m],o.evolve).lower().get_d() < 0) {
                        try {
                            auto approximator = NonlinearCandidateValidationInnerApproximator(ParallelLinearisationContractor(NativeSimplex(),2,1));
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

BoundingBoxType bounding_box(Orbit<LabelledEnclosure> const& orbit) {
    BoundingBoxType result = orbit.initial().euclidean_set().bounding_box();
    for (auto const& e : orbit.reach()) {
        result = hull(result,e.euclidean_set().bounding_box());
    }
    return result;
}

ConstrainedEvolutionResult constrained_evolution(VectorField const& dynamics, RealExpressionBoundedConstraintSet const& initial_set, Real const& evolution_time,
                                                 Configuration<VectorFieldEvolver> const& configuration, List<RealExpression> const& constraints) {
    auto satisfaction = ConstraintSatisfaction(constraints,dynamics.state_space());
    return constrained_evolution(dynamics,initial_set,evolution_time,configuration,constraints,satisfaction);
}

ConstrainedEvolutionResult constrained_evolution(VectorField const& dynamics, RealExpressionBoundedConstraintSet const& initial_set, Real const& evolution_time,
                                                 Configuration<VectorFieldEvolver> const& configuration, List<RealExpression> const& constraints, ConstraintSatisfaction const& constraint_satisfaction) {
    CONCLOG_SCOPE_CREATE

    auto satisfaction = constraint_satisfaction;
    auto h = satisfaction.indeterminate_constraints_function();

    Orbit<LabelledEnclosure> approximate_orbit({});
    Orbit<LabelledEnclosure> rigorous_orbit({});
    Orbit<LabelledEnclosure> controlled_orbit({});

    size_t num_indeterminates = constraints.size();

    for (size_t i=0; ; ++i) {
        CONCLOG_PRINTLN("Round #" << i)

        for (size_t j=0; ; ++j) {
            CONCLOG_PRINTLN_AT(1,"Singleton evolution try #" << j)
            auto analysis_point = configuration.search_space().initial_point();
            CONCLOG_PRINTLN_AT(1,"Using point " << analysis_point << " for singleton analyses.")

            auto singleton_configuration = make_singleton(configuration, analysis_point);
            singleton_configuration.set_enable_clobbering(true);

            VectorFieldEvolver approximate_evolver(dynamics, singleton_configuration);

            CONCLOG_PRINTLN_AT(1,"Computing approximate singleton evolution...")
            approximate_orbit = approximate_evolver.orbit(initial_set, evolution_time, Semantics::LOWER);

            CONCLOG_PRINTLN_AT(1,"Computing rigorous singleton evolution...")
            auto approximate_bounding_box = bounding_box(approximate_orbit);
            singleton_configuration.set_enable_clobbering(false).set_enable_premature_termination(true).set_maximum_enclosure_radius(approximate_bounding_box.radius().get_d());

            VectorFieldEvolver rigorous_evolver(dynamics, singleton_configuration);

            h = satisfaction.indeterminate_constraints_function();
            auto uncontrolled_task_constraints = build_uncontrolled_task_constraints(h);
            rigorous_evolver.set_constraints(uncontrolled_task_constraints);
            rigorous_orbit = rigorous_evolver.orbit(initial_set, evolution_time, Semantics::LOWER);

            satisfaction.merge_from_uncontrolled(rigorous_evolver.constraining_state());

            if (not rigorous_orbit.final().empty()) break;
            else { CONCLOG_PRINTLN_AT(1,"Early terminated due to set radius, retrying...") }
        }

        num_indeterminates = satisfaction.indeterminate_indexes().size();

        CONCLOG_PRINTLN_VAR(satisfaction)

        if (num_indeterminates > 0) {

            CONCLOG_PRINTLN("Processing singleton evolutions for constraints... ")

            h = satisfaction.indeterminate_constraints_function();
            auto analysis = evaluate_singleton_orbits(approximate_orbit, rigorous_orbit, h);

            auto controlled_task_constraints = build_controlled_task_constraints(h, analysis);

            auto indeterminate_idxs = satisfaction.indeterminate_indexes();
            for (size_t m = 0; m < analysis.number_of_constraints(); ++m) {
                CONCLOG_PRINTLN(constraints.at(indeterminate_idxs.at(m)) << " >= 0 : " << analysis.usage(m))
            }

            auto controlled_configuration = configuration;
            auto approximate_bounding_box = bounding_box(approximate_orbit);
            controlled_configuration.set_enable_premature_termination(true).set_maximum_enclosure_radius(
                    approximate_bounding_box.radius().get_d());
            VectorFieldEvolver controlled_evolver(dynamics, controlled_configuration);

            CONCLOG_PRINTLN("Computing controlled evolution... ")
            controlled_evolver.set_constraints(controlled_task_constraints);

            controlled_orbit = controlled_evolver.orbit(initial_set, evolution_time, Semantics::LOWER);

            satisfaction.merge_from_controlled(controlled_evolver.constraining_state(), analysis);

            if (satisfaction.indeterminate_indexes().size() == num_indeterminates) {
                CONCLOG_PRINTLN("No improvement in this round, aborting.")
                break;
            }

            num_indeterminates = satisfaction.indeterminate_indexes().size();

            if (num_indeterminates == 0) {
                CONCLOG_PRINTLN("All constraints satisfiability determined, terminating.")
                break;
            }
        } else {
            CONCLOG_PRINTLN("Completed during singleton analysis, terminating earlier (controlled orbit will be the previous one).")
            break;
        }
    }

    return {approximate_orbit,rigorous_orbit,controlled_orbit,satisfaction};
}

} // namespace Ariadne
