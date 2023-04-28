/***************************************************************************
 *            constrained_vanderpol.cpp
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

#include "ariadne_main.hpp"
#include "dynamics/inner_approximation.hpp"
#include "helper/stopwatch.hpp"
#include "pexplore/task_runner.tpl.hpp"

using std::ostream;

using pExplore::ConstrainingSpecification;
using pExplore::ConstraintBuilder;
using pExplore::ConstraintObjectiveImpact;
using pExplore::ConstraintSuccessAction;
using pExplore::ConstraintFailureKind;
using pExplore::TaskInput;
using pExplore::TaskOutput;
using ConstraintIndexType = size_t;

//! \brief The prescription to the satisfaction of the constraint
//! \details TRUE : always true for all points (required over-approximation)
//!          FALSE_FOR_ALL : sometimes false for all points (can use over-approximation)
//!          FALSE_FOR_SOME : sometimes false only for some points (must use under-approximation)
enum class SatisfactionPrescription { TRUE, FALSE_FOR_ALL, FALSE_FOR_SOME };
std::ostream& operator<<(std::ostream& os, const SatisfactionPrescription prescription) {
    switch(prescription) {
        case SatisfactionPrescription::TRUE: os << "TRUE"; return os;
        case SatisfactionPrescription::FALSE_FOR_ALL: os << "FALSE_FOR_ALL"; return os;
        case SatisfactionPrescription::FALSE_FOR_SOME: os << "FALSE_FOR_SOME"; return os;
        default: HELPER_FAIL_MSG("Unhandled SatisfactionPrescription value for printing")
    }
}

FloatDPBounds evaluate_from_function(EffectiveScalarMultivariateFunction const& function, LabelledEnclosure const& enclosure);
double get_chi(Vector<FloatDPBounds>const& bnds, EffectiveScalarMultivariateFunction const& constraint, SatisfactionPrescription prescription);
double get_rho(double chi, size_t N, double volume, SatisfactionPrescription prescription);
double volume(Vector<FloatDPBounds> const& bnds);

struct TimedBoxEvaluation {
    TimedBoxEvaluation(double time_, Vector<FloatDPBounds> const& box_, Vector<FloatDPBounds> const& evaluation_) : time(time_), box(box_), evaluation(evaluation_) { }
    double time;
    Vector<FloatDPBounds> box;
    Vector<FloatDPBounds> evaluation;
    friend ostream& operator<<(ostream& os, TimedBoxEvaluation& tbe) { return os << tbe.time << ": box=" << tbe.box << ", eval=" << tbe.evaluation; }
};

struct HardConstraintPrescription {
  HardConstraintPrescription(SatisfactionPrescription sigma_, double t_star_, double alpha_) : sigma(sigma_), t_star(t_star_), alpha(alpha_) { }
  SatisfactionPrescription sigma;
  double t_star;
  double alpha;
  friend ostream& operator<<(ostream& os, HardConstraintPrescription const& hcc) { return os << "{sigma=" << hcc.sigma << ",T*=" << hcc.t_star << ",alpha=" << hcc.alpha << "}"; }
};

struct TimedVolume {
  TimedVolume(double time_, double volume_) : time(time_), volume(volume_) { }
  double time;
  double volume;
  friend ostream& operator<<(ostream& os, TimedVolume const& tv) { return os << tv.time << ":" << tv.volume; }
};

class EvaluationSequenceBuilder;

class EvaluationSequence {
    friend class EvaluationSequenceBuilder;
  protected:
    EvaluationSequence(List<TimedVolume> const& tv, Vector<HardConstraintPrescription> const& usages) : _sequence(tv), _prescriptions(usages) { }
  public:

    size_t number_of_constraints() const { return _prescriptions.size(); }

    //! \brief timed volume accessor by index
    TimedVolume const& at(size_t idx) const { return _sequence.at(idx); }

    //! \brief Get the volume with time closest to \a time
    double const& near(double time) const {
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

    //! \brief The usage given the constraint index \a idx
    HardConstraintPrescription const& usage(ConstraintIndexType idx) const { return _prescriptions.at(idx); }

    //! \brief The number of elements in the sequence
    size_t size() const { return _sequence.size(); }

    friend ostream& operator<<(ostream& os, EvaluationSequence const& es) {
        os << "timed_volumes:{";
        for (size_t i=0; i<es.size()-1; ++i) os << es.at(i) << ",";
        return os << es.at(es.size()-1) << "}, usages:" << es._prescriptions;
    }

  private:
    List<TimedVolume> _sequence;
    Vector<HardConstraintPrescription> _prescriptions;
};

class EvaluationSequenceBuilder {
  public:

    EvaluationSequenceBuilder(size_t N, List<EffectiveScalarMultivariateFunction> const& hs) : _N(N), _M(hs.size()), _hs(hs), _prescriptions(hs.size(),SatisfactionPrescription::TRUE),
                                          _max_robustness_false(hs.size(),{{SatisfactionPrescription::FALSE_FOR_ALL,0.0},{SatisfactionPrescription::FALSE_FOR_SOME,0.0}}),
                                          _t_false(hs.size(),{{SatisfactionPrescription::FALSE_FOR_ALL,0.0},{SatisfactionPrescription::FALSE_FOR_SOME,0.0}}) { }

    void add_from(LabelledEnclosure const& e) {
        Vector<FloatDPBounds> eval(_M,DoublePrecision());
        for (size_t m=0; m<_M; ++m) eval.at(m) = evaluate_from_function(_hs.at(m),e);

        auto bb = e.euclidean_set().bounding_box();
        add({e.time_function().range().midpoint().get_d(),reinterpret_cast<Vector<FloatDPBounds>const&>(bb),eval});
    }

    void add(TimedBoxEvaluation const& tbe) {
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
                auto chi = get_chi(tbe.box,_hs[m],prescription_to_check);
                auto robustness = get_rho(chi,_N,volume(tbe.box),prescription_to_check);
                if (_max_robustness_false[m].get(prescription_to_check) < robustness) {
                    _max_robustness_false[m].at(prescription_to_check) = robustness;
                    _t_false[m].at(prescription_to_check) = tbe.time;
                }
            }
        }
    }

    EvaluationSequence build() const {
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
                        double chi = get_chi(tbe.box,_hs[m],_prescriptions[m]);
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

  private:

    size_t _list_index(double time) const {
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

  private:
    size_t const _N;
    size_t const _M;
    List<EffectiveScalarMultivariateFunction> const _hs;
    List<TimedBoxEvaluation> _timed_box_evaluations;

    Vector<SatisfactionPrescription> _prescriptions;
    Vector<Map<SatisfactionPrescription,double>> _max_robustness_false;
    Vector<Map<SatisfactionPrescription,double>> _t_false;
};

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

List<EffectiveScalarMultivariateFunction> convert_to_functions(List<RealExpression> const& constraints, RealSpace const& spc) {
    List<EffectiveScalarMultivariateFunction> result;
    for (auto const& c : constraints)
        result.push_back(make_function(spc,c));
    return result;
}

FloatDPBounds evaluate_from_function(EffectiveScalarMultivariateFunction const& function, LabelledEnclosure const& enclosure) {
    auto bb = enclosure.euclidean_set().bounding_box();
    return function(reinterpret_cast<Vector<FloatDPBounds>const&>(bb));
}

EvaluationSequence evaluate_approximate_orbit(Orbit<LabelledEnclosure> const& orbit, List<EffectiveScalarMultivariateFunction> const& hs) {
    CONCLOG_SCOPE_CREATE

    EvaluationSequenceBuilder sb(orbit.initial().dimension(),hs);
    sb.add_from(orbit.initial());
    for (auto const& e : orbit.intermediate()) { sb.add_from(e); }
    return sb.build();
}

ConstrainingSpecification<VectorFieldEvolver> build_constraining_specification(EvaluationSequence const& evaluation, List<EffectiveScalarMultivariateFunction> const& hs) {
    using A = VectorFieldEvolver;
    using I = TaskInput<A>;
    using O = TaskOutput<A>;

    List<pExplore::Constraint<A>> constraints;

    for (ConstraintIndexType m=0; m<evaluation.number_of_constraints(); ++m) {
        constraints.push_back(ConstraintBuilder<A>([evaluation,m](I const&, O const& o) {
            auto v0 = evaluation.at(0).volume;
            auto v_approx = evaluation.near(o.time.get_d());
            auto alpha = evaluation.usage(m).alpha;
            return (alpha+1)*v_approx - alpha*v0 - o.evolve.euclidean_set().bounding_box().volume().get_d();
        }).set_name("objective_soft"+to_string(m)).set_group_id(m).set_objective_impact(ConstraintObjectiveImpact::UNSIGNED).set_failure_kind(ConstraintFailureKind::SOFT).build()
        );
        constraints.push_back(ConstraintBuilder<A>([evaluation,hs,m](I const&, O const& o) {
            if (evaluation.usage(m).sigma == SatisfactionPrescription::TRUE)
                return evaluate_from_function(hs.at(m),o.reach).lower().get_d();
            else return evaluation.usage(m).t_star - o.time.get_d();
        }).set_name("hard"+to_string(m)).set_group_id(m).set_failure_kind(ConstraintFailureKind::HARD).build()
        );
        if (evaluation.usage(m).sigma != SatisfactionPrescription::TRUE) {
            constraints.push_back(ConstraintBuilder<A>([evaluation,hs,m](I const&, O const& o) {
                if (evaluation.usage(m).sigma == SatisfactionPrescription::FALSE_FOR_ALL) {
                    return -evaluate_from_function(hs.at(m), o.evolve).upper().get_d();
                } else {
                    if (evaluate_from_function(hs.at(m),o.evolve).lower().get_d() < 0) {
                        try {
                            auto approximator = NonlinearCandidateValidationInnerApproximator(ParallelLinearisationContractor(GLPKSimplex(),2,1));
                            auto inner_evolve = approximator.compute_from(o.evolve);
                            return -evaluate_from_function(hs.at(m),inner_evolve).lower().get_d();
                        } catch (std::exception&) { }
                    }
                    return -1.0;
                }
            }).set_name("falsify_success"+to_string(m)).set_group_id(m).set_success_action(ConstraintSuccessAction::DEACTIVATE).build()
            );
        }
    }

    return {constraints};
}

void ariadne_main()
{
    RealConstant mu("mu",1);
    RealVariable x("x"), y("y");

    VectorField dynamics({dot(x)=y, dot(y)= mu*y*(1-sqr(x))-x});

    RealConstant ymax("ymax",2.8_x);
    RealConstant ymin("ymin",-3.0_x);
    RealConstant xmin("xmin",-1.0_x);
    RealConstant xmax("xmax",2.5_x);
    RealConstant rsqr("r^2",2.0_x);
    List<RealExpression> constraints = {y - ymin, x - xmin, ymax - y, xmax - x, sqr(x) + sqr(y) - rsqr};
    //List<RealExpression> constraints = {x-xmin};

    auto approximate_configuration = Configuration<VectorFieldEvolver>().
        set_maximum_step_size(0.1);

    AffineIntegrator approximate_integrator(2,2);
    VectorFieldEvolver approximate_evolver(dynamics,approximate_configuration,approximate_integrator);
    CONCLOG_PRINTLN_VAR(approximate_evolver.configuration())

    auto rigorous_configuration = Configuration<VectorFieldEvolver>().
            set_maximum_enclosure_radius(1.0).
            set_maximum_step_size(0.02).
            set_maximum_spacial_error(1e-6);

    StepMaximumError max_err=1e-6;
    TaylorPicardIntegrator rigorous_integrator(max_err);
    VectorFieldEvolver rigorous_evolver(dynamics,rigorous_configuration,rigorous_integrator);
    CONCLOG_PRINTLN_VAR(rigorous_evolver.configuration())

    Real x0 = 1.40_dec;
    Real y0 = 2.40_dec;
    Real eps_x0 = 0.15_dec;
    Real eps_y0 = 0.05_dec;

    RealExpressionBoundedConstraintSet initial_set({x0-eps_x0<=x<=x0+eps_x0,y0-eps_y0<=y<=y0+eps_y0});

    CONCLOG_PRINTLN("Initial set: " << initial_set)
    Real evolution_time = 7;

    Helper::Stopwatch<std::chrono::microseconds> sw;
    CONCLOG_PRINTLN("Computing approximate evolution...")
    auto approximate_orbit = approximate_evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    sw.click();
    CONCLOG_PRINTLN_AT(0,"Done in " << sw.elapsed_seconds() << " seconds.")

    sw.restart();
    CONCLOG_PRINTLN("Processing approximate evolution for constraints... ")

    auto hs = convert_to_functions(constraints, dynamics.state_space());
    CONCLOG_PRINTLN_VAR(hs)

    auto analysis = evaluate_approximate_orbit(approximate_orbit, hs);

    auto constraining = build_constraining_specification(analysis, hs);

    sw.click();
    CONCLOG_PRINTLN_AT(0,"Done in " << sw.elapsed_seconds() << " seconds.")

    CONCLOG_PRINTLN_VAR_AT(1,analysis)

    CONCLOG_PRINTLN("Plotting...")
    LabelledFigure fig=LabelledFigure({-3<=x<=3,-3<=y<=3});
    fig.draw(approximate_orbit);
    CONCLOG_RUN_MUTED(fig.write("vanderpol_approximate"))

    sw.restart();
    CONCLOG_PRINTLN("Computing rigorous evolution... ")
    rigorous_evolver.set_constraining(constraining);
    auto rigorous_orbit = rigorous_evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    sw.click();
    CONCLOG_PRINTLN_AT(0,"Done in " << sw.elapsed_seconds() << " seconds.")

    CONCLOG_PRINTLN("Plotting...")
    fig.clear();
    fig.draw(rigorous_orbit);
    CONCLOG_RUN_MUTED(fig.write("vanderpol_rigorous"))
}
