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

#include "helper/stopwatch.hpp"
#include "ariadne_main.hpp"

using std::ostream;

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

struct TimedBounds {
  TimedBounds(double time_, Vector<FloatDPBounds> const& bounds_) : time(time_), bounds(bounds_) { }
  double time;
  Vector<FloatDPBounds> bounds;
  friend ostream& operator<<(ostream& os, TimedBounds& tb) { return os << tb.time << ":" << tb.bounds; }
};

struct HardConstraintConstants {
  HardConstraintConstants(SatisfactionPrescription mu_, double t_star_, double alpha_) : mu(mu_), t_star(t_star_), alpha(alpha_) { }
  SatisfactionPrescription mu;
  double t_star;
  double alpha;
  friend ostream& operator<<(ostream& os, HardConstraintConstants const& hcc) { return os << "{mu=" << hcc.mu << ",T*=" << hcc.t_star << ",alpha=" << hcc.alpha << "}"; }
};

struct TimedWidths {
  TimedWidths(double time_, Map<ConstraintIndexType,double> const& widths_) : time(time_), widths(widths_) { }
  double time;
  Map<ConstraintIndexType,double> widths;
  friend ostream& operator<<(ostream& os, TimedWidths const& te) { return os << te.time << ":" << te.widths; }
};

class EvaluationSequenceBuilder;

class EvaluationSequence {
    friend class EvaluationSequenceBuilder;
  protected:
    EvaluationSequence(List<TimedWidths> const& te, Vector<HardConstraintConstants> const& usages) : _sequence(te), _usages(usages) { }
  public:
    TimedWidths const& at(size_t idx) const { return _sequence.at(idx); }
    //! \brief Get the widths with time closest to \a time
    Map<ConstraintIndexType,double> const& near(double time) const {
        size_t lower = 0;
        size_t upper = size()-1;

        if (_sequence.at(lower).time >= time) return _sequence.at(lower).widths;
        if (_sequence.at(upper).time <= time) return _sequence.at(upper).widths;

        while (upper - lower > 1) {
            auto mid = (lower+upper)/2;
            if (_sequence.at(mid).time < time) lower = mid;
            else upper = mid;
        }
        if (time - _sequence.at(lower).time > _sequence.at(upper).time - time) return _sequence.at(upper).widths;
        else return _sequence.at(lower).widths;
    }

    size_t size() const { return _sequence.size(); }

    friend ostream& operator<<(ostream& os, EvaluationSequence const& es) {
        os << "{";
        for (size_t i=0; i<es.size()-1; ++i) os << es.at(i) << ",";
        return os << es.at(es.size()-1) << "}";
    }

  private:
    List<TimedWidths> _sequence;
    Vector<HardConstraintConstants> _usages;
};

class EvaluationSequenceBuilder {
  public:

    EvaluationSequenceBuilder(size_t M) : _M(M), _prescriptions(M,SatisfactionPrescription::TRUE),
                                          _max_robustness_false(M,{{SatisfactionPrescription::FALSE_FOR_ALL,0.0},{SatisfactionPrescription::FALSE_FOR_SOME,0.0}}),
                                          _t_false(M,{{SatisfactionPrescription::FALSE_FOR_ALL,0.0},{SatisfactionPrescription::FALSE_FOR_SOME,0.0}}) { }

    void add(TimedBounds const& tb) {
        HELPER_PRECONDITION(_timed_bounds_list.empty() or _timed_bounds_list.at(_timed_bounds_list.size()-1).time < tb.time)
        _timed_bounds_list.push_back(tb);
        for (ConstraintIndexType m=0; m<_M; ++m) {
            auto const& bounds = tb.bounds[m];
            auto upper = bounds.upper().get_d();
            auto lower = bounds.lower().get_d();
            if (_prescriptions[m] != SatisfactionPrescription::FALSE_FOR_ALL) {
                if (upper < 0) _prescriptions[m] = SatisfactionPrescription::FALSE_FOR_ALL;
                else if (lower < 0) _prescriptions[m] = SatisfactionPrescription::FALSE_FOR_SOME;
            }
            if (upper < 0) {
                if (_max_robustness_false[m].get(SatisfactionPrescription::FALSE_FOR_ALL) < -upper) {
                    _max_robustness_false[m].at(SatisfactionPrescription::FALSE_FOR_ALL) = -upper;
                    _t_false[m].at(SatisfactionPrescription::FALSE_FOR_ALL) = tb.time;
                }
            } else if (lower < 0) {
                if (_max_robustness_false[m].get(SatisfactionPrescription::FALSE_FOR_SOME) < -lower) {
                    _max_robustness_false[m].at(SatisfactionPrescription::FALSE_FOR_SOME) = -lower;
                    _t_false[m].at(SatisfactionPrescription::FALSE_FOR_SOME) = tb.time;
                }
            }
        }
    }

    EvaluationSequence build() const {
        HELPER_PRECONDITION(not _timed_bounds_list.empty())

        Vector<double> T_star(_M, 0.0);
        for (ConstraintIndexType m=0; m<_M; ++m) {
            switch (_prescriptions[m]) {
                case SatisfactionPrescription::TRUE :           T_star[m] = _timed_bounds_list.at(_timed_bounds_list.size()-1).time; break;
                case SatisfactionPrescription::FALSE_FOR_ALL :  T_star[m] = _t_false[m].get(SatisfactionPrescription::FALSE_FOR_ALL); break;
                case SatisfactionPrescription::FALSE_FOR_SOME : T_star[m] = _t_false[m].get(SatisfactionPrescription::FALSE_FOR_SOME); break;
                default: HELPER_FAIL_MSG("Unhandled SatisfactionPrescription value")
            }
        }

        List<TimedWidths> timed_widths;

        Vector<double> width0(_M,0.0);
        for (ConstraintIndexType m=0; m<_M; ++m)
            width0[m] = _timed_bounds_list.at(0).bounds[m].upper().get_d()-_timed_bounds_list.at(0).bounds[m].lower().get_d();

        Vector<double> alpha(_M, std::numeric_limits<double>::max());
        for (size_t i=1; i < _timed_bounds_list.size(); ++i) {
            auto const& tb = _timed_bounds_list.at(i);
            Map<ConstraintIndexType,double> widths;
            for (ConstraintIndexType m=0; m<_M; ++m) {
                if (tb.time <= T_star[m]) {
                    double width = tb.bounds[m].upper().get_d()-tb.bounds[m].lower().get_d();
                    widths.insert(m,width);
                    if (_prescriptions[m] == SatisfactionPrescription::TRUE or tb.time == T_star[m]) {
                        double rho = 0;
                        switch (_prescriptions[m]) {
                            case SatisfactionPrescription::TRUE : rho = tb.bounds[m].lower().get_d(); break;
                            case SatisfactionPrescription::FALSE_FOR_ALL : rho = -tb.bounds[m].upper().get_d(); break;
                            case SatisfactionPrescription::FALSE_FOR_SOME : rho = -tb.bounds[m].lower().get_d(); break;
                            default: HELPER_FAIL_MSG("Unhandled SatisfactionPrescription value")
                        }
                        //CONCLOG_PRINTLN("i="<< i <<",m="<< m << ",rho="<< rho << ",w-w0=" << width-width0[m])
                        alpha[m] = std::min(alpha[m],rho/(width-width0[m]));
                    }
                }
            }
            if (widths.empty()) break;
            timed_widths.push_back({tb.time,widths});
        }

        Vector<HardConstraintConstants> constants(_M, {SatisfactionPrescription::TRUE, 0.0, 0.0});
        for (ConstraintIndexType m=0; m<_M; ++m)
            constants[m] = {_prescriptions[m],T_star[m],alpha[m]};
        CONCLOG_PRINTLN_VAR(constants)

        return {timed_widths, constants};
    }

  private:

    size_t _list_index(double time) const {
        auto initial_time = _timed_bounds_list.at(0).time;
        auto final_time = _timed_bounds_list.at(_timed_bounds_list.size()-1).time;
        HELPER_PRECONDITION(initial_time <= time)
        HELPER_PRECONDITION(final_time >= time)

        size_t lower = 0;
        size_t upper = _timed_bounds_list.size()-1;

        if (time == initial_time) return 0;
        if (time == final_time) return _timed_bounds_list.size()-1;

        while (upper - lower > 0) {
            auto mid = (lower+upper)/2;
            auto midtime = _timed_bounds_list.at(mid).time;
            if (midtime == time) return mid;
            else if (midtime < time) lower = mid;
            else upper = mid;
        }
        HELPER_ASSERT_MSG(_timed_bounds_list.at(lower).time == time, "The time " << time << " could not be found in the list")
        return lower;
    }

    double _find_latest_false_time(ConstraintIndexType m, bool false_for_all) const {
        double starting_time = _t_false[m].get(false_for_all ? SatisfactionPrescription::FALSE_FOR_ALL : SatisfactionPrescription::FALSE_FOR_SOME);
        size_t idx = _list_index(starting_time)+1;
        size_t max_idx = _timed_bounds_list.size()-1;
        for ( ; idx<=max_idx; ++idx) {
            double bound = (false_for_all ? _timed_bounds_list.at(idx).bounds[m].upper().get_d() : _timed_bounds_list.at(idx).bounds[m].lower().get_d());
            if (bound >= 0) {
                break;
            }
        }
        return _timed_bounds_list.at(--idx).time;
    }

  private:
    size_t const _M;
    List<TimedBounds> _timed_bounds_list;

    Vector<SatisfactionPrescription> _prescriptions;
    Vector<Map<SatisfactionPrescription,double>> _max_robustness_false;
    Vector<Map<SatisfactionPrescription,double>> _t_false;
};

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

void evaluate_approximate_orbit(Orbit<LabelledEnclosure> const& orbit, List<EffectiveScalarMultivariateFunction> const& hs) {
    CONCLOG_SCOPE_CREATE

    auto M = hs.size();

    EvaluationSequenceBuilder sb(M);
    Vector<FloatDPBounds> initial_eval(M,DoublePrecision());
    for (size_t m=0; m<M; ++m)
        initial_eval.at(m) = evaluate_from_function(hs.at(m),orbit.initial());
    sb.add({orbit.initial().time_function().range().midpoint().get_d(),initial_eval});
    for (auto const& e : orbit.intermediate()) {
        Vector<FloatDPBounds> eval(M,DoublePrecision());
        for (size_t m=0; m<M; ++m)
            eval.at(m) = evaluate_from_function(hs.at(m),e);
        sb.add({e.time_function().range().midpoint().get_d(),eval});
    }
    auto es = sb.build();

    CONCLOG_PRINTLN_VAR_AT(1,es)
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
    List<RealExpression> constraints = {y - ymin, x - xmin, ymax - y, xmax - x};

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
    CONCLOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.")

    sw.restart();
    CONCLOG_PRINTLN("Processing approximate evolution for constraints... ")

    auto hs = convert_to_functions(constraints, dynamics.state_space());
    CONCLOG_PRINTLN_VAR(hs)

    evaluate_approximate_orbit(approximate_orbit, hs);
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.")

    CONCLOG_PRINTLN("Plotting...")
    LabelledFigure fig=LabelledFigure({-2.5<=x<=2.5,-3<=y<=3});
    fig.draw(approximate_orbit);
    CONCLOG_RUN_MUTED(fig.write("vanderpol_approximate"))

    sw.restart();
    CONCLOG_PRINTLN("Computing rigorous evolution... ")
    auto rigorous_orbit = rigorous_evolver.orbit(initial_set,evolution_time,Semantics::UPPER);
    sw.click();
    CONCLOG_PRINTLN_AT(1,"Done in " << sw.elapsed_seconds() << " seconds.")

    CONCLOG_PRINTLN("Plotting...")
    fig.clear();
    fig.draw(rigorous_orbit);
    CONCLOG_RUN_MUTED(fig.write("vanderpol_rigorous"))
}
