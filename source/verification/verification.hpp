/***************************************************************************
 *            verification.hpp
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

#ifndef ARIADNE_VERIFICATION
#define ARIADNE_VERIFICATION

#include "numeric/float_bounds.hpp"
#include "function/function.hpp"
#include "dynamics/vector_field_evolver.hpp"
#include "dynamics/inner_approximation.hpp"
#include "helper/stopwatch.hpp"
#include "pexplore/task_runner.tpl.hpp"

namespace Ariadne {

using std::ostream;
using Helper::Milliseconds;
using Helper::Stopwatch;
using pExplore::ConstrainingState;
using pExplore::ConstraintBuilder;
using pExplore::ConstraintObjectiveImpact;
using pExplore::ConstraintSuccessAction;
using pExplore::ConstraintFailureKind;
using pExplore::TaskInput;
using pExplore::TaskOutput;
using ConstraintIndexType = size_t;

double nthroot(double value, size_t n);

//! \brief The prescription to the satisfaction of the constraint
//! \details TRUE : always true for all points (required over-approximation)
//!          FALSE_FOR_ALL : sometimes false for all points (can use over-approximation)
//!          FALSE_FOR_SOME : sometimes false only for some points (must use under-approximation)
enum class SatisfactionPrescriptionKind { TRUE, FALSE_FOR_ALL, FALSE_FOR_SOME };
std::ostream& operator<<(std::ostream& os, const SatisfactionPrescriptionKind prescription) {
    switch(prescription) {
        case SatisfactionPrescriptionKind::TRUE: os << "TRUE"; return os;
        case SatisfactionPrescriptionKind::FALSE_FOR_ALL: os << "FALSE_FOR_ALL"; return os;
        case SatisfactionPrescriptionKind::FALSE_FOR_SOME: os << "FALSE_FOR_SOME"; return os;
        default: HELPER_FAIL_MSG("Unhandled SatisfactionPrescriptionKind value for printing")
    }
}

class EvaluationSequence;
class ConstrainedEvolutionResult;
class ConstraintSatisfaction;

FloatDPBounds evaluate_from_function(EffectiveScalarMultivariateFunction const& function, LabelledEnclosure const& enclosure);
Vector<FloatDPBounds> evaluate_from_function(EffectiveVectorMultivariateFunction const& function, LabelledEnclosure const& enclosure);
Vector<FloatDPBounds> widen(Vector<FloatDPBounds> const& bx, double chi);
Vector<FloatDPBounds> shrink(Vector<FloatDPBounds> const& bx, double chi);

List<pExplore::Constraint<VectorFieldEvolver>> build_uncontrolled_task_constraints(EffectiveVectorMultivariateFunction const& h);
//! \brief Build the controlled constraints; we \a exclude_truth if the rigorous orbit from preanalysis was cut short due to the set being too large
List<pExplore::Constraint<VectorFieldEvolver>> build_controlled_task_constraints(EffectiveVectorMultivariateFunction const& h, EvaluationSequence const& evaluation, bool exclude_truth);
ConstrainedEvolutionResult constrained_evolution(VectorField const& dynamics, RealExpressionBoundedConstraintSet const& initial_set, Real const& evolution_time,
                                                 Configuration<VectorFieldEvolver> const& configuration, List<RealExpression> const& constraints);

double get_chi(Vector<FloatDPBounds>const& bnds, EffectiveScalarMultivariateFunction const& constraint, SatisfactionPrescriptionKind prescription);
double get_rho(double chi, double beta, SatisfactionPrescriptionKind prescription);
double get_beta(Vector<FloatDPBounds> const& bnds, size_t n);
Vector<FloatDPBounds> widen(Vector<FloatDPBounds> const& bx, double chi);
Vector<FloatDPBounds> shrink(Vector<FloatDPBounds> const& bx, double chi);

struct TimedBoxEvaluation {
    TimedBoxEvaluation(double time_, Vector<FloatDPBounds> const& approximate_box_, Vector<FloatDPBounds> const& rigorous_box_, Vector<FloatDPBounds> const& approximate_evaluation_) :
        time(time_), approximate_box(approximate_box_), rigorous_box(rigorous_box_), approximate_evaluation(approximate_evaluation_) { }
    double time;
    Vector<FloatDPBounds> approximate_box;
    Vector<FloatDPBounds> rigorous_box;
    Vector<FloatDPBounds> approximate_evaluation;
    friend ostream& operator<<(ostream& os, TimedBoxEvaluation& tbe) { return os << tbe.time << ": a_box=" << tbe.approximate_box << ", r_box=" << tbe.rigorous_box << ", eval=" << tbe.approximate_evaluation; }
};

struct ControlSpecification {
    ControlSpecification(SatisfactionPrescriptionKind sigma_, double t_star_, double alpha_) : sigma(sigma_), t_star(t_star_), alpha(alpha_) { }
    SatisfactionPrescriptionKind sigma;
    double t_star;
    double alpha;
    friend ostream& operator<<(ostream& os, ControlSpecification const& hcc) { return os << "{sigma=" << hcc.sigma << ",T*=" << hcc.t_star << ",alpha=" << hcc.alpha << "}"; }
};

//! \brief A timed measurement in terms of beta (normalised volume) both for approximate and rigorous reachability
struct TimedMeasurement {
    TimedMeasurement(double time_, double approximate_beta_, double rigorous_beta_) : time(time_), approximate_beta(approximate_beta_), rigorous_beta(rigorous_beta_) { }
    double time;
    double approximate_beta;
    double rigorous_beta;
    friend ostream& operator<<(ostream& os, TimedMeasurement const& tm) { return os << tm.time << ":" << tm.approximate_beta << "(a)," << tm.rigorous_beta << "(r)"; }
};

struct ConstraintPrescription {
    ConstraintPrescription(RealExpression const& expression_, SatisfactionPrescriptionKind const& prescription_) : expression(expression_), prescription(prescription_) { }
    RealExpression expression;
    SatisfactionPrescriptionKind prescription;
    friend ostream& operator<<(ostream& os, ConstraintPrescription const& cp) { return os << cp.expression << " >= 0 : " << cp.prescription; }
};

class ConstraintSatisfaction;

class ConstraintSatisfactionSnapshot {
    friend ConstraintSatisfaction;
  public:
    ConstraintSatisfactionSnapshot(double time, Vector<SatisfactionPrescriptionKind> const& prescriptions, Vector<LogicalValue> const& outcomes);

    double time() const;
    Vector<SatisfactionPrescriptionKind> const& prescriptions() const;
    Vector<LogicalValue> const& outcomes() const;

    Map<SatisfactionPrescriptionKind,double> prescription_ratios() const;
    Map<SatisfactionPrescriptionKind,double> success_ratios() const;
    double global_success_ratio() const;
    List<size_t> indeterminate_indexes() const;
    bool completed() const;
    friend ostream& operator<<(ostream& os, ConstraintSatisfactionSnapshot const& s);

  protected:
    void set_prescription(size_t m, SatisfactionPrescriptionKind const& prescription);
    void set_outcome(size_t m, bool outcome);

  private:
    double _time;
    Vector<SatisfactionPrescriptionKind> _prescriptions;
    Vector<LogicalValue> _outcomes;
};

class ConstraintSatisfaction {
  public:
    ConstraintSatisfaction(List<RealExpression> const& cs, RealSpace const& spc);

    //! \brief The dimension of the satisfaction, i.e., the number of constraints
    size_t dimension() const;

    //! \brief Construct a vector function from the constraints that are still indeterminate
    EffectiveVectorMultivariateFunction indeterminate_constraints_function() const;

    void merge_from_uncontrolled(ConstrainingState<VectorFieldEvolver> const& state, bool exclude_truth);
    void merge_from_controlled(ConstrainingState<VectorFieldEvolver> const& state, EvaluationSequence const& preanalysis, bool exclude_truth);
    List<size_t> indeterminate_indexes() const;
    bool completed() const;

    Map<SatisfactionPrescriptionKind,double> prescription_ratios() const;
    Map<SatisfactionPrescriptionKind,double> success_ratios() const;
    double global_success_ratio() const;

    //! \brief The integral across snapshots of the global success ratio represents the cost
    double cost() const;

    RealExpression const& expression(size_t m) const;
    SatisfactionPrescriptionKind const& prescription(size_t m) const;
    LogicalValue const& outcome(size_t m) const;

    List<ConstraintSatisfactionSnapshot> const& snapshots() const;

    friend ostream& operator<<(ostream& os, ConstraintSatisfaction const& cs);

  private:
    ConstraintSatisfactionSnapshot const& _last() const;
    ConstraintSatisfactionSnapshot& _last();
    void _add_snapshot();
    void _set_last_outcome(size_t m, bool outcome);
    void _set_last_prescription(size_t m, SatisfactionPrescriptionKind prescription);

  private:
    Vector<RealExpression> const _cs;
    RealSpace const _space;
    List<ConstraintSatisfactionSnapshot> _snapshots;
    Stopwatch<Milliseconds> _sw;
};

class ConstrainedEvolutionResult {
  public:
    ConstrainedEvolutionResult(Orbit<LabelledEnclosure> const& approximate, Orbit<LabelledEnclosure> const& rigorous, Orbit<LabelledEnclosure> const& constrained, ConstraintSatisfaction const& satisfaction) :
        _approximate(approximate), _rigorous(rigorous), _constrained(constrained), _satisfaction(satisfaction) { }


    Orbit<LabelledEnclosure> const& approximate() const { return _approximate; }
    Orbit<LabelledEnclosure> const& rigorous() const { return _rigorous; }
    Orbit<LabelledEnclosure> const& constrained() const { return _constrained; }
    ConstraintSatisfaction const& satisfaction() const { return _satisfaction; }

  private:
    Orbit<LabelledEnclosure> const _approximate;
    Orbit<LabelledEnclosure> const _rigorous;
    Orbit<LabelledEnclosure> const _constrained;
    ConstraintSatisfaction const _satisfaction;
};

class EvaluationSequenceBuilder;

class EvaluationSequence {
    friend class EvaluationSequenceBuilder;
  protected:
    EvaluationSequence(List<TimedMeasurement> const& tv, Vector<ControlSpecification> const& usages);
  public:

    //! \brief The number of constraints inserted during creation
    size_t number_of_constraints() const;

    //! \brief The maximum time registered by the sequence
    double maximum_time() const;

    //! \brief Timed measurement accessor by index
    TimedMeasurement const& at(size_t idx) const;

    //! \brief Get the time measurement with time closest to \a time
    TimedMeasurement const& near(double time) const;

    //! \brief The usage given the constraint index \a idx
    ControlSpecification const& usage(ConstraintIndexType idx) const;

    //! \brief The number of elements in the sequence
    size_t size() const;

    friend ostream& operator<<(ostream& os, EvaluationSequence const& es);

  private:
    List<TimedMeasurement> _sequence;
    Vector<ControlSpecification> _prescriptions;
};

class EvaluationSequenceBuilder {
  public:

    EvaluationSequenceBuilder(size_t N, EffectiveVectorMultivariateFunction const& h);

    //! \brief Add an element of the sequence from \a approximate and \a rigorous enclosures
    //! \details Assumes that the time has already been chosen the closest between the two and will use the approximate one
    void add_from(LabelledEnclosure const& approximate, LabelledEnclosure const& rigorous);

    void add(TimedBoxEvaluation const& tbe);

    EvaluationSequence build() const;

  private:
    size_t _list_index(double time) const;

  private:
    size_t const _N;
    size_t const _M;
    EffectiveVectorMultivariateFunction const _h;
    List<TimedBoxEvaluation> _timed_box_evaluations;

    Vector<SatisfactionPrescriptionKind> _prescriptions;
    Vector<Map<SatisfactionPrescriptionKind,double>> _max_robustness_false;
    Vector<Map<SatisfactionPrescriptionKind,double>> _t_false;
};

} // namespace Ariadne

#endif // ARIADNE_VERIFICATION