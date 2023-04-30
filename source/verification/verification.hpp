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

using pExplore::ConstrainingState;
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

class EvaluationSequence;

FloatDPBounds evaluate_from_function(EffectiveScalarMultivariateFunction const& function, LabelledEnclosure const& enclosure);
Vector<FloatDPBounds> widen(Vector<FloatDPBounds> const& bx, double chi);
Vector<FloatDPBounds> shrink(Vector<FloatDPBounds> const& bx, double chi);

List<pExplore::Constraint<VectorFieldEvolver>> build_task_constraints(EvaluationSequence const& evaluation, Vector<EffectiveScalarMultivariateFunction> const& hs);
Vector<Kleenean> synthesise_outcomes(EvaluationSequence const& preanalysis, ConstrainingState<VectorFieldEvolver> const& constraining);
Tuple<Orbit<LabelledEnclosure>,Orbit<LabelledEnclosure>,Vector<Kleenean>> constrained_evolution(VectorField const& dynamics, RealExpressionBoundedConstraintSet const& initial_set, Real const& evolution_time,
                                                                                                    List<RealExpression> const& constraints, Configuration<VectorFieldEvolver> const& configuration);

double get_chi(Vector<FloatDPBounds>const& bnds, EffectiveScalarMultivariateFunction const& constraint, SatisfactionPrescription prescription);
double get_rho(double chi, size_t N, double volume, SatisfactionPrescription prescription);
double volume(Vector<FloatDPBounds> const& bnds);
Vector<FloatDPBounds> widen(Vector<FloatDPBounds> const& bx, double chi);
Vector<FloatDPBounds> shrink(Vector<FloatDPBounds> const& bx, double chi);

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
    EvaluationSequence(List<TimedVolume> const& tv, Vector<HardConstraintPrescription> const& usages);
  public:

    //! \brief The number of constraints inserted during creation
    size_t number_of_constraints() const;

    //! \brief timed volume accessor by index
    TimedVolume const& at(size_t idx) const;

    //! \brief Get the volume with time closest to \a time
    double const& near(double time) const;

    //! \brief The usage given the constraint index \a idx
    HardConstraintPrescription const& usage(ConstraintIndexType idx) const;

    //! \brief The number of elements in the sequence
    size_t size() const;

    friend ostream& operator<<(ostream& os, EvaluationSequence const& es);

  private:
    List<TimedVolume> _sequence;
    Vector<HardConstraintPrescription> _prescriptions;
};

class EvaluationSequenceBuilder {
  public:

    EvaluationSequenceBuilder(size_t N, Vector<EffectiveScalarMultivariateFunction> const& hs);

    void add_from(LabelledEnclosure const& e);

    void add(TimedBoxEvaluation const& tbe);

    EvaluationSequence build() const;

  private:
    size_t _list_index(double time) const;

  private:
    size_t const _N;
    size_t const _M;
    Vector<EffectiveScalarMultivariateFunction> const _hs;
    List<TimedBoxEvaluation> _timed_box_evaluations;

    Vector<SatisfactionPrescription> _prescriptions;
    Vector<Map<SatisfactionPrescription,double>> _max_robustness_false;
    Vector<Map<SatisfactionPrescription,double>> _t_false;
};

} // namespace Ariadne

#endif // ARIADNE_VERIFICATION