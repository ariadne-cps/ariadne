/***************************************************************************
 *            inner_approximation.hpp
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

#ifndef ARIADNE_INNER_APPROXIMATION
#define ARIADNE_INNER_APPROXIMATION

#include "numeric/floatdp.hpp"
#include "function/constraint.hpp"
#include "dynamics/enclosure.hpp"
#include "dynamics/orbit.hpp"
#include "solvers/linear_programming.hpp"
#include "solvers/nonlinear_programming.hpp"
#include "glpk.h"

namespace Ariadne {

class LinearSolverInterface {
  public:
    virtual FloatDP minimise(SizeType i, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const = 0;
    virtual FloatDP maximise(SizeType i, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const = 0;

    virtual LinearSolverInterface* clone() const = 0;
    virtual ~LinearSolverInterface() = default;
};

class NativeLinearSolver : public LinearSolverInterface {
  public:
    FloatDP minimise(SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const;
    FloatDP maximise(SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const;
  protected:
    virtual FloatDP solve(SizeType k, Vector<FloatDP> const& c, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const = 0;
};

class NativeSimplex : public NativeLinearSolver {
  protected:
    FloatDP solve(SizeType k, Vector<FloatDP> const& c, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const override;
  public:
    LinearSolverInterface* clone() const override;
};

class NativeIPM : public NativeLinearSolver {
  protected:
    FloatDP solve(SizeType k, Vector<FloatDP> const& c, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const override;
  public:
    LinearSolverInterface* clone() const override;
};

class GLPKSolver : public LinearSolverInterface {
  public:
    FloatDP minimise(SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const override;
    FloatDP maximise(SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const override;
  protected:
    virtual void optimisation_method(glp_prob* lp) const = 0;
  private:
    FloatDP _solve(int dir, SizeType k, Matrix<FloatDP> const& A, Vector<FloatDP> const& b, Vector<FloatDP> const& xl, Vector<FloatDP> const& xu) const;
};

class GLPKSimplex : public GLPKSolver {
  protected:
    void optimisation_method(glp_prob* lp) const override;
  public:
    LinearSolverInterface* clone() const override;
};

class GLPKIPM : public GLPKSolver {
  protected:
    void optimisation_method(glp_prob* lp) const override;
  public:
    LinearSolverInterface* clone() const override;
};

List<LabelledEnclosure> boundary(LabelledEnclosure const& enclosure);

Tuple<Matrix<FloatDP>,Vector<FloatDP>,Vector<FloatDP>,Vector<FloatDP>> construct_parallel_linearisation_problem(ValidatedVectorMultivariateFunction const& f, ExactBoxType const& d);

ApproximateOptimiser::ValidatedNumericType nonlinear_intersection_bound(ValidatedVectorMultivariateFunctionPatch const& h, ExactBoxType const& d, SizeType var_idx, bool maximise);

ExactBoxType nonlinear_nonintersection_domain(ApproximateOptimiser::ValidatedNumericType const& bound, ExactBoxType const& d, SizeType var_idx, bool maximise, double scaling);

ExactBoxType inner_difference(ExactBoxType const& bx1, ExactBoxType const& bx2);

//! \brief Interface for searching an interval where upper values are expected to have a false outcome, lower values
//! to have a true outcome, with separated value sets
template<class T> class IntervalSearchInterface {
  public:
    //! \brief Return the current point
    virtual T const& current() const = 0;
    //! \brief Move to the next point according to the current outcome
    //! \return Whether a next point is available
    virtual bool move_next(bool current_outcome) = 0;
    //! \brief Return whether a next point is available
    virtual bool ended() const = 0;
    //! \brief Return whether one solution has been found
    //! \details Even if one solution has been found, the search might not have ended, in case we want to refine it
    virtual bool solution_found() const = 0;
};

template<class T> class DownwardSearch : public IntervalSearchInterface<T> {
  public:
    DownwardSearch(T const& lower, T const& upper, T const& decrement) : _current(upper), _lower(lower), _decrement(decrement), _solution_found(false) { }
    T const& current() const override { return _current; }
    bool move_next(bool current_outcome) override {
        if (current_outcome) {
            _solution_found = true;
            _current = _lower;
            return false;
        } else if (not ended()) {
            _current -= _decrement;
            return true;
        } else {
            return false;
        }
    }

    bool ended() const override { return _current-_lower < _decrement; }
    bool solution_found() const override { return _solution_found; }
  private:
    T _current;
    T _lower;
    T const _decrement;
    bool _solution_found;
};

template<class T> class BisectionSearch : public IntervalSearchInterface<T> {
public:
    BisectionSearch(T const& lower, T const& upper, T const& threshold) :
        _current(upper), _lower(lower), _lower_is_set(false), _upper(upper), _threshold(threshold), _solution_found(false) { }
    T const& current() const override { return _current; }
    bool move_next(bool current_outcome) override {
        if (not ended()) {
            if (current_outcome) {
                _solution_found = true;
                _lower = _current;
                _lower_is_set = true;
                if (_current == _upper) {
                    return false;
                }
            } else {
                _upper = _current;
                if (_current == _lower) {
                    return false;
                }
            }
            _current = (_lower_is_set ? _lower + (_upper-_lower)/2 : _lower);
            _lower_is_set = true;
            return true;
        } else {
            return false;
        }
    }
    bool ended() const override { return _upper-_lower <= _threshold; }
    bool solution_found() const override { return _solution_found; }
private:
    T _current;
    T _lower;
    bool _lower_is_set;
    T _upper;
    T const _threshold;
    bool _solution_found;
};

//! \brief A class for contracting the domain of a function
class ContractorInterface {
  public:
    //! \brief Contract the domain \a d of a function \a f
    //! \return Returns the contracted box, the empty box if no feasible point is found anymore
    //! \details The result is guaranteed to be a subset of \a d
    virtual ExactBoxType contract(ValidatedVectorMultivariateFunctionPatch const& f, ExactBoxType const& d) const = 0;
    virtual ContractorInterface* clone() const = 0;
    virtual ~ContractorInterface() = default;
};

class ParallelLinearisationContractor : public ContractorInterface {
  private:
    ParallelLinearisationContractor(ParallelLinearisationContractor const& other);
  public:
    ParallelLinearisationContractor(LinearSolverInterface const& solver, SizeType num_iterations, SizeType split_depth = 0);
    ExactBoxType contract(ValidatedVectorMultivariateFunctionPatch const& f, ExactBoxType const& d) const override;
    ContractorInterface* clone() const override;
  private:
    std::shared_ptr<LinearSolverInterface> _solver_ptr;
    SizeType const _num_iterations;
    SizeType const _split_depth;
};

class InnerApproximatorInterface {
  public:
    virtual LabelledEnclosure compute_from(LabelledEnclosure const& outer) const = 0;
    virtual ListSet<LabelledEnclosure> compute_from(ListSet<LabelledEnclosure> const& outer_list) const = 0;
    virtual Orbit<LabelledEnclosure> compute_from(Orbit<LabelledEnclosure> const& outer_orbit) const = 0;
};

class InnerApproximatorBase : public InnerApproximatorInterface {
  public:
    InnerApproximatorBase(ContractorInterface const& contractor);
    virtual LabelledEnclosure compute_from(LabelledEnclosure const& outer) const override = 0;
  protected:
    ListSet<LabelledEnclosure> compute_from(ListSet<LabelledEnclosure> const& outer_list) const override;
    Orbit<LabelledEnclosure> compute_from(Orbit<LabelledEnclosure> const& outer_orbit) const override;
  protected:
    std::shared_ptr<ContractorInterface> _contractor_ptr;
};

class NonlinearCandidateValidationInnerApproximator : public InnerApproximatorBase {
  private:
    SizeType const _max_rounds;
  public:
    //! \brief Constructs from a \a contractor and a \a max_rounds for refinement, subject to the minimum given by the number of boundaries
    NonlinearCandidateValidationInnerApproximator(ContractorInterface const& contractor, SizeType max_rounds = std::numeric_limits<SizeType>::max());

    LabelledEnclosure compute_from(LabelledEnclosure const& outer) const override final;
    ListSet<LabelledEnclosure> compute_from(ListSet<LabelledEnclosure> const& outer_list) const override final { return InnerApproximatorBase::compute_from(outer_list); }
    Orbit<LabelledEnclosure> compute_from(Orbit<LabelledEnclosure> const& outer_orbit) const override final { return InnerApproximatorBase::compute_from(outer_orbit); }
};

} // namespace Ariadne

#endif // ARIADNE_INNER_APPROXIMATION