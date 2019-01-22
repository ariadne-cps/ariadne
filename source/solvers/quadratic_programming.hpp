/***************************************************************************
 *            quadratic programming.hpp
 *
 *          Copyright  2018  Nicola Dessì
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

/*! \file quadratic_programming.hpp
 *  \brief Quadratic programming.
 */

#ifndef ARIADNE_QUADRATIC_PROGRAMMING_HPP
#define ARIADNE_QUADRATIC_PROGRAMMING_HPP

#include "../algebra/matrix.hpp"
#include "../algebra/vector.hpp"
#include "../numeric/numeric.hpp"
#include "../output/logging.hpp"
#include "../utility/tuple.hpp"
#include "linear_programming.hpp"

#include "AlgebraAddOn.hpp"

namespace Ariadne {

template <class X> class Vector;
template <class X> class Matrix;
template <class X> class Affine;

enum class QuadraticProgramStatus : std::uint8_t {
  INDETERMINATE_FEASIBILITY = 0,
  PRIMAL_FEASIBLE = 1,
  DUAL_FEASIBLE = 2,
  PRIMAL_DUAL_FEASIBLE = 3,
  DEGENERATE_FEASIBILITY = 4
};

// class DegenerateFeasibilityProblemException : public std::runtime_error
// {
//   public:
//     DegenerateFeasibilityProblemException() : std::runtime_error("") {}
//     DegenerateFeasibilityProblemException(const StringType &what)
//         : std::runtime_error(what)
//     {
//     }
// };

class InfeasibleStartingPoint : std::runtime_error {
public:
  InfeasibleStartingPoint(const StringType &what) : std::runtime_error(what) {}
};

struct SingularQuadraticProgram : std::runtime_error {
  SingularQuadraticProgram(const StringType &what) : std::runtime_error(what){};
};

struct UnboundedQuadraticProgram : std::runtime_error {
  UnboundedQuadraticProgram(const StringType &what)
      : std::runtime_error(what) {}
};

struct InfeasibleQuadraticProgram : std::runtime_error {
  InfeasibleQuadraticProgram(const StringType &what)
      : std::runtime_error(what) {}
};

//! \ingroup OptimisationModule
//! Solver for quadratic programming problems using active set methods.
class ASMQPSolver : public Loggable {
  //! \brief Structure used in algorithm to store highly used variables and
  //! constants
  struct StepData;

public:
  //! \brief Find approximate optimal solution of \f$\min 1/2x^T Q x + x^T d
  //! \text{ s.t. } l<=Ax<=u; x\geq0\f$. Returns the pair (x,y) where x is the
  //! optimal point, and y the corresponding dual feasible point.
  //! @param Q RawFloatMatrix
  //! @param d RawFloatVector
  //! @param xl Lower bounds on x. This is considered only if @param x_bounds
  //! is set to true
  //! @param xu Upper bounds on x. This is considered only if @param x_bounds
  //! is set to true
  //! @param A Matrix<FloatDP>
  //! @param l RawFloatVector
  //! @param u RawFloatVector
  //! @param status Determines the status of solver when exit. It is a
  //! reference used on an upper level i.e. SQP non linear programming. Only 2
  //! status are catched here: -2 = QP did not converged in k_max steps; -3 =
  //! QP failed due to singularity.
  //! @param x_bounds Determines if xl or xu are used as bounds. If x_bounds =
  //! false the bounds are not set and starting point x is equal to xl (in
  //! case x_bounds = false, xl=xu. If x_bounds = true then starting point x
  //! is the middle between xl and xu, and they are used as bounds (on
  //! constraints)
  Tuple<FloatDP, Vector<FloatDP>, Vector<FloatDP>>
  minimise(const RawFloatMatrix &Q, const RawFloatVector &d,
           const RawFloatVector &xl, const RawFloatVector &xu,
           const Matrix<FloatDP> &A, const RawFloatVector &l,
           const RawFloatVector &u, int &status,
           const RawFloatVector x0 = Vector<FloatDP>()) const;

  //! \brief Verify if the problem is feasible for a such x in Ax=a and Bx>=b
  //! with a rtol as tolerance
  bool feasible(const RawFloatMatrix &A, const RawFloatVector &a,
                const RawFloatMatrix &B, const RawFloatVector &b,
                const RawFloatVector &x, const FloatDP &rtol) const;

  //! \brief Implement the phase I where a feasible starting point is found. It
  //! uses glpk library.
  void feasible_hotstart(Vector<FloatDP> &x, const Matrix<FloatDP> &A,
                         const Vector<FloatDP> &a, const Matrix<FloatDP> &B,
                         const Vector<FloatDP> &b, const FloatDP &rtol) const;
  //! \brief Normalize the problem to Ax=a and Bx>=b unifying all the bounds
  Tuple<Matrix<FloatDP>, Vector<FloatDP>, Matrix<FloatDP>, Vector<FloatDP>>
  normalize_problem(const Matrix<FloatDP> &A, const Vector<FloatDP> &A_lb,
                    const Vector<FloatDP> &A_ub, const Vector<FloatDP> &x,
                    const Vector<FloatDP> &x_lb,
                    const Vector<FloatDP> &x_ub) const;

private:
  const Vector<FloatDP> EMPTY_VEC = Vector<FloatDP>();

  //! \brief Internal minimise function wich take serialized problem inside v
  //! StepData structure
  Tuple<FloatDP, Vector<FloatDP>, Vector<FloatDP>> _minimise(struct StepData &v,
                                                             int &status) const;

  //! Initialize step data and include phase I using feasible_hotstart function
  void initialize_step_data(const Matrix<FloatDP> &H, const Vector<FloatDP> &d,
                            const Matrix<FloatDP> &A, const Vector<FloatDP> &a,
                            const Matrix<FloatDP> &B, const Vector<FloatDP> &b,
                            const Vector<FloatDP> &x, struct StepData &v) const;

  //! Linesearch to find
  FloatDP linesearch(struct StepData &v) const;

  //! \brief Perform a step of the optimization of \f$\min 1/2x^T Q x + x^T d
  //! \text{ s.t. } Ax>=l; x_l \leq x\leq x_u\f$. Returns true if a full ASMQP
  //! step is taken. In this case, the problem remains feasible (up to
  //! roundoff error).
  //! @param v Serialized StepData structure
  QuadraticProgramStatus _minimisation_step(struct StepData &v) const;

  //! \brief null space method to solve the eq problem relative to KKT condition
  //! system
  void null_space(struct StepData &v) const;
};

} // namespace Ariadne

#endif
