/***************************************************************************
 *            nonlinear_programming.hpp
 *
 *  Copyright 2009-17  Pieter Collins
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

/*! \file nonlinear_programming.hpp
 *  \brief Nonlinear programming.
 */

#ifndef ARIADNE_NONLINEAR_PROGRAMMING_HPP
#define ARIADNE_NONLINEAR_PROGRAMMING_HPP

#include "../utility/declarations.hpp"

#include "../numeric/numeric.hpp"
#include "../output/logging.hpp"
#include "../utility/tuple.hpp"

#include "../solvers/quadratic_programming.hpp"
namespace Ariadne {

template<class X, class R>
class Constraint;

class InfeasibleProblemException : public std::runtime_error
{
public:
  InfeasibleProblemException()
    : std::runtime_error("InfeasibleProblemException")
  {}
};
class IndeterminateFeasibilityException : public std::runtime_error
{
public:
  IndeterminateFeasibilityException()
    : std::runtime_error("IndeterminateFeasibilityException")
  {}
};
class DegenerateNonlinearFeasibilityProblemException : public std::runtime_error
{
public:
  DegenerateNonlinearFeasibilityProblemException()
    : std::runtime_error("DegenerateNonlinearFeasibilityProblemException")
  {}
};
class NearBoundaryOfFeasibleDomainException : public std::runtime_error
{
public:
  NearBoundaryOfFeasibleDomainException()
    : std::runtime_error("NearBoundaryOfFeasibleDomainException")
  {}
};

// ND -- 22/12/2019 -- Sequential quadratic programming -- BEGIN
#if defined HAVE_EIGEN3_H && defined HAVE_GLPK_H

class SQPStatusLeq1Exception : public std::runtime_error
{
public:
  int st;
  SQPStatusLeq1Exception(int status)
    : std::runtime_error("SQPStatusLeq1Exception")
  {
    st = status;
  }
};
class BFGSException : public std::runtime_error
{
public:
  BFGSException()
    : std::runtime_error("BFGS failed")
  {}
};
// ND -- 28/12/2018 -- Sequential quadratic programming -- END

#endif
//! \ingroup OptimisationModule EvaluationModule
//! Interface for nonlinear programming solvers.
class OptimiserInterface
{
public:
  // ND -- 28/12/2018 -- Developing dynamic solver -- BEGIN
  std::string problem_name;
  unsigned complexity_order;
  bool active_set_high = false;
  Vector<FloatDP> initial_guess;
  bool use_initial_guess = false;
  // ND -- 28/12/2018 -- Developing dynamic solver -- END

  //! \brief Virtual destructor.
  virtual ~OptimiserInterface() = default;
  //! \brief Create a dynamically-allocated copy.
  virtual OptimiserInterface* clone() const = 0;

  //! \brief Solve the general nonlinear programming problem \f$\min f(x) \text{
  //! such that } x\in D \text{ and } g(x)\in C\f$.
  virtual Vector<ValidatedNumericType> minimise(
    ValidatedScalarMultivariateFunction f,
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C) const = 0;
  //! \brief Solve the standard nonlinear programming problem \f$\min f(x)
  //! \text{ such that } x\in ,D\ g(x)\leq 0  \text{ and } h(x) = 0\f$.
  virtual Vector<ValidatedNumericType> minimise(
    ValidatedScalarMultivariateFunction f,
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ValidatedVectorMultivariateFunction h) const = 0;

  //! \brief Tests is the general nonlinear feasibility problem \f$x\in D \text{
  //! and } g(x)\in C\f$ is feasible.
  virtual ValidatedKleenean feasible(ExactBoxType D,
                                     ValidatedVectorMultivariateFunction g,
                                     ExactBoxType C) const = 0;
  //! \brief Tests is the standard nonlinear feasibility problem \f$x\in D,\
  //! g(x)\leq 0 \text{ and } h(x) = 0\f$ is feasible. Assumes \fD\f$ is
  //! singleton with nonempty interior. \internal This is one of the simplest
  //! nonlinear programming problems, and is a good test case for new
  //! algorithms.
  virtual ValidatedKleenean feasible(
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ValidatedVectorMultivariateFunction h) const = 0;

  //! \brief Tests if the point \a x is feasible, in that \f$x\in D\f$ and
  //! \f$g(x)\in N_\epsilon(C)\f$.
  virtual Bool almost_feasible_point(ExactBoxType D,
                                     ValidatedVectorMultivariateFunction g,
                                     ExactBoxType C,
                                     FloatApproximationVector x,
                                     FloatDPApproximation eps) const = 0;
  //! \brief Tests if the point \a x is feasible.
  virtual Bool is_feasible_point(ExactBoxType D,
                                 ValidatedVectorMultivariateFunction g,
                                 ExactBoxType C,
                                 ExactFloatVector x) const = 0;
  //! \brief Tests if the point \a x is near feasible.
  virtual Bool validate_feasibility(ExactBoxType D,
                                    ValidatedVectorMultivariateFunction g,
                                    ExactBoxType C,
                                    ExactVector x) const = 0;
  //! \brief Tests if the point \a x is near feasible, using approximate
  //! multipliers \a y to guide the search.
  virtual Bool validate_feasibility(ExactBoxType D,
                                    ValidatedVectorMultivariateFunction g,
                                    ExactBoxType C,
                                    ExactVector x,
                                    ExactVector y) const = 0;
  //! \brief Tests if the feasibility problem is definitely unsolvable, using
  //! multipliers \a y and local centering point \a x.
  virtual Bool validate_infeasibility(ExactBoxType D,
                                      ValidatedVectorMultivariateFunction g,
                                      ExactBoxType C,
                                      ExactVector x,
                                      ExactVector y) const = 0;
  //! \brief Tests if the box \a X definitely containss a feasible point.
  virtual ValidatedKleenean contains_feasible_point(
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C,
    FloatBoundsVector X) const = 0;
  //! \brief Tests if the Lagrange multipliers \a y are a certificate of
  //! infeasiblity.
  virtual Bool is_infeasibility_certificate(
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C,
    ExactVector y) const = 0;
};

//! \ingroup OptimisationModule
//! Common routines for nonlinear minimisation
class OptimiserBase
  : public OptimiserInterface
  , public Loggable
{
protected:
  static const FloatDPValue zero;
  static const FloatDPValue one;

public:
  virtual Vector<ValidatedNumericType> minimise(
    ValidatedScalarMultivariateFunction f,
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C) const = 0;
  virtual Vector<ValidatedNumericType> minimise(
    ValidatedScalarMultivariateFunction f,
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ValidatedVectorMultivariateFunction h) const;

  virtual ValidatedKleenean feasible(ExactBoxType D,
                                     ValidatedVectorMultivariateFunction g,
                                     ExactBoxType C) const = 0;
  virtual ValidatedKleenean feasible(
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ValidatedVectorMultivariateFunction h) const;

  virtual Bool almost_feasible_point(ExactBoxType D,
                                     ValidatedVectorMultivariateFunction g,
                                     ExactBoxType C,
                                     ApproximateVector x,
                                     FloatDPApproximation error) const;
  virtual Bool is_feasible_point(ExactBoxType D,
                                 ValidatedVectorMultivariateFunction g,
                                 ExactBoxType C,
                                 ExactVector x) const;
  virtual Bool validate_feasibility(ExactBoxType D,
                                    ValidatedVectorMultivariateFunction g,
                                    ExactBoxType C,
                                    ExactVector x) const;
  virtual Bool validate_feasibility(ExactBoxType D,
                                    ValidatedVectorMultivariateFunction g,
                                    ExactBoxType C,
                                    ExactVector x,
                                    ExactVector y) const;
  virtual Bool validate_infeasibility(ExactBoxType D,
                                      ValidatedVectorMultivariateFunction g,
                                      ExactBoxType C,
                                      ExactVector x,
                                      ExactVector y) const;
  virtual ValidatedKleenean contains_feasible_point(
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C,
    ValidatedVector X) const;
  virtual Bool is_infeasibility_certificate(
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C,
    ExactVector lambda) const;
};

//! \ingroup OptimisationModule
//! \brief Solver for feasibility problems based on a penalty-function approach
//!
//! For the feasibility problem \f$x\in D,\ g(x)\in C,\ h(x)=0\f$ where \f$D\f$
//! is singleton and \f$D,C\f$ have nonempty interiors. Introduce slack variable
//! \f$w=g(x)\f$, and minimise \f[ \sum_{j} (g_j(x)-w_j)^2 + \sum_k h_k(x)^2 \f]
//! with \f$x\in D\f$ and \f$w\in C\f$. Since the minimiser is not unique, we
//! need to add penalty terms \f$-\mu/2(\log(x_u-x)+log(x-x_l))\f$ and
//! \f$-\nu/2(\log(w_u-w)+log(w-w_l))\f$. It suffices to add a penalty in
//! \f$x\f$.
class PenaltyFunctionOptimiser : public OptimiserBase
{
public:
  using OptimiserBase::feasible;
  using OptimiserBase::minimise;

public:
  virtual PenaltyFunctionOptimiser* clone() const;
  virtual ValidatedKleenean check_feasibility(
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C,
    ExactVector x,
    ExactVector y) const;
  virtual Vector<ValidatedNumericType> minimise(
    ValidatedScalarMultivariateFunction f,
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C) const;
  virtual ValidatedKleenean feasible(ExactBoxType D,
                                     ValidatedVectorMultivariateFunction g,
                                     ExactBoxType C) const;
  virtual Void feasibility_step(const ExactBoxType& D,
                                const ApproximateVectorMultivariateFunction& g,
                                const ExactBoxType& C,
                                FloatApproximationVector& x,
                                FloatApproximationVector& w,
                                FloatDPApproximation& mu) const;
  virtual Void feasibility_step(const ExactBoxType& D,
                                const ValidatedVectorMultivariateFunction& g,
                                const ExactBoxType& C,
                                FloatBoundsVector& x,
                                FloatBoundsVector& w) const;
  virtual Void feasibility_step(const ExactBoxType& D,
                                const ApproximateVectorMultivariateFunction& g,
                                const ExactBoxType& C,
                                FloatApproximationVector& x,
                                FloatApproximationVector& y,
                                FloatApproximationVector& z) const;
};

//! \ingroup OptimisationModule
//! \brief Solver for linear programming problems using invalid interior point
//! methods. \details Introduces variables \f$w\f$ and attempts to find \f$x\in
//! D\f$ and \f$w\in C\f$ such that \f$g(x)=w\f$.
//!   The dual variables \f$y\f$ are unconstrained Lagrange multipliers for
//!   \f$y\cdot(g(x)-w)=0\f$.
class NonlinearInfeasibleInteriorPointOptimiser : public OptimiserBase
{
public:
  virtual NonlinearInfeasibleInteriorPointOptimiser* clone() const
  {
    return new NonlinearInfeasibleInteriorPointOptimiser(*this);
  }

  using OptimiserBase::feasible;
  using OptimiserBase::minimise;

  struct PrimalDualData;
  struct StepData;

  //! \brief Compute a \em local optimum of linear programming problem \f$\max
  //! f(x) \text{ such that } x\in D, g(x)\in C \text{ and } h(x)=0.\f$.
  //! \precondition The domain \f$D\f$ is singleton and has nonempty interior,
  //! and the codomain \f$C\f$ is nonempty. \return A box \f$X\f$ which
  //! definitely contains a feasible point, and contains a local optimum.
  virtual Vector<ValidatedNumericType> minimise(
    ValidatedScalarMultivariateFunction f,
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C) const;
  //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and }
  //! g(x)\in C\f$ is feasible.
  virtual ValidatedKleenean feasible(ExactBoxType D,
                                     ValidatedVectorMultivariateFunction g,
                                     ExactBoxType C) const;

  //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in
  //! D\f$ using a nonlinear feasibility test, hotstarting the method with the
  //! overall primal and dual variables.
  Pair<ValidatedKleenean, FloatApproximationVector> feasible_hotstarted(
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C,
    const PrimalDualData& wxy0) const;

  Void setup_feasibility(const ExactBoxType& D,
                         const ApproximateVectorMultivariateFunction& g,
                         const ExactBoxType& C,
                         StepData& stp) const;
  Void step(const ApproximateScalarMultivariateFunction& f,
            const ExactBoxType& D,
            const ApproximateVectorMultivariateFunction& g,
            const ExactBoxType& C,
            StepData& stp) const;
  //    FloatDPApproximation compute_mu(const ExactBoxType& D, const
  //    ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
  //                     FloatApproximationVector& w, FloatApproximationVector&
  //                     x, FloatApproximationVector& y) const;
};

//! \ingroup OptimisationModule
//! Solver for linear programming problems using interior point methods.
class NonlinearInteriorPointOptimiser : public OptimiserBase
{
public:
  virtual NonlinearInteriorPointOptimiser* clone() const
  {
    return new NonlinearInteriorPointOptimiser(*this);
  }

  using OptimiserBase::feasible;
  using OptimiserBase::minimise;

  //! \brief Compute a \em local optimum of linear programming problem \f$\max
  //! f(x) \text{ such that } x\in D, g(x)\in C \text{ and } h(x)=0.\f$.
  //! \precondition The domain \f$D\f$ is singleton and has nonempty interior,
  //! and the codomain \f$C\f$ is nonempty. \return A box \f$X\f$ which
  //! definitely contains a feasible point, and contains a local optimum.
  virtual Vector<ValidatedNumericType> minimise(
    ValidatedScalarMultivariateFunction f,
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C) const;
  //! \brief Tests is the nonlinear programming problem \f$x\in D, g(x)\in C
  //! \text{ and } h(x)= 0 \f$ is feasible.
  virtual ValidatedKleenean feasible(ExactBoxType D,
                                     ValidatedVectorMultivariateFunction g,
                                     ExactBoxType C) const;

  //! \brief Test if the constraints \f$g(y)\in C\f$ are solvable for \f$y\in
  //! D\f$ using a nonlinear feasibility test, hotstarting the method with the
  //! overall constraint violation, primal and dual variables.
  Pair<ValidatedKleenean, FloatApproximationVector> feasible_hotstarted(
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C,
    const FloatApproximationVector& x0,
    const FloatApproximationVector& lambda0,
    const FloatDPApproximation& violation0) const;

  //! \brief Test if the constraints \f$g(y)\in C\f$ are solvable for \f$y\in
  //! D\f$ using a nonlinear feasibility test, hotstarting the method with the
  //! primal and dual.
  Pair<ValidatedKleenean, FloatApproximationVector> feasible_hotstarted(
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C,
    const FloatApproximationVector& x0,
    const FloatApproximationVector& lambda0) const;

  Void minimisation_step(const ApproximateScalarMultivariateFunction& f,
                         const ExactBoxType& D,
                         const ApproximateVectorMultivariateFunction& g,
                         const ExactBoxType& C,
                         const ApproximateVectorMultivariateFunction& h,
                         FloatApproximationVector& x,
                         FloatApproximationVector& w,
                         FloatApproximationVector& kappa,
                         FloatApproximationVector& lambda,
                         const FloatDPApproximation& mu) const;
  Void setup_feasibility(const ExactBoxType& D,
                         const ApproximateVectorMultivariateFunction& g,
                         const ExactBoxType& C,
                         FloatApproximationVector& x,
                         FloatApproximationVector& lambda) const;
  Void setup_feasibility(const ExactBoxType& D,
                         const ApproximateVectorMultivariateFunction& g,
                         const ExactBoxType& C,
                         FloatApproximationVector& x,
                         FloatApproximationVector& lambda,
                         FloatDPApproximation& t) const;
  Void feasibility_step(const ExactBoxType& D,
                        const ApproximateVectorMultivariateFunction& g,
                        const ExactBoxType& C,
                        FloatApproximationVector& x,
                        FloatApproximationVector& lambda) const;
  Void feasibility_step(const ExactBoxType& D,
                        const ApproximateVectorMultivariateFunction& g,
                        const ExactBoxType& C,
                        FloatApproximationVector& x,
                        FloatApproximationVector& lambda,
                        FloatDPApproximation& violation) const;
  Void initialise_lagrange_multipliers(
    const ExactBoxType& D,
    const ApproximateVectorMultivariateFunction& g,
    const ExactBoxType& C,
    const FloatApproximationVector& x,
    FloatApproximationVector& lambda) const;
  FloatDPApproximation compute_mu(
    const ExactBoxType& D,
    const ApproximateVectorMultivariateFunction& g,
    const ExactBoxType& C,
    const FloatApproximationVector& x,
    const FloatApproximationVector& lambda) const;

public: // Deprecated
  Void compute_tz(const ExactBoxType& D,
                  const ApproximateVectorMultivariateFunction& g,
                  const ExactBoxType& C,
                  FloatApproximationVector& x,
                  FloatDPApproximation& t,
                  FloatApproximationVector& z) const
  {}
  Void feasibility_step(const ExactBoxType& D,
                        const ApproximateVectorMultivariateFunction& g,
                        const ExactBoxType& C,
                        FloatApproximationVector& x,
                        FloatApproximationVector& y,
                        FloatApproximationVector& z,
                        FloatDPApproximation& violation) const {};
  Void linearised_feasibility_step(
    const ExactBoxType& D,
    const ApproximateVectorMultivariateFunction& g,
    const ExactBoxType& C,
    FloatDPApproximation& slack,
    FloatApproximationVector& x,
    FloatApproximationVector& lambda) const {};
  Void linearised_feasibility_step(
    const ExactBoxType& D,
    const ApproximateVectorMultivariateFunction& g,
    const ExactBoxType& C,
    FloatApproximationVector& x,
    FloatApproximationVector& y,
    FloatApproximationVector& z,
    FloatDPApproximation& t) const {};

private:
  FloatDPApproximation compute_mu(
    const ApproximateScalarMultivariateFunction& f,
    const ExactBoxType& D,
    const ValidatedVectorMultivariateFunction& g,
    const ExactBoxType& C,
    const FloatApproximationVector& x,
    const FloatApproximationVector& y) const;
  Void compute_violation(const ExactBoxType& D,
                         const ApproximateVectorMultivariateFunction& g,
                         const ExactBoxType& C,
                         FloatApproximationVector& x,
                         FloatDPApproximation& t) const;
};

class IntervalOptimiser : public NonlinearInteriorPointOptimiser
{
  virtual IntervalOptimiser* clone() const
  {
    return new IntervalOptimiser(*this);
  }
  virtual ValidatedKleenean feasible_zero(
    ExactBoxType D,
    ValidatedVectorMultivariateFunction h) const;
  Void feasibility_step(const ExactFloatVector& xl,
                        const ExactFloatVector& xu,
                        const ValidatedVectorMultivariateFunction& h,
                        FloatBoundsVector& x,
                        FloatBoundsVector& y,
                        FloatBoundsVector& zl,
                        FloatBoundsVector zu,
                        FloatDPBounds& mu) const;
};

class ApproximateOptimiser : public NonlinearInteriorPointOptimiser
{
  virtual ApproximateOptimiser* clone() const
  {
    return new ApproximateOptimiser(*this);
  }
  virtual ValidatedKleenean feasible_zero(
    ExactBoxType D,
    ValidatedVectorMultivariateFunction h) const;
  Void feasibility_step(const ExactBoxType& D,
                        const ApproximateVectorMultivariateFunction& h,
                        FloatApproximationVector& X,
                        FloatApproximationVector& Lambda) const;
};

/*//! \ingroup OptimisationModule
//! Solver for linear programming problems using interior point methods.
//! WARNING: This class currently does not work; maybe there is a problem with
the algorithms. class KrawczykOptimiser : public OptimiserBase
{

  public:
    virtual KrawczykOptimiser* clone() const { return new
KrawczykOptimiser(*this); }

    //! \brief Solve the linear programming problem \f$\max f(x) \text{ such
that } x\in D \text{ and } g(x)\in C\f$. virtual Vector<ValidatedNumericType>
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D,
ValidatedVectorMultivariateFunction g, ExactBoxType C) const;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and }
g(x)\in C\f$ is feasible. virtual ValidatedKleenean feasible(ExactBoxType D,
ValidatedVectorMultivariateFunction g, ExactBoxType C) const;

  public:
    //! \brief Try to solve the nonlinear constraint problem by applying the
Krawczyk contractor to the Kuhn-Tucker conditions,
    //! hotstarting the iteration with the primal and dual variables.
    ValidatedKleenean minimise(ValidatedScalarMultivariateFunction f,
ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C, const
FloatDPBounds& t0, const FloatBoundsVector& x0, const FloatBoundsVector& y0,
const FloatBoundsVector& z0) const;

    //! \brief A primal-dual feasibility step for the problem \f$g(y)\in C;\
y\in D\f$. Void minimisation_step(const ExactBoxType& D, const
ValidatedVectorMultivariateFunction& g, const ExactBoxType& C,
                           FloatBoundsVector& x, FloatBoundsVector& y,
FloatBoundsVector& z, FloatDPBounds& t) const;
    //! \brief A primal-dual feasibility step for the problem \f$g(y)\in C;\
y\in D\f$. Void feasibility_step(const ExactBoxType& D, const
ValidatedVectorMultivariateFunction& g, const ExactBoxType& C,
                          FloatBoundsVector& x, FloatBoundsVector& y,
FloatBoundsVector& z, FloatDPBounds& t) const;

    //! \brief A primal feasibility step for the problem \f$g(y)\in C;\ y\in
D\f$. \deprecated Void feasibility_step(const ExactBoxType& D, const
ValidatedVectorMultivariateFunction& g, const ExactBoxType& C,
                          FloatBoundsVector& y, FloatDPBounds& t) const;
    //! \brief A feasibility step for the problem \f$g(y)\leq 0\f$. \deprecated
    Void feasibility_step(const ValidatedVectorMultivariateFunction& g,
                          FloatBoundsVector& x, FloatBoundsVector& y,
FloatBoundsVector& z, FloatDPBounds& t) const;
    //! \brief An optimization step for the problem \f$\max f(y) \text{ s.t. }
g(y)\leq 0\f$. \deprecated Void minimisation_step(const
ValidatedScalarMultivariateFunction& f, const
ValidatedVectorMultivariateFunction& g, FloatBoundsVector& x, FloatBoundsVector&
y, FloatBoundsVector& z) const; protected: Void setup_feasibility(const
ExactBoxType& D, const ValidatedVectorMultivariateFunction& g, const
ExactBoxType& C, FloatBoundsVector& x, FloatBoundsVector& y, FloatBoundsVector&
z, FloatDPBounds& t) const; protected: Void compute_tz(const ExactBoxType& D,
const ValidatedVectorMultivariateFunction& g, const ExactBoxType& C, const
FloatBoundsVector& y, FloatDPBounds& t, FloatBoundsVector& z) const;
};

*/

//-------------------------------------------------------------------------------

// ND -- 22/12/2019 -- Sequential quadratic programming -- BEGIN
// need eigen3 to run
#if defined HAVE_EIGEN3_H && defined HAVE_GLPK_H

class NonlinearMixedOptimiser;
//! \ingroup OptimisationModule
//! Solver for non linear programming problems using sequential quadratic
//! programming.
class NonlinearSQPOptimiser : public OptimiserBase
{
  friend NonlinearMixedOptimiser;

public:
  virtual NonlinearSQPOptimiser* clone() const
  {
    return new NonlinearSQPOptimiser(*this);
  }

  NonlinearSQPOptimiser();

  using OptimiserBase::feasible;
  using OptimiserBase::minimise;

  //! Struct to serialize data used by algorithm
  struct StepData;
  //! Enum for exit status
  enum class SQPStatus : int8_t
  {
    OK = 1,
    NOTHING = 0,
    NULL_STEP = -1,
    QP_NOT_CONV = -2,
    QP_SINGULAR = -3,
    QP_INFEASIBLE = -4,
    NULL_GRADIENT = -5,
    UNSAFE_ZONE = -6,
    B_SINGULAR = -7,
    BFGS_FAILED = -8,
    UNSAFE_AND_NULL_GRADIENT = -9
  };

  //! \brief Compute a \em local optimum of linear programming problem \f$\max
  //! f(x) \text{ such that } x\in D, g(x)\in C \text{ and } h(x)=0.\f$.
  //! \precondition The domain \f$D\f$ is singleton and has nonempty interior,
  //! and the codomain \f$C\f$ is nonempty. \return A box \f$X\f$ which
  //! definitely contains a feasible point, and contains a local optimum.
  //! @param f objective function
  //! @param D boundaries on x
  //! @param g constraints function
  //! @param C boundaries on g
  virtual Vector<ValidatedNumericType> minimise(
    ValidatedScalarMultivariateFunction f,
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C) const;

  virtual ValidatedKleenean feasible(ExactBoxType D,
                                     ValidatedVectorMultivariateFunction g,
                                     ExactBoxType C) const;

  virtual ValidatedKleenean check_feasibility(
    const ExactBoxType& d,
    const ValidatedVectorMultivariateFunction& f,
    const ExactBoxType& c,
    const ExactPoint& y) const;

  bool feasible_point(const ExactBoxType domain,
                      const ValidatedVectorMultivariateFunction g,
                      const ExactBoxType codomain,
                      RawFloatVector& point) const;

  //! Take a step on a descendent direction of function f s.t. g
  Void feasibility_step(const ExactBoxType& D,
                        const ApproximateVectorMultivariateFunction& g,
                        const ExactBoxType& C,
                        StepData& stp) const;

  //! Take a step on a descendent direction of function f s.t. g
  Void step(const ApproximateScalarMultivariateFunction& f,
            const ExactBoxType& D,
            const ApproximateVectorMultivariateFunction& g,
            const ExactBoxType& C,
            StepData& stp) const;

private:
  //! Compute Quasi-Newton formula to update an approximation of the Hessian
  //! @param grd_f_cap
  //! @param grd_g_cap
  //! @param s direction vector s (computed as x_new - x)
  //! @param v Algorithm data serialized
  Void BFGS(const Vector<FloatDP>& x_cap,
            const Vector<FloatDP>& grd_f_cap,
            const Matrix<FloatDP>& grd_g_cap,
            const struct StepData& v,
            Matrix<FloatDP>& B) const;

  //! Method to get the best step length in the right direction. To find the
  //! best step length alpha is used a merit function.
  //! @see merit_function_l1
  //! @param x current position
  //! @param p vector direction used as x_new = x + (step-length)*p
  //! @param grd_f gradient of f on x computed last cycle
  //! @param g constraints function
  //! @param g_l lower boundaries on g
  //! @param g_u upper boundaries on g
  //! @param lambda Lambda (full-)multipliers computed considering all
  //! constraints and boundaries
  FloatDP linesearch(
    struct StepData& v,
    const std::function<FloatDP(const Vector<FloatDP>&, const FloatDP&)>&
      phi_l1) const;

  void initialize_step_data(ApproximateScalarMultivariateFunction f,
                            ExactBoxType D,
                            ApproximateVectorMultivariateFunction g,
                            ExactBoxType C,
                            struct StepData& v) const;

  //! Element x Element product between vectors
  Matrix<FloatDP> prodVec(const Vector<FloatDP>& v1,
                          const Vector<FloatDP>& v2T) const;

  const Vector<FloatDP> EMPTY_VEC = Vector<FloatDP>();
  std::shared_ptr<ASMQPSolver> qpsolver_ptr; //< Sub-solver QP as shared pointer
};
// ND -- 22/12/2019 -- Sequential quadratic programming -- END
//-------------------------------------------------------------------------------

// ND -- 22/12/2019 -- Mixed quadratic programming -- BEGIN
//! \ingroup OptimisationModule
//! Solver for non linear programming problems using sequential quadratic
//! programming.
class NonlinearMixedOptimiser : public OptimiserBase
{
public:
  virtual NonlinearMixedOptimiser* clone() const
  {
    return new NonlinearMixedOptimiser(*this);
  }

  NonlinearMixedOptimiser();

  using OptimiserBase::feasible;
  using OptimiserBase::minimise;

  //! Step data to store infos
  struct StepData;

  //! Type of strategy to use
  enum class Strategy : uint8_t
  {
    SQP = 0,
    IPM = 1,
    BOTH = 2
  };

  //! \brief Compute a \em local optimum of linear programming problem \f$\max
  //! f(x) \text{ such that } x\in D, g(x)\in C \text{ and } h(x)=0.\f$.
  //! \precondition The domain \f$D\f$ is singleton and has nonempty interior,
  //! and the codomain \f$C\f$ is nonempty. \return A box \f$X\f$ which
  //! definitely contains a feasible point, and contains a local optimum.
  virtual Vector<ValidatedNumericType> minimise(
    ValidatedScalarMultivariateFunction f,
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C) const;

  //! \brief Compute a \em local optimum of nonlinear programming problem using
  //! static policy
  virtual Vector<ValidatedNumericType> minimise_static(
    ValidatedScalarMultivariateFunction f,
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C) const;

  //! \brief Compute a \em local optimum of nonlinear programming problem using
  //! dynamic policy
  virtual Vector<ValidatedNumericType> minimise_dynamic(
    ValidatedScalarMultivariateFunction f,
    ExactBoxType D,
    ValidatedVectorMultivariateFunction g,
    ExactBoxType C) const;

  //! \brief Search a feasible point
  virtual ValidatedKleenean feasible(ExactBoxType D,
                                     ValidatedVectorMultivariateFunction g,
                                     ExactBoxType C) const;

private:
  //! \brief Minimise using a strategy (should be used in dynamic minimise)
  Void minimise_with_strategy(const ApproximateScalarMultivariateFunction& f,
                              const ExactBoxType& d,
                              const ApproximateVectorMultivariateFunction& g,
                              const ExactBoxType& c,
                              const ApproximateVectorMultivariateFunction& h,
                              struct StepData& v) const;

  //! \brief Single step of minimisation using strategy
  Void step(const ApproximateScalarMultivariateFunction& f,
            const ExactBoxType& d,
            const ApproximateVectorMultivariateFunction& g,
            const ExactBoxType& c,
            const ApproximateVectorMultivariateFunction& h,
            struct StepData& v) const;

  //! \brief Check if a change of strategy is needed and update the strategy
  //! variable inside SteData parameter
  Void change_strategy(struct StepData&) const;

  //! \brief Function to decide which policy is correct
  Strategy lookup_policy(struct StepData&) const;

  Void initialize_step_data(struct StepData& stepData,
                            const Strategy& strategy,
                            const ValidatedScalarMultivariateFunction f,
                            const ExactBoxType D,
                            const ValidatedVectorMultivariateFunction g,
                            const ExactBoxType C) const;

  std::shared_ptr<NonlinearInteriorPointOptimiser> nlipm_ptr;
  std::shared_ptr<NonlinearSQPOptimiser> nlsqp_ptr;
};

#endif
// ND -- 22/12/2019 -- Mixed quadratic programming -- END

} // namespace Ariadne

#endif
