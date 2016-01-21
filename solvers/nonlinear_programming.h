/***************************************************************************
 *            nonlinear_programming.h
 *
 *  Copyright 2009  Pieter Collins
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

/*! \file nonlinear_programming.h
 *  \brief Nonlinear programming.
 */

#ifndef ARIADNE_NONLINEAR_PROGRAMMING_H
#define ARIADNE_NONLINEAR_PROGRAMMING_H

#include "utility/declarations.h"

#include "utility/logging.h"
#include "numeric/numeric.h"
#include "utility/tuple.h"


namespace Ariadne {

template<class X, class R> class Constraint;
typedef Constraint<EffectiveScalarFunction,EffectiveNumericType> EffectiveConstraint;
typedef Constraint<ValidatedScalarFunction,ValidatedNumericType> ValidatedConstraint;

class InfeasibleProblemException : public std::runtime_error {
  public: InfeasibleProblemException() : std::runtime_error("InfeasibleProblemException") { }
};
class IndeterminateFeasibilityException : public std::runtime_error {
  public: IndeterminateFeasibilityException() : std::runtime_error("IndeterminateFeasibilityException") { }
};
class DegenerateNonlinearFeasibilityProblemException : public std::runtime_error {
  public: DegenerateNonlinearFeasibilityProblemException() : std::runtime_error("DegenerateNonlinearFeasibilityProblemException") { }
};
class NearBoundaryOfFeasibleDomainException : public std::runtime_error {
  public: NearBoundaryOfFeasibleDomainException() : std::runtime_error("NearBoundaryOfFeasibleDomainException") { }
};

//! \ingroup OptimisationModule EvaluationModule
//! Interface for nonlinear programming solvers.
class OptimiserInterface {
  public:
    //! \brief Virtual destructor.
    virtual ~OptimiserInterface() { }
    //! \brief Create a dynamically-allocated copy.
    virtual OptimiserInterface* clone() const = 0;

    //! \brief Solve the general nonlinear programming problem \f$\min f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarFunction f, ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C) const = 0;
    //! \brief Solve the standard nonlinear programming problem \f$\min f(x) \text{ such that } x\in ,D\ g(x)\leq 0  \text{ and } h(x) = 0\f$.
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarFunction f, ExactBoxType D, ValidatedVectorFunction g, ValidatedVectorFunction h) const = 0;

    //! \brief Tests is the general nonlinear feasibility problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C) const = 0;
    //! \brief Tests is the standard nonlinear feasibility problem \f$x\in D,\ g(x)\leq 0 \text{ and } h(x) = 0\f$ is feasible. Assumes \fD\f$ is singleton with nonempty interior.
    //! \internal This is one of the simplest nonlinear programming problems, and is a good test case for new algorithms.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorFunction g, ValidatedVectorFunction h) const = 0;

    //! \brief Tests if the point \a x is feasible, in that \f$x\in D\f$ and \f$g(x)\in N_\epsilon(C)\f$.
    virtual Bool almost_feasible_point(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                       FloatApproximationVector x, Float64Approximation eps) const = 0;
    //! \brief Tests if the point \a x is feasible.
    virtual Bool is_feasible_point(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                   ExactFloatVector x) const = 0;
    //! \brief Tests if the point \a x is near feasible.
    virtual Bool validate_feasibility(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                      ExactVector x) const = 0;
    //! \brief Tests if the point \a x is near feasible, using approximate multipliers \a y to guide the search.
    virtual Bool validate_feasibility(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                      ExactVector x, ExactVector y) const = 0;
    //! \brief Tests if the feasibility problem is definitely unsolvable, using multipliers \a y and local centering point \a x.
    virtual Bool validate_infeasibility(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                      ExactVector x, ExactVector y) const = 0;
    //! \brief Tests if the box \a X definitely containss a feasible point.
    virtual ValidatedKleenean contains_feasible_point(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                            FloatBoundsVector X) const = 0;
    //! \brief Tests if the Lagrange multipliers \a y are a certificate of infeasiblity.
    virtual Bool is_infeasibility_certificate(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                              ExactVector y) const = 0;
};

//! \ingroup OptimisationModule
//! Common routines for nonlinear minimisation
class OptimiserBase
    : public OptimiserInterface
    , public Loggable
{
  public:
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarFunction f, ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C) const = 0;
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarFunction f, ExactBoxType D, ValidatedVectorFunction g, ValidatedVectorFunction h) const;

    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C) const = 0;
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorFunction g, ValidatedVectorFunction h) const;

    virtual Bool almost_feasible_point(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                       ApproximateVector x, Float64Approximation error) const;
    virtual Bool is_feasible_point(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                   ExactVector x) const;
    virtual Bool validate_feasibility(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                      ExactVector x) const;
    virtual Bool validate_feasibility(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                      ExactVector x, ExactVector y) const;
    virtual Bool validate_infeasibility(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                        ExactVector x, ExactVector y) const;
    virtual ValidatedKleenean contains_feasible_point(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                            ValidatedVector X) const;
    virtual Bool is_infeasibility_certificate(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                              ExactVector lambda) const;
};

//! \ingroup OptimisationModule
//! \brief Solver for feasibility problems based on a penalty-function approach
//!
//! For the feasibility problem \f$x\in D,\ g(x)\in C,\ h(x)=0\f$ where \f$D\f$ is singleton and \f$D,C\f$ have nonempty interiors.
//! Introduce slack variable \f$w=g(x)\f$, and minimise \f[ \sum_{j} (g_j(x)-w_j)^2 + \sum_k h_k(x)^2 \f] with \f$x\in D\f$ and \f$w\in C\f$.
//! Since the minimiser is not unique, we need to add penalty terms \f$-\mu/2(\log(x_u-x)+log(x-x_l))\f$ and \f$-\nu/2(\log(w_u-w)+log(w-w_l))\f$.
//! It suffices to add a penalty in \f$x\f$.
class PenaltyFunctionOptimiser
    : public OptimiserBase
{
  public:
    virtual PenaltyFunctionOptimiser* clone() const;
    virtual ValidatedKleenean check_feasibility(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C, ExactVector x, ExactVector y) const;
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarFunction f, ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C) const;
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C) const;
    virtual Void feasibility_step(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                                  FloatApproximationVector& x, FloatApproximationVector& w, Float64Approximation& mu) const;
    virtual Void feasibility_step(const ExactBoxType& D, const ValidatedVectorFunction& g, const ExactBoxType& C,
                                  FloatBoundsVector& x, FloatBoundsVector& w) const;
    virtual Void feasibility_step(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                                  FloatApproximationVector& x, FloatApproximationVector& y, FloatApproximationVector& z) const;
};




//! \ingroup OptimisationModule
//! \brief Solver for linear programming problems using invalid interior point methods.
//! \details Introduces variables \f$w\f$ and attempts to find \f$x\in D\f$ and \f$w\in C\f$ such that \f$g(x)=w\f$.
//!   The dual variables \f$y\f$ are unconstrained Lagrange multipliers for \f$y\cdot(g(x)-w)=0\f$.
class NonlinearInfeasibleInteriorPointOptimiser
    : public OptimiserBase
{
  public:
    virtual NonlinearInfeasibleInteriorPointOptimiser* clone() const { return new NonlinearInfeasibleInteriorPointOptimiser(*this); }

    using OptimiserBase::minimise;
    using OptimiserBase::feasible;

    struct PrimalDualData;
    struct StepData;

    //! \brief Compute a \em local optimum of linear programming problem \f$\max f(x) \text{ such that } x\in D, g(x)\in C \text{ and } h(x)=0.\f$.
    //! \precondition The domain \f$D\f$ is singleton and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A box \f$X\f$ which definitely contains a feasible point, and contains a local optimum.
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarFunction f, ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C) const;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C) const;

    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the overall primal and dual variables.
    Pair<ValidatedKleenean,FloatApproximationVector> feasible_hotstarted(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                                                const PrimalDualData& wxy0) const;

    Void setup_feasibility(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                           StepData& stp) const;
    Void step(const ApproximateScalarFunction& f, const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
              StepData& stp) const;
//    Float64Approximation compute_mu(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
//                     FloatApproximationVector& w, FloatApproximationVector& x, FloatApproximationVector& y) const;
};


//! \ingroup OptimisationModule
//! Solver for linear programming problems using interior point methods.
class NonlinearInteriorPointOptimiser
    : public OptimiserBase
{
  public:
    virtual NonlinearInteriorPointOptimiser* clone() const { return new NonlinearInteriorPointOptimiser(*this); }

    using OptimiserBase::minimise;
    using OptimiserBase::feasible;

    //! \brief Compute a \em local optimum of linear programming problem \f$\max f(x) \text{ such that } x\in D, g(x)\in C \text{ and } h(x)=0.\f$.
    //! \precondition The domain \f$D\f$ is singleton and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A box \f$X\f$ which definitely contains a feasible point, and contains a local optimum.
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarFunction f, ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C) const;
    //! \brief Tests is the nonlinear programming problem \f$x\in D, g(x)\in C \text{ and } h(x)= 0 \f$ is feasible.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C) const;

    //! \brief Test if the constraints \f$g(y)\in C\f$ are solvable for \f$y\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the overall constraint violation, primal and dual variables.
    Pair<ValidatedKleenean,FloatApproximationVector> feasible_hotstarted(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                                  const FloatApproximationVector& x0, const FloatApproximationVector& lambda0, const Float64Approximation& violation0) const;

    //! \brief Test if the constraints \f$g(y)\in C\f$ are solvable for \f$y\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal and dual.
    Pair<ValidatedKleenean,FloatApproximationVector> feasible_hotstarted(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                                                  const FloatApproximationVector& x0, const FloatApproximationVector& lambda0) const;

    Void minimisation_step(const ApproximateScalarFunction& f, const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C, const ApproximateVectorFunction& h,
                           FloatApproximationVector& x, FloatApproximationVector& w, FloatApproximationVector& kappa, FloatApproximationVector& lambda, const Float64Approximation& mu) const;
    Void setup_feasibility(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                           FloatApproximationVector& x, FloatApproximationVector& lambda) const;
    Void setup_feasibility(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                           FloatApproximationVector& x, FloatApproximationVector& lambda, Float64Approximation& t) const;
    Void feasibility_step(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                          FloatApproximationVector& x, FloatApproximationVector& lambda) const;
    Void feasibility_step(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                          FloatApproximationVector& x, FloatApproximationVector& lambda, Float64Approximation& violation) const;
    Void initialise_lagrange_multipliers(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                                         const FloatApproximationVector& x, FloatApproximationVector& lambda) const;
    Float64Approximation compute_mu(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                     const FloatApproximationVector& x, const FloatApproximationVector& lambda) const;

  public: // Deprecated
    Void compute_tz(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                          FloatApproximationVector& x, Float64Approximation& t, FloatApproximationVector& z) const { }
    Void feasibility_step(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                          FloatApproximationVector& x, FloatApproximationVector& y, FloatApproximationVector& z, Float64Approximation& violation) const { };
    Void linearised_feasibility_step(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                                     Float64Approximation& slack, FloatApproximationVector& x, FloatApproximationVector& lambda) const { };
    Void linearised_feasibility_step(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                                     FloatApproximationVector& x, FloatApproximationVector& y, FloatApproximationVector& z, Float64Approximation& t) const { };
  private:
    Float64Approximation compute_mu(const ApproximateScalarFunction& f, const ExactBoxType& D, const ValidatedVectorFunction& g, const ExactBoxType& C,
                     const FloatApproximationVector& x, const FloatApproximationVector& y) const;
    Void compute_violation(const ExactBoxType& D, const ApproximateVectorFunction& g, const ExactBoxType& C,
                           FloatApproximationVector& x, Float64Approximation& t) const;
};




class IntervalOptimiser
    : public NonlinearInteriorPointOptimiser
{
    virtual IntervalOptimiser* clone() const { return new IntervalOptimiser(*this); }
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorFunction h) const;
    Void feasibility_step(const ExactFloatVector& xl, const ExactFloatVector& xu, const ValidatedVectorFunction& h,
                          FloatBoundsVector& x, FloatBoundsVector& y, FloatBoundsVector& zl, FloatBoundsVector zu, Float64Bounds& mu) const;
};


class ApproximateOptimiser
    : public NonlinearInteriorPointOptimiser
{
    virtual ApproximateOptimiser* clone() const { return new ApproximateOptimiser(*this); }
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorFunction h) const;
    Void feasibility_step(const ExactBoxType& D, const ApproximateVectorFunction& h,
                          FloatApproximationVector& X, FloatApproximationVector& Lambda) const;
};


/*//! \ingroup OptimisationModule
//! Solver for linear programming problems using interior point methods.
//! WARNING: This class currently does not work; maybe there is a problem with the algorithms.
class KrawczykOptimiser
    : public OptimiserBase
{

  public:
    virtual KrawczykOptimiser* clone() const { return new KrawczykOptimiser(*this); }

    //! \brief Solve the linear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarFunction f, ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C) const;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C) const;

  public:
    //! \brief Try to solve the nonlinear constraint problem by applying the Krawczyk contractor to the Kuhn-Tucker conditions,
    //! hotstarting the iteration with the primal and dual variables.
    ValidatedKleenean minimise(ValidatedScalarFunction f, ExactBoxType D, ValidatedVectorFunction g, ExactBoxType C,
                     const Float64Bounds& t0, const FloatBoundsVector& x0, const FloatBoundsVector& y0, const FloatBoundsVector& z0) const;

    //! \brief A primal-dual feasibility step for the problem \f$g(y)\in C;\ y\in D\f$.
    Void minimisation_step(const ExactBoxType& D, const ValidatedVectorFunction& g, const ExactBoxType& C,
                           FloatBoundsVector& x, FloatBoundsVector& y, FloatBoundsVector& z, Float64Bounds& t) const;
    //! \brief A primal-dual feasibility step for the problem \f$g(y)\in C;\ y\in D\f$.
    Void feasibility_step(const ExactBoxType& D, const ValidatedVectorFunction& g, const ExactBoxType& C,
                          FloatBoundsVector& x, FloatBoundsVector& y, FloatBoundsVector& z, Float64Bounds& t) const;

    //! \brief A primal feasibility step for the problem \f$g(y)\in C;\ y\in D\f$. \deprecated
    Void feasibility_step(const ExactBoxType& D, const ValidatedVectorFunction& g, const ExactBoxType& C,
                          FloatBoundsVector& y, Float64Bounds& t) const;
    //! \brief A feasibility step for the problem \f$g(y)\leq 0\f$. \deprecated
    Void feasibility_step(const ValidatedVectorFunction& g,
                          FloatBoundsVector& x, FloatBoundsVector& y, FloatBoundsVector& z, Float64Bounds& t) const;
    //! \brief An optimization step for the problem \f$\max f(y) \text{ s.t. } g(y)\leq 0\f$. \deprecated
    Void minimisation_step(const ValidatedScalarFunction& f, const ValidatedVectorFunction& g,
                           FloatBoundsVector& x, FloatBoundsVector& y, FloatBoundsVector& z) const;
  protected:
    Void setup_feasibility(const ExactBoxType& D, const ValidatedVectorFunction& g, const ExactBoxType& C,
                           FloatBoundsVector& x, FloatBoundsVector& y, FloatBoundsVector& z, Float64Bounds& t) const;
    protected:
    Void compute_tz(const ExactBoxType& D, const ValidatedVectorFunction& g, const ExactBoxType& C, const FloatBoundsVector& y, Float64Bounds& t, FloatBoundsVector& z) const;
};

*/


} // namespace Ariadne

#endif
