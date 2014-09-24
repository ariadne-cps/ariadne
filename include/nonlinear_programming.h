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

#include "declarations.h"

#include "logging.h"
#include "numeric.h"
#include "tuple.h"


namespace Ariadne {

template<class X, class R> class Constraint;
typedef Constraint<EffectiveScalarFunction,EffectiveNumberType> EffectiveConstraint;
typedef Constraint<ValidatedScalarFunction,ExactNumberType> ValidatedConstraint;

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
    virtual Vector<ValidatedNumberType> minimise(ValidatedScalarFunction f, Box D, ValidatedVectorFunction g, Box C) const = 0;
    //! \brief Solve the standard nonlinear programming problem \f$\min f(x) \text{ such that } x\in ,D\ g(x)\leq 0  \text{ and } h(x) = 0\f$.
    virtual Vector<ValidatedNumberType> minimise(ValidatedScalarFunction f, Box D, ValidatedVectorFunction g, ValidatedVectorFunction h) const = 0;

    //! \brief Tests is the general nonlinear feasibility problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual Tribool feasible(Box D, ValidatedVectorFunction g, Box C) const = 0;
    //! \brief Tests is the standard nonlinear feasibility problem \f$x\in D,\ g(x)\leq 0 \text{ and } h(x) = 0\f$ is feasible. Assumes \fD\f$ is bounded with nonempty interior.
    //! \internal This is one of the simplest nonlinear programming problems, and is a good test case for new algorithms.
    virtual Tribool feasible(Box D, ValidatedVectorFunction g, ValidatedVectorFunction h) const = 0;

    //! \brief Tests if the point \a x is feasible, in that \f$x\in D\f$ and \f$g(x)\in N_\epsilon(C)\f$.
    virtual Bool almost_feasible_point(Box D, ValidatedVectorFunction g, Box C,
                                       ApproximateFloatVectorType x, ApproximateFloatType eps) const = 0;
    //! \brief Tests if the point \a x is feasible.
    virtual Bool is_feasible_point(Box D, ValidatedVectorFunction g, Box C,
                                   ExactFloatVectorType x) const = 0;
    //! \brief Tests if the point \a x is near feasible.
    virtual Bool validate_feasibility(Box D, ValidatedVectorFunction g, Box C,
                                      ExactPointType x) const = 0;
    //! \brief Tests if the point \a x is near feasible, using approximate multipliers \a y to guide the search.
    virtual Bool validate_feasibility(Box D, ValidatedVectorFunction g, Box C,
                                      ExactPointType x, ExactPointType y) const = 0;
    //! \brief Tests if the feasibility problem is definitely unsolvable, using multipliers \a y and local centering point \a x.
    virtual Bool validate_infeasibility(Box D, ValidatedVectorFunction g, Box C,
                                      ExactPointType x, ExactPointType y) const = 0;
    //! \brief Tests if the box \a X definitely containss a feasible point.
    virtual Tribool contains_feasible_point(Box D, ValidatedVectorFunction g, Box C,
                                            ValidatedFloatVectorType X) const = 0;
    //! \brief Tests if the Lagrange multipliers \a y are a certificate of infeasiblity.
    virtual Bool is_infeasibility_certificate(Box D, ValidatedVectorFunction g, Box C,
                                              ExactPointType y) const = 0;
};

//! \ingroup OptimisationModule
//! Common routines for nonlinear minimisation
class OptimiserBase
    : public OptimiserInterface
    , public Loggable
{
  public:
    virtual ValidatedPointType minimise(ValidatedScalarFunction f, Box D, ValidatedVectorFunction g, Box C) const = 0;
    virtual ValidatedPointType minimise(ValidatedScalarFunction f, Box D, ValidatedVectorFunction g, ValidatedVectorFunction h) const;

    virtual Tribool feasible(Box D, ValidatedVectorFunction g, Box C) const = 0;
    virtual Tribool feasible(Box D, ValidatedVectorFunction g, ValidatedVectorFunction h) const;

    virtual Bool almost_feasible_point(Box D, ValidatedVectorFunction g, Box C,
                                       ApproximatePointType x, ApproximateFloatType error) const;
    virtual Bool is_feasible_point(Box D, ValidatedVectorFunction g, Box C,
                                   ExactPointType x) const;
    virtual Bool validate_feasibility(Box D, ValidatedVectorFunction g, Box C,
                                      ExactPointType x) const;
    virtual Bool validate_feasibility(Box D, ValidatedVectorFunction g, Box C,
                                      ExactPointType x, ExactPointType y) const;
    virtual Bool validate_infeasibility(Box D, ValidatedVectorFunction g, Box C,
                                        ExactPointType x, ExactVectorType y) const;
    virtual Tribool contains_feasible_point(Box D, ValidatedVectorFunction g, Box C,
                                            ValidatedPointType X) const;
    virtual Bool is_infeasibility_certificate(Box D, ValidatedVectorFunction g, Box C,
                                              ExactVectorType lambda) const;
};

//! \ingroup OptimisationModule
//! \brief Solver for feasibility problems based on a penalty-function approach
//!
//! For the feasibility problem \f$x\in D,\ g(x)\in C,\ h(x)=0\f$ where \f$D\f$ is bounded and \f$D,C\f$ have nonempty interiors.
//! Introduce slack variable \f$w=g(x)\f$, and minimise \f[ \sum_{j} (g_j(x)-w_j)^2 + \sum_k h_k(x)^2 \f] with \f$x\in D\f$ and \f$w\in C\f$.
//! Since the minimiser is not unique, we need to add penalty terms \f$-\mu/2(\log(x_u-x)+log(x-x_l))\f$ and \f$-\nu/2(\log(w_u-w)+log(w-w_l))\f$.
//! It suffices to add a penalty in \f$x\f$.
class PenaltyFunctionOptimiser
    : public OptimiserBase
{
  public:
    virtual PenaltyFunctionOptimiser* clone() const;
    virtual Tribool check_feasibility(Box D, ValidatedVectorFunction g, Box C, ExactPointType x, ExactVectorType y) const;
    virtual ValidatedPointType minimise(ValidatedScalarFunction f, Box D, ValidatedVectorFunction g, Box C) const;
    virtual Tribool feasible(Box D, ValidatedVectorFunction g, Box C) const;
    virtual Void feasibility_step(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                                  ApproximateFloatVectorType& x, ApproximateFloatVectorType& w, ApproximateFloatType& mu) const;
    virtual Void feasibility_step(const Box& D, const ValidatedVectorFunction& g, const Box& C,
                                  ValidatedFloatVectorType& x, ValidatedFloatVectorType& w) const;
    virtual Void feasibility_step(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                                  ApproximateFloatVectorType& x, ApproximateFloatVectorType& y, ApproximateFloatVectorType& z) const;
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
    //! \precondition The domain \f$D\f$ is bounded and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A box \f$X\f$ which definitely contains a feasible point, and contains a local optimum.
    virtual ValidatedPointType minimise(ValidatedScalarFunction f, Box D, ValidatedVectorFunction g, Box C) const;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual Tribool feasible(Box D, ValidatedVectorFunction g, Box C) const;

    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the overall primal and dual variables.
    Pair<Tribool,ApproximateFloatVectorType> feasible_hotstarted(Box D, ValidatedVectorFunction g, Box C,
                                                  const PrimalDualData& wxy0) const;

    Void setup_feasibility(const Box& D, const ApproximateVectorFunctionInterface& g, const Box& C,
                           StepData& stp) const;
    Void step(const ApproximateScalarFunctionInterface& f, const Box& D, const ApproximateVectorFunctionInterface& g, const Box& C,
              StepData& stp) const;
//    ApproximateFloatType compute_mu(const Box& D, const ApproximateVectorFunction& g, const Box& C,
//                     ApproximateFloatVectorType& w, ApproximateFloatVectorType& x, ApproximateFloatVectorType& y) const;
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
    //! \precondition The domain \f$D\f$ is bounded and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A box \f$X\f$ which definitely contains a feasible point, and contains a local optimum.
    virtual ValidatedPointType minimise(ValidatedScalarFunction f, Box D, ValidatedVectorFunction g, Box C) const;
    //! \brief Tests is the nonlinear programming problem \f$x\in D, g(x)\in C \text{ and } h(x)= 0 \f$ is feasible.
    virtual Tribool feasible(Box D, ValidatedVectorFunction g, Box C) const;

    //! \brief Test if the constraints \f$g(y)\in C\f$ are solvable for \f$y\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the overall constraint violation, primal and dual variables.
    Pair<Tribool,ApproximateFloatVectorType> feasible_hotstarted(Box D, ValidatedVectorFunction g, Box C,
                                                  const ApproximateFloatVectorType& x0, const ApproximateFloatVectorType& lambda0, const ApproximateFloatType& violation0) const;

    //! \brief Test if the constraints \f$g(y)\in C\f$ are solvable for \f$y\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal and dual.
    Pair<Tribool,ApproximateFloatVectorType> feasible_hotstarted(Box D, ValidatedVectorFunction g, Box C,
                                                  const ApproximateFloatVectorType& x0, const ApproximateFloatVectorType& lambda0) const;

    Void minimisation_step(const ApproximateScalarFunction& f, const Box& D, const ApproximateVectorFunction& g, const Box& C, const ApproximateVectorFunction& h,
                           ApproximateFloatVectorType& x, ApproximateFloatVectorType& w, ApproximateFloatVectorType& kappa, ApproximateFloatVectorType& lambda, const ApproximateFloatType& mu) const;
    Void setup_feasibility(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                           ApproximateFloatVectorType& x, ApproximateFloatVectorType& lambda) const;
    Void setup_feasibility(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                           ApproximateFloatVectorType& x, ApproximateFloatVectorType& lambda, ApproximateFloatType& t) const;
    Void feasibility_step(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                          ApproximateFloatVectorType& x, ApproximateFloatVectorType& lambda) const;
    Void feasibility_step(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                          ApproximateFloatVectorType& x, ApproximateFloatVectorType& lambda, ApproximateFloatType& violation) const;
    Void initialise_lagrange_multipliers(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                                         const ApproximateFloatVectorType& x, ApproximateFloatVectorType& lambda) const;
    ApproximateFloatType compute_mu(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                     const ApproximateFloatVectorType& x, const ApproximateFloatVectorType& lambda) const;

  public: // Deprecated
    Void compute_tz(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                          ApproximateFloatVectorType& x, ApproximateFloatType& t, ApproximateFloatVectorType& z) const { }
    Void feasibility_step(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                          ApproximateFloatVectorType& x, ApproximateFloatVectorType& y, ApproximateFloatVectorType& z, ApproximateFloatType& violation) const { };
    Void linearised_feasibility_step(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                                     ApproximateFloatType& slack, ApproximateFloatVectorType& x, ApproximateFloatVectorType& lambda) const { };
    Void linearised_feasibility_step(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                                     ApproximateFloatVectorType& x, ApproximateFloatVectorType& y, ApproximateFloatVectorType& z, ApproximateFloatType& t) const { };
  private:
    ApproximateFloatType compute_mu(const ApproximateScalarFunction& f, const Box& D, const ValidatedVectorFunction& g, const Box& C,
                     const ApproximateFloatVectorType& x, const ApproximateFloatVectorType& y) const;
    Void compute_violation(const Box& D, const ApproximateVectorFunction& g, const Box& C,
                           ApproximateFloatVectorType& x, ApproximateFloatType& t) const;
};




class IntervalOptimiser
    : public NonlinearInteriorPointOptimiser
{
    virtual IntervalOptimiser* clone() const { return new IntervalOptimiser(*this); }
    virtual Tribool feasible(Box D, ValidatedVectorFunction h) const;
    Void feasibility_step(const ApproximateFloatVectorType& xl, const ApproximateFloatVectorType& xu, const ValidatedVectorFunction& h,
                          ValidatedFloatVectorType& x, ValidatedFloatVectorType& y, ValidatedFloatVectorType& zl, ValidatedFloatVectorType zu, ValidatedFloatType& mu) const;
};


class ApproximateOptimiser
    : public NonlinearInteriorPointOptimiser
{
    virtual ApproximateOptimiser* clone() const { return new ApproximateOptimiser(*this); }
    virtual Tribool feasible(Box D, ValidatedVectorFunction h) const;
    Void feasibility_step(const Box& D, const ApproximateVectorFunction& h,
                          ApproximateFloatVectorType& X, ApproximateFloatVectorType& Lambda) const;
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
    virtual Vector<ValidatedNumberType> minimise(ValidatedScalarFunction f, Box D, ValidatedVectorFunction g, Box C) const;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual tribool feasible(Box D, ValidatedVectorFunction g, Box C) const;

  public:
    //! \brief Try to solve the nonlinear constraint problem by applying the Krawczyk contractor to the Kuhn-Tucker conditions,
    //! hotstarting the iteration with the primal and dual variables.
    tribool minimise(ValidatedScalarFunction f, Box D, ValidatedVectorFunction g, Box C,
                     const ValidatedFloatType& t0, const ValidatedFloatVectorType& x0, const ValidatedFloatVectorType& y0, const ValidatedFloatVectorType& z0) const;

    //! \brief A primal-dual feasibility step for the problem \f$g(y)\in C;\ y\in D\f$.
    void minimisation_step(const Box& D, const ValidatedVectorFunction& g, const Box& C,
                           ValidatedFloatVectorType& x, ValidatedFloatVectorType& y, ValidatedFloatVectorType& z, ValidatedFloatType& t) const;
    //! \brief A primal-dual feasibility step for the problem \f$g(y)\in C;\ y\in D\f$.
    void feasibility_step(const Box& D, const ValidatedVectorFunction& g, const Box& C,
                          ValidatedFloatVectorType& x, ValidatedFloatVectorType& y, ValidatedFloatVectorType& z, ValidatedFloatType& t) const;

    //! \brief A primal feasibility step for the problem \f$g(y)\in C;\ y\in D\f$. \deprecated
    void feasibility_step(const Box& D, const ValidatedVectorFunction& g, const Box& C,
                          ValidatedFloatVectorType& y, ValidatedFloatType& t) const;
    //! \brief A feasibility step for the problem \f$g(y)\leq 0\f$. \deprecated
    void feasibility_step(const ValidatedVectorFunction& g,
                          ValidatedFloatVectorType& x, ValidatedFloatVectorType& y, ValidatedFloatVectorType& z, ValidatedFloatType& t) const;
    //! \brief An optimization step for the problem \f$\max f(y) \text{ s.t. } g(y)\leq 0\f$. \deprecated
    void minimisation_step(const ValidatedScalarFunction& f, const ValidatedVectorFunction& g,
                           ValidatedFloatVectorType& x, ValidatedFloatVectorType& y, ValidatedFloatVectorType& z) const;
  protected:
    void setup_feasibility(const Box& D, const ValidatedVectorFunction& g, const Box& C,
                           ValidatedFloatVectorType& x, ValidatedFloatVectorType& y, ValidatedFloatVectorType& z, ValidatedFloatType& t) const;
    protected:
    void compute_tz(const Box& D, const ValidatedVectorFunction& g, const Box& C, const ValidatedFloatVectorType& y, ValidatedFloatType& t, ValidatedFloatVectorType& z) const;
};

*/


} // namespace Ariadne

#endif
