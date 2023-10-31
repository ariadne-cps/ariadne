/***************************************************************************
 *            solvers/nonlinear_programming.hpp
 *
 *  Copyright  2009-20  Pieter Collins
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

/*! \file solvers/nonlinear_programming.hpp
 *  \brief Nonlinear programming.
 */

#ifndef ARIADNE_NONLINEAR_PROGRAMMING_HPP
#define ARIADNE_NONLINEAR_PROGRAMMING_HPP

#include "utility/declarations.hpp"

#include "conclog/logging.hpp"
#include "numeric/numeric.hpp"
#include "utility/tuple.hpp"

using namespace ConcLog;

namespace Ariadne {

template<class X, class R> class Constraint;

class ProblemException : public std::runtime_error {
  public: ProblemException() : std::runtime_error("ProblemException") { }
};
class IndeterminateFeasibilityException : public std::runtime_error {
  public: IndeterminateFeasibilityException() : std::runtime_error("IndeterminateFeasibilityException") { }
};
class DegenerateFeasibilityProblemException : public std::runtime_error {
  public: DegenerateFeasibilityProblemException() : std::runtime_error("DegenerateFeasibilityProblemException") { }
};
class NearBoundaryOfFeasibleDomainException : public std::runtime_error {
  public: NearBoundaryOfFeasibleDomainException() : std::runtime_error("NearBoundaryOfFeasibleDomainException") { }
};

template<class P> struct FeasibilityProblem {
  private:
    template<class PP> struct BoxTrait;
    template<> struct BoxTrait<ApproximateTag> { typedef ApproximateBoxType Type; };
    template<> struct BoxTrait<ValidatedTag> { typedef ExactBoxType Type; };
  public:
    template<class PP> using BoxType = typename BoxTrait<PP>::Type;
    BoxType<P> D;
    VectorMultivariateFunction<P> g;
    BoxType<P> C;
  public:
    FeasibilityProblem(BoxType<P> D, VectorMultivariateFunction<P> g, BoxType<P> C)
        : D(D), g(g), C(C) { }
    template<class PP> requires std::is_convertible<PP,P>::value FeasibilityProblem(const FeasibilityProblem<PP>& p)
        : FeasibilityProblem(p.D,p.g,p.C) { }
};
template<class P> OutputStream& operator<<(OutputStream& os, FeasibilityProblem<P> const& p);
using ApproximateFeasibilityProblem = FeasibilityProblem<ApproximateTag>;
using ValidatedFeasibilityProblem = FeasibilityProblem<ValidatedTag>;

template<class P> struct OptimisationProblem : public FeasibilityProblem<P> {
  public:
    template<class PP> using BoxType = typename FeasibilityProblem<PP>::template BoxType<PP>;
    ScalarMultivariateFunction<P> f;
  public:
    OptimisationProblem(ScalarMultivariateFunction<P> f, BoxType<P> D, VectorMultivariateFunction<P> g, BoxType<P> C)
        : FeasibilityProblem<P>{D,g,C}, f{f} { }
    template<class PP> requires std::is_convertible<PP,P>::value OptimisationProblem(const OptimisationProblem<PP>& p)
        : OptimisationProblem(p.f,p.D,p.g,p.C) { }
};
template<class P> OutputStream& operator<<(OutputStream& os, OptimisationProblem<P> const& p);
using ApproximateOptimisationProblem = OptimisationProblem<ApproximateTag>;
using ValidatedOptimisationProblem = OptimisationProblem<ValidatedTag>;

//! \brief Data for the solution of an optimization problem, with primal variables of type \a X.
template<class X> struct PrimalData {
    X x;
    //! \brief .
    X const& primal() const { return this->x; }
    //! \brief Convert to the primal variables.
    operator X () const { return this->x; }
};

//! \brief Data for the solution of an optimization problem, with primal variables of type \a X and dual variables of type \a Y.
template<class X, class Y> struct PrimalDualData : public PrimalData<X> {
    Y y;
    //! \brief .
    PrimalDualData(X x, Y y) : PrimalData<X>(x), y(y) { }
    //! \brief .
    Y const& dual() const { return this->y; }
};
using ApproximatePrimalDualData = PrimalDualData<Vector<FloatDPApproximation>,Vector<FloatDPApproximation>>;
using ValidatedPrimalDualData = PrimalDualData<Vector<FloatDPBounds>,Vector<FloatDPBounds>>;

//! \brief Data for the solution of an optimization problem, with primal variables of type \a X, dual variables of type \a Y, and complementary variables of type \a Z.
template<class X, class Y, class Z=X> struct PrimalDualComplementaryData : public PrimalDualData<X,Y> {
    Z z;
    //! \brief .
    PrimalDualComplementaryData(X x, Y y, Z z) : PrimalDualData<X,Y>(x,y), z(z) { }
    //! \brief .
    Z const& complementary() const { return this->z; }
};
using ApproximatePrimalDualComplementaryData = PrimalDualComplementaryData<Vector<FloatDPApproximation>,Vector<FloatDPApproximation>>;
using ValidatedPrimalDualComplementaryData = PrimalDualComplementaryData<Vector<FloatDPBounds>,Vector<FloatDPBounds>>;

//! \brief Data for the solution of an optimization problem, with primal variables of type \a X and slack variables of type \a W.
template<class X, class W=X> struct SlackPrimalData : public PrimalData<X> {
    W w;
    //! \brief .
    SlackPrimalData(W w, X x) : PrimalData<X>(x), w(w) { }
    //! \brief .
    W const& slack() const { return this->w; }
};
using ApproximateSlackPrimalData = SlackPrimalData<Vector<FloatDPApproximation>,Vector<FloatDPApproximation>>;
using ValidatedSlackPrimalData = SlackPrimalData<Vector<FloatDPBounds>,Vector<FloatDPBounds>>;

//! \brief Data for the solution of an optimization problem, with primal variables of type \a X dual variables of type \a Y, slack variables of type \a W and complementary variables of type \a Z.
template<class W, class X=W, class Y=W, class Z=X> struct SlackPrimalDualComplementaryData : public PrimalDualComplementaryData<X,Y,Z> {
    //! \brief The slack variables.
    W w;
    //! \brief .
    SlackPrimalDualComplementaryData(W w, X x, Y y, Z z) : PrimalDualComplementaryData<X,Y,Z>(x,y,z), w(w) { }
    //! \brief .
    W const& slack() const { return this->w; }
};
using ApproximateSlackPrimalDualComplementaryData = SlackPrimalDualComplementaryData<Vector<FloatDPApproximation>,Vector<FloatDPApproximation>>;
using ValidatedSlackPrimalDualComplementaryData = SlackPrimalDualComplementaryData<Vector<FloatDPBounds>,Vector<FloatDPBounds>>;


//! \ingroup OptimisationSubModule EvaluationModule
//! Interface for nonlinear programming solvers.
class OptimiserInterface {
    using FLT=FloatDP;
  public:
    typedef Bounds<FLT> ValidatedNumericType;
    typedef Approximation<FLT> ApproximateNumericType;
    typedef Vector<FLT> ExactVectorType;
    typedef Vector<Bounds<FLT>> ValidatedVectorType;
    typedef Vector<Approximation<FLT>> ApproximateVectorType;
  public:
    //! \brief Virtual destructor.
    virtual ~OptimiserInterface() = default;
    //! \brief Create a dynamically-allocated copy.
    virtual OptimiserInterface* clone() const = 0;

    //! \brief Solve the general nonlinear programming problem \f$\min f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    virtual Vector<ValidatedNumericType> minimise(ValidatedOptimisationProblem p) const = 0;
    //! \brief Solve the general nonlinear programming problem \f$\min f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    //! \pre The domain \f$D\f$ is bounded and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A point \f$x\f$ is approximately satisfies the feasibility constraints (it may not be which definitely contains a feasible point, and contains a local optimum.
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const = 0;
    //! \brief Solve the standard nonlinear programming problem \f$\min f(x) \text{ such that } x\in D\ g(x)\leq 0  \text{ and } h(x) = 0\f$.
    //! \pre The domain \f$D\f$ is bounded and has nonempty interior.
    //! \return A point \f$x\f$ is approximately satisfies the feasibility constraints (it may not be which definitely contains a feasible point, and contains a local optimum.
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const = 0;

    //! \brief Solve the general nonlinear programming problem \f$\min f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    virtual Vector<ApproximateNumericType> minimise(ApproximateOptimisationProblem p) const = 0;
    //! \brief Approximatedly solve the general nonlinear programming problem \f$\min f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    //! \return A point \f$x\f$ is approximately satisfies the feasibility constraints (it may not be close to a genuinely feasible point) and is approximately locally optimal.
    virtual Vector<ApproximateNumericType> minimise(ApproximateScalarMultivariateFunction f, ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateBoxType C) const = 0;
    //! \brief Approximately solve the standard nonlinear programming problem \f$\min f(x) \text{ such that } x\in D\ g(x)\leq 0  \text{ and } h(x) = 0\f$.
    virtual Vector<ApproximateNumericType> minimise(ApproximateScalarMultivariateFunction f, ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateVectorMultivariateFunction h) const = 0;

    //! \brief Tests is the general nonlinear feasibility problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ValidatedFeasibilityProblem p) const = 0;
    //! \brief Tests is the general nonlinear feasibility problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const = 0;
    //! \brief Tests is the standard nonlinear feasibility problem \f$x\in D,\ g(x)\leq 0 \text{ and } h(x) = 0\f$ is feasible. Assumes \f$D\f$ is bounded with nonempty interior.
    //! \internal This is one of the simplest nonlinear programming problems, and is a good test case for new algorithms.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const = 0;

    //! \brief Tests if the point \a x is almost feasible, in that \f$x\in D\f$ and \f$g(x)\in N_\epsilon(C)\f$.
    virtual ApproximateKleenean almost_feasible_point(ValidatedFeasibilityProblem p,
                                                      ApproximateVectorType x, FloatDPApproximation eps) const = 0;
    //! \brief Tests whether the point \a x is feasible.
    virtual ValidatedKleenean is_feasible_point(ValidatedFeasibilityProblem p,
                                                ExactVectorType x) const = 0;
    //! \brief Tests whether the box \a X contains a feasible point.
    virtual ValidatedKleenean contains_feasible_point(ValidatedFeasibilityProblem p,
                                                      ValidatedVectorType x) const = 0;
    //! \brief Tests whether the point \a x is a feasible, or if the Lagrange multipliers \a y are a certificate of infeasibility.
    virtual ValidatedKleenean check_feasibility(ValidatedFeasibilityProblem p,
                                                ExactVectorType x, ExactVectorType y) const = 0;

    //! \brief Checks if the point \a x is definitely feasible.
    virtual Bool validate_feasibility(ValidatedFeasibilityProblem p,
                                      ExactVectorType x) const = 0;
    //! \brief Tests if the feasibility problem is definitely unsolvable, using multipliers \a y and local centering point \a x.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        ExactVectorType x, ExactVectorType y) const = 0;
    //! \brief Tests if the feasibility problem is definitely unsolvable over the box \a X, using multipliers \a y.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        UpperBoxType X, ExactVectorType y) const = 0;
    //! \brief Tests if the Lagrange multipliers \a y are a certificate of infeasiblity.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        ExactVectorType y) const = 0;
};

//! \ingroup OptimisationSubModule
//! Common routines for nonlinear programming
class OptimiserBase
    : public OptimiserInterface
{
  protected:
    static const FloatDP zero;
    static const FloatDP one;
  public:
    virtual Vector<ValidatedNumericType> minimise(ValidatedOptimisationProblem p) const override = 0;
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const override;
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const override;

    virtual Vector<ApproximateNumericType> minimise(ApproximateOptimisationProblem p) const override = 0;
    virtual Vector<ApproximateNumericType> minimise(ApproximateScalarMultivariateFunction f, ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateBoxType C) const override;
    virtual Vector<ApproximateNumericType> minimise(ApproximateScalarMultivariateFunction f, ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateVectorMultivariateFunction h) const override;

    virtual ValidatedKleenean feasible(ValidatedFeasibilityProblem p) const override = 0;
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const override;
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const override;

    virtual ApproximateKleenean almost_feasible_point(ValidatedFeasibilityProblem p,
                                                      ApproximateVectorType x, FloatDPApproximation error) const override;
    virtual ValidatedKleenean is_feasible_point(ValidatedFeasibilityProblem p,
                                                ExactVectorType x) const override;
    virtual ValidatedKleenean contains_feasible_point(ValidatedFeasibilityProblem p,
                                                      ValidatedVectorType X) const override;
    virtual ValidatedKleenean check_feasibility(ValidatedFeasibilityProblem p,
                                                ExactVectorType x, ExactVectorType y) const override;

    virtual Bool validate_feasibility(ValidatedFeasibilityProblem p,
                                      ExactVectorType x) const override;
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        ExactVectorType x, ExactVectorType y) const override;
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        UpperBoxType X, ExactVectorType y) const override;
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        ExactVectorType y) const override;
};

//! \ingroup OptimisationSubModule
//! \brief Solver for nonlinear programming problems based on a penalty-function approach
//!
//! For the feasibility problem \f$x\in D,\ g(x)\in C,\ h(x)=0\f$ where \f$D\f$ is bounded and \f$D,C\f$ have nonempty interiors.
//! Introduce slack variable \f$w=g(x)\f$, and minimise \f[ \sum_{j} (g_j(x)-w_j)^2 + \sum_k h_k(x)^2 \f] with \f$x\in D\f$ and \f$w\in C\f$.
//! Since the minimiser is not unique, we need to add penalty terms \f$-\mu/2(\log(x_u-x)+log(x-x_l))\f$ and \f$-\nu/2(\log(w_u-w)+log(w-w_l))\f$.
//! It suffices to add a penalty in \f$x\f$.
class PenaltyFunctionOptimiser
    : public OptimiserBase
{
  public:
    using OptimiserBase::minimise;
    using OptimiserBase::feasible;
  public:
    virtual PenaltyFunctionOptimiser* clone() const override;

    virtual Vector<ValidatedNumericType> minimise(ValidatedOptimisationProblem p) const override;
    virtual Vector<ApproximateNumericType> minimise(ApproximateOptimisationProblem p) const override;
    virtual ValidatedKleenean feasible(ValidatedFeasibilityProblem p) const override;

    virtual ValidatedKleenean check_feasibility(ValidatedFeasibilityProblem p, ExactVectorType x, ExactVectorType y) const override;
  public:
    virtual Void feasibility_step(ValidatedFeasibilityProblem p,
                                  ValidatedVectorType& w, ValidatedVectorType& x) const;
    virtual Void feasibility_step(ApproximateFeasibilityProblem p,
                                  ApproximateVectorType& w, ApproximateVectorType& x, FloatDPApproximation& mu) const;
    virtual Void feasibility_step(ApproximateFeasibilityProblem p,
                                  ApproximateVectorType& w, ApproximateVectorType& x, ApproximateVectorType& y) const;
};




//! \ingroup OptimisationSubModule
//! Solver for nonlinear programming problems using the infeasible interior point method with primal, dual and slack variables.
//! \details Relies on an engine minimising \f$f(x)\f$ for \f$x\in D^\circ\f$ subject to \f$g(x)=w\f$ with \f$w\in C^\circ\f$ and \f$h(x)=0\f$.
class InteriorPointOptimiser
    : public OptimiserBase
{
  public:
    virtual InteriorPointOptimiser* clone() const override;

    using OptimiserBase::minimise;
    using OptimiserBase::feasible;

  public:
    //! \brief Construct with default accuracy parameters.
    InteriorPointOptimiser();

    //! \brief Compute a \em local optimum of nonlinear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C .\f$.
    //! \pre The domain \f$D\f$ is bounded and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A box \f$X\f$ which definitely contains a feasible point, and contains a local optimum.
    virtual Vector<ValidatedNumericType> minimise(ValidatedOptimisationProblem p) const override;
    //! \brief Compute an approximate \em local optimum of nonlinear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C.\f$.
    virtual Vector<ApproximateNumericType> minimise(ApproximateOptimisationProblem p) const override;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ValidatedFeasibilityProblem) const override;

    //! \brief Test if the constraints \f$g(y)\in C\f$ are solvable for \f$y\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the overall constraint violation, primal and dual variables.
    Pair<ValidatedKleenean,FloatDPApproximationVector> feasible_hotstarted(ValidatedFeasibilityProblem p,
                                                  const FloatDPApproximationVector& x0, const FloatDPApproximationVector& lambda0, const FloatDPApproximation& t0) const;

    //! \brief Test if the constraints \f$g(y)\in C\f$ are solvable for \f$y\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal and dual.
    Pair<ValidatedKleenean,FloatDPApproximationVector> feasible_hotstarted(ValidatedFeasibilityProblem p,
                                                  const FloatDPApproximationVector& x0, const FloatDPApproximationVector& lambda0) const;

    Void minimisation_step(const ApproximateOptimisationProblem& p,
                           const ApproximateVectorMultivariateFunction& h,
                           FloatDPApproximationVector& w, FloatDPApproximationVector& x,
                           FloatDPApproximationVector& kappa, FloatDPApproximationVector& lambda, const FloatDPApproximation& mu) const;
    Void setup_feasibility(const ApproximateFeasibilityProblem& p,
                           FloatDPApproximationVector& x, FloatDPApproximationVector& y) const;
    Void setup_feasibility(const ApproximateFeasibilityProblem& p,
                           FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximation& t) const;
    Void feasibility_step(const ApproximateFeasibilityProblem& p,
                          FloatDPApproximationVector& x, FloatDPApproximationVector& y) const;
    Void feasibility_step(const ApproximateFeasibilityProblem& p,
                          FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximation& t) const;
    Void initialise_lagrange_multipliers(const ApproximateFeasibilityProblem& p,
                                         const FloatDPApproximationVector& x, FloatDPApproximationVector& y) const;
    FloatDPApproximation compute_mu(const ApproximateFeasibilityProblem& p,
                                    const FloatDPApproximationVector& x, const FloatDPApproximationVector& y) const;

    friend OutputStream& operator<<(OutputStream& os, InteriorPointOptimiser const& opt);
  public: // Deprecated
    Void compute_tz(const ApproximateBoxType& D, const ApproximateVectorMultivariateFunction& g, const ApproximateBoxType& C,
                    FloatDPApproximationVector& x, FloatDPApproximation& t, FloatDPApproximationVector& z) const;
    Void feasibility_step(const ApproximateBoxType& D, const ApproximateVectorMultivariateFunction& g, const ApproximateBoxType& C,
                          FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximationVector& z, FloatDPApproximation& t) const;
    Void linearised_feasibility_step(const ApproximateBoxType& D, const ApproximateVectorMultivariateFunction& g, const ApproximateBoxType& C,
                                     FloatDPApproximation& w, FloatDPApproximationVector& x, FloatDPApproximationVector& lambda) const;
    Void linearised_feasibility_step(const ApproximateBoxType& D, const ApproximateVectorMultivariateFunction& g, const ApproximateBoxType& C,
                                     FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximationVector& z, FloatDPApproximation& t) const;
  private:
    FloatDPApproximation compute_mu(const ApproximateScalarMultivariateFunction& f, const ApproximateBoxType& D, const ValidatedVectorMultivariateFunction& g, const ApproximateBoxType& C,
                     const FloatDPApproximationVector& x, const FloatDPApproximationVector& y) const;
    Void compute_violation(const ApproximateFeasibilityProblem& p,
                           FloatDPApproximationVector& x, FloatDPApproximation& t) const;
};


//! \ingroup OptimisationSubModule
//! \brief Solver for nonlinear programming problems using infeasible interior point methods.
//! \details Introduces variables \f$w\f$ and attempts to find \f$x\in D\f$ and \f$w\in C\f$ such that \f$g(x)=w\f$.
//!   The dual variables \f$y\f$ are unconstrained Lagrange multipliers for \f$y\cdot(g(x)-w)=0\f$.
class InfeasibleInteriorPointOptimiser
    : public OptimiserBase
{
  public:
    virtual InfeasibleInteriorPointOptimiser* clone() const override;

    using OptimiserBase::minimise;
    using OptimiserBase::feasible;

    struct PrimalDualData;
    struct StepData;

  public:
    //! \brief Construct with default accuracy parameters.
    InfeasibleInteriorPointOptimiser();

    //! \brief Compute a \em local optimum of nonlinear programming problem \f$\max f(x) \text{ such that } x\in D, g(x)\in C \text{ and } h(x)=0.\f$.
    //! \pre The domain \f$D\f$ is bounded and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A box \f$X\f$ which definitely contains a feasible point, and contains a local optimum.
    virtual Vector<ValidatedNumericType> minimise(ValidatedOptimisationProblem p) const override;
    //! \brief Compute an approximate \em local optimum of nonlinear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C.\f$.
    virtual Vector<ApproximateNumericType> minimise(ApproximateOptimisationProblem p) const override;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ValidatedFeasibilityProblem) const override;

    friend OutputStream& operator<<(OutputStream& os, InfeasibleInteriorPointOptimiser const& opt);
  public:
    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the overall primal and dual variables.
    Pair<ValidatedKleenean,FloatDPApproximationVector> feasible_hotstarted(ValidatedFeasibilityProblem p,
                                                                const PrimalDualData& wxy0) const;

    Void setup_feasibility(const ApproximateFeasibilityProblem& p,
                           StepData& stp) const;
    Void step(const ApproximateOptimisationProblem& p,
              StepData& stp) const;




//    FloatDPApproximation compute_mu(const ApproximateFeasibilityProblem& p,
//                     FloatDPApproximationVector& w, FloatDPApproximationVector& x, FloatDPApproximationVector& y) const;
};





class IntervalOptimiser
    : public InteriorPointOptimiser
{
    virtual IntervalOptimiser* clone() const override;
    virtual ValidatedKleenean feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const;
    Void feasibility_step(const FloatDPVector& xl, const FloatDPVector& xu, const ValidatedVectorMultivariateFunction& h,
                          FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& zl, FloatDPBoundsVector zu, FloatDPBounds& mu) const;
};


class ApproximateOptimiser
    : public InteriorPointOptimiser
{
    virtual ApproximateOptimiser* clone() const override;
    virtual ValidatedKleenean feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const;
    Void feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& h,
                          FloatDPApproximationVector& X, FloatDPApproximationVector& Lambda) const;
};


/*//! \ingroup OptimisationSubModule
//! Solver for nonlinear programming problems using interior point methods.
//! WARNING: This class currently does not work; maybe there is a problem with the algorithms.
class KrawczykOptimiser
    : public OptimiserBase
{

  public:
    virtual KrawczykOptimiser* clone() const;

    //! \brief Solve the linear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    virtual Vector<ValidatedNumericType> minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const;

  public:
    //! \brief Try to solve the nonlinear constraint problem by applying the Krawczyk contractor to the Kuhn-Tucker conditions,
    //! hotstarting the iteration with the primal and dual variables.
    ValidatedKleenean minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                     const FloatDPBounds& t0, const FloatDPBoundsVector& x0, const FloatDPBoundsVector& y0, const FloatDPBoundsVector& z0) const;

    //! \brief A primal-dual feasibility step for the problem \f$g(y)\in C;\ y\in D\f$.
    Void minimisation_step(const ExactBoxType& D, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& C,
                           FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& z, FloatDPBounds& t) const;
    //! \brief A primal-dual feasibility step for the problem \f$g(y)\in C;\ y\in D\f$.
    Void feasibility_step(const ExactBoxType& D, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& C,
                          FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& z, FloatDPBounds& t) const;

    //! \brief A primal feasibility step for the problem \f$g(y)\in C;\ y\in D\f$. \deprecated
    Void feasibility_step(const ExactBoxType& D, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& C,
                          FloatDPBoundsVector& y, FloatDPBounds& t) const;
    //! \brief A feasibility step for the problem \f$g(y)\leq 0\f$. \deprecated
    Void feasibility_step(const ValidatedVectorMultivariateFunction& g,
                          FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& z, FloatDPBounds& t) const;
    //! \brief An optimization step for the problem \f$\max f(y) \text{ s.t. } g(y)\leq 0\f$. \deprecated
    Void minimisation_step(const ValidatedScalarMultivariateFunction& f, const ValidatedVectorMultivariateFunction& g,
                           FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& z) const;
  protected:
    Void setup_feasibility(const ExactBoxType& D, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& C,
                           FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& z, FloatDPBounds& t) const;
    protected:
    Void compute_tz(const ExactBoxType& D, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& C, const FloatDPBoundsVector& y, FloatDPBounds& t, FloatDPBoundsVector& z) const;
};

*/


} // namespace Ariadne

#endif
