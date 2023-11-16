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

class InfeasibleProblemException : public std::runtime_error {
  public: InfeasibleProblemException() : std::runtime_error("InfeasibleProblemException") { }
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

template<class P> struct BoxTrait;
template<> struct BoxTrait<ApproximateTag> { typedef ApproximateBoxType Type; };
template<> struct BoxTrait<ValidatedTag> { typedef ExactBoxType Type; };

//! \ingroup OptimisationSubModule
//! \brief The type used as a box domain or codomain for a nonlinear programming problem with information \a P.
template<class P> using BoxType = typename BoxTrait<P>::Type;

//! \ingroup OptimisationSubModule
//! \brief Data for the feasibility problem \f$x\in D \wedge g(x) \in C\f$, where \f$D\subset\R^n\f$ and \f$C\subset \R^m\f$ are boxes, and \f$g:\R^n\to\R^m\f$ is continuous.
template<class P> struct FeasibilityProblem {
  public:
    BoxType<P> D; //!< <p/>
    VectorMultivariateFunction<P> g; //!< <p/>
    BoxType<P> C; //!< <p/>
  public:
    //! <p/>
    FeasibilityProblem(BoxType<P> D_, VectorMultivariateFunction<P> g_, BoxType<P> C_)
        : D(D_), g(g_), C(C_) { }
    //! <p/>
    template<class PP> requires Convertible<PP,P> FeasibilityProblem(const FeasibilityProblem<PP>& p)
        : FeasibilityProblem(p.D,p.g,p.C) { }
};
template<class P> OutputStream& operator<<(OutputStream& os, FeasibilityProblem<P> const& p);
//! \relates FeasibilityProblem
//! \ingroup OptimisationSubModule
//! \name Type synonyms
//!@{
using ApproximateFeasibilityProblem = FeasibilityProblem<ApproximateTag>; //!< <p/>
using ValidatedFeasibilityProblem = FeasibilityProblem<ValidatedTag>; //!< <p/>
//!@}

//! \ingroup OptimisationSubModule
//! \brief Data for the feasibility problem \f$\text{minimise } f(x) \text{ subject to } x\in D \text{ and } g(x) \in C\f$, where \f$D\subset\R^n\f$ and \f$C\subset \R^m\f$ are boxes, and \f$f:\R^n\to\R\f$ and \f$g:\R^n\to\R^m\f$ are continuous.
template<class P> struct OptimisationProblem : public FeasibilityProblem<P> {
  public:
    ScalarMultivariateFunction<P> f;//!< <p/>
  public:
    //! <p/>
    OptimisationProblem(ScalarMultivariateFunction<P> f_, BoxType<P> D_, VectorMultivariateFunction<P> g_, BoxType<P> C_)
        : FeasibilityProblem<P>(D_,g_,C_), f(f_) { }
    //! <p/>
    template<class PP> requires Convertible<PP,P> OptimisationProblem(const OptimisationProblem<PP>& p)
        : OptimisationProblem(p.f,p.D,p.g,p.C) { }
};
template<class P> OutputStream& operator<<(OutputStream& os, OptimisationProblem<P> const& p);
//! \relates OptimisationProblem
//! \ingroup OptimisationSubModule
//! \name Type synonyms
//!@{
using ApproximateOptimisationProblem = OptimisationProblem<ApproximateTag>; //!< <p/>
using ValidatedOptimisationProblem = OptimisationProblem<ValidatedTag>; //!< <p/>
//!@}

//! \ingroup OptimisationSubModule
//! \name Data structures storing primal, dual, complementary and/or slack variables.
//!@{

//! \brief Data for the solution of an optimization problem, with primal variables of type \a X.
template<class X> struct PrimalData {
    Vector<X> x;
    //! \brief .
    PrimalData(Vector<X> x_) : x(x_) { }
    //! \brief .
    Vector<X> const& primal() const { return this->x; }
    //! \brief Convert to the primal variables.
    operator Vector<X> () const { return this->x; }
};
using ApproximatePrimalData = PrimalData<FloatDPApproximation>;
using ValidatedPrimalData = PrimalData<FloatDPBounds>;

//! \brief Data for the solution of an optimization problem, with primal variables of type \a X and dual variables of type \a Y.
template<class X> struct PrimalDualData : public PrimalData<X> {
    Vector<X> y;
    //! \brief .
    PrimalDualData(Vector<X> x_, Vector<X> y_) : PrimalData<X>(x_), y(y_) { }
    //! \brief .
    Vector<X> const& dual() const { return this->y; }
};
using ApproximatePrimalDualData = PrimalDualData<FloatDPApproximation>;
using ValidatedPrimalDualData = PrimalDualData<FloatDPBounds>;

//! \brief Data for the solution of an optimization problem, with primal variables of type \a Vector<X>, dual variables of type \a Vector<X>, and complementary variables of type \a Vector<X>.
template<class X> struct PrimalDualComplementaryData : public PrimalDualData<X> {
    Vector<X> z;
    //! \brief .
    PrimalDualComplementaryData(Vector<X> x_, Vector<X> y_, Vector<X> z_) : PrimalDualData<X>(x_,y_), z(z_) { }
    //! \brief .
    Vector<X> const& complementary() const { return this->z; }
};
using ApproximatePrimalDualComplementaryData = PrimalDualComplementaryData<FloatDPApproximation>;
using ValidatedPrimalDualComplementaryData = PrimalDualComplementaryData<FloatDPBounds>;

//! \brief Data for the solution of an optimization problem, with primal variables of type \a Vector<X> and slack variables of type \a Vector<X>.
template<class X> struct SlackPrimalData : public PrimalData<X> {
    Vector<X> w;
    //! \brief .
    SlackPrimalData(Vector<X> w_, Vector<X> x_) : PrimalData<X>(x_), w(w_) { }
    //! \brief .
    Vector<X> const& slack() const { return this->w; }
};
using ApproximateSlackPrimalData = SlackPrimalData<FloatDPApproximation>;
using ValidatedSlackPrimalData = SlackPrimalData<FloatDPBounds>;

//! \brief Data for the solution of an optimization problem, with primal variables of type \a Vector<X> dual variables of type \a Vector<X>, slack variables of type \a Vector<X> and complementary variables of type \a Vector<X>.
template<class X> struct SlackPrimalDualComplementaryData : public PrimalDualComplementaryData<X> {
    //! \brief The slack variables.
    Vector<X> w;
    //! \brief .
    SlackPrimalDualComplementaryData(Vector<X> w_, Vector<X> x_, Vector<X> y_, Vector<X> z_) : PrimalDualComplementaryData<X>(x_,y_,z_), w(w_) { }
    //! \brief .
    Vector<X> const& slack() const { return this->w; }
};
using ApproximateSlackPrimalDualComplementaryData = SlackPrimalDualComplementaryData<FloatDPApproximation>;
using ValidatedSlackPrimalDualComplementaryData = SlackPrimalDualComplementaryData<FloatDPBounds>;

//!@}


//! \ingroup OptimisationSubModule EvaluationModule
//! Interface for methods to check feasibility of constraint systems.
class FeasibilityCheckerInterface {
  public:
    typedef Vector<ExactNumber> ExactVectorType;
    typedef Vector<ValidatedNumber> ValidatedVectorType;
    typedef Vector<ApproximateNumber> ApproximateVectorType;
  public:
    //! \brief Virtual destructor.
    virtual ~FeasibilityCheckerInterface() = default;
    //! \brief Create a dynamically-allocated copy.
    virtual FeasibilityCheckerInterface* clone() const = 0;

    //! \brief Tests if the point \a x is almost feasible, in that \f$x\in D\f$ and \f$g(x)\in N_\epsilon(C)\f$.
    virtual ApproximateKleenean almost_feasible_point(ValidatedFeasibilityProblem p,
                                                      ApproximateVectorType x, ApproximateNumber eps) const = 0;
    //! \brief Tests whether the point \a x is feasible.
    //! \details If there are non-algebraic equality constraints, then it is unlikely to find \f$x\f$ exactly atisfying constraints.
    //! For this reason, this method should only be used if all constraints are inequalities, so the feasible set has nonempty interior.
    virtual ValidatedKleenean is_feasible_point(ValidatedFeasibilityProblem p,
                                                ExactVectorType x) const = 0;
    //! \brief Tests whether the box \a X contains a feasible point.
    virtual ValidatedKleenean contains_feasible_point(ValidatedFeasibilityProblem p,
                                                      UpperBoxType X) const = 0;
    //! \brief Tests whether the validated point \a x contains an exact feasible point, or if the Lagrange multipliers \a y are a certificate of infeasibility.
    virtual ValidatedKleenean check_feasibility(ValidatedFeasibilityProblem p,
                                                ValidatedVectorType x, ExactVectorType y) const = 0;
    //! \brief Checks if the point \a x is definitely feasible.
    virtual Bool validate_feasibility(ValidatedFeasibilityProblem p,
                                      ValidatedVectorType x) const = 0;
    //! \brief Check whether the system \f$h(x)=0\f$ definitely has a solution consistent with \f$x\f$.
    virtual Bool validate_feasibility(ValidatedVectorMultivariateFunction h,
                                      ValidatedVectorType x) const = 0;
    //! \brief Tests if the feasibility problem is definitely unsolvable, using multipliers \a y and local centering point \a xa.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        ApproximateVectorType xa, ExactVectorType y) const = 0;
    //! \brief Tests if the feasibility problem is definitely unsolvable over the box \a X, using multipliers \a y.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        UpperBoxType X, ExactVectorType y) const = 0;
    //! \brief Tests if the feasibility problem is definitely unsolvable, using multipliers \a y.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        ExactVectorType y) const = 0;

    //! \brief Tests if the feasibility problem is definitely unsolvable.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p) const = 0;

};


//! \ingroup OptimisationSubModule EvaluationModule
//! Interface for nonlinear programming solvers.
class OptimiserInterface {
    using FLT=FloatDP;
  public:
    typedef FLT ExactNumericType;
    typedef Bounds<FLT> ValidatedNumericType;
    typedef Approximation<FLT> ApproximateNumericType;
    typedef Vector<FLT> ExactVectorType;
    typedef Vector<Bounds<FLT>> ValidatedVectorType;
    typedef Scalar<Bounds<FLT>> ValidatedScalarType;
    typedef Vector<Approximation<FLT>> ApproximateVectorType;
    typedef Scalar<Approximation<FLT>> ApproximateScalarType;
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

    //! \brief Tests is the general nonlinear feasibility problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ApproximateKleenean feasible(ApproximateFeasibilityProblem p) const = 0;
};

//! \ingroup OptimisationSubModule
//! Common routines for nonlinear programming
class FeasibilityChecker
    : public virtual FeasibilityCheckerInterface
{
    template<class FLT> using Exact = FLT;
    using FLT=FloatDP;
    static constexpr DP pr=dp;
  public:
    virtual FeasibilityChecker* clone() const override;

    //! \brief Tests if the point \a x is almost feasible, in that \f$x\in D\f$ and \f$g(x)\in N_\epsilon(C)\f$.
    virtual ApproximateKleenean almost_feasible_point(ValidatedFeasibilityProblem p,
                                                      ApproximateVectorType x, ApproximateNumber eps) const override;
    //! \brief Tests whether the point \a x is feasible.
    virtual ValidatedKleenean is_feasible_point(ValidatedFeasibilityProblem p,
                                                ExactVectorType x) const override;
    //! \brief Tests whether the validated point \a x contains an exact feasible point.
    //! \details First tests whether \f$[x]\f$ is an element of \f$D\f$.
    //! If definitely not, returns \a false.
    //! If overlap, restricts \f$[x]\f$ to \f$D\f$.
    //! Then computes \f$[w]=g([x])\f$ and tests whther this is an element of \f$C\f$.
    //! If definitely not, returns \a false.
    //! For any other component for which \f$[w]_j\f$ is not a subset of \f$C_j\f$, introduce an equality constraint which needs to be satisfied.
    //! Check these equality constraints using an interval Newton contractor.
    virtual ValidatedKleenean contains_feasible_point(ValidatedFeasibilityProblem p,
                                                      UpperBoxType x) const override;
    //! \brief Tests whether the validated point \a x contains an exact feasible point, or if the Lagrange multipliers \a y are a certificate of infeasibility.
    //! \details First checks feasibility of \f$[x]\f$ using \ref is_feasible_point.
    //! If not definitely feasible, uses Lagrange multipliers \f$y\f$ to make a linear combination of constraint functions
//     //! and tests \f$\sum y_j g_j(D)\f$ overlaps \f$\sum y_j C_j\f$, possibly using power series or linearisation to increase accuracy.
    virtual ValidatedKleenean check_feasibility(ValidatedFeasibilityProblem p,
                                                ValidatedVectorType x, ExactVectorType y) const override;

    //! \brief Checks if the point \a x is definitely feasible.
    virtual Bool validate_feasibility(ValidatedFeasibilityProblem p,
                                      ValidatedVectorType x) const override;
    //! \brief Check whether the system \f$h(x)=0\f$ definitely has a solution consistent with \f$x\f$.
    virtual Bool validate_feasibility(ValidatedVectorMultivariateFunction h,
                                      ValidatedVectorType x) const override;
    //! \brief Tests if the feasibility problem is definitely unsolvable, using multipliers \a y and local centering point \a xa.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        ApproximateVectorType xa, ExactVectorType y) const override;
    //! \brief Tests if the feasibility problem is definitely unsolvable over the box \a X, using multipliers \a y.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        UpperBoxType X, ExactVectorType y) const override;
    //! \brief Tests if the feasibility problem is definitely unsolvable, using multipliers \a y.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        ExactVectorType y) const override;
    //! \brief Tests if the feasibility problem is definitely unsolvable.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p) const override;

    friend OutputStream& operator<<(OutputStream& os, FeasibilityChecker const& fc);
  private: public:
    ValidatedKleenean check_feasibility(ValidatedFeasibilityProblem p,
                                        Vector<Bounds<FLT>> x, Vector<Bounds<FLT>> y) const;

    Bool validate_feasibility(ValidatedFeasibilityProblem p,
                                      Vector<Bounds<FLT>> x) const;
    Bool validate_feasibility(ValidatedVectorMultivariateFunction h,
                                      Vector<Bounds<FLT>> x) const;

    Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        UpperBoxType X, Vector<Bounds<FLT>> y) const;
    Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        Vector<Approximation<FLT>> xa, Vector<Bounds<FLT>> y) const;
    Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        Vector<Bounds<FLT>> y) const;

    Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        UpperBoxType X, Vector<Exact<FLT>> y) const;
    Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        Vector<Approximation<FLT>> xa, Vector<Exact<FLT>> y) const;
    Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        Vector<Exact<FLT>> y) const;
};


//! \ingroup OptimisationSubModule
//! Common routines for nonlinear programming
class OptimiserBase
    : public virtual OptimiserInterface
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

    virtual ApproximateKleenean feasible(ApproximateFeasibilityProblem p) const override = 0;

};


//! \ingroup OptimisationSubModule
//! Common routines for nonlinear programming
class ApproximateOptimiserBase
    : public OptimiserBase
{
};

struct InteriorPointOptimiserProperties {
    const double VALUE_TOLERANCE=1e-8;
    const double STATE_TOLERANCE=1e-8;
    const CounterType MAXIMUM_STEPS=24;
};


class InteriorPointOptimiserBase
    : public OptimiserBase
{
    InteriorPointOptimiserProperties _properties;
  public:
    //! \brief Compute a \em local optimum of nonlinear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C .\f$.
    //! \pre The domain \f$D\f$ is bounded and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A box \f$X\f$ which definitely contains a feasible point, and contains a local optimum.
    virtual Vector<ValidatedNumericType> minimise(ValidatedOptimisationProblem p) const override;
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ValidatedFeasibilityProblem) const override;
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
    virtual ApproximateKleenean feasible(ApproximateFeasibilityProblem p) const override;
  public:
    virtual Void feasibility_step(ValidatedFeasibilityProblem p,
                                  ValidatedVectorType& w, ValidatedVectorType& x) const;
    virtual Void feasibility_step(ApproximateFeasibilityProblem p,
                                  ApproximateVectorType& w, ApproximateVectorType& x, FloatDPApproximation& mu) const;
    virtual Void feasibility_step(ApproximateFeasibilityProblem p,
                                  ApproximateVectorType& w, ApproximateVectorType& x, ApproximateVectorType& y) const;
};


struct FeasibilityCertificate
{
    using ValidatedVectorType = OptimiserInterface::ValidatedVectorType;
    using ExactVectorType = OptimiserInterface::ExactVectorType;

    ValidatedVectorType x; ExactVectorType y;

    FeasibilityCertificate(ValidatedVectorType x_, ExactVectorType y_) : x(x_), y(y_) { }
    FeasibilityCertificate(Pair<ValidatedVectorType,ExactVectorType>const& pr)
        : FeasibilityCertificate(std::get<0>(pr),std::get<1>(pr)) { }
    friend OutputStream& operator<<(OutputStream& os, FeasibilityCertificate& fc) {
        return os << "FeasibilityCertificate(x=" << fc.x << ", y=" << fc.y << ")"; }
};

struct LocalInfeasibilityCertificate
{
    using ExactVectorType = OptimiserInterface::ExactVectorType;
    ExactBoxType B; ExactVectorType y;
    ExactBoxType bounds() const { return B; }
    ExactVectorType dual() const { return y; }

    LocalInfeasibilityCertificate(ExactBoxType B_, ExactVectorType y_) : B(B_), y(y_) { }
    friend OutputStream& operator<<(OutputStream& os, const LocalInfeasibilityCertificate& lifc) {
        return os << "{B=" << lifc.B << ", y=" << lifc.y << "}"; }
};

struct InfeasibilityCertificate
    : public List<LocalInfeasibilityCertificate>
{
};

struct LocalOptimalityCertificate
{
    using ValidatedNumericType = OptimiserInterface::ValidatedNumericType;
    using ValidatedVectorType = OptimiserInterface::ValidatedVectorType;
    using ExactVectorType = OptimiserInterface::ExactVectorType;

    ValidatedNumericType v; ValidatedVectorType x; ExactVectorType y;

    ValidatedNumericType value() const { return v; }
    ValidatedVectorType primal() const { return x; }
    ExactVectorType dual() const { return y; }

    LocalOptimalityCertificate(ValidatedNumericType v_, ValidatedVectorType x_, ExactVectorType y_)
        : v(v_), x(x_), y(y_) { }
    LocalOptimalityCertificate(Tuple<ValidatedNumericType,ValidatedVectorType,ExactVectorType>const& tup)
        : LocalOptimalityCertificate(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup)) { }
    friend OutputStream& operator<<(OutputStream& os, const LocalOptimalityCertificate& loc) {
        return os << "{v=" << loc.v << ", x=" << loc.x << ", y=" << loc.y << "}"; }
};

struct OptimalityCertificate
    : public List<LocalOptimalityCertificate>
{
};




//! \ingroup OptimisationSubModule
//! Solver for nonlinear programming problems using the infeasible interior point method with primal (\f$x\f$), dual (\f$y\f$) and complementary (\f$z\f$)  variables.
//! \details Relies on an engine minimising \f$f(x)\f$ for \f$x\in D^\circ\f$ subject to \f$g(x)=w\f$ with \f$w\in C^\circ\f$.
//!
//! To check feasibility, we maximise \f$\mu \sum_{i}\bigl(\log(x_i-\unl{d}_i)+\log(\ovl{d}_i-x_i)\bigr) + \mu \sum_{j}\bigl(\log(w_j-\unl{c}_j)+\log(\ovl{c}_j-w_j)\bigr)\f$ under the constraints \f$g(x)-w=0\f$.
//! Introducing the Lagrange multipliers \f$y_j = \bigl(1/(w_j-\unl{c}_j)-1/(\ovl{c}_j-x_j)\bigr)\f$ and \f$z_i = 1/(x_i-\unl{d}_i)-1/(\ovl{d}_i-x_i)\f$, we find optimality conditions \f$ \sum_{j} y_j \nabla{g_j}(x) + z = 0\f$, \f$g(x)-w=0\f$.
//! These are exactly the central path equations with \f$\mu=1\f$ for the problem \f$\minimise f(x)\equiv 0\f$ subject to \f$x\in D\f$ and \f$g(x)\in C\f$, so we can use an optimisation step to improve feasibility.
class InteriorPointOptimiser
    : public OptimiserBase
{
    InteriorPointOptimiserProperties _properties;
  public:
    template<class T> struct SlackPrimalDualComplementaryData;
    template<class T> struct SlackPrimalDualComplementaryMatrix;
    struct StepData;

  public:
    virtual InteriorPointOptimiser* clone() const override;

    using OptimiserBase::minimise;
    using OptimiserBase::feasible;

  public:
    //! \brief Construct with default accuracy parameters.
    InteriorPointOptimiser();

    virtual Vector<ValidatedNumericType> minimise(ValidatedOptimisationProblem p) const override;
    virtual ValidatedKleenean feasible(ValidatedFeasibilityProblem p) const override;

    //! \brief Compute an approximate \em local optimum of nonlinear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C.\f$.
    virtual Vector<ApproximateNumericType> minimise(ApproximateOptimisationProblem p) const override;
    //! \brief Compute an approximate \em local optimum of nonlinear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C.\f$.
    virtual ApproximateKleenean feasible(ApproximateFeasibilityProblem p) const override;

        //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
        //! hotstarting the method with the primal and dual variables.
        Tuple<ValidatedScalarType,ValidatedVectorType,ValidatedVectorType>
        minimise_hotstarted(const ValidatedOptimisationProblem& p,
                            const ApproximateVectorType& x0, const ApproximateVectorType& y0) const;

        //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
        //! hotstarting the method with the primal and dual variables.
        Tuple<ValidatedKleenean,ValidatedVectorType,ValidatedVectorType>
        feasible_hotstarted(const ValidatedFeasibilityProblem& p,
                            const ApproximateVectorType& x0, const ApproximateVectorType& y0) const;


    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal and dual variables.
    virtual Tuple<ApproximateScalarType,ApproximateVectorType,ApproximateVectorType>
    minimise_hotstarted(const ApproximateOptimisationProblem& p,
                        const ApproximateVectorType& x0, const ApproximateVectorType& y0) const;
    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal and dual variables.
    virtual Tuple<ApproximateKleenean,ApproximateVectorType,ApproximateVectorType>
    feasible_hotstarted(const ApproximateFeasibilityProblem& p,
                        const ApproximateVectorType& x0, const ApproximateVectorType& y0) const;


    //! \brief .
    SlackPrimalDualComplementaryData<FloatDPApproximation>
    minimisation_update(const ApproximateOptimisationProblem& p,
                        FloatDPApproximationVector& w, FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximationVector& z,
                        FloatDPApproximation& mu) const;

    //! \brief .
    Void minimisation_step(const ApproximateOptimisationProblem& p,
                           StepData& d) const;
    //! \brief .
    Void minimisation_step(const ApproximateOptimisationProblem& p,
                           FloatDPApproximationVector& w, FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximationVector& z,
                           FloatDPApproximation& mu) const;
    //! \brief .
    Void feasibility_step(const ApproximateFeasibilityProblem& p,
                          FloatDPApproximationVector& w, FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximationVector& z) const;

    //! \brief <p/>
    StepData initial_step_data(const ApproximateFeasibilityProblem& p) const;
    //! \brief <p/>
    StepData initial_step_data_hotstarted(const ApproximateFeasibilityProblem& p,
                                          const FloatDPApproximationVector& x, const FloatDPApproximationVector& y) const;

    //! \brief .
    FloatDPApproximationVector compute_dual(const ApproximateBoxType& D, const FloatDPApproximationVector& x, const FloatDPApproximation& mu) const;
    //! \brief .
    FloatDPApproximationVector compute_x(const ApproximateFeasibilityProblem& p) const;
    //! \brief .
    FloatDPApproximationVector compute_y(const ApproximateFeasibilityProblem& p, const FloatDPApproximationVector& x, const FloatDPApproximation& mu) const;
    //! \brief .
    FloatDPApproximationVector compute_w(const ApproximateFeasibilityProblem& p,
                                         const FloatDPApproximationVector& x, const FloatDPApproximationVector& y, const FloatDPApproximation& mu) const;
    //! \brief .
    FloatDPApproximationVector compute_w(const ApproximateFeasibilityProblem& p,
                                         const FloatDPApproximationVector& x, const FloatDPApproximation& mu) const;
    //! \brief .
    FloatDPApproximationVector compute_z(const ApproximateFeasibilityProblem& p,
                                         const FloatDPApproximationVector& x, const FloatDPApproximation& mu) const;
    //! \brief .
    FloatDPApproximation compute_t(const ApproximateFeasibilityProblem& p,
                                   const FloatDPApproximationVector& x) const;
    //! \brief .
    FloatDPApproximation compute_mu(const ApproximateFeasibilityProblem& p,
                                    const FloatDPApproximationVector& x, const FloatDPApproximationVector& y) const;
    Void compute_tz(const ApproximateFeasibilityProblem& p,
                    FloatDPApproximationVector& x, FloatDPApproximation& t, FloatDPApproximationVector& z) const;

    friend OutputStream& operator<<(OutputStream& os, InteriorPointOptimiser const& opt);
  private: public:
    //! \brief Checks whether the point \f$x\f$ is feasible, or the multipliers \f$y\f$ provide a certificate of infeasibility based around \f$x\f$.
    Tuple<ValidatedKleenean,ValidatedVectorType,ValidatedVectorType>
    check_feasibility(const ValidatedFeasibilityProblem& p,
                      const ApproximateVectorType& x, const ApproximateVectorType& y) const;

    //! \brief Tests whether the validated point \a x contains a locally optimal point, using Lagrange multipliers \a y.
    Tuple<ValidatedKleenean,ValidatedVectorType,ValidatedVectorType>
    check_minimality(const ValidatedOptimisationProblem& p,
                     const FloatDPApproximationVector& w, const FloatDPApproximationVector& x, const FloatDPApproximationVector& y, const FloatDPApproximationVector& z) const;

    //! \brief Tests whether the validated point \a x contains a locally optimal point, using Lagrange multipliers \a y.
    Tuple<ValidatedKleenean,ValidatedVectorType,ValidatedVectorType>
    check_minimality(const ValidatedOptimisationProblem& p,
                     const FloatDPBoundsVector& w, const FloatDPBoundsVector& x, const FloatDPBoundsVector& y, const FloatDPBoundsVector& z) const;

    Tuple<ValidatedScalarType,ValidatedVectorType,ValidatedVectorType>
    nonsplitting_minimise_hotstarted(const ValidatedOptimisationProblem& p,
                                     const ApproximateVectorType& x0, const ApproximateVectorType& y0) const;

    Tuple<ValidatedKleenean,ValidatedVectorType,ValidatedVectorType>
    nonsplitting_feasible_hotstarted(const ValidatedFeasibilityProblem& p,
                                     const ApproximateVectorType& x0, const ApproximateVectorType& y0) const;


    Variant<OptimalityCertificate,InfeasibilityCertificate>
    splitting_minimise_hotstarted(const ValidatedOptimisationProblem& p,
                                  UpperBoxType B, FloatDPApproximationVector x, FloatDPApproximationVector y) const;

    Variant<FeasibilityCertificate,InfeasibilityCertificate>
    splitting_feasible_hotstarted(const ValidatedFeasibilityProblem& p,
                                  UpperBoxType B, FloatDPApproximationVector x, FloatDPApproximationVector y) const;
};


struct InteriorPointOptimiser::StepData {
    FloatDPApproximationVector w,x,y,z; FloatDPApproximation mu;

    StepData(SizeType m, SizeType n, DP pr)
        : w(m,pr), x(n,pr), y(m,pr), z(n,pr), mu(pr) { }
    StepData(FloatDPApproximationVector w_, FloatDPApproximationVector x_, FloatDPApproximationVector y_, FloatDPApproximationVector z_, FloatDPApproximation mu_)
        : w(w_), x(x_), y(y_), z(z_), mu(mu_) { }
    friend OutputStream& operator<<(OutputStream& os, StepData const& d) {
        return os << "InteriorPointOptimiser::StepData(w=" << d.w << ", x=" << d.x << ", y=" << d.y << ", z=" << d.z << ", mu=" << d.mu << ")"; }
};


//! \ingroup OptimisationSubModule
//! Solver for nonlinear programming problems using the infeasible interior point method with primal (\f$x\f$), dual (\f$y\f$), slack (\f$w\f$) and complementary (\f$z\f$) variables.
//! The dual and complementary variables are split into lower and upper versions, corresponding to lower and upper constraint bounds.
//! \details The complementary variables satisfy the central path equations \f$(x-\unl{d})\unl{z} = \mu\f$ and \f$(\ovl{d}-x)\ovl{z} = \mu\f$ with \f$\unl{z},\ovl{z}>0\f$, and can be combined into \f$z = \ovl{z}-\unl{z}\f$ satisfying \f$z = \mu\bigl(1/(\ovl{d}-x)-1/(x-\unl{d})\bigr)\f$, or equivalently \f$(x-\unl{d})(\ovl{d}-x)z + (\unl{d}+\ovl{d}-2x) \mu = 0\f$. Similarly, \f$(w-\unl{c})\unl{y} = \mu\f$ and \f$(\ovl{c}-w)\ovl{y} = \mu\f$ where \f$w=g(x)\f$.
class PrimalSplitDualInteriorPointOptimiser
    : public InteriorPointOptimiser
{
    virtual PrimalSplitDualInteriorPointOptimiser* clone() const override;
  public:
    struct StepData;

    //! \brief .
    Tuple<FloatDPApproximation,FloatDPApproximationVector,FloatDPApproximationVector>
    virtual minimise_hotstarted(const ApproximateOptimisationProblem& p,
                                const FloatDPApproximationVector& x0, const FloatDPApproximationVector& y0) const override;

    //! \brief .
    Void minimisation_step(const ApproximateOptimisationProblem& p,
                           FloatDPApproximationVector& w, FloatDPApproximationVector& x,
                           FloatDPApproximationVector& yl, FloatDPApproximationVector& yu, FloatDPApproximationVector& zl, FloatDPApproximationVector& zu,
                           FloatDPApproximation mu) const;
};

//! \relates PrimalSplitDualInteriorPointOptimiser
struct PrimalSplitDualInteriorPointOptimiser::StepData {
    FloatDPApproximationVector w,x,yl,yu,zl,zu; FloatDPApproximation mu;

    StepData(SizeType m, SizeType n, DP pr)
        : w(m,pr), x(n,pr), yl(m,pr), yu(m,pr), zl(n,pr), zu(n,pr), mu(pr) { }
    StepData(FloatDPApproximationVector w_, FloatDPApproximationVector x_, FloatDPApproximationVector yl_, FloatDPApproximationVector yu_, FloatDPApproximationVector zl_, FloatDPApproximationVector zu_, FloatDPApproximation mu_)
        : w(w_), x(x_), yl(yl_), yu(yu_), zl(zl_), zu(zu_), mu(mu_) { }
};

//! \ingroup OptimisationSubModule
//! An interior-point optimiser using only primal (\f$x\f$) and dual (\f$y\f$) variables.
//! In particular, there are no complementary variables dual to the constraints \f$\unl{d}_i\leq x_i \leq \ovl{d}_i\f$, or slack variables \f$w=g(x)\f$.
class PrimalDualOnlyInteriorPointOptimiser
    : public InteriorPointOptimiser
{
    virtual PrimalDualOnlyInteriorPointOptimiser* clone() const override;
  public:
    Tuple<ApproximateKleenean,FloatDPApproximationVector,FloatDPApproximationVector>
    virtual feasible_hotstarted(const ApproximateFeasibilityProblem& p,
                                const FloatDPApproximationVector& x0, const FloatDPApproximationVector& y0) const override;

    Void feasibility_step(const ApproximateFeasibilityProblem& p,
                          FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximation& t) const;
};

//! \ingroup OptimisationSubModule
//! \brief Solver for nonlinear programming problems using infeasible interior point methods.
//! Uses primal variables \f$x\f$, slack variables \f$w=g(x)\f$, dual variables \f$y\f$ for the constraints \f$w \in C\f$ and complementary variables \f$z\f$ for the constraints \f$x\in D\f$.
//! \details The dual variables \f$y\f$ are unconstrained Lagrange multipliers for \f$y\cdot(g(x)-w)=0\f$.
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
    //! \brief Tests is the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ValidatedFeasibilityProblem p) const override;

    //! \brief Compute an approximate \em local optimum of nonlinear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C.\f$.
    virtual Vector<ApproximateNumericType> minimise(ApproximateOptimisationProblem p) const override;
    //! \brief Tests is the general nonlinear feasibility problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ApproximateKleenean feasible(ApproximateFeasibilityProblem p) const override;

    friend OutputStream& operator<<(OutputStream& os, InfeasibleInteriorPointOptimiser const& opt);
  public:
    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the overall primal and dual variables.
    Pair<ValidatedKleenean,FloatDPApproximationVector> feasible_hotstarted(ValidatedFeasibilityProblem p,
                                                                           const PrimalDualData& wxy0) const;

    Void setup_feasibility(const ApproximateFeasibilityProblem& p,
                           StepData& stp) const;
    Void minimisation_step(const ApproximateOptimisationProblem& p,
                           StepData& stp) const;
};

struct InfeasibleInteriorPointOptimiser::PrimalDualData {
    PrimalDualData() : PrimalDualData(0u,0u,dp) { }
    PrimalDualData(SizeType m, SizeType n, DP pr) : w(m,pr), x(n,pr), y(m,pr) { }
    FloatDPApproximationVector w,x,y;
};

//! \relates InfeasibleInteriorPointOptimiser
struct InfeasibleInteriorPointOptimiser::StepData : public InfeasibleInteriorPointOptimiser::PrimalDualData {
    StepData() : StepData(0u,0u,dp) { }
    StepData(SizeType m, SizeType n, DP pr)
        : PrimalDualData(m,n,pr), vl(m,pr), wl(m,pr), xl(n,pr), zl(n,pr), vu(m,pr), wu(m,pr), xu(n,pr), zu(n,pr), mu(pr) { }
    FloatDPApproximationVector vl,wl,xl,zl,vu,wu,xu,zu; FloatDPApproximation mu;
};





//! \ingroup OptimisationSubModule
//! \brief An optimiser using rigorous interval-arithmetic methods based on the John conditions for feasibility problems.
//! <br>Currently only provides rigorous solvers for feasibility problems; optimisation problems are not supported.
class IntervalOptimiser
    : public InteriorPointOptimiser
{
  private: public:
    virtual IntervalOptimiser* clone() const override;
    //! \brief .
    virtual ValidatedKleenean feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const;
    //! \brief .
    Void feasibility_step(const FloatDPVector& xl, const FloatDPVector& xu, const ValidatedVectorMultivariateFunction& h,
                          FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& zl, FloatDPBoundsVector zu, FloatDPBounds& mu) const;
};


//! \ingroup OptimisationSubModule
//! \brief An optimiser based on the InteriorPointOptimiser. \deprecated
class ApproximateOptimiser
    : public InteriorPointOptimiser
{
  private: public:
    virtual ApproximateOptimiser* clone() const override;
    //! \brief .
    virtual ValidatedKleenean feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const;
    //! \brief .
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
