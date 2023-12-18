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
template<class P> using ConstraintType = Constraint<ScalarMultivariateFunction<P>,Number<P>>; //!< <p/>

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
    FeasibilityProblem(BoxType<P> D_, VectorMultivariateFunction<P> g_, BoxType<P> C_);
    //! <p/>
    FeasibilityProblem(BoxType<P> D_, const List<ConstraintType<P>>& cl);
    //! <p/>
    template<class PP> requires Convertible<PP,P> FeasibilityProblem(const FeasibilityProblem<PP>& p)
        : FeasibilityProblem(p.D,p.g,p.C) { }

    SizeType number_of_variables() const { return this->g.argument_size(); }
    SizeType number_of_constraints() const { return this->g.result_size(); }
  private:
    static VectorMultivariateFunction<P> _function(SizeType as, const List<ConstraintType<P>>& cl);
    static BoxType<P> _bounds(const List<ConstraintType<P>>& cl);
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
//! \brief Data for the optimisation problem \f$\text{minimise } f(x) \text{ subject to } x\in D \text{ and } g(x) \in C\f$, where \f$D\subset\R^n\f$ and \f$C\subset \R^m\f$ are boxes, and \f$f:\R^n\to\R\f$ and \f$g:\R^n\to\R^m\f$ are continuous.
template<class P> struct OptimisationProblem : public FeasibilityProblem<P> {
  public:
    ScalarMultivariateFunction<P> f;//!< <p/>
  public:
    //! <p/>
    OptimisationProblem(ScalarMultivariateFunction<P> f_, BoxType<P> D_, VectorMultivariateFunction<P> g_, BoxType<P> C_)
        : FeasibilityProblem<P>(D_,g_,C_), f(f_) { ARIADNE_PRECONDITION(f_.argument_size()==D_.dimension()); }
    //! <p/>
    OptimisationProblem(ScalarMultivariateFunction<P> f_, BoxType<P> D_, const List<ConstraintType<P>>& cl)
        : FeasibilityProblem<P>(D_,cl), f(f_) { ARIADNE_PRECONDITION(f_.argument_size()==D_.dimension()); }
    //! <p/>
    OptimisationProblem(ScalarMultivariateFunction<P> f_, BoxType<P> D_)
        : OptimisationProblem(f_,D_, List<ConstraintType<P>>()) { }
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

//! \brief Data for the solution of an optimization problem, with primal variables \a x of type \a X.
template<class X> struct PrimalData {
    Vector<X> x;
    //! \brief <p/>
    PrimalData(Vector<X> x_) : x(x_) { }
    //! \brief <p/>
    Vector<X> const& primal() const { return this->x; }
    //! \brief Convert to the primal variables.
    operator Vector<X> () const { return this->x; }
};
using ApproximatePrimalData = PrimalData<Approximation<FloatDP>>;
using ValidatedPrimalData = PrimalData<Bounds<FloatDP>>;


//! \brief Data for the solution of an optimization problem, with primal variables of type \a X and dual variables of type \a Y.
template<class X, class Y=X> struct PrimalDualData : public PrimalData<X> {
    Vector<Y> y;
    template<class PR> requires Constructible<X,PR> PrimalDualData(SizeType m, SizeType n, PR pr)
        : PrimalDualData(Vector<X>(n,pr),Vector<X>(m,pr)) { }
    PrimalDualData(SizeType m, SizeType n, Vector<X> r)
        : PrimalDualData(r[range(0,n)],r[range(n,m+n)]) { ARIADNE_DEBUG_PRECONDITION(r.size()==m+n); }
    //! \brief <p/>
    PrimalDualData(Vector<X> x_, Vector<Y> y_) : PrimalData<X>(x_), y(y_) { }
    //! \brief <p/>
    Vector<Y> const& dual() const { return this->y; }
    //! \brief <p/>
    friend OutputStream& operator<<(OutputStream& os, const PrimalDualData<X,Y>& d) {
        return os << "PrimalDualData(x=" << d.x << ", y=" << d.y << ")"; }
    friend PrimalDualData<X,Y> operator-(PrimalDualData<X,Y> d) {
        return PrimalDualData<X,Y>(-d.x,-d.y); }
};
using ApproximatePrimalDualData = PrimalDualData<Approximation<FloatDP>>;
using ValidatedPrimalDualData = PrimalDualData<Bounds<FloatDP>>;

//! \brief Data for the solution of an optimization problem, with primal variables \a x, dual variables \a y, and complementary variables \a z dual to bounds on \a x, all of type \a Vector<X>.
template<class X> struct PrimalDualComplementaryData : public PrimalDualData<X> {
    Vector<X> z;
    //! \brief <p/>
    PrimalDualComplementaryData(SizeType m, SizeType n, PrecisionType<X> pr)
        : PrimalDualComplementaryData(Vector<X>(n,pr), Vector<X>(m,pr), Vector<X>(n,pr)) { }
    PrimalDualComplementaryData(SizeType m, SizeType n, Vector<X> r)
        : PrimalDualComplementaryData(r[range(0,n)], r[range(n,m+n)], r[range(m+n,m+2*n)]) { ARIADNE_DEBUG_PRECONDITION(r.size()==m+2*n); }
    //! \brief <p/>
    PrimalDualComplementaryData(Vector<X> x_, Vector<X> y_, Vector<X> z_) : PrimalDualData<X>(x_,y_), z(z_) { }
    //! \brief <p/>
    Vector<X> const& complementary() const { return this->z; }

    static PrimalDualComplementaryData<X> disassemble(SizeType m, SizeType n, Vector<X> r);
    Vector<X> assemble() const;

    //! \brief <p/>
    friend OutputStream& operator<<(OutputStream& os, const PrimalDualComplementaryData<X>& d) {
        return os << "PrimalDualComplementaryData(x=" << d.x << ", y=" << d.y << ", z=" << d.z << ")"; }
};
using ApproximatePrimalDualComplementaryData = PrimalDualComplementaryData<Approximation<FloatDP>>;
using ValidatedPrimalDualComplementaryData = PrimalDualComplementaryData<Bounds<FloatDP>>;


//! \brief Data for the solution of an optimization problem, with primal variables \a x and slack variables \a w of type \a X.
template<class X> struct SlackPrimalData : public PrimalData<X> {
    Vector<X> w;
    //! \brief <p/>
    SlackPrimalData(Vector<X> w_, Vector<X> x_) : PrimalData<X>(x_), w(w_) { }
    //! \brief <p/>
    Vector<X> const& slack() const { return this->w; }
};
using ApproximateSlackPrimalData = SlackPrimalData<Approximation<FloatDP>>;
using ValidatedSlackPrimalData = SlackPrimalData<Bounds<FloatDP>>;

//! \brief Data for the solution of an optimization problem, with primal variables \a x, slack variables \a w approximating \a g(x), and dual variables \a y of type \a X.
template<class X> struct SlackPrimalDualData : public PrimalDualData<X> {
    Vector<X> w;
    SlackPrimalDualData(SizeType m, SizeType n, PrecisionType<X> pr)
        : PrimalDualData<X>(Vector<X>(n,pr),Vector<X>(m,pr)), w(m,pr) { }
    //! \brief <p/>
    SlackPrimalDualData(Vector<X> w_, Vector<X> x_, Vector<X> y_) : PrimalDualData<X>(x_,y_), w(w_) { }
    //! \brief <p/>
    Vector<X> const& slack() const { return this->w; }
};
using ApproximateSlackPrimalDualData = SlackPrimalDualData<Approximation<FloatDP>>;
using ValidatedSlackPrimalDualData = SlackPrimalDualData<Bounds<FloatDP>>;

//! \brief Data for the solution of an optimization problem, with primal variables \a x, dual variables \a y, slack variables \a w and complementary variables \a z, all of type \a Vector<X>.
template<class X> struct SlackPrimalDualComplementaryData : public PrimalDualComplementaryData<X> {
    Vector<X> w;

    //! \brief <p/>
    SlackPrimalDualComplementaryData(SizeType m, SizeType n, DP pr)
        : SlackPrimalDualComplementaryData(Vector<X>(m,pr),Vector<X>(n,pr),Vector<X>(m,pr),Vector<X>(n,pr)) { }
    //! \brief <p/>
    SlackPrimalDualComplementaryData(Vector<X> w_, Vector<X> x_, Vector<X> y_, Vector<X> z_)
        : PrimalDualComplementaryData<X>(x_,y_,z_), w(w_) { }
    //! \brief <p/>
    template<class XX> requires Convertible<XX,X> SlackPrimalDualComplementaryData(SlackPrimalDualComplementaryData<XX> const& d)
        : SlackPrimalDualComplementaryData(d.w,d.x,d.y,d.z) { }

    //! \brief <p/>
    static SlackPrimalDualComplementaryData<X> disassemble(SizeType m, SizeType n, Vector<X> r);
    //! \brief <p/>
    Vector<X> assemble() const;

    //! \brief <p/>
    friend OutputStream& operator<<(OutputStream& os, SlackPrimalDualComplementaryData<X>& d) {
        return os << "SlackPrimalDualComplementaryData(w=" << d.w << ", x=" << d.x << ", y=" << d.y << ", z=" << d.z << ")"; }

    friend SlackPrimalDualComplementaryData<X> operator-(SlackPrimalDualComplementaryData<X> const& d) {
        return SlackPrimalDualComplementaryData<X>(-d.w, -d.x, -d.y, -d.z); }
    friend SlackPrimalDualComplementaryData<X> operator-(SlackPrimalDualComplementaryData<X> const& d1, SlackPrimalDualComplementaryData<X> const& d2) {
        return SlackPrimalDualComplementaryData<X>(d1.w-d2.w, d1.x-d2.x, d1.y-d2.y, d1.z-d2.z); }
};
using ApproximateSlackPrimalDualComplementaryData = SlackPrimalDualComplementaryData<Approximation<FloatDP>>;
using ValidatedSlackPrimalDualComplementaryData = SlackPrimalDualComplementaryData<Bounds<FloatDP>>;

//! \brief Data for the solution of an optimization problem, with primal variables \a x, slack variables \a w approximating \a g(x), dual variables \a yl and \a yu for lower- and upper-bounds of \a w, and complementary variables \a zl, \a zu for lower- and upper-bounds of \a x.
template<class X> struct SlackPrimalSplitDualComplementaryData : public SlackPrimalData<X> {
    Vector<X> yl,yu,zl,zu;

    //! \brief <p/>
    SlackPrimalSplitDualComplementaryData(SizeType m, SizeType n, DP pr)
        : SlackPrimalSplitDualComplementaryData(Vector<X>(m,pr),Vector<X>(n,pr),Vector<X>(m,pr),Vector<X>(m,pr),Vector<X>(n,pr),Vector<X>(n,pr)) { }
    //! \brief <p/>
    SlackPrimalSplitDualComplementaryData(Vector<X> w_, Vector<X> x_, Vector<X> yl_, Vector<X> yu_, Vector<X> zl_, Vector<X> zu_)
        : SlackPrimalData<X>(w_,x_), yl(yl_), yu(yu_), zl(zl_), zu(zu_) { }
    //! \brief <p/>
    template<class XX> requires Convertible<XX,X> SlackPrimalSplitDualComplementaryData(SlackPrimalSplitDualComplementaryData<XX> const& d)
        : SlackPrimalSplitDualComplementaryData(d.w,d.x,d.yl,d.yu,d.zl,d.zu) { }

    //! \brief <p/>
    static SlackPrimalSplitDualComplementaryData<X> disassemble(SizeType m, SizeType n, Vector<X> r);
    //! \brief <p/>
    Vector<X> assemble() const;

    //! \brief <p/>
    friend OutputStream& operator<<(OutputStream& os, SlackPrimalSplitDualComplementaryData<X>& d) {
        return os << "SlackPrimalSplitDualComplementaryData(w=" << d.w << ", x=" << d.x << ", yl=" << d.yl << ", yu=" << d.yu << ", zl=" << d.zl << ", zu=" << d.zu << ")"; }
};
using ApproximateSlackPrimalSplitDualComplementaryData = SlackPrimalSplitDualComplementaryData<Approximation<FloatDP>>;
using ValidatedSlackPrimalSplitDualComplementaryData = SlackPrimalSplitDualComplementaryData<Bounds<FloatDP>>;



template<class X> struct FeasiblePrimalDualData : PrimalDualData<X> {
    ComparisonType<X> r;
    FeasiblePrimalDualData(ComparisonType<X> r_, Vector<X> x_, Vector<X> y_) : PrimalDualData<X>(x_,y_), r(r_) { }
    ComparisonType<X> const& is_feasible() const { return this->r; }

};

template<class X> struct ValuePrimalDualData : PrimalDualData<X> {
    X v;
    ValuePrimalDualData(X v_, Vector<X> x_, Vector<X> y_) : PrimalDualData<X>(x_,y_), v(v_) { }
    X const& value() const { return this->v; }
};

template<class X> struct SlackPrimalDualComplementaryMatrix;

//! \ingroup OptimisationSubModule
//! \brief Information providing evidence of feasibility.
struct FeasibilityCertificate
{
    using ValidatedVector = Vector<ValidatedNumber>;
    using ExactVector = Vector<ExactNumber>;

    ValidatedVector x; ExactVector y;

    FeasibilityCertificate(ValidatedVector x_, ExactVector y_) : x(x_), y(y_) { }
    FeasibilityCertificate(Pair<ValidatedVector,ExactVector>const& pr)
        : FeasibilityCertificate(std::get<0>(pr),std::get<1>(pr)) { }
    friend OutputStream& operator<<(OutputStream& os, FeasibilityCertificate& fc) {
        return os << "FeasibilityCertificate(x=" << fc.x << ", y=" << fc.y << ")"; }
};

//! \ingroup OptimisationSubModule
//! \brief Information providing evidence of local infeasibility.
struct LocalInfeasibilityCertificate
{
    using ExactVector = Vector<ExactNumber>;

    //! A subdomain \f$B\subset D\f$ such that \f$g(x)\not\in C\f$ for all \f$x\in C\f$. i.e. For which all (primal) variables are infeasible.
    ExactBoxType B;
    //! Dual variables (Lagrange multipliers) such that \f$\sum_{j=1} y_j g_j(x) \not\in \sum_{j} y_j C_j\f$ for all \f$x\in\B\f$.
    ExactVector y;

    ExactBoxType bounds() const { return B; }
    ExactVector dual() const { return y; }

    LocalInfeasibilityCertificate(ExactBoxType B_, ExactVector y_) : B(B_), y(y_) { }
    friend OutputStream& operator<<(OutputStream& os, const LocalInfeasibilityCertificate& lifc) {
        return os << "{B=" << lifc.B << ", y=" << lifc.y << "}"; }
};

//! \ingroup OptimisationSubModule
//! \brief Information providing evidence of global infeasibility.
struct InfeasibilityCertificate
    : public List<LocalInfeasibilityCertificate>
{
    using List<LocalInfeasibilityCertificate>::List;
};

//! \ingroup OptimisationSubModule
//! \brief Information providing evidence of local optimality (minimality).
struct LocalOptimalityCertificate
{
    using ValidatedVector = Vector<ValidatedNumber>;
    using ExactVector = Vector<ExactNumber>;

    ValidatedNumber v; //!< The optimal value.
    ValidatedVector x; //!< The optimal (primal) variables.
    ExactVector y; //!< The optimal dual variables (Lagrange multipliers).

    ValidatedNumber value() const { return v; }
    ValidatedVector primal() const { return x; }
    ExactVector dual() const { return y; }

    LocalOptimalityCertificate(ValidatedNumber v_, ValidatedVector x_, ExactVector y_)
        : v(v_), x(x_), y(y_) { }
    LocalOptimalityCertificate(Bounds<FloatDP> v_, Vector<Bounds<FloatDP>> x_, Vector<FloatDP> y_)
        : v(v_), x(x_), y(y_) { }
    LocalOptimalityCertificate(Tuple<ValidatedNumber,ValidatedVector,ExactVector>const& tup)
        : LocalOptimalityCertificate(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup)) { }
    friend OutputStream& operator<<(OutputStream& os, const LocalOptimalityCertificate& loc) {
        return os << "{v=" << loc.v << ", x=" << loc.x << ", y=" << loc.y << "}"; }
};

//! \ingroup OptimisationSubModule
//! \brief Information providing evidence of global optimality (minimality).
struct OptimalityCertificate
    : public List<LocalOptimalityCertificate>
{
    using List<LocalOptimalityCertificate>::List;
};

//!@}


//! \ingroup OptimisationSubModule
//! Interface for methods to check feasibility of constraint systems.
class FeasibilityCheckerInterface {
  public:
    typedef Vector<ExactNumber> ExactVector;
    typedef Vector<ValidatedNumber> ValidatedVector;
    typedef Vector<ApproximateNumber> ApproximateVector;
  public:
    //! \brief Virtual destructor.
    virtual ~FeasibilityCheckerInterface() = default;
    //! \brief Create a dynamically-allocated copy.
    virtual FeasibilityCheckerInterface* clone() const = 0;

    //! \brief Tests if the point \a x is almost feasible, in that \f$x\in D\f$ and \f$g(x)\in N_\epsilon(C)\f$.
    virtual ApproximateKleenean almost_feasible_point(ValidatedFeasibilityProblem p,
                                                      ApproximateVector x, ApproximateNumber eps) const = 0;
    //! \brief Tests whether the point \a x is feasible.
    //! \details If there are non-algebraic equality constraints, then it is unlikely to find \f$x\f$ exactly atisfying constraints.
    //! For this reason, this method should only be used if all constraints are inequalities, so the feasible set has nonempty interior.
    virtual ValidatedKleenean is_feasible_point(ValidatedFeasibilityProblem p,
                                                ExactVector x) const = 0;
    //! \brief Tests whether the box \a X contains a feasible point.
    virtual ValidatedKleenean contains_feasible_point(ValidatedFeasibilityProblem p,
                                                      UpperBoxType X) const = 0;
    //! \brief Tests whether the validated point \a x contains an exact feasible point, or if the Lagrange multipliers \a y are a certificate of infeasibility.
    virtual ValidatedKleenean check_feasibility(ValidatedFeasibilityProblem p,
                                                ValidatedVector x, ExactVector y) const = 0;
    //! \brief Checks if the point \a x is definitely feasible.
    virtual Bool validate_feasibility(ValidatedFeasibilityProblem p,
                                      ValidatedVector x) const = 0;
    //! \brief Check whether the system \f$h(x)=0\f$ definitely has a solution consistent with \f$x\f$.
    virtual Bool validate_feasibility(ValidatedVectorMultivariateFunction h,
                                      ValidatedVector x) const = 0;
    //! \brief Tests if the feasibility problem is definitely unsolvable, using multipliers \a y and local centering point \a xa.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        ApproximateVector xa, ExactVector y) const = 0;
    //! \brief Tests if the feasibility problem is definitely unsolvable over the box \a X, using multipliers \a y.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        UpperBoxType X, ExactVector y) const = 0;
    //! \brief Tests if the feasibility problem is definitely unsolvable, using multipliers \a y.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        ExactVector y) const = 0;

    //! \brief Tests if the feasibility problem is definitely unsolvable.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p) const = 0;

};


//! \ingroup OptimisationSubModule
//! \brief Interfaces for nonlinear programming solvers.
//! \details Specialised by \ref OptimiserInterface<ApproximateTag> and \ref OptimiserInterface<ValidatedTag>.
template<class P> class OptimiserInterface;

//! \ingroup OptimisationSubModule
//! \relates OptimiserInterface<P>
//! \name Type synonyms
//!@{
using ApproximateOptimiserInterface = OptimiserInterface<ApproximateTag>; //!< <p/>
using ValidatedOptimiserInterface = OptimiserInterface<ValidatedTag>; //!< <p/>
//!@}


//! \relates OptimiserInterface<P>
//! \brief Interface for approximate nonlinear programming solvers.
template<> class OptimiserInterface<ApproximateTag> {
  public:
    typedef Vector<ApproximateNumber> ApproximateVector;
  public:
    //! \brief Virtual destructor.
    virtual ~OptimiserInterface() = default;
    //! \brief Create a dynamically-allocated copy.
    virtual OptimiserInterface* clone() const = 0;

    //! \brief Solve the general nonlinear programming problem \f$\text{minimise } f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    virtual ApproximateVector minimise(ApproximateOptimisationProblem p) const = 0;
    //! \brief Approximatedly solve the general nonlinear programming problem \f$\text{minimise } f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    //! \return A point \f$x\f$ is approximately satisfies the feasibility constraints (it may not be close to a genuinely feasible point) and is approximately locally optimal.
    virtual ApproximateVector minimise(ApproximateScalarMultivariateFunction f, ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateBoxType C) const = 0;
    //! \brief Approximately solve the standard nonlinear programming problem \f$\text{minimise } f(x) \text{ such that } x\in D\ g(x)\leq 0  \text{ and } h(x) = 0\f$.
    virtual ApproximateVector minimise(ApproximateScalarMultivariateFunction f, ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateVectorMultivariateFunction h) const = 0;

    //! \brief Tests if the general nonlinear feasibility problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ApproximateKleenean feasible(ApproximateFeasibilityProblem p) const = 0;
    //! \brief Tests if the general nonlinear feasibility problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ApproximateKleenean feasible(ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateBoxType C) const = 0;
};

//! \relates OptimiserInterface<P>
//! \brief Interface for validated nonlinear programming solvers.
template<> class OptimiserInterface<ValidatedTag> {
  public:
    typedef Vector<ExactNumber> ExactVector;
    typedef Vector<ValidatedNumber> ValidatedVector;
  public:
    //! \brief Virtual destructor.
    virtual ~OptimiserInterface() = default;
    //! \brief Create a dynamically-allocated copy.
    virtual OptimiserInterface* clone() const = 0;

    //! \brief Solve the general nonlinear programming problem \f$\text{minimise } f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    virtual ValidatedVector minimise(ValidatedOptimisationProblem p) const = 0;
    //! \brief Solve the general nonlinear programming problem \f$\text{minimise } f(x) \text{ such that } x\in D \text{ and } g(x)\in C\f$.
    //! \pre The domain \f$D\f$ is bounded and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A point \f$x\f$ is approximately satisfies the feasibility constraints (it may not be which definitely contains a feasible point, and contains a local optimum.
    virtual ValidatedVector minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const = 0;
    //! \brief Solve the standard nonlinear programming problem \f$\text{minimise } f(x) \text{ such that } x\in D,\ g(x)\leq 0  \text{ and } h(x) = 0\f$.
    //! \pre The domain \f$D\f$ is bounded and has nonempty interior.
    //! \return A point \f$x\f$ is approximately satisfies the feasibility constraints (it may not be which definitely contains a feasible point, and contains a local optimum.
    virtual ValidatedVector minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const = 0;

    //! \brief Tests if the general nonlinear feasibility problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ValidatedFeasibilityProblem p) const = 0;
    //! \brief Tests if the general nonlinear feasibility problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const = 0;
    //! \brief Tests if the standard nonlinear feasibility problem \f$x\in D,\ g(x)\leq 0 \text{ and } h(x) = 0\f$ is feasible. Assumes \f$D\f$ is bounded with nonempty interior.
    //! \internal This is one of the simplest nonlinear programming problems, and is a good test case for new algorithms.
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const = 0;
};

//! \ingroup OptimisationSubModule
//! General-purpose routines for validated feasibility checking.
class FeasibilityChecker
    : public virtual FeasibilityCheckerInterface
{
    template<class FLT> using Exact = FLT;
    using FLT=FloatDP;
    using PR=DP;
    static constexpr PR pr=dp;
  public:
    virtual FeasibilityChecker* clone() const override;

    //! \brief Tests if the point \a x is almost feasible, in that \f$x\in D\f$ and \f$g(x)\in N_\epsilon(C)\f$.
    virtual ApproximateKleenean almost_feasible_point(ValidatedFeasibilityProblem p,
                                                      ApproximateVector x, ApproximateNumber eps) const override;
    //! \brief Tests whether the point \a x is feasible.
    virtual ValidatedKleenean is_feasible_point(ValidatedFeasibilityProblem p,
                                                ExactVector x) const override;
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
                                                ValidatedVector x, ExactVector y) const override;

    //! \brief Checks if the point \a x is definitely feasible.
    virtual Bool validate_feasibility(ValidatedFeasibilityProblem p,
                                      ValidatedVector x) const override;
    //! \brief Check whether the system \f$h(x)=0\f$ definitely has a solution consistent with \f$x\f$.
    virtual Bool validate_feasibility(ValidatedVectorMultivariateFunction h,
                                      ValidatedVector x) const override;
    //! \brief Tests if the feasibility problem is definitely unsolvable, using multipliers \a y and local centering point \a xa.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        ApproximateVector xa, ExactVector y) const override;
    //! \brief Tests if the feasibility problem is definitely unsolvable over the box \a X, using multipliers \a y.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        UpperBoxType X, ExactVector y) const override;
    //! \brief Tests if the feasibility problem is definitely unsolvable, using multipliers \a y.
    virtual Bool validate_infeasibility(ValidatedFeasibilityProblem p,
                                        ExactVector y) const override;
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
//! \brief Common routines for nonlinear programming
//! \details Specialised by \ref OptimiserBase<ApproximateTag> and \ref OptimiserBase<ValidatedTag>.
template<class P> class OptimiserBase;

//! \ingroup OptimisationSubModule
//! \relates OptimiserBase<P>
//! \name Type synonyms
//!@{
using ApproximateOptimiserBase = OptimiserBase<ApproximateTag>; //!< <p/>
using ValidatedOptimiserBase = OptimiserBase<ValidatedTag>; //!< <p/>
//!@}

//! \relates OptimiserBase<P>
//! \brief Common routines for approximate nonlinear programming
template<> class OptimiserBase<ApproximateTag>
    : public virtual OptimiserInterface<ApproximateTag>
{
  protected:
    template<class FLT> using Exact = FLT;
    using FLT=FloatDP;
    using PR=DP;
    static constexpr PR pr=dp;
  public:
    typedef Approximation<FLT> ApproximateNumberType;
    typedef Vector<Approximation<FLT>> ApproximateVectorType;
  protected:
    static const FLT zero;
    static const FLT one;
  public:
    virtual ApproximateVector minimise(ApproximateOptimisationProblem p) const override = 0;
    virtual ApproximateVector minimise(ApproximateScalarMultivariateFunction f, ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateBoxType C) const override;
    virtual ApproximateVector minimise(ApproximateScalarMultivariateFunction f, ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateVectorMultivariateFunction h) const override;

    virtual ApproximateKleenean feasible(ApproximateFeasibilityProblem p) const override = 0;
    virtual ApproximateKleenean feasible(ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateBoxType C) const override;
};

//! \relates OptimiserBase<P>
//! \brief Common routines for validated nonlinear programming
template<> class OptimiserBase<ValidatedTag>
    : public virtual OptimiserInterface<ValidatedTag>
{
  protected:
    template<class FLT> using Exact = FLT;
    using FLT=FloatDP;
    using PR=DP;
    static constexpr PR pr=dp;
  public:
    typedef FLT ExactNumberType;
    typedef Bounds<FLT> ValidatedNumberType;
    typedef Approximation<FLT> ApproximateNumberType;
    typedef Vector<FLT> ExactVectorType;
    typedef Vector<Bounds<FLT>> ValidatedVectorType;
    typedef Vector<Approximation<FLT>> ApproximateVectorType;
  protected:
    static const FLT zero;
    static const FLT one;
  public:
    virtual ValidatedVector minimise(ValidatedOptimisationProblem p) const override = 0;
    virtual ValidatedVector minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const override;
    virtual ValidatedVector minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const override;

    virtual ValidatedKleenean feasible(ValidatedFeasibilityProblem p) const override = 0;
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const override;
    virtual ValidatedKleenean feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const override;
};

using ApproximateOptimiserBase =  OptimiserBase<ApproximateTag>;
using ValidatedOptimiserBase =  OptimiserBase<ValidatedTag>;



struct InteriorPointOptimiserProperties {
    const double VALUE_TOLERANCE=1e-8;
    const double STATE_TOLERANCE=1e-8;
    const CounterType MAXIMUM_STEPS=24;
};

//! \ingroup OptimisationSubModule
//! \brief Solver for nonlinear programming problems based on a penalty-function approach
//!
//! For the feasibility problem \f$x\in D,\ g(x)\in C,\ h(x)=0\f$ where \f$D\f$ is bounded and \f$D,C\f$ have nonempty interiors.
//! Introduce slack variable \f$w=g(x)\f$, and minimise \f[ \sum_{j} (g_j(x)-w_j)^2 + \sum_k h_k(x)^2 \f] with \f$x\in D\f$ and \f$w\in C\f$.
//! Since the minimiser is not unique, we need to add penalty terms \f$-\mu/2(\log(x_u-x)+log(x-x_l))\f$ and \f$-\nu/2(\log(w_u-w)+log(w-w_l))\f$.
//! It suffices to add a penalty in \f$x\f$.
class PenaltyFunctionOptimiser
    : public ApproximateOptimiserBase
{
  public:
    using ApproximateOptimiserBase::minimise;
    using ApproximateOptimiserBase::feasible;
  public:
    virtual PenaltyFunctionOptimiser* clone() const override;
    virtual ApproximateVector minimise(ApproximateOptimisationProblem p) const override;
    virtual ApproximateKleenean feasible(ApproximateFeasibilityProblem p) const override;
  public:
    virtual Void feasibility_step(ApproximateFeasibilityProblem p,
                                  Vector<Approximation<FLT>>& w, Vector<Approximation<FLT>>& x, Approximation<FLT>& mu) const;
    virtual Void feasibility_step(ApproximateFeasibilityProblem p,
                                  Vector<Approximation<FLT>>& w, Vector<Approximation<FLT>>& x, Vector<Approximation<FLT>>& y) const;
  private:
    ValidatedKleenean feasible(ValidatedFeasibilityProblem p) const;
};


//! \brief Common functionality for optimisers using interior-point methods.
class InteriorPointOptimiserBase
    : public ApproximateOptimiserBase
{
  protected:
    struct StepData;
  protected:
    InteriorPointOptimiserProperties _properties;
  public:
    //! \brief Compute an approximate \em local optimum of nonlinear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C.\f$.
    virtual ApproximateVector minimise(ApproximateOptimisationProblem p) const override;
    //! \brief Compute an approximate \em local optimum of nonlinear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C.\f$.
    virtual ApproximateKleenean feasible(ApproximateFeasibilityProblem p) const override;


    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal and dual variables.
    virtual ValuePrimalDualData<ApproximateNumber>
    minimise_hotstarted(const ApproximateOptimisationProblem& p,
                        const Vector<Approximation<FLT>>& x0, const Vector<Approximation<FLT>>& y0) const;

    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal and dual variables.
    virtual FeasiblePrimalDualData<ApproximateNumber>
    feasible_hotstarted(const ApproximateFeasibilityProblem& p,
                        const Vector<Approximation<FLT>>& x0, const Vector<Approximation<FLT>>& y0) const;

    //! \brief <p/>
    StepData* initial_step_data(const ApproximateFeasibilityProblem& p) const;

    //! \brief <p/>
    StepData* initial_step_data_hotstarted(const ApproximateFeasibilityProblem& p,
                                           const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y) const;

    Void minimisation_step(const ApproximateOptimisationProblem& p, StepData& d) const;

    Void feasibility_step(const ApproximateFeasibilityProblem& p, StepData& d) const;


  public:
    //! \brief <p/>
    Vector<Approximation<FLT>> compute_dual(const ApproximateBoxType& D, const Vector<Approximation<FLT>>& x, const Approximation<FLT>& mu) const;
    //! \brief <p/>
    Vector<Approximation<FLT>> compute_x(const ApproximateFeasibilityProblem& p) const;
    //! \brief <p/>
    Vector<Approximation<FLT>> compute_y(const ApproximateFeasibilityProblem& p, const Vector<Approximation<FLT>>& x, const Approximation<FLT>& mu) const;
    //! \brief <p/>
    Vector<Approximation<FLT>> compute_w(const ApproximateFeasibilityProblem& p,
                                         const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y, const Approximation<FLT>& mu) const;
    //! \brief <p/>
    Vector<Approximation<FLT>> compute_w(const ApproximateFeasibilityProblem& p,
                                         const Vector<Approximation<FLT>>& x, const Approximation<FLT>& mu) const;
    //! \brief <p/>
    Vector<Approximation<FLT>> compute_z(const ApproximateFeasibilityProblem& p,
                                         const Vector<Approximation<FLT>>& x, const Approximation<FLT>& mu) const;
    //! \brief <p/>
    Approximation<FLT> compute_mu(const ApproximateFeasibilityProblem& p,
                                    const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y) const;
    Void compute_tz(const ApproximateFeasibilityProblem& p,
                    Vector<Approximation<FLT>>& x, Approximation<FLT>& t, Vector<Approximation<FLT>>& z) const;

    //! \brief <p/>
    Approximation<FLT> compute_t(const ApproximateFeasibilityProblem& p,
                                 const Vector<Approximation<FLT>>& x) const;
  protected:
    virtual StepData* _initial_step_data_hotstarted(const ApproximateFeasibilityProblem& p,
                                                    const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y) const = 0;
    //! \brief <p/>
    virtual Void _minimisation_step(const ApproximateOptimisationProblem& p, StepData& d) const = 0;

    //! \brief <p/>
    virtual Void _feasibility_step(const ApproximateFeasibilityProblem& p, StepData& d) const;

};

struct InteriorPointOptimiserBase::StepData {
    using T=Approximation<FLT>;
    T mu;
    virtual ~StepData() = default;
    explicit StepData(T mu_) : mu(mu_) { }
    explicit StepData(DP pr) : StepData(T(pr)) { }
    virtual Vector<T> primal() const = 0;
    virtual Vector<T> dual() const = 0;
    friend OutputStream& operator<<(OutputStream& os, const StepData& d) {
        return d._write(os); }
  private:
    virtual OutputStream& _write(OutputStream& os) const = 0;
};



//! \ingroup OptimisationSubModule
//! An interior-point optimiser using only primal (\f$x\f$) and dual (\f$y\f$) variables.
//! In particular, there are no complementary variables dual to the constraints \f$\unl{d}_i\leq x_i \leq \ovl{d}_i\f$, or slack variables \f$w=g(x)\f$.
//! The initial point \f$x_0\f$ for the optimisation iteration must satisfy \f$x_0\in D^\circ\f$ and \f$g(x_0)\in C^\circ\f$.
class PrimalDualInteriorPointOptimiser
    : public InteriorPointOptimiserBase
{
    using StepDataBase = InteriorPointOptimiserBase::StepData;

    virtual PrimalDualInteriorPointOptimiser* clone() const override;
  public:
    struct StepData;

    StepData* initial_step_data(const ApproximateFeasibilityProblem& p) const;
    StepData* initial_step_data_hotstarted(const ApproximateFeasibilityProblem& p,
                                           const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y) const;
    Void minimisation_step(const ApproximateOptimisationProblem& p,
                           StepData& d) const;

    friend OutputStream& operator<<(OutputStream& os, const PrimalDualInteriorPointOptimiser& opt);
  public:
    Void error_feasibility_step(const ApproximateFeasibilityProblem& p,
                                ApproximateVectorType& x, ApproximateVectorType& y, ApproximateNumberType& t) const;
  private:
    virtual StepDataBase* _initial_step_data_hotstarted(const ApproximateFeasibilityProblem& p,
                                                        const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y) const override;

    virtual Void _minimisation_step(const ApproximateOptimisationProblem& p, StepDataBase& d) const override;
};

//! \relates PrimalDualInteriorPointOptimiser
struct PrimalDualInteriorPointOptimiser::StepData
    : virtual InteriorPointOptimiserBase::StepData, PrimalDualData<Approximation<FLT>>
{
    StepData(Vector<T> x_, Vector<T> y_)
        : InteriorPointOptimiserBase::StepData(T(1,pr)), PrimalDualData<T>(x_,y_) { }
    virtual Vector<Approximation<FLT>> primal() const override { return this->x; }
    virtual Vector<Approximation<FLT>> dual() const override { return this->y; }
    friend OutputStream& operator<<(OutputStream& os, const StepData& d) {
      return os << "(x=" << d.x << ", y=" << d.y << ", mu=" << d.mu << ")"; }
  private:
    virtual OutputStream& _write(OutputStream& os) const override { return os << *this; }
};


//! \ingroup OptimisationSubModule
//! Solver for nonlinear programming problems using the infeasible interior point method with primal (\f$x\f$), dual (\f$y\f$) and complementary (\f$z\f$)  variables.
//! The initial point \f$x_0\f$ for the optimisation iteration must satisfy \f$x_0\in D^\circ\f$ and \f$g(x_0)\in C^\circ\f$.
//! \details Relies on an engine minimising \f$f(x)\f$ for \f$x\in D^\circ\f$ subject to \f$g(x)=w\f$ with \f$w\in C^\circ\f$.
//!
//! To check feasibility, we maximise \f$\mu \sum_{i}\bigl(\log(x_i-\unl{d}_i)+\log(\ovl{d}_i-x_i)\bigr) + \mu \sum_{j}\bigl(\log(w_j-\unl{c}_j)+\log(\ovl{c}_j-w_j)\bigr)\f$ where \f$w=g(x)\f$.
//! Introducing the Lagrange multipliers \f$y_j = \bigl(1/(w_j-\unl{c}_j)-1/(\ovl{c}_j-x_j)\bigr)\f$ and \f$z_i = 1/(x_i-\unl{d}_i)-1/(\ovl{d}_i-x_i)\f$, we find optimality conditions \f$ \sum_{j} y_j \nabla{g_j}(x) + z = 0\f$.
//! These are exactly the central path equations with \f$\mu=1\f$ for the problem \f$\mathrm{minimise} f(x)\equiv 0\f$ subject to \f$x\in D\f$ and \f$g(x)\in C\f$, so we can use an optimisation step to improve feasibility.
class PrimalDualComplementaryInteriorPointOptimiser
    : public InteriorPointOptimiserBase
{
  public:
    using StepDataBase = InteriorPointOptimiserBase::StepData;

    struct StepData;

  public:
    virtual PrimalDualComplementaryInteriorPointOptimiser* clone() const override;

    using OptimiserBase::minimise;
    using OptimiserBase::feasible;

  public:
    //! \brief Construct with default accuracy parameters.
    PrimalDualComplementaryInteriorPointOptimiser();

    //! \brief <p/>
    StepData* initial_step_data(const ApproximateFeasibilityProblem& p) const;

    //! \brief <p/>
    StepData* initial_step_data_hotstarted(const ApproximateFeasibilityProblem& p,
                                           const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y) const;
    //! \brief <p/>
    Void minimisation_step(const ApproximateOptimisationProblem& p,
                           StepData&) const;

    //! \brief <p/>
    friend OutputStream& operator<<(OutputStream& os, PrimalDualComplementaryInteriorPointOptimiser const& opt);
  public:
    //! \brief <p/>
    Void minimisation_step(const ApproximateOptimisationProblem& p,
                           Vector<Approximation<FLT>>& x, Vector<Approximation<FLT>>& y, Vector<Approximation<FLT>>& z,
                           Approximation<FLT>& mu) const;
    //! \brief <p/>
    Void feasibility_step(const ApproximateFeasibilityProblem& p,
                          Vector<Approximation<FLT>>& x, Vector<Approximation<FLT>>& y, Vector<Approximation<FLT>>& z) const;

    PrimalDualComplementaryData<Approximation<FLT>> minimisation_update(
        const ApproximateOptimisationProblem& p,
        ApproximateVectorType& x, ApproximateVectorType& y, ApproximateVectorType& z, Approximation<FLT>& mu) const;
  public:
    Void error_feasibility_step(const ApproximateFeasibilityProblem& p,
                                ApproximateVectorType& x, ApproximateVectorType& yl, ApproximateVectorType& yu, ApproximateVectorType& z,
                                ApproximateNumberType& t) const;
  private:
    virtual StepDataBase* _initial_step_data_hotstarted(const ApproximateFeasibilityProblem& p,
                                                        const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y) const override;
    virtual Void _minimisation_step(const ApproximateOptimisationProblem& p,
                                    StepDataBase&) const override;
};

//! \relates PrimalDualComplementaryInteriorPointOptimiser
struct PrimalDualComplementaryInteriorPointOptimiser::StepData
    : public virtual InteriorPointOptimiserBase::StepData, public PrimalDualComplementaryData<Approximation<FLT>>
{
    StepData(SizeType m, SizeType n, DP pr)
        : InteriorPointOptimiserBase::StepData(pr), PrimalDualComplementaryData<T>(m,n,pr) { }
    StepData(Vector<T> x_, Vector<T> y_, Vector<T> z_, T mu_)
        : InteriorPointOptimiserBase::StepData(mu_), PrimalDualComplementaryData<T>(x_,y_,z_) { }
    virtual Vector<T> primal() const override { return this->x; }
    virtual Vector<T> dual() const override { return this->y; }
    friend OutputStream& operator<<(OutputStream& os, StepData const& d) {
        return os << "(x=" << d.x << ", y=" << d.y << ", z=" << d.z << ", mu=" << d.mu << ")"; }
  private:
    virtual OutputStream& _write(OutputStream& os) const override { return os << *this; }
};




//! \ingroup OptimisationSubModule
//! Solver for nonlinear programming problems using the infeasible interior point method with primal (\f$x\f$), dual (\f$y\f$) and complementary (\f$z\f$)  variables.
//! \details Relies on an engine minimising \f$f(x)\f$ for \f$x\in D^\circ\f$ subject to \f$g(x)=w\f$ with \f$w\in C\f$.
//!
//! To check feasibility, we maximise \f$\mu \sum_{i}\bigl(\log(x_i-\unl{d}_i)+\log(\ovl{d}_i-x_i)\bigr) + \mu \sum_{j}\bigl(\log(w_j-\unl{c}_j)+\log(\ovl{c}_j-w_j)\bigr)\f$ under the constraints \f$g(x)-w=0\f$.
//! Introducing the Lagrange multipliers \f$y_j = \bigl(1/(w_j-\unl{c}_j)-1/(\ovl{c}_j-x_j)\bigr)\f$ and \f$z_i = 1/(x_i-\unl{d}_i)-1/(\ovl{d}_i-x_i)\f$, we find optimality conditions \f$ \sum_{j} y_j \nabla{g_j}(x) + z = 0\f$, \f$g(x)-w=0\f$.
//! These are exactly the central path equations with \f$\mu=1\f$ for the problem \f$\mathrm{minimise} f(x)\equiv 0\f$ subject to \f$x\in D\f$ and \f$g(x)\in C\f$, so we can use an optimisation step to improve feasibility.
class SlackPrimalDualComplementaryInteriorPointOptimiser
    : public InteriorPointOptimiserBase
{
  public:
    using StepDataBase = InteriorPointOptimiserBase::StepData;

    struct StepData;

  public:
    virtual SlackPrimalDualComplementaryInteriorPointOptimiser* clone() const override;

    using OptimiserBase::minimise;
    using OptimiserBase::feasible;

  public:
    //! \brief Construct with default accuracy parameters.
    SlackPrimalDualComplementaryInteriorPointOptimiser();

    //! \brief <p/>
    StepData* initial_step_data(const ApproximateFeasibilityProblem& p) const;
    //! \brief <p/>
    StepData* initial_step_data_hotstarted(const ApproximateFeasibilityProblem& p,
                                           const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y) const;
    //! \brief <p/>
    Approximation<FLT> minimisation_step(const ApproximateOptimisationProblem& p,
                                         StepData&) const;

    //! \brief <p/>
    friend OutputStream& operator<<(OutputStream& os, SlackPrimalDualComplementaryInteriorPointOptimiser const& opt);
  public:
    //! \brief <p/>
    Approximation<FLT> minimisation_step(const ApproximateOptimisationProblem& p,
                                         Vector<Approximation<FLT>>& w, Vector<Approximation<FLT>>& x,
                                         Vector<Approximation<FLT>>& y, Vector<Approximation<FLT>>& z,
                                         Approximation<FLT>& mu) const;
    //! \brief <p/>
    Void feasibility_step(const ApproximateFeasibilityProblem& p,
                          Vector<Approximation<FLT>>& w, Vector<Approximation<FLT>>& x, Vector<Approximation<FLT>>& y, Vector<Approximation<FLT>>& z) const;

    Approximation<FLT> minimisation_step_size(
        const ApproximateOptimisationProblem& p,
        const ApproximateVectorType& w, const ApproximateVectorType& x, const ApproximateVectorType& y, const ApproximateVectorType& z,
        ApproximateVectorType& dw, ApproximateVectorType& dx, ApproximateVectorType& dy, ApproximateVectorType& dz,
        const Approximation<FLT>& mu) const;

    SlackPrimalDualComplementaryData<Approximation<FLT>> minimisation_update(
        const ApproximateOptimisationProblem& p,
        ApproximateVectorType& w, ApproximateVectorType& x, ApproximateVectorType& y, ApproximateVectorType& z, Approximation<FLT>& mu) const;

  private:
    virtual StepDataBase* _initial_step_data_hotstarted(const ApproximateFeasibilityProblem& p,
                                                        const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y) const override;
    virtual Void _minimisation_step(const ApproximateOptimisationProblem& p,
                                    StepDataBase&) const override;
};

//! \relates SlackPrimalDualComplementaryInteriorPointOptimiser
struct SlackPrimalDualComplementaryInteriorPointOptimiser::StepData
    : public virtual InteriorPointOptimiserBase::StepData, public SlackPrimalDualComplementaryData<Approximation<FLT>>
{
    StepData(SizeType m, SizeType n, DP pr)
        : InteriorPointOptimiserBase::StepData(pr), SlackPrimalDualComplementaryData<T>(m,n,pr) { }
    StepData(Vector<T> w_, Vector<T> x_, Vector<T> y_, Vector<T> z_, T mu_)
        : InteriorPointOptimiserBase::StepData(mu_), SlackPrimalDualComplementaryData<T>(w_,x_,y_,z_) { }
    virtual Vector<T> primal() const override { return this->x; }
    virtual Vector<T> dual() const override { return this->y; }
    friend OutputStream& operator<<(OutputStream& os, StepData const& d) {
        return os << "(w=" << d.w << ", x=" << d.x << ", y=" << d.y << ", z=" << d.z << ", mu=" << d.mu << ")"; }
  private:
    virtual OutputStream& _write(OutputStream& os) const override { return os << *this; }
};




//! \ingroup OptimisationSubModule
//! Solver for nonlinear programming problems using the infeasible interior point method with primal (\f$x\f$), dual (\f$\unl{y},\ovl{y}\f$), slack (\f$w\f$) and complementary (\f$\unl{z},\ovl{z}\f$) variables.
//! The dual and complementary variables are split into lower and upper versions, corresponding to lower and upper constraint bounds.
//! \details The complementary variables satisfy the central path equations \f$(x-\unl{d})\unl{z} = \mu\f$ and \f$(\ovl{d}-x)\ovl{z} = \mu\f$ with \f$\unl{z},\ovl{z}>0\f$, and can be combined into \f$z = \ovl{z}-\unl{z}\f$ satisfying \f$z = \mu\bigl(1/(\ovl{d}-x)-1/(x-\unl{d})\bigr)\f$, or equivalently \f$(x-\unl{d})(\ovl{d}-x)z + (\unl{d}+\ovl{d}-2x) \mu = 0\f$. Similarly, \f$(w-\unl{c})\unl{y} = \mu\f$ and \f$(\ovl{c}-w)\ovl{y} = \mu\f$ where \f$w=g(x)\f$.
class SlackPrimalSplitDualComplementaryInteriorPointOptimiser
    : public InteriorPointOptimiserBase
{
    using StepDataBase = InteriorPointOptimiserBase::StepData;

    virtual SlackPrimalSplitDualComplementaryInteriorPointOptimiser* clone() const override;
  public:
    struct StepData;

    //! \brief <p/>
    StepData* initial_step_data(const ApproximateFeasibilityProblem& p) const;
    //! \brief <p/>
    StepData* initial_step_data_hotstarted(const ApproximateFeasibilityProblem& p,
                                           const Vector<Approximation<FLT>>& x0, const Vector<Approximation<FLT>>& y0) const;
    //! \brief Performs one step of the interior point algorithm on \a d in-place. Returns the step-size parameter \f$\alpha\f$ used.
    Approximation<FLT> minimisation_step(const ApproximateOptimisationProblem& p,
                                         StepData& d) const;
    //! \brief <p/>
    friend OutputStream& operator<<(OutputStream& os, SlackPrimalSplitDualComplementaryInteriorPointOptimiser const& opt);
  private:
    virtual StepDataBase* _initial_step_data_hotstarted(const ApproximateFeasibilityProblem& p,
                                                        const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y) const override;

    virtual Void _minimisation_step(const ApproximateOptimisationProblem& p, StepDataBase& d) const override;
};

//! \relates SlackPrimalSplitDualComplementaryInteriorPointOptimiser
struct SlackPrimalSplitDualComplementaryInteriorPointOptimiser::StepData
    : virtual InteriorPointOptimiserBase::StepData, SlackPrimalSplitDualComplementaryData<Approximation<FLT>>
{
    StepData(SizeType m, SizeType n, PR pr)
        : InteriorPointOptimiserBase::StepData(T(pr)), SlackPrimalSplitDualComplementaryData<T>(m,n,pr) { }
    StepData(Vector<T> w_, Vector<T> x_, Vector<T> yl_, Vector<T> yu_, Vector<T> zl_, Vector<T> zu_, T mu_)
        : InteriorPointOptimiserBase::StepData(mu_), SlackPrimalSplitDualComplementaryData<T>(w_,x_,yl_,yu_,zl_,zu_) { }
    virtual Vector<T> primal() const override { return this->x; }
    virtual Vector<T> dual() const override { return this->yu-this->yl; }
    friend OutputStream& operator<<(OutputStream& os, const StepData& d) {
      return os << "(w=" << d.w << ", x=" << d.x << ", yl=" << d.yl << ", yu=" << d.yu << ", zl=" << d.zl << ", zu=" << d.zu << ", mu=" << d.mu << ")"; }
  private:
    virtual OutputStream& _write(OutputStream& os) const override { return os << *this; }
};


using InteriorPointOptimiser = PrimalDualComplementaryInteriorPointOptimiser;
using InfeasibleInteriorPointOptimiser = SlackPrimalDualComplementaryInteriorPointOptimiser;
using SplitInfeasibleInteriorPointOptimiser = SlackPrimalSplitDualComplementaryInteriorPointOptimiser;


class KarushKuhnTuckerOptimiser
    : public ValidatedOptimiserBase
{
  public:
    KarushKuhnTuckerOptimiser* clone() const override;

    //! \brief Compute a \em local optimum of nonlinear programming problem \f$\max f(x) \text{ such that } x\in D \text{ and } g(x)\in C .\f$.
    //! \pre The domain \f$D\f$ is bounded and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A box \f$X\f$ which definitely contains a feasible point, and contains a local optimum.
    virtual ValidatedVector minimise(ValidatedOptimisationProblem p) const override;
    //! \brief Tests if the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ValidatedFeasibilityProblem) const override;

    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal and dual variables.
    ValuePrimalDualData<ValidatedNumber>
    minimise_hotstarted(const ValidatedOptimisationProblem& p,
                        const Vector<Approximation<FLT>>& x0, const Vector<Approximation<FLT>>& y0) const;

    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal and dual variables.
    FeasiblePrimalDualData<ValidatedNumber>
    feasible_hotstarted(const ValidatedFeasibilityProblem& p,
                        const Vector<Approximation<FLT>>& x0, const Vector<Approximation<FLT>>& y0) const;

    //! \brief Checks whether the point \f$x\f$ is feasible, or the multipliers \f$y\f$ provide a certificate of infeasibility based around \f$x\f$.
    FeasiblePrimalDualData<ValidatedNumber>
    check_feasibility(const ValidatedFeasibilityProblem& p,
                      const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y) const;

    //! \brief Tests whether the validated point \a x contains a locally optimal point, using Lagrange multipliers \a y.
    FeasiblePrimalDualData<ValidatedNumber>
    check_minimality(const ValidatedOptimisationProblem& p,
                     const Vector<Approximation<FLT>>& w, const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y, const Vector<Approximation<FLT>>& z) const;

    //! \brief Tests whether the validated point \a x contains a locally optimal point, using Lagrange multipliers \a y.
    FeasiblePrimalDualData<ValidatedNumber>
    check_minimality(const ValidatedOptimisationProblem& p,
                     const Vector<Bounds<FLT>>& w, const Vector<Bounds<FLT>>& x, const Vector<Bounds<FLT>>& y, const Vector<Bounds<FLT>>& z) const;

    ValuePrimalDualData<ValidatedNumber>
    nonsplitting_minimise_hotstarted(const ValidatedOptimisationProblem& p,
                                     const Vector<Approximation<FLT>>& x0, const Vector<Approximation<FLT>>& y0) const;

    FeasiblePrimalDualData<ValidatedNumber>
    nonsplitting_feasible_hotstarted(const ValidatedFeasibilityProblem& p,
                                     const Vector<Approximation<FLT>>& x0, const Vector<Approximation<FLT>>& y0) const;


    Variant<OptimalityCertificate,InfeasibilityCertificate>
    splitting_minimise_hotstarted(const ValidatedOptimisationProblem& p,
                                  UpperBoxType B, Vector<Approximation<FLT>> x, Vector<Approximation<FLT>> y) const;

    Variant<FeasibilityCertificate,InfeasibilityCertificate>
    splitting_feasible_hotstarted(const ValidatedFeasibilityProblem& p,
                                  UpperBoxType B, Vector<Approximation<FLT>> x, Vector<Approximation<FLT>> y) const;
};

class InfeasibleKarushKuhnTuckerOptimiser
    : public ValidatedOptimiserBase
{
  public:
    using ValidatedVector = ValidatedOptimiserInterface::ValidatedVector;

    virtual InfeasibleKarushKuhnTuckerOptimiser* clone() const override;

    using ValidatedOptimiserBase::minimise;
    using ValidatedOptimiserBase::feasible;

    //! \brief Compute a \em local optimum of nonlinear programming problem \f$\max f(x) \text{ such that } x\in D, g(x)\in C \text{ and } h(x)=0.\f$.
    //! \pre The domain \f$D\f$ is bounded and has nonempty interior, and the codomain \f$C\f$ is nonempty.
    //! \return A box \f$X\f$ which definitely contains a feasible point, and contains a local optimum.
    virtual ValidatedVector minimise(ValidatedOptimisationProblem p) const override;
    //! \brief Tests if the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
    virtual ValidatedKleenean feasible(ValidatedFeasibilityProblem p) const override;

    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the primal and dual variables.
    ValuePrimalDualData<ValidatedNumber>
    minimise_hotstarted(const ValidatedOptimisationProblem& p,
                        const Vector<Approximation<FLT>>& x0, const Vector<Approximation<FLT>>& y0) const;

    //! \brief Test if the constraints \f$g(x)\in C\f$ are solvable for \f$x\in D\f$ using a nonlinear feasibility test,
    //! hotstarting the method with the overall primal and dual variables.
    FeasiblePrimalDualData<ValidatedNumber>
    feasible_hotstarted(const ValidatedFeasibilityProblem& p,
        const SlackPrimalDualData<Approximation<FLT>>& wxy0) const;

    Variant<OptimalityCertificate,InfeasibilityCertificate>
    splitting_minimise_hotstarted(const ValidatedOptimisationProblem& p,
                                  UpperBoxType B, Vector<Approximation<FLT>> x, Vector<Approximation<FLT>> y) const;

    //! \brief Tests whether the validated point \a x contains a locally optimal point, using Lagrange multipliers \a y.
    Tuple<ValidatedKleenean,ValidatedVector,ValidatedVector>
    check_minimality(const ValidatedOptimisationProblem& p,
                     const Vector<Approximation<FLT>>& w, const Vector<Approximation<FLT>>& x, const Vector<Approximation<FLT>>& y, const Vector<Approximation<FLT>>& z) const;

    //! \brief Tests whether the validated point \a x contains a locally optimal point, using Lagrange multipliers \a y.
    Tuple<ValidatedKleenean,ValidatedVector,ValidatedVector>
    check_minimality(const ValidatedOptimisationProblem& p,
                     const Vector<Bounds<FLT>>& w, const Vector<Bounds<FLT>>& x, const Vector<Bounds<FLT>>& y, const Vector<Bounds<FLT>>& z) const;
  private:
    Variant<OptimalityCertificate,InfeasibilityCertificate>
    _splitting_minimise_subdivide(const ValidatedOptimisationProblem& p,
                                  UpperBoxType B, Vector<Approximation<FLT>> x, Vector<Approximation<FLT>> y) const;

    auto minimise_old(ValidatedOptimisationProblem p) const -> ValidatedVector;
};


//! \ingroup OptimisationSubModule
//! \brief An optimiser using rigorous interval-arithmetic methods based on the John conditions for feasibility problems.
//! <br>Currently only provides rigorous solvers for feasibility problems; optimisation problems are not supported.
class IntervalOptimiser
    : public KarushKuhnTuckerOptimiser
{
  private: public:
    virtual IntervalOptimiser* clone() const override;
    //! \brief <p/>
    virtual ValidatedKleenean feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const;
    //! \brief <p/>
    Void feasibility_step(const FloatDPVector& xl, const FloatDPVector& xu, const ValidatedVectorMultivariateFunction& h,
                          FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& zl, FloatDPBoundsVector zu, FloatDPBounds& mu) const;
};


// \ingroup OptimisationSubModule
// \brief An optimiser based on the PrimalDualComplementaryInteriorPointOptimiser. \deprecated
class ApproximateOptimiser
    : public PrimalDualComplementaryInteriorPointOptimiser
{
  private: public:
    virtual ApproximateOptimiser* clone() const override;
    // \brief <p/>
    virtual ValidatedKleenean feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const;
    // \brief <p/>
    Void feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& h,
                          Vector<Approximation<FLT>>& X, Vector<Approximation<FLT>>& Lambda) const;
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
    //! \brief Tests if the nonlinear programming problem \f$x\in D \text{ and } g(x)\in C\f$ is feasible.
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
