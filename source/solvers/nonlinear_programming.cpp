/***************************************************************************
 *            solvers/nonlinear_programming.cpp
 *
 *  Copyright  2010-20  Pieter Collins
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

// See Hande Y. Benson, David F. Shanno, And Robert J. Vanderbei,
// "Interior-point methods for nonconvex nonlinear programming: Jamming and comparative numerical testing"
// For some of the terminology used

#include "function/functional.hpp"
#include "config.hpp"

#include <limits>

#include "utility/macros.hpp"
#include "conclog/logging.hpp"
#include "utility/tuple.hpp"
#include "utility/tribool.hpp"
#include "numeric/numeric.hpp"
#include "algebra/linear_algebra.decl.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/diagonal_matrix.hpp"
#include "algebra/symmetric_matrix.hpp"
#include "algebra/differential.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/function_mixin.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.hpp"
#include "function/procedure.hpp"

#include "solvers/nonlinear_programming.hpp"
#include "solvers/solver.hpp"
#include "algebra/multi_index-noaliasing.hpp"
#include "solvers/constraint_solver.hpp"

#include "algebra/expansion.inl.hpp"
#include "solvers/nonlinear_programming.tpl.hpp"

using namespace ConcLog;

// CONCLOG_NUM_ARGS extracts the number of arguments of a variadic argument pack
#define _CONCLOG_NUM_ARGS(X8, X7, X6, X5, X4, X3, X2, X1, N, ...) N
#define CONCLOG_NUM_ARGS(...) _CONCLOG_NUM_ARGS(__VA_ARGS__, 8, 7, 6, 5, 4, 3, 2, 1)

#define CONCLOG_EXPAND(x)  x
#define CONCLOG_FIRSTARG(x, ...)  (x)
#define CONCLOG_RESTARGS(x, ...)  (__VA_ARGS__)
#define CONCLOG_FOREACH(m, xs)  CONCLOG_FOREACH_(CONCLOG_NUM_ARGS xs, m, xs)
#define CONCLOG_FOREACH_(n, m, xs)  CONCLOG_FOREACH__(n, m, xs)
#define CONCLOG_FOREACH__(n, m, xs)  CONCLOG_FOREACH_##n(m, xs)
#define CONCLOG_FOREACH_1(m, xs)  m xs
#define CONCLOG_FOREACH_2(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_1(m, CONCLOG_RESTARGS xs)
#define CONCLOG_FOREACH_3(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_2(m, CONCLOG_RESTARGS xs)
#define CONCLOG_FOREACH_4(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_3(m, CONCLOG_RESTARGS xs)
#define CONCLOG_FOREACH_5(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_4(m, CONCLOG_RESTARGS xs)
#define CONCLOG_FOREACH_6(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_5(m, CONCLOG_RESTARGS xs)
#define CONCLOG_FOREACH_7(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_6(m, CONCLOG_RESTARGS xs)
#define CONCLOG_FOREACH_8(m, xs)  CONCLOG_EXPAND(m CONCLOG_FIRSTARG xs) << ", " CONCLOG_FOREACH_7(m, CONCLOG_RESTARGS xs)

#define CONCLOG_WRITE_VAR(x)    << #x << "=" << (x)
#define CONCLOG_WRITE_VARS(...) CONCLOG_FOREACH(CONCLOG_WRITE_VAR, (__VA_ARGS__))

#define CONCLOG_INITIALISE_AT(level) if (!Logger::instance().is_muted_at(level)) { std::ostringstream logger_stream; logger_stream << std::boolalpha
#define CONCLOG_FINALISE_AT(level) ; Logger::instance().println(level,logger_stream.str()); }

#define CONCLOG_PRINTLN_VARS_AT(level,...) { CONCLOG_INITIALISE_AT(level) CONCLOG_WRITE_VARS(__VA_ARGS__) CONCLOG_FINALISE_AT(level) }
#define CONCLOG_PRINTLN_VARS(...) { CONCLOG_PRINTLN_VARS_AT(0,__VA_ARGS__) }


namespace Ariadne {

namespace {

template<class T> using UniquePointer = std::unique_ptr<T>;

template<class FLT> using Exact = FLT;

template<class FLT> requires ARawFloat<FLT> inline Scalar<ExactNumber> to_generic(Scalar<Exact<FLT>> x) { return Scalar<ExactNumber>(x); }
template<class FLT> inline Scalar<ValidatedNumber> to_generic(Scalar<Bounds<FLT>> x) { return Scalar<ValidatedNumber>(x); }
template<class FLT> inline Scalar<ApproximateNumber> to_generic(Scalar<Approximation<FLT>> x) { return Scalar<ApproximateNumber>(x); }

template<class FLT> inline Vector<ExactNumber> to_generic(Vector<Exact<FLT>> v) { return Vector<ExactNumber>(v); }
template<class FLT> inline Vector<ValidatedNumber> to_generic(Vector<Bounds<FLT>> v) { return Vector<ValidatedNumber>(v); }
template<class FLT> inline Vector<ApproximateNumber> to_generic(Vector<Approximation<FLT>> v) { return Vector<ApproximateNumber>(v); }

inline Sweeper<FloatDP> default_sweeper() { return Sweeper<FloatDP>(); }

} // namespace

typedef DiagonalMatrix<FloatDPBounds> FloatDPBoundsDiagonalMatrix;

typedef Vector<FloatDPApproximation> ApproximateVectorType;
typedef VectorRange<ApproximateVectorType> ApproximateVectorTypeRange;
typedef Matrix<FloatDPApproximation> FloatDPApproximationMatrix;
typedef DiagonalMatrix<FloatDPApproximation> FloatDPApproximationDiagonalMatrix;
typedef Differential<FloatDPApproximation> FloatDPApproximationDifferential;

typedef Vector<FloatDPBounds> FloatDPBoundsVector;
typedef Matrix<FloatDPBounds> FloatDPBoundsMatrix;
typedef Differential<FloatDPBounds> FloatDPBoundsDifferential;
typedef Vector<FloatDP> ExactFloatDPVectorType;

typedef Vector<UpperIntervalType> UpperIntervalVectorType;
typedef Matrix<UpperIntervalType> UpperIntervalMatrixType;

typedef FloatDPApproximation ApproximateNumericType;

template<class T> struct Pretty;

template<> struct Pretty<Matrix<FloatDPApproximation>> {
    static const uint prec=4; static const uint wdth=prec+6;
    Matrix<FloatDPApproximation> const& ref;
    friend OutputStream& operator<<(OutputStream& os, Pretty<Matrix<FloatDPApproximation>> const& pA) {
        auto& A=pA.ref;
        os << "\n";
        for(SizeType i=0; i!=A.row_size(); ++i) {
            os << "[";
            for(SizeType j=0; j!=A.column_size(); ++j) {
                if (j!=0) { os << ","; } os << std::setprecision(prec) << std::setw(wdth) << A[i][j].get_d();
            }
            os << "]\n";
        }
        return os;
    }
};
template<class T> Pretty<T> pretty(T const& t) { return Pretty<T>{t};}



inline FloatDPApproximation operator*(ApproximateDouble x1, FloatDPApproximation x2) {
    return FloatDPApproximation(x1,x2.precision()) * x2;
}
inline FloatDPApproximation operator*(FloatDPApproximation x1, ApproximateDouble x2) {
    return x1 * FloatDPApproximation(x2,x1.precision());
}
inline ApproximateKleenean operator<(FloatDPApproximation x1, ApproximateDouble x2) {
    return x1 < FloatDPApproximation(x2,x1.precision());
}


template OutputStream& operator<<(OutputStream&, FeasibilityProblem<ApproximateTag> const&);
template OutputStream& operator<<(OutputStream&, FeasibilityProblem<ValidatedTag> const&);

template<class P> OutputStream& operator<<(OutputStream& os, OptimisationProblem<P> const& p) {
    return os << "OptimisationProblem( f=" << p.f << ", D=" << p.D << ", g=" << p.g << ", C=" << p.C << " )";
}
template OutputStream& operator<<(OutputStream&, OptimisationProblem<ApproximateTag> const&);
template OutputStream& operator<<(OutputStream&, OptimisationProblem<ValidatedTag> const&);


inline Box<Interval<FloatDP>> cast_exact_widen(Box<Interval<FloatDP>> const& bx, FloatDP e) {
    Box<Interval<FloatDP>> r(bx);
    for(SizeType i=0; i!=bx.size(); ++i) {
        r[i]=Interval<FloatDP>(sub(down,bx[i].lower_bound(),e),add(up,bx[i].upper_bound(),e));
    }
    return r;
}



//------- FeasibilityProblem -----------------------------------//

template<class P> FeasibilityProblem<P>::FeasibilityProblem(BoxType<P> D_, VectorMultivariateFunction<P> g_, BoxType<P> C_)
    : D(D_), g(g_), C(C_)
{
    ARIADNE_PRECONDITION(g.argument_size()==D.dimension());
    ARIADNE_PRECONDITION(g.result_size()==C.dimension());
}

template<class P> FeasibilityProblem<P>::FeasibilityProblem(BoxType<P> D_, const List<ConstraintType<P>>& cl)
    : D(D_), g(_function(D.dimension(),cl)), C(_bounds(cl)) { }

template<class P> auto FeasibilityProblem<P>::_function(SizeType as, const List<ConstraintType<P>>& cl) -> VectorMultivariateFunction<P> {
    if (cl.size()==0) {
        return VectorMultivariateFunction<P>(0u, ScalarMultivariateFunction<P>(as));
    } else {
        List<ScalarMultivariateFunction<P>> lst;
        for (SizeType i=0; i!=cl.size(); ++i) {
            lst.append(cl[i].function());
        }
        return VectorMultivariateFunction<P>(lst);
    }
}

template<class P> auto FeasibilityProblem<P>::_bounds(const List<ConstraintType<P>>& cl) -> BoxType<P> {
    return BoxType<P>(cl.size(),[&cl](SizeType i){return cl[i].cast_exact_bounds();});
}

template struct FeasibilityProblem<ValidatedTag>;
template struct FeasibilityProblem<ApproximateTag>;


//------- FeasibilityChecker -----------------------------------//

FeasibilityChecker* FeasibilityChecker::
clone() const
{
    return new FeasibilityChecker(*this);
}


ApproximateKleenean FeasibilityChecker::
almost_feasible_point(ValidatedFeasibilityProblem p, ApproximateVector ax, ApproximateNumber e) const
{
    Vector<Approximation<FLT>> axf(ax,pr);
    Approximation<FLT> ef(e,pr);

    if (!probably(contains(p.D,axf))) { return false; }
    Vector<Approximation<FLT>> gxf=p.g(axf);
    return contains(cast_exact_widen(p.C,cast_exact(ef)),gxf);
}


ValidatedKleenean FeasibilityChecker::
is_feasible_point(ValidatedFeasibilityProblem p, ExactVector x) const
{
    Vector<Bounds<FLT>> xf(x,pr);

    if (not possibly(contains(p.D,xf))) { return false; }
    Vector<Bounds<FLT>> gxf=p.g(xf);
    return contains(p.C,gxf);
}

ValidatedKleenean FeasibilityChecker::
contains_feasible_point(ValidatedFeasibilityProblem p, UpperBoxType X) const
{
    if (this->validate_feasibility(p,cast_singleton(X))) {
        return true;
    } else {
        p.D=intersection(p.D,cast_exact_box(X));
        if (this->validate_infeasibility(p)) {
            return false;
        } else {
            return indeterminate;
        }
    }
}


ValidatedKleenean FeasibilityChecker::
check_feasibility(ValidatedFeasibilityProblem p,
                  ValidatedVector x, ExactVector y) const
{
    return this->check_feasibility(p, Vector<Bounds<FLT>>(x,dp), Vector<Bounds<FLT>>(y,pr));
}

ValidatedKleenean FeasibilityChecker::
check_feasibility(ValidatedFeasibilityProblem p,
                  Vector<Bounds<FLT>> x, Vector<Bounds<FLT>> y) const
{
    if (this->validate_feasibility(p,x)) {
        return true;
    } else if (this->validate_infeasibility(p,x,y)) {
        return false;
    } else {
        return indeterminate;
    }
}


Bool FeasibilityChecker::
validate_feasibility(ValidatedVectorMultivariateFunction h,
                     ValidatedVector x) const
{
    return this->validate_feasibility(h, Vector<Bounds<FLT>>(x,pr));
}


Bool FeasibilityChecker::
validate_feasibility(ValidatedVectorMultivariateFunction h,
                     Vector<Bounds<FLT>> X) const
{
    // Let h:R^n->R^m with m<n, and X define a box in R^n.
    // Attempt to solve h(c+R*B*w)=0 for w, where c=mid(X), R is a diagonal scaling matrix, and B is an n*m matrix.
    // Let z lie in the unit box Z=[-1:+1]^n, and take x=R*z+c where R=diag(rad(X)) is a scaling matrix.
    // Then the function z->h(c+R*z) has Jacobian Dh(X)*R = [A]*R over Z.
    // Then for an interval Newton step, we have [A]*R*B dw = -h(c+B*w)
    // Take A to be an approximation to [A], such as A=mid([A]) or A=Dh(C), and set B=(A*R)^T=R*A^T
    // Then we solve h(c+R*R*AT*w)=0 using Newton's method centred at w=0, yielding W' = -([A]*R*R*A^T)\h(c)

    Vector<Approximation<FLT>> ca = midpoint(X);
    Vector<FLT> c = cast_exact(ca);
    Matrix<FLT> AT = transpose(cast_exact(h.jacobian(ca)));
    CONCLOG_PRINTLN("A="<<transpose(AT));


    Vector<Bounds<FLT>> W=h(X);
    Matrix<Bounds<FLT>> A = h.jacobian(X);
    DiagonalMatrix<FloatDP> R(Array<FloatDP>(X.size(), [&X](SizeType i){return cast_exact(X[i].error());}));


    // Could take w to be centre of set W
    // Vector<FloatDP> w=cast_exact(W);
    // new_W = w - gs_solve(A*AT,h(x+AT*w));
    // Easier to take w = 0, which should be an element of W
    Vector<Bounds<FLT>> new_W = - gs_solve(A*(R*R*AT),h(c));
    Vector<Bounds<FLT>> new_X = c + (R * R) * (AT * new_W);

    if (refines(new_X,X)) {
        return true;
    } else {
        if (refines(new_W,W)) {
            ARIADNE_WARN("Did not verify feasible point in "<<X<<", but one may exist.");
        }
        return false;
    }
}

Bool FeasibilityChecker::
validate_feasibility(ValidatedFeasibilityProblem p,
                     ValidatedVector x) const
{
    return this->validate_feasibility(p,Vector<Bounds<FLT>>(x,pr));
}

Bool FeasibilityChecker::
validate_feasibility(ValidatedFeasibilityProblem p,
                     Vector<Bounds<FLT>> x) const
{
    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C);
    CONCLOG_PRINTLN("x="<<x);

    ARIADNE_PRECONDITION(x.size()==D.size());


    for(SizeType i=0; i!=D.size(); ++i) {
        CONCLOG_PRINTLN_AT(1,"x["<<i<<"]="<<x[i]<<", D["<<i<<"]="<<D[i]);
        if (!possibly(contains(D[i],x[i]))) {
            return false;
        } else if (definitely(contains(D[i],x[i]))) {
        } else {
            x[i]=cast_singleton(intersection(UpperIntervalType(x[i]),D[i]));
        }
    }

    Vector<Bounds<FLT>> w=g(x);
    CONCLOG_PRINTLN_AT(1,"w=g(x)="<<w);

    List<SizeType> equalities;
    for(SizeType j=0; j!=C.size(); ++j) {
        CONCLOG_PRINTLN_AT(1,"w["<<j<<"]="<<w[j]<<", C["<<j<<"]="<<C[j]);
        if (!possibly(contains(C[j],w[j]))) {
            return false;
        } else if (definitely(contains(C[j],w[j]))) {
        } else {
            // NOTE: It is safe to try solving a constraint as an equality
            if ( true || decide(C[j].lower_bound()==C[j].upper_bound()) ) {
                equalities.append(j);
            }
        }
    }
    if(equalities.empty()) {
        return true;
    }

    ValidatedVectorMultivariateFunction h(equalities.size(),g.domain());
    for(SizeType i=0; i!=equalities.size(); ++i) {
        SizeType j=equalities[i];
        h[i] = g[j]-static_cast<ExactNumber>(intersection(UpperIntervalType(w[j]),C[j]).midpoint());
    }
    CONCLOG_PRINTLN_AT(1,"h="<<h);

    return this->validate_feasibility(h,x);
}


Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p,
                       UpperBoxType B, ExactVector y) const
{
    return validate_infeasibility(p,B,Vector<Bounds<FLT>>(y,pr));
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p,
                       UpperBoxType B, Vector<Exact<FLT>> y) const
{
    return validate_infeasibility(p,B,Vector<Bounds<FLT>>(y));
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p,
                       UpperBoxType B, Vector<Bounds<FLT>> y) const
{
    auto& D=p.D; auto& g=p.g; auto& C=p.C;
    ExactBoxType DB = intersection(D,cast_exact_box(B));
    Vector<Approximation<FLT>> xa = midpoint(B);
    return this->validate_infeasibility(ValidatedFeasibilityProblem(DB,g,C),xa,y);
}


Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p,
                       ApproximateVector xa, ExactVector y) const
{
    return validate_infeasibility(p, Vector<Approximation<FLT>>(xa,dp), Vector<Bounds<FLT>>(y,pr));
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p,
                       Vector<Approximation<FLT>> xa, Vector<Exact<FLT>> y) const
{
    return validate_infeasibility(p, xa, Vector<Bounds<FLT>>(y));
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p,
                       Vector<Approximation<FLT>> xa, Vector<Bounds<FLT>> y) const
{
    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C<<", xa="<<xa<<", y="<<y);

    ARIADNE_PRECONDITION(xa.size()==D.size());
    ARIADNE_PRECONDITION(y.size()==C.size());

    if(y.size()==0) { return D.is_empty(); }

    UpperIntervalType yC = dot(y,C);

    // Compute Taylor estimate of y g(X)
    ValidatedVectorMultivariateTaylorFunctionModelDP tg(D,g,default_sweeper());
    ValidatedScalarMultivariateTaylorFunctionModelDP tyg(D,default_sweeper());
        for(SizeType j=0; j!=y.size(); ++j) { tyg += y[j]*tg[j];
    }
    UpperIntervalType ygD = apply(tyg.function(),D);
    // UpperIntervalType ygD = dot(y,apply(g,D));

    if(definitely(disjoint(yC,ygD))) {
        CONCLOG_PRINTLN("Infeasible");
        return true;
    }

    UpperIntervalMatrixType dgD = jacobian_range(g,cast_vector(D));
    UpperIntervalVectorType ydgD = transpose(dgD)*y;

    Vector<Bounds<FLT>> x = cast_exact(xa);
    Bounds<FLT> ygx = tyg(x);
    // ValidatedNumericType ygx = dot(y,g(x));
    UpperIntervalType ygDx = UpperIntervalType(ygx);
    for(SizeType i=0; i!=x.size(); ++i) {
        ygDx += ydgD[i] * (D[i]-x[i]);
    }

    CONCLOG_PRINTLN("yC="<<yC<<", ygD="<<ygD<<", ygx="<<ygx<<", ydgD="<<ydgD<<", ygDx="<<ygDx);

    if(definitely(disjoint(yC,intersection(ygD,ygDx)))) {
        CONCLOG_PRINTLN("Infeasible"); return true; }
    else { return false; }
}


Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p, ExactVector y) const
{
    return this->validate_infeasibility(p,Vector<Bounds<FLT>>(y,pr));
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p, Vector<Exact<FLT>> y) const
{
    return this->validate_infeasibility(p,Vector<Bounds<FLT>>(y));
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p, Vector<Bounds<FLT>> y) const
{
    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C<<", y="<<y);

    ARIADNE_PRECONDITION(y.size()==C.size());

    if(y.size()==0) { return D.is_empty(); }

    UpperIntervalType yC = dot(y,C);

    // Compute Taylor estimate of y g(X)
    ValidatedVectorMultivariateTaylorFunctionModelDP tg(D,g,default_sweeper());
    ValidatedScalarMultivariateTaylorFunctionModelDP tyg(D,default_sweeper());
    for(SizeType j=0; j!=y.size(); ++j) {
        tyg+=y[j]*tg[j];
    }
    UpperIntervalType ygD = apply(tyg.function(),D);
    // UpperIntervalType ygD = dot(y,apply(g,D));

    if(definitely(disjoint(yC,ygD))) {
        CONCLOG_PRINTLN("Infeasible");
        return true;
    } else {
        return false;
    }
}

Bool FeasibilityChecker::
validate_infeasibility(ValidatedFeasibilityProblem p) const
{
    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C);

    if (D.is_empty()) {
        return true;
    }

    UpperBoxType W=apply(g,D);
    return (definitely(disjoint(W,C)));
}


OutputStream& operator<<(OutputStream& os, FeasibilityChecker const& fc) {
    return os << "FeasibilityChecker()";
}



//------- OptimiserBase -----------------------------------//

const OptimiserBase<ApproximateTag>::FLT OptimiserBase<ApproximateTag>::zero = FLT(0,pr);
const OptimiserBase<ApproximateTag>::FLT OptimiserBase<ApproximateTag>::one = FLT(1,pr);
const OptimiserBase<ValidatedTag>::FLT OptimiserBase<ValidatedTag>::zero = FLT(0,pr);
const OptimiserBase<ValidatedTag>::FLT OptimiserBase<ValidatedTag>::one = FLT(1,pr);


auto OptimiserBase<ValidatedTag>::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const -> ValidatedVector
{
    CONCLOG_SCOPE_CREATE;
    return this->minimise(ValidatedOptimisationProblem{f,D,g,C});
}

auto OptimiserBase<ValidatedTag>::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const -> ValidatedVector
{
    CONCLOG_SCOPE_CREATE;
    ValidatedVectorMultivariateFunction gh=join(g,h);
    ExactBoxType C(gh.result_size(),ExactIntervalType(0,0));
    for(SizeType i=0; i!=g.result_size(); ++i) { C[i]=ExactIntervalType(-inf,0); }
    return this->minimise(f,D,gh,C);
}

auto OptimiserBase<ApproximateTag>::
minimise(ApproximateScalarMultivariateFunction f, ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateBoxType C) const -> ApproximateVector
{
    CONCLOG_SCOPE_CREATE;
    return this->minimise(ApproximateOptimisationProblem{f,D,g,C});
}

auto OptimiserBase<ApproximateTag>::
minimise(ApproximateScalarMultivariateFunction f, ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateVectorMultivariateFunction h) const -> ApproximateVector
{
    CONCLOG_SCOPE_CREATE;
    ApproximateVectorMultivariateFunction gh=join(g,h);
    ApproximateBoxType C(gh.result_size(),ExactIntervalType(0,0));
    for(SizeType i=0; i!=g.result_size(); ++i) { C[i]=ExactIntervalType(-inf,0); }
    return this->minimise(f,D,gh,C);
}



auto OptimiserBase<ValidatedTag>::
feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const -> ValidatedKleenean
{
    CONCLOG_SCOPE_CREATE;
    return this->feasible(ValidatedFeasibilityProblem{D,g,C});
}

auto OptimiserBase<ValidatedTag>::
feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const -> ValidatedKleenean
{
    CONCLOG_SCOPE_CREATE;
    ValidatedVectorMultivariateFunction gh=join(g,h);
    ExactBoxType C(gh.result_size(),ExactIntervalType(0,0));
    for(SizeType i=0; i!=g.result_size(); ++i) { C[i]=ExactIntervalType(-inf,0); }
    return this->feasible(D,gh,C);
}

auto OptimiserBase<ApproximateTag>::
feasible(ApproximateBoxType D, ApproximateVectorMultivariateFunction g, ApproximateBoxType C) const -> ApproximateKleenean
{
    CONCLOG_SCOPE_CREATE;
    return this->feasible(ApproximateFeasibilityProblem{D,g,C});
}


//------- PenaltyFunctionOptimiser ------------------------------------------//

PenaltyFunctionOptimiser* PenaltyFunctionOptimiser::
clone() const
{
    return new PenaltyFunctionOptimiser(*this);
}

auto PenaltyFunctionOptimiser::
minimise(ApproximateOptimisationProblem p) const ->ApproximateVector
{
    ARIADNE_NOT_IMPLEMENTED;
}

auto PenaltyFunctionOptimiser::
feasible(ValidatedFeasibilityProblem p) const -> ValidatedKleenean
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<p.D<<" g="<<p.g<<" C="<<p.C<<" ");

    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    ApproximateVectorType x=midpoint(D);

    ApproximateVectorType w=midpoint(C);
    for(SizeType i=0; i!=C.size(); ++i) {
        if(C[i].upper_bound()==+infty) { w[i]=C[i].lower_bound()+one; }
        else if(C[i].lower_bound()==-infty) { w[i]=C[i].upper_bound()-one; }
    }

    ApproximateVectorType y(C.size(),zero);

    CONCLOG_PRINTLN("x="<<x<<" w="<<w<<" y="<<y);

    for(SizeType i=0; i!=10; ++i) {
        this->feasibility_step(ApproximateFeasibilityProblem(D,g,C),w,x,y);
    }
    return FeasibilityChecker().check_feasibility(p,cast_exact(x),cast_exact(y));
}

auto PenaltyFunctionOptimiser::
feasible(ApproximateFeasibilityProblem p) const -> ApproximateKleenean
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<p.D<<" g="<<p.g<<" C="<<p.C<<" ");

    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    ApproximateVectorType x=midpoint(D);

    ApproximateVectorType w=midpoint(C);
    for(SizeType i=0; i!=C.size(); ++i) {
        if (decide(C[i].upper_bound()==+infty)) { w[i]=C[i].lower_bound()+one; }
        else if (decide(C[i].lower_bound()==-infty)) { w[i]=C[i].upper_bound()-one; }
    }

    ApproximateVectorType y(C.size(),zero);

    CONCLOG_PRINTLN("x="<<x<<" w="<<w<<" y="<<y);

    for(SizeType i=0; i!=10; ++i) {
        this->feasibility_step(p,w,x,y);
    }
    return C.contains(g(x));
}

Void PenaltyFunctionOptimiser::
feasibility_step(ApproximateFeasibilityProblem p,
                 ApproximateVectorType& w, ApproximateVectorType& x, ApproximateNumberType& mu) const
{
    CONCLOG_SCOPE_CREATE;
    auto& d=p.D; auto& g=p.g; auto& c=p.C;

    ApproximateVectorMultivariateFunction h(0u,d.dimension());
    const SizeType n=d.size();
    const SizeType m=c.size();
    const SizeType l=h.result_size();

    CONCLOG_PRINTLN("x="<<x);
    CONCLOG_PRINTLN("w="<<w);

    Vector<FloatDPApproximationDifferential> ddgx=g.evaluate(FloatDPApproximationDifferential::variables(2,x));
    Vector<FloatDPApproximationDifferential> ddhx=h.evaluate(FloatDPApproximationDifferential::variables(2,x));

    mu *= 0.5;
    CONCLOG_PRINTLN("mu="<<mu);

    // G is the constraint value vector
    ApproximateVectorType gx = ddgx.value();
    ApproximateVectorType hx = ddhx.value();
    CONCLOG_PRINTLN("g(x)="<<gx);
    CONCLOG_PRINTLN("h(x)="<<hx);

    // A is the transpose derivative matrix aij=dgi/dxj
    FloatDPApproximationMatrix A = transpose(ddgx.jacobian());
    CONCLOG_PRINTLN("A=Dg(x)="<<A);
    FloatDPApproximationMatrix B = transpose(ddhx.jacobian());
    // FIXME: Due to problems with zero-element differential, need to resize matrix if no h
    if(l==0) { B.resize(n,0); }
    CONCLOG_PRINTLN("B=Dh(x)="<<B);

    // H is the Hessian matrix H[i1,i2] = df/dx[i1]dx[i2] + Sum_[j] lambda[j]*dg[j]/dx[i1]dx[i2]
    FloatDPApproximationMatrix H(n,n,dp);
    for(SizeType j=0; j!=m; ++j) { H += (gx[j]-w[j]) * ddgx[j].hessian(); }
    for(SizeType k=0; k!=l; ++k) { H += (hx[k]) * ddhx[k].hessian(); }
    CONCLOG_PRINTLN("H="<<H);

    FloatDPApproximationDiagonalMatrix E(n,dp);
    FloatDPApproximationDiagonalMatrix D(m,dp);
    for(SizeType i=0; i!=n; ++i) { E[i] = rec(sqr(x[i]-d[i].lower_bound())) + rec(sqr(d[i].upper_bound()-x[i])); }
    for(SizeType j=0; j!=m; ++j) { D[j] = rec(sqr(w[j]-c[j].lower_bound())) + rec(sqr(c[j].upper_bound()-w[j])); }
    CONCLOG_PRINTLN("E="<<E);
    CONCLOG_PRINTLN("D="<<D);

    FloatDPApproximationMatrix S = H + B * transpose(B);
    S += E;
    CONCLOG_PRINTLN("S="<<S);

    FloatDPApproximationMatrix R=inverse(S);
    CONCLOG_PRINTLN("inverse(S)="<<R);

    // Compute residuals
    ApproximateVectorType rx = A*gx + B * hx ; // + 1/(x.upper_bound()-x) + 1/x.lower_bound()-x if no regularisation
    ApproximateVectorType rw = w-gx;

    CONCLOG_PRINTLN("rx="<<rx);
    CONCLOG_PRINTLN("rw="<<rw);

    ApproximateVectorType dx = R * (rx + A * rw);
    ApproximateVectorType dw = rw + transpose(A)*dx;
    CONCLOG_PRINTLN("dx="<<dx);
    CONCLOG_PRINTLN("dw="<<dw);


    ApproximateVectorType newx(n,dp);
    ApproximateVectorType neww(m,dp);

    static const FloatDPApproximation ALPHA_SCALE_FACTOR = 0.75_approx;

    FloatDPApproximation alpha = 1.0_approx;
    do {
        newx = x - alpha * dx;
        neww = w - alpha * dw;
        alpha *= ALPHA_SCALE_FACTOR;
    } while ( ! decide(contains(d,newx)) || ! decide(contains(c,neww)) );
    alpha /= ALPHA_SCALE_FACTOR;

    CONCLOG_PRINTLN("alpha="<<alpha);
    CONCLOG_PRINTLN("newx="<<newx<<", neww="<<neww);

    x=newx;
    w=neww;

    return;
}


// Use a penalty approach without multipliers on the constraint functions
// Solve g(x)=w, x in D, w in C; Lagrangian y.(g(x)-w)
Void PenaltyFunctionOptimiser::
feasibility_step(ApproximateFeasibilityProblem p,
                 ApproximateVectorType& w, ApproximateVectorType& x, ApproximateVectorType& y) const
{
    CONCLOG_SCOPE_CREATE;
    auto& D=p.D; auto& g=p.g; auto& C=p.C;
    auto m=y.size(); auto n=x.size();

    Vector<ApproximateNumberType> cl=lower_bounds(C);
    Vector<ApproximateNumberType> cu=upper_bounds(C);
    Vector<ApproximateNumberType> dl=lower_bounds(D);
    Vector<ApproximateNumberType> du=upper_bounds(D);

    CONCLOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C);
    CONCLOG_PRINTLN("dl="<<dl<<", du="<<du);
    CONCLOG_PRINTLN("cl="<<cl<<", cu="<<cu);
    CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y);

    ARIADNE_ASSERT_MSG(g.argument_size()==D.size(),"D="<<D<<", g="<<g<<", C="<<C);
    ARIADNE_ASSERT_MSG(g.result_size()==C.size(),  "D="<<D<<", g="<<g<<", C="<<C);
    ARIADNE_ASSERT(w.size()==m);
    ARIADNE_ASSERT(x.size()==n);
    ARIADNE_ASSERT(y.size()==m);

    // Solve the problem
    //   minimise Sum -log(w-cl)-log(cu-w)-log(x-dl)-log(du-x)
    //   subject to g(x)-w=0

    // Lagrangian
    //   -log(w-cl)-log(cu-w)-log(x-dl)-log(du-x) - y.(g(x)-w)

    // Conditions for a constrained minimum
    // 1/(cu-w)-1/(w-cu) + y       = 0
    // 1/(du-x)-1/(x-dl) - y.Dg(x) = 0
    //          w - g(x)           = 0

    // Second-order conditions
    // (1/(w-cl)^2 + 1/(cu-w)^2) dw                 +   dy = - ( 1/(cu-w) - 1/(w-cl) + y       )
    //   (1/(x-dl)^2 + 1/(du-x)^2 - y.D^2x) dx - Dg'(x) dy = - ( 1/(du-x) - 1/(x-dl) - y.Dg(x) )
    //                           dw - Dg(x) dx             = - ( w - g(x) )

    Vector<Differential<ApproximateNumberType>> ddgx=g.evaluate(Differential<ApproximateNumberType>::variables(2,x));
    CONCLOG_PRINTLN("ddgx="<<ddgx);

    Vector<ApproximateNumberType> gx = ddgx.value();
    CONCLOG_PRINTLN("g(x)="<<gx);
    Matrix<ApproximateNumberType> A = ddgx.jacobian();
    CONCLOG_PRINTLN("Dg(x)="<<A);

    Vector<ApproximateNumberType> yA=transpose(A)*y;

    // H is the Hessian matrix H of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) $
    Matrix<ApproximateNumberType> YH(x.size(),x.size(),dp);
    for(SizeType i=0; i!=y.size(); ++i) {
        YH+=y[i]*ddgx[i].hessian();
    }
    CONCLOG_PRINTLN("Y.D2g(x)="<<YH);

    Vector<ApproximateNumberType> recwu=cu-w; recwu=erec(recwu);
    Vector<ApproximateNumberType> recwl=w-cl; recwl=erec(recwl);
    Vector<ApproximateNumberType> recxu=du-x; recxu=erec(recxu);
    Vector<ApproximateNumberType> recxl=x-dl; recxl=erec(recxl);

    Vector<ApproximateNumberType> diagDw=esqr(recwu)+esqr(recwl);
    Matrix<ApproximateNumberType> Dw(m,m,dp); for(SizeType i=0; i!=m; ++i) { Dw[i][i]=diagDw[i]; }
    DiagonalMatrix<ApproximateNumberType> Dx(esqr(recxu)+esqr(recxl));


    for(SizeType i=0; i!=n; ++i) { YH[i][i]-=Dx[i]; }

    Matrix<ApproximateNumberType> AT=transpose(A);
    Matrix<ApproximateNumberType> Znm(n,m,dp);
    Matrix<ApproximateNumberType> Zmn(m,n,dp);
    Matrix<ApproximateNumberType> Zmm(m,m,dp);
    Matrix<ApproximateNumberType> Im=Matrix<ApproximateNumberType>::identity(m,dp);


    Matrix<ApproximateNumberType> S=cojoin(join(Dw,Zmn,Im),join(Znm,-YH,-AT),join(Im,-A,Zmm));
    Vector<ApproximateNumberType> r=join(recwu-recwl+y,recxu-recxl-yA,w-gx);

    for(SizeType j=0; j!=m; ++j) {
        if( decide(C[j].lower_bound()==C[j].upper_bound()) ) {
            S[j][j]=1;
            S[j][m+n+j]=0;
            S[m+n+j][j]=0;
            r[j]=0;
        }
    }

    Vector<ApproximateNumberType> swxy = -solve(S,r);

    Vector<ApproximateNumberType> sw(m,dp),sx(n,dp),sy(m,dp);
    sw = project(swxy,range(0,m));
    sx = project(swxy,range(m,m+n));
    sy = project(swxy,range(m+n,m+n+m));

    ApproximateNumberType al=one;
    Vector<ApproximateNumberType> nw=w+al*sw;
    Vector<ApproximateNumberType> nx=x+al*sx;
    Vector<ApproximateNumberType> ny(m,dp);
    CONCLOG_PRINTLN("sx="<<sx);
    CONCLOG_PRINTLN("sw="<<sw);
    while( ! decide(contains(C,nw)) || ! decide(contains(D,nx)) ) {
        al*=0.75;
        nw=w+al*sw;
        nx=x+al*sx;
    }
    CONCLOG_PRINTLN("al="<<sw);
    ny=y+al*sy;

    w=nw; x=nx; y=ny;
}


//------- InteriorPointOptimiserBase ------------------------------------------//


// Given the constraint w in C, dual variables y, and centering parameter mu
// Compute the residual r(w,x), and the partial derivatives -dw=dr/dw and cdy=dr/dy
template<class FLT>
Void compute_linear_equation(const ApproximateIntervalType& C, const Approximation<FLT>& w, const Approximation<FLT>& y, const Approximation<FLT>& mu,
                             Approximation<FLT>& cdw, Approximation<FLT>& cdy, Approximation<FLT>& r)
{
    auto& Cl=C.lower_bound(); auto& Cu=C.upper_bound();
    if (decide(Cl==Cu)) {
        // w-Clu=0
        r=w-Cl; cdw=1; cdy=0;
    } else if (decide(Cu==+inf)) {
        // (w-Cl)*y+mu=0
        r=(w-Cl)*y+mu; cdw=y; cdy=w-Cl;
    } else if (decide(Cl==-inf)) {
        // (Cu-w)*y-mu=0
        r=(Cu-w)*y-mu; cdw=-y; cdy=Cu-w;
    } else {
        // (w-Cl)*(Cu-w)*y+mu*(Cl+Cu-2*w)=0
        r=(w-Cl)*(Cu-w)*y+mu*(Cl+Cu-2*w); cdw=(Cl+Cu-2*w)*y-2*mu; cdy=(w-Cl)*(Cu-w);
    }
}

template<class FLT>
Approximation<FLT> compute_complementarity_residual(const ApproximateIntervalType& C, const Approximation<FLT>& w, const Approximation<FLT>& y, const Approximation<FLT>& mu)
{
    auto& Cl=C.lower_bound(); auto& Cu=C.upper_bound();
    if (decide(Cl==Cu)) {
        return w-Cl; // w-Clu=0
    } else if (decide(Cu==+inf)) {
        return (w-Cl)*y+mu; // (w-Cl)*y+mu=0
    } else if (decide(Cl==-inf)) {
        return (Cu-w)*y-mu; // (Cu-w)*y-mu=0
    } else {
        return (w-Cl)*(Cu-w)*y+mu*(Cl+Cu-2*w); // (w-Cl)*(Cu-w)*y+mu*(Cl+Cu-2w)=0
    }
}

template<class X>
Void compute_linear_equation(const ExactIntervalType& C, const X& w, const X& y,
                             ArithmeticType<X>& cdw, ArithmeticType<X>& cdy, ArithmeticType<X>& r)
{
    auto& Cl=C.lower_bound(); auto& Cu=C.upper_bound();
    if (Cl==Cu) {
        // w-Clu=0
        r=w-Cl; cdw=1; cdy=0;
    } else if (Cu==+inf) {
        // (w-Cl)*y+mu=0
        r=(w-Cl)*y; cdw=y; cdy=w-Cl;
    } else if (Cl==-inf) {
        // (Cu-w)*y-mu=0
        r=(Cu-w)*y; cdw=-y; cdy=Cu-w;
    } else {
        // (w-Cl)*(Cu-w)*y-mu*(Cl+Cu-2*w)=0
        r=(w-Cl)*(Cu-w)*y; cdw=(Cl+Cu-2*w)*y; cdy=(w-Cl)*(Cu-w);
    }
}

template<class X, class BX>
Tuple<DiagonalMatrix<X>,DiagonalMatrix<X>,Vector<X>>
compute_linear_equations(const BX& C, const Vector<X>& w, const Vector<X>& y, const X& mu)
{
    const SizeType m=C.dimension(); const X zero(mu.precision());
    DiagonalMatrix<X> Y(m,zero), W(m,zero); Vector<X> r(m,zero);
    for (SizeType j=0; j!=m; ++j) {
        compute_linear_equation(C[j],w[j],y[j],mu, Y._at(j), W._at(j), r.at(j));
    }
    return make_tuple(Y,W,r);
}

template<class X, class BX>
Vector<X> compute_complementarity_residuals(const BX& C, const Vector<X>& w, const Vector<X>& y, const X& mu)
{
    const SizeType m=C.dimension(); const X zero(mu.precision());
    return Vector<X>(C.dimension(), [&](SizeType j){return compute_complementarity_residual(C[j],w[j],y[j],mu);});
}

template<class X, class BX>
Tuple<DiagonalMatrix<ArithmeticType<X>>,DiagonalMatrix<ArithmeticType<X>>,Vector<ArithmeticType<X>>>
compute_linear_equations(const BX& C, const Vector<X>& w, const Vector<X>& y)
{
    const SizeType m=C.dimension(); const X zero(w.zero_element().precision());
    DiagonalMatrix<ArithmeticType<X>> Y(m,zero), W(m,zero); Vector<ArithmeticType<X>> r(m,zero);
    for (SizeType j=0; j!=m; ++j) {
        compute_linear_equation(C[j],w[j],y[j], Y._at(j), W._at(j), r.at(j));
    }
    return make_tuple(Y,W,r);
}

auto
InteriorPointOptimiserBase::compute_t(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x) const
        -> Approximation<FLT>
{
    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    // Compute a the minimum constraint satisfaction
    Approximation<FLT> t=zero; t=+infty;
    ApproximateVectorType gx = g(x);

    for(SizeType i=0; i!=D.size(); ++i) {
        t = min( t, min(x[i]-D[i].lower_bound(), D[i].upper_bound()-x[i]) );
    }
    for(SizeType j=0; j!=C.size(); ++j) {
        t = min( t, min(gx[j]-C[j].lower_bound(), C[j].upper_bound()-gx[j]) );
    }
    return t;
}

auto
InteriorPointOptimiserBase::minimise(ApproximateOptimisationProblem p) const -> ApproximateVector
{
    CONCLOG_SCOPE_CREATE;

    Approximation<FLT> mu0(1,pr);
    auto x0 = midpoint(p.D);
    auto y0 = this->compute_y(p,x0,mu0);
    return minimise_hotstarted(p,x0,y0).primal();
}



auto
InteriorPointOptimiserBase::minimise_hotstarted(
    const ApproximateOptimisationProblem& p,
    const Vector<Approximation<FLT>>& x0, const Vector<Approximation<FLT>>& y0) const
        -> ValuePrimalDualData<ApproximateNumber>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;

    static const ApproximateDouble VALUE_TOLERANCE=1e-8;
    static const ApproximateDouble STATE_TOLERANCE=1e-8;
    static const ApproximateDouble MU_MIN = 1e-12;

    static const CounterType MAXIMUM_STEPS=24;

    UniquePointer<StepData> d_ptr(this->initial_step_data_hotstarted(p,x0,y0));
    StepData& d=*d_ptr;

    Vector<Approximation<FLT>> x=d.primal();

    Vector<Approximation<FLT>> oldx(x.size(),x.zero_element());

    for(SizeType i=0; i!=MAXIMUM_STEPS; ++i) {
        CONCLOG_PRINTLN_AT(1,"step="<<i);
        CONCLOG_PRINTLN_AT(1,d);
        oldx=x;
        this->_minimisation_step(p, d);
        x=d.primal();
        if(probably(mag(f(x)-f(x))<VALUE_TOLERANCE) && probably(norm(oldx-x)<STATE_TOLERANCE)) {
            break;
        }
        if(probably(d.mu<MU_MIN)) {
            break;
        }
        if (2<=i && i+2<=MAXIMUM_STEPS) {
            d.mu *= 0.25_exact;
        }
    }
    CONCLOG_PRINTLN(d);

    if (probably(D.contains(x)) && probably(C.contains(g(x)))) {
        Vector<Approximation<FLT>> y=d.dual();
        CONCLOG_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        return ValuePrimalDualData<ApproximateNumber>(f(x),x,y);
    } else {
        CONCLOG_PRINTLN("indeterminate_feasibility")
        throw IndeterminateFeasibilityException();
    }
}


auto
InteriorPointOptimiserBase::feasible(ApproximateFeasibilityProblem p) const -> ApproximateKleenean
{
    ApproximateVectorType x0 = p.D.midpoint();
    ApproximateVectorType y0(p.g.result_size(),zero);
    return this->feasible_hotstarted(p,x0,y0).is_feasible();
}



auto
InteriorPointOptimiserBase::feasible_hotstarted(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> FeasiblePrimalDualData<ApproximateNumber>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& g=p.g;

    ARIADNE_ASSERT(x0.size()==p.number_of_variables());
    ARIADNE_ASSERT(y0.size()==p.number_of_constraints());

    //Vector<Approximation<FLT>> x(x0,pr);
    //Vector<Approximation<FLT>> y(y0,pr);

    UniquePointer<StepData> d_ptr(this->initial_step_data_hotstarted(p,x0,y0));
    StepData& d=*d_ptr;

    // FIXME: Allow more steps
    Vector<Approximation<FLT>> x=d.primal();
    Vector<Approximation<FLT>> y=d.dual();
    Approximation<FLT> t=this->compute_t(p,x);

    for (SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"t="<<t<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        this->feasibility_step(p, d);
        x=d.primal();
        y=d.dual();
        t=this->compute_t(p,x);
        if (probably(LogicalValue(t>0))) {
            return FeasiblePrimalDualData<ApproximateNumber>(true,x,y);
        }
    }
    CONCLOG_PRINTLN("x="<<x<<", y="<<y<<", g(x)="<<g(x));
    return FeasiblePrimalDualData<ApproximateNumber>(false,x,y);
}


auto
InteriorPointOptimiserBase::initial_step_data(const ApproximateFeasibilityProblem& p) const -> StepData*
{
    CONCLOG_PRINTLN("InteriorPointOptimiserBase::initial_step_data(ApproximateFeasibilityProblem p)");
    CONCLOG_PRINTLN("p="<<p);

    Vector<Approximation<FLT>> x0=midpoint(p.D);
    Vector<Approximation<FLT>> y0(p.number_of_constraints(),zero);
    for (SizeType i=0; i!=p.number_of_constraints(); ++i) {
        if (decide(p.C[i].lower_bound()==inf)) {
            y0[i]=1;
        } else if (decide(p.C[i].upper_bound()==inf)) {
            y0[i]=-1;
        }
    }
    return this->_initial_step_data_hotstarted(p,x0,y0);
}

auto
InteriorPointOptimiserBase::initial_step_data_hotstarted(const ApproximateFeasibilityProblem& p, const ApproximateVectorType& x0, const ApproximateVectorType& y0) const -> StepData*
{
    return this->_initial_step_data_hotstarted(p,x0,y0);
}

auto
InteriorPointOptimiserBase::minimisation_step(const ApproximateOptimisationProblem& p, StepData& d) const -> Void
{
    this->_minimisation_step(p,d);
}

auto
InteriorPointOptimiserBase::feasibility_step(const ApproximateFeasibilityProblem& p, StepData& d) const -> Void
{
    this->_feasibility_step(p,d);
}

auto
InteriorPointOptimiserBase::_feasibility_step(const ApproximateFeasibilityProblem& p, StepData& d) const -> Void
{
    auto f=ApproximateScalarMultivariateFunction::constant(p.number_of_variables(),0);
    ApproximateOptimisationProblem optp(f,p.D,p.g,p.C);
    this->_minimisation_step(optp,d);
}


auto
InteriorPointOptimiserBase::compute_dual(
    const ApproximateBoxType& D,
    const ApproximateVectorType& x, const Approximation<FLT>& mu) const
        -> ApproximateVectorType
{
    ApproximateVectorType z(x.size(),x.zero_element());
    for(SizeType i=0; i!=D.size(); ++i) {
        if (decide(D[i].lower_bound() == D[i].upper_bound())) { }
        else if (decide(D[i].lower_bound()==-inf)) { z[i] = mu / (D[i].upper_bound()-x[i]); }
        else if (decide(D[i].upper_bound()==+inf)) { z[i] = -mu / (x[i]-D[i].lower_bound()); }
        else { z[i] = mu * ( rec(D[i].upper_bound()-x[i]) - rec(x[i]-D[i].lower_bound()) ); }
    }
    return z;
}


auto
InteriorPointOptimiserBase::compute_x(const ApproximateFeasibilityProblem& p) const
        -> ApproximateVectorType
{
    return p.D.midpoint();
}

auto
InteriorPointOptimiserBase::compute_y(const ApproximateFeasibilityProblem& p, const ApproximateVectorType& x, const Approximation<FLT>& mu) const
        -> ApproximateVectorType
{
    return this->compute_dual(p.C,p.g(x),mu);
}

auto
InteriorPointOptimiserBase::compute_w(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x, const ApproximateVectorType& y, const Approximation<FLT>& mu) const
        -> ApproximateVectorType
{
    auto& C=p.C;
    ApproximateVectorType w=midpoint(C);
    for (SizeType j=0; j!=w.size(); ++j) {
        if (decide(C[j].lower_bound()==C[j].upper_bound())) {
            w[j]=C[j].lower_bound();
        } else if (decide(C[j].lower_bound()==-inf)) {
            ARIADNE_ASSERT (decide(y[j]>0));
            w[j]=C[j].upper_bound()-mu/y[j];
        } else if (decide(C[j].upper_bound()==+inf)) {
            ARIADNE_ASSERT (decide(y[j]<0));
            w[j]=C[j].lower_bound()-mu/y[j];
        } else {
            auto m=C[j].midpoint();
            auto r=C[j].radius();
            auto k=r*y[j]/mu;
            w[j]=m+r*(k/(sqrt(1+sqr(k))+1));
        }
    }
    return w;
}


auto
InteriorPointOptimiserBase::compute_w(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x, const Approximation<FLT>& mu) const
        -> ApproximateVectorType
{
    static const ApproximateDouble GAMMA=0.0009765625; // 1.0/1024;

    auto& C=p.C;
    auto gx=p.g(x);

    ApproximateVectorType w=midpoint(C);
    for (SizeType j=0; j!=w.size(); ++j) {
        if (decide(C[j].lower_bound()==C[j].upper_bound())) {
            w[j]=C[j].lower_bound();
        } else if (decide(C[j].lower_bound()==-inf)) {
            w[j]=min(gx[j],C[j].upper_bound()-mu*GAMMA);
        } else if (decide(C[j].upper_bound()==+inf)) {
            max(gx[j], C[j].lower_bound()+mu*GAMMA);
        } else {
            if (decide(C[j].radius()<=mu*GAMMA)) {
                w[j]=C[j].midpoint();
            } else {
                w[j]=min( max( gx[j], C[j].lower_bound()+mu*GAMMA ), C[j].upper_bound()-mu*GAMMA );
            }
        }
    }
    return w;
}

auto
InteriorPointOptimiserBase::compute_z(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x, const Approximation<FLT>& mu) const
        -> ApproximateVectorType
{
    return this->compute_dual(p.D,x,mu);
}

auto
InteriorPointOptimiserBase::compute_mu(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x, const ApproximateVectorType& y) const
        -> ApproximateNumberType
{
    auto& g=p.g; auto& C=p.C;

    // Compute the relaxation parameter mu as the average of the product of the Lyapunov exponents and constraint satisfactions
    Approximation<FLT> mu=zero;
    ApproximateVectorType gx = g(x);

    for(SizeType i=0; i!=C.size(); ++i) {
        if (decide(C[i].lower_bound()==C[i].upper_bound())) { }
        else if (decide(C[i].lower_bound()==-infty)) { mu += y[i] * (gx[i] - C[i].upper_bound()); }
        else if (decide(C[i].upper_bound()==+infty)) { mu += y[i] * (gx[i] - C[i].lower_bound()); }
        else {
            if ( decide(y[i] <=0.0) ) { mu += y[i] * (gx[i] - C[i].upper_bound()); }
            else { mu += y[i] * (gx[i] - C[i].lower_bound()); }
        }
    }
    mu /= C.size();
    return mu;
}


Void InteriorPointOptimiserBase::compute_tz(
    const ApproximateFeasibilityProblem& p,
    ApproximateVectorType& x, ApproximateNumberType& t, ApproximateVectorType& z) const
{
    t = this->compute_t(p,x);
    z = this->compute_z(p,x,one);
}


//------- ApproximateOptimiser ------------------------------------------//

ApproximateOptimiser* ApproximateOptimiser::clone() const {
    return new ApproximateOptimiser(*this);
}

ValidatedKleenean ApproximateOptimiser::
feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<D<<", h="<<h);
    ApproximateVectorType x=midpoint(D);
    ApproximateVectorType y(h.result_size(),zero);

    for(SizeType i=0; i!=8; ++i) {
        this->feasibility_step(D,h,x,y);
    }

    if( decide(norm(h(x))<1e-10) ) { return true; }

    if(!possibly(contains(UpperIntervalType(dot(UpperIntervalVectorType(cast_exact(y)),apply(h,D))),zero))) { return false; }

    return indeterminate;
}

Void ApproximateOptimiser::
feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& h,
                 ApproximateVectorType& x, ApproximateVectorType& y) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("x="<<x<<" y="<<y);
    static const double SCALE_FACTOR = 0.75;
    const SizeType n=x.size();
    const SizeType m=y.size();
    // Solve equations y Dh(x) - 1/(x-xl) + 1/(xu-x) = 0; h(x) = 0
    Vector<FloatDPApproximationDifferential> ddhx=h.evaluate(FloatDPApproximationDifferential::variables(2,x));
    FloatDPApproximationMatrix A = ddhx.jacobian();
    CONCLOG_PRINTLN_AT(1,"A="<<A<<" b="<<ddhx.value());

    FloatDPApproximationMatrix H(n,n,dp);
    for(SizeType i=0; i!=m; ++i) { H += y[i] * ddhx[i].hessian(); }
    for(SizeType j=0; j!=n; ++j) {
        H[j][j] += rec(sqr(x[j]-D[j].lower_bound()));
        H[j][j] += rec(sqr(D[j].upper_bound()-x[j]));
    }

    ApproximateVectorType rx = transpose(A) * y;
    for(SizeType j=0; j!=n; ++j) {
        rx[j] -= rec(x[j]-D[j].lower_bound());
        rx[j] += rec(D[j].upper_bound()-x[j]);
    }
    ApproximateVectorType ry = ddhx.value();
    CONCLOG_PRINTLN("rx="<<rx<<" ry="<<ry);

    // S = A Hinv AT
    // H dx + AT dy = rx; A dx = ry;
    //  dx = Hinv ( rx - AT dy )
    //  dy = Sinv ( A Hinv rx - ry )
    FloatDPApproximationMatrix Hinv=inverse(H);
    CONCLOG_PRINTLN_AT(1,"H="<<H<<" Hinv="<<Hinv);
    FloatDPApproximationMatrix S=A*Hinv*transpose(A);
    FloatDPApproximationMatrix Sinv=inverse(S);
    CONCLOG_PRINTLN_AT(1,"S="<<S<<" Sinv="<<Sinv);
    ApproximateVectorType dy = Sinv * ( A*(Hinv*rx) - ry );
    ApproximateVectorType dx = Hinv * ( rx - transpose(A) * dy);
    CONCLOG_PRINTLN("dx="<<dx<<" dy="<<dy);

    FloatDPApproximation ax = one;
    ApproximateVectorType nx = x-ax*dx;
    while(!contains(D,cast_exact(nx))) {
        ax*=SCALE_FACTOR;
        nx = x - ax * dx;
    }
    ApproximateVectorType ny = y-ax*dy;
    CONCLOG_PRINTLN("nx="<<nx<<" ax="<<ax<<" ny="<<ny);
    CONCLOG_PRINTLN_AT(1,"h(x)="<<h(nx));

    x=nx; y=ny;
}





//------- ApproximateOptimiserBase -----------------------------------//

//------- PrimalDualInteriorPointOptimiser -------------------------//

auto PrimalDualInteriorPointOptimiser::clone() const -> PrimalDualInteriorPointOptimiser* {
    return new PrimalDualInteriorPointOptimiser(*this);
}

OutputStream& operator<<(OutputStream& os, PrimalDualInteriorPointOptimiser const& opt) {
    return os << "PrimalDualInteriorPointOptimiser()";
}

auto
PrimalDualInteriorPointOptimiser::initial_step_data(const ApproximateFeasibilityProblem& p) const -> StepData*
{
    return dynamic_cast<StepData*>(this->InteriorPointOptimiserBase::initial_step_data(p));
}

auto
PrimalDualInteriorPointOptimiser::initial_step_data_hotstarted(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const -> StepData*
{
    CONCLOG_SCOPE_CREATE;

    return new StepData(x0,y0);
}

auto
PrimalDualInteriorPointOptimiser::_initial_step_data_hotstarted(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const -> StepDataBase*
{
    return this->initial_step_data_hotstarted(p,x0,y0);
}

Void
PrimalDualInteriorPointOptimiser::_minimisation_step(
    const ApproximateOptimisationProblem& p,
    StepDataBase& d) const
{
    this->minimisation_step(p,dynamic_cast<StepData&>(d));
}

Void
PrimalDualInteriorPointOptimiser::minimisation_step(
    const ApproximateOptimisationProblem& p,
    StepData& d) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("PrimalDualInteriorPointOptimiser::minimisation_step(p,d)");

    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;
    CONCLOG_PRINTLN("d="<<d);

    auto Dl=lower_bounds(D);
    auto Du=upper_bounds(D);
    auto Cl=lower_bounds(C);
    auto Cu=upper_bounds(C);

    static const ApproximateDouble GAMMA=0.0009765625; // 1.0/1024;
    static const ApproximateDouble SIGMA=0.125;
    static const ApproximateDouble SCALE=0.75;

    const SizeType m=p.number_of_constraints();
    const SizeType n=p.number_of_variables();

    // Consider central path for min f(x) - mu * Sum(log(g(x)-Cl)+log(Cu-g(x))) - mu * Sum(log(x-Dl)+log(Du-x))
    // Criticality conditions:
    //   grad{f}(x) - mu * Sum(1/(Cu-g(x))-1/(g(x)-Cl))*grad{g}(x) - mu * Sum(1/(Du-x)-1/(x-Dl)) = 0
    // Set:
    //   y = -mu * (1/(Cu-g(x))-1/(g(x)-Cl))
    // Optimality conditions
    //   Let w=g(x), z=-mu * Sum(1/(Du-x)-1/(x-Dl)) in
    //     grad{f}(x) + y*grad{g}(x) + z = 0
    //     (Cu-w)*(w-Cl)*y + mu*(2*w-Cl-Cu) = 0


    Vector<Approximation<FLT>>& x=d.x;
    Vector<Approximation<FLT>>& y=d.y;
    Approximation<FLT>& mu=d.mu;

    ARIADNE_ASSERT(x.size()==n);
    ARIADNE_ASSERT(y.size()==m);

    FloatDPApproximationDifferential ddfx=f.evaluate(FloatDPApproximationDifferential::variables(2,x));
    CONCLOG_PRINTLN("ddfx="<<ddfx);
    Vector<FloatDPApproximationDifferential> ddgx=g.evaluate(FloatDPApproximationDifferential::variables(2,x));
    CONCLOG_PRINTLN("ddgx="<<ddgx);

    Vector<Approximation<FLT>> w = ddgx.value();
    CONCLOG_PRINTLN("w=g(x)="<<w<<" ");
    Covector<FloatDPApproximation> c = ddfx.gradient();
    CONCLOG_PRINTLN("c="<<c<<" ");
    Matrix<FloatDPApproximation> A = ddgx.jacobian();
    CONCLOG_PRINTLN("A="<<A<<" ");

    Vector<FloatDPApproximation> z(n,dp);
    for(SizeType i=0; i!=n; ++i) {
        z[i]=-mu*(rec(Du[i]-x[i])-rec(x[i]-Dl[i]));
    }

    // Q is the Hessian matrix of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) \lambda_k$
    Matrix<FloatDPApproximation> Q=ddfx.hessian();
    for(SizeType i=0; i!=m; ++i) {
        Q+=y[i]*ddgx[i].hessian();
    }
    // Add correction for D to diagonal elements of Hessian
    for(SizeType i=0; i!=n; ++i) {
        Q[i][i]-=mu*(rec(sqr(Du[i]-x[i]))+rec(sqr(x[i]-Dl[i])));
    }
    CONCLOG_PRINTLN("Q="<<Q<<" ");

    DiagonalMatrix<FloatDPApproximation> Y(m,dp);
    DiagonalMatrix<FloatDPApproximation> W(m,dp);
    for (SizeType j=0; j!=m; ++j) {
        Y._at(j)=(Cl[j]+Cu[j]-2*w[j])*y[j]+2*w[j];
        W._at(j)=(Cu[j]-w[j])*(w[j]-Cl[j]);
    }
    CONCLOG_PRINTLN("W="<<W<<", Y="<<Y);


    Vector<FloatDPApproximation> rx=transpose(c)+transpose(A)*y+z;
    Vector<FloatDPApproximation> ry=emul(Cu-w,w-Cl,y)+mu*(2*w-Cl-Cu);


    FloatDPApproximation sigma(SIGMA,dp);
    if(!egtr(emul(w,y),GAMMA*mu)) {
        CONCLOG_PRINTLN("WARNING: near-degeneracy in Lyapunov multipliers in interior-point solver:");
        CONCLOG_PRINTLN("w=g(x)="<<w<<", x="<<x<<", y="<<y<<", z="<<z);
    }

    // Solve the system Q*dx+AT*dy=-rx; Y*A*dx+W*dy=-ry
    // Symmetrise A*dx+W/Y*dy=-ry/Y
    // Simplify to (Q-AT*(Y/W)*A*dx=-(rx-AT/W*ry)
    FloatDPApproximationMatrix S=qatda(Q,A,-Y/W);
    CONCLOG_PRINTLN("S="<<S);
    FloatDPApproximationMatrix Sinv=inverse(S);
    CONCLOG_PRINTLN("Sinv="<<Sinv);

    // FIXME: What if S is not invertible?
    ApproximateVectorType rr=rx - transpose(A)*(ry/W);

    // Compute the differences
    ApproximateVectorType dx=-solve(S,rr);
    ApproximateVectorType dy=-(ry+Y*(A*dx))/W;
    CONCLOG_PRINTLN("dx="<<dx<<" dy="<<dy);

    ApproximateVectorType nw,nx,ny;

    // Since we need to keep the point feasible, but the updates are linear
    // we need to validate feasibility directly rather than assuming the
    // linear update of y and z are good enough.
    Bool allfeasible=false;
    FloatDPApproximation alpha=1/FloatDPApproximation(SCALE,dp);
    if(!egtr(emul(x,z) , GAMMA*mu/16)) {
        CONCLOG_PRINTLN("WARNING: x="<<x<<", z="<<z<< ", x.z="<<emul(x,z)<<"<"<<GAMMA*mu / 16);
        throw NearBoundaryOfFeasibleDomainException();
    }
    while(!allfeasible) {
        alpha=alpha*SCALE;
        nx=x+alpha*dx;
        nw=g(nx);
        allfeasible = egtr(nx-Dl,GAMMA*mu) && egtr(Du-nx,GAMMA*mu) && egtr(nw-Cl,GAMMA*mu) && egtr(Cu-nw,GAMMA*mu);
    }
    CONCLOG_PRINTLN("alpha="<<alpha);
    CONCLOG_PRINTLN("nx="<<nx<<" ny="<<ny<<" nw="<<nw);

    x=nx; y=ny;
}


// Perform one step of the optmisation problem
//   minimise -t subject to x in D,
//     gj(x)-t>=Cjl, gj(x)+t<=Cju.
//     gj(x)-t-Cjl>=0, Cju-gj(x)-t>=0
//   Assume Cl+t<=g(x)<=Cu-t

// Solve the central path equations with w=g(x)
//   (yu-yl).Dg(x)+z=0
//   sum (yl[j]+yu[j])=1
//   (w-t-Cl).*(Cu-w-t)*(yu-yl)-mu*(2*w-Cl-Cu)=0
//   (x-Dl).*(Du-x)*z-mu*(2*x-Dl-Du)=0
Void
PrimalDualComplementaryInteriorPointOptimiser::error_feasibility_step(
    const ApproximateFeasibilityProblem& p,
    Vector<Approximation<FLT>>& x, Vector<Approximation<FLT>>& yl, Vector<Approximation<FLT>>& yu, Vector<Approximation<FLT>>& z,
    Approximation<FLT>& t) const
{

    CONCLOG_SCOPE_CREATE;

    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    auto Cl=lower_bounds(C); auto Cu=upper_bounds(C);

    CONCLOG_PRINTLN("x="<<x<<", t="<<t<<", yl="<<yl<<", yu="<<yu<<", z="<<z);

    static const ApproximateDouble GAMMA=0.0009765625; // 1.0/1024;
    static const ApproximateDouble SIGMA=0.125;
    static const ApproximateDouble SCALE=0.75;

    const SizeType m=p.number_of_constraints();
    const SizeType n=p.number_of_variables();
    CONCLOG_PRINTLN("m="<<m<<" n="<<n);

    ARIADNE_ASSERT(x.size()==n);
    ARIADNE_ASSERT(yl.size()==m);
    ARIADNE_ASSERT(yu.size()==m);
    ARIADNE_ASSERT(z.size()==n);

    Vector<Approximation<FLT>> y=yu-yl;
    CONCLOG_PRINTLN("y=yu-yl"<<y);

    Vector<Differential<Approximation<FLT>>> ddgx=g.evaluate(Differential<Approximation<FLT>>::variables(2,x));
    CONCLOG_PRINTLN("ddgx="<<ddgx);

    Vector<Approximation<FLT>> w = ddgx.value();
    CONCLOG_PRINTLN("w=g(x)="<<w<<" ");
    Matrix<FloatDPApproximation> A = ddgx.jacobian();
    CONCLOG_PRINTLN("A="<<A<<" ");

    // Q is the Hessian matrix of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) \lambda_k$
    SymmetricMatrix<Approximation<FLT>> Q(n,dp);
    for(SizeType i=0; i!=m; ++i) {
        Q+=y[i]*ddgx[i].hessian();
    }
    CONCLOG_PRINTLN("Q="<<Q<<" ");

    DiagonalMatrix<Approximation<FLT>> W(m,dp),Y(m,dp);
    // Compute diagonal entries of KKT Hessian
    for(SizeType j=0; j!=m; ++j) {
        if (decide(C[j].lower_bound()==C[j].upper_bound())) {

        } else if (decide(C[j].upper_bound()==+inf)) {
        } else if (decide(C[j].lower_bound()==-inf)) {
        } else {
            ARIADNE_DEBUG_ASSERT(decide(-infty<C[j].lower_bound() && C[j].lower_bound()<C[j].upper_bound() && C[j].upper_bound()<+infty));
            #warning
            //Q[i][i]-=mu*(rec(sqr(Du[i]-t-x[i]))+rec(sqr(x[i]-Dl-t[i])));
        }
    }

    FloatDPApproximation sigma(SIGMA,dp);
    FloatDPApproximation mu=dot(x,z)/n;
    if(!egtr(emul(x,z),GAMMA*mu)) {
        CONCLOG_PRINTLN("WARNING: near-degeneracy in Lyapunov multipliers in interior-point solver:");
        CONCLOG_PRINTLN("x="<<x<<", y="<<y<<", z="<<z);
        x=(1-sigma)*x+ApproximateVectorType(n,sigma/n);
        mu=dot(x,z)/n;
    }

    ApproximateVectorType yt=join(y,t);
    CONCLOG_PRINTLN("x="<<x<<" yt="<<yt<<" z="<<z);


    // Construct diagonal matrices
    ApproximateVectorType DE=ediv(x,z);
    CONCLOG_PRINTLN("DE="<<DE);

    // Construct the extended valuation GY=(gy-cu+te,cl-gy+te,y-bu+te,bl-y+te)
    ApproximateVectorType gye(2*(m+n),dp);
    //for(SizeType j=0; j!=n; ++j) { gxe[j]=gy[j]-C[j].upper_bound()+t; gye[n+j]=C[j].lower_bound()-gy[j]+t; }
    //for(SizeType i=0; i!=m; ++i) { gye[2*n+i]=y[i]-D[i].upper_bound()+t; gye[2*n+m+i]=D[i].lower_bound()-y[i]+t; }
    CONCLOG_PRINTLN("GE="<<gye);

    // Construct the extended matrix AE=(A -A I -I \\ e e 0 0)
    FloatDPApproximationMatrix AE(m+1,2*(m+n),dp);
    //for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { AE[i][j]=A[i][j]; AE[i][n+j]=-A[i][j]; } }
    //for(SizeType i=0; i!=m; ++i) { AE[i][2*n+i]=1; AE[i][2*n+m+i]=-1; }
    //for(SizeType k=0; k!=o; ++k) { AE[m][k]=1; }
    FloatDPApproximationMatrix AET=transpose(AE);

    // Construct the symmetric matrix and its inverse
    //FloatDPMatrix S(m+1,m+1); adat(S,AE,DE);
    //CONCLOG_PRINTLN("S="<<S);
    //S=FloatDPMatrix(m+1,m+1); simple_adat(S,AE,DE);
    //CONCLOG_PRINTLN("S="<<S);
    FloatDPApproximationMatrix S=feasibility_adat(Q,A,DiagonalMatrix<FloatDPApproximation>(DE));
    CONCLOG_PRINTLN("S="<<S);
    FloatDPApproximationMatrix Sinv=inverse(S);
    CONCLOG_PRINTLN("Sinv="<<Sinv);

    // FIXME: What if S is not invertible?

    // Construct the residuals
    ApproximateVectorType rx=esub(emul(x,z),mu*sigma);
    //RawFloatDPVector ryt=-prod(AE,x); ryt[m]+=1; // FIXME: Need hessian
    ApproximateVectorType ryt=-feasibility_mul(A,x); ryt[m]+=1; // FIXME: Need hessian
    ApproximateVectorType rz=gye+z;
    CONCLOG_PRINTLN("rx="<<rx<<" ryt="<<ryt<<" rz="<<rz);

    //RawFloatDPVector rr=prod(AE,ediv(RawFloatDPVector(rx-emul(x,rz)),z))-ryt;
    ApproximateVectorType rr=ryt + AE*ediv(ApproximateVectorType(rx-emul(x,rz)),z) - ryt;


    // Compute the differences
    ApproximateVectorType dyt=solve(S,rr);
    //RawFloatDPVector dz=-rz-prod(AET,dyt);
    ApproximateVectorType dz=-rz-feasibility_trmul(A,dyt);
    ApproximateVectorType dx=-ediv(ApproximateVectorType(rx+emul(x,dz)),z);
    CONCLOG_PRINTLN("dx="<<dx<<" dyt="<<dyt<<" dz="<<dz);

    ApproximateVectorType nx,ny,nyt,nz; FloatDPApproximation nt(dp);

    // Since we need to keep the point feasible, but the updates are linear
    // we need to validate feasibility directly rather than assuming the
    // linear update of y and z are good enough.
    Bool allpositive=false;
    FloatDPApproximation alpha=1/FloatDPApproximation(SCALE,dp);
    if(!egtr(emul(x,z) , GAMMA*mu/16)) {
        CONCLOG_PRINTLN("WARNING: x="<<x<<", z="<<z<< ", x.z="<<emul(x,z)<<"<"<<GAMMA*mu / 16);
        throw NearBoundaryOfFeasibleDomainException();
    }
    while(!allpositive) {
        alpha=alpha*SCALE;
        nx=x+alpha*dx;
        nyt=yt+alpha*dyt;
        ny=project(nyt,range(0,m));
        nt=nyt[m];
        //PrimalDualComplementaryInteriorPointOptimiser::compute_z(D,g,C,ny,nt,nz);
        allpositive = egtr(nx,0.0) && egtr(nz,0.0) && egtr(emul(nx,nz),GAMMA*mu);
    }
    CONCLOG_PRINTLN("alpha="<<alpha);
    CONCLOG_PRINTLN("nx="<<nx<<" nyt="<<nyt<<" nz="<<nz<<" nxz="<<emul(nx,nz));

    x=nx; y=project(nyt,range(0,m)); z=nz; t=nyt[m];
}


//------- PrimalDualComplementaryInteriorPointOptimiser -----------------------------------//

PrimalDualComplementaryInteriorPointOptimiser::PrimalDualComplementaryInteriorPointOptimiser() {
}

auto
PrimalDualComplementaryInteriorPointOptimiser::clone() const -> PrimalDualComplementaryInteriorPointOptimiser* {
    return new PrimalDualComplementaryInteriorPointOptimiser(*this);
}

OutputStream& operator<<(OutputStream& os, PrimalDualComplementaryInteriorPointOptimiser const& opt) {
    return os << "PrimalDualComplementaryInteriorPointOptimiser()";
}









// See Hande Y. Benson, David F. Shanno, And Robert J. Vanderbei,
// "Interior-point methods for nonconvex nonlinear programming: Jamming and comparative numerical testing"
// For some of the terminology used


Void
PrimalDualComplementaryInteriorPointOptimiser::_minimisation_step(const ApproximateOptimisationProblem& p, StepDataBase& d) const
{
    this->minimisation_step(p, dynamic_cast<StepData&>(d));
}

// min f(x) | x\in D & w\in C | g(x) = w
// Lagrange multipliers y d(g(x)-w); z dx
Void
PrimalDualComplementaryInteriorPointOptimiser::minimisation_step(const ApproximateOptimisationProblem& p, StepData& d) const
{
    auto& x=d.x; auto& y=d.y; auto& z=d.z; auto& mu=d.mu;
    this->minimisation_step(p, x,y,z,mu);
}

auto
PrimalDualComplementaryInteriorPointOptimiser::minimisation_update(
    const ApproximateOptimisationProblem& p,
    ApproximateVectorType& x, ApproximateVectorType& y, ApproximateVectorType& z,
    Approximation<FLT>& mu) const -> PrimalDualComplementaryData<Approximation<FLT>>
{
    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;

    const SizeType m=g.result_size();
    const SizeType n=g.argument_size();

    ARIADNE_DEBUG_PRECONDITION(decide(contains(D,x)));
    ARIADNE_DEBUG_PRECONDITION(decide(mu>0));

    ApproximateVectorType w=g(x);
    ARIADNE_DEBUG_PRECONDITION(decide(contains(C,w)));

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("x="<<x);
    CONCLOG_PRINTLN("w="<<w);
    CONCLOG_PRINTLN("y="<<y);
    CONCLOG_PRINTLN("z="<<z);
    CONCLOG_PRINTLN("mu="<<mu);


    Differential<Approximation<FLT>> ddfx=f.evaluate(Differential<Approximation<FLT>>::variables(2,x));
    Vector<Differential<Approximation<FLT>>> ddgx=g.evaluate(Differential<Approximation<FLT>>::variables(2,x));

    // G is the constraint value vector
    Approximation<FLT> fx = ddfx.value();
    Vector<Approximation<FLT>> gx = ddgx.value();
    CONCLOG_PRINTLN("f(x)="<<fx);
    CONCLOG_PRINTLN("g(x)="<<gx);
    CONCLOG_PRINTLN_AT(1,"g(x)-w="<<(gx-w));

    // A, B are the derivative matrices aij=dgi/dxj
    // HACK: Need to explicitly set size of Jacobian if g or h have result_size of zero
    Vector<Approximation<FLT>> c = transpose(ddfx.gradient());
    CONCLOG_PRINTLN_AT(1,"c=df(x)="<<c);
    Matrix<Approximation<FLT>> A = ddgx.jacobian();
    if(m==0) { A=Matrix<Approximation<FLT>>(m,n,dp); }
    CONCLOG_PRINTLN("A=Dg(x)="<<A);

    // Q is the Hessian matrix Q[i1,i2] = df/dx[i1]dx[i2] + Sum_[j]y[j]*dg[j]/dx[i1]dx[i2]
    SymmetricMatrix<Approximation<FLT>> Q = ddfx.hessian();
    for(SizeType j=0; j!=m; ++j) { Q += y[j] * ddgx[j].hessian(); }
    CONCLOG_PRINTLN("Q="<<Q);

    // Determines the weighting to give to the relaxation parameter mu
    // for equality constraints relative to other constraints
    // Currently unused
    //    static const double EQUALITY_RELAXATION_MULTIPLIER = 1.0;

    // Set-up equations
    //   (w-cl)*(cu-w)*y + (cl+cu-2*w) * mu = 0
    //   Df(x) + sum y_j Dg_j(x) + z = 0
    //   g(x) - w = 0
    //   (x-dl)*(du-x)*z + (dl+du-2*w) * mu = 0

    auto Cl = lower_bounds(C);
    auto Cu = upper_bounds(C);
    auto Dl = lower_bounds(D);
    auto Du = upper_bounds(D);

    // Compute the residuals and contributions from slack in x and w
    //   rx[i] = df/dx[i] + Sum[j] dg[j]/dx[i] * y[j] + z[i]
    Vector<Approximation<FLT>> rw = gx - w;
    Vector<Approximation<FLT>> rx = c + transpose(A) * y + z;
    Vector<Approximation<FLT>> ry = emul(w-Cl,Cu-w,y) + (Cl+Cu-2*w)*mu;
    Vector<Approximation<FLT>> rz = emul(x-Dl,Du-x,z) + (Dl+Du-2*x)*mu;

    DiagonalMatrix<Approximation<FLT>> W(emul(w-Cl,Cu-w));
    DiagonalMatrix<Approximation<FLT>> Y(esub(emul(Cl+Cu-2*w,y),2*mu));
    DiagonalMatrix<Approximation<FLT>> X(emul(x-Dl,Du-x));
    DiagonalMatrix<Approximation<FLT>> Z(esub(emul(Dl+Du-2*x,z),2*mu));

    // FIXME: Make this more consistent
    for (SizeType j=0; j!=m; ++j) {
        if (decide(Cl[j]==-inf)) { ry[j]=(Cu[j]-w[j])*y[j]-mu; W._at(j)=Cu[j]-w[j]; Y._at(j)=y[j]; }
        else if (decide(Cu[j]==+inf)) { ry[j]=(w[j]-Cl[j])*y[j]-mu; W._at(j)=w[j]-Cl[j]; Y._at(j)=y[j]; }
    }

    CONCLOG_PRINTLN("rw="<<rw);
    CONCLOG_PRINTLN("rx="<<rx);
    CONCLOG_PRINTLN("ry="<<ry);
    CONCLOG_PRINTLN("rz="<<rz);

    SlackPrimalDualComplementaryMatrix<Approximation<FLT>> S(Q,A,W,X,Y,Z);
    SlackPrimalDualComplementaryData<Approximation<FLT>> r(rw,rx,ry,rz);

    SlackPrimalDualComplementaryData<Approximation<FLT>> dwr=S.solve(r);

    PrimalDualComplementaryData<Approximation<FLT>> dr(dwr.x,dwr.y,dwr.z);

    return dr;
}

Void
PrimalDualComplementaryInteriorPointOptimiser::minimisation_step(
    const ApproximateOptimisationProblem& p,
    ApproximateVectorType& x, ApproximateVectorType& y, ApproximateVectorType& z,
    Approximation<FLT>& mu) const
{
    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;

    ARIADNE_DEBUG_PRECONDITION(x.size()==f.argument_size());
    ARIADNE_DEBUG_PRECONDITION(z.size()==x.size());

    ARIADNE_DEBUG_PRECONDITION(decide(contains(D,x)));
    ARIADNE_DEBUG_PRECONDITION(decide(contains(C,g(x))));
    ARIADNE_DEBUG_PRECONDITION(decide(mu>0));

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("w=g(x)="<<g(x)<<", x="<<x<<", y="<<y<<", z="<<z<<", mu="<<mu);

    PrimalDualComplementaryData<Approximation<FLT>> dr = this->minimisation_update(p, x,y,z, mu);

    static const Approximation<FLT> ALPHA_SCALE_FACTOR = 0.75_approx;
    static const Approximation<FLT> MINIMUM_ALPHA = 1e-16_approx;

    // Compute distance to move variables preserving feasibility
    // FIXME: Current implementation might fail due to getting too close to boundary!
    ApproximateVectorType newx(x.size(),dp);
    ApproximateVectorType neww(y.size(),dp);
    Approximation<FLT> alpha = 1.0_approx;
    Bool success = false;
    do {
        newx = x - alpha * dr.x;
        neww = g(newx);
        if (probably(contains(D,newx)) && probably(contains(C,neww))) { success = true; }
        else { alpha *= ALPHA_SCALE_FACTOR; }
        if (probably(alpha<MINIMUM_ALPHA)) { throw NearBoundaryOfFeasibleDomainException(); }
    } while (!success);
    CONCLOG_PRINTLN("alpha="<<alpha);

    ApproximateVectorType newz = z - alpha * dr.z;
    ApproximateVectorType newy = y - alpha * dr.y;

    CONCLOG_PRINTLN("neww="<<neww<<", newx="<<newx<<", newy="<<newy<<", newz="<<newz);

    x=newx; y=newy; z=newz;
}


Void
PrimalDualComplementaryInteriorPointOptimiser::feasibility_step(
    const ApproximateFeasibilityProblem& p,
    ApproximateVectorType& x, ApproximateVectorType& y, ApproximateVectorType& z) const
{
    CONCLOG_SCOPE_CREATE;

    EffectiveScalarMultivariateFunction f(p.g.argument_size());
    ApproximateOptimisationProblem optp(f,p.D,p.g,p.C);

    Approximation<FLT> mu=one;
    this->minimisation_step(optp, x,y,z, mu);
}



auto
PrimalDualComplementaryInteriorPointOptimiser::initial_step_data(
    const ApproximateFeasibilityProblem& p) const
        -> StepData*
{
    return dynamic_cast<StepData*>(this->InteriorPointOptimiserBase::initial_step_data(p));
}

auto
PrimalDualComplementaryInteriorPointOptimiser::initial_step_data_hotstarted(
    const ApproximateFeasibilityProblem& p, const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> StepData*
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("PrimalDualComplementaryInteriorPointOptimiser::initial_step_data_hotstarted(p,x0,y0)");
    CONCLOG_PRINTLN("p="<<p);
    CONCLOG_PRINTLN("x0="<<x0<<", y0="<<y0);
    Approximation<FLT> mu=one;
    CONCLOG_PRINTLN("g(x0)="<<p.g(x0));
    ARIADNE_PRECONDITION(decide(contains(p.C,p.g(x0))));
    ApproximateVectorType z0=this->compute_z(p,x0,mu);
    CONCLOG_PRINTLN("z0="<<z0);
    return new StepData(x0,y0,z0,mu);
}

auto
PrimalDualComplementaryInteriorPointOptimiser::_initial_step_data_hotstarted(
    const ApproximateFeasibilityProblem& p, const ApproximateVectorType& x, const ApproximateVectorType& y) const
        -> StepDataBase*
{
    return this->initial_step_data_hotstarted(p,x,y);
}






//------- SlackPrimalDualComplementaryInteriorPointOptimiser -----------------------------------//

SlackPrimalDualComplementaryInteriorPointOptimiser::SlackPrimalDualComplementaryInteriorPointOptimiser() {
}

auto
SlackPrimalDualComplementaryInteriorPointOptimiser::clone() const -> SlackPrimalDualComplementaryInteriorPointOptimiser* {
    return new SlackPrimalDualComplementaryInteriorPointOptimiser(*this);
}

OutputStream& operator<<(OutputStream& os, SlackPrimalDualComplementaryInteriorPointOptimiser const& opt) {
    return os << "SlackPrimalDualComplementaryInteriorPointOptimiser()";
}









// See Hande Y. Benson, David F. Shanno, And Robert J. Vanderbei,
// "Interior-point methods for nonconvex nonlinear programming: Jamming and comparative numerical testing"
// For some of the terminology used


Void
SlackPrimalDualComplementaryInteriorPointOptimiser::_minimisation_step(const ApproximateOptimisationProblem& p, StepDataBase& d) const
{
    this->minimisation_step(p, dynamic_cast<StepData&>(d));
}

// min f(x) | x\in D & w\in C | g(x) = w
// Lagrange multipliers y d(g(x)-w); z dx
auto
SlackPrimalDualComplementaryInteriorPointOptimiser::minimisation_step(const ApproximateOptimisationProblem& p, StepData& d) const -> Approximation<FLT>
{
    auto& w=d.w; auto& x=d.x; auto& y=d.y; auto& z=d.z; auto& mu=d.mu;
    return this->minimisation_step(p, w,x,y,z,mu);
}

auto
SlackPrimalDualComplementaryInteriorPointOptimiser::minimisation_update(
    const ApproximateOptimisationProblem& p,
    ApproximateVectorType& w, ApproximateVectorType& x, ApproximateVectorType& y, ApproximateVectorType& z,
    Approximation<FLT>& mu) const -> SlackPrimalDualComplementaryData<Approximation<FLT>>
{
    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;

    const SizeType m=g.result_size();
    const SizeType n=g.argument_size();

    ARIADNE_DEBUG_PRECONDITION(w.size()==y.size());
    ARIADNE_DEBUG_PRECONDITION(x.size()==z.size());
    ARIADNE_DEBUG_PRECONDITION(w.size()==m);
    ARIADNE_DEBUG_PRECONDITION(x.size()==n);

    ARIADNE_DEBUG_PRECONDITION(decide(contains(D,x)));
    ARIADNE_DEBUG_PRECONDITION(decide(contains(C,w)));
    ARIADNE_DEBUG_PRECONDITION(decide(mu>0));

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("x="<<x);
    CONCLOG_PRINTLN("w="<<w);
    CONCLOG_PRINTLN("y="<<y);
    CONCLOG_PRINTLN("z="<<z);
    CONCLOG_PRINTLN("mu="<<mu);


    Differential<Approximation<FLT>> ddfx=f.evaluate(Differential<Approximation<FLT>>::variables(2,x));
    Vector<Differential<Approximation<FLT>>> ddgx=g.evaluate(Differential<Approximation<FLT>>::variables(2,x));

    // G is the constraint value vector
    Approximation<FLT> fx = ddfx.value();
    Vector<Approximation<FLT>> gx = ddgx.value();
    CONCLOG_PRINTLN("f(x)="<<fx);
    CONCLOG_PRINTLN("g(x)="<<gx);
    CONCLOG_PRINTLN_AT(1,"g(x)-w="<<(gx-w));

    // A, B are the derivative matrices aij=dgi/dxj
    // HACK: Need to explicitly set size of Jacobian if g or h have result_size of zero
    Vector<Approximation<FLT>> c = transpose(ddfx.gradient());
    CONCLOG_PRINTLN_AT(1,"c=df(x)="<<c);
    Matrix<Approximation<FLT>> A = ddgx.jacobian();
    if(m==0) { A=Matrix<Approximation<FLT>>(m,n,dp); }
    CONCLOG_PRINTLN("A=Dg(x)="<<A);

    // Q is the Hessian matrix Q[i1,i2] = df/dx[i1]dx[i2] + Sum_[j]y[j]*dg[j]/dx[i1]dx[i2]
    SymmetricMatrix<Approximation<FLT>> Q = ddfx.hessian();
    for(SizeType j=0; j!=m; ++j) { Q += y[j] * ddgx[j].hessian(); }
    CONCLOG_PRINTLN("Q="<<Q);

    // Determines the weighting to give to the relaxation parameter mu
    // for equality constraints relative to other constraints
    // Currently unused
    //    static const double EQUALITY_RELAXATION_MULTIPLIER = 1.0;

    // Set-up equations
    //   (w-cl)*(cu-w)*y + (cl+cu-2*w) * mu = 0
    //   Df(x) + sum y_j Dg_j(x) + z = 0
    //   g(x) - w = 0
    //   (x-dl)*(du-x)*z + (dl+du-2*w) * mu = 0

    auto Cl = lower_bounds(C);
    auto Cu = upper_bounds(C);
    auto Dl = lower_bounds(D);
    auto Du = upper_bounds(D);

    // Compute the residuals and contributions from slack in x and w
    //   rx[i] = df/dx[i] + Sum[j] dg[j]/dx[i] * y[j] + z[i]
    Vector<Approximation<FLT>> rw = gx - w;
    Vector<Approximation<FLT>> rx = c + transpose(A) * y + z;

    Vector<Approximation<FLT>> ry = compute_complementarity_residuals(C,w,y,mu);
    Vector<Approximation<FLT>> rz = compute_complementarity_residuals(D,x,z,mu);

    auto YW=compute_linear_equations(C,w,y,mu);
    auto ZX=compute_linear_equations(D,x,z,mu);

    DiagonalMatrix<Approximation<FLT>>& W=std::get<1>(YW);
    DiagonalMatrix<Approximation<FLT>>& Y=std::get<0>(YW);
    DiagonalMatrix<Approximation<FLT>>& X=std::get<1>(ZX);
    DiagonalMatrix<Approximation<FLT>>& Z=std::get<0>(ZX);

    CONCLOG_PRINTLN("rw="<<rw);
    CONCLOG_PRINTLN("rx="<<rx);
    CONCLOG_PRINTLN("ry="<<ry);
    CONCLOG_PRINTLN("rz="<<rz);

    SlackPrimalDualComplementaryMatrix<Approximation<FLT>> S(Q,A,W,X,Y,Z);
    SlackPrimalDualComplementaryData<Approximation<FLT>> r(rw,rx,ry,rz);

    CONCLOG_PRINTLN("S="<<pretty(S.assemble()));
    SlackPrimalDualComplementaryData<Approximation<FLT>> d=S.solve(r);
    CONCLOG_PRINTLN("dw="<<d.w);
    CONCLOG_PRINTLN("dx="<<d.x);
    CONCLOG_PRINTLN("dy="<<d.y);
    CONCLOG_PRINTLN("dz="<<d.z);
    return d;
}

auto
SlackPrimalDualComplementaryInteriorPointOptimiser::minimisation_step_size(
    const ApproximateOptimisationProblem& p,
    const ApproximateVectorType& w, const ApproximateVectorType& x, const ApproximateVectorType& y, const ApproximateVectorType& z,
    ApproximateVectorType& dw, ApproximateVectorType& dx, ApproximateVectorType& dy, ApproximateVectorType& dz,
    const Approximation<FLT>& mu) const -> Approximation<FLT>
{
    Approximation<FLT> alpha = 1.0_approx;
    if constexpr (true) {
        static const Approximation<FLT> SIGMA = 1.0_approx/1024;
        Approximation<FLT> e=mu*SIGMA;
        CONCLOG_PRINTLN("e="<<e);
        for (SizeType i=0; i!=x.size(); ++i) {
            if (probably(dx[i]<0)) {
                if (decide(x[i]<p.D[i].upper_bound()-e)) { alpha=min(alpha,(x[i]-p.D[i].upper_bound()+e)/dx[i]); }
                else { dx[i]=0; }
            } else if (probably(dx[i]>0)) {
                if (decide(x[i]>p.D[i].lower_bound()+e)) { alpha=min(alpha,(x[i]-p.D[i].lower_bound()-e)/dx[i]); }
                else { dx[i]=0; }
            }
        }
        for (SizeType j=0; j!=w.size(); ++j) {
            if (probably(dw[j]<0)) {
                if (decide(w[j]<p.C[j].upper_bound()-e)) { alpha=min(alpha,(w[j]-p.C[j].upper_bound()+e)/dw[j]); }
                else { dw[j]=0; }
            } else if (probably(dw[j]>0)) {
                if (decide(w[j]>p.C[j].lower_bound()+e)) { alpha=min(alpha,(w[j]-p.C[j].lower_bound()-e)/dw[j]); }
                else { dw[j]=0; }
            }
        }
        CONCLOG_PRINTLN("dw="<<dw<<", dx="<<dx<<", dy="<<dy<<", dz="<<dz);
        CONCLOG_PRINTLN("dw="<<cast_exact(dw)<<", dx="<<cast_exact(dx)<<", dy="<<cast_exact(dy)<<", dz="<<cast_exact(dz));
    }
    if constexpr (false) {
        static const Approximation<FLT> ALPHA_SCALE_FACTOR = 0.75_approx;
        static const Approximation<FLT> MINIMUM_ALPHA = 1e-16_approx;

        // Compute distance to move variables preserving feasibility
        // FIXME: Current implementation might fail due to getting too close to boundary!
        ApproximateVectorType newx(x.size(),dp);
        ApproximateVectorType neww(w.size(),dp);
        Approximation<FLT> alpha = 1.0_approx;
        Bool success = false;
        do {
            CONCLOG_PRINTLN("al="<<alpha<<"\nx="<<x<<", dx="<<dx<<"\ny="<<y<<", dy="<<dy)
            newx = x - alpha * dx;
            neww = w - alpha * dw;
            if (probably(contains(p.D,newx)) && probably(contains(p.C,neww))) { success = true; }
            else { alpha *= ALPHA_SCALE_FACTOR; }
            if (probably(alpha<MINIMUM_ALPHA)) { throw NearBoundaryOfFeasibleDomainException(); }
        } while (!success);
    }
    CONCLOG_PRINTLN("alpha="<<cast_exact(alpha));

    return alpha;
}

auto
SlackPrimalDualComplementaryInteriorPointOptimiser::minimisation_step(
    const ApproximateOptimisationProblem& p,
    ApproximateVectorType& w, ApproximateVectorType& x, ApproximateVectorType& y, ApproximateVectorType& z,
    Approximation<FLT>& mu) const -> Approximation<FLT>
{
    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;

    ARIADNE_DEBUG_PRECONDITION(w.size()==g.result_size());
    ARIADNE_DEBUG_PRECONDITION(x.size()==f.argument_size());
    ARIADNE_DEBUG_PRECONDITION(y.size()==w.size());
    ARIADNE_DEBUG_PRECONDITION(z.size()==x.size());

    ARIADNE_DEBUG_PRECONDITION(decide(contains(D,x)));
    ARIADNE_DEBUG_PRECONDITION(decide(contains(C,w)));
    ARIADNE_DEBUG_PRECONDITION(decide(mu>0));

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("SlackPrimalDualComplementaryInteriorPointOptimiser::minimisation_step(p, w,x,y,z,my)");
    CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y<<", z="<<z<<", mu="<<mu);

    SlackPrimalDualComplementaryData<Approximation<FLT>> dr = this->minimisation_update(p, w,x,y,z, mu);
    CONCLOG_PRINTLN("dw="<<dr.w<<", dx="<<dr.x<<", dy="<<dr.y<<", dz="<<dr.z);

    ApproximateNumericType alpha = minimisation_step_size(p, w,x,y,z, dr.w,dr.x,dr.y,dr.z, mu);

    ApproximateVectorType newx = x - alpha * dr.x;
    ApproximateVectorType neww = w - alpha * dr.w;
    ApproximateVectorType newz = z - alpha * dr.z;
    ApproximateVectorType newy = y - alpha * dr.y;

    CONCLOG_PRINTLN("neww="<<neww<<", newx="<<newx<<", newy="<<newy<<", newz="<<newz);

    x=newx; w=neww; y=newy; z=newz;

    return alpha;
}


Void
SlackPrimalDualComplementaryInteriorPointOptimiser::feasibility_step(
    const ApproximateFeasibilityProblem& p,
    ApproximateVectorType& w, ApproximateVectorType& x, ApproximateVectorType& y, ApproximateVectorType& z) const
{
    CONCLOG_SCOPE_CREATE;

    EffectiveScalarMultivariateFunction f(p.g.argument_size());
    ApproximateOptimisationProblem optp(f,p.D,p.g,p.C);

    Approximation<FLT> mu=one;
    this->minimisation_step(optp, w,x,y,z, mu);
}



auto
SlackPrimalDualComplementaryInteriorPointOptimiser::initial_step_data(const ApproximateFeasibilityProblem& p) const -> StepData*
{
    return dynamic_cast<StepData*>(this->InteriorPointOptimiserBase::initial_step_data(p));
}


auto
SlackPrimalDualComplementaryInteriorPointOptimiser::initial_step_data_hotstarted(
    const ApproximateFeasibilityProblem& p, const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> StepData*
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("SlackPrimalDualComplementaryInteriorPointOptimiser::initial_step_data_hotstarted(p,x0,y0)");
    CONCLOG_PRINTLN("p="<<p);
    CONCLOG_PRINTLN("x0="<<x0<<", y0="<<y0);
    Approximation<FLT> mu0=one;
    ApproximateVectorType w0=this->compute_w(p,x0,y0,mu0);
    ApproximateVectorType z0=this->compute_z(p,x0,mu0);
    CONCLOG_PRINTLN("w0="<<w0<<", z0="<<z0<<", mu0="<<mu0);
    return new StepData(w0,x0,y0,z0,mu0);
}

auto
SlackPrimalDualComplementaryInteriorPointOptimiser::_initial_step_data_hotstarted(
    const ApproximateFeasibilityProblem& p, const ApproximateVectorType& x, const ApproximateVectorType& y) const
        -> StepDataBase*
{
    return this->initial_step_data_hotstarted(p,x,y);
}






//------- SlackPrimalSplitDualComplementaryInteriorPointOptimiser -------------------------//


namespace {
template<class X> inline Vector<X> operator-(Vector<X> v1, X const& s2) { return esub(v1,s2); }
template<class X> inline Vector<X> operator-(X const& s1, Vector<X> v2) { return esub(s1,v2); }
template<class X> inline Vector<X> operator*(Vector<X> v1, Vector<X> const& v2) { return emul(v1,v2); }
template<class X> inline Vector<X> operator/(Vector<X> v1, Vector<X> const& v2) { return ediv(v1,v2); }
template<class X> inline Vector<X> operator/(X const& s1, Vector<X> v2) { return ediv(s1,v2); }
}

OutputStream& operator<<(OutputStream& os, SlackPrimalSplitDualComplementaryInteriorPointOptimiser const& opt) {
    return os << "SlackPrimalSplitDualComplementaryInteriorPointOptimiser1()";
}

auto
SlackPrimalSplitDualComplementaryInteriorPointOptimiser::initial_step_data(const ApproximateFeasibilityProblem& p) const -> StepData*
{
    return dynamic_cast<StepData*>(this->InteriorPointOptimiserBase::initial_step_data(p));
}

auto
SlackPrimalSplitDualComplementaryInteriorPointOptimiser::initial_step_data_hotstarted(
    const ApproximateFeasibilityProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const -> StepData*
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("SlackPrimalSplitDualComplementaryInteriorPointOptimiser::initial_step_data_hotstarted(p,x0,y0)");
    CONCLOG_PRINTLN("p="<<p);
    CONCLOG_PRINTLN("x0="<<x0<<", y0="<<y0);
    auto& D=p.D; auto& C=p.C;

    ApproximateNumericType mu = one;

    auto x = x0;
    CONCLOG_PRINTLN("x="<<x);

    auto w = this->compute_w(p,x0,y0,mu);
    CONCLOG_PRINTLN("w="<<w);

    auto yl = mu/(w-lower_bounds(C));
    auto yu = mu/(upper_bounds(C)-w);
    auto zl = mu/(x-lower_bounds(D));
    auto zu = mu/(upper_bounds(D)-x);
    CONCLOG_PRINTLN("yl="<<yl<<", yu="<<yu<<", zl="<<zl<<", zu="<<zu);

    CONCLOG_PRINTLN("d="<<StepData(w,x,yl,yu,zl,zu, mu));

    return new StepData(w,x,yl,yu,zl,zu, mu);
}


auto
SlackPrimalSplitDualComplementaryInteriorPointOptimiser::clone() const -> SlackPrimalSplitDualComplementaryInteriorPointOptimiser*
{
    return new SlackPrimalSplitDualComplementaryInteriorPointOptimiser(*this);
}

auto
SlackPrimalSplitDualComplementaryInteriorPointOptimiser::_initial_step_data_hotstarted(
    const ApproximateFeasibilityProblem& p, const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> StepDataBase*
{
    return this->initial_step_data_hotstarted(p,x0,y0);
}

Void
SlackPrimalSplitDualComplementaryInteriorPointOptimiser::_minimisation_step(
    const ApproximateOptimisationProblem& p,
    StepDataBase& d) const
{
    this->minimisation_step(p,dynamic_cast<StepData&>(d));
}

auto
SlackPrimalSplitDualComplementaryInteriorPointOptimiser::minimisation_step(
    const ApproximateOptimisationProblem& p,
    StepData& d) const -> ApproximateNumericType
{
    using XA = ApproximateNumberType;

    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("SlackPrimalSplitDualComplementaryInteriorPointOptimiser::minimisation_step(p,d)");
    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;
    CONCLOG_PRINTLN("p="<<p);
    CONCLOG_PRINTLN("d="<<d);
    auto& w=d.w; auto& x=d.x; auto& yl=d.yl; auto& yu=d.yu; auto& zl=d.zl; auto& zu=d.zu; auto& mu=d.mu;

    auto Cl = lower_bounds(C);
    auto Cu = upper_bounds(C);
    auto Dl = lower_bounds(D);
    auto Du = upper_bounds(D);

    auto y=yu-yl;
    auto z=zu-zl;

    auto ddfx = differential(f,x,2);
    auto ddgx = differential(g,x,2);

    auto c=ddfx.gradient();
    auto A=ddgx.jacobian();
    auto Q=ddfx.hessian();
    for (SizeType i=0; i!=y.size(); ++i) {
        Q+=y[i]*ddgx[i].hessian();
    }

// Equations (w-Cl)*yl=mu, (Cu-w)*yu=mu
// y=yu-yl
// (w-Cl)*dyl + yl * dw = -((w-Cl)*yl-mu) =: -rwl
// (Cu-w)*dyu - yu * dw = -((Cu-w)*yu-mu) =: -rwu
// (w-Cl)*(Cu-w)*dyl + (Cu-w)*yl * dw = - (Cu-w)*((w-Cl)*yl-mu)
// (w-Cl)*(Cu-w)*dyu - (w-Cl)*yu * dw = - (w-Cl)*((Cu-w)*yu-mu)
// (w-Cl)*(Cu-w)*dy  - ((w-Cl)*yu+(Cu-w)*yl) * dw = - ((w-Cl)*rwu-(Cu-w)*rwl) =: -rw
//      = -((w-Cl)*((Cu-w)*yu-mu)-(Cu-w)*((w-Cl)*yl-mu))
//      = -((w-Cl)*(Cu-w)*y-(2*w-Cl-Cu)*mu) := rw
//
// If Cu=+inf, yu=0, so (w-Cl)*(-dyl) - yl*dw = - (mu-(w-Cl)*yl)
// If Cl=-inf, yl=0, so (Cu-w)*(+dyu) - yu*dw = + (mu-(Cu-w)*yu)
// If Cl=Cu=0, -w^2*dy - (w*yu - w*yl) * dw = w*(mu + w*yu)+w*(mu-w*yl)
//    -w*dy - (yu-yl) * dw = 2*mu + w * (yu-yl)

    auto ryl = (w-Cl)*yl-mu;
    auto ryu = (Cu-w)*yu-mu;
    auto rzl = (x-Dl)*zl-mu;
    auto rzu = (Du-x)*zu-mu;

    auto rw = g(x)-w;
    auto rx = transpose(c)+transpose(A)*y+z;
    auto ry = (w-Cl)*ryu-(Cu-w)*ryl;
    auto rz = (x-Dl)*rzu-(Du-x)*rzl;

    auto W = DiagonalMatrix<XA>((w-Cl)*(Cu-w));
    auto X = DiagonalMatrix<XA>((x-Dl)*(Du-x));
    auto Y = DiagonalMatrix<XA>(-((w-Cl)*yu+(Cu-w)*yl));
    auto Z = DiagonalMatrix<XA>(-((x-Dl)*zu+(Du-x)*zl));

    CONCLOG_PRINTLN("ryl="<<ryl<<", ryu="<<ryu<<", rzl="<<rzl<<", rzu="<<rzu);
    CONCLOG_PRINTLN("rw="<<rw<<", rx="<<rx<<", ry="<<ry<<", rz="<<rz);

    SlackPrimalDualComplementaryMatrix<XA> S(Q,A,W,X,Y,Z);
    SlackPrimalDualComplementaryData<XA> r(rw,rx,ry,rz);
    SlackPrimalDualComplementaryData<XA> dd=-S.solve(r);

    auto& dw=dd.w; auto& dx=dd.x; auto& dy=dd.y; auto& dz=dd.z;


    // Use relations for dw,dy (similary for dx,dz)
    // Note must have dy = dyu-dyl
    // (cu-w) * dyu - yu * dw = (mu - (cu-w)*yu)
    // (w-cl) * dyl + yl * dw = (mu - (w-cl)*yl)

    auto dyl = - (ryl+yl*dw) / (w-Cl);
    auto dyu = - (ryu-yu*dw) / (Cu-w);
    auto dzl = - (rzl+zl*dx) / (x-Dl);
    auto dzu = - (rzu-zu*dx) / (Du-x);

    CONCLOG_PRINTLN("dw="<<dw<<", dx="<<dx<<", dy="<<dy<<", dz="<<dz);
    CONCLOG_PRINTLN("dyl="<<dyl<<", dyu="<<dyu<<", dzl="<<dzl<<", dzu="<<dzu)
    CONCLOG_PRINTLN("dyu-dyl="<<dyu-dyl<<", dyu-dyl-dy="<<(dyu-dyl-dy)<<", norm(dyu-dyl-dy)="<<norm(dyu-dyl-dy));
    CONCLOG_PRINTLN("dzu-dzl="<<dzu-dzl<<", dzu-dzl-dz="<<(dzu-dzl-dz)<<", norm(dzu-dzl-dz)="<<norm(dzu-dzl-dz));

    ARIADNE_DEBUG_ASSERT(decide(norm(dyu-dyl-dy)<1e-8));
    ARIADNE_DEBUG_ASSERT(decide(norm(dzu-dzl-dz)<1e-8));

    XA alpha=one;

    static const ExactDouble GAMMA(1.0/1024);

    auto max_vector_step = [&](Vector<XA> const& v, Vector<XA> const& dv) {
        XA m(1,v.zero_element().precision()); XA e = mu*GAMMA; for (SizeType i=0; i!=v.size(); ++i) { if (decide(dv[i]>0)) { m=min((v[i]-e)/dv[i],m); } } return m; };

    alpha = min( one, min( min(max_vector_step(yl,dyl),max_vector_step(yu,dyu)), min(max_vector_step(zl,dzl),max_vector_step(zu,dzu)) ) );

    CONCLOG_PRINTLN("alpha="<<alpha);

    w += alpha * dw;
    x += alpha * dx;
    yl += alpha * dyl;
    yu += alpha * dyu;
    zl += alpha * dzl;
    zu += alpha * dzu;

    CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", yl="<<yl<<", yu="<<yu<<", zl="<<zl<<", zu="<<zu);

    rx=transpose(gradient(f,x))+transpose(jacobian(g,x))*y+z;
    CONCLOG_PRINTLN("rwl="<<(w-Cl)*yl-mu<<", rwu="<<(Cu-w)*yu-mu<<", rx=df(x)+y*dg(x)+z="<<rx<<", ry=g(x)-w="<<g(x)-w<<", rzl="<<(x-Dl)*zl-mu<<", rzu="<<(Du-x)*zu-mu);

    return alpha;
}


//------- KarushKuhnTuckerOptimiser -----------------------------------//

auto
KarushKuhnTuckerOptimiser::clone() const -> KarushKuhnTuckerOptimiser*
{
    return new KarushKuhnTuckerOptimiser(*this);
}


auto
KarushKuhnTuckerOptimiser::minimise(ValidatedOptimisationProblem p) const -> ValidatedVector
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("KarushKuhnTuckerOptimiser::minimise(ValidatedOptimisationProblem p) const -> ValidatedVector");

    ARIADNE_PRECONDITION(p.D.has_nonempty_interior());
    ARIADNE_PRECONDITION(p.C.has_nonempty_interior());
    ApproximateVectorType x0=p.D.midpoint();
    ApproximateVectorType y0(p.C.dimension(),zero);
    return this->minimise_hotstarted(p,x0,y0).primal();
}

auto
KarushKuhnTuckerOptimiser::feasible(ValidatedFeasibilityProblem p) const -> ValidatedKleenean
{
    CONCLOG_SCOPE_CREATE;

    ARIADNE_PRECONDITION(p.D.has_nonempty_interior());
    ARIADNE_PRECONDITION(p.C.has_nonempty_interior());

    ApproximateVectorType x0=p.D.midpoint();
    ApproximateVectorType y0(p.C.dimension(),zero);

    return this->feasible_hotstarted(p,x0,y0).is_feasible();
}

auto
KarushKuhnTuckerOptimiser::minimise_hotstarted(
    const ValidatedOptimisationProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> ValuePrimalDualData<ValidatedNumber>
{
    UpperBoxType B=p.D;

    auto r = this->splitting_minimise_hotstarted(p,B,x0,y0);

    auto ro0 = std::get<OptimalityCertificate>(r)[0];
    return ValuePrimalDualData<ValidatedNumber>(ro0.value(),ro0.primal(),ro0.dual());
}


auto
KarushKuhnTuckerOptimiser::nonsplitting_minimise_hotstarted(
    const ValidatedOptimisationProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> ValuePrimalDualData<ValidatedNumber>
{
    ARIADNE_NOT_IMPLEMENTED;

    PrimalDualComplementaryInteriorPointOptimiser approximate_optimiser;

    auto r = approximate_optimiser.minimise_hotstarted(ApproximateOptimisationProblem(p),x0,y0);
    // FIXME: Find a better way of making the answer exact
    const Scalar<ApproximateNumber>& av=r.value();
    const Vector<ApproximateNumber>& ax=r.primal();
    const Vector<ApproximateNumber>& ay=r.dual();

    Scalar<ValidatedNumber> v=to_generic(cast_exact(Scalar<Approximation<FLT>>(av,pr)));
    Vector<ValidatedNumber> x=to_generic(cast_exact(Vector<Approximation<FLT>>(ax,pr)));
    Vector<ValidatedNumber> y=to_generic(cast_exact(Vector<Approximation<FLT>>(ay,pr)));

    ARIADNE_WARN("KarushKuhnTuckerOptimiser::nonsplitting_minimise_hotstarted(p,x0a,y0a): Result may not be reliable.");
    return ValuePrimalDualData<ValidatedNumber>(v,x,y);
}


auto
KarushKuhnTuckerOptimiser::splitting_minimise_hotstarted(
    const ValidatedOptimisationProblem& p,
    UpperBoxType B, ApproximateVectorType x, ApproximateVectorType y) const
        -> Variant<OptimalityCertificate,InfeasibilityCertificate>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;

    PrimalDualComplementaryInteriorPointOptimiser approximate_optimiser;

    // Extra variables for use of optimization solver
    Approximation<FLT> mu = one;
    ApproximateVectorType z=approximate_optimiser.compute_z(p,x,mu);

    for(SizeType i=0; i!=18; ++i) {
        CONCLOG_PRINTLN_AT(1,"x="<<x<<", y="<<y<<", z="<<z<<", f(x)="<<f(x)<<", g(x)="<<g(x))
        approximate_optimiser.minimisation_step(p, x,y,z, mu);
        if (i>=2 && i<17) { mu=mu/4; }
    }
    CONCLOG_PRINTLN("x="<<x<<", y="<<y<<", z="<<z<<", mu="<<mu<<", f(x)="<<f(x)<<", g(x)="<<g(x))

    FeasiblePrimalDualData<ValidatedNumber> r=this->check_minimality(p,g(x),x,y,z);
    if (definitely(r.is_feasible())) {
        CONCLOG_PRINTLN("r="<<r);
        FloatDPBounds v=f(cast_exact(x));
        return OptimalityCertificate({LocalOptimalityCertificate(v,cast_exact(x),cast_exact(y))});
    } else {
        ConstraintSolver cslvr;
        UpperBoxType oldB=B;
        cslvr.reduce(B,g,C);
        CONCLOG_PRINTLN("r="<<r<<", D="<<D<<", B="<<oldB<<"="<<B);
        if (definitely(B.is_empty())) {
            return InfeasibilityCertificate({LocalInfeasibilityCertificate(cast_exact_box(oldB),cast_exact(y))});
        } else {
            if (decide(measure(B)<measure(oldB)/2)) {
                return splitting_minimise_hotstarted(p,B,midpoint(B),y);
            } else {
                auto sB=split(B); auto& B1=sB.first; auto& B2=sB.second;
                auto r1=splitting_minimise_hotstarted(p,B1,midpoint(B),y);
                auto r2=splitting_minimise_hotstarted(p,B2,midpoint(B),y);

                if (std::holds_alternative<OptimalityCertificate>(r1) && std::holds_alternative<OptimalityCertificate>(r2)) {
                    // Both sub-boxes are feasible. Filter out clearly non-optimal points
                    auto o=catenate(std::get<OptimalityCertificate>(r1),std::get<OptimalityCertificate>(r2));
                    ARIADNE_DEBUG_ASSERT(not o.empty());
                    auto v=o[0].v;
                    for (auto lo : o) { v=min(v,lo.v); }
                    OptimalityCertificate ro;
                    for (auto lo : o) { if (possibly(lo.v<=v)) { ro.append(lo); } }
                    return ro;
                } else if (std::holds_alternative<OptimalityCertificate>(r1)) {
                    // B1 is feasible, but B2 is not
                    return r1;
                } else if (std::holds_alternative<OptimalityCertificate>(r2)) {
                    return r2;
                } else {
                    return InfeasibilityCertificate(catenate(std::get<InfeasibilityCertificate>(r1),std::get<InfeasibilityCertificate>(r2)));
                }
            }
        }
    }

}


auto
KarushKuhnTuckerOptimiser::splitting_feasible_hotstarted(
    const ValidatedFeasibilityProblem& p,
    UpperBoxType B, ApproximateVectorType x, ApproximateVectorType y) const
        -> Variant<FeasibilityCertificate,InfeasibilityCertificate>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    PrimalDualComplementaryInteriorPointOptimiser approximate_optimiser;

    // Extra variables for use of optimization solver
    ApproximateVectorType w=midpoint(C);
    Approximation<FLT> mu = one;
    ApproximateVectorType z=approximate_optimiser.compute_z(p,x,mu);

    // FIXME: Allow more steps
    Approximation<FLT> t=approximate_optimiser.compute_t(p,x);

    for(SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"t="<<t<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        approximate_optimiser.feasibility_step(p, x,y,z);
        t=approximate_optimiser.compute_t(p,x);
        if(probably(LogicalValue(t>0))) {
            break;
        }
    }
    CONCLOG_PRINTLN("t="<<t<<", x="<<x<<", y="<<y<<", g(x)="<<g(x))

    FeasiblePrimalDualData<ValidatedNumber> r=this->check_feasibility(p,x,y);
    if (definitely(r.is_feasible())) {
        CONCLOG_PRINTLN("r="<<r);
        return FeasibilityCertificate(Pair<ValidatedVector,ExactVector>(cast_exact(x),cast_exact(y)));
    } else {
        ConstraintSolver cslvr;
        UpperBoxType oldB=B;
        cslvr.reduce(B,g,C);
        CONCLOG_PRINTLN("r="<<r<<", D="<<D<<", B="<<oldB<<"="<<B);
        if (definitely(B.is_empty())) {
            return InfeasibilityCertificate({LocalInfeasibilityCertificate(cast_exact_box(oldB),cast_exact(y))});
        } else {
            if (decide(measure(B)<measure(oldB)/2)) {
                return splitting_feasible_hotstarted(p,B,midpoint(B),y);
            } else {
                auto sB=split(B);
                auto r1=splitting_feasible_hotstarted(p,sB.first,midpoint(B),y);
                if (std::holds_alternative<FeasibilityCertificate>(r1)) { return r1; }
                auto r2=splitting_feasible_hotstarted(p,sB.second,midpoint(B),y);
                if (std::holds_alternative<FeasibilityCertificate>(r2)) { return r2; }
                return InfeasibilityCertificate(catenate(std::get<InfeasibilityCertificate>(r1),std::get<InfeasibilityCertificate>(r2)));
            }
        }
    }

}


auto
KarushKuhnTuckerOptimiser::feasible_hotstarted(
    const ValidatedFeasibilityProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> FeasiblePrimalDualData<ValidatedNumber>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p<<", x0="<<x0<<", y0="<<y0)

    auto r=splitting_feasible_hotstarted(p,p.D,x0,y0);
    if (std::holds_alternative<FeasibilityCertificate>(r)) {
        FeasibilityCertificate fc=std::get<FeasibilityCertificate>(r);
        CONCLOG_PRINTLN("fc="<<fc)
        return FeasiblePrimalDualData<ValidatedNumber>(true,fc.x,fc.y);
    } else {
        InfeasibilityCertificate ifc=std::get<InfeasibilityCertificate>(r);
        CONCLOG_PRINTLN("ifc="<<ifc)
        return FeasiblePrimalDualData<ValidatedNumber>(false,cast_singleton(ifc[0].B),ifc[0].y);
    }
}

auto
KarushKuhnTuckerOptimiser::nonsplitting_feasible_hotstarted(
    const ValidatedFeasibilityProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> FeasiblePrimalDualData<ValidatedNumber>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("p="<<p);

    InteriorPointOptimiser approximate_optimiser;

    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    ARIADNE_ASSERT(x0.size()==g.argument_size());
    ARIADNE_ASSERT(y0.size()==g.result_size());

    Vector<Approximation<FLT>> x=x0;
    Vector<Approximation<FLT>> y=y0;

    // Extra variables for use of optimization solver
    Vector<Approximation<FLT>> w=midpoint(C);
    Approximation<FLT> mu = one;
    Vector<Approximation<FLT>> z=approximate_optimiser.compute_z(p,x0,mu);

    // FIXME: Allow more steps
    Approximation<FLT> t=approximate_optimiser.compute_t(p,x);

    CONCLOG_PRINTLN_AT(1,"D="<<D<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
    for(SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"t="<<t<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        CONCLOG_PRINTLN_AT(2,"w="<<w<<", z="<<z);
        approximate_optimiser.feasibility_step(p, x,y,z);
        t=approximate_optimiser.compute_t(p,x);
        if(probably(LogicalValue(t>0))) {
            return this->check_feasibility(p,x,y);
        }
    }
    CONCLOG_PRINTLN("t="<<t<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
    return this->check_feasibility(p,x,y);
}


auto
KarushKuhnTuckerOptimiser::check_feasibility(
    const ValidatedFeasibilityProblem& p,
    const ApproximateVectorType& xa, const ApproximateVectorType& ya) const
        -> FeasiblePrimalDualData<ValidatedNumber>
{
    CONCLOG_SCOPE_CREATE;

    auto& g=p.g;
    const Vector<Exact<FLT>> x=cast_exact(xa);
    const Vector<Exact<FLT>> y=cast_exact(ya);

    CONCLOG_PRINTLN("p="<<p)
    CONCLOG_PRINTLN("x="<<x<<", y="<<y<<", g(x)="<<g(x))

    FeasibilityChecker fc;
    ValidatedKleenean k=fc.check_feasibility(p,x,y);
    return FeasiblePrimalDualData<ValidatedNumber>(k,to_generic(x),to_generic(y));
}


auto
KarushKuhnTuckerOptimiser::check_minimality(
    const ValidatedOptimisationProblem& p,
    const ApproximateVectorType& wa, const ApproximateVectorType& xa, const ApproximateVectorType& ya, const ApproximateVectorType& za) const
        -> FeasiblePrimalDualData<ValidatedNumber>
{
    //CONCLOG_SCOPE_CREATE;
    auto& f=p.f; auto& g=p.g;

    CONCLOG_PRINTLN("KarushKuhnTuckerOptimiser::check_minimality(ValidatedOptimisationProblem p, ApproximateVectorType wa, xa, ya, za)");
    CONCLOG_PRINTLN("wa="<<wa<<", xa="<<xa<<", ya="<<ya<<", za="<<za<<", f(xa)="<<f(xa)<<", g(xa)="<<g(xa));

    ExactDouble se(std::numeric_limits<float>::epsilon());
    se=ExactDouble(1.0/1024);
    se=ExactDouble(1.0/1048576);
    FloatDPBounds e=se*FloatDPBounds(-1,+1,dp);
    auto cast_exact_widen = [&e](auto va){return eadd(cast_exact(va)*(1+e),e);};

    FloatDPBoundsVector w=cast_exact_widen(wa);
    FloatDPBoundsVector x=cast_exact_widen(xa);
    FloatDPBoundsVector y=cast_exact_widen(ya);
    FloatDPBoundsVector z=cast_exact_widen(za);

    // auto restrict_number = [](FloatDPBounds xb, UpperIntervalType ivl) {
    //     return FloatDPBounds(max(xb.lower(),ivl.lower_bound()),min(xb.upper(),ivl.upper_bound())); };
    // auto restrict_vector = [&](Vector<FloatDPBounds> v, ExactBoxType bx) {
    //     return Vector<FloatDPBounds>( v.size(), [&](SizeType i){return restrict_number(v[i],bx[i]);} ); };
    // w=restrict_vector(w,p.C);
    // x=restrict_vector(x,p.D);

    CONCLOG_PRINTLN_AT(0, "w="<<w<<", x="<<x<<", y="<<y<<", z="<<z<<", f(x)="<<f(x)<<", g(x)="<<g(x));
    return this->check_minimality(p,w,x,y,z);
}

auto
KarushKuhnTuckerOptimiser::check_minimality(
    const ValidatedOptimisationProblem& p,
    const FloatDPBoundsVector& w_, const FloatDPBoundsVector& x_, const FloatDPBoundsVector& y_, const FloatDPBoundsVector& z_) const
        -> FeasiblePrimalDualData<ValidatedNumber>
{
    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;
    Vector<Bounds<FLT>> w(w_), x(x_), y(y_), z(z_);

    //CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("PrimalDualComplementaryInteriorPointOptimiser::check_minimality(ValidatedOptimisationProblem p, Vector<Bounds<FLT>> w, x, y, z)");
    CONCLOG_PRINTLN_AT(1, "w="<<w<<", x="<<x<<", y="<<y<<", z="<<z<<", f(x)="<<f(x)<<", g(x)="<<g(x));

    // Set-up interval Newton contractor
    // Want to solve KKT conditions:
    // g(x)-w=0
    // df(x)+dg(x)'*y+z=0
    // (w-cl)(cu-w)y=0
    // (x-dl)(du-x)z=0

    // [-I  A  0  0]
    // [ 0  Q AT  I]
    // [ Y  0  W  0]
    // [ 0  Z  0  X]

    auto Cl=lower_bounds(C);
    auto Cu=upper_bounds(C);
    auto Dl=lower_bounds(D);
    auto Du=upper_bounds(D);

    auto midpoint = [](ExactIntervalType const& ivl) {
        return FloatDPBounds(hlf(ivl.lower_bound()+ivl.upper_bound())); };
    auto midpoints = [&midpoint](ExactBoxType const& bx) {
        return Vector<FloatDPBounds>( bx.dimension(), [&midpoint,&bx](SizeType i){return midpoint(bx[i]);} ); };

    bool is_refinement = false;
    int step=0;

    while (not is_refinement && step<8) {
        auto ddfx=differential(f,x,2);
        auto ddgx=differential(g,x,2);

        auto A=ddgx.jacobian();
        auto Q=ddfx.hessian();
        for (SizeType i=0; i!=y.size(); ++i) { Q+=y[i]*ddgx[i].hessian(); }

        auto W=DiagonalMatrix<Bounds<FLT>>(emul(w-Cl,Cu-w));
        auto X=DiagonalMatrix<Bounds<FLT>>(emul(x-Dl,Du-x));
        auto Y=DiagonalMatrix<Bounds<FLT>>(emul(2,midpoints(C)-w,y));
        auto Z=DiagonalMatrix<Bounds<FLT>>(emul(2,midpoints(D)-x,z));

        Vector<Exact<FLT>> we=cast_exact(w);
        Vector<Exact<FLT>> xe=cast_exact(x);
        Vector<Exact<FLT>> ye=cast_exact(y);
        Vector<Exact<FLT>> ze=cast_exact(z);
        Covector<Bounds<FLT>> ce=f.gradient(xe);
        Matrix<Bounds<FLT>> Ae=g.jacobian(xe);

        Vector<Bounds<FLT>> rw=g(xe)-we;
        Vector<Bounds<FLT>> rx=transpose(ce)+transpose(Ae)*ye+ze;
        Vector<Bounds<FLT>> ry=emul(we-Cl,Cu-we,ye);
        Vector<Bounds<FLT>> rz=emul(xe-Dl,Du-xe,ze);

        SlackPrimalDualComplementaryMatrix<Bounds<FLT>> S(Q,A,W,X,Y,Z);
        SlackPrimalDualComplementaryData<Bounds<FLT>> r(rw,rx,ry,rz);
        SlackPrimalDualComplementaryData<Exact<FLT>> q(we,xe,ye,ze);

        auto Sa=Matrix<Approximation<FLT>>(S.assemble());
        CONCLOG_PRINTLN("step="<<step);
        CONCLOG_PRINTLN("C="<<C<<", D="<<D);
        CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y<<", z="<<z);
        CONCLOG_PRINTLN("Sa="<<pretty(Sa))

        try {
            SlackPrimalDualComplementaryData<Bounds<FLT>> nd=q-S.solve(r);
            auto& nw=nd.w; auto& nx=nd.x; auto& ny=nd.y; auto& nz=nd.z;
            CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y<<", z="<<z)
            CONCLOG_PRINTLN("nw="<<nw<<", nx="<<nx<<", ny="<<ny<<", nz="<<nz)
            CONCLOG_PRINTLN("refines(nw,w)="<<refines(nw,w)<<", refines(nx,x)="<<refines(nx,x)<<", refines(ny,y)="<<refines(ny,y)<<", refines(nz,z)="<<refines(nz,z));
            CONCLOG_PRINTLN("inconsistent(nw,w)="<<inconsistent(nw,w)<<", inconsistent(nx,x)="<<inconsistent(nx,x)<<", inconsistent(ny,y)="<<inconsistent(ny,y)<<", inconsistent(nz,z)="<<inconsistent(nz,z));
            if (refines(nw,w) && refines(nx,x) && refines(ny,y) && refines(nz,z)) {
                return FeasiblePrimalDualData<ValidatedNumber>(true,nx,ny);
            } else if (inconsistent(nw,w) || inconsistent(nx,x) || inconsistent(ny,y) || inconsistent(nz,z)) {
                return FeasiblePrimalDualData<ValidatedNumber>(false,x,ny);
            } else {
                w=refinement(w,nw); x=refinement(x,nx); y=refinement(y,ny); z=refinement(z,nz);
                ++step;
            }
        } catch (const SingularMatrixException& e) {
            std::cerr<<"S="<<S.assemble()<<"\nSa="<<pretty(Sa)<<"\n" << std::flush;
            //std::cerr<<"inv(Sa)="<<pretty(inverse(Sa))<<"\n";
            throw(e);
        }
    }
    return FeasiblePrimalDualData<ValidatedNumber>(indeterminate,x,y);
}

//------- InfeasibleKarushKuhnTuckerOptimiser -----------------------------------//

auto
InfeasibleKarushKuhnTuckerOptimiser::clone() const -> InfeasibleKarushKuhnTuckerOptimiser*
{
    return new InfeasibleKarushKuhnTuckerOptimiser(*this);
}


auto
InfeasibleKarushKuhnTuckerOptimiser::minimise(ValidatedOptimisationProblem p) const -> ValidatedVector
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("InfeasibleKarushKuhnTuckerOptimiser::minimise(ValidatedOptimisationProblem p)");
    CONCLOG_PRINTLN("p="<<p);

    ARIADNE_PRECONDITION(p.D.has_nonempty_interior());
    ARIADNE_PRECONDITION(p.C.has_nonempty_interior());
    ApproximateVectorType x0=p.D.midpoint();
    ApproximateNumericType mu0(1,x0.zero_element().precision());
    InfeasibleInteriorPointOptimiser approximate_optimiser;
    ApproximateVectorType y0=approximate_optimiser.compute_y(p,x0,mu0);
    return this->minimise_hotstarted(p,x0,y0).primal();
}


auto
InfeasibleKarushKuhnTuckerOptimiser::minimise_hotstarted(
    const ValidatedOptimisationProblem& p,
    const ApproximateVectorType& x0, const ApproximateVectorType& y0) const
        -> ValuePrimalDualData<ValidatedNumber>
{
    CONCLOG_PRINTLN("InfeasibleKarushKuhnTuckerOptimiser::minimise_hotstarted(ValidatedOptimisationProblem p, ApproximateVectorType xa0, ya0)");

    UpperBoxType B=p.D;

    auto r = this->splitting_minimise_hotstarted(p,B,x0,y0);

    auto ro0 = std::get<OptimalityCertificate>(r)[0];
    return ValuePrimalDualData<ValidatedNumber>(ro0.value(),ro0.primal(),ro0.dual());
}



auto
InfeasibleKarushKuhnTuckerOptimiser::_splitting_minimise_subdivide(
    const ValidatedOptimisationProblem& p,
    UpperBoxType B, ApproximateVectorType x, ApproximateVectorType y) const
        -> Variant<OptimalityCertificate,InfeasibilityCertificate>
{
    auto sB=split(B); auto& B1=sB.first; auto& B2=sB.second;
    auto r1=splitting_minimise_hotstarted(p,B1,midpoint(B),y);
    auto r2=splitting_minimise_hotstarted(p,B2,midpoint(B),y);

    if (std::holds_alternative<OptimalityCertificate>(r1) && std::holds_alternative<OptimalityCertificate>(r2)) {
        // Both sub-boxes are feasible. Filter out clearly non-optimal points
        auto o=catenate(std::get<OptimalityCertificate>(r1),std::get<OptimalityCertificate>(r2));
        ARIADNE_DEBUG_ASSERT(not o.empty());
        auto v=o[0].v;
        for (auto lo : o) { v=min(v,lo.v); }
        OptimalityCertificate ro;
        for (auto lo : o) { if (possibly(lo.v<=v)) { ro.append(lo); } }
        return ro;
    } else if (std::holds_alternative<OptimalityCertificate>(r1)) {
        // B1 is feasible, but B2 is not
        return r1;
    } else if (std::holds_alternative<OptimalityCertificate>(r2)) {
        return r2;
    } else {
        return InfeasibilityCertificate(catenate(std::get<InfeasibilityCertificate>(r1),std::get<InfeasibilityCertificate>(r2)));
    }
}

auto
InfeasibleKarushKuhnTuckerOptimiser::splitting_minimise_hotstarted(
    const ValidatedOptimisationProblem& p,
    UpperBoxType B, ApproximateVectorType x, ApproximateVectorType y) const
        -> Variant<OptimalityCertificate,InfeasibilityCertificate>
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("")
    CONCLOG_PRINTLN("InfeasibleKarushKuhnTuckerOptimiser::splitting_minimise_hotstarted(p, B0,xa0,ya0)");
    CONCLOG_PRINTLN("D="<<p.D<<", C="<<p.C);
    CONCLOG_PRINTLN("B0="<<B<<", xa0="<<x<<", ya0="<<y);

    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;

    SlackPrimalDualComplementaryInteriorPointOptimiser approximate_optimiser;

    // Extra variables for use of optimization solver
    Approximation<FLT> mu = one;
    ApproximateVectorType z=approximate_optimiser.compute_z(p,x,mu);
    ApproximateVectorType w=approximate_optimiser.compute_w(p,x,y,mu);

    CONCLOG_PRINTLN("w0="<<w<<", x0="<<x<<", y0="<<y<<", z0="<<z<<", mu0="<<mu<<", f(x0)="<<f(x)<<", g(x0)="<<g(x));

#warning FIXME: Change to pure mimimisation. Remove verbosity switch
    for(SizeType i=0; i!=18; ++i) {
        try {
            auto v=ConcLog::Logger::instance().configuration().verbosity();
            ConcLog::Logger::instance().configuration().set_verbosity(0);
            approximate_optimiser.minimisation_step(p, w,x,y,z, mu);
            ConcLog::Logger::instance().configuration().set_verbosity(v);
        } catch (SingularMatrixException e) {
            return this->_splitting_minimise_subdivide(p,B,x,y);
        }
        if (i>=2 && i<17) { mu=mu/4; }
        CONCLOG_PRINTLN_AT(1,"w"<<i+1<<"="<<w<<", x="<<x<<", y="<<y<<", z="<<z<<", mu="<<mu<<", f(x)="<<f(x)<<", g(x)="<<g(x));
    }
    CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y<<", z="<<z<<", mu="<<mu<<", f(x)="<<f(x)<<", g(x)="<<g(x));

    Tuple<ValidatedKleenean,ValidatedVector,ValidatedVector> r=this->check_minimality(p,w,x,y,z);

    if (definitely(std::get<0>(r))) {
        CONCLOG_PRINTLN("r="<<r);
        FloatDPBounds v=f(cast_exact(x));
        return OptimalityCertificate({LocalOptimalityCertificate(v,cast_exact(x),cast_exact(y))});
    } else {
        ConstraintSolver cslvr;
        UpperBoxType oldB=B;
        cslvr.reduce(B,g,C);
        CONCLOG_PRINTLN("r="<<r<<", D="<<D<<", B="<<oldB<<"="<<B);
        if (definitely(B.is_empty())) {
            return InfeasibilityCertificate({LocalInfeasibilityCertificate(cast_exact_box(oldB),cast_exact(y))});
        } else {
            if (decide(measure(B)<measure(oldB)/2)) {
                return splitting_minimise_hotstarted(p,B,midpoint(B),y);
            } else {
                return this->_splitting_minimise_subdivide(p,B,x,y);
            }
        }
    }

}

auto InfeasibleKarushKuhnTuckerOptimiser::
feasible(ValidatedFeasibilityProblem p) const -> ValidatedKleenean
{
    CONCLOG_SCOPE_CREATE
    CONCLOG_PRINTLN("p="<<p);

    auto& D=p.D; auto& g=p.g; auto& C=p.C;

    ApproximateScalarMultivariateFunction f(EuclideanDomain(D.dimension()));
    auto R=intersection(cast_exact_box(widen(apply(g,D),1)),C);

    ApproximateOptimisationProblem optp(f,D,g,R);

    InfeasibleInteriorPointOptimiser approximate_optimiser;

    UniquePointer<InfeasibleInteriorPointOptimiser::StepData> d_ptr(approximate_optimiser.initial_step_data(optp));
    InfeasibleInteriorPointOptimiser::StepData& d=*d_ptr;

    ApproximateVectorType& x=d.x;
    ApproximateVectorType& y=d.y;

    static const ExactDouble MU_MIN = 1e-12_pr;

    // FIXME: Allow more steps
    for(SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        approximate_optimiser.minimisation_step(optp,d);
        if(FeasibilityChecker().validate_feasibility(p,cast_exact(x))) {
            CONCLOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
            CONCLOG_PRINTLN("Feasible");
            return true;
        }
        if(FeasibilityChecker().validate_infeasibility(p,cast_exact(y))) {
            CONCLOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
            CONCLOG_PRINTLN("Infeasible");
            return false;
        }
        if(d.mu.raw()<MU_MIN) {
            break;
        }
    }
    CONCLOG_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
    CONCLOG_PRINTLN("Indeterminate");
    return indeterminate;
}



auto InfeasibleKarushKuhnTuckerOptimiser::
feasible_hotstarted(ValidatedFeasibilityProblem p,
                    const SlackPrimalDualData<ApproximateNumericType>& wxy0) const
                        -> Tuple<ValidatedKleenean,ValidatedVector,ValidatedVector>
{
    ARIADNE_NOT_IMPLEMENTED;
}


auto
InfeasibleKarushKuhnTuckerOptimiser::check_minimality(
    const ValidatedOptimisationProblem& p,
    const ApproximateVectorType& wa, const ApproximateVectorType& xa, const ApproximateVectorType& ya, const ApproximateVectorType& za) const
        -> Tuple<ValidatedKleenean,ValidatedVector,ValidatedVector>
{
    auto& f=p.f; auto& g=p.g;

    CONCLOG_PRINTLN("")
    CONCLOG_PRINTLN("InfeasibleKarushKuhnTuckerOptimiser::check_minimality(ValidatedOptimisationProblem p, ApproximateVectorType wa, xa, ya, za)");
    CONCLOG_PRINTLN("wa="<<wa<<", xa="<<xa<<", ya="<<ya<<", za="<<za<<", f(xa)="<<f(xa)<<", g(xa)="<<g(xa));

    Approximation<FLT> sa=sqrt(FLT::eps(xa.zero_element().precision()));
    FLT se=cast_exact(sa);
    CONCLOG_PRINTLN("se="<<se);
#warning
    //se=ExactDouble(1.0/1024);
    //se=ExactDouble(1.0/1048576);
    FloatDPBounds e=se*FloatDPBounds(-1,+1,dp);
    auto cast_exact_widen = [&e](auto va){return eadd(cast_exact(va)*(1+e),e);};

    FloatDPBoundsVector w=cast_exact_widen(wa);
    FloatDPBoundsVector x=cast_exact_widen(xa);
    FloatDPBoundsVector y=cast_exact_widen(ya);
    FloatDPBoundsVector z=cast_exact_widen(za);

    // auto restrict_number = [](FloatDPBounds xb, UpperIntervalType ivl) {
    //     return FloatDPBounds(max(xb.lower(),ivl.lower_bound()),min(xb.upper(),ivl.upper_bound())); };
    // auto restrict_vector = [&](Vector<FloatDPBounds> v, ExactBoxType bx) {
    //     return Vector<FloatDPBounds>( v.size(), [&](SizeType i){return restrict_number(v[i],bx[i]);} ); };
    // w=restrict_vector(w,p.C);
    // x=restrict_vector(x,p.D);

    CONCLOG_PRINTLN_AT(0, "w="<<w<<", x="<<x<<", y="<<y<<", z="<<z<<", f(x)="<<f(x)<<", g(x)="<<g(x));
    return this->check_minimality(p,w,x,y,z);
}

auto
InfeasibleKarushKuhnTuckerOptimiser::check_minimality(
    const ValidatedOptimisationProblem& p,
    const ValidatedVectorType& w_, const ValidatedVectorType& x_, const ValidatedVectorType& y_, const ValidatedVectorType& z_) const
        -> Tuple<ValidatedKleenean,ValidatedVector,ValidatedVector>
{
    auto& f=p.f; auto& D=p.D; auto& g=p.g; auto& C=p.C;
    Vector<Bounds<FLT>> w(w_), x(x_), y(y_), z(z_);

    auto restrict_number = [](FloatDPBounds xb, UpperIntervalType ivl) {
        return FloatDPBounds(max(xb.lower(),ivl.lower_bound()),min(xb.upper(),ivl.upper_bound())); };
    auto restrict_vector = [&](Vector<FloatDPBounds> v, ExactBoxType bx) {
        return Vector<FloatDPBounds>( v.size(), [&](SizeType i){return restrict_number(v[i],bx[i]);} ); };

    CONCLOG_PRINTLN("\n")
    CONCLOG_PRINTLN("InfeasibleKarushKuhnTuckerOptimiser::check_minimality(ValidatedOptimisationProblem p, Vector<Bounds<FLT>> w, x, y, z)");
    CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y<<", z="<<z);
    CONCLOG_PRINTLN("f(x)="<<f(x)<<", g(x)="<<g(x));

#warning Do we want to restrict?
//    w=restrict_vector(w,p.C);
//    x=restrict_vector(x,p.D);

    CONCLOG_PRINTLN("w="<<w<<", x="<<x);
    CONCLOG_PRINTLN("f(x)="<<f(x)<<", g(x)="<<g(x));


    // Set-up interval Newton contractor
    // Want to solve KKT conditions:
    // g(x)-w=0
    // df(x)+dg(x)'*y+z=0
    // (w-cl)(cu-w)y=0
    // (x-dl)(du-x)z=0

    // [-I  A  0  0]
    // [ 0  Q AT  I]
    // [ Y  0  W  0]
    // [ 0  Z  0  X]

    auto Cl=lower_bounds(C);
    auto Cu=upper_bounds(C);
    auto Dl=lower_bounds(D);
    auto Du=upper_bounds(D);

    auto midpoint = [](ExactIntervalType const& ivl) {
        return FloatDPBounds(hlf(ivl.lower_bound()+ivl.upper_bound())); };
    auto midpoints = [&midpoint](ExactBoxType const& bx) {
        return Vector<FloatDPBounds>( bx.dimension(), [&midpoint,&bx](SizeType i){return midpoint(bx[i]);} ); };

    bool is_refinement = false;
    int step=0;

    while (not is_refinement && step<1) {
        auto ddfx=differential(f,x,2);
        auto ddgx=differential(g,x,2);

        auto A=ddgx.jacobian();
        auto Q=ddfx.hessian();
        for (SizeType i=0; i!=y.size(); ++i) { Q+=y[i]*ddgx[i].hessian(); }


#warning
//        Tuple<DiagonalMatrix<X>,DiagonalMatrix<X>,Vector<X>> compute_linear_equations(const BX& C, const Vector<X>& w, const Vector<X>& y, const X& mu)

        auto YWrw = compute_linear_equations(C,w,y);
        DiagonalMatrix<Bounds<FLT>>& Y=std::get<0>(YWrw);
        DiagonalMatrix<Bounds<FLT>>& W=std::get<1>(YWrw);
        auto ZXrx = compute_linear_equations(D,x,z);
        DiagonalMatrix<Bounds<FLT>>& Z=std::get<0>(ZXrx);
        DiagonalMatrix<Bounds<FLT>>& X=std::get<1>(ZXrx);

        Vector<Exact<FLT>> we=cast_exact(w);
        Vector<Exact<FLT>> xe=cast_exact(x);
        Vector<Exact<FLT>> ye=cast_exact(y);
        Vector<Exact<FLT>> ze=cast_exact(z);
        Covector<Bounds<FLT>> ce=f.gradient(xe);
        Matrix<Bounds<FLT>> Ae=g.jacobian(xe);

        Vector<Bounds<FLT>> rwe=g(xe)-we;
        Vector<Bounds<FLT>> rxe=transpose(ce)+transpose(Ae)*ye+ze;
        Vector<Bounds<FLT>> rye = std::get<2>(compute_linear_equations(C,we,ye));
        Vector<Bounds<FLT>> rze = std::get<2>(compute_linear_equations(D,xe,ze));

        SlackPrimalDualComplementaryMatrix<Bounds<FLT>> S(Q,A,W,X,Y,Z);
#warning Check order of variables in matrix
        SlackPrimalDualComplementaryData<Bounds<FLT>> re(rwe,rxe,rye,rze);
        SlackPrimalDualComplementaryData<Exact<FLT>> qe(we,xe,ye,ze);

        auto Sa=Matrix<Approximation<FLT>>(S.assemble());
        CONCLOG_PRINTLN("step="<<step);
        CONCLOG_PRINTLN("C="<<C<<", D="<<D);
        CONCLOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y<<", z="<<z);
        CONCLOG_PRINTLN_AT(1,"Q="<<Q<<", A="<<A);
        CONCLOG_PRINTLN_AT(1,"W="<<W<<", X="<<X<<", Y="<<Y<<", Z="<<Z);
        CONCLOG_PRINTLN("we="<<we<<", xe="<<xe<<", ye="<<ye<<", ze="<<ze);
        CONCLOG_PRINTLN("rwe="<<rwe<<", rxe="<<rxe<<", rye="<<rye<<", rze="<<rze);
        CONCLOG_PRINTLN("Sa="<<pretty(Sa))
        CONCLOG_PRINTLN("re="<<re)

        try {
            SlackPrimalDualComplementaryData<Bounds<FLT>> nd=qe-S.solve(re);
            auto& nw=nd.w; auto& nx=nd.x; auto& ny=nd.y; auto& nz=nd.z;
            CONCLOG_PRINTLN("w="<<w)
            CONCLOG_PRINTLN("nw="<<nw)
            CONCLOG_PRINTLN("x="<<x)
            CONCLOG_PRINTLN("nx="<<nx)
            CONCLOG_PRINTLN("y="<<y)
            CONCLOG_PRINTLN("ny="<<ny)
            CONCLOG_PRINTLN("z="<<z)
            CONCLOG_PRINTLN("nz="<<nz)
            CONCLOG_PRINTLN("refines(nw,w)="<<refines(nw,w)<<", refines(nx,x)="<<refines(nx,x)<<", refines(ny,y)="<<refines(ny,y)<<", refines(nz,z)="<<refines(nz,z));
            CONCLOG_PRINTLN("inconsistent(nw,w)="<<inconsistent(nw,w)<<", inconsistent(nx,x)="<<inconsistent(nx,x)<<", inconsistent(ny,y)="<<inconsistent(ny,y)<<", inconsistent(nz,z)="<<inconsistent(nz,z));
            if (refines(nw,w) && refines(nx,x) && refines(ny,y) && refines(nz,z)) {
                return make_tuple(true,nx,ny);
            } else if (inconsistent(nw,w) || inconsistent(nx,x) || inconsistent(ny,y) || inconsistent(nz,z)) {
                return make_tuple(false,x,ny);
            } else {
                w=refinement(w,nw); x=refinement(x,nx); y=refinement(y,ny); z=refinement(z,nz);
                ++step;
            }
        } catch (const SingularMatrixException& e) {
            std::cerr<<"S="<<(S.assemble())<<"\nSa="<<pretty(Sa)<<"\n"<<std::flush;
            //std::cerr<<"inv(Sa)="<<pretty(lu_inverse(Sa))<<"\n"<<std::flush;
            throw(e);
        }
    }

    return make_tuple(indeterminate,x,y);
}





/*

// Solve max log(x-xl) + log(xu-x) + log(zu-z) + log(z-zl) such that g(x)=z
//   if zl[i]=zu[i] then z=zc is hard constraint
//   alternatively, use the penalty (z[j]-zc[j])^2/2 instead

// KKT conditions
//     1/(x-xl) - 1/(xu-x) + y Dg(x) = 0
//     1/(z-zl) - 1/(zu-z) - y = 0
//     g(x) - z = 0
//   If zu[j]=inf, then 1/(z[j-zl[j]) - y[j] = 0, so y[j]>=0
//   If zl[j]=zu[j], then use the equation z[j]=zc[j] instead
//
// Re-write KKT conditions for x,z as
//     (xu-xl) + y g(x) (x-xl)(xu-x) = 0
//     (zu-zl) - y(z-zl)(zu-z) = 0
//  Or 1 - y(z-zl)(zu-z)/(zu-zl) = 0
//   If zu[j] = inf, then 1 - y(z-zl) = 0
//      zl[j]=zu[j]=zc[j], then y(z-zc)^2 = 0
//
// Derivative matrix
//     - 1/(x-xl)^2 - 1/(xu-x)^2 + y D^2g = 0
//
// PROBLEM:
//   Dg can be singular, even at intermediate points
//
// FJ conditions
//     mu/(x-xl) - mu/(xu-x) + y D g(x) = 0
//     mu/(z-zl) - mu/(zu-z) - y = 0
//     g(x) - z = 0
//     sum y^2 - mu = 0
//   If zu[j]=inf, then mu/(z[j-zl[j]) - y[j] = 0, so y[j]>=0
//   If zl[j]=zu[j], then -(z-zc)/mu - y = 0 instead
//
// Re-write FJ conditions for x,z as
// Derivative matrix
//     - mu/(x-xl)^2 dx - mu/(xu-x)^2 dx + y D^2g dx + Dg^T dy + (1/(x-xl) - (1/xu-x)) dmu
//     - mu/(z-zl)^2 dz - mu/(zu-z)^2 dz + I dy + (1/(z-zl) - (1/zu-z)) dmu
//     Dg dx - I dz
//    2y . dy - dmu
Void PenaltyFunctionOptimiser::
feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                 RawFloatDPVector& x, RawFloatDPVector& y, RawFloatDPVector& z) const
{
    CONCLOG_PRINTLN_AT(2,"feasibility_step");
    RawFloatDPVector xl=lower_bounds(D); RawFloatDPVector xu=upper_bounds(D);
    RawFloatDPVector zl=lower_bounds(C); RawFloatDPVector zu=upper_bounds(C);

    const SizeType n=x.size();
    const SizeType m=y.size();

    CONCLOG_PRINTLN_AT(4,"x="<<x<<" y="<<y<<" z="<<z);
    Vector<FloatDPDifferential> ddx = FloatDPDifferential::variables(2,x);
    Vector<FloatDPDifferential> ddgx = g.evaluate(ddx);

    FloatDPMatrix A=ddgx.jacobian();
    CONCLOG_PRINTLN_AT(6,"A="<<A);
    RawFloatDPVector v = join(join(x,z),y);

    RawFloatDPVector r(n+2*m,n+2*m);
    project(r,range(0,n)) = y * A;
    for(SizeType i=0; i!=n; ++i) {
        r[i] += ( rec(x[i]-xl[i]) - rec(xu[i]-x[i]) );
    }
    for(SizeType j=0; j!=m; ++j) {
        if(zl[j]==zu[j]) { assert(zu[j]==zl[j]); r[n+j] = z[j]-zl[j]; }
        else { r[n+j] = ( rec(z[j]-zl[j]) - rec(zu[j]-z[j]) - y[j] ); }
    }
    project(r,range(n+m,n+2*m)) = ddgx.value() - z;
    r[n+2*m]=0.0;
    CONCLOG_PRINTLN_AT(5,"r="<<r);

    FloatDPMatrix S(n+2*m+1,n+2*m+1);
    for(SizeType j=0; j!=m; ++j) {
        FloatDPMatrix H=ddgx[j].hessian();
        for(SizeType i1=0; i1!=n; ++i1) {
            for(SizeType i2=0; i2!=n; ++i2) {
                S[i1][i2]+=y[j]*H[i1][i2];
            }
        }
    }
    for(SizeType j=0; j!=m; ++j) {
        for(SizeType i=0; i!=n; ++i) {
            S[i][j+m+n]=A[j][i];
            S[j+m+n][i]=A[j][i];
        }
    }
    for(SizeType j=0; j!=m; ++j) {
        S[n+j][n+m+j] = -1.0;
        S[n+m+j][n+j] = -1.0;
        //if(zl[j]==zu[j]) { S[n+j][n+j] = -1.0; S[n+j][n+m+j] = 0.0; }
        if(zl[j]==zu[j]) { S[n+j][n+j] = +inf; }
        else { S[n+j][n+j] = - rec(sqr(z[j]-zl[j])) - rec(sqr(zu[j]-z[j])); }
    }
    for(SizeType i=0; i!=n; ++i) {
        S[i][i]-= rec(sqr(xu[i]-x[i]));
        S[i][i]-= rec(sqr(x[i]-xl[i]));
    }

    for(SizeType i=0; i!=n; ++i) {
        S[i][n+2*m] -= rec(xu[i]-x[i]);
        S[i][n+2*m] += rec(x[i]-xl[i]);
    }

    for(SizeType j=0; j!=n; ++j) {
        //S[n+m+j][n+m+j] = -1.0/1024;
    }

    CONCLOG_PRINTLN_AT(5,"S="<<S);
    //CONCLOG_PRINTLN_AT(5,"S="<<std::fixed<<pretty(S));

    FloatDPMatrix Sinv = inverse(S);
    //CONCLOG_PRINTLN_AT(9,"Sinv="<<Sinv);
    //CONCLOG_PRINTLN_AT(5,"Sinv="<<std::fixed<<pretty(Sinv));

    RawFloatDPVector dv = Sinv * r;
    CONCLOG_PRINTLN_AT(5,"dv="<<dv);

    FloatDP alpha = 1.0;
    RawFloatDPVector nv = v-dv;
    while(!contains(D,RawFloatDPVector(project(nv,range(0,n)))) || !contains(C,RawFloatDPVector(project(nv,range(n,n+m)))) ) {
        alpha *= 0.75;
        nv = v-alpha*dv;
    }

    CONCLOG_PRINTLN_AT(4,"nv="<<nv<<" a="<<alpha);

    x=project(nv,range(0,n));
    z=project(nv,range(n,n+m));
    y=project(nv,range(n+m,n+2*m));
    CONCLOG_PRINTLN_AT(4,"g(x)-z="<<g(x)-z);

}
*/

//------- IntervalOptimiser -----------------------------------//

IntervalOptimiser* IntervalOptimiser::clone() const {
    return new IntervalOptimiser(*this);
}


// Solve equations y Dh(x) - 1/(x-xl) + 1/(xu-x) = 0; h(x) = 0
ValidatedKleenean IntervalOptimiser::
feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("D="<<D<<", h="<<h);

    const SizeType n=D.size();

    FloatDPBoundsVector zl(n,dp), zu(n,dp);
    ExactFloatDPVectorType xl = Ariadne::lower_bounds(D);
    ExactFloatDPVectorType xu = Ariadne::upper_bounds(D);

    FloatDPBoundsVector x=cast_singleton(D);
    FloatDPBoundsVector y(h.result_size(),FloatDPBounds(-1,+1,dp));
    FloatDPBounds mu(0,1,dp);

    for(SizeType i=0; i!=8; ++i) {
        this->feasibility_step(xl,xu,h,x,y,zl,zu,mu);
    }

    return indeterminate;
}

Void IntervalOptimiser::
feasibility_step(const ExactFloatDPVectorType& xl, const ExactFloatDPVectorType& xu, const ValidatedVectorMultivariateFunction& h,
                 FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& zl, FloatDPBoundsVector zu, FloatDPBounds& mu) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("[x]="<<x<<" [lambda]="<<y<<", [zl]="<<zl<<", [zu]="<<zu<<" [mu]="<<mu);

    const SizeType n=x.size();
    const SizeType m=y.size();

    FloatDPBoundsVector mx=midpoint(x);
    FloatDPBoundsVector my=midpoint(y);
    FloatDPBoundsVector mzl=midpoint(zl);
    FloatDPBoundsVector mzu=midpoint(zu);
    FloatDPBounds mmu(midpoint(mu));
    CONCLOG_PRINTLN_AT(1,"x~"<<x<<" lambda~="<<y<<", mu~"<<mu);

    // Solve equations y Dh(x) - zl + zu = 0; h(x) = 0; (x-xl).zl - mu = 0;  (xu-x).zu - mu = 0; Sum_j y_j^2 - mu = 0
    Vector<FloatDPBoundsDifferential> ddhx=h.evaluate(FloatDPBoundsDifferential::variables(2,x));
    Vector<FloatDPBoundsDifferential> dhmx=h.evaluate(FloatDPBoundsDifferential::variables(1,mx));
    FloatDPBoundsMatrix A = ddhx.jacobian();
    FloatDPBoundsMatrix mA = dhmx.jacobian();
    CONCLOG_PRINTLN_AT(1,"A="<<A<<" b="<<ddhx.value());

    FloatDPBoundsVector rx = transpose(mA) * my;
    for(SizeType j=0; j!=n; ++j) {
        rx[j] -= mmu*rec(mx[j]-xl[j]);
        rx[j] += mmu*rec(xu[j]-mx[j]);
    }
    FloatDPBoundsVector ry = dhmx.value();
    FloatDPBoundsVector rzl = esub(emul(FloatDPBoundsVector(mx-xl),mzl),mmu);
    FloatDPBoundsVector rzu = esub(emul(FloatDPBoundsVector(xu-mx),mzu),mmu);
    CONCLOG_PRINTLN("rx="<<rx<<" ry="<<ry<<" rzl="<<rzl<<" rzu="<<rzu);

    FloatDPBoundsMatrix H(n,n,dp);
    for(SizeType i=0; i!=m; ++i) { H += y[i] * ddhx[i].hessian(); }
    for(SizeType j=0; j!=n; ++j) {
        H[j][j] += mu*rec(sqr(x[j]-xl[j]));
        H[j][j] += mu*rec(sqr(xu[j]-x[j]));
    }

    // S = A Hinv AT
    // H dx + AT dy = rx; A dx = ry;
    //  dx = Hinv ( rx - AT dy )
    //  dy = Sinv ( A Hinv rx - ry )
    FloatDPBoundsMatrix Hinv=inverse(H);
    CONCLOG_PRINTLN_AT(1,"H="<<H<<" Hinv="<<Hinv);
    FloatDPBoundsMatrix S=A*Hinv*transpose(A);
    FloatDPBoundsMatrix Sinv=inverse(S);
    CONCLOG_PRINTLN_AT(1,"S="<<S<<" Sinv="<<Sinv);
    FloatDPBoundsVector dy = Sinv * ( A*(Hinv*rx) - ry );
    FloatDPBoundsVector dx = Hinv * ( rx - transpose(A) * dy);
    CONCLOG_PRINTLN("dx="<<dx<<" dy="<<dy);

    FloatDPBoundsVector nx = x-dx;
    FloatDPBoundsVector ny = y-dy;
    CONCLOG_PRINTLN("nx="<<nx<<" ny="<<ny);
    CONCLOG_PRINTLN_AT(1,"h(x)="<<h(nx));

    x = refinement(x,nx); y=refinement(y,ny);
    FloatDPBounds nmu = zero;
    for(SizeType i=0; i!=m; ++i) { nmu += sqr(y[i]); }
    mu=refinement(mu,nmu);
}


//------- Optimality condition functions -----------------------------------//

/*

struct KuhnTuckerFunctionBody : VectorMultivariateFunctionMixin<KuhnTuckerFunctionBody,ExactIntervalType>
{
    ValidatedScalarMultivariateFunction f;
    Array<ValidatedScalarMultivariateFunction> g;
    Array<ValidatedScalarMultivariateFunction> df;
    Array<Array<ValidatedScalarMultivariateFunction> > dg;

    KuhnTuckerFunctionBody(ValidatedScalarMultivariateFunction _f, ValidatedVectorMultivariateFunction _g) {
        ARIADNE_ASSERT(_f.argument_size()==_g.argument_size());
        const SizeType m=_g.argument_size();
        const SizeType n=_g.result_size();
        g.resize(n); df.resize(m); dg.resize(n); for(SizeType j=0; j!=n; ++j) { dg[j].resize(m); }
        f=_f;
        for(SizeType j=0; j!=n; ++j) { g[j]=_g[j]; }
        for(SizeType i=0; i!=m; ++i) { df[i]=f.derivative(i); }
        for(SizeType j=0; j!=n; ++j) { for(SizeType i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    SizeType result_size() const { return g.size()*2+f.argument_size(); }
    SizeType argument_size() const { return g.size()*2+f.argument_size(); }
    ValidatedScalarMultivariateFunction operator[](SizeType) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream&) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> Void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const SizeType m=f.argument_size();
        const SizeType n=g.size();
        Vector<X> x(project(arg,range(0,n)));
        Vector<X> y(project(arg,range(n,n+m)));
        Vector<X> z(project(arg,range(n+m,n+m+n)));
        Vector<X> rx(m), rz(n), rs(n);
        for(SizeType i=0; i!=m; ++i) { rx[i]=df[i].evaluate(y); for(SizeType j=0; j!=n; ++j) { rx[i]=rx[i]-x[j]*dg[j][i].evaluate(y); } }
        for(SizeType j=0; j!=n; ++j) { rz[j]=g[j].evaluate(y) + z[j]; }
        for(SizeType j=0; j!=n; ++j) { rs[j]=x[j]*z[j]; }
        project(res,range(0,n))=rz;
        project(res,range(n,n+m))=rx;
        project(res,range(n+m,n+m+n))=rs;
    }
};

struct FeasibilityKuhnTuckerFunctionBody : VectorMultivariateFunctionMixin<FeasibilityKuhnTuckerFunctionBody,ExactIntervalType>
{
    Array<ValidatedScalarMultivariateFunction> g;
    Array<Array<ValidatedScalarMultivariateFunction> > dg;

    FeasibilityKuhnTuckerFunctionBody(ValidatedVectorMultivariateFunction _g) {
        const SizeType m=_g.argument_size();
        const SizeType n=_g.result_size();
        g.resize(n); dg.resize(n); for(SizeType j=0; j!=n; ++j) { dg[j].resize(m); }
        for(SizeType j=0; j!=n; ++j) { g[j]=_g[j]; for(SizeType i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    SizeType result_size() const { return g.size()*2+g[0].argument_size()+1; }
    SizeType argument_size() const { return g.size()*2+g[0].argument_size()+1; }
    ValidatedScalarMultivariateFunction operator[](SizeType) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream&) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> Void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const SizeType m=g[0].argument_size();
        const SizeType n=g.size();
        Vector<X> x(project(arg,range(0,n)));
        Vector<X> y(project(arg,range(n,n+m)));
        Vector<X> z(project(arg,range(n+m,n+m+n)));
        X t(arg[n+m+n]);
        Vector<X> rx(m), rz(n), rs(n); X rt;
        for(SizeType i=0; i!=m; ++i) { rx[i]=x[0]*dg[0][i].evaluate(y); for(SizeType j=1; j!=n; ++j) { rx[i]=rx[i]+x[j]*dg[j][i].evaluate(y); } }
        for(SizeType j=0; j!=n; ++j) { rz[j]=g[j].evaluate(y) + t + z[j]; }
        for(SizeType j=0; j!=n; ++j) { rs[j]=x[j]*z[j]; }
        rt=1-x[0]; for(SizeType j=1; j!=n; ++j) { rt=rt-x[j]; }
        project(res,range(0,n))=rz;
        project(res,range(n,n+m))=rx;
        project(res,range(n+m,n+m+n))=rs;
        res[n+m+n]=rt;
    }
};



struct ConstrainedFeasibilityKuhnTuckerFunctionBody : VectorMultivariateFunctionMixin<FeasibilityKuhnTuckerFunctionBody,ExactIntervalType>
{
    SizeType m;
    SizeType n;
    ExactIntervalVectorType d;
    Array<ValidatedScalarMultivariateFunction> g;
    ExactIntervalVectorType c;
    Array<Array<ValidatedScalarMultivariateFunction> > dg;

    ConstrainedFeasibilityKuhnTuckerFunctionBody(ExactBoxType D, ValidatedVectorMultivariateFunction _g, ExactBoxType C) {
        m=_g.argument_size();
        n=_g.result_size();
        d=D; c=C;
        g.resize(n); dg.resize(n); for(SizeType j=0; j!=n; ++j) { dg[j].resize(m); }
        for(SizeType j=0; j!=n; ++j) { g[j]=_g[j]; for(SizeType i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    SizeType result_size() const { return 5*m+4*n+1u; }
    SizeType argument_size() const { return 5*m+4*n+1u; }
    ValidatedScalarMultivariateFunction operator[](SizeType) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream& os) const { return os << "KuhnTuckerFunctionBody"; }

    template<class X> Void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const X zero=arg[0].zero_element();
        const SizeType l=2*(m+n);
        assert(arg.size()==l+m+l+1);
        Vector<X> x(project(arg,range(0u,l)));
        Vector<X> y(project(arg,range(l,l+m)));
        Vector<X> z(project(arg,range(l+m,l+m+l)));
        X t(arg[l+m+l]);
        Vector<X> rx(m,zero), rz(l,zero), rs(l,zero); X rt(zero);
        Vector<X> gy(n);
        for(SizeType j=0; j!=n; ++j) { gy[j]=g[j].evaluate(y); }
        Matrix<X> dgy(n,m);
        for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { dgy[j][i]=dg[j][i].evaluate(y); } }

        for(SizeType i=0; i!=m; ++i) {
            for(SizeType j=0; j!=n; ++j) { rx[i]+=x[j]*(dgy[j][i]-c[j].upper_bound()); rx[i]+=x[n+j]*(c[j].lower_bound()-dgy[j][i]); }
            rx[i]+=x[2*n+i]*(y[i]-d[i].upper_bound())-x[2*n+m+i]*(d[i].lower_bound()-y[i]);
        }
        for(SizeType j=0; j!=n; ++j) { rz[j]=gy[j] + t + z[j]; rz[n+j]=t+z[n+j]-gy[j]; }
        for(SizeType i=0; i!=m; ++i) { rz[2*n+i]=y[i]+t+z[2*n+i]; rz[2*n+m+i]=y[i]+t+z[2*n+m+i]; }
        for(SizeType k=0; k!=l; ++k) { rs[k]=x[k]*z[k]; }
        rt+=1.0; for(SizeType j=0; j!=2*n; ++j) { rt=rt-x[j]; }
        project(res,range(0,l))=rz;
        project(res,range(l,l+m))=rx;
        project(res,range(l+m,l+m+l))=rs;
        res[l+m+l]=rt;
    }
};

*/


//------- KrawczykOptimiser -----------------------------------//

/*

ValidatedVector KrawczykOptimiser::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


ValidatedKleenean KrawczykOptimiser::
feasible(ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c) const
{
    CONCLOG_PRINTLN("KrawczykOptimiser::feasible(ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c)");
    CONCLOG_PRINTLN("d="<<d<<", g="<<g<<", c="<<c);

    ARIADNE_ASSERT(g.argument_size()==d.size());
    ARIADNE_ASSERT(g.result_size()==c.size());

    ExactIntervalType t; ExactIntervalVectorType x,y,z;
    setup_feasibility(d,g,c,x,y,z,t);


    // FIXME: Allow more steps
    for(SizeType i=0; i!=12; ++i) {
        CONCLOG_PRINTLN_AT(1,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", x="<<x<<", z="<<z);
        try {
            this->feasibility_step(d,g,c,x,y,z,t);
        }
        catch(const SingularMatrixException& e) {
            return indeterminate;
        }
        if(t.lower_bound()>t.upper_bound()) {
            CONCLOG_PRINTLN_AT(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
            return indeterminate;
        }
        if(t.lower_bound()>0.0) {
            CONCLOG_PRINTLN_AT(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
            if(this->is_feasible_point(d,g,c,midpoint(y))) {
                return true;
            }
        }
        if(t.upper_bound()<0.0) {
            CONCLOG_PRINTLN_AT(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
            return false;
        }
    }
    CONCLOG_PRINTLN_AT(1,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
    if(this->is_infeasibility(d,g,c,midpoint(x))) {
        return false;
    }
    return indeterminate;
}



Void KrawczykOptimiser::setup_feasibility(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                          ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z, ExactIntervalType& t) const
{
    const SizeType m=g.argument_size();
    const SizeType n=g.result_size();
    const SizeType l=2*(m+n);
    x=ExactIntervalVectorType(l, ExactIntervalType(0,1)/l);
    y=d;
    z.resize(2*(m+n));
    compute_tz(d,g,c,y,t,z);
}


Void KrawczykOptimiser::compute_tz(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                   const ExactIntervalVectorType& y, ExactIntervalType& t, ExactIntervalVectorType& z) const
{
    ARIADNE_ASSERT(d.size()>0u);
    //static const double EPS=1.0/8;
    static const float min_float=std::numeric_limits<float>::min();

    const SizeType m=g.argument_size();
    const SizeType n=g.result_size();

    // Compute the image of y under the constraint function
    ExactIntervalVectorType gy=g(y);
    gy+=ExactIntervalVectorType(gy.size(),ExactIntervalType(-min_float,+min_float));
    ExactIntervalVectorType my=midpoint(y);
    ExactIntervalVectorType mgy=g(my);

    // Find the range of possible values of the optimal t
    // This range is too pessimistic
    t=ExactIntervalType(+inf,+inf);
    for(SizeType j=0; j!=n; ++j) {
        t=min(t,c[j]-gy[j]);
        t=min(t,gy[j]-c[j]);
    }
    for(SizeType i=0; i!=m; ++i) {
        t=min(t,d[i]-y[i]);
        t=min(t,y[i]-d[i]);
    }

    // Find the range of possible values of the optimal t
    FloatDP tmin=+inf;
    FloatDP tmax=+inf;
    for(SizeType j=0; j!=n; ++j) {
        tmax=min(tmax,sub(up,c[j].upper_bound(),gy[j].lower_bound()));
        tmax=min(tmax,sub(up,gy[j].upper_bound(),c[j].lower_bound()));
        tmin=min(tmin,sub(down,c[j].upper_bound(),mgy[j].upper_bound()));
        tmin=min(tmin,sub(down,mgy[j].lower_bound(),c[j].lower_bound()));
    }
    for(SizeType i=0; i!=m; ++i) {
        tmin=min(tmin,sub(up,d[i].upper_bound(),y[i].lower_bound()));
        tmax=min(tmax,sub(up,y[i].upper_bound(),d[i].lower_bound()));
        tmin=min(tmin,sub(down,d[i].upper_bound(),my[i].upper_bound()));
        tmax=min(tmax,sub(down,my[i].lower_bound(),d[i].lower_bound()));
    }
    tmin-=0.0625;
    t=ExactIntervalType(tmin,tmax);


    // Find the range of possible values of the optimal z
    // This range is too pessimistic
    for(SizeType j=0; j!=n; ++j) {
        z[j]=max(c[j].upper_bound()-gy[j]-t,0.0);
        z[n+j]=max(gy[j]-c[j].lower_bound()-t,0.0);
    }
    for(SizeType i=0; i!=m; ++i) {
        z[2*n+i]=max(d[i].upper_bound()-y[i]-t,0.0);
        z[2*n+m+i]=max(y[i]-d[i].lower_bound()-t,0.0);
    }

    // Find the range of possible values of the optimal z
    // This range is too pessimistic
    for(SizeType j=0; j!=n; ++j) {
        z[j]=ExactIntervalType(0.0,c[j].upper_bound()-mgy[j].lower_bound()-tmin);
        z[n+j]=ExactIntervalType(0.0,mgy[j].upper_bound()-c[j].lower_bound()-tmin);
    }
    for(SizeType i=0; i!=m; ++i) {
        z[2*n+i]=ExactIntervalType(0.0,d[i].upper_bound()-my[i].lower_bound()-tmin);
        z[2*n+m+i]=ExactIntervalType(0.0,my[i].upper_bound()-d[i].lower_bound()-tmin);
    }

    CONCLOG_PRINTLN_AT(9,"  d="<<d<<", c="<<c<<", y="<<y<<", g(y)="<<gy<<", t="<<t<<", z="<<z);

}


Void KrawczykOptimiser::
minimisation_step(const ValidatedScalarMultivariateFunction& f, const ValidatedVectorMultivariateFunction& g,
                  ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z) const
{
    const SizeType m=f.argument_size();
    const SizeType n=g.result_size();

    Differential<UpperIntervalType> ddf=f.evaluate(Differential<UpperIntervalType>::variables(2,y));
    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));

    ExactIntervalMatrixType H(m,m);
    set_hessian(H,ddf);
    for(SizeType j=0; j!=n; ++j) { add_hessian(H,-x[j],ddg[j]); }

    ExactIntervalMatrixType A(m,n);
    set_jacobian_transpose(A,ddg);

    CONCLOG_PRINTLN_AT(9,"f="<<f<<",g="<<g);
    CONCLOG_PRINTLN_AT(9,"x="<<x<<" y="<<y<<" z="<<z);
    CONCLOG_PRINTLN_AT(9,"A="<<A);
    CONCLOG_PRINTLN_AT(9,"H="<<H);

    ARIADNE_NOT_IMPLEMENTED;

}



Void KrawczykOptimiser::feasibility_step(const ValidatedVectorMultivariateFunction& g,
                                         ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z, ExactIntervalType& t) const
{
    ARIADNE_NOT_IMPLEMENTED;
    const SizeType m=y.size();
    const SizeType n=x.size();

    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));

    // A is the transpose derivative matrix aij=dgj/dyi
    ExactIntervalMatrixType A(m,n);
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    CONCLOG_PRINTLN_AT(9,"A="<<A);

    // H is the Hessian matrix Hik = xj*dgj/dyidyk
    ExactIntervalMatrixType H(m,m);
    for(SizeType j=0; j!=n; ++j) {
        add_hessian(H,x[j],ddg[j]);
    }
    CONCLOG_PRINTLN_AT(9," H="<<H);

    FloatDPMatrix mA=midpoint(A);
    CONCLOG_PRINTLN_AT(9," mA="<<mA);
    FloatDPMatrix mH=midpoint(H);
    CONCLOG_PRINTLN_AT(9," mH="<<mH);

    RawFloatDPVector mD(n);
    for(SizeType j=0; j!=n; ++j) { mD[j]=midpoint(x[j])/midpoint(z[j]); }
    CONCLOG_PRINTLN_AT(9," mD="<<mD);

    FloatDPMatrix& mS=mH;
    adat(mS,mA,mD);
    CONCLOG_PRINTLN_AT(9,"mS="<<mS);
    FloatDPMatrix mSinv=inverse(mS);
    CONCLOG_PRINTLN_AT(9,"mSinv="<<mSinv);
}

// Feasibility step for dual (inequality constrained) problem without using slack variables
// FIXME: Do we need a slackness parameter mu? Probably not; hopefully the infinities are kept in check...
// This method has the advantage of not needing to update the primal variables
Void KrawczykOptimiser::feasibility_step(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                         ExactIntervalVectorType& y, ExactIntervalType& t) const
{
    const SizeType m=d.size();
    const SizeType n=c.size();

    // Compute function values
    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));

    // gy is the vector of values of g(y)
    ExactIntervalVectorType gy(n);
    for(SizeType j=0; j!=n; ++j) { gy[j]=ddg[j].value(); }

    // z is the vector of slack variables z[k]=cu[k]-gy[k]-t or z[k]=gy[k]-cl[k]-t
    ExactIntervalVectorType z(2*(m+n));
    for(SizeType j=0; j!=n; ++j) { z[j]=d[j].upper_bound()-gy[j]-t; z[n+j]=gy[j]-d[j].lower_bound()-t; }
    for(SizeType i=0; i!=m; ++i) { z[i]=c[2*n+i].upper_bound()-y[i]-t; z[2*n+m+i]=y[i]-c[i].lower_bound()-t; }

    ExactIntervalVectorType zr(2*(m+n));
    for(SizeType k=0; k!=2*(m+n); ++k) { zr[k]=1.0/z[k]; }

    ExactIntervalVectorType D(2*(m+n));
    for(SizeType k=0; k!=2*(m+n); ++k) { D[k]=zr[k]*zr[k]; }

    // A is the transpose derivative matrix aij=dgj/dyi
    ExactIntervalMatrixType A(m,n);
    for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { A[i][j]=ddg[j][i]; } }

    // A is the sum of scaled Hessian matrices hi1i2=zj*ddgj/dyi1yi2
    ExactIntervalMatrixType H(m,m);

    ExactIntervalMatrixType SE(m+1,m+1);
    // SE[0:m][0:m] is the matrix H+/-A(D1+D2)AT+(D3+D4) where D=z.z
    for(SizeType i1=0; i1!=m; ++i1) { for(SizeType i2=0; i2!=m; ++i2) { SE[i1][i2]=H[i1][i2];
        for(SizeType j=0; j!=n; ++j) { SE[i1][i2]+=A[i1][j]*(D[j]+D[n+j])*A[i2][j]; }
    } }
    for(SizeType i=0; i!=m; ++i) { SE[i][i]+=(D[2*n+i]+D[2*n+m+i]); }
    // SE[m][0:m]=SE[0:m][m] is the vector A(D1-D2)e+(D3-D4)e
    for(SizeType i=0; i!=m; ++i) { SE[i][m]=D[2*n+i]-D[2*n+m+i];
        for(SizeType j=0; j!=n; ++j) { SE[i][m]+=A[i][j]*(D[j]-D[n+j]); }
        SE[m][i]=SE[i][m];
    }
    // SE[m][m] is the scalar eT(D1+D2)e+eT(D3+D4)e
    for(SizeType k=0; k!=2*(m+n); ++k) { SE[m][m]+=D[k]; }

    // Vector of residuals
    ExactIntervalVectorType re(m+1);
    for(SizeType i=0; i!=m; ++i) { re[i]+=(zr[2*n+i]-zr[2*n+m+i]);
        for(SizeType j=0; j!=n; ++j) { re[i]+=A[i][j]*(zr[j]-zr[n+j]); }
    }
    for(SizeType j=0; j!=n; ++j) { re[m]+=(zr[j]+zr[n+n]); }
    for(SizeType i=0; i!=m; ++i) { re[m]+=(zr[2*n+i]+zr[2*n+m+i]); }

    // Compute inverse Jacobian matrix
    ExactIntervalMatrixType JE;
    try {
        JE=inverse(midpoint(SE));
    }
    catch(const SingularMatrixException& e) {
        ARIADNE_WARN("Matrix S="<<midpoint(SE)<<" is not invertible");
        CONCLOG_PRINTLN_AT(1,"WARNING: Matrix S="<<midpoint(SE)<<" is not invertible");
        throw e;
    }

    // Krawczyk step
    ExactIntervalVectorType dyt=prod(JE,ExactIntervalVectorType(midpoint(re)))+prod(ExactIntervalMatrixType::identity(m+1)-prod(JE,SE),re-midpoint(re));

    // Extract y and t
    ExactIntervalVectorType yt=join(y,t);
    ExactIntervalVectorType nyt=yt+dyt;

    yt=intersection(yt,nyt);
    y=project(yt,range(0,m));
    t=yt[m];
}




Void KrawczykOptimiser::feasibility_step(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                         ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z, ExactIntervalType& t) const
{
    const SizeType m=d.size();
    const SizeType n=c.size();
    const SizeType o=2*(m+n);

    ARIADNE_ASSERT_MSG(g.argument_size()==m,"d="<<d<<" g="<<g);
    ARIADNE_ASSERT_MSG(g.result_size()==n,"d="<<d<<" g="<<g<<" c="<<c);
    ARIADNE_ASSERT(x.size()==o);
    ARIADNE_ASSERT(y.size()==m);
    ARIADNE_ASSERT(z.size()==o);

    ExactIntervalVectorType yt=join(y,t);
    CONCLOG_PRINTLN_AT(9,"m="<<m<<" n="<<n);
    CONCLOG_PRINTLN_AT(9,"x="<<x<<" yt="<<yt<<" z="<<z);

    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));
    CONCLOG_PRINTLN_AT(9,"  ddg="<<ddg);

    // gy is the vector of values of g(y)
    UpperIntervalVectorType gy(n); for(SizeType j=0; j!=n; ++j) { gy[j]=ddg[j].value(); }
    CONCLOG_PRINTLN_AT(9,"  g(y)="<<gy<<" ");

    // A is the transpose derivative matrix aij=dgj/dyi, extended with a column of ones
    UpperIntervalMatrixType A(m,n);
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    CONCLOG_PRINTLN_AT(9," A="<<A<<" ");

    // H is the Hessian matrix Hik = (xcuj-xclj)*dgj/dyidyk
    UpperIntervalMatrixType H(m,m);
    for(SizeType j=0; j!=n; ++j) {
        add_hessian(H,x[j]-x[n+j],ddg[j]);
    }
    CONCLOG_PRINTLN_AT(9," H="<<H);

    // Construct the extended valuation GY=(gy-cu+te,cl-gy+te,y-bu+te,bl-y+te)
    UpperIntervalVectorType gye(o);
    for(SizeType j=0; j!=n; ++j) { gye[j]=gy[j]-c[j].upper_bound()+t; gye[n+j]=c[j].lower_bound()-gy[j]+t; }
    for(SizeType i=0; i!=m; ++i) { gye[2*n+i]=y[i]-d[i].upper_bound()+t; gye[2*n+m+i]=d[i].lower_bound()-y[i]+t; }
    CONCLOG_PRINTLN_AT(9,"  GE="<<gye);

    // Construct the extended matrix AE=(A -A I -I \\ e e 0 0)
    UpperIntervalMatrixType AE(m+1,o);
    for(SizeType i=0; i!=m; ++i) { for(SizeType j=0; j!=n; ++j) { AE[i][j]=A[i][j]; AE[i][n+j]=-A[i][j]; } }
    for(SizeType i=0; i!=m; ++i) { AE[i][2*n+i]=1; AE[i][2*n+m+i]=-1; }
    for(SizeType k=0; k!=o; ++k) { AE[m][k]=1; }
    UpperIntervalMatrixType AET=transpose(AE);

    FloatDPMatrix mA=midpoint(A);
    FloatDPMatrix mAE=midpoint(AE);
    FloatDPMatrix mAET=midpoint(AET);
    FloatDPMatrix mH=midpoint(H);
    RawFloatDPVector mx=midpoint(x);
    RawFloatDPVector myt=midpoint(yt);
    RawFloatDPVector mz=midpoint(z);
    RawFloatDPVector mDE=ediv(mx,mz);


    // Construct the symmetric matrix and its inverse
    //FloatDPMatrix S(m+1,m+1); adat(S,AE,DE);
    //CONCLOG_PRINTLN_AT(9,"S="<<S);
    //S=FloatDPMatrix(m+1,m+1); simple_adat(S,AE,DE);
    //CONCLOG_PRINTLN_AT(9,"S="<<S);
    FloatDPMatrix mS=feasibility_adat(mH,mA,mDE);
    CONCLOG_PRINTLN_AT(9,"mS="<<mS);
    FloatDPMatrix mSinv=inverse(mS);
    CONCLOG_PRINTLN_AT(9,"mSinv="<<mSinv);

    // FIXME: What if S is not invertible?

    // Construct the residuals
    UpperIntervalVectorType rx=emul(mx,mz);
    //RawFloatDPVector ryt=-prod(AE,x); ryt[m]+=1; // FIXME: Need hessian
    UpperIntervalVectorType ryt=-feasibility_mul(mA,mx); ryt[m]+=1; // FIXME: Need hessian
    UpperIntervalVectorType rz=midpoint(gye)+mz;
    CONCLOG_PRINTLN_AT(9,"rx="<<rx<<" ryt="<<ryt<<" rz="<<rz);

    // Construct the errors on the residuals ([M]-M)([x]-x)
    UpperIntervalVectorType ex=x-mx;
    UpperIntervalVectorType eyt=yt-myt;
    UpperIntervalVectorType ez=z-mz;
    UpperIntervalMatrixType eA=A-mA;
    UpperIntervalMatrixType eH=H-mH;

    UpperIntervalVectorType erx=2.0*emul(ex,ez);
    UpperIntervalVectorType eryt=UpperIntervalMatrixType(AE-mAE)*ex;
    UpperIntervalVectorType erz=UpperIntervalMatrixType(AET-mAET)*eyt;
    CONCLOG_PRINTLN_AT(9,"erx="<<erx<<" eryt="<<eryt<<" erz="<<erz);

    rx+=2.0*emul(ex,ez);
    ryt+=UpperIntervalMatrixType(AE-mAE)*ex;
    rz+=UpperIntervalMatrixType(AET-mAET)*eyt;
    CONCLOG_PRINTLN_AT(9,"rx="<<rx<<" ryt="<<ryt<<" rz="<<rz);

    //RawFloatDPVector rr=prod(AE,ediv(RawFloatDPVector(rx-emul(x,rz)),z))-ryt;

    // Compute the error differences
    UpperIntervalVectorType erxdz=ediv(erx,mz);
    UpperIntervalVectorType edyt=(mSinv*mAE)*erxdz + mSinv*eyt - (mSinv*(mAE*DiagonalMatrix<FloatDP>(mDE))) * ez;
    UpperIntervalVectorType edz=-erz-feasibility_trmul(mA,edyt);
    UpperIntervalVectorType edx=-ediv(UpperIntervalVectorType(erx+emul(mx,edz)),mz);
    CONCLOG_PRINTLN_AT(9,"edx="<<edx<<" edyt="<<edyt<<" edz="<<edz);

    // Compute the error differences
    UpperIntervalVectorType eerr=prod(mAE,ediv(esub(erx,emul(mx,erz)),mz))-eryt;
    CONCLOG_PRINTLN_AT(9,"  eerr="<<eerr);
    UpperIntervalVectorType eedyt=prod(mSinv,eerr);
    UpperIntervalVectorType eedz=-erz-feasibility_trmul(mA,eedyt);
    UpperIntervalVectorType eedx=-ediv(UpperIntervalVectorType(erx+emul(mx,eedz)),mz);
    CONCLOG_PRINTLN_AT(9,"eedx="<<eedx<<" eedyt="<<eedyt<<" eedz="<<eedz);


    // Compute the differences
    UpperIntervalVectorType rr=prod(mAE,ediv(esub(rx,emul(mx,rz)),mz))-ryt;
    UpperIntervalVectorType dyt=prod(mSinv,rr);
    UpperIntervalVectorType dz=-rz-feasibility_trmul(mA,dyt);
    UpperIntervalVectorType dx=-ediv(UpperIntervalVectorType(rx+emul(mx,dz)),mz);
    CONCLOG_PRINTLN_AT(9,"dx="<<dx<<" dyt="<<dyt<<" dz="<<dz);

    UpperIntervalVectorType nx,ny,nyt,nz; FloatDP nt;
    nx=mx+dx;
    nyt=myt+dyt;
    nz=mz+dz;

    x=intersection(x,nx);
    yt=intersection(yt,nyt);
    z=intersection(z,nz);

    y=project(yt,range(0,m));
    t=yt[m];
}

*/


} // namespace Ariadne
