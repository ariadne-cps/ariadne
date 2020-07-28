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

#include "../function/functional.hpp"
#include "../config.hpp"

#include <limits>

#include "../utility/macros.hpp"
#include "../output/logging.hpp"
#include "../utility/tuple.hpp"
#include "../utility/tribool.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/linear_algebra.decl.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/diagonal_matrix.hpp"
#include "../algebra/differential.hpp"
#include "../algebra/algebra.hpp"
#include "../function/function.hpp"
#include "../function/function_mixin.hpp"
#include "../function/taylor_function.hpp"
#include "../function/formula.hpp"
#include "../function/procedure.hpp"

#include "../solvers/nonlinear_programming.hpp"
#include "../solvers/solver.hpp"
#include "../algebra/multi_index-noaliasing.hpp"
#include "../solvers/constraint_solver.hpp"

#include "../algebra/expansion.inl.hpp"

namespace Ariadne {

inline Sweeper<FloatDP> default_sweeper() { return Sweeper<FloatDP>(); }

typedef DiagonalMatrix<FloatDPBounds> FloatDPBoundsDiagonalMatrix;

typedef Vector<FloatDPApproximation> FloatDPApproximationVector;
typedef VectorRange<FloatDPApproximationVector> FloatDPApproximationVectorRange;
typedef Matrix<FloatDPApproximation> FloatDPApproximationMatrix;
typedef DiagonalMatrix<FloatDPApproximation> FloatDPApproximationDiagonalMatrix;
typedef Differential<FloatDPApproximation> FloatDPApproximationDifferential;

typedef Vector<FloatDPBounds> FloatDPBoundsVector;
typedef Matrix<FloatDPBounds> FloatDPBoundsMatrix;
typedef Differential<FloatDPBounds> FloatDPBoundsDifferential;
typedef Vector<FloatDPValue> ExactFloatDPVectorType;

typedef Vector<UpperIntervalType> UpperIntervalVectorType;
typedef Matrix<UpperIntervalType> UpperIntervalMatrixType;

typedef FloatDPApproximation ApproximateNumericType;

Matrix<ApproximateNumericType> join(Matrix<ApproximateNumericType> const&, Matrix<ApproximateNumericType> const&, Matrix<ApproximateNumericType> const&);

inline Vector<Differential<RawFloatDP>>const& cast_raw(Vector<Differential<FloatDPApproximation>>const& v) {
    return reinterpret_cast<Vector<Differential<RawFloatDP>>const&>(v);
}

inline Vector<FloatDPApproximation>& cast_approximate(Vector<RawFloatDP>& v) {
    return reinterpret_cast<Vector<FloatDPApproximation>&>(v);
}

inline Vector<FloatDPApproximation>& cast_approximate(Vector<FloatDPApproximation>& v) {
    return v;
}
inline Vector<Differential<FloatDPApproximation>>const& cast_approximate(Vector<Differential<RawFloatDP>>const& v) {
    return reinterpret_cast<Vector<Differential<FloatDPApproximation>>const&>(v);
}

inline UpperIntervalType dot(Vector<UpperIntervalType> const& bx1, Vector<ExactIntervalType> const& bx2) {
    return dot(bx1,Vector<UpperIntervalType>(bx2));
}

template<class X1, class IVL2> inline decltype(auto) dot(Vector<X1> const& x1, Box<IVL2> const& bx2) {
    return dot(x1,cast_vector(bx2)); }

template<class X> inline
DiagonalMatrix<X> diagonal_matrix(const Vector<X>& v) {
    return DiagonalMatrix<X>(v.array());
}

template<class X> inline
Bool epos(const Vector<X>& x) {
    for(Nat i=0; i!=x.size(); ++i) { if(x[i]<=0) { return false; } } return true;
}

template<class X> inline
Bool eneg(const Vector<X>& x) {
    for(Nat i=0; i!=x.size(); ++i) { if(x[i]>=0) { return false; } } return true;
}

template<class X, class XX> inline
Bool egtr(const Vector<X>& x, const XX& s) {
    for(Nat i=0; i!=x.size(); ++i) { if(decide(x[i]<=s)) { return false; } } return true;
}

template<class X, class XX> inline
Bool elss(const Vector<X>& x, const XX& s) {
    for(Nat i=0; i!=x.size(); ++i) { if(decide(x[i]>=s)) { return false; } } return true;
}

template<class X> inline
Vector<X> eadd(const Vector<X>& x, const Vector<X>& y) {
    Vector<X> r(x.size(),dp); for(Nat i=0; i!=r.size(); ++i) { r[i]=x[i]+y[i]; } return r;
}

template<class X> inline
Vector<X> esub(const Vector<X>& x, const X& s) {
    Vector<X> r(x.size(),dp); for(Nat i=0; i!=r.size(); ++i) { r[i]=x[i]-s; } return r;
}

template<class X> inline
Vector<X> esub(const Vector<X>& x, const Vector<X>& y) {
    Vector<X> r(x.size(),dp); for(Nat i=0; i!=r.size(); ++i) { r[i]=x[i]-y[i]; } return r;
}

template<class X> inline
Vector<X> emul(const Vector<X>& x, const Vector<X>& z) {
    Vector<X> r(x.size(),dp); for(Nat i=0; i!=r.size(); ++i) { r[i]=x[i]*z[i]; } return r;
}

template<class X, class XX> inline
Vector<X> ediv(const Vector<X>& x, const Vector<XX>& z) {
    Vector<X> r(x.size(),dp); for(Nat i=0; i!=r.size(); ++i) { r[i]=x[i]/z[i]; } return r;
}

template<class X> inline
Vector<X> ediv(const X& s, const Vector<X>& z) {
    Vector<X> r(z.size(),dp); for(Nat i=0; i!=r.size(); ++i) { r[i]=s/z[i]; } return r;
}

template<class X> inline
Vector<X> erec(const Vector<X>& z) {
    Vector<X> r(z.size(),dp); for(Nat i=0; i!=r.size(); ++i) { r[i]=rec(z[i]); } return r;
}

template<class X> inline
Vector<X> esqr(const Vector<X>& z) {
    Vector<X> r(z.size(),dp); for(Nat i=0; i!=r.size(); ++i) { r[i]=sqr(z[i]); } return r;
}

inline
ExactIntervalType eivl(const RawFloatDPVector& x) {
    ARIADNE_ASSERT(x.size()>0); ExactIntervalType r=ExactIntervalType(FloatDPValue(x[0]));
    for(Nat i=1; i!=x.size(); ++i) { r=hull(r,FloatDPValue(x[i])); } return r;
}

Matrix<ApproximateNumericType> join(Matrix<ApproximateNumericType> const& A1, Matrix<ApproximateNumericType> const& A2, Matrix<ApproximateNumericType> const& A3) {
    Nat m=A1.row_size(); Nat n1=A1.column_size(); Nat n2=A2.column_size(); Nat n3=A3.column_size();
    Matrix<ApproximateNumericType> A123(m,n1+n2+n3,(A1.zero_element()+A2.zero_element()+A3.zero_element()));
    project(A123,range(0,m),range(0,n1))=A1;
    project(A123,range(0,m),range(n1,n1+n2))=A2;
    project(A123,range(0,m),range(n1+n2,n1+n2+n3))=A3;
    return A123;
}

template<class X> Matrix<X> cojoin(Matrix<X> const& A1, Matrix<X> const& A2, Matrix<X> const& A3) {
    Nat n=A1.column_size(); Nat m1=A1.row_size(); Nat m2=A2.row_size(); Nat m3=A3.row_size();
    Matrix<X> A123(m1+m2+m3,n,(A1.zero_element()+A2.zero_element()+A3.zero_element()));
    project(A123,range(0,m1),range(0,n))=A1;
    project(A123,range(m1,m1+m2),range(0,n))=A2;
    project(A123,range(m1+m2,m1+m2+m3),range(0,n))=A3;
    return A123;
}


// Compute S+=ADA^T, where D is diagonal and S is symmetric.
template<class X>
Void adat(Matrix<X>& S, const Matrix<X>& A, const Vector<X>& D)
{
    const Nat m=A.row_size();
    const Nat n=A.column_size();
    for(Nat i1=0; i1!=m; ++i1) {
        for(Nat j=0; j!=n; ++j) {
            X ADij=A[i1][j]*D[j];
            for(Nat i2=i1; i2!=m; ++i2) {
                S[i1][i2]+=ADij*A[i2][j];
            }
        }
    }
    for(Nat i1=1; i1!=m; ++i1) {
        for(Nat i2=0; i2!=i1; ++i2) {
            S[i1][i2]=S[i2][i1];
        }
    }
}

// Compute S+=A^TDA, where D is diagonal and S is symmetric.
template<class X>
Void atda(Matrix<X>& S, const Matrix<X>& A, const Vector<X>& D)
{
    assert(S.row_size()==S.column_size());
    assert(S.column_size()==A.column_size());
    assert(D.size()==A.row_size());

    const Nat m=A.column_size();
    const Nat n=A.row_size();
    for(Nat i1=0; i1!=m; ++i1) {
        for(Nat j=0; j!=n; ++j) {
            X ATDij=A[j][i1]*D[j];
            for(Nat i2=i1; i2!=m; ++i2) {
                S[i1][i2]+=ATDij*A[j][i2];
            }
        }
    }
    for(Nat i1=1; i1!=m; ++i1) {
        for(Nat i2=0; i2!=i1; ++i2) {
            S[i1][i2]=S[i2][i1];
        }
    }
}

// Compute S+=A^TDA, where D is diagonal and S is symmetric.
template<class X>
Void atda(Matrix<X>& S, const Matrix<X>& A, const DiagonalMatrix<X>& D)
{
    atda(S,A,D.diagonal());
}

// Compute S=ADA^T, where D is diagonal.
template<class X>
Matrix<X> adat(const Matrix<X>& A, const Vector<X>& D)
{
    const Nat m=A.row_size();
    Matrix<X> S=Matrix<X>::zero(m,m);
    adat(S,A,D);
    return S;
}

// Compute S+=AA^T
template<class X>
Matrix<X> amulat(const Matrix<X>& A)
{
    const Nat m=A.row_size();
    const Nat n=A.column_size();
    Matrix<X> S(m,m);
    for(Nat i1=0; i1!=m; ++i1) {
        for(Nat j=0; j!=n; ++j) {
            for(Nat i2=i1; i2!=m; ++i2) {
                S[i1][i2]+=A[i1][j]*A[i2][j];
            }
        }
    }
    for(Nat i1=1; i1!=m; ++i1) {
        for(Nat i2=0; i2!=i1; ++i2) {
            S[i1][i2]=S[i2][i1];
        }
    }
    return S;
}

template<class X> inline Bool all_greater(const Vector<X>& x, const X& e) {
    for(Nat i=0; i!=x.size(); ++i) { if(x[i]<=e) { return false; } } return true;
}






template<class X> Vector< Differential<X> > second_derivative(const ValidatedVectorMultivariateFunction& f, const Vector<X>& x) {
    Vector< Differential<X> > d=Differential<X>::variables(f.result_size(),f.argument_size(),2);
    return f.evaluate(d);
}

template<class Vec, class Diff> Void set_gradient(Vec& g, const Diff& D) {
    Nat i=0;
    typename Diff::ConstIterator iter=D.begin();
    if(iter!=D.end() && iter->index().degree()==0) { ++iter; }
    while(iter!=D.end() && iter->index().degree()<=2) {
        while(iter->index()[i]==0) { ++i; }
        g[i]=iter->coefficient();
        ++iter;
    }
}

template<class Mx, class Diff> Void set_jacobian_transpose(Mx& A, const Vector<Diff>& D) {
    for(Nat j=0; j!=A.column_size(); ++j) {
        for(Nat i=0; i!=A.row_size(); ++i) {
            A[i][j]=D[j][i];
        }
    }
}

template<class Mx, class Diff> Void set_hessian(Mx& H, const Diff& D) {
    typedef typename Diff::ValueType X;
    Nat i=0; Nat j=1;
    typename Diff::ConstIterator iter=D.begin();
    while(iter!=D.end() && iter->index().degree()<=1) { ++iter; }
    while(iter!=D.end() && iter->index().degree()<=2) {
        UniformConstReference<MultiIndex> a=iter->index();
        UniformConstReference<X> c=iter->coefficient();
        while(a[i]==0) { ++i; j=i+1; }
        if(a[i]==2) { H[i][i]=c; }
        else { while(a[j]==0) { ++j; } H[i][j]=c; H[j][i]=c; }
        ++iter;
    }
}

template<class Mx, class S, class Diff> Void add_hessian(Mx& H, const S& s, const Diff& D) {
    typedef typename Diff::ValueType X;
    typename Diff::ConstIterator iter=D.begin();
    while(iter!=D.end() && iter->index().degree()<=1) { ++iter; }
    while(iter!=D.end() && iter->index().degree()==2) {
        UniformConstReference<MultiIndex> a=iter->index();
        UniformConstReference<X> c=iter->coefficient();
        Nat i=0;
        while(a[i]==0) { ++i; }
        if(a[i]==2) { H[i][i]+=s*c; }
        else { Nat j=i+1; while(a[j]==0) { ++j; } H[i][j]+=s*c; H[j][i]+=s*c; }
        ++iter;
    }
}

// Compute the product (A -A I -I ; 1 1 1 1) v
template<class X, class XX> Vector<X> feasibility_mul(const Matrix<XX>& A, const Vector<X>& v)
{
    const Nat m=A.row_size();
    const Nat n=A.column_size();
    ARIADNE_ASSERT(v.size()==2*(m+n));
    Vector<X> r(m+1u,v.zero_element());
    for(Nat i=0; i!=m; ++i) {
        r[i]=v[2*n+i]-v[2*n+m+i];
        for(Nat j=0; j!=n; ++j) {
            r[i]+=A[i][j]*(v[j]-v[n+j]);
        }
    }
    for(Nat k=0; k!=2*(m+n); ++k) {
        r[m]+=v[k];
    }
    return r;
}

// Compute the product (AT 1 \\ -AT 1 \\ I 1 \\ -I 1) I -I ; 1 1 1 1) v
template<class X, class XX> Vector<X> feasibility_trmul(const Matrix<XX>& A, const Vector<X>& w)
{
    const Nat m=A.row_size();
    const Nat n=A.column_size();
    ARIADNE_ASSERT(w.size()==m+1);
    Vector<X> r(2*(m+n),w.zero_element());
    for(Nat j=0; j!=n; ++j) {
        r[j]=0;
        for(Nat i=0; i!=m; ++i) {
            r[j]+=A[i][j]*w[i];
        }
        r[n+j]=-r[j];
        r[j]+=w[m];
        r[n+j]+=w[m];
    }
    for(Nat i=0; i!=m; ++i) {
        r[2*n+i]=w[i]+w[m];
        r[2*n+m+i]=-w[i]+w[m];
    }
    return r;
}


// Compute the product \f$\hat{A}^T \hat{D} \hat{A} + \hat{H}\f$ where \f$\hat{A}=\left(\begin{matrix}A&-A&I&-I\\1&1&1&1\end{matrix}\right)\f$ and \f$\hat{D}=D\f$ is diagonal.
template<class X> Matrix<X> feasibility_adat(const Matrix<X>& H, const Matrix<X>& A, const Vector<X>& D)
{
    const Nat m=A.row_size();
    const Nat n=A.column_size();
    ARIADNE_ASSERT(H.row_size()==m);
    ARIADNE_ASSERT(H.column_size()==m);
    ARIADNE_ASSERT(D.size()==2*(m+n));
    Matrix<X> S(m+1,m+1,H.zero_element());

    for(Nat i=0; i!=m; ++i) { for(Nat j=0; j!=m; ++j) { S[i][j] = H[i][j]; } }
    for(Nat i=0; i!=m; ++i) { S[i][m]=0; S[m][i]=0; } S[m][m]=0;

    for(Nat i1=0; i1!=m; ++i1) {
        for(Nat j=0; j!=n; ++j) {
            X ADij=A[i1][j]*(D[j]+D[n+j]);
            for(Nat i2=0; i2!=m; ++i2) {
                S[i1][i2]+=ADij*A[i2][j];
            }
        }
    }
    for(Nat i=0; i!=m; ++i) {
        S[i][i]+=(D[2*n+i]+D[2*n+m+i]);
    }
    for(Nat i=0; i!=m; ++i) {
        for(Nat j=0; j!=n; ++j) {
            S[i][m]+=A[i][j]*(D[j]-D[n+j]);
        }
        S[i][m]+=(D[2*n+i]-D[2*n+m+i]);
        S[m][i]=S[i][m];
    }
    for(Nat k=0; k!=2*(m+n); ++k) {
        S[m][m]+=D[k];
    }

    return S;
}






template<class R>
class ConstrainedFeasibilityMatrix {
    ConstrainedFeasibilityMatrix(const Vector<R>& x, const Vector<R>& z, const Matrix<R>& a, const Matrix<R>& h)
        : X(x), Z(z), D(ediv(x,z)), A(a), H(h) { }

    template<class RR> Tuple< Vector<RR>,Vector<RR>,Vector<RR> >
    mul(const Vector<RR>& x, const Vector<RR>& yt, const Vector<RR>& z) const {
        Vector<RR> nx=Z*x+X*x;
        Vector<RR> nyt=H*yt-A*x;
        Vector<RR> nz=H*yt-A*x;
        return make_tuple(nx,nyt,nz);
    }

    template<class RR> Tuple< Vector<RR>,Vector<RR>,Vector<RR> >
    solve(const Vector<RR>& x, const Vector<RR>& yt, const Vector<RR>& z) const {
        Sinv=inverse(feasibility_adat(H,A,D));
        Vector<RR> rx=Z.solve(x);
        Vector<RR> ryt=yt+feasibility_mul(A,rx-D*z);
        Vector<RR> rz=z;
        ryt=Sinv*ryt;
        rz=rz-feasibility_trmul(A,ryt);
        rx=rz-D*rz;
        return make_tuple(rx,ryt,rz);
    }

    DiagonalMatrix<R> X;
    DiagonalMatrix<R> Z;
    DiagonalMatrix<R> D;
    const Matrix<R>& A;
    const Matrix<R>& H;
    Matrix<R> Sinv;
};


inline ExactBoxType widen(ExactBoxType bx, RawFloatDP e) {
    for(Nat i=0; i!=bx.size(); ++i) {
        bx[i]=ExactIntervalType(sub(down,bx[i].lower().raw(),e),add(up,bx[i].upper().raw(),e));
    }
    return bx;
}


const FloatDPValue OptimiserBase::zero = FloatDPValue(0,dp);
const FloatDPValue OptimiserBase::one = FloatDPValue(1,dp);

Bool OptimiserBase::
almost_feasible_point(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C, ApproximateVectorType ax, FloatDPApproximation error) const
{
    ExactVectorType ex=cast_exact(ax);
    if(!contains(D,ex)) { return false; }
    ApproximateVectorType gx=g(ax);
    return probably(contains(widen(C,cast_exact(error)),gx));
}


Bool OptimiserBase::
is_feasible_point(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C, ExactVectorType x) const
{
    if(!contains(D,x)) { return false; }
    Vector<FloatDPBounds> gx=g(x);
    return definitely(contains(C,gx));
}


ValidatedKleenean OptimiserBase::
contains_feasible_point(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C, ValidatedVectorType X) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C<<", X="<<X);

    // Now test if the (reduced) box X satisfies other constraints
    if(definitely(disjoint(Box<UpperIntervalType>(X),D))) { return false; }
    if(definitely(not subset(Box<UpperIntervalType>(X),D))) { return indeterminate; }

    // Test inequality constraints
    ValidatedKleenean result = true;
    Vector<FloatDPBounds> gx=g(X);
    ARIADNE_LOG_PRINTLN("g(X)="<<gx);
    for(Nat i=0; i!=C.size(); ++i) {
        if(definitely(disjoint(UpperIntervalType(gx[i]),C[i]))) {
            return false;
        }
        if(!C[i].is_singleton()) {
            if(definitely(not subset(UpperIntervalType(gx[i]),C[i]))) { result = indeterminate; }
        }
    }

    // Break if some inequality constraints indefinite
    if(!definitely(result)) { return result; }

    // Extract the equality constraints
    List<Nat> equality_constraints;
    equality_constraints.reserve(C.size());
    for(Nat i=0; i!=C.size(); ++i) {
        if(C[i].is_singleton()) { equality_constraints.append(i); }
    }

    // Construct the function g_e(x) = g_{e_i}(x)
    ARIADNE_ASSERT(g.result_size()>0);
    ValidatedVectorMultivariateFunction ge(equality_constraints.size(),g.domain());
    ExactBoxType ce(equality_constraints.size());
    for(Nat i=0; i!=ge.result_size(); ++i) {
        ge[i]=g[equality_constraints[i]];
        ce[i]=C[equality_constraints[i]];
    }

    ARIADNE_LOG_PRINTLN("ge="<<ge<<", ce="<<ce);

    // FIXME: Carefully change this code!
    FloatDPBoundsMatrix ivlA=jacobian(ge,X);
    ARIADNE_LOG_PRINTLN("ivlA="<<ivlA);
    FloatDPApproximationVector fltD(X.size(),dp);
    for(Nat i=0; i!=X.size(); ++i) { fltD[i]=rec(sqr(X[i].error())); }
    FloatDPApproximationMatrix fltA=midpoint(ivlA);
    ARIADNE_LOG_PRINTLN("A="<<fltA);
    ARIADNE_LOG_PRINTLN("D="<<fltD);
    FloatDPApproximationMatrix fltL = FloatDPApproximationDiagonalMatrix(fltD.array())*transpose(fltA);
    ARIADNE_LOG_PRINTLN("L="<<fltL);

    FloatDPBoundsMatrix ivlS = ivlA * cast_exact(fltL);
    ARIADNE_LOG_PRINTLN("ivlS="<<ivlS);

    FloatDPBoundsMatrix ivlR = inverse(ivlS);
    try {
        ivlR=inverse(ivlS);
    }
    catch (SingularMatrixException e) {
        return indeterminate;
    }

    ARIADNE_LOG_PRINTLN("ivlR="<<ivlR);
    FloatDPBoundsMatrix& valR=reinterpret_cast<FloatDPBoundsMatrix&>(ivlR);

    // Projected interval Newton step. For h:R^n->R^m; Dh mxn, take L nxm.
    // ExactIntervalType Newton update X' = x - L * (Dh(X)*L)^{-1} * h(x)
    // Choose L = rad(X)^2 Dh(x)^T where rad(X) is the diagonal matrix of radii of X
    Vector<FloatDPBounds> x=midpoint(X);
    Vector<FloatDPBounds> new_X = x - cast_exact(fltL) * (valR * (ge(x)-cast_singleton(ce)) );
    ARIADNE_LOG_PRINTLN("old_X="<<X);
    ARIADNE_LOG_PRINTLN("new_X="<<new_X);
    Vector<FloatDPBounds> reduced_X = refinement(X,new_X);
    ARIADNE_LOG_PRINTLN("reduced_X="<<reduced_X);

    if(refines(new_X,X)) { return true; }
    else { return indeterminate; }
}




Bool OptimiserBase::
validate_feasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                     ExactVectorType x0, ExactVectorType y0) const
{
    return this->validate_feasibility(D,g,C,x0);
}

Bool OptimiserBase::
validate_feasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                     ExactVectorType x0) const
{
    ARIADNE_PRECONDITION(D.size()==g.argument_size());
    ARIADNE_PRECONDITION(C.size()==g.result_size());
    ARIADNE_PRECONDITION(x0.size()==D.size());
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C);
    ARIADNE_LOG_PRINTLN("x0="<<x0);

    Vector<FloatDPBounds> x(x0);
    ARIADNE_LOG_PRINTLN("x="<<x);

    Vector<FloatDPBounds> gx=g(x);
    ARIADNE_LOG_PRINTLN("gx="<<gx);

    List<Nat> equalities, inequalities;
    for(Nat i=0; i!=C.size(); ++i) {
        if(C[i].lower()==C[i].upper()) {
            equalities.append(i);
        } else {
            inequalities.append(i);
            if(!definitely(contains(C[i],gx[i]))) {
                ARIADNE_LOG_PRINTLN_AT(1,"g["<<i<<"](x)="<<gx[i]<<", C["<<i<<"]="<<C[i]);
                return false; }
        }
    }

    if(equalities.empty()) { ARIADNE_LOG_PRINTLN("Feasible"); return true; }

    Nat k=equalities.size();
    Nat n=D.size();
    ValidatedVectorMultivariateFunction h(equalities.size(),g.domain());
    ExactFloatDPVectorType c(equalities.size(),dp);
    for(Nat i=0; i!=equalities.size(); ++i) {
        h[i] = g[equalities[i]];
        c[i] = C[equalities[i]].lower();
    }
    ARIADNE_LOG_PRINTLN("h="<<h<<" c="<<c<<" h(x)-c="<<(h(x0)-c));

    // Attempt to solve h(x0+AT*w)=0
    // TODO: Change to use validated numbers
    Matrix<FloatDPBounds> AT = transpose(h.jacobian(x0));
    ARIADNE_LOG_PRINTLN("A="<<transpose(AT));
    Vector<FloatDPBounds> w0(k,FloatDPBounds(0,dp));

    Bool found_solution=false;
    Bool validated_solution=false;

    Vector<FloatDPBounds> w(k,dp), mw(k,dp), nw(k,dp);
    Vector<FloatDPBounds> mx(n,dp);

    for(Nat ii=0; ii!=12; ++ii) {
        mw=midpoint(w);
        x=x0+AT*w;
        mx=x0+AT*mw;
        nw = mw - solve(h.jacobian(x)*AT,Vector<FloatDPBounds>(h(mx)-c));
        ARIADNE_LOG_PRINTLN_AT(1,"w="<<w<<", h(x0+AT*w)="<<h(x)<<", nw="<<nw<<", refines="<<refines(nw,w));

        if(!found_solution) {
            if(refines(nw,w)) {
                found_solution=true;
                w=w+FloatDPBounds(0,1,dp)*(w0-w);
            } else {
                w=w+FloatDPBounds(0,1,dp)*(FloatDPBounds(2,dp)*nw-w);
            }
        } else {
            if(refines(nw,w)) {
                validated_solution=true;
            } else if(validated_solution) {
                // Not a contraction, so presumably accurate enough
                break;
            }
            w=refinement(nw,w);
        }

    }
    ARIADNE_LOG_PRINTLN("w="<<w<<", validated="<<validated_solution);

    if(!validated_solution) { return false; }

    // Compute x value
    x=x0+AT*w;
    ARIADNE_LOG_PRINTLN("x="<<x);
    gx=g(x);
    ARIADNE_LOG_PRINTLN("g(x)="<<gx);

    // Check that equality constraints are plausible
    ARIADNE_DEBUG_ASSERT(models(h(x)-c,ExactFloatDPVectorType(k,dp)));

    // Check inequality constraints once more
    for(Nat i=0; i!=C.size(); ++i) {
        if(C[i].lower()==C[i].upper()) {
            ARIADNE_DEBUG_ASSERT(models(gx[i],C[i].midpoint()));
        } else {
            if(!definitely(element(gx[i],C[i]))) {
                return false;
            }
        }
    }
    return true;
}


Bool OptimiserBase::
validate_infeasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                       ExactVectorType x, ExactVectorType y) const
{
    ARIADNE_PRECONDITION(D.size()==g.argument_size());
    ARIADNE_PRECONDITION(C.size()==g.result_size());
    ARIADNE_PRECONDITION(x.size()==D.size());
    ARIADNE_PRECONDITION(y.size()==C.size());
    ARIADNE_LOG_SCOPE_CREATE;
    // Compute first-order approximation to g(D) centred at x.
    // For feasibilty, have yg(D) cap yC nonempty.
    // Estimate y g(X) = y g(x) + y Dg(X).(X-x)

    // Compute y.C
    UpperIntervalType yC = dot(y,cast_vector(C));

    // Compute Taylor estimate of y g(X)
    ValidatedVectorMultivariateTaylorFunctionModelDP tg(D,g,default_sweeper());
    ValidatedScalarMultivariateTaylorFunctionModelDP tyg(D,default_sweeper());
    for(Nat j=0; j!=y.size(); ++j) { tyg += y[j]*tg[j]; }
    UpperIntervalType tygD = apply(tyg,D);

    UpperIntervalMatrixType dgD = jacobian_range(g,cast_vector(D));
    UpperIntervalVectorType ydgD = transpose(dgD) * y;

    UpperIntervalType ygx = dot(y,apply(g,UpperIntervalVectorType(x)));

    UpperIntervalType ygD = ygx;
    for(Nat i=0; i!=x.size(); ++i) {
        ygD += ydgD[i] * (D[i]-x[i]);
    }

    ARIADNE_LOG_PRINTLN("yC="<<yC<<" tygD="<<tygD<<" ygD="<<ygD);

    if(definitely(intersection(yC,ygD).is_empty())) { ARIADNE_LOG_PRINTLN("Infeasible"); return true; }
    else { return false; }
}

// FIXME: Look at this code again, especially relating to generalised Lagrange multipliers
Bool OptimiserBase::
is_infeasibility_certificate(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C, ExactFloatDPVectorType y) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C<<", y="<<y);

    if(y.size()==0) { return D.is_empty(); }

    // Try to prove lambda.(g(y)-c) != 0
    const Nat n=C.size();

    ValidatedScalarMultivariateTaylorFunctionModelDP tyg(D,default_sweeper());
    for(Nat i=0; i!=n; ++i) {
        tyg+=y[i]*ValidatedScalarMultivariateTaylorFunctionModelDP(D,g[i],default_sweeper());
    }
    ValidatedNumericType iygx = tyg(cast_singleton(D));

    UpperIntervalType iyC(0,0);
    for(Nat i=0; i!=n; ++i) {
        iyC+=y[i]*C[i];
    }

    if(definitely(disjoint(iyC,UpperIntervalType(iygx)))) {
        return true;
    } else {
        return false;
    }
}




typedef OptimiserBase::ValidatedVectorType ValidatedVectorType;

ValidatedVectorType OptimiserBase::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ValidatedVectorMultivariateFunction gh=join(g,h);
    ExactBoxType C(gh.result_size(),ExactIntervalType(0,0));
    for(Nat i=0; i!=g.result_size(); ++i) { C[i]=ExactIntervalType(-inf,0); }
    return this->minimise(f,D,gh,C);
}



ValidatedKleenean OptimiserBase::
feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ValidatedVectorMultivariateFunction h) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ValidatedVectorMultivariateFunction gh=join(g,h);
    ExactBoxType C(gh.result_size(),ExactIntervalType(0,0));
    for(Nat i=0; i!=g.result_size(); ++i) { C[i]=ExactIntervalType(-inf,0); }
    return this->feasible(D,gh,C);
}


//------- NonlinearInfeasibleInteriorPointOptimiser -------------------------//

struct NonlinearInfeasibleInteriorPointOptimiser::PrimalDualData {
    PrimalDualData() : PrimalDualData(0u,0u,dp) { }
    PrimalDualData(SizeType m, SizeType n, DP pr) : w(m,pr), x(n,pr), y(m,pr) { }
    FloatDPApproximationVector w,x,y;
};

struct NonlinearInfeasibleInteriorPointOptimiser::StepData : public PrimalDualData {
    StepData() : StepData(0u,0u,dp) { }
    StepData(SizeType m, SizeType n, DP pr)
        : PrimalDualData(m,n,pr), vl(m,pr), wl(m,pr), xl(n,pr), zl(n,pr), vu(m,pr), wu(m,pr), xu(n,pr), zu(n,pr), mu(pr) { }
    FloatDPApproximationVector vl,wl,xl,zl,vu,wu,xu,zu; FloatDPApproximation mu;
};

ValidatedVectorType NonlinearInfeasibleInteriorPointOptimiser::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("f="<<f<<", D="<<D<<", g="<<g<<", C="<<C);

    static const double VALUE_TOLERANCE=1e-8;
    static const double STATE_TOLERANCE=1e-8;
    static const Nat MAXIMUM_STEPS=24;

    ARIADNE_ASSERT(f.argument_size()==D.size());
    ARIADNE_ASSERT(g.argument_size()==D.size());
    ARIADNE_ASSERT(g.result_size()==C.size());
    StepData v;
    FloatDPApproximationVector& x=cast_approximate(v.x);
    FloatDPApproximationVector& y=cast_approximate(v.y);
    C=intersection(C,cast_exact_box(apply(g,D)+UpperIntervalVectorType(C.size(),UpperIntervalType(-1,+1))));
    this->setup_feasibility(D,g,C,v);
    FloatDPApproximationVector oldx=x;

    static const float MU_MIN = 1e-12;

    // FIXME: Allow more steps
    for(Nat i=0; i!=MAXIMUM_STEPS; ++i) {
        ARIADNE_LOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        oldx=x;
        FloatDPApproximation oldfx=f(oldx);
        this->step(f,D,g,C,v);
        if(this->is_infeasibility_certificate(D,g,C,cast_exact(y))) {
            ARIADNE_LOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
            ARIADNE_LOG_PRINTLN_AT(1,"Infeasible");
            std::cerr<<"EXCEPTION: "<<InfeasibleProblemException().what()<<"\n";
            throw InfeasibleProblemException();
        }
        FloatDPApproximation fx=f(x);
        if(probably(mag(fx-oldfx)<VALUE_TOLERANCE) && probably(norm(oldx-x)<STATE_TOLERANCE)) {
            break;
        }
        if(v.mu.raw()<MU_MIN) {
            break;
        }
    }
    ARIADNE_LOG_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));

    if(this->validate_feasibility(D,g,C,cast_exact(x))) {
        ARIADNE_LOG_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        return cast_exact(x);
    }
    ARIADNE_LOG_PRINTLN("indeterminate_feasibility");
    throw IndeterminateFeasibilityException();
}

ValidatedKleenean NonlinearInfeasibleInteriorPointOptimiser::
feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const
{
    ARIADNE_LOG_SCOPE_CREATE
    ARIADNE_LOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C);

    ARIADNE_ASSERT(g.argument_size()==D.size());
    ARIADNE_ASSERT(g.result_size()==C.size());

    StepData v;
    FloatDPApproximationVector& x=cast_approximate(v.x);
    FloatDPApproximationVector& y=cast_approximate(v.y);

    ApproximateScalarMultivariateFunction f(EuclideanDomain(D.dimension()));
    ExactBoxType R=intersection(cast_exact_box(widen(apply(g,D),1)),C);
    this->setup_feasibility(D,g,R,v);

    static const float MU_MIN = 1e-12;

    // FIXME: Allow more steps
    for(Nat i=0; i!=12; ++i) {
        ARIADNE_LOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
        this->step(f,D,g,R,v);
        if(this->validate_feasibility(D,g,C,cast_exact(x))) {
            ARIADNE_LOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
            ARIADNE_LOG_PRINTLN("Feasible");
            return true;
        }
        if(this->is_infeasibility_certificate(D,g,C,cast_exact(y))) {
            ARIADNE_LOG_PRINTLN_AT(1,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
            ARIADNE_LOG_PRINTLN("Infeasible");
            return false;
        }
        if(v.mu.raw()<MU_MIN) {
            break;
        }
    }
    ARIADNE_LOG_PRINTLN("f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x));
    ARIADNE_LOG_PRINTLN("Indeterminate");
    return indeterminate;
}

Void NonlinearInfeasibleInteriorPointOptimiser::
setup_feasibility(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
                  StepData& v) const
{
    ExactIntervalType I(-1,+1);
    Nat m=C.size(); Nat n=D.size();

    v.x=midpoint(D);
    v.y=FloatDPApproximationVector(m,zero);
    v.w=midpoint(C);

    //stp.xl=lower(D)-x;
    v.wl=Vector(m,-one);
    v.wu=Vector(m,+one);
    v.xl=lower_bounds(D)-v.x;
    v.xu=upper_bounds(D)-v.x;
    v.vl=Vector(m,-one);
    v.vu=Vector(m,+one);
    v.zl=Vector(n,-one);
    v.zu=Vector(n,+one);
    // FIXME: What should relaxation parameter be?
    v.mu=1.0_x;
}


Void
NonlinearInfeasibleInteriorPointOptimiser::step(
    const ApproximateScalarMultivariateFunction& f, const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c,
    StepData& v) const
{
    FloatDPApproximationVector& w=v.w; FloatDPApproximationVector& x=v.x; FloatDPApproximationVector& y=v.y;
    FloatDPApproximation& mu=v.mu;
    FloatDPApproximationVector& wl=v.wl; FloatDPApproximationVector& wu=v.wu;
    FloatDPApproximationVector& xl=v.xl; FloatDPApproximationVector& xu=v.xu;
    FloatDPApproximationVector& vl=v.vl; FloatDPApproximationVector& vu=v.vu;
    FloatDPApproximationVector& zl=v.zl; FloatDPApproximationVector& zu=v.zu;
    FloatDPApproximationVector cl=lower_bounds(c); FloatDPApproximationVector cu=upper_bounds(c);
    FloatDPApproximationVector dl=lower_bounds(d); FloatDPApproximationVector du=upper_bounds(d);

    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_LOG_PRINTLN("w ="<<w<<",  x ="<<x<<", y ="<<y<<" mu="<<mu);
    ARIADNE_LOG_PRINTLN("wl="<<wl<<", wu="<<wu<<", xl="<<xl<<", xu="<<xu);
    ARIADNE_LOG_PRINTLN("vl="<<vl<<", vu="<<vu<<", zl="<<zl<<", zu="<<zu);
    ARIADNE_LOG_PRINTLN("cl-wl="<<cl-wl<<", dl-xl="<<dl-xl);
    ARIADNE_LOG_PRINTLN("  w  ="<<w<<",   x  ="<<x);
    ARIADNE_LOG_PRINTLN("cu-wu="<<cu-wu<<", du-xu="<<du-xu);
    static const ExactDouble gamma=0.0009765625_x;
    static const ExactDouble sigma=0.125_x;
    static const ExactDouble scale=0.75_x;

    const Nat n=d.size();
    const Nat m=c.size();

    ARIADNE_ASSERT_MSG(f.argument_size()==d.size(),"f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_ASSERT_MSG(g.argument_size()==d.size(),"f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_ASSERT_MSG(g.result_size()==c.size(),  "f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_ASSERT(w.size()==m);
    ARIADNE_ASSERT(x.size()==n);
    ARIADNE_ASSERT(y.size()==m);

    mu = mu * sigma;

    FloatDPApproximationVector ax(x);
    FloatDPApproximationDifferential ddfx=f.differential(ax,2u);
    ARIADNE_LOG_PRINTLN("ddfx="<<ddfx);
    Vector<FloatDPApproximationDifferential> ddgx=g.differential(ax,2u);
    ARIADNE_LOG_PRINTLN("ddgx="<<ddgx);

    FloatDPApproximation fx = ddfx.value();
    Vector<FloatDPApproximation> gx = ddgx.value();
    ARIADNE_LOG_PRINTLN("f(x)="<<fx);
    ARIADNE_LOG_PRINTLN("g(x)="<<gx);
    Vector<FloatDPApproximation> Jfx = transpose(ddfx.gradient());
    Matrix<FloatDPApproximation> A = ddgx.jacobian();
    Matrix<FloatDPApproximation>& Jgx = A;
    ARIADNE_LOG_PRINTLN("Df(x)="<<Jfx);
    ARIADNE_LOG_PRINTLN("Dg(x)="<<Jgx);

    // H is the Hessian matrix H of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) $
    Matrix<FloatDPApproximation> YH = ddfx.hessian();
    for(Nat i=0; i!=m; ++i) {
        YH+=y[i]*ddgx[i].hessian();
    }
    ARIADNE_LOG_PRINTLN("D2f(x)="<<ddfx.hessian());
    ARIADNE_LOG_PRINTLN("D2f(x)+Y.D2g(x)="<<YH);

    // Set up the system of equations
    // (A^TDA + E - Y.H) dx = A^T(r_w-Dr_y)+r_x
    // dw = A \delta x + r_y
    // dy = r_w - D dw

    FloatDPApproximationDiagonalMatrix const& Vl=diagonal_matrix(vl);
    FloatDPApproximationDiagonalMatrix const& Vu=diagonal_matrix(vu);
    FloatDPApproximationDiagonalMatrix const& Wl=diagonal_matrix(wl);
    FloatDPApproximationDiagonalMatrix const& Wu=diagonal_matrix(wu);
    FloatDPApproximationDiagonalMatrix const& Xl=diagonal_matrix(xl);
    FloatDPApproximationDiagonalMatrix const& Xu=diagonal_matrix(xu);
    FloatDPApproximationDiagonalMatrix const& Zl=diagonal_matrix(zl);
    FloatDPApproximationDiagonalMatrix const& Zu=diagonal_matrix(zu);

    // Compute the diagonal matrices
    //   D=XL/ZL+XU/ZU  E=WL/VL+WU/VU
    FloatDPApproximationDiagonalMatrix Dl=Vl/Wl;
    FloatDPApproximationDiagonalMatrix Du=Vu/Wu;
    FloatDPApproximationDiagonalMatrix D=Dl+Du;
    ARIADNE_LOG_PRINTLN("D="<<D);
    FloatDPApproximationDiagonalMatrix El=Zl/Xl;
    FloatDPApproximationDiagonalMatrix Eu=Zu/Xu;
    FloatDPApproximationDiagonalMatrix E=El+Eu;
    ARIADNE_LOG_PRINTLN("E="<<E);

    // normal equation matrix
    FloatDPApproximationMatrix S=YH;
    atda(S,A,D);
    S+=E;

    //FloatDPApproximationMatrix EE(n,n); for(Nat j=0; j!=n; ++j) { EE[j][j]=E[j]; }
    //FloatDPApproximationMatrix DD(m,m); for(Nat i=0; i!=m; ++i) { DD[i][i]=E[i]; }

    ARIADNE_LOG_PRINTLN("S="<<S);
    ARIADNE_DEBUG_ASSERT(decide(norm(FloatDPApproximationMatrix(S-(YH+E+transpose(A)*(D*A))))/norm(S)<1e-8));
    FloatDPApproximationMatrix Sinv=inverse(S);
    ARIADNE_LOG_PRINTLN("Sinv="<<Sinv);

    // Construct the residuals
    // The residual for the slack variable xl is given by the duality condition xl.zl=mu as mu/xl-zl
    // The residual for the dual variable zl is given by the slackness condition x-xl-cl
    // The residual for the auxiliary variable w is given by y-(vu-vl)
    // The residual for the dual variable y is given by g(x)-w
    FloatDPApproximationVector ew=(vl+vu)-y;
    FloatDPApproximationVector ex=Jfx+transpose(Jgx)*y+(zl+zu);
    FloatDPApproximationVector ey=gx-w;
    FloatDPApproximationVector ewl=esub(vl,ediv(mu,wl));
    FloatDPApproximationVector ewu=esub(vu,ediv(mu,wu));
    FloatDPApproximationVector exl=esub(zl,ediv(mu,xl));
    FloatDPApproximationVector exu=esub(zu,ediv(mu,xu));
    FloatDPApproximationVector evl=w+wl-cl;
    FloatDPApproximationVector evu=w+wu-cu;
    FloatDPApproximationVector ezl=x+xl-dl;
    FloatDPApproximationVector ezu=x+xu-du;

    ARIADNE_LOG_PRINTLN("ew="<<ew<<", ex="<<ex<<", ey="<<ey);
    ARIADNE_LOG_PRINTLN("ewl="<<ewl<<", ewu="<<ewu<<", exl="<<exl<<" exu="<<exu);
    ARIADNE_LOG_PRINTLN("evl="<<evl<<", evu="<<evu<<", ezl="<<ezl<<" ezu="<<ezu);

    FloatDPApproximationVector rw = ew - (ewl+ewu) + Dl*evl + Du*evu;
    FloatDPApproximationVector rx = ex - (exl+exu) + El*ezl + Eu*ezu;
    FloatDPApproximationVector& ry = ey;

    // Solve linear system
    // ( D   0  -I ) (dw)   (rw)
    // ( 0  H+E A^T) (dx) = (rx)
    // (-I   A   0 ) (dy) = (ry)

    FloatDPApproximationVector r = transpose(A)*(rw+D*ry)+rx;
    ARIADNE_LOG_PRINTLN("rw="<<rw<<" rx="<<rx<<" ry="<<ry);
    ARIADNE_LOG_PRINTLN("r="<<r);

    // Compute the differences
    FloatDPApproximationVector dx = solve(S,r);
    ARIADNE_LOG_PRINTLN("S*dx="<<S*dx<<" r="<<r);
    ARIADNE_LOG_PRINTLN("S*inverse(S)-I="<<S*inverse(S)-FloatDPApproximationMatrix::identity(n,dp));
    ARIADNE_DEBUG_ASSERT(decide(norm(S*dx - r)/max(1.0_x,norm(r))<1e-4));

    FloatDPApproximationVector dw = A*dx-ry;
    FloatDPApproximationVector dy = D*dw-rw;
    ARIADNE_LOG_PRINTLN("dw="<<dw<<" dx="<<dx<<" dy="<<dy);

    ARIADNE_LOG_PRINTLN("YH*dx+E*dx+dy*A="<<(YH*dx+E*dx+transpose(A)*dy)<<", rx="<<rx);

    // Check solution of linear system for residuals
    ARIADNE_DEBUG_ASSERT(decide(norm(D*dw-dy-rw)/max(one,norm(rw))<1e-4));
    ARIADNE_DEBUG_ASSERT(decide(norm(YH*dx+E*dx+transpose(A)*dy-rx)/max(one,norm(rx))<1e-2));
    ARIADNE_DEBUG_ASSERT(decide(norm(-dw+A*dx-ry)/max(one,norm(ry))<1e-4));

    FloatDPApproximationVector dwl = evl-dw;
    FloatDPApproximationVector dwu = evu-dw;
    FloatDPApproximationVector dxl = ezl-dx;
    FloatDPApproximationVector dxu = ezu-dx;
    FloatDPApproximationVector dvl = ewl-Dl*dwl;
    FloatDPApproximationVector dvu = ewu-Du*dwu;
    FloatDPApproximationVector dzl = exl-El*dxl;
    FloatDPApproximationVector dzu = exu-Eu*dxu;

    ARIADNE_LOG_PRINTLN("dwl="<<dwl<<", dwu="<<dwu<<", dxl="<<dxl<<" dxu="<<dxu);
    ARIADNE_LOG_PRINTLN("dvl="<<dvl<<", dvu="<<dvu<<", dzl="<<dzl<<" dzu="<<dzu);

    ARIADNE_LOG_PRINTLN("YH*dx+dy*A+dzl+dzu="<<(YH*dx+transpose(A)*dy+dzl+dzu)<<", ex="<<ex);
    // Check solution of linear system
/*
    ARIADNE_DEBUG_ASSERT(norm(-dy+dvl+dvu - ew)/max(1.0,norm(ew))<1e-4);
    ARIADNE_DEBUG_ASSERT(norm(YH*dx+dy*A+dzl+dzu - ex)/max(1.0,norm(ex))<1e-2);
    ARIADNE_DEBUG_ASSERT(norm(-dw+A*dx - ey)/max(1.0,norm(ey))<1e-4);
    ARIADNE_DEBUG_ASSERT(norm(Dl*dwl+dvl - ewl)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(Du*dwu+dvu - ewu)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(El*dxl+dzl - exl)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(Eu*dxu+dzu - exu)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(dw+dwl - evl)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(dw+dwu - evu)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(dx+dxl - ezl)<1e-12);
    ARIADNE_DEBUG_ASSERT(norm(dx+dxu - ezu)<1e-12);
*/

    FloatDPApproximationVector nw; FloatDPApproximationVector nx; FloatDPApproximationVector ny;
    FloatDPApproximationVector nwl; FloatDPApproximationVector nwu; FloatDPApproximationVector nxl; FloatDPApproximationVector nxu;
    FloatDPApproximationVector nvl; FloatDPApproximationVector nvu; FloatDPApproximationVector nzl; FloatDPApproximationVector nzu;


    FloatDPApproximation alpha=one;
    nx = x-alpha*dx;
    // Pick an update value which minimises the objective function
    FloatDPApproximation fxmin=f(nx);
    FloatDPApproximation alphamin=one;
    static const Nat REDUCTION_STEPS=4;
    for(Nat i=0; i!=REDUCTION_STEPS; ++i) {
        alpha*=scale;
        nx = x-alpha*dx;
        FloatDPApproximation fnx=f(nx);
        if(decide(fnx<fxmin*scale)) {
            fxmin=fnx;
            alphamin=alpha;
        }
    }
    //alpha=alphamin;
    alpha=one;

    // Since we need to keep the point feasible, but the updates are linear
    // we need to validate feasibility directly.
    static const double MINIMUM_ALPHA=1e-16;
    Bool allfeasible=false;
    while(decide(alpha>MINIMUM_ALPHA) && !allfeasible) {
        nwl=wl-alpha*dwl;
        nwu=wu-alpha*dwu;
        nxl=xl-alpha*dxl;
        nxu=xu-alpha*dxu;
        nvl=vl-alpha*dvl;
        nvu=vu-alpha*dvu;
        nzl=zl-alpha*dzl;
        nzu=zu-alpha*dzu;
        allfeasible = elss(nwl,mu*gamma) && egtr(nwu,mu*gamma) && elss(nxl,mu*gamma) && egtr(nxu,mu*gamma)
                           && elss(nvl,mu*gamma)&& egtr(nvu,mu*gamma) && elss(nzl,mu*gamma) && egtr(nzu,mu*gamma);
        //allfeasible = eneg(nwl) && epos(nwu) && eneg(nxl) && epos(nxu) && elss(nvl,mu*gamma) && egtr(nvu,mu*gamma) && elss(nzl,mu*gamma) && egtr(nzu,mu*gamma);
        if(!allfeasible) { alpha*=scale; }
    }
    nw=w-alpha*dw;
    nx=x-alpha*dx;
    ny=y-alpha*dy;
    if(decide(alpha<=MINIMUM_ALPHA)) {
        ARIADNE_LOG_PRINTLN_AT(1,"w="<<w<<"  x="<<x<<"  y="<<y);
        ARIADNE_LOG_PRINTLN_AT(1,"nw="<<nw<<"  nx="<<nx<<"  ny="<<ny);
        throw NearBoundaryOfFeasibleDomainException(); }
    ARIADNE_LOG_PRINTLN("alpha="<<alpha);
    ARIADNE_LOG_PRINTLN("nw="<<nw<<" nx="<<nx<<" ny="<<ny);
    ARIADNE_LOG_PRINTLN("nwl="<<nwl<<", nwu="<<nwu<<", nxl="<<nxl<<" nxu="<<nxu);
    ARIADNE_LOG_PRINTLN("nvl="<<nvl<<", nvu="<<nvu<<", nzl="<<nzl<<" nzu="<<nzu);

    w=nw; x=nx; y=ny;
    wl=nwl; wu=nwu; xl=nxl; xu=nxu;
    vl=nvl; vu=nvu; zl=nzl; zu=nzu;

    FloatDPApproximation nmu = zero;
    for(Nat i=0; i!=m; ++i) {
        nmu = nmu + wl[i]*vl[i] + wu[i]*vu[i];
    }
    for(Nat j=0; j!=n; ++j) {
        nmu = nmu + xl[j]*zl[j] + xu[j]*zu[j];
    }
    nmu /= (2*(m+n));
    mu = nmu;

    ARIADNE_LOG_PRINTLN("nmu="<<nmu);

}





//------- NonlinearInteriorPointOptimiser -----------------------------------//

ValidatedVectorType NonlinearInteriorPointOptimiser::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("f="<<f<<" D="<<D<<" g="<<g<<" C="<<C);
    ValidatedVectorMultivariateFunction h(0,D.dimension());

    UpperBoxType gD = apply(g,D);
    if(definitely(disjoint(gD,C))) { throw InfeasibleProblemException(); }

    FloatDPApproximationVector x = midpoint(D);
    FloatDPApproximationVector w = midpoint(intersection(UpperBoxType(gD),C));

    FloatDPApproximationVector kappa(g.result_size(),zero);
    FloatDPApproximationVector lambda(h.result_size(),zero);
    FloatDPApproximation mu = one;


    for(Nat i=0; i!=12; ++i) {
        this->minimisation_step(f,D,g,C,h, x,w, kappa,lambda, mu);
        if(i%3==0 && i<=10) { mu *= 0.25_exact; }
    }

    return ValidatedVectorType(cast_exact(x));
}



// See Hande Y. Benson, David F. Shanno, And Robert J. Vanderbei,
// "Interior-point methods for nonconvex nonlinear programming: Jamming and comparative numerical testing"
// For some of the terminology used


// min f(x) | x\in D & w\in C | g(x) = w & h(x) = 0
// Lagrange multipliers kappa d(g(x)-w); lambda dh(x)
Void NonlinearInteriorPointOptimiser::
minimisation_step(const ApproximateScalarMultivariateFunction& f, const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c, const ApproximateVectorMultivariateFunction& h,
                  FloatDPApproximationVector& x, FloatDPApproximationVector& w,
                  FloatDPApproximationVector& kappa, FloatDPApproximationVector& lambda, const FloatDPApproximation& mu) const
{
    const Nat n=x.size();
    const Nat m=kappa.size();
    const Nat l=lambda.size();

    ARIADNE_DEBUG_PRECONDITION(w.size()==kappa.size());
    ARIADNE_DEBUG_PRECONDITION(f.argument_size()==n);
    ARIADNE_DEBUG_PRECONDITION(g.argument_size()==n);
    ARIADNE_DEBUG_PRECONDITION(h.argument_size()==n);
    ARIADNE_DEBUG_PRECONDITION(g.result_size()==m);
    ARIADNE_DEBUG_PRECONDITION(h.result_size()==l);
    ARIADNE_DEBUG_PRECONDITION(contains(d,cast_exact(x)));
    ARIADNE_DEBUG_PRECONDITION(contains(c,cast_exact(w)));
    ARIADNE_DEBUG_PRECONDITION(mu.raw()>0);

    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("x="<<x);
    ARIADNE_LOG_PRINTLN("w="<<w);
    ARIADNE_LOG_PRINTLN_AT(1,"kappa="<<kappa);
    ARIADNE_LOG_PRINTLN_AT(1,"lambda="<<lambda);
    ARIADNE_LOG_PRINTLN_AT(1,"mu="<<mu);

    FloatDPApproximationVector slack(2*n,dp);
    FloatDPApproximationVectorRange slackl(slack,range(0,n));
    FloatDPApproximationVectorRange slacku(slack,range(n,2*n));

    FloatDPApproximationDifferential ddfx=f.evaluate(FloatDPApproximationDifferential::variables(2,x));
    Vector<FloatDPApproximationDifferential> ddgx=g.evaluate(FloatDPApproximationDifferential::variables(2,x));
    Vector<FloatDPApproximationDifferential> ddhx=h.evaluate(FloatDPApproximationDifferential::variables(2,x));

    // G is the constraint value vector
    FloatDPApproximation fx = ddfx.value();
    FloatDPApproximationVector gx = ddgx.value();
    FloatDPApproximationVector hx = ddhx.value();
    ARIADNE_LOG_PRINTLN("f(x)="<<fx);
    ARIADNE_LOG_PRINTLN("g(x)="<<gx);
    ARIADNE_LOG_PRINTLN("h(x)="<<hx);
    ARIADNE_LOG_PRINTLN_AT(1,"g(x)-w="<<(gx-w));

    // A, B are the derivative matrices aij=dgi/dxj
    // HACK: Need to explicitly set size of Jacobian if g or h have result_size of zero
    FloatDPApproximationVector df = transpose(ddfx.gradient());
    ARIADNE_LOG_PRINTLN_AT(1,"df(x)="<<df);
    FloatDPApproximationMatrix A = ddgx.jacobian();
    if(m==0) { A=FloatDPApproximationMatrix(m,n,dp); }
    ARIADNE_LOG_PRINTLN("A="<<A);
    FloatDPApproximationMatrix B = ddhx.jacobian();
    if(l==0) { B=FloatDPApproximationMatrix(l,n,dp); }
    ARIADNE_LOG_PRINTLN("B="<<B);



    // H is the Hessian matrix H[i1,i2] = df/dx[i1]dx[i2] + Sum_[j]kappa[j]*dg[j]/dx[i1]dx[i2] + Sum[k]lambda[k]*dh[k]/dx[i1]dx[i2]
    FloatDPApproximationMatrix H = ddfx.hessian();
    for(Nat j=0; j!=m; ++j) { H += kappa[j] * ddgx[j].hessian(); }
    for(Nat k=0; k!=l; ++k) { H += lambda[k] * ddhx[k].hessian(); }
    ARIADNE_LOG_PRINTLN("H="<<H);

    // Determines the weighting to give to the relaxation parameter mu
    // for equality constraints relative to other constraints
    static const double EQUALITY_RELAXATION_MULTIPLIER = 1.0;

    // Compute the residuals and contributions from slack in x and w
    //   rx = df/dx[i] + Sum[j] dg[j]/dx[i] * kappa[j] + Sum[k] dh[k]/dx[i] * lambda[j] + mu *( 1/(xu[i]-x[i]) - 1/(x[i]-xl[i]) )
    FloatDPApproximationVector rx = df + transpose(A) * kappa + transpose(B)* lambda;
    FloatDPApproximationDiagonalMatrix D(n,dp);
    for(Nat i=0; i!=n; ++i) {
        FloatDPApproximation nuu = rec(d[i].upper()-x[i]);
        FloatDPApproximation nul = rec(x[i]-d[i].lower());
        rx[i] += mu * ( nuu - nul );
        D[i] = mu * ( nuu*nuu + nul*nul );
    }

    //   rw = - kappa[j] + mu *( 1/(wu[i]-w[i]) - 1/(w[i]-wl[i]) )
    FloatDPApproximationVector rw = -kappa;
    FloatDPApproximationDiagonalMatrix C(m,dp);
    for(Nat j=0; j!=m; ++j) {
        FloatDPApproximation nuu = rec(c[j].upper()-w[j]);
        FloatDPApproximation nul = rec(w[j]-c[j].lower());
        rw[j] += (mu*EQUALITY_RELAXATION_MULTIPLIER) * ( nuu - nul );
        C[j] = (mu*EQUALITY_RELAXATION_MULTIPLIER) * ( nuu*nuu + nul*nul );
    }

    //   rkappa = g(x) - w
    FloatDPApproximationVector rkappa = gx - w;

    //   rlambda = h(x)
    FloatDPApproximationVector const& rlambda = hx;

    ARIADNE_LOG_PRINTLN("rx="<<rx);
    ARIADNE_LOG_PRINTLN("rw="<<rw);
    ARIADNE_LOG_PRINTLN("rkappa="<<rkappa);
    ARIADNE_LOG_PRINTLN("rlambda="<<rlambda);

    // Solve the equations
    //   H+D dx        + AT dk + BT dl = rx
    //            C dw -  I dk         = rw
    //    A  dx - I dw                 = rk
    //    B  dx                        = rl

    // Eliminate dw, dk without fill-in to obtain
    //   (H+D+ATCA) dx + BT dl = rx + AT rw + ATC rk
    //         B    dx         = rl

    // Set S=(H+D+ATCA); invert, and eliminate dx
    //   dx = Sinv * (rx + AT rw + ATC rk - BT dl)
    //   (B * Sinv * BT) dl = B * Sinv * (rx + AT rw + ATC rk) - rl
    FloatDPApproximationMatrix& S=H;
    S+=D;
    S+=FloatDPApproximationMatrix(transpose(A))*C*A;
    ARIADNE_LOG_PRINTLN("S="<<S);

    FloatDPApproximationMatrix Sinv=inverse(S);
    ARIADNE_LOG_PRINTLN("R=Sinv="<<Sinv);

    FloatDPApproximationMatrix BSinvBT = (B*Sinv)*FloatDPApproximationMatrix( transpose(B) );
    ARIADNE_LOG_PRINTLN("B*inverse(S)*BT="<<BSinvBT);
    ARIADNE_LOG_PRINTLN("inverse(B*inverse(S)*BT)="<<inverse(BSinvBT));

    FloatDPApproximationVector rr = Sinv * (rx + transpose(A) * (rkappa * C + rw));
    FloatDPApproximationVector dlambda = inverse(BSinvBT) * (B * rr - rlambda);
    FloatDPApproximationVector dx = rr - transpose(B*Sinv) * dlambda;
    FloatDPApproximationVector dw = A * dx - rkappa;
    FloatDPApproximationVector dkappa = rw - C * dw;

    static const FloatDPApproximation ALPHA_SCALE_FACTOR = 0.75_approx;
    static const FloatDPApproximation MINIMUM_ALPHA = 1e-16_approx;

    // Compute distance to move variables preserving feasibility
    // FIXME: Current implementation might fail due to getting too close to boundary!
    FloatDPApproximationVector newx(n,dp);
    FloatDPApproximationVector neww(m,dp);
    FloatDPApproximation alpha = 1.0_approx;
    Bool success = false;
    do {
        newx = x - alpha * dx;
        neww = w - alpha * dw;
        if (contains(d,cast_exact(newx)) && contains(c,cast_exact(neww))) { success = true; }
        else { alpha *= ALPHA_SCALE_FACTOR; }
        if(probably(alpha<MINIMUM_ALPHA)) { throw NearBoundaryOfFeasibleDomainException(); }
    } while(!success);
    ARIADNE_LOG_PRINTLN("alpha="<<alpha);

    FloatDPApproximationVector newlambda = lambda - alpha * dlambda;
    FloatDPApproximationVector newkappa = kappa - alpha * dkappa;

    ARIADNE_LOG_PRINTLN("newx="<<newx);
    ARIADNE_LOG_PRINTLN("neww="<<neww);
    ARIADNE_LOG_PRINTLN("newkappa="<<newkappa);
    ARIADNE_LOG_PRINTLN("newlambda="<<newlambda);

    x=newx; w=neww; kappa=newkappa; lambda=newlambda;
}



ValidatedKleenean NonlinearInteriorPointOptimiser::
feasible(ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("d="<<d<<", g="<<g<<", c="<<c);

    ARIADNE_ASSERT(g.argument_size()==d.size());
    ARIADNE_ASSERT(g.result_size()==c.size());
    FloatDPApproximation t(dp);
    FloatDPApproximationVector x,y,z;

    this->setup_feasibility(d,g,c,x,y);

    // FIXME: Allow more steps
    for(Nat i=0; i!=12; ++i) {
        ARIADNE_LOG_PRINTLN_AT(1,"t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", x="<<x<<", z="<<z);
        this->feasibility_step(d,g,c,x,y);
        if(probably(LogicalValue(t>0))) {
            ARIADNE_LOG_PRINTLN_AT(1,"y="<<y<<", g(y)="<<g(y));
            if(this->is_feasible_point(d,g,c,cast_exact(y))) {
                return true;
            }
        }
    }
    ARIADNE_LOG_PRINTLN("t="<<t<<", y="<<y<<", g(y)="<<g(y));
    if(this->is_infeasibility_certificate(d,g,c,cast_exact(x))) {
        return false;
    }
    return indeterminate;
}


Void
NonlinearInteriorPointOptimiser::feasibility_step(
    const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c,
    FloatDPApproximationVector& x, FloatDPApproximationVector& y) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


Void
NonlinearInteriorPointOptimiser::feasibility_step(
    const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c,
    FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximation& t) const
{
    ARIADNE_LOG_SCOPE_CREATE;

    static const ExactDouble GAMMA=0.0009765625_x; // 1.0/1024;
    static const ExactDouble SIGMA=0.125_x;
    static const ExactDouble SCALE=0.75_x;

    const Nat m=d.size();
    const Nat n=c.size();

    FloatDPApproximationVector z(n,dp);

    ARIADNE_ASSERT_MSG(g.argument_size()==m,"d="<<d<<" g="<<g);
    ARIADNE_ASSERT_MSG(g.result_size()==n,"d="<<d<<" g="<<g<<" c="<<c);
    ARIADNE_ASSERT(x.size()==m);
    ARIADNE_ASSERT(y.size()==n);

    Vector<FloatDPApproximationDifferential> ddgx=g.evaluate(FloatDPApproximationDifferential::variables(2,x));
    ARIADNE_LOG_PRINTLN("ddgx="<<ddgx);

    Vector<FloatDPApproximation> gx = ddgx.value();
    ARIADNE_LOG_PRINTLN("g(x)="<<gx<<" ");
    Matrix<FloatDPApproximation> A = transpose(ddgx.jacobian());
    ARIADNE_LOG_PRINTLN("A="<<A<<" ");

    // H is the Hessian matrix H of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) \lambda_k$
    Matrix<FloatDPApproximation> H(m,m,dp);
    for(Nat i=0; i!=m; ++i) {
        H+=y[i]*ddgx[i].hessian();
    }
    ARIADNE_LOG_PRINTLN("H="<<H<<" ");


    // Add correction for singleton domain to diagonal elements of Hessian
    for(Nat i=0; i!=m; ++i) {
    }

    // Compute diagonal entries of KKT Hessian
    Vector<FloatDPApproximation> D(n,dp);
    for(Nat j=0; j!=n; ++j) {
        if(c[j].lower()==c[j].upper()) {
        } else if(c[j].upper().raw()==+inf) {
        } else if(c[j].lower().raw()==-inf) {
        } else {
            ARIADNE_DEBUG_ASSERT(definitely(-infty<c[j].lower() && c[j].lower()<c[j].upper() && c[j].upper()<+infty));
        }
    }

    FloatDPApproximation sigma(SIGMA,dp);
    FloatDPApproximation mu=dot(x,z)/m;
    if(!egtr(emul(x,z),GAMMA*mu)) {
        ARIADNE_LOG_PRINTLN("WARNING: near-degeneracy in Lyapunov multipliers in interior-point solver:\n  x="<<x<<", y="<<y<<", z="<<z);
        x=(1-sigma)*x+FloatDPApproximationVector(x.size(),sigma/x.size());
        mu=dot(x,z)/m;
    }

    FloatDPApproximationVector yt=join(y,t);
    ARIADNE_LOG_PRINTLN("m="<<m<<" n="<<n);
    ARIADNE_LOG_PRINTLN("x="<<x<<" yt="<<yt<<" z="<<z);


    // Construct diagonal matrices
    FloatDPApproximationVector DE=ediv(x,z);
    ARIADNE_LOG_PRINTLN("D="<<DE);

    // Construct the extended valuation GY=(gy-cu+te,cl-gy+te,y-bu+te,bl-y+te)
    FloatDPApproximationVector gye(2*(m+n),dp);
    //for(Nat j=0; j!=n; ++j) { gxe[j]=gy[j]-c[j].upper()+t; gye[n+j]=c[j].lower()-gy[j]+t; }
    //for(Nat i=0; i!=m; ++i) { gye[2*n+i]=y[i]-d[i].upper()+t; gye[2*n+m+i]=d[i].lower()-y[i]+t; }
    ARIADNE_LOG_PRINTLN("GE="<<gye);

    // Construct the extended matrix AE=(A -A I -I \\ e e 0 0)
    FloatDPApproximationMatrix AE(m+1,2*(m+n),dp);
    //for(Nat i=0; i!=m; ++i) { for(Nat j=0; j!=n; ++j) { AE[i][j]=A[i][j]; AE[i][n+j]=-A[i][j]; } }
    //for(Nat i=0; i!=m; ++i) { AE[i][2*n+i]=1; AE[i][2*n+m+i]=-1; }
    //for(Nat k=0; k!=o; ++k) { AE[m][k]=1; }
    FloatDPApproximationMatrix AET=transpose(AE);

    // Construct the symmetric matrix and its inverse
    //FloatDPMatrix S(m+1,m+1); adat(S,AE,DE);
    //ARIADNE_LOG_PRINTLN("S="<<S);
    //S=FloatDPMatrix(m+1,m+1); simple_adat(S,AE,DE);
    //ARIADNE_LOG_PRINTLN("S="<<S);
    FloatDPApproximationMatrix S=feasibility_adat(H,A,DE);
    ARIADNE_LOG_PRINTLN("S="<<S);
    FloatDPApproximationMatrix Sinv=inverse(S);
    ARIADNE_LOG_PRINTLN("Sinv="<<Sinv);

    // FIXME: What if S is not invertible?

    // Construct the residuals
    FloatDPApproximationVector rx=esub(emul(x,z),mu*sigma);
    //RawFloatDPVector ryt=-prod(AE,x); ryt[m]+=1; // FIXME: Need hessian
    FloatDPApproximationVector ryt=-feasibility_mul(A,x); ryt[m]+=1; // FIXME: Need hessian
    FloatDPApproximationVector rz=gye+z;
    ARIADNE_LOG_PRINTLN("rx="<<rx<<" ryt="<<ryt<<" rz="<<rz);

    //RawFloatDPVector rr=prod(AE,ediv(RawFloatDPVector(rx-emul(x,rz)),z))-ryt;
    FloatDPApproximationVector rr=ryt + AE*ediv(FloatDPApproximationVector(rx-emul(x,rz)),z) - ryt;


    // Compute the differences
    FloatDPApproximationVector dyt=Sinv*rr;
    //RawFloatDPVector dz=-rz-prod(AET,dyt);
    FloatDPApproximationVector dz=-rz-feasibility_trmul(A,dyt);
    FloatDPApproximationVector dx=-ediv(FloatDPApproximationVector(rx+emul(x,dz)),z);
    ARIADNE_LOG_PRINTLN("dx="<<dx<<" dyt="<<dyt<<" dz="<<dz);

    FloatDPApproximationVector nx,ny,nyt,nz; FloatDPApproximation nt(dp);

    // Since we need to keep the point feasible, but the updates are linear
    // we need to validate feasibility directly rather than assuming the
    // linear update of y and z are good enough.
    Bool allpositive=false;
    FloatDPApproximation alpha=1/FloatDPApproximation(SCALE,dp);
    if(!egtr(emul(x,z) , GAMMA*mu/16)) {
        ARIADNE_LOG_PRINTLN("WARNING: x="<<x<<", z="<<z<< ", x.z="<<emul(x,z)<<"<"<<GAMMA*mu / 16);
        throw NearBoundaryOfFeasibleDomainException();
    }
    while(!allpositive) {
        alpha=alpha*SCALE;
        nx=x+alpha*dx;
        nyt=yt+alpha*dyt;
        ny=project(nyt,range(0,m));
        nt=nyt[m];
        //NonlinearInteriorPointOptimiser::compute_z(d,g,c,ny,nt,nz);
        allpositive = egtr(nx,0.0) && egtr(nz,0.0) && egtr(emul(nx,nz),GAMMA*mu);
    }
    ARIADNE_LOG_PRINTLN("alpha="<<alpha);
    ARIADNE_LOG_PRINTLN("nx="<<nx<<" nyt="<<nyt<<" nz="<<nz<<" nxz="<<emul(nx,nz));

    x=nx; y=project(nyt,range(0,m)); z=nz; t=nyt[m];
}

FloatDPApproximation NonlinearInteriorPointOptimiser::
compute_mu(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& C,
           const FloatDPApproximationVector& x, const FloatDPApproximationVector& lambda) const
{
    // Compute the relaxation parameter mu as the average of the product of the Lyapunov exponents and constraint satisfactions
    FloatDPApproximation mu=zero;
    FloatDPApproximationVector gx = g(x);

    for(Nat i=0; i!=C.size(); ++i) {
        if(C[i].lower()==C[i].upper()) { }
        else if(C[i].lower()==-infty) { mu += lambda[i] * (gx[i] - C[i].upper()); }
        else if(C[i].upper()==+infty) { mu += lambda[i] * (gx[i] - C[i].lower()); }
        else { // std::cerr<<"FIXME: Compute mu for singleton constraint\n";
            if ( decide(lambda[i] <=0.0) ) { mu += lambda[i] * (gx[i] - C[i].upper()); }
            else { mu += lambda[i] * (gx[i] - C[i].lower()); }
        }
    }
    mu /= C.size();
    return mu;
}


Void NonlinearInteriorPointOptimiser::
setup_feasibility(const ExactBoxType& d, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& c,
                  FloatDPApproximationVector& x, FloatDPApproximationVector& y) const
{
    const Nat l=2*(d.size()+c.size());
    y=midpoint(d);
    x=FloatDPApproximationVector(l,one/l);
    //compute_tz(d,g,c,y,t,z);
}




//------- PenaltyFunctionOptimiser ------------------------------------------//

PenaltyFunctionOptimiser* PenaltyFunctionOptimiser::
clone() const
{
    return new PenaltyFunctionOptimiser(*this);
}

ValidatedVectorType PenaltyFunctionOptimiser::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const
{
    ARIADNE_NOT_IMPLEMENTED;
}

ValidatedKleenean PenaltyFunctionOptimiser::
feasible(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("D="<<D<<" g="<<g<<" C="<<C<<" ");

    FloatDPApproximationVector x=midpoint(D);

    FloatDPApproximationVector w=midpoint(C);
    for(Nat i=0; i!=C.size(); ++i) {
        if(C[i].upper()==+infty) { w[i]=C[i].lower()+one; }
        else if(C[i].lower()==-infty) { w[i]=C[i].upper()-one; }
    }

    FloatDPApproximationVector y(C.size(),zero);

    ARIADNE_LOG_PRINTLN("x="<<x<<" w="<<w<<" y="<<y);

    for(Nat i=0; i!=10; ++i) {
        this->feasibility_step(D,g,C,x,y,w);
    }
    return this->check_feasibility(D,g,C,cast_exact(x),cast_exact(y));
}

Void PenaltyFunctionOptimiser::
feasibility_step(const ExactBoxType& X, const ApproximateVectorMultivariateFunction& g, const ExactBoxType& W,
                 FloatDPApproximationVector& x, FloatDPApproximationVector& w, FloatDPApproximation& mu) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ApproximateVectorMultivariateFunction h(0u,X.dimension());
    const Nat n=X.size();
    const Nat m=W.size();
    const Nat l=h.result_size();

    ARIADNE_LOG_PRINTLN("x="<<x);
    ARIADNE_LOG_PRINTLN("w="<<w);

    Vector<FloatDPApproximationDifferential> ddgx=g.evaluate(FloatDPApproximationDifferential::variables(2,x));
    Vector<FloatDPApproximationDifferential> ddhx=h.evaluate(FloatDPApproximationDifferential::variables(2,x));

    mu *= 0.5;
    ARIADNE_LOG_PRINTLN("mu="<<mu);

    // G is the constraint value vector
    FloatDPApproximationVector gx = ddgx.value();
    FloatDPApproximationVector hx = ddhx.value();
    ARIADNE_LOG_PRINTLN("g(x)="<<gx);
    ARIADNE_LOG_PRINTLN("h(x)="<<hx);

    // A is the transpose derivative matrix aij=dgi/dxj
    FloatDPApproximationMatrix A = transpose(ddgx.jacobian());
    ARIADNE_LOG_PRINTLN("A=Dg(x)="<<A);
    FloatDPApproximationMatrix B = transpose(ddhx.jacobian());
    // FIXME: Due to problems with zero-element differential, need to resize matrix if no h
    if(l==0) { B.resize(n,0); }
    ARIADNE_LOG_PRINTLN("B=Dh(x)="<<B);

    // H is the Hessian matrix H[i1,i2] = df/dx[i1]dx[i2] + Sum_[j] lambda[j]*dg[j]/dx[i1]dx[i2]
    FloatDPApproximationMatrix H(n,n,dp);
    for(Nat j=0; j!=m; ++j) { H += (gx[j]-w[j]) * ddgx[j].hessian(); }
    for(Nat k=0; k!=l; ++k) { H += (hx[k]) * ddhx[k].hessian(); }
    ARIADNE_LOG_PRINTLN("H="<<H);

    FloatDPApproximationDiagonalMatrix D(n,dp);
    FloatDPApproximationDiagonalMatrix E(m,dp);
    for(Nat i=0; i!=n; ++i) { D[i] = rec(sqr(x[i]-X[i].lower())) + rec(sqr(X[i].upper()-x[i])); }
    for(Nat j=0; j!=m; ++j) { E[j] = rec(sqr(w[j]-W[j].lower())) + rec(sqr(W[j].upper()-w[j])); }
    ARIADNE_LOG_PRINTLN("D="<<D);
    ARIADNE_LOG_PRINTLN("E="<<E);

    FloatDPApproximationMatrix S = H + B * transpose(B);
    S += D;
    ARIADNE_LOG_PRINTLN("S="<<S);

    FloatDPApproximationMatrix R=inverse(S);
    ARIADNE_LOG_PRINTLN("inverse(S)="<<R);

    // Compute residuals
    FloatDPApproximationVector rx = A*gx + B * hx ; // + 1/(x.upper()-x) + 1/x.lower()-x if no regularisation
    FloatDPApproximationVector rw = w-gx;

    ARIADNE_LOG_PRINTLN("rx="<<rx);
    ARIADNE_LOG_PRINTLN("rw="<<rw);

    FloatDPApproximationVector dx = R * (rx + A * rw);
    FloatDPApproximationVector dw = rw + transpose(A)*dx;
    ARIADNE_LOG_PRINTLN("dx="<<dx);
    ARIADNE_LOG_PRINTLN("dw="<<dw);


    FloatDPApproximationVector newx(n,dp);
    FloatDPApproximationVector neww(m,dp);

    static const FloatDPApproximation ALPHA_SCALE_FACTOR = 0.75_approx;

    FloatDPApproximation alpha = 1.0_approx;
    do {
        newx = x - alpha * dx;
        neww = w - alpha * dw;
        alpha *= ALPHA_SCALE_FACTOR;
    } while ( !contains(X,cast_exact(newx)) || !contains(W,cast_exact(neww)) );
    alpha /= ALPHA_SCALE_FACTOR;

    ARIADNE_LOG_PRINTLN("alpha="<<alpha);

    ARIADNE_LOG_PRINTLN("newx="<<newx);
    ARIADNE_LOG_PRINTLN("neww="<<neww<<"\n");

    x=newx;
    w=neww;

    return;
}


Void PenaltyFunctionOptimiser::
feasibility_step(const ExactBoxType& D, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& C,
                 FloatDPBoundsVector& x, FloatDPBoundsVector& w) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


// Use a penalty approach without multipliers on the constraint functions
// Solve g(x)=w, x in D, w in C; Lagrangian y.(g(x)-w)
Void PenaltyFunctionOptimiser::
feasibility_step(ExactBoxType const& D, ApproximateVectorMultivariateFunction const& g, ExactBoxType const& C,
                 FloatDPApproximationVector& x, FloatDPApproximationVector& y, FloatDPApproximationVector& w) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    auto m=y.size(); auto n=x.size();

    FloatDPApproximationVector cl=lower_bounds(C);
    FloatDPApproximationVector cu=upper_bounds(C);
    FloatDPApproximationVector dl=lower_bounds(D);
    FloatDPApproximationVector du=upper_bounds(D);

    ARIADNE_LOG_PRINTLN("D="<<D<<", g="<<g<<", C="<<C);
    ARIADNE_LOG_PRINTLN("dl="<<dl<<", du="<<du);
    ARIADNE_LOG_PRINTLN("cl="<<cl<<", cu="<<cu);
    ARIADNE_LOG_PRINTLN("w="<<w<<", x="<<x<<", y="<<y);

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

    Vector<Differential<ApproximateNumericType>> ddgx=g.evaluate(Differential<ApproximateNumericType>::variables(2,x));
    ARIADNE_LOG_PRINTLN("ddgx="<<ddgx);

    Vector<ApproximateNumericType> gx = ddgx.value();
    ARIADNE_LOG_PRINTLN("g(x)="<<gx);
    Matrix<ApproximateNumericType> A = ddgx.jacobian();
    ARIADNE_LOG_PRINTLN("Dg(x)="<<A);

    Vector<ApproximateNumericType> yA=transpose(A)*y;

    // H is the Hessian matrix H of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) $
    Matrix<ApproximateNumericType> YH(x.size(),x.size(),dp);
    for(Nat i=0; i!=y.size(); ++i) {
        YH+=y[i]*ddgx[i].hessian();
    }
    ARIADNE_LOG_PRINTLN("Y.D2g(x)="<<YH);

    Vector<ApproximateNumericType> recwu=cu-w; recwu=erec(recwu);
    Vector<ApproximateNumericType> recwl=w-cl; recwl=erec(recwl);
    Vector<ApproximateNumericType> recxu=du-x; recxu=erec(recxu);
    Vector<ApproximateNumericType> recxl=x-dl; recxl=erec(recxl);

    Vector<ApproximateNumericType> diagDw=esqr(recwu)+esqr(recwl);
    Matrix<ApproximateNumericType> Dw(m,m,dp); for(Nat i=0; i!=m; ++i) { Dw[i][i]=diagDw[i]; }
    DiagonalMatrix<ApproximateNumericType> Dx(esqr(recxu)+esqr(recxl));


    for(Nat i=0; i!=n; ++i) { YH[i][i]-=Dx[i]; }

    Matrix<ApproximateNumericType> AT=transpose(A);
    Matrix<ApproximateNumericType> Znm(n,m,dp);
    Matrix<ApproximateNumericType> Zmn(m,n,dp);
    Matrix<ApproximateNumericType> Zmm(m,m,dp);
    Matrix<ApproximateNumericType> Im=Matrix<ApproximateNumericType>::identity(m,dp);


    Matrix<ApproximateNumericType> S=cojoin(join(Dw,Zmn,Im),join(Znm,-YH,-AT),join(Im,-A,Zmm));
    Vector<ApproximateNumericType> r=join(recwu-recwl+y,recxu-recxl-yA,w-gx);

    for(Nat j=0; j!=m; ++j) {
        if(C[j].lower()==C[j].upper()) {
            S[j][j]=1;
            S[j][m+n+j]=0;
            S[m+n+j][j]=0;
            r[j]=0;
        }
    }

    Vector<ApproximateNumericType> swxy = -solve(S,r);

    Vector<ApproximateNumericType> sw(m,dp),sx(n,dp),sy(m,dp);
    sw = project(swxy,range(0,m));
    sx = project(swxy,range(m,m+n));
    sy = project(swxy,range(m+n,m+n+m));

    ApproximateNumericType al=one;
    ApproximateVectorType nw=w+al*sw;
    ApproximateVectorType nx=x+al*sx;
    ApproximateVectorType ny(m,dp);
    ARIADNE_LOG_PRINTLN("sx="<<sx);
    ARIADNE_LOG_PRINTLN("sw="<<sw);
    while(!contains(C,cast_exact(nw)) || !contains(D,cast_exact(nx))) {
        al*=0.75;
        nw=w+al*sw;
        nx=x+al*sx;
    }
    ARIADNE_LOG_PRINTLN("al="<<sw);
    ny=y+al*sy;

    w=nw; x=nx; y=ny;
}


ValidatedKleenean ApproximateOptimiser::
feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("D="<<D<<", h="<<h);
    FloatDPApproximationVector x=midpoint(D);
    FloatDPApproximationVector y(h.result_size(),zero);

    for(Nat i=0; i!=8; ++i) {
        this->feasibility_step(D,h,x,y);
    }

    if( decide(norm(h(x))<1e-10) ) { return true; }

    if(!possibly(contains(UpperIntervalType(dot(UpperIntervalVectorType(cast_exact(y)),apply(h,D))),zero))) { return false; }

    return indeterminate;
}

Void ApproximateOptimiser::
feasibility_step(const ExactBoxType& D, const ApproximateVectorMultivariateFunction& h,
                 FloatDPApproximationVector& x, FloatDPApproximationVector& y) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("x="<<x<<" y="<<y);
    static const double SCALE_FACTOR = 0.75;
    const Nat n=x.size();
    const Nat m=y.size();
    // Solve equations y Dh(x) - 1/(x-xl) + 1/(xu-x) = 0; h(x) = 0
    Vector<FloatDPApproximationDifferential> ddhx=h.evaluate(FloatDPApproximationDifferential::variables(2,x));
    FloatDPApproximationMatrix A = ddhx.jacobian();
    ARIADNE_LOG_PRINTLN_AT(1,"A="<<A<<" b="<<ddhx.value());

    FloatDPApproximationMatrix H(n,n,dp);
    for(Nat i=0; i!=m; ++i) { H += y[i] * ddhx[i].hessian(); }
    for(Nat j=0; j!=n; ++j) {
        H[j][j] += rec(sqr(x[j]-D[j].lower()));
        H[j][j] += rec(sqr(D[j].upper()-x[j]));
    }

    FloatDPApproximationVector rx = transpose(A) * y;
    for(Nat j=0; j!=n; ++j) {
        rx[j] -= rec(x[j]-D[j].lower());
        rx[j] += rec(D[j].upper()-x[j]);
    }
    FloatDPApproximationVector ry = ddhx.value();
    ARIADNE_LOG_PRINTLN("rx="<<rx<<" ry="<<ry);

    // S = A Hinv AT
    // H dx + AT dy = rx; A dx = ry;
    //  dx = Hinv ( rx - AT dy )
    //  dy = Sinv ( A Hinv rx - ry )
    FloatDPApproximationMatrix Hinv=inverse(H);
    ARIADNE_LOG_PRINTLN_AT(1,"H="<<H<<" Hinv="<<Hinv);
    FloatDPApproximationMatrix S=A*Hinv*transpose(A);
    FloatDPApproximationMatrix Sinv=inverse(S);
    ARIADNE_LOG_PRINTLN_AT(1,"S="<<S<<" Sinv="<<Sinv);
    FloatDPApproximationVector dy = Sinv * ( A*(Hinv*rx) - ry );
    FloatDPApproximationVector dx = Hinv * ( rx - transpose(A) * dy);
    ARIADNE_LOG_PRINTLN("dx="<<dx<<" dy="<<dy);

    FloatDPApproximation ax = one;
    FloatDPApproximationVector nx = x-ax*dx;
    while(!contains(D,cast_exact(nx))) {
        ax*=SCALE_FACTOR;
        nx = x - ax * dx;
    }
    FloatDPApproximationVector ny = y-ax*dy;
    ARIADNE_LOG_PRINTLN("nx="<<nx<<" ax="<<ax<<" ny="<<ny);
    ARIADNE_LOG_PRINTLN_AT(1,"h(x)="<<h(nx));

    x=nx; y=ny;
}


ValidatedKleenean PenaltyFunctionOptimiser::
check_feasibility(ExactBoxType D, ValidatedVectorMultivariateFunction g, ExactBoxType C,
                     ExactFloatDPVectorType x, ExactFloatDPVectorType y) const
{
    ARIADNE_PRECONDITION(D.size()==g.argument_size());
    ARIADNE_PRECONDITION(C.size()==g.result_size());
    ARIADNE_PRECONDITION(x.size()==D.size());
    ARIADNE_PRECONDITION(y.size()==C.size());

    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN_AT(1,"D="<<D<<" C="<<C);

    FloatDPBoundsVector gx=g(x);
    ARIADNE_LOG_PRINTLN_AT(1,"x="<<x<<" y="<<y<<" g(x)="<<gx);

    ValidatedKleenean result = true;

    List<Nat> equalities;
    for(Nat j=0; j!=C.size(); ++j) {
        if( definitely(gx[j].upper()<C[j].lower() || gx[j].lower()>C[j].upper()) ) {
            return false;
        }
        if(C[j].lower()==C[j].upper()) {
            equalities.append(j);
        } else {
            if(!definitely(contains(C[j],gx[j]))) { result = indeterminate; }
        }
    }

    if(definitely(result)) {
        if(equalities.empty()) { ARIADNE_LOG_PRINTLN("Feasible"); return true; }

        ValidatedVectorMultivariateFunction h(equalities.size(),g.domain());
        FloatDPBoundsVector c(equalities.size(),dp);
        for(Nat i=0; i!=equalities.size(); ++i) {
            h[i] = g[equalities[i]];
            c[i] = C[equalities[i]].lower();
        }
        ARIADNE_LOG_PRINTLN_AT(1,"g="<<g);
        ARIADNE_LOG_PRINTLN_AT(1,"h="<<h<<" c="<<c<<" h(x)-c="<<(h(x)-c));

        static const ExactDouble SEARCH_TOLERANCE=1e-8_pr;
        FloatDPBoundsVector W(h.result_size(),FloatDPBounds(-1,+1,dp)*SEARCH_TOLERANCE);
        FloatDPBoundsMatrix AT = transpose(midpoint(h.jacobian(x)));
        FloatDPBoundsVector B = x+AT*W;
        FloatDPBoundsMatrix IA = h.jacobian(B);
        ARIADNE_LOG_PRINTLN_AT(1,"AT="<<AT<<" IA="<<IA);
        ARIADNE_LOG_PRINTLN_AT(1,"B="<<B);

        // Perform an interval Newton step to try to attain feasibility
        FloatDPBoundsVector nW = inverse(IA*AT) * FloatDPBoundsVector(h(x)-cast_exact(c));
        ARIADNE_LOG_PRINTLN_AT(1,"W="<<W<<"\nnew_W="<<nW);
        if(definitely(subset(UpperBoxType(B),D)) && refines(nW,W)) { ARIADNE_LOG_PRINTLN("Feasible"); return true; }
        else { result=indeterminate; }
    }

    // Compute first-order approximation to g(D) centred at x.
    // For feasibilty, have yg(D) cap yC nonempty.
    // Estimate y g(X) = y g(x) + y Dg(X).(X-x)

    // Compute y.C
    UpperIntervalVectorType iy(y);
    UpperIntervalType yC = dot(iy,C);

    // Compute Taylor estimate of y g(X)
    ValidatedVectorMultivariateTaylorFunctionModelDP tg(D,g,default_sweeper());
    ValidatedScalarMultivariateTaylorFunctionModelDP tyg(D,default_sweeper());
    for(Nat j=0; j!=y.size(); ++j) { tyg += y[j]*tg[j]; }
    UpperIntervalType tygD = UpperIntervalType(tyg(cast_singleton(D)));

    UpperIntervalMatrixType dgD = jacobian_range(g,cast_vector(D));
    UpperIntervalVectorType ydgD = transpose(dgD) * y;

    FloatDPBounds ygx = dot(y,gx);

    UpperIntervalType ygD = UpperIntervalType(ygx);
    for(Nat i=0; i!=x.size(); ++i) {
        ygD += ydgD[i] * (D[i]-UpperIntervalType(x[i]));
    }

    ARIADNE_LOG_PRINTLN_AT(1,"yC="<<yC<<" tygD="<<tygD<<" ygD="<<ygD);

    if(definitely(is_empty(intersection(yC,ygD)))) { ARIADNE_LOG_PRINTLN("Infeasible"); return false; }
    else { return indeterminate; }
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
    ARIADNE_LOG(2,"feasibility_step");
    RawFloatDPVector xl=lower_bounds(D); RawFloatDPVector xu=upper_bounds(D);
    RawFloatDPVector zl=lower_bounds(C); RawFloatDPVector zu=upper_bounds(C);

    const Nat n=x.size();
    const Nat m=y.size();

    ARIADNE_LOG(4,"x="<<x<<" y="<<y<<" z="<<z);
    Vector<FloatDPDifferential> ddx = FloatDPDifferential::variables(2,x);
    Vector<FloatDPDifferential> ddgx = g.evaluate(ddx);

    FloatDPMatrix A=ddgx.jacobian();
    ARIADNE_LOG(6,"A="<<A);
    RawFloatDPVector v = join(join(x,z),y);

    RawFloatDPVector r(n+2*m,n+2*m);
    project(r,range(0,n)) = y * A;
    for(Nat i=0; i!=n; ++i) {
        r[i] += ( rec(x[i]-xl[i]) - rec(xu[i]-x[i]) );
    }
    for(Nat j=0; j!=m; ++j) {
        if(zl[j]==zu[j]) { assert(zu[j]==zl[j]); r[n+j] = z[j]-zl[j]; }
        else { r[n+j] = ( rec(z[j]-zl[j]) - rec(zu[j]-z[j]) - y[j] ); }
    }
    project(r,range(n+m,n+2*m)) = ddgx.value() - z;
    r[n+2*m]=0.0;
    ARIADNE_LOG(5,"r="<<r);

    FloatDPMatrix S(n+2*m+1,n+2*m+1);
    for(Nat j=0; j!=m; ++j) {
        FloatDPMatrix H=ddgx[j].hessian();
        for(Nat i1=0; i1!=n; ++i1) {
            for(Nat i2=0; i2!=n; ++i2) {
                S[i1][i2]+=y[j]*H[i1][i2];
            }
        }
    }
    for(Nat j=0; j!=m; ++j) {
        for(Nat i=0; i!=n; ++i) {
            S[i][j+m+n]=A[j][i];
            S[j+m+n][i]=A[j][i];
        }
    }
    for(Nat j=0; j!=m; ++j) {
        S[n+j][n+m+j] = -1.0;
        S[n+m+j][n+j] = -1.0;
        //if(zl[j]==zu[j]) { S[n+j][n+j] = -1.0; S[n+j][n+m+j] = 0.0; }
        if(zl[j]==zu[j]) { S[n+j][n+j] = +inf; }
        else { S[n+j][n+j] = - rec(sqr(z[j]-zl[j])) - rec(sqr(zu[j]-z[j])); }
    }
    for(Nat i=0; i!=n; ++i) {
        S[i][i]-= rec(sqr(xu[i]-x[i]));
        S[i][i]-= rec(sqr(x[i]-xl[i]));
    }

    for(Nat i=0; i!=n; ++i) {
        S[i][n+2*m] -= rec(xu[i]-x[i]);
        S[i][n+2*m] += rec(x[i]-xl[i]);
    }

    for(Nat j=0; j!=n; ++j) {
        //S[n+m+j][n+m+j] = -1.0/1024;
    }

    ARIADNE_LOG(5,"S="<<S);
    //ARIADNE_LOG(5,"S="<<std::fixed<<pretty(S));

    FloatDPMatrix Sinv = inverse(S);
    //ARIADNE_LOG(9,"Sinv="<<Sinv<<"\n);
    //ARIADNE_LOG(5,"Sinv="<<std::fixed<<pretty(Sinv));

    RawFloatDPVector dv = Sinv * r;
    ARIADNE_LOG(5,"dv="<<dv);

    FloatDP alpha = 1.0;
    RawFloatDPVector nv = v-dv;
    while(!contains(D,RawFloatDPVector(project(nv,range(0,n)))) || !contains(C,RawFloatDPVector(project(nv,range(n,n+m)))) ) {
        alpha *= 0.75;
        nv = v-alpha*dv;
    }

    ARIADNE_LOG(4,"nv="<<nv<<" a="<<alpha);

    x=project(nv,range(0,n));
    z=project(nv,range(n,n+m));
    y=project(nv,range(n+m,n+2*m));
    ARIADNE_LOG(4,"g(x)-z="<<g(x)-z);

}
*/

// Solve equations y Dh(x) - 1/(x-xl) + 1/(xu-x) = 0; h(x) = 0
ValidatedKleenean IntervalOptimiser::
feasible_zero(ExactBoxType D, ValidatedVectorMultivariateFunction h) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("D="<<D<<", h="<<h);

    const Nat n=D.size();

    FloatDPBoundsVector zl(n,dp), zu(n,dp);
    ExactFloatDPVectorType xl = Ariadne::lower_bounds(D);
    ExactFloatDPVectorType xu = Ariadne::upper_bounds(D);

    FloatDPBoundsVector x=cast_singleton(D);
    FloatDPBoundsVector y(h.result_size(),FloatDPBounds(-1,+1,dp));
    FloatDPBounds mu(0,1,dp);

    for(Nat i=0; i!=8; ++i) {
        this->feasibility_step(xl,xu,h,x,y,zl,zu,mu);
    }

    return indeterminate;
}

Void IntervalOptimiser::
feasibility_step(const ExactFloatDPVectorType& xl, const ExactFloatDPVectorType& xu, const ValidatedVectorMultivariateFunction& h,
                 FloatDPBoundsVector& x, FloatDPBoundsVector& y, FloatDPBoundsVector& zl, FloatDPBoundsVector zu, FloatDPBounds& mu) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    ARIADNE_LOG_PRINTLN("[x]="<<x<<" [lambda]="<<y<<", [zl]="<<zl<<", [zu]="<<zu<<" [mu]="<<mu);

    const Nat n=x.size();
    const Nat m=y.size();

    FloatDPBoundsVector mx=midpoint(x);
    FloatDPBoundsVector my=midpoint(y);
    FloatDPBoundsVector mzl=midpoint(zl);
    FloatDPBoundsVector mzu=midpoint(zu);
    FloatDPBounds mmu(midpoint(mu));
    ARIADNE_LOG_PRINTLN_AT(1,"x~"<<x<<" lambda~="<<y<<", mu~"<<mu);

    // Solve equations y Dh(x) - zl + zu = 0; h(x) = 0; (x-xl).zl - mu = 0;  (xu-x).zu - mu = 0; Sum_j y_j^2 - mu = 0
    Vector<FloatDPBoundsDifferential> ddhx=h.evaluate(FloatDPBoundsDifferential::variables(2,x));
    Vector<FloatDPBoundsDifferential> dhmx=h.evaluate(FloatDPBoundsDifferential::variables(1,mx));
    FloatDPBoundsMatrix A = ddhx.jacobian();
    FloatDPBoundsMatrix mA = dhmx.jacobian();
    ARIADNE_LOG_PRINTLN_AT(1,"A="<<A<<" b="<<ddhx.value());

    FloatDPBoundsVector rx = transpose(mA) * my;
    for(Nat j=0; j!=n; ++j) {
        rx[j] -= mmu*rec(mx[j]-xl[j]);
        rx[j] += mmu*rec(xu[j]-mx[j]);
    }
    FloatDPBoundsVector ry = dhmx.value();
    FloatDPBoundsVector rzl = esub(emul(FloatDPBoundsVector(mx-xl),mzl),mmu);
    FloatDPBoundsVector rzu = esub(emul(FloatDPBoundsVector(xu-mx),mzu),mmu);
    ARIADNE_LOG_PRINTLN("rx="<<rx<<" ry="<<ry<<" rzl="<<rzl<<" rzu="<<rzu);

    FloatDPBoundsMatrix H(n,n,dp);
    for(Nat i=0; i!=m; ++i) { H += y[i] * ddhx[i].hessian(); }
    for(Nat j=0; j!=n; ++j) {
        H[j][j] += mu*rec(sqr(x[j]-xl[j]));
        H[j][j] += mu*rec(sqr(xu[j]-x[j]));
    }

    // S = A Hinv AT
    // H dx + AT dy = rx; A dx = ry;
    //  dx = Hinv ( rx - AT dy )
    //  dy = Sinv ( A Hinv rx - ry )
    FloatDPBoundsMatrix Hinv=inverse(H);
    ARIADNE_LOG_PRINTLN_AT(1,"H="<<H<<" Hinv="<<Hinv);
    FloatDPBoundsMatrix S=A*Hinv*transpose(A);
    FloatDPBoundsMatrix Sinv=inverse(S);
    ARIADNE_LOG_PRINTLN_AT(1,"S="<<S<<" Sinv="<<Sinv);
    FloatDPBoundsVector dy = Sinv * ( A*(Hinv*rx) - ry );
    FloatDPBoundsVector dx = Hinv * ( rx - transpose(A) * dy);
    ARIADNE_LOG_PRINTLN("dx="<<dx<<" dy="<<dy);

    FloatDPBoundsVector nx = x-dx;
    FloatDPBoundsVector ny = y-dy;
    ARIADNE_LOG_PRINTLN("nx="<<nx<<" ny="<<ny);
    ARIADNE_LOG_PRINTLN_AT(1,"h(x)="<<h(nx));

    x = refinement(x,nx); y=refinement(y,ny);
    FloatDPBounds nmu = zero;
    for(Nat i=0; i!=m; ++i) { nmu += sqr(y[i]); }
    mu=refinement(mu,nmu);
}

/*

struct KuhnTuckerFunctionBody : VectorMultivariateFunctionMixin<KuhnTuckerFunctionBody,ExactIntervalType>
{
    ValidatedScalarMultivariateFunction f;
    Array<ValidatedScalarMultivariateFunction> g;
    Array<ValidatedScalarMultivariateFunction> df;
    Array<Array<ValidatedScalarMultivariateFunction> > dg;

    KuhnTuckerFunctionBody(ValidatedScalarMultivariateFunction _f, ValidatedVectorMultivariateFunction _g) {
        ARIADNE_ASSERT(_f.argument_size()==_g.argument_size());
        const Nat m=_g.argument_size();
        const Nat n=_g.result_size();
        g.resize(n); df.resize(m); dg.resize(n); for(Nat j=0; j!=n; ++j) { dg[j].resize(m); }
        f=_f;
        for(Nat j=0; j!=n; ++j) { g[j]=_g[j]; }
        for(Nat i=0; i!=m; ++i) { df[i]=f.derivative(i); }
        for(Nat j=0; j!=n; ++j) { for(Nat i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    Nat result_size() const { return g.size()*2+f.argument_size(); }
    Nat argument_size() const { return g.size()*2+f.argument_size(); }
    ValidatedScalarMultivariateFunction operator[](Nat) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream&) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> Void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const Nat m=f.argument_size();
        const Nat n=g.size();
        Vector<X> x(project(arg,range(0,n)));
        Vector<X> y(project(arg,range(n,n+m)));
        Vector<X> z(project(arg,range(n+m,n+m+n)));
        Vector<X> rx(m), rz(n), rs(n);
        for(Nat i=0; i!=m; ++i) { rx[i]=df[i].evaluate(y); for(Nat j=0; j!=n; ++j) { rx[i]=rx[i]-x[j]*dg[j][i].evaluate(y); } }
        for(Nat j=0; j!=n; ++j) { rz[j]=g[j].evaluate(y) + z[j]; }
        for(Nat j=0; j!=n; ++j) { rs[j]=x[j]*z[j]; }
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
        const Nat m=_g.argument_size();
        const Nat n=_g.result_size();
        g.resize(n); dg.resize(n); for(Nat j=0; j!=n; ++j) { dg[j].resize(m); }
        for(Nat j=0; j!=n; ++j) { g[j]=_g[j]; for(Nat i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    Nat result_size() const { return g.size()*2+g[0].argument_size()+1; }
    Nat argument_size() const { return g.size()*2+g[0].argument_size()+1; }
    ValidatedScalarMultivariateFunction operator[](Nat) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream&) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> Void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const Nat m=g[0].argument_size();
        const Nat n=g.size();
        Vector<X> x(project(arg,range(0,n)));
        Vector<X> y(project(arg,range(n,n+m)));
        Vector<X> z(project(arg,range(n+m,n+m+n)));
        X t(arg[n+m+n]);
        Vector<X> rx(m), rz(n), rs(n); X rt;
        for(Nat i=0; i!=m; ++i) { rx[i]=x[0]*dg[0][i].evaluate(y); for(Nat j=1; j!=n; ++j) { rx[i]=rx[i]+x[j]*dg[j][i].evaluate(y); } }
        for(Nat j=0; j!=n; ++j) { rz[j]=g[j].evaluate(y) + t + z[j]; }
        for(Nat j=0; j!=n; ++j) { rs[j]=x[j]*z[j]; }
        rt=1-x[0]; for(Nat j=1; j!=n; ++j) { rt=rt-x[j]; }
        project(res,range(0,n))=rz;
        project(res,range(n,n+m))=rx;
        project(res,range(n+m,n+m+n))=rs;
        res[n+m+n]=rt;
    }
};



struct ConstrainedFeasibilityKuhnTuckerFunctionBody : VectorMultivariateFunctionMixin<FeasibilityKuhnTuckerFunctionBody,ExactIntervalType>
{
    Nat m;
    Nat n;
    ExactIntervalVectorType d;
    Array<ValidatedScalarMultivariateFunction> g;
    ExactIntervalVectorType c;
    Array<Array<ValidatedScalarMultivariateFunction> > dg;

    ConstrainedFeasibilityKuhnTuckerFunctionBody(ExactBoxType D, ValidatedVectorMultivariateFunction _g, ExactBoxType C) {
        m=_g.argument_size();
        n=_g.result_size();
        d=D; c=C;
        g.resize(n); dg.resize(n); for(Nat j=0; j!=n; ++j) { dg[j].resize(m); }
        for(Nat j=0; j!=n; ++j) { g[j]=_g[j]; for(Nat i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    Nat result_size() const { return 5*m+4*n+1u; }
    Nat argument_size() const { return 5*m+4*n+1u; }
    ValidatedScalarMultivariateFunction operator[](Nat) const { ARIADNE_NOT_IMPLEMENTED; }
    OutputStream& _write(OutputStream& os) const { return os << "KuhnTuckerFunctionBody"; }

    template<class X> Void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const X zero=arg[0].zero_element();
        const Nat l=2*(m+n);
        assert(arg.size()==l+m+l+1);
        Vector<X> x(project(arg,range(0u,l)));
        Vector<X> y(project(arg,range(l,l+m)));
        Vector<X> z(project(arg,range(l+m,l+m+l)));
        X t(arg[l+m+l]);
        Vector<X> rx(m,zero), rz(l,zero), rs(l,zero); X rt(zero);
        Vector<X> gy(n);
        for(Nat j=0; j!=n; ++j) { gy[j]=g[j].evaluate(y); }
        Matrix<X> dgy(n,m);
        for(Nat i=0; i!=m; ++i) { for(Nat j=0; j!=n; ++j) { dgy[j][i]=dg[j][i].evaluate(y); } }

        for(Nat i=0; i!=m; ++i) {
            for(Nat j=0; j!=n; ++j) { rx[i]+=x[j]*(dgy[j][i]-c[j].upper()); rx[i]+=x[n+j]*(c[j].lower()-dgy[j][i]); }
            rx[i]+=x[2*n+i]*(y[i]-d[i].upper())-x[2*n+m+i]*(d[i].lower()-y[i]);
        }
        for(Nat j=0; j!=n; ++j) { rz[j]=gy[j] + t + z[j]; rz[n+j]=t+z[n+j]-gy[j]; }
        for(Nat i=0; i!=m; ++i) { rz[2*n+i]=y[i]+t+z[2*n+i]; rz[2*n+m+i]=y[i]+t+z[2*n+m+i]; }
        for(Nat k=0; k!=l; ++k) { rs[k]=x[k]*z[k]; }
        rt+=1.0; for(Nat j=0; j!=2*n; ++j) { rt=rt-x[j]; }
        project(res,range(0,l))=rz;
        project(res,range(l,l+m))=rx;
        project(res,range(l+m,l+m+l))=rs;
        res[l+m+l]=rt;
    }
};



ValidatedVectorType KrawczykOptimiser::
minimise(ValidatedScalarMultivariateFunction f, ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


ValidatedKleenean KrawczykOptimiser::
feasible(ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c) const
{
    ARIADNE_LOG(2,"KrawczykOptimiser::feasible(ExactBoxType d, ValidatedVectorMultivariateFunction g, ExactBoxType c)");
    ARIADNE_LOG(2,"  d="<<d<<", g="<<g<<", c="<<c);

    ARIADNE_ASSERT(g.argument_size()==d.size());
    ARIADNE_ASSERT(g.result_size()==c.size());

    ExactIntervalType t; ExactIntervalVectorType x,y,z;
    setup_feasibility(d,g,c,x,y,z,t);


    // FIXME: Allow more steps
    for(Nat i=0; i!=12; ++i) {
        ARIADNE_LOG(4,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", x="<<x<<", z="<<z);
        try {
            this->feasibility_step(d,g,c,x,y,z,t);
        }
        catch(const SingularMatrixException& e) {
            return indeterminate;
        }
        if(t.lower()>t.upper()) {
            ARIADNE_LOG(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
            return indeterminate;
        }
        if(t.lower()>0.0) {
            ARIADNE_LOG(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
            if(this->is_feasible_point(d,g,c,midpoint(y))) {
                return true;
            }
        }
        if(t.upper()<0.0) {
            ARIADNE_LOG(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
            return false;
        }
    }
    ARIADNE_LOG(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c);
    if(this->is_infeasibility_certificate(d,g,c,midpoint(x))) {
        return false;
    }
    return indeterminate;
}



Void KrawczykOptimiser::setup_feasibility(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                          ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z, ExactIntervalType& t) const
{
    const Nat m=g.argument_size();
    const Nat n=g.result_size();
    const Nat l=2*(m+n);
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

    const Nat m=g.argument_size();
    const Nat n=g.result_size();

    // Compute the image of y under the constraint function
    ExactIntervalVectorType gy=g(y);
    gy+=ExactIntervalVectorType(gy.size(),ExactIntervalType(-min_float,+min_float));
    ExactIntervalVectorType my=midpoint(y);
    ExactIntervalVectorType mgy=g(my);

    // Find the range of possible values of the optimal t
    // This range is too pessimistic
    t=ExactIntervalType(+inf,+inf);
    for(Nat j=0; j!=n; ++j) {
        t=min(t,c[j]-gy[j]);
        t=min(t,gy[j]-c[j]);
    }
    for(Nat i=0; i!=m; ++i) {
        t=min(t,d[i]-y[i]);
        t=min(t,y[i]-d[i]);
    }

    // Find the range of possible values of the optimal t
    FloatDP tmin=+inf;
    FloatDP tmax=+inf;
    for(Nat j=0; j!=n; ++j) {
        tmax=min(tmax,sub(up,c[j].upper(),gy[j].lower()));
        tmax=min(tmax,sub(up,gy[j].upper(),c[j].lower()));
        tmin=min(tmin,sub(down,c[j].upper(),mgy[j].upper()));
        tmin=min(tmin,sub(down,mgy[j].lower(),c[j].lower()));
    }
    for(Nat i=0; i!=m; ++i) {
        tmin=min(tmin,sub(up,d[i].upper(),y[i].lower()));
        tmax=min(tmax,sub(up,y[i].upper(),d[i].lower()));
        tmin=min(tmin,sub(down,d[i].upper(),my[i].upper()));
        tmax=min(tmax,sub(down,my[i].lower(),d[i].lower()));
    }
    tmin-=0.0625;
    t=ExactIntervalType(tmin,tmax);


    // Find the range of possible values of the optimal z
    // This range is too pessimistic
    for(Nat j=0; j!=n; ++j) {
        z[j]=max(c[j].upper()-gy[j]-t,0.0);
        z[n+j]=max(gy[j]-c[j].lower()-t,0.0);
    }
    for(Nat i=0; i!=m; ++i) {
        z[2*n+i]=max(d[i].upper()-y[i]-t,0.0);
        z[2*n+m+i]=max(y[i]-d[i].lower()-t,0.0);
    }

    // Find the range of possible values of the optimal z
    // This range is too pessimistic
    for(Nat j=0; j!=n; ++j) {
        z[j]=ExactIntervalType(0.0,c[j].upper()-mgy[j].lower()-tmin);
        z[n+j]=ExactIntervalType(0.0,mgy[j].upper()-c[j].lower()-tmin);
    }
    for(Nat i=0; i!=m; ++i) {
        z[2*n+i]=ExactIntervalType(0.0,d[i].upper()-my[i].lower()-tmin);
        z[2*n+m+i]=ExactIntervalType(0.0,my[i].upper()-d[i].lower()-tmin);
    }

    ARIADNE_LOG(9,"  d="<<d<<", c="<<c<<", y="<<y<<", g(y)="<<gy<<", t="<<t<<", z="<<z);

}


Void KrawczykOptimiser::
minimisation_step(const ValidatedScalarMultivariateFunction& f, const ValidatedVectorMultivariateFunction& g,
                  ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z) const
{
    const Nat m=f.argument_size();
    const Nat n=g.result_size();

    Differential<UpperIntervalType> ddf=f.evaluate(Differential<UpperIntervalType>::variables(2,y));
    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));

    ExactIntervalMatrixType H(m,m);
    set_hessian(H,ddf);
    for(Nat j=0; j!=n; ++j) { add_hessian(H,-x[j],ddg[j]); }

    ExactIntervalMatrixType A(m,n);
    set_jacobian_transpose(A,ddg);

    ARIADNE_LOG(9,"f="<<f<<"\ng="<<g<<"\nx="<<x<<" y="<<y<<" z="<<z);
    ARIADNE_LOG(9,"A="<<A<<"\nH="<<H);

    ARIADNE_NOT_IMPLEMENTED;

}



Void KrawczykOptimiser::feasibility_step(const ValidatedVectorMultivariateFunction& g,
                                         ExactIntervalVectorType& x, ExactIntervalVectorType& y, ExactIntervalVectorType& z, ExactIntervalType& t) const
{
    ARIADNE_NOT_IMPLEMENTED;
    const Nat m=y.size();
    const Nat n=x.size();

    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));

    // A is the transpose derivative matrix aij=dgj/dyi
    ExactIntervalMatrixType A(m,n);
    for(Nat i=0; i!=m; ++i) {
        for(Nat j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    ARIADNE_LOG(9,"A="<<A);

    // H is the Hessian matrix Hik = xj*dgj/dyidyk
    ExactIntervalMatrixType H(m,m);
    for(Nat j=0; j!=n; ++j) {
        add_hessian(H,x[j],ddg[j]);
    }
    ARIADNE_LOG(9," H="<<H);

    FloatDPMatrix mA=midpoint(A);
    ARIADNE_LOG(9," mA="<<mA);
    FloatDPMatrix mH=midpoint(H);
    ARIADNE_LOG(9," mH="<<mH);

    RawFloatDPVector mD(n);
    for(Nat j=0; j!=n; ++j) { mD[j]=midpoint(x[j])/midpoint(z[j]); }
    ARIADNE_LOG(9," mD="<<mD);

    FloatDPMatrix& mS=mH;
    adat(mS,mA,mD);
    ARIADNE_LOG(9,"mS="<<mS);
    FloatDPMatrix mSinv=inverse(mS);
    ARIADNE_LOG(9,"mSinv="<<mSinv);
}

// Feasibility step for dual (inequality constrained) problem without using slack variables
// FIXME: Do we need a slackness parameter mu? Probably not; hopefully the infinities are kept in check...
// This method has the advantage of not needing to update the primal variables
Void KrawczykOptimiser::feasibility_step(const ExactBoxType& d, const ValidatedVectorMultivariateFunction& g, const ExactBoxType& c,
                                         ExactIntervalVectorType& y, ExactIntervalType& t) const
{
    const Nat m=d.size();
    const Nat n=c.size();

    // Compute function values
    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));

    // gy is the vector of values of g(y)
    ExactIntervalVectorType gy(n);
    for(Nat j=0; j!=n; ++j) { gy[j]=ddg[j].value(); }

    // z is the vector of slack variables z[k]=cu[k]-gy[k]-t or z[k]=gy[k]-cl[k]-t
    ExactIntervalVectorType z(2*(m+n));
    for(Nat j=0; j!=n; ++j) { z[j]=d[j].upper()-gy[j]-t; z[n+j]=gy[j]-d[j].lower()-t; }
    for(Nat i=0; i!=m; ++i) { z[i]=c[2*n+i].upper()-y[i]-t; z[2*n+m+i]=y[i]-c[i].lower()-t; }

    ExactIntervalVectorType zr(2*(m+n));
    for(Nat k=0; k!=2*(m+n); ++k) { zr[k]=1.0/z[k]; }

    ExactIntervalVectorType D(2*(m+n));
    for(Nat k=0; k!=2*(m+n); ++k) { D[k]=zr[k]*zr[k]; }

    // A is the transpose derivative matrix aij=dgj/dyi
    ExactIntervalMatrixType A(m,n);
    for(Nat i=0; i!=m; ++i) { for(Nat j=0; j!=n; ++j) { A[i][j]=ddg[j][i]; } }

    // A is the sum of scaled Hessian matrices hi1i2=zj*ddgj/dyi1yi2
    ExactIntervalMatrixType H(m,m);

    ExactIntervalMatrixType SE(m+1,m+1);
    // SE[0:m][0:m] is the matrix H+/-A(D1+D2)AT+(D3+D4) where D=z.z
    for(Nat i1=0; i1!=m; ++i1) { for(Nat i2=0; i2!=m; ++i2) { SE[i1][i2]=H[i1][i2];
        for(Nat j=0; j!=n; ++j) { SE[i1][i2]+=A[i1][j]*(D[j]+D[n+j])*A[i2][j]; }
    } }
    for(Nat i=0; i!=m; ++i) { SE[i][i]+=(D[2*n+i]+D[2*n+m+i]); }
    // SE[m][0:m]=SE[0:m][m] is the vector A(D1-D2)e+(D3-D4)e
    for(Nat i=0; i!=m; ++i) { SE[i][m]=D[2*n+i]-D[2*n+m+i];
        for(Nat j=0; j!=n; ++j) { SE[i][m]+=A[i][j]*(D[j]-D[n+j]); }
        SE[m][i]=SE[i][m];
    }
    // SE[m][m] is the scalar eT(D1+D2)e+eT(D3+D4)e
    for(Nat k=0; k!=2*(m+n); ++k) { SE[m][m]+=D[k]; }

    // Vector of residuals
    ExactIntervalVectorType re(m+1);
    for(Nat i=0; i!=m; ++i) { re[i]+=(zr[2*n+i]-zr[2*n+m+i]);
        for(Nat j=0; j!=n; ++j) { re[i]+=A[i][j]*(zr[j]-zr[n+j]); }
    }
    for(Nat j=0; j!=n; ++j) { re[m]+=(zr[j]+zr[n+n]); }
    for(Nat i=0; i!=m; ++i) { re[m]+=(zr[2*n+i]+zr[2*n+m+i]); }

    // Compute inverse Jacobian matrix
    ExactIntervalMatrixType JE;
    try {
        JE=inverse(midpoint(SE));
    }
    catch(const SingularMatrixException& e) {
        ARIADNE_WARN("Matrix S="<<midpoint(SE)<<" is not invertible");
        ARIADNE_LOG(1,"WARNING: Matrix S="<<midpoint(SE)<<" is not invertible");
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
    const Nat m=d.size();
    const Nat n=c.size();
    const Nat o=2*(m+n);

    ARIADNE_ASSERT_MSG(g.argument_size()==m,"d="<<d<<" g="<<g);
    ARIADNE_ASSERT_MSG(g.result_size()==n,"d="<<d<<" g="<<g<<" c="<<c);
    ARIADNE_ASSERT(x.size()==o);
    ARIADNE_ASSERT(y.size()==m);
    ARIADNE_ASSERT(z.size()==o);

    ExactIntervalVectorType yt=join(y,t);
    ARIADNE_LOG(9,"m="<<m<<" n="<<n);
    ARIADNE_LOG(9,"x="<<x<<" yt="<<yt<<" z="<<z);

    Vector< Differential<UpperIntervalType> > ddg=g.evaluate(Differential<UpperIntervalType>::variables(2,y));
    ARIADNE_LOG(9,"  ddg="<<ddg);

    // gy is the vector of values of g(y)
    UpperIntervalVectorType gy(n); for(Nat j=0; j!=n; ++j) { gy[j]=ddg[j].value(); }
    ARIADNE_LOG(9,"  g(y)="<<gy<<" ");

    // A is the transpose derivative matrix aij=dgj/dyi, extended with a column of ones
    UpperIntervalMatrixType A(m,n);
    for(Nat i=0; i!=m; ++i) {
        for(Nat j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    ARIADNE_LOG(9," A="<<A<<" ");

    // H is the Hessian matrix Hik = (xcuj-xclj)*dgj/dyidyk
    UpperIntervalMatrixType H(m,m);
    for(Nat j=0; j!=n; ++j) {
        add_hessian(H,x[j]-x[n+j],ddg[j]);
    }
    ARIADNE_LOG(9," H="<<H);

    // Construct the extended valuation GY=(gy-cu+te,cl-gy+te,y-bu+te,bl-y+te)
    UpperIntervalVectorType gye(o);
    for(Nat j=0; j!=n; ++j) { gye[j]=gy[j]-c[j].upper()+t; gye[n+j]=c[j].lower()-gy[j]+t; }
    for(Nat i=0; i!=m; ++i) { gye[2*n+i]=y[i]-d[i].upper()+t; gye[2*n+m+i]=d[i].lower()-y[i]+t; }
    ARIADNE_LOG(9,"  GE="<<gye);

    // Construct the extended matrix AE=(A -A I -I \\ e e 0 0)
    UpperIntervalMatrixType AE(m+1,o);
    for(Nat i=0; i!=m; ++i) { for(Nat j=0; j!=n; ++j) { AE[i][j]=A[i][j]; AE[i][n+j]=-A[i][j]; } }
    for(Nat i=0; i!=m; ++i) { AE[i][2*n+i]=1; AE[i][2*n+m+i]=-1; }
    for(Nat k=0; k!=o; ++k) { AE[m][k]=1; }
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
    //ARIADNE_LOG(9,"S="<<S);
    //S=FloatDPMatrix(m+1,m+1); simple_adat(S,AE,DE);
    //ARIADNE_LOG(9,"S="<<S);
    FloatDPMatrix mS=feasibility_adat(mH,mA,mDE);
    ARIADNE_LOG(9,"mS="<<mS);
    FloatDPMatrix mSinv=inverse(mS);
    ARIADNE_LOG(9,"mSinv="<<mSinv);

    // FIXME: What if S is not invertible?

    // Construct the residuals
    UpperIntervalVectorType rx=emul(mx,mz);
    //RawFloatDPVector ryt=-prod(AE,x); ryt[m]+=1; // FIXME: Need hessian
    UpperIntervalVectorType ryt=-feasibility_mul(mA,mx); ryt[m]+=1; // FIXME: Need hessian
    UpperIntervalVectorType rz=midpoint(gye)+mz;
    ARIADNE_LOG(9,"rx="<<rx<<" ryt="<<ryt<<" rz="<<rz);

    // Construct the errors on the residuals ([M]-M)([x]-x)
    UpperIntervalVectorType ex=x-mx;
    UpperIntervalVectorType eyt=yt-myt;
    UpperIntervalVectorType ez=z-mz;
    UpperIntervalMatrixType eA=A-mA;
    UpperIntervalMatrixType eH=H-mH;

    UpperIntervalVectorType erx=2.0*emul(ex,ez);
    UpperIntervalVectorType eryt=UpperIntervalMatrixType(AE-mAE)*ex;
    UpperIntervalVectorType erz=UpperIntervalMatrixType(AET-mAET)*eyt;
    ARIADNE_LOG(9,"erx="<<erx<<" eryt="<<eryt<<" erz="<<erz);

    rx+=2.0*emul(ex,ez);
    ryt+=UpperIntervalMatrixType(AE-mAE)*ex;
    rz+=UpperIntervalMatrixType(AET-mAET)*eyt;
    ARIADNE_LOG(9,"rx="<<rx<<" ryt="<<ryt<<" rz="<<rz);

    //RawFloatDPVector rr=prod(AE,ediv(RawFloatDPVector(rx-emul(x,rz)),z))-ryt;

    // Compute the error differences
    UpperIntervalVectorType erxdz=ediv(erx,mz);
    UpperIntervalVectorType edyt=(mSinv*mAE)*erxdz + mSinv*eyt - (mSinv*(mAE*DiagonalMatrix<FloatDP>(mDE))) * ez;
    UpperIntervalVectorType edz=-erz-feasibility_trmul(mA,edyt);
    UpperIntervalVectorType edx=-ediv(UpperIntervalVectorType(erx+emul(mx,edz)),mz);
    ARIADNE_LOG(9,"edx="<<edx<<" edyt="<<edyt<<" edz="<<edz);

    // Compute the error differences
    UpperIntervalVectorType eerr=prod(mAE,ediv(esub(erx,emul(mx,erz)),mz))-eryt;
    ARIADNE_LOG(9,"  eerr="<<eerr);
    UpperIntervalVectorType eedyt=prod(mSinv,eerr);
    UpperIntervalVectorType eedz=-erz-feasibility_trmul(mA,eedyt);
    UpperIntervalVectorType eedx=-ediv(UpperIntervalVectorType(erx+emul(mx,eedz)),mz);
    ARIADNE_LOG(9,"eedx="<<eedx<<" eedyt="<<eedyt<<" eedz="<<eedz);


    // Compute the differences
    UpperIntervalVectorType rr=prod(mAE,ediv(esub(rx,emul(mx,rz)),mz))-ryt;
    UpperIntervalVectorType dyt=prod(mSinv,rr);
    UpperIntervalVectorType dz=-rz-feasibility_trmul(mA,dyt);
    UpperIntervalVectorType dx=-ediv(UpperIntervalVectorType(rx+emul(mx,dz)),mz);
    ARIADNE_LOG(9,"dx="<<dx<<" dyt="<<dyt<<" dz="<<dz<<"\n");

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
