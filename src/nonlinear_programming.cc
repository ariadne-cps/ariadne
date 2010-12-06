/***************************************************************************
 *            nonlinear_programming.cc
 *
 *  Copyright 2010  Pieter Collins
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

#include <limits>

#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "boost/numeric/ublas/storage.hpp"
#include "boost/numeric/ublas/vector.hpp"

#include "macros.h"
#include "logging.h"
#include "tuple.h"
#include "tribool.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "differential.h"
#include "function.h"
#include "function_mixin.h"
#include "taylor_function.h"

#include "nonlinear_programming.h"
#include "solver.h"

namespace Ariadne {

static const double error =  1e-2;

typedef Vector<Float> FloatVector;
typedef Matrix<Float> FloatMatrix;
typedef boost::numeric::ublas::vector_range<FloatVector> FloatVectorRange;
typedef DiagonalMatrix<Float> FloatDiagonalMatrix;

typedef Vector<Interval> IntervalVector;
typedef Matrix<Interval> IntervalMatrix;
typedef boost::numeric::ublas::vector_range<IntervalVector> IntervalVectorRange;
typedef DiagonalMatrix<Interval> IntervalDiagonalMatrix;

inline Vector<Float> lower_bounds(const Vector<Interval>& iv) {
    Vector<Float> lv(iv.size()); for(uint i=0; i!=lv.size(); ++i) { lv[i]=iv[i].lower(); } return lv;
}

template<class X, class XX> inline
bool egtr(const Vector<X>& x, const XX& s) {
    for(uint i=0; i!=x.size(); ++i) { if(x[i]<=s) { return false; } } return true;
}

template<class X> inline
Vector<X> esub(const Vector<X>& x, const X& s) {
    Vector<X> r(x.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=x[i]-s; } return r;
}

template<class X> inline
Vector<X> esub(const Vector<X>& x, const Vector<X>& y) {
    Vector<X> r(x.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=x[i]-y[i]; } return r;
}

template<class X> inline
Vector<X> emul(const Vector<X>& x, const Vector<X>& z) {
    Vector<X> r(x.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=x[i]*z[i]; } return r;
}

inline
Vector<Interval> emul(const Vector<Interval>& x, const Vector<Float>& z) {
    Vector<Interval> r(x.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=x[i]*z[i]; } return r;
}

inline
Vector<Interval> emul(const Vector<Float>& x, const Vector<Interval>& z) {
    Vector<Interval> r(x.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=x[i]*z[i]; } return r;
}

template<class X, class XX> inline
Vector<X> ediv(const Vector<X>& x, const Vector<XX>& z) {
    Vector<X> r(x.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=x[i]/z[i]; } return r;
}

inline
Interval eivl(const FloatVector& x) {
    Interval r(x[0]); for(uint i=0; i!=x.size(); ++i) { r=hull(r,x[i]); } return r;
}

// Compute S+=ADA^T, where D is diagonal and S is symmetric.
template<class X>
void adat(Matrix<X>& S, const Matrix<X>& A, const Vector<X>& D)
{
    const uint m=A.row_size();
    const uint n=A.column_size();
    for(uint i1=0; i1!=m; ++i1) {
        for(uint j=0; j!=n; ++j) {
            X ADij=A[i1][j]*D[j];
            for(uint i2=i1; i2!=m; ++i2) {
                S[i1][i2]+=ADij*A[i2][j];
            }
        }
    }
    for(uint i1=1; i1!=m; ++i1) {
        for(uint i2=0; i2!=i1; ++i2) {
            S[i1][i2]=S[i2][i1];
        }
    }
}

// Compute S+=A^TDA, where D is diagonal and S is symmetric.
template<class X>
void atda(Matrix<X>& S, const Matrix<X>& A, const Vector<X>& D)
{
    assert(S.row_size()==S.column_size());
    assert(S.column_size()==A.column_size());
    assert(D.size()==A.row_size());

    const uint m=A.column_size();
    const uint n=A.row_size();
    for(uint i1=0; i1!=m; ++i1) {
        for(uint j=0; j!=n; ++j) {
            X ATDij=A[j][i1]*D[j];
            for(uint i2=i1; i2!=m; ++i2) {
                S[i1][i2]+=ATDij*A[j][i2];
            }
        }
    }
    for(uint i1=1; i1!=m; ++i1) {
        for(uint i2=0; i2!=i1; ++i2) {
            S[i1][i2]=S[i2][i1];
        }
    }
}

// Compute L=DA^T, where D is diagonal and S is symmetric.
template<class X>
Matrix<X> dat(const Matrix<X>& A, const Vector<X>& D)
{
    Matrix<X> S(A.column_size(),A.row_size());
    for(uint i=0; i!=S.row_size(); ++i) {
        for(uint j=0; j!=S.column_size(); ++j) {
            S[i][j] = D[i]*A[j][i];
        }
    }
    return S;
}

// Compute S=ADA^T, where D is diagonal.
template<class X>
Matrix<X> adat(const Matrix<X>& A, const Vector<X>& D)
{
    const uint m=A.row_size();
    Matrix<X> S=Matrix<X>::zero(m,m);
    adat(S,A,D);
    return S;
}

// Compute S+=AA^T
template<class X>
Matrix<X> amulat(const Matrix<X>& A)
{
    const uint m=A.row_size();
    const uint n=A.column_size();
    Matrix<X> S(m,m);
    for(uint i1=0; i1!=m; ++i1) {
        for(uint j=0; j!=n; ++j) {
            for(uint i2=i1; i2!=m; ++i2) {
                S[i1][i2]+=A[i1][j]*A[i2][j];
            }
        }
    }
    for(uint i1=1; i1!=m; ++i1) {
        for(uint i2=0; i2!=i1; ++i2) {
            S[i1][i2]=S[i2][i1];
        }
    }
    return S;
}

template<class X> inline bool all_greater(const Vector<X>& x, const X& e) {
    for(uint i=0; i!=x.size(); ++i) { if(x[i]<=e) { return false; } } return true;
}

template<class X> inline Vector<X> operator*(const Matrix<X>& A, const Vector<X>& b) {
    return prod(A,b);
}

template<class X> inline Matrix<X> operator*(const Matrix<X>& A, const Matrix<X>& B) {
    return prod(A,B);
}

template<class X> inline Matrix<X> operator*(const Matrix<X>& A, const DiagonalMatrix<X>& B) {
    Matrix<X> R(A.row_size(),A.column_size());
    for(uint i=0; i!=A.row_size(); ++i) { for(uint j=0; j!=A.column_size(); ++j) { R[i][j]=A[i][j]*B.diagonal()[j]; } }
    return R;
}

template<class X> inline Vector<X> prod(const Vector<X>& v, const DiagonalMatrix<X>& D) {
    Vector<X> r(v.size()); for(uint i=0; i!=r.size(); ++i) { r[i] = v[i] * D[i]; } return r;
}

inline Vector<Float> operator*(const Vector<Float>& v1, const DiagonalMatrix<Float>& D2) { return prod(v1,D2); }

inline Matrix<Float> operator*(const Matrix<Float>& A1, const Matrix<Float>& A2) { return prod(A1,A2); }
inline Vector<Float> operator*(const Matrix<Float>& A1, const Vector<Float>& v2) { return prod(A1,v2); }
inline Vector<Float> operator*(const Vector<Float>& v1, const Matrix<Float>& A2) { return prod(v1,A2); }
inline Matrix<Interval> operator*(const Matrix<Float>& A1, const Matrix<Interval>& A2) { return prod(A1,A2); }
inline Vector<Interval> operator*(const Matrix<Float>& A1, const Vector<Interval>& v2) { return prod(A1,v2); }
inline Vector<Interval> operator*(const Vector<Float>& v1, const Matrix<Interval>& A2) { return prod(v1,A2); }
inline Matrix<Interval> operator*(const Matrix<Interval>& A1, const Matrix<Float>& A2) { return prod(A1,A2); }
inline Vector<Interval> operator*(const Matrix<Interval>& A1, const Vector<Float>& v2) { return prod(A1,v2); }
inline Vector<Interval> operator*(const Vector<Interval>& v1, const Matrix<Float>& A2) { return prod(v1,A2); }
inline Matrix<Interval> operator*(const Matrix<Interval>& A1, const Matrix<Interval>& A2) { return prod(A1,A2); }
inline Vector<Interval> operator*(const Matrix<Interval>& A1, const Vector<Interval>& v2) { return prod(A1,v2); }
inline Vector<Interval> operator*(const Vector<Interval>& v1, const Matrix<Interval>& A2) { return prod(v1,A2); }

template<class X> inline Matrix<X>& operator+=(Matrix<X>& A, const DiagonalMatrix<X>& D) {
    for(uint i=0; i!=D.size(); ++i) { A[i][i]+=D[i]; } return A;
}


template<class X> inline Vector<X> operator+(const Vector<X>& x, const Vector<X>& y) {
    Vector<X> r(x.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=x[i]+y[i]; } return r;
}

template<class X> inline Vector<X> operator-(const Vector<X>& x, const Vector<X>& y) {
    Vector<X> r(x.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=x[i]-y[i]; } return r;
}

template<class X> Vector<X> join(const Vector<X>& v1, const Vector<X>& v2, const Vector<X>& v3) {
    Vector<X> r(v1.size()+v2.size()+v3.size());
    for(uint i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(uint i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    for(uint i=0; i!=v3.size(); ++i) { r[v1.size()+v2.size()+i]=v3[i]; }
    return r;
}

template<class X> Vector<X> join(const Vector<X>& v1, const Vector<X>& v2, const Vector<X>& v3, const X& s4) {
    Vector<X> r(v1.size()+v2.size()+v3.size()+1u);
    for(uint i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(uint i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    for(uint i=0; i!=v3.size(); ++i) { r[v1.size()+v2.size()+i]=v3[i]; }
    r[r.size()-1]=s4;
    return r;
}


template<class X> Vector< Differential<X> > second_derivative(const IntervalVectorFunction& f, const Vector<X>& x) {
    Vector< Differential<X> > d=Differential<X>::variables(f.result_size(),f.argument_size(),2);
    return f.evaluate(d);
}

template<class Vec, class Diff> void set_gradient(Vec& g, const Diff& D) {
    typedef typename Diff::ValueType X;
    uint i=0;
    typename Diff::const_iterator iter=D.begin();
    if(iter!=D.end() && iter->key().degree()==0) { ++iter; }
    while(iter!=D.end() && iter->key().degree()<=2) {
        while(iter->key()[i]==0) { ++i; }
        g[i]=iter->data();
        ++iter;
    }
}

template<class Mx, class Diff> void set_jacobian_transpose(Mx& A, const Vector<Diff>& D) {
    for(uint j=0; j!=A.column_size(); ++j) {
        for(uint i=0; i!=A.row_size(); ++i) {
            A[i][j]=D[j][i];
        }
    }
}

template<class Mx, class Diff> void set_hessian(Mx& H, const Diff& D) {
    typedef typename Diff::ValueType X;
    uint i=0; uint j=1;
    typename Diff::const_iterator iter=D.begin();
    while(iter!=D.end() && iter->key().degree()<=1) { ++iter; }
    while(iter!=D.end() && iter->key().degree()<=2) {
        const MultiIndex& a=iter->key();
        const X& c=iter->data();
        while(a[i]==0) { ++i; j=i+1; }
        if(a[i]==2) { H[i][i]=c; }
        else { while(a[j]==0) { ++j; } H[i][j]=c; H[j][i]=c; }
        ++iter;
    }
}

template<class Mx, class S, class Diff> void add_hessian(Mx& H, const S& s, const Diff& D) {
    typedef typename Diff::ValueType X;
    typename Diff::const_iterator iter=D.begin();
    while(iter!=D.end() && iter->key().degree()<=1) { ++iter; }
    while(iter!=D.end() && iter->key().degree()==2) {
        const MultiIndex& a=iter->key();
        const X& c=iter->data();
        uint i=0;
        while(a[i]==0) { ++i; }
        if(a[i]==2) { H[i][i]+=s*c; }
        else { uint j=i+1; while(a[j]==0) { ++j; } H[i][j]+=s*c; H[j][i]+=s*c; }
        ++iter;
    }
}

// Compute the product (A -A I -I ; 1 1 1 1) v
template<class X, class XX> Vector<X> feasibility_mul(const Matrix<XX>& A, const Vector<X>& v)
{
    const uint m=A.row_size();
    const uint n=A.column_size();
    ARIADNE_ASSERT(v.size()==2*(m+n));
    Vector<X> r(m+1u);
    for(uint i=0; i!=m; ++i) {
        r[i]=v[2*n+i]-v[2*n+m+i];
        for(uint j=0; j!=n; ++j) {
            r[i]+=A[i][j]*(v[j]-v[n+j]);
        }
    }
    for(uint k=0; k!=2*(m+n); ++k) {
        r[m]+=v[k];
    }
    return r;
}

// Compute the product (AT 1 \\ -AT 1 \\ I 1 \\ -I 1) I -I ; 1 1 1 1) v
template<class X, class XX> Vector<X> feasibility_trmul(const Matrix<XX>& A, const Vector<X>& w)
{
    const uint m=A.row_size();
    const uint n=A.column_size();
    ARIADNE_ASSERT(w.size()==m+1);
    Vector<X> r(2*(m+n));
    for(uint j=0; j!=n; ++j) {
        r[j]=0;
        for(uint i=0; i!=m; ++i) {
            r[j]+=A[i][j]*w[i];
        }
        r[n+j]=-r[j];
        r[j]+=w[m];
        r[n+j]+=w[m];
    }
    for(uint i=0; i!=m; ++i) {
        r[2*n+i]=w[i]+w[m];
        r[2*n+m+i]=-w[i]+w[m];
    }
    return r;
}


// Compute the product \f$\hat{A}^T \hat{D} \hat{A} + \hat{H}\f$ where \f$\hat{A}=\left(\begin{matrix}A&-A&I&-I\\1&1&1&1\end{matrix}\right)\f$ and \f$\hat{D}=D\f$ is diagonal.
template<class X> Matrix<X> feasibility_adat(const Matrix<X>& H, const Matrix<X>& A, const Vector<X>& D)
{
    const uint m=A.row_size();
    const uint n=A.column_size();
    ARIADNE_ASSERT(H.row_size()==m);
    ARIADNE_ASSERT(H.column_size()==m);
    ARIADNE_ASSERT(D.size()==2*(m+n));
    Matrix<X> S(m+1,m+1);

    for(uint i=0; i!=m; ++i) { for(uint j=0; j!=m; ++j) { S[i][j] = H[i][j]; } }
    for(uint i=0; i!=m; ++i) { S[i][m]=0; S[m][i]=0; } S[m][m]=0;

    for(uint i1=0; i1!=m; ++i1) {
        for(uint j=0; j!=n; ++j) {
            X ADij=A[i1][j]*(D[j]+D[n+j]);
            for(uint i2=0; i2!=m; ++i2) {
                S[i1][i2]+=ADij*A[i2][j];
            }
        }
    }
    for(uint i=0; i!=m; ++i) {
        S[i][i]+=(D[2*n+i]+D[2*n+m+i]);
    }
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=n; ++j) {
            S[i][m]+=A[i][j]*(D[j]-D[n+j]);
        }
        S[i][m]+=(D[2*n+i]-D[2*n+m+i]);
        S[m][i]=S[i][m];
    }
    for(uint k=0; k!=2*(m+n); ++k) {
        S[m][m]+=D[k];
    }

    return S;
}






template<class R>
class ConstrainedFeasibilityMatrix {
    ConstrainedFeasibilityMatrix(const Vector<R>& x, const Vector<R>& z, const Matrix<R>& a, const Matrix<R>& h)
        : X(x), Z(z), D(ediv(x,z)), A(a), H(h) { }

    template<class RR> tuple< Vector<RR>,Vector<RR>,Vector<RR> >
    mul(const Vector<RR>& x, const Vector<RR>& yt, const Vector<RR>& z) const {
        Vector<RR> nx=Z*x+X*x;
        Vector<RR> nyt=H*yt-A*x;
        Vector<RR> nz=H*yt-A*x;
        return make_tuple(nx,nyt,nz);
    }

    template<class RR> tuple< Vector<RR>,Vector<RR>,Vector<RR> >
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



Bool OptimiserBase::
is_feasible_point(IntervalVector D, IntervalVectorFunction g, IntervalVector C, FloatVector x) const
{
    if(!contains(D,x)) { return false; }
    IntervalVector gx=g(IntervalVector(x));
    return subset(gx,C);
}


Tribool OptimiserBase::
contains_feasible_point(IntervalVector D, IntervalVectorFunction g, IntervalVector C, IntervalVector X) const
{
    ARIADNE_LOG(4,"OptimiserBase::contains_feasible_point(D,g,C,X):\n");
    ARIADNE_LOG(5,"  D="<<D<<", g="<<g<<", C="<<C<<", X="<<X<<"\n");

    // Now test if the (reduced) box X satisfies other constraints
    if(disjoint(X,D)) { return false; }
    if(!subset(X,D)) { return indeterminate; }

    // Test inequality constraints
    Tribool result = true;
    IntervalVector gx=g(X);
    ARIADNE_LOG(7,"g(X)="<<gx<<"\n");
    for(uint i=0; i!=C.size(); ++i) {
        if(disjoint(gx[i],C[i])) {
            return false;
        }
        if(!C[i].singleton()) {
            if(!subset(gx[i],C[i])) { result = indeterminate; }
        }
    }

    // Break if some inequality constraints indefinite
    if(!definitely(result)) { return result; }

    // Extract the equality constraints
    List<uint> equality_constraints;
    equality_constraints.reserve(C.size());
    for(uint i=0; i!=C.size(); ++i) {
        if(C[i].singleton()) { equality_constraints.append(i); }
    }

    // Construct the function g_e(x) = g_{e_i}(x)
    IntervalVectorFunction ge(equality_constraints.size(),g[0]);
    IntervalVector ce(equality_constraints.size());
    for(uint i=0; i!=ge.result_size(); ++i) {
        ge[i]=g[equality_constraints[i]];
        ce[i]=C[equality_constraints[i]];
    }

    ARIADNE_LOG(7,"ge="<<ge<<", ce="<<ce<<"\n");

    IntervalMatrix ivlA=ge.jacobian(X);
    ARIADNE_LOG(7,"ivlA="<<ivlA<<"\n");
    FloatVector fltD(X.size());
    for(uint i=0; i!=X.size(); ++i) { fltD[i]=rec(sqr(rad(X[i]))); }
    FloatMatrix fltA=midpoint(ivlA);
    ARIADNE_LOG(7,"A="<<fltA<<"\n");
    ARIADNE_LOG(7,"D="<<fltD<<"\n");
    FloatMatrix fltL = dat(fltA,fltD);
    ARIADNE_LOG(7,"L="<<fltL<<"\n");

    IntervalMatrix ivlS = ivlA * fltL;
    ARIADNE_LOG(7,"ivlS="<<ivlS<<"\n");

    IntervalMatrix ivlR = inverse(ivlS);
    try {
        ivlR=inverse(ivlS);
    }
    catch (SingularMatrixException e) {
        return indeterminate;
    }

    ARIADNE_LOG(7,"ivlR="<<ivlR<<"\n");


    // Projected interval Newton step. For h:R^n->R^m; Dh mxn, take L nxm.
    // Interval Newton update X' = x - L * (Dh(X)*L)^{-1} * h(x)
    // Choose L = rad(X)^2 Dh(x)^T where rad(X) is the diagonal matrix of radii of X
    IntervalVector x=midpoint(X);
    IntervalVector new_X = x - fltL * (ivlR * (ge(x)-ce) );
    ARIADNE_LOG(5,"old_X="<<X<<"\n");
    ARIADNE_LOG(5,"new_X="<<new_X<<"\n");
    IntervalVector reduced_X = intersection(X,new_X);
    ARIADNE_LOG(5,"reduced_X="<<reduced_X<<"\n");

    if(subset(new_X,X)) { return true; }
    else { return indeterminate; }
}



// FIXME: Look at this code again, especially relating to generalised Lagrange multipliers
Bool OptimiserBase::
is_infeasibility_certificate(IntervalVector d, IntervalVectorFunction g, IntervalVector c, FloatVector lambda) const
{
    // Try to prove lambda.g(y) > 0
    const uint m=d.size();
    const uint n=c.size();
    VectorTaylorFunction tg(d,g);
    VectorTaylorFunction ti=VectorTaylorFunction::identity(d);
    ScalarTaylorFunction ts(d);
    for(uint i=0; i!=n; ++i) {
        ts+=lambda[i]*(tg[i]-c[i].upper())+lambda[i+n]*(c[i].lower()-tg[i]);
    }
    for(uint i=0; i!=m; ++i) {
        ts+=lambda[2*n+i]*(ti[i]-d[i].upper())+lambda[2*n+m+i]*(ti[i]-d[i].lower());
    }
    Interval lambdagd=ts(d);
    ARIADNE_LOG(2,"  lambda="<<lambda<<"  lambda.g(D)="<<lambdagd<<"\n");
    if(lambdagd.lower()>0.0) {
        return true;
    } else {
        return false;
    }
}

IntervalVector OptimiserBase::
optimise(IntervalScalarFunction f, IntervalVector D, IntervalVectorFunction g, IntervalVector C) const
{
    IntervalVectorFunction h(0u,RealScalarFunction::constant(D.size(),0));
    return static_cast<const OptimiserInterface*>(this)->optimise(f,D,g,C,h);
}

IntervalVector OptimiserBase::
optimise(IntervalScalarFunction f, IntervalVector D, IntervalVectorFunction h) const
{
    IntervalVectorFunction g(0u,RealScalarFunction::constant(D.size(),0));
    IntervalVector C(0u,Interval(-infty,+infty));
    return static_cast<const OptimiserInterface*>(this)->optimise(f,D,g,C,h);
}



enum ConstraintKind { EQUALITY, UPPER_BOUNDED, LOWER_BOUNDED, BOUNDED };

inline ConstraintKind constraint_kind(Interval C) {
    if(C.lower()==C.upper()) { return EQUALITY; }
    else if(C.lower()==-infty) { return UPPER_BOUNDED; }
    else if(C.upper()==+infty) { return LOWER_BOUNDED; }
    else { return BOUNDED; }
}


Void NonlinearInteriorPointOptimiser::
initialise_lagrange_multipliers(const IntervalVector& D, const FloatVectorFunction& g, const IntervalVector& C,
                                const FloatVector& x, FloatVector& lambda) const
{
    lambda.resize(C.size());
    for(uint i=0; i!=C.size(); ++i) {
        if(C[i].lower()==C[i].upper()) { lambda[i]=0.0; }
        else if(C[i].lower()==-infty) { lambda[i]=-1.0; }
        else if(C[i].upper()==+infty) { lambda[i]=+1.0; }
        else { lambda[i] = 0.0; }
    }
}

Float NonlinearInteriorPointOptimiser::
compute_mu(const IntervalVector& D, const FloatVectorFunction& g, const IntervalVector& C,
           const FloatVector& x, const FloatVector& lambda) const
{
    // Compute the relaxation parameter mu as the average of the product of the Lyapunov exponents and constraint satisfactions

    Float mu = 0.0;

    FloatVector gx = g(x);

    for(uint i=0; i!=C.size(); ++i) {
        if(C[i].lower()==C[i].upper()) { }
        else if(C[i].lower()==-infty) { mu += lambda[i] * (gx[i] - C[i].upper()); }
        else if(C[i].upper()==+infty) { mu += lambda[i] * (gx[i] - C[i].lower()); }
        else { std::cerr<<"FIXME: Compute mu for bounded constraint\n";
            if (lambda[i] <=0.0) { mu += lambda[i] * (gx[i] - C[i].upper()); }
            else { mu += lambda[i] * (gx[i] - C[i].lower()); }
        }
    }

    mu /= C.size();

    return mu;
}


IntervalVector NonlinearInteriorPointOptimiser::
optimise(IntervalScalarFunction f, IntervalVector D, IntervalVectorFunction g, IntervalVector C, IntervalVectorFunction h) const
{
    assert(h.result_size()==0u);

    FloatVector x = midpoint(D);

    std::cerr<<"x="<<x<<" g(x)="<<g(IntervalVector(x))<<" C="<<C<<"\n";
    // FIXME: Temporarily assume x is feasible
    ARIADNE_ASSERT(subset(g(IntervalVector(x)),C));

    Vector<Float> lambda;
    this->initialise_lagrange_multipliers(D,g,C,x,lambda);
    std::cerr << "x="<<x<<" lambda="<<lambda<<"\n";

    for(uint i=0; i!=10; ++i) {
        this->optimisation_step(f,D,g,C,x,lambda);
    }

    std::cerr << "x="<<x<<" lambda="<<lambda<<" mu="<<compute_mu(D,g,C,x,lambda)<<"\ng(x)="<<g(x)<<" C="<<C<<"\nf(x)="<<f(x)<<"\n";
    return IntervalVector(x);
}


Tribool NonlinearInteriorPointOptimiser::
feasible(IntervalVector d, IntervalVectorFunction g, IntervalVector c) const
{
    ARIADNE_LOG(2,"NonlinearInteriorPointOptimiser::feasible(IntervalVector d, IntervalVectorFunction g, IntervalVector c)\n");
    ARIADNE_LOG(2,"  d="<<d<<", g="<<g<<", c="<<c<<"\n");

    ARIADNE_ASSERT(g.argument_size()==d.size());
    ARIADNE_ASSERT(g.result_size()==c.size());
    Float t;
    FloatVector x,y,z;

    this->setup_feasibility(d,g,c,x,y);

    // FIXME: Allow more steps
    for(uint i=0; i!=12; ++i) {
        ARIADNE_LOG(4,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", x="<<x<<", z="<<z<<"\n");
        this->feasibility_step(d,g,c,x,y);
        if(t>0) {
            ARIADNE_LOG(2,"  y="<<y<<", g(y)="<<g(y)<<"\n");
            if(this->is_feasible_point(d,g,c,y)) {
                return true;
            }
        }
    }
    ARIADNE_LOG(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<"\n");
    if(this->is_infeasibility_certificate(d,g,c,x)) {
        return false;
    }
    return indeterminate;
}

// See Hande Y. Benson, David F. Shanno, And Robert J. Vanderbei,
// "Interior-point methods for nonconvex nonlinear programming: Jamming and comparative numerical testing"
// For some of the terminology used


Void NonlinearInteriorPointOptimiser::
optimisation_step(const FloatScalarFunction& f, const IntervalVector& d, const FloatVectorFunction& g, const IntervalVector& c,
                  FloatVector& x, FloatVector& lambda) const
{
    const uint m=x.size();
    const uint n=lambda.size();

    ARIADNE_LOG(7,"x="<<x<<"\n");
    ARIADNE_LOG(7,"lambda="<<lambda<<"\n");

    FloatVector slack(2*n);
    FloatVectorRange slackl(slack,range(0,n));
    FloatVectorRange slacku(slack,range(n,2*n));

    Differential<Float> ddfx=f.evaluate(Differential<Float>::variables(2,x));
    Vector< Differential<Float> > ddgx=g.evaluate(Differential<Float>::variables(2,x));

    Float mu = this->compute_mu(d,g,c,x,lambda);
    ARIADNE_LOG(9,"mu="<<mu<<"\n");

    // G is the constraint value vector
    FloatVector gx = ddgx.value();
    ARIADNE_LOG(9,"g(x)="<<gx<<"\n");

    // A is the derivative matrix aij=dgi/dxj
    FloatVector gradfx = ddfx.gradient();
    ARIADNE_LOG(9,"Df(x)="<<gradfx<<"\n");

    // A is the derivative matrix aij=dgi/dxj
    FloatMatrix A = ddgx.jacobian();
    ARIADNE_LOG(9,"A="<<A<<"\n");

    // H is the Hessian matrix H[i1,i2] = df/dx[i1]dx[i2] + Sum_[j] lambda[j]*dg[j]/dx[i1]dx[i2]
    FloatMatrix H = ddfx.hessian();
    for(uint j=0; j!=n; ++j) { H += lambda[j] * ddgx[j].hessian(); }
    ARIADNE_LOG(9,"H="<<H<<"\n");


    // Determines the weighting to give to the relaxation parameter mu
    // for equality constraints relative to other constraints
    static const double EQUALITY_RELAXATION_MULTIPLIER = 1.0;

    FloatVector primal_residuals = gradfx + lambda * A; // The residuals df/dx[i] + Sum[j] dg[j]/dx[i] * l[j]
    ARIADNE_LOG(9,"rx="<<primal_residuals<<"\n");

    FloatVector dual_residuals(n);        // The residuals s(g(x),lambda,mu)
    FloatVector dual_inverse_hessian(n);  // The matrix (ds/dw) / (ds/dl) where w[j]=g[j](x)
    for(uint j=0; j!=n; ++j) {
        if(c[j].lower()==c[j].upper()) {
            // s(w,l,mu) = (w-c) + l * mu * erm
            dual_residuals[j] = (gx[j]-c[j].midpoint()) + lambda[j] * (mu * EQUALITY_RELAXATION_MULTIPLIER);
            dual_inverse_hessian[j] = 1.0 / (mu * EQUALITY_RELAXATION_MULTIPLIER);
        } else if(c[j].lower()==-infty) {
            // s(w,l,mu) = (w-cu) * l - mu
            dual_residuals[j] = (gx[j]-c[j].upper()) * lambda[j] - mu;
            dual_inverse_hessian[j] = lambda[j] / (c[j].upper() - gx[j]);
        }
        else if(c[j].upper()==+infty) {
            // s(w,l,mu) = (w-cl) * l - mu
            dual_residuals[j] = (gx[j]-c[j].lower()) * lambda[j] - mu;
            dual_inverse_hessian[j] = lambda[j] / (c[j].lower() - gx[j]);
        }
        else {
            // s(w,l,mu) = (w-cl) * (cu-w) * l - (2*w-cl-cu) * mu
            Float zr = (gx[j]-c[j].lower()) * (c[j].upper()-gx[j]);
            Float zm = 2 * gx[j] - ( c[j].lower() + c[j].upper() );
            dual_residuals[j] = zr * lambda[j] + zm * mu;
            dual_inverse_hessian[j] = ( -zm * lambda[j] + 2 * mu) / zr;
        }
    }
    ARIADNE_LOG(9,"rl="<<dual_residuals<<"\n");
    ARIADNE_LOG(9,"D="<<dual_inverse_hessian<<"\n");

    // Add corrections to H due to state constraints
    FloatVector xnu(m);
    for(uint i=0; i!=m; ++i) {
        // D = -2*mu*(r^2+(x-c)^2)/(r^2-(x-c)^2)^2
        Float rsqr = sqr(d[i].radius());
        Float xmincsqr = sqr(x[i] - d[i].midpoint());
        xnu[i] = -2*mu*(rsqr+xmincsqr)/sqr(rsqr-xmincsqr);
        H[i][i] += xnu[i];
    }
    ARIADNE_LOG(9,"xnu="<<xnu<<"\n");

    FloatMatrix S=H;
    FloatVector D=-dual_inverse_hessian;
    atda(S,A,D);
    ARIADNE_LOG(9,"S="<<S<<"\n");
    FloatMatrix Sinv=inverse(S);
    ARIADNE_LOG(9,"Sinv="<<Sinv<<"\n");

    FloatVector& rx=primal_residuals;
    FloatVector& rlambda=dual_residuals;

    // Compute directions to move primal and dual variables
    FloatVector dx = Sinv * FloatVector( rx - ( rlambda * FloatDiagonalMatrix(D) ) * A );
    FloatVector dlambda = FloatDiagonalMatrix(D) * FloatVector(rlambda - A * dx);

    // Compute distance to move variables preserving feasibility

    FloatVector newx(m);
    FloatVector newlambda(n);
    FloatVector newgx(n);

    Float alpha = 1.0;
    bool success = false;
    do {
        newx = x - alpha * dx;
        newlambda = lambda - alpha * dlambda;
        if(egtr(newlambda,0.0)) {
            newgx = g(newx);
            if(!egtr(newgx-lower_bounds(c),0.0)) { success = true; }
        }
        if(!success) { alpha*=0.75; }
    } while(!success);
    ARIADNE_LOG(9,"alpha="<<alpha<<"\n");

    x = newx;
    lambda = newlambda;

    ARIADNE_LOG(2,"\n");
}


IntervalVector NonlinearInteriorPointOptimiser::
optimise(IntervalScalarFunction f, IntervalVector d, IntervalVectorFunction h) const
{
    ARIADNE_PRECONDITION(f.argument_size()==d.size());
    ARIADNE_PRECONDITION(h.argument_size()==d.size());
    ARIADNE_PRECONDITION(!empty(d));
    FloatVector x = midpoint(d);
    FloatVector lambda(h.result_size(),0.0);

    ARIADNE_LOG(4,"NonlinearInteriorPointOptimiser::optimise(f,D,h)\n");
    ARIADNE_LOG(4,"f="<<f<<"\n");
    ARIADNE_LOG(4,"D="<<d<<"\n");
    ARIADNE_LOG(4,"h="<<h<<"\n");

    static const double MU_REDUCTION_FACTOR = 0.25;

    Float mu=MU_REDUCTION_FACTOR;
    for(uint i=0; i!=5; ++i) {
        for(uint i=0; i!=2; ++i) {
            this -> optimisation_step(f,d,h,x,lambda, mu);
        }
        mu *= MU_REDUCTION_FACTOR;
    }

    return x;
}

Void NonlinearInteriorPointOptimiser::
optimisation_step(const FloatScalarFunction& f, const IntervalVector& d, const FloatVectorFunction& h,
                  FloatVector& x, FloatVector& lambda, Float& mu) const
{
    ARIADNE_LOG(4,"NonlinearInteriorPointOptimiser::feasibility_step(...)\n");

    ARIADNE_DEBUG_ASSERT(x.size()==d.size());
    ARIADNE_DEBUG_ASSERT(f.argument_size()==x.size());
    ARIADNE_DEBUG_ASSERT(h.argument_size()==x.size());
    ARIADNE_DEBUG_ASSERT(h.result_size()==lambda.size());

    static const Float ALPHA_SCALE_FACTOR = 0.75;

    const uint n=x.size();
    const uint m=lambda.size();

    ARIADNE_LOG(7,"x="<<x<<"\n");
    ARIADNE_LOG(7,"lambda="<<lambda<<"\n");

    FloatDifferential ddfx=f.evaluate(FloatDifferential::variables(2,x));
    Vector<FloatDifferential> ddhx=h.evaluate(FloatDifferential::variables(2,x));

    FloatVector hx = ddhx.value();
    ARIADNE_LOG(9,"h(x)="<<hx<<"\n");

    FloatVector gradfx = ddfx.gradient();
    ARIADNE_LOG(9,"Df(x)="<<gradfx<<"\n");

    // B is the derivative matrix aij=dhi/dxj
    FloatMatrix B = ddhx.jacobian();
    ARIADNE_LOG(9,"B=Dh(x)="<<B<<"\n");

    ARIADNE_LOG(9,"DDf(x)="<<ddfx.hessian()<<"\n");
    ARIADNE_LOG(9,"DDh[0](x)="<<ddhx[0].hessian()<<"\n");
    // H is the Hessian matrix H[i1,i2] = df/dx[i1]dx[i2] + Sum_[j] lambda[j]*dg[j]/dx[i1]dx[i2]
    FloatMatrix H = ddfx.hessian();
    for(uint j=0; j!=m; ++j) { H += lambda[j] * ddhx[j].hessian(); }
    ARIADNE_LOG(9,"H="<<H<<"\n");

    ARIADNE_LOG(7,"mu="<<mu<<"\n");

    // Using a barriers log(x-xl)+log(xu-x) and penalties h(x)^2, we obtain the
    // max f(x) + mu * sum_i log(x_i-xl_i)+log(xu_i-x_i) - sum_k h_k(x)^2 / 2nu
    //
    // Differentiating gives d_if(x) + mu/(x_i-xl_i) - mu/(xu_i-x_i) - sum_k dh_k/nu = 0
    // Set lambda_k = -h_k/nu and then take nu->0
    //
    // The KKT conditions are
    //    d_if(x) + sum_k kappa_k d_ih_k + mu s_i(x_i) = 0
    //       h_k(x)   = 0
    //   xl_i <= x_i <= xu_i
    //
    // The Jacobian of the KKT conditions is
    //   ( H + ds(X)    dh^T )
    //   (  dh           0   )

    // Compute residual rx = df + lambda * dg + mu*(1/xu-x - 1/x-xl)
    FloatMatrix& S=H;
    const FloatVector& rlambda = hx;
    FloatVector rx = gradfx + lambda * B;

    ARIADNE_LOG(9,"df+lambda*dh="<<rx<<"\n");
    for(uint i=0; i!=n; ++i) {
        Float reczil = rec(x[i]-d[i].lower());
        Float recziu = rec(d[i].upper() - x[i]);
        rx[i] += mu * ( recziu - reczil);
        S[i][i] += mu * (reczil*reczil + recziu*recziu);
    }

    ARIADNE_LOG(9,"S="<<S<<"\n");
    ARIADNE_LOG(9,"B="<<B<<"\n");
    ARIADNE_LOG(9,"rx="<<rx<<"\n");
    ARIADNE_LOG(9,"rlambda="<<rlambda<<"\n");

    FloatMatrix Sinv = inverse(S);
    ARIADNE_LOG(9,"inverse(S)="<<Sinv<<"\n");
    FloatMatrix BSinvBT = (B*Sinv)*transpose(B);
    ARIADNE_LOG(9,"BSinvBT="<<BSinvBT<<"\n");
    ARIADNE_LOG(9,"inverse(BSinvBT)="<<inverse(BSinvBT)<<"\n");

    ARIADNE_LOG(9,"BSinvBT*(B*rx)="<<(BSinvBT*(B*rx))<<"\n");
    FloatVector dlambda = inverse(BSinvBT) * (B*(Sinv*rx)-rlambda);
    ARIADNE_LOG(9,"dlambda="<<dlambda<<"\n");
    FloatVector dx = Sinv * (rx - dlambda*B);
    ARIADNE_LOG(9,"dx="<<dx<<"\n");

    ARIADNE_LOG(9,"S*dx+B^T*dlambda="<<S*dx+dlambda*B<<"; rx="<<rx<<"; ex="<<S*dx+dlambda*B-rx<<"\n");
    ARIADNE_LOG(9,"B*dx="<<B*dx<<"; rlambda="<<rlambda<<"; elambda="<<B*dx-rlambda<<"\n");


    Float alpha = 1.0;
    FloatVector newx(x.size());
    do {
        newx = x - alpha * dx;
        alpha *= ALPHA_SCALE_FACTOR;
    } while(!contains(d,newx));
    alpha /= ALPHA_SCALE_FACTOR;

    newx = x - alpha * dx;
    FloatVector newlambda = lambda - alpha * dlambda;

    ARIADNE_LOG(9,"alpha="<<alpha<<"\n");
    ARIADNE_LOG(9,"new x="<<newx<<"\n");
    ARIADNE_LOG(9,"new lambda="<<newlambda<<"\n\n");

    ARIADNE_LOG(9,"old f(x)="<<f(x)<<"\n");
    ARIADNE_LOG(9,"old h(x)="<<h(x)<<"\n");
    ARIADNE_LOG(9,"new f(x)="<<f(newx)<<"\n");
    ARIADNE_LOG(9,"new h(x)="<<h(newx)<<"\n\n");

    x = newx;
    lambda = newlambda;
}


void
NonlinearInteriorPointOptimiser::feasibility_step(
    const IntervalVector& d, const FloatVectorFunction& g, const IntervalVector& c,
    FloatVector& x, FloatVector& y) const
{
}


void
NonlinearInteriorPointOptimiser::feasibility_step(
    const IntervalVector& d, const FloatVectorFunction& g, const IntervalVector& c,
    FloatVector& x, FloatVector& y, Float& t) const
{
    static const double infty = std::numeric_limits<double>::infinity();

    static const double gamma=1.0/1024;
    static const double sigma=1.0/8;
    static const double scale=0.75;

    const uint m=d.size();
    const uint n=c.size();

    FloatVector z(n);

    ARIADNE_ASSERT_MSG(g.argument_size()==m,"d="<<d<<" g="<<g);
    ARIADNE_ASSERT_MSG(g.result_size()==n,"d="<<d<<" g="<<g<<" c="<<c);
    ARIADNE_ASSERT(x.size()==m);
    ARIADNE_ASSERT(y.size()==n);

    Vector< Differential<Float> > ddgx=g.evaluate(Differential<Float>::variables(2,x));
    ARIADNE_LOG(9,"  ddgx="<<ddgx<<"\n");

    Vector<Float> gx = ddgx.value();
    ARIADNE_LOG(7," g(x)="<<gx<<" ");
    Matrix<Float> A = transpose(ddgx.jacobian());
    ARIADNE_LOG(7," A="<<A<<" ");

    // H is the Hessian matrix H of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) \lambda_k$
    Matrix<Float> H(m,m);
    for(uint i=0; i!=m; ++i) {
        H+=y[i]*ddgx[i].hessian();
    }
    ARIADNE_LOG(7," H="<<H<<" ");




    // Add correction for bounded domain to diagonal elements of Hessian
    for(uint i=0; i!=m; ++i) {
    }

    // Compute diagonal entries of KKT Hessian
    Vector<Float> D(n);
    for(uint j=0; j!=n; ++j) {
        if(c[j].lower()==c[j].upper()) {
        } else if(c[j].upper()==+infty) {
        } else if(c[j].lower()==-infty) {
        } else {
            ARIADNE_DEBUG_ASSERT(-infty<c[j].lower() && c[j].lower()<c[j].upper() && c[j].upper()<+infty);
        }
    }

    Float mu=dot(x,z)/m;
    if(!egtr(emul(x,z),gamma*mu)) {
        if(verbosity>=1) { ARIADNE_WARN("Near-degeneracy in Lyapunov multipliers in interior-point solver:\n  x="<<x<<", y="<<y<<", z="<<z<<"\n"); }
        x=(1-sigma)*x+FloatVector(x.size(),sigma/x.size());
        mu=dot(x,z)/m;
    }

    FloatVector yt=join(y,t);
    ARIADNE_LOG(9,"m="<<m<<" n="<<n<<"\n");
    ARIADNE_LOG(9,"x="<<x<<" yt="<<yt<<" z="<<z<<"\n");


    // Construct diagonal matrices
    FloatVector DE=ediv(x,z);
    ARIADNE_LOG(9,"  D="<<DE<<"\n");

    // Construct the extended valuation GY=(gy-cu+te,cl-gy+te,y-bu+te,bl-y+te)
    FloatVector gye(2*(m+n));
    //for(uint j=0; j!=n; ++j) { gxe[j]=gy[j]-c[j].upper()+t; gye[n+j]=c[j].lower()-gy[j]+t; }
    //for(uint i=0; i!=m; ++i) { gye[2*n+i]=y[i]-d[i].upper()+t; gye[2*n+m+i]=d[i].lower()-y[i]+t; }
    ARIADNE_LOG(9,"  GE="<<gye<<"\n");

    // Construct the extended matrix AE=(A -A I -I \\ e e 0 0)
    FloatMatrix AE(m+1,2*(m+n));
    //for(uint i=0; i!=m; ++i) { for(uint j=0; j!=n; ++j) { AE[i][j]=A[i][j]; AE[i][n+j]=-A[i][j]; } }
    //for(uint i=0; i!=m; ++i) { AE[i][2*n+i]=1; AE[i][2*n+m+i]=-1; }
    //for(uint k=0; k!=o; ++k) { AE[m][k]=1; }
    FloatMatrix AET=transpose(AE);

    // Construct the symmetric matrix and its inverse
    //FloatMatrix S(m+1,m+1); adat(S,AE,DE);
    //ARIADNE_LOG(9,"S="<<S<<"\n");
    //S=FloatMatrix(m+1,m+1); simple_adat(S,AE,DE);
    //ARIADNE_LOG(9,"S="<<S<<"\n");
    FloatMatrix S=feasibility_adat(H,A,DE);
    ARIADNE_LOG(9,"S="<<S<<"\n");
    FloatMatrix Sinv=inverse(S);
    ARIADNE_LOG(9,"Sinv="<<Sinv<<"\n");

    // FIXME: What if S is not invertible?

    // Construct the residuals
    FloatVector rx=esub(emul(x,z),mu*sigma);
    //FloatVector ryt=-prod(AE,x); ryt[m]+=1; // FIXME: Need hessian
    FloatVector ryt=-feasibility_mul(A,x); ryt[m]+=1; // FIXME: Need hessian
    FloatVector rz=gye+z;
    ARIADNE_LOG(9,"rx="<<rx<<" ryt="<<ryt<<" rz="<<rz<<"\n");

    //FloatVector rr=prod(AE,ediv(FloatVector(rx-emul(x,rz)),z))-ryt;
    FloatVector rr=ryt+prod(AE,ediv(FloatVector(rx-emul(x,rz)),z))-ryt;


    // Compute the differences
    FloatVector dyt=prod(Sinv,rr);
    //FloatVector dz=-rz-prod(AET,dyt);
    FloatVector dz=-rz-feasibility_trmul(A,dyt);
    FloatVector dx=-ediv(FloatVector(rx+emul(x,dz)),z);
    ARIADNE_LOG(9,"dx="<<dx<<" dyt="<<dyt<<" dz="<<dz<<"\n");

    FloatVector nx,ny,nyt,nz; Float nt;

    // Since we need to keep the point feasible, but the updates are linear
    // we need to validate feasibility directly rather than assuming the
    // linear update of y and z are good enough.
    bool allpositive=false;
    Float alpha=1/scale;
    if(!egtr(emul(x,z) , gamma*mu/16)) {
        ARIADNE_LOG(1,"WARNING: x="<<x<<", z="<<z<< ", x.z="<<emul(x,z)<<"<"<<gamma*mu / 16);
        throw NearBoundaryOfFeasibleDomainException();
    }
    while(!allpositive) {
        alpha=alpha*scale;
        nx=x+alpha*dx;
        nyt=yt+alpha*dyt;
        ny=project(nyt,range(0,m));
        nt=nyt[m];
        //NonlinearInteriorPointOptimiser::compute_z(d,g,c,ny,nt,nz);
        allpositive = egtr(nx,0.0) && egtr(nz,0.0) && egtr(emul(nx,nz),gamma*mu);
    }
    ARIADNE_LOG(9,"alpha="<<alpha<<"\n");
    ARIADNE_LOG(9,"nx="<<nx<<" nyt="<<nyt<<" nz="<<nz<<" nxz="<<eivl(emul(nx,nz))<<"\n");

    x=nx; y=project(nyt,range(0,m)); z=nz; t=nyt[m];
}

/*
Void NonlinearInteriorPointOptimiser::linearised_feasibility_step(
    const IntervalVector& d, const FloatVectorFunction& g, const IntervalVector& c,
    Float& t, FloatVector& x, FloatVector& y) const
{
    static const double gamma=1.0/1024;
    static const double sigma=1.0/8;
    static const double scale=0.75;

    const uint m=d.size();
    const uint n=c.size();
    const uint o=2*(m+n);

    ARIADNE_ASSERT_MSG(g.argument_size()==m,"d="<<d<<" g="<<g);
    ARIADNE_ASSERT_MSG(g.result_size()==n,"d="<<d<<" g="<<g<<" c="<<c);
    ARIADNE_ASSERT(x.size()==o);
    ARIADNE_ASSERT(y.size()==m);

    FloatVector z(o);

    FloatVector yt=join(y,t);
    ARIADNE_LOG(9,"m="<<m<<" n="<<n<<"\n");
    ARIADNE_LOG(9,"x="<<x<<" yt="<<yt<<" z="<<z<<"\n");

    Float mu=dot(x,z)/o;

    Vector< Differential<Float> > dg=g.evaluate(Differential<Float>::variables(1,y));
    ARIADNE_LOG(9,"  dg="<<dg<<"\n");

    // gy is the vector of values of g(y)
    FloatVector gy(n); for(uint j=0; j!=n; ++j) { gy[j]=dg[j].value(); }
    ARIADNE_LOG(9,"  g(y)="<<gy<<" ");

    // A is the transpose derivative matrix aij=dgj/dyi, extended with a column of ones
    FloatMatrix A(m,n);
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=n; ++j) {
            A[i][j]=dg[j][i];
        }
    }
    ARIADNE_LOG(9," A="<<A<<" ");

    // H is the Hessian matrix Hik = (xcuj-xclj)*dgj/dyidyk
    FloatMatrix H=FloatMatrix::zero(m,m);
    ARIADNE_LOG(9," H="<<H);

    // Construct diagonal matrices
    FloatVector DE=ediv(x,z);
    ARIADNE_LOG(9,"  D="<<DE<<"\n");

    // Construct the extended valuation GY=(gy-cu+te,cl-gy+te,y-bu+te,bl-y+te)
    FloatVector gye(o);
    for(uint j=0; j!=n; ++j) { gye[j]=gy[j]-c[j].upper()+t; gye[n+j]=c[j].lower()-gy[j]+t; }
    for(uint i=0; i!=m; ++i) { gye[2*n+i]=y[i]-d[i].upper()+t; gye[2*n+m+i]=d[i].lower()-y[i]+t; }
    ARIADNE_LOG(9,"  GE="<<gye<<"\n");

    // Construct the extended matrix AE=(A -A I -I \\ e e 0 0)
    FloatMatrix AE(m+1,o);
    for(uint i=0; i!=m; ++i) { for(uint j=0; j!=n; ++j) { AE[i][j]=A[i][j]; AE[i][n+j]=-A[i][j]; } }
    for(uint i=0; i!=m; ++i) { AE[i][2*n+i]=1; AE[i][2*n+m+i]=-1; }
    for(uint k=0; k!=o; ++k) { AE[m][k]=1; }
    for(uint k=0; k!=2*(m+n); ++k) { AE[m][k]=1; }
    FloatMatrix AET=transpose(AE);

    // Construct the symmetric matrix and its inverse
    //FloatMatrix S(m+1,m+1); adat(S,AE,DE);
    //ARIADNE_LOG(9,"S="<<S<<"\n");
    //S=FloatMatrix(m+1,m+1); simple_adat(S,AE,DE);
    //ARIADNE_LOG(9,"S="<<S<<"\n");
    FloatMatrix S=feasibility_adat(H,A,DE);
    ARIADNE_LOG(9,"S="<<S<<"\n");
    FloatMatrix Sinv=inverse(S);
    ARIADNE_LOG(9,"Sinv="<<Sinv<<"\n");

    // FIXME: What if S is not invertible?

    // Construct the residuals
    FloatVector rx=esub(emul(x,z),mu*sigma);
    //FloatVector ryt=-prod(AE,x); ryt[m]+=1; // FIXME: Need hessian
    FloatVector ryt=-feasibility_mul(A,x); ryt[m]+=1; // FIXME: Need hessian
    FloatVector rz=gye+z;
    ARIADNE_LOG(9,"rx="<<rx<<" ryt="<<ryt<<" rz="<<rz<<"\n");

    //FloatVector rr=prod(AE,ediv(FloatVector(rx-emul(x,rz)),z))-ryt;
    FloatVector rr=ryt+prod(AE,ediv(FloatVector(rx-emul(x,rz)),z))-ryt;


    // Compute the differences
    FloatVector dyt=prod(Sinv,rr);
    //FloatVector dz=-rz-prod(AET,dyt);
    FloatVector dz=-rz-feasibility_trmul(A,dyt);
    FloatVector dx=-ediv(FloatVector(rx+emul(x,dz)),z);
    ARIADNE_LOG(9,"dx="<<dx<<" dyt="<<dyt<<" dz="<<dz<<"\n");

    FloatVector nx,ny,nyt,nz; Float nt;

    // Since we need to keep the point feasible, but the updates are linear
    // we need to validate feasibility directly rather than assuming the
    // linear update of y and z are good enough.
    bool allpositive=false;
    Float alpha=1/scale;
    ARIADNE_ASSERT_MSG(egtr(emul(x,z),gamma*mu),emul(x,z)<<"<"<<gamma*mu);
    while(!allpositive) {
        alpha=alpha*scale;
        nx=x+alpha*dx;
        nyt=yt+alpha*dyt;
        ny=project(nyt,range(0,m));
        nt=nyt[m];
        //NonlinearInteriorPointOptimiser::compute_z(d,g,c,ny,nt,nz);
        allpositive = egtr(nx,0.0) && egtr(nz,0.0) && egtr(emul(nx,nz),gamma*mu);
    }
    ARIADNE_LOG(9,"alpha="<<alpha<<"\n");
    ARIADNE_LOG(9,"nx="<<nx<<" nyt="<<nyt<<" nz="<<nz<<" nxz="<<eivl(emul(nx,nz))<<"\n");

    x=nx; y=project(nyt,range(0,m)); z=nz; t=nyt[m];

}
*/


void NonlinearInteriorPointOptimiser::
setup_feasibility(const IntervalVector& d, const FloatVectorFunction& g, const IntervalVector& c,
                  FloatVector& x, FloatVector& y) const
{
    const uint l=2*(d.size()+c.size());
    y=midpoint(d);
    x=FloatVector(l,1.0/l);
    //compute_tz(d,g,c,y,t,z);
}


PenaltyFunctionOptimiser* PenaltyFunctionOptimiser::
clone() const
{
    return new PenaltyFunctionOptimiser(*this);
}

IntervalVector PenaltyFunctionOptimiser::
optimise(IntervalScalarFunction f, IntervalVector D, IntervalVectorFunction g, IntervalVector C, IntervalVectorFunction h) const
{
    ARIADNE_NOT_IMPLEMENTED;
    return D;
}

Tribool PenaltyFunctionOptimiser::
feasible(IntervalVector D, IntervalVectorFunction g, IntervalVector C) const
{
    ARIADNE_LOG(3,"D="<<D<<" g="<<g<<" C="<<C<<" \n");
    RealVectorFunction h(0u,D.size());
    return this->feasible(D,g,C,h);
}

Tribool PenaltyFunctionOptimiser::
feasible(IntervalVector D, IntervalVectorFunction g, IntervalVector C, IntervalVectorFunction h) const
{
    ARIADNE_LOG(3,"D="<<D<<" g="<<g<<" C="<<C<<" h="<<h<<"\n");

    FloatVector x=midpoint(D);
    FloatVector w=midpoint(C);
    Float mu=1.0;

    ARIADNE_LOG(5,"x="<<x<<" w="<<w<<" mu="<<mu<<"\n");

    for(uint i=0; i!=10; ++i) {
        this->feasibility_step(D,g,C,h,x,w,mu);
    }
    return indeterminate;
}

Void PenaltyFunctionOptimiser::
feasibility_step(const IntervalVector& X, const FloatVectorFunction& g, const IntervalVector& W, const FloatVectorFunction& h,
                 FloatVector& x, FloatVector& w, Float& mu) const
{
    const uint n=X.size();
    const uint m=W.size();
    const uint l=h.result_size();

    ARIADNE_LOG(4,"PenaltyFunctionOptimiser::feasibility_step(...)\n");
    ARIADNE_LOG(5,"x="<<x<<"\n");
    ARIADNE_LOG(5,"w="<<w<<"\n");

    Vector<FloatDifferential> ddgx=g.evaluate(FloatDifferential::variables(2,x));
    Vector<FloatDifferential> ddhx=h.evaluate(FloatDifferential::variables(2,x));

    mu *= 0.5;
    ARIADNE_LOG(9,"mu="<<mu<<"\n");

    // G is the constraint value vector
    FloatVector gx = ddgx.value();
    FloatVector hx = ddhx.value();
    ARIADNE_LOG(9,"g(x)="<<gx<<"\n");
    ARIADNE_LOG(9,"h(x)="<<hx<<"\n");

    // A is the transpose derivative matrix aij=dgi/dxj
    FloatMatrix A = transpose(ddgx.jacobian());
    ARIADNE_LOG(9,"A=Dg(x)="<<A<<"\n");
    FloatMatrix B = transpose(ddhx.jacobian());
    // FIXME: Due to problems with zero-element differential, need to resize matrix if no h
    if(l==0) { B.resize(n,0); }
    ARIADNE_LOG(9,"B=Dh(x)="<<B<<"\n");

    // H is the Hessian matrix H[i1,i2] = df/dx[i1]dx[i2] + Sum_[j] lambda[j]*dg[j]/dx[i1]dx[i2]
    FloatMatrix H(n,n);
    for(uint j=0; j!=m; ++j) { H += (gx[j]-w[j]) * ddgx[j].hessian(); }
    for(uint k=0; k!=l; ++k) { H += (hx[k]) * ddhx[k].hessian(); }
    ARIADNE_LOG(9,"H="<<H<<"\n");

    FloatDiagonalMatrix D(n);
    FloatDiagonalMatrix E(m);
    for(uint i=0; i!=n; ++i) { D[i] = rec(sqr(x[i]-X[i].lower())) + rec(sqr(X[i].upper()-x[i])); }
    for(uint j=0; j!=m; ++j) { E[j] = rec(sqr(w[j]-W[j].lower())) + rec(sqr(W[j].upper()-w[j])); }
    ARIADNE_LOG(9,"D="<<D<<"\n");
    ARIADNE_LOG(9,"E="<<E<<"\n");

    FloatMatrix S = H + B * transpose(B);
    S += D;
    ARIADNE_LOG(9,"S="<<S<<"\n");

    FloatMatrix R=inverse(S);
    ARIADNE_LOG(9,"inverse(S)="<<R<<"\n");

    // Compute residuals
    FloatVector rx = A*gx + B * hx ; // + 1/(x.upper()-x) + 1/x.lower()-x if no regularisation
    FloatVector rw = w-gx;

    ARIADNE_LOG(9,"rx="<<rx<<"\n");
    ARIADNE_LOG(9,"rw="<<rw<<"\n");

    FloatVector dx = R * (rx + A * rw);
    FloatVector dw = rw + dx*A;
    ARIADNE_LOG(9,"dx="<<dx<<"\n");
    ARIADNE_LOG(9,"dw="<<dw<<"\n");


    FloatVector newx(n);
    FloatVector neww(m);

    static const Float ALPHA_SCALE_FACTOR = 0.75;

    Float alpha = 1.0;
    do {
        newx = x - alpha * dx;
        neww = w - alpha * dw;
        alpha *= ALPHA_SCALE_FACTOR;
    } while ( !contains(X,newx) || !contains(W,neww) );
    alpha /= ALPHA_SCALE_FACTOR;

    ARIADNE_LOG(9,"alpha="<<alpha<<"\n");

    ARIADNE_LOG(9,"newx="<<newx<<"\n");
    ARIADNE_LOG(9,"neww="<<neww<<"\n\n");

    x=newx;
    w=neww;

    return;
}



Void PenaltyFunctionOptimiser::
feasibility_step(const IntervalVector& X, const IntervalVectorFunction& g, const IntervalVector& W, const IntervalVectorFunction& h,
                 IntervalVector& x, IntervalVector& w) const
{
}

/*
void NonlinearInteriorPointOptimiser::
compute_tz(const IntervalVector& d, const FloatVectorFunction& g, const IntervalVector& b,
           const FloatVector& y, Float& t, FloatVector& z) const
{
    static const double ZMIN=0.5;

    const uint m=g.argument_size();
    const uint n=g.result_size();

    FloatVector gy=g(y);

    t=+inf<Float>();
    for(uint j=0; j!=n; ++j) {
        t=min(t,b[j].upper()-gy[j]);
        t=min(t,gy[j]-b[j].lower());
    }
    for(uint i=0; i!=m; ++i) {
        t=min(t,d[i].upper()-y[i]);
        t=min(t,y[i]-d[i].lower());
    }

    // Ensures all z start out strictly positive for interior point method
    // TODO: Find a good initialization for t
    if(t>0.0) { t/=2; }
    else { t-=ZMIN; }

    z.resize(2*(m+n));
    for(uint j=0; j!=n; ++j) {
        z[j]=b[j].upper()-gy[j]-t;
        z[n+j]=gy[j]-b[j].lower()-t;
    }
    for(uint i=0; i!=m; ++i) {
        z[2*n+i]=d[i].upper()-y[i]-t;
        z[2*n+m+i]=y[i]-d[i].lower()-t;
    }
}

void NonlinearInteriorPointOptimiser::compute_z(const IntervalVector& d, const FloatVectorFunction& g, const IntervalVector& b,
                                                const FloatVector& y, const Float& t, FloatVector& z) const
{
    const uint m=g.argument_size();
    const uint n=g.result_size();

    FloatVector gy=g(y);

    z.resize(2*(m+n));
    for(uint j=0; j!=n; ++j) {
        z[j]=b[j].upper()-gy[j]-t;
        z[n+j]=gy[j]-b[j].lower()-t;
    }
    for(uint i=0; i!=m; ++i) {
        z[2*n+i]=d[i].upper()-y[i]-t;
        z[2*n+m+i]=y[i]-d[i].lower()-t;
    }
}
*/







/*

struct KuhnTuckerFunctionBody : VectorFunctionMixin<KuhnTuckerFunctionBody,Interval>
{
    IntervalScalarFunction f;
    Array<IntervalScalarFunction> g;
    Array<IntervalScalarFunction> df;
    Array<Array<IntervalScalarFunction> > dg;

    KuhnTuckerFunctionBody(IntervalScalarFunction _f, IntervalVectorFunction _g) {
        ARIADNE_ASSERT(_f.argument_size()==_g.argument_size());
        const uint m=_g.argument_size();
        const uint n=_g.result_size();
        g.resize(n); df.resize(m); dg.resize(n); for(uint j=0; j!=n; ++j) { dg[j].resize(m); }
        f=_f;
        for(uint j=0; j!=n; ++j) { g[j]=_g[j]; }
        for(uint i=0; i!=m; ++i) { df[i]=f.derivative(i); }
        for(uint j=0; j!=n; ++j) { for(uint i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    uint result_size() const { return g.size()*2+f.argument_size(); }
    uint argument_size() const { return g.size()*2+f.argument_size(); }
    IntervalScalarFunction operator[](uint) const { ARIADNE_NOT_IMPLEMENTED; }
    std::ostream& write(std::ostream&) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const uint m=f.argument_size();
        const uint n=g.size();
        Vector<X> x(project(arg,range(0,n)));
        Vector<X> y(project(arg,range(n,n+m)));
        Vector<X> z(project(arg,range(n+m,n+m+n)));
        Vector<X> rx(m), rz(n), rs(n);
        for(uint i=0; i!=m; ++i) { rx[i]=df[i].evaluate(y); for(uint j=0; j!=n; ++j) { rx[i]=rx[i]-x[j]*dg[j][i].evaluate(y); } }
        for(uint j=0; j!=n; ++j) { rz[j]=g[j].evaluate(y) + z[j]; }
        for(uint j=0; j!=n; ++j) { rs[j]=x[j]*z[j]; }
        project(res,range(0,n))=rz;
        project(res,range(n,n+m))=rx;
        project(res,range(n+m,n+m+n))=rs;
    }
};

struct FeasibilityKuhnTuckerFunctionBody : VectorFunctionMixin<FeasibilityKuhnTuckerFunctionBody,Interval>
{
    Array<IntervalScalarFunction> g;
    Array<Array<IntervalScalarFunction> > dg;

    FeasibilityKuhnTuckerFunctionBody(IntervalVectorFunction _g) {
        const uint m=_g.argument_size();
        const uint n=_g.result_size();
        g.resize(n); dg.resize(n); for(uint j=0; j!=n; ++j) { dg[j].resize(m); }
        for(uint j=0; j!=n; ++j) { g[j]=_g[j]; for(uint i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    uint result_size() const { return g.size()*2+g[0].argument_size()+1; }
    uint argument_size() const { return g.size()*2+g[0].argument_size()+1; }
    IntervalScalarFunction operator[](uint) const { ARIADNE_NOT_IMPLEMENTED; }
    std::ostream& write(std::ostream&) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class X> void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const uint m=g[0].argument_size();
        const uint n=g.size();
        Vector<X> x(project(arg,range(0,n)));
        Vector<X> y(project(arg,range(n,n+m)));
        Vector<X> z(project(arg,range(n+m,n+m+n)));
        X t(arg[n+m+n]);
        Vector<X> rx(m), rz(n), rs(n); X rt;
        for(uint i=0; i!=m; ++i) { rx[i]=x[0]*dg[0][i].evaluate(y); for(uint j=1; j!=n; ++j) { rx[i]=rx[i]+x[j]*dg[j][i].evaluate(y); } }
        for(uint j=0; j!=n; ++j) { rz[j]=g[j].evaluate(y) + t + z[j]; }
        for(uint j=0; j!=n; ++j) { rs[j]=x[j]*z[j]; }
        rt=1-x[0]; for(uint j=1; j!=n; ++j) { rt=rt-x[j]; }
        project(res,range(0,n))=rz;
        project(res,range(n,n+m))=rx;
        project(res,range(n+m,n+m+n))=rs;
        res[n+m+n]=rt;
    }
};



struct ConstrainedFeasibilityKuhnTuckerFunctionBody : VectorFunctionMixin<FeasibilityKuhnTuckerFunctionBody,Interval>
{
    uint m;
    uint n;
    IntervalVector d;
    Array<IntervalScalarFunction> g;
    IntervalVector c;
    Array<Array<IntervalScalarFunction> > dg;

    ConstrainedFeasibilityKuhnTuckerFunctionBody(IntervalVector D, IntervalVectorFunction _g, IntervalVector C) {
        m=_g.argument_size();
        n=_g.result_size();
        d=D; c=C;
        g.resize(n); dg.resize(n); for(uint j=0; j!=n; ++j) { dg[j].resize(m); }
        for(uint j=0; j!=n; ++j) { g[j]=_g[j]; for(uint i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    uint result_size() const { return 5*m+4*n+1u; }
    uint argument_size() const { return 5*m+4*n+1u; }
    IntervalScalarFunction operator[](uint) const { ARIADNE_NOT_IMPLEMENTED; }
    std::ostream& write(std::ostream& os) const { return os << "KuhnTuckerFunctionBody"; }

    template<class X> void _compute(Vector<X>& res, const Vector<X>& arg) const {
        const X zero=arg[0]*0.0;
        const uint l=2*(m+n);
        assert(arg.size()==l+m+l+1);
        Vector<X> x(project(arg,range(0u,l)));
        Vector<X> y(project(arg,range(l,l+m)));
        Vector<X> z(project(arg,range(l+m,l+m+l)));
        X t(arg[l+m+l]);
        Vector<X> rx(m,zero), rz(l,zero), rs(l,zero); X rt(zero);
        Vector<X> gy(n);
        for(uint j=0; j!=n; ++j) { gy[j]=g[j].evaluate(y); }
        Matrix<X> dgy(n,m);
        for(uint i=0; i!=m; ++i) { for(uint j=0; j!=n; ++j) { dgy[j][i]=dg[j][i].evaluate(y); } }

        for(uint i=0; i!=m; ++i) {
            for(uint j=0; j!=n; ++j) { rx[i]+=x[j]*(dgy[j][i]-c[j].upper()); rx[i]+=x[n+j]*(c[j].lower()-dgy[j][i]); }
            rx[i]+=x[2*n+i]*(y[i]-d[i].upper())-x[2*n+m+i]*(d[i].lower()-y[i]);
        }
        for(uint j=0; j!=n; ++j) { rz[j]=gy[j] + t + z[j]; rz[n+j]=t+z[n+j]-gy[j]; }
        for(uint i=0; i!=m; ++i) { rz[2*n+i]=y[i]+t+z[2*n+i]; rz[2*n+m+i]=y[i]+t+z[2*n+m+i]; }
        for(uint k=0; k!=l; ++k) { rs[k]=x[k]*z[k]; }
        rt+=1.0; for(uint j=0; j!=2*n; ++j) { rt=rt-x[j]; }
        project(res,range(0,l))=rz;
        project(res,range(l,l+m))=rx;
        project(res,range(l+m,l+m+l))=rs;
        res[l+m+l]=rt;
    }
};



IntervalVector KrawczykOptimiser::
optimise(IntervalScalarFunction f, IntervalVector d, IntervalVectorFunction g, IntervalVector c) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


Tribool KrawczykOptimiser::
feasible(IntervalVector d, IntervalVectorFunction g, IntervalVector c) const
{
    ARIADNE_LOG(2,"KrawczykOptimiser::feasible(IntervalVector d, IntervalVectorFunction g, IntervalVector c)\n");
    ARIADNE_LOG(2,"  d="<<d<<", g="<<g<<", c="<<c<<"\n");

    ARIADNE_ASSERT(g.argument_size()==d.size());
    ARIADNE_ASSERT(g.result_size()==c.size());

    Interval t; IntervalVector x,y,z;
    setup_feasibility(d,g,c,x,y,z,t);


    // FIXME: Allow more steps
    for(uint i=0; i!=12; ++i) {
        ARIADNE_LOG(4,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", x="<<x<<", z="<<z<<"\n");
        try {
            this->feasibility_step(d,g,c,x,y,z,t);
        }
        catch(SingularMatrixException) {
            return indeterminate;
        }
        if(t.lower()>t.upper()) {
            ARIADNE_LOG(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c<<"\n");
            return indeterminate;
        }
        if(t.lower()>0.0) {
            ARIADNE_LOG(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c<<"\n");
            if(this->is_feasible_point(d,g,c,midpoint(y))) {
                return true;
            }
        }
        if(t.upper()<0.0) {
            ARIADNE_LOG(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c<<"\n");
            return false;
        }
    }
    ARIADNE_LOG(2,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", d="<<d<<", c="<<c<<"\n");
    if(this->is_infeasibility_certificate(d,g,c,midpoint(x))) {
        return false;
    }
    return indeterminate;
}



void KrawczykOptimiser::setup_feasibility(const IntervalVector& d, const IntervalVectorFunction& g, const IntervalVector& c,
                                          IntervalVector& x, IntervalVector& y, IntervalVector& z, Interval& t) const
{
    const uint m=g.argument_size();
    const uint n=g.result_size();
    const uint l=2*(m+n);
    x=IntervalVector(l, Interval(0,1)/l);
    y=d;
    z.resize(2*(m+n));
    compute_tz(d,g,c,y,t,z);
}


void KrawczykOptimiser::compute_tz(const IntervalVector& d, const IntervalVectorFunction& g, const IntervalVector& c,
                                   const IntervalVector& y, Interval& t, IntervalVector& z) const
{
    ARIADNE_ASSERT(d.size()>0u);
    //static const double EPS=1.0/8;
    static const float min_float=std::numeric_limits<float>::min();

    const uint m=g.argument_size();
    const uint n=g.result_size();

    // Compute the image of y under the constraint function
    IntervalVector gy=g(y);
    gy+=IntervalVector(gy.size(),Interval(-min_float,+min_float));
    IntervalVector my=midpoint(y);
    IntervalVector mgy=g(my);

    // Find the range of possible values of the optimal t
    // This range is too pessimistic
    t=Interval(+inf<Float>(),+inf<Float>());
    for(uint j=0; j!=n; ++j) {
        t=min(t,c[j]-gy[j]);
        t=min(t,gy[j]-c[j]);
    }
    for(uint i=0; i!=m; ++i) {
        t=min(t,d[i]-y[i]);
        t=min(t,y[i]-d[i]);
    }

    // Find the range of possible values of the optimal t
    Float tmin=+inf<Float>();
    Float tmax=+inf<Float>();
    for(uint j=0; j!=n; ++j) {
        tmax=min(tmax,sub_up(c[j].upper(),gy[j].lower()));
        tmax=min(tmax,sub_up(gy[j].upper(),c[j].lower()));
        tmin=min(tmin,sub_down(c[j].upper(),mgy[j].upper()));
        tmin=min(tmin,sub_down(mgy[j].lower(),c[j].lower()));
    }
    for(uint i=0; i!=m; ++i) {
        tmin=min(tmin,sub_up(d[i].upper(),y[i].lower()));
        tmax=min(tmax,sub_up(y[i].upper(),d[i].lower()));
        tmin=min(tmin,sub_down(d[i].upper(),my[i].upper()));
        tmax=min(tmax,sub_down(my[i].lower(),d[i].lower()));
    }
    tmin-=0.0625;
    t=Interval(tmin,tmax);


    // Find the range of possible values of the optimal z
    // This range is too pessimistic
    for(uint j=0; j!=n; ++j) {
        z[j]=max(c[j].upper()-gy[j]-t,0.0);
        z[n+j]=max(gy[j]-c[j].lower()-t,0.0);
    }
    for(uint i=0; i!=m; ++i) {
        z[2*n+i]=max(d[i].upper()-y[i]-t,0.0);
        z[2*n+m+i]=max(y[i]-d[i].lower()-t,0.0);
    }

    // Find the range of possible values of the optimal z
    // This range is too pessimistic
    for(uint j=0; j!=n; ++j) {
        z[j]=Interval(0.0,c[j].upper()-mgy[j].lower()-tmin);
        z[n+j]=Interval(0.0,mgy[j].upper()-c[j].lower()-tmin);
    }
    for(uint i=0; i!=m; ++i) {
        z[2*n+i]=Interval(0.0,d[i].upper()-my[i].lower()-tmin);
        z[2*n+m+i]=Interval(0.0,my[i].upper()-d[i].lower()-tmin);
    }

    ARIADNE_LOG(9,"  d="<<d<<", c="<<c<<", y="<<y<<", g(y)="<<gy<<", t="<<t<<", z="<<z<<"\n");

}


void KrawczykOptimiser::
optimisation_step(const IntervalScalarFunction& f, const IntervalVectorFunction& g,
                  IntervalVector& x, IntervalVector& y, IntervalVector& z) const
{
    const uint m=f.argument_size();
    const uint n=g.result_size();

    Differential<Interval> ddf=f.evaluate(Differential<Interval>::variables(2,y));
    Vector< Differential<Interval> > ddg=g.evaluate(Differential<Interval>::variables(2,y));

    IntervalMatrix H(m,m);
    set_hessian(H,ddf);
    for(uint j=0; j!=n; ++j) { add_hessian(H,-x[j],ddg[j]); }

    IntervalMatrix A(m,n);
    set_jacobian_transpose(A,ddg);

    ARIADNE_LOG(9,"f="<<f<<"\ng="<<g<<"\nx="<<x<<" y="<<y<<" z="<<z<<"\n");
    ARIADNE_LOG(9,"A="<<A<<"\nH="<<H<<"\n");

    ARIADNE_NOT_IMPLEMENTED;

}



void KrawczykOptimiser::feasibility_step(const IntervalVectorFunction& g,
                                         IntervalVector& x, IntervalVector& y, IntervalVector& z, Interval& t) const
{
    ARIADNE_NOT_IMPLEMENTED;
    const uint m=y.size();
    const uint n=x.size();

    Vector< Differential<Interval> > ddg=g.evaluate(Differential<Interval>::variables(2,y));

    // A is the transpose derivative matrix aij=dgj/dyi
    IntervalMatrix A(m,n);
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    ARIADNE_LOG(9,"A="<<A<<"\n");

    // H is the Hessian matrix Hik = xj*dgj/dyidyk
    IntervalMatrix H(m,m);
    for(uint j=0; j!=n; ++j) {
        add_hessian(H,x[j],ddg[j]);
    }
    ARIADNE_LOG(9," H="<<H<<"\n");

    FloatMatrix mA=midpoint(A);
    ARIADNE_LOG(9," mA="<<mA<<"\n");
    FloatMatrix mH=midpoint(H);
    ARIADNE_LOG(9," mH="<<mH<<"\n");

    FloatVector mD(n);
    for(uint j=0; j!=n; ++j) { mD[j]=midpoint(x[j])/midpoint(z[j]); }
    ARIADNE_LOG(9," mD="<<mD<<"\n");

    FloatMatrix& mS=mH;
    adat(mS,mA,mD);
    ARIADNE_LOG(9,"mS="<<mS<<"\n");
    FloatMatrix mSinv=inverse(mS);
    ARIADNE_LOG(9,"mSinv="<<mSinv<<"\n");
}

// Feasibility step for dual (inequality constrained) problem without using slack variables
// FIXME: Do we need a slackness parameter mu? Probably not; hopefully the infinities are kept in check...
// This method has the advantage of not needing to update the primal variables
void KrawczykOptimiser::feasibility_step(const IntervalVector& d, const IntervalVectorFunction& g, const IntervalVector& c,
                                         IntervalVector& y, Interval& t) const
{
    const uint m=d.size();
    const uint n=c.size();

    // Compute function values
    Vector< Differential<Interval> > ddg=g.evaluate(Differential<Interval>::variables(2,y));

    // gy is the vector of values of g(y)
    IntervalVector gy(n);
    for(uint j=0; j!=n; ++j) { gy[j]=ddg[j].value(); }

    // z is the vector of slack variables z[k]=cu[k]-gy[k]-t or z[k]=gy[k]-cl[k]-t
    IntervalVector z(2*(m+n));
    for(uint j=0; j!=n; ++j) { z[j]=d[j].upper()-gy[j]-t; z[n+j]=gy[j]-d[j].lower()-t; }
    for(uint i=0; i!=m; ++i) { z[i]=c[2*n+i].upper()-y[i]-t; z[2*n+m+i]=y[i]-c[i].lower()-t; }

    IntervalVector zr(2*(m+n));
    for(uint k=0; k!=2*(m+n); ++k) { zr[k]=1.0/z[k]; }

    IntervalVector D(2*(m+n));
    for(uint k=0; k!=2*(m+n); ++k) { D[k]=zr[k]*zr[k]; }

    // A is the transpose derivative matrix aij=dgj/dyi
    IntervalMatrix A(m,n);
    for(uint i=0; i!=m; ++i) { for(uint j=0; j!=n; ++j) { A[i][j]=ddg[j][i]; } }

    // A is the sum of scaled Hessian matrices hi1i2=zj*ddgj/dyi1yi2
    IntervalMatrix H(m,m);

    IntervalMatrix SE(m+1,m+1);
    // SE[0:m][0:m] is the matrix H+/-A(D1+D2)AT+(D3+D4) where D=z.z
    for(uint i1=0; i1!=m; ++i1) { for(uint i2=0; i2!=m; ++i2) { SE[i1][i2]=H[i1][i2];
        for(uint j=0; j!=n; ++j) { SE[i1][i2]+=A[i1][j]*(D[j]+D[n+j])*A[i2][j]; }
    } }
    for(uint i=0; i!=m; ++i) { SE[i][i]+=(D[2*n+i]+D[2*n+m+i]); }
    // SE[m][0:m]=SE[0:m][m] is the vector A(D1-D2)e+(D3-D4)e
    for(uint i=0; i!=m; ++i) { SE[i][m]=D[2*n+i]-D[2*n+m+i];
        for(uint j=0; j!=n; ++j) { SE[i][m]+=A[i][j]*(D[j]-D[n+j]); }
        SE[m][i]=SE[i][m];
    }
    // SE[m][m] is the scalar eT(D1+D2)e+eT(D3+D4)e
    for(uint k=0; k!=2*(m+n); ++k) { SE[m][m]+=D[k]; }

    // Vector of residuals
    IntervalVector re(m+1);
    for(uint i=0; i!=m; ++i) { re[i]+=(zr[2*n+i]-zr[2*n+m+i]);
        for(uint j=0; j!=n; ++j) { re[i]+=A[i][j]*(zr[j]-zr[n+j]); }
    }
    for(uint j=0; j!=n; ++j) { re[m]+=(zr[j]+zr[n+n]); }
    for(uint i=0; i!=m; ++i) { re[m]+=(zr[2*n+i]+zr[2*n+m+i]); }

    // Compute inverse Jacobian matrix
    IntervalMatrix JE;
    try {
        JE=inverse(midpoint(SE));
    }
    catch(const SingularMatrixException& e) {
        ARIADNE_WARN("Matrix S="<<midpoint(SE)<<" is not invertible");
        ARIADNE_LOG(1,"WARNING: Matrix S="<<midpoint(SE)<<" is not invertible");
        throw e;
    }

    // Krawczyk step
    IntervalVector dyt=prod(JE,IntervalVector(midpoint(re)))+prod(IntervalMatrix::identity(m+1)-prod(JE,SE),re-midpoint(re));

    // Extract y and t
    IntervalVector yt=join(y,t);
    IntervalVector nyt=yt+dyt;

    yt=intersection(yt,nyt);
    y=project(yt,range(0,m));
    t=yt[m];
}




void KrawczykOptimiser::feasibility_step(const IntervalVector& d, const IntervalVectorFunction& g, const IntervalVector& c,
                                         IntervalVector& x, IntervalVector& y, IntervalVector& z, Interval& t) const
{
    const uint m=d.size();
    const uint n=c.size();
    const uint o=2*(m+n);

    ARIADNE_ASSERT_MSG(g.argument_size()==m,"d="<<d<<" g="<<g);
    ARIADNE_ASSERT_MSG(g.result_size()==n,"d="<<d<<" g="<<g<<" c="<<c);
    ARIADNE_ASSERT(x.size()==o);
    ARIADNE_ASSERT(y.size()==m);
    ARIADNE_ASSERT(z.size()==o);

    IntervalVector yt=join(y,t);
    ARIADNE_LOG(9,"m="<<m<<" n="<<n<<"\n");
    ARIADNE_LOG(9,"x="<<x<<" yt="<<yt<<" z="<<z<<"\n");

    Vector< Differential<Interval> > ddg=g.evaluate(Differential<Interval>::variables(2,y));
    ARIADNE_LOG(9,"  ddg="<<ddg<<"\n");

    // gy is the vector of values of g(y)
    IntervalVector gy(n); for(uint j=0; j!=n; ++j) { gy[j]=ddg[j].value(); }
    ARIADNE_LOG(9,"  g(y)="<<gy<<" ");

    // A is the transpose derivative matrix aij=dgj/dyi, extended with a column of ones
    IntervalMatrix A(m,n);
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    ARIADNE_LOG(9," A="<<A<<" ");

    // H is the Hessian matrix Hik = (xcuj-xclj)*dgj/dyidyk
    IntervalMatrix H(m,m);
    for(uint j=0; j!=n; ++j) {
        add_hessian(H,x[j]-x[n+j],ddg[j]);
    }
    ARIADNE_LOG(9," H="<<H);

    // Construct the extended valuation GY=(gy-cu+te,cl-gy+te,y-bu+te,bl-y+te)
    IntervalVector gye(o);
    for(uint j=0; j!=n; ++j) { gye[j]=gy[j]-c[j].upper()+t; gye[n+j]=c[j].lower()-gy[j]+t; }
    for(uint i=0; i!=m; ++i) { gye[2*n+i]=y[i]-d[i].upper()+t; gye[2*n+m+i]=d[i].lower()-y[i]+t; }
    ARIADNE_LOG(9,"  GE="<<gye<<"\n");

    // Construct the extended matrix AE=(A -A I -I \\ e e 0 0)
    IntervalMatrix AE(m+1,o);
    for(uint i=0; i!=m; ++i) { for(uint j=0; j!=n; ++j) { AE[i][j]=A[i][j]; AE[i][n+j]=-A[i][j]; } }
    for(uint i=0; i!=m; ++i) { AE[i][2*n+i]=1; AE[i][2*n+m+i]=-1; }
    for(uint k=0; k!=o; ++k) { AE[m][k]=1; }
    IntervalMatrix AET=transpose(AE);

    FloatMatrix mA=midpoint(A);
    FloatMatrix mAE=midpoint(AE);
    FloatMatrix mAET=midpoint(AET);
    FloatMatrix mH=midpoint(H);
    FloatVector mx=midpoint(x);
    FloatVector myt=midpoint(yt);
    FloatVector mz=midpoint(z);
    FloatVector mDE=ediv(mx,mz);


    // Construct the symmetric matrix and its inverse
    //FloatMatrix S(m+1,m+1); adat(S,AE,DE);
    //ARIADNE_LOG(9,"S="<<S<<"\n");
    //S=FloatMatrix(m+1,m+1); simple_adat(S,AE,DE);
    //ARIADNE_LOG(9,"S="<<S<<"\n");
    FloatMatrix mS=feasibility_adat(mH,mA,mDE);
    ARIADNE_LOG(9,"mS="<<mS<<"\n");
    FloatMatrix mSinv=inverse(mS);
    ARIADNE_LOG(9,"mSinv="<<mSinv<<"\n");

    // FIXME: What if S is not invertible?

    // Construct the residuals
    IntervalVector rx=emul(mx,mz);
    //FloatVector ryt=-prod(AE,x); ryt[m]+=1; // FIXME: Need hessian
    IntervalVector ryt=-feasibility_mul(mA,mx); ryt[m]+=1; // FIXME: Need hessian
    IntervalVector rz=midpoint(gye)+mz;
    ARIADNE_LOG(9,"rx="<<rx<<" ryt="<<ryt<<" rz="<<rz<<"\n");

    // Construct the errors on the residuals ([M]-M)([x]-x)
    IntervalVector ex=x-mx;
    IntervalVector eyt=yt-myt;
    IntervalVector ez=z-mz;
    IntervalMatrix eA=A-mA;
    IntervalMatrix eH=H-mH;

    IntervalVector erx=2.0*emul(ex,ez);
    IntervalVector eryt=IntervalMatrix(AE-mAE)*ex;
    IntervalVector erz=IntervalMatrix(AET-mAET)*eyt;
    ARIADNE_LOG(9,"erx="<<erx<<" eryt="<<eryt<<" erz="<<erz<<"\n");

    rx+=2.0*emul(ex,ez);
    ryt+=IntervalMatrix(AE-mAE)*ex;
    rz+=IntervalMatrix(AET-mAET)*eyt;
    ARIADNE_LOG(9,"rx="<<rx<<" ryt="<<ryt<<" rz="<<rz<<"\n");

    //FloatVector rr=prod(AE,ediv(FloatVector(rx-emul(x,rz)),z))-ryt;

    // Compute the error differences
    IntervalVector erxdz=ediv(erx,mz);
    IntervalVector edyt=(mSinv*mAE)*erxdz + mSinv*eyt - (mSinv*(mAE*DiagonalMatrix<Float>(mDE))) * ez;
    IntervalVector edz=-erz-feasibility_trmul(mA,edyt);
    IntervalVector edx=-ediv(IntervalVector(erx+emul(mx,edz)),mz);
    ARIADNE_LOG(9,"edx="<<edx<<" edyt="<<edyt<<" edz="<<edz<<"\n");

    // Compute the error differences
    IntervalVector eerr=prod(mAE,ediv(esub(erx,emul(mx,erz)),mz))-eryt;
    ARIADNE_LOG(9,"  eerr="<<eerr<<"\n");
    IntervalVector eedyt=prod(mSinv,eerr);
    IntervalVector eedz=-erz-feasibility_trmul(mA,eedyt);
    IntervalVector eedx=-ediv(IntervalVector(erx+emul(mx,eedz)),mz);
    ARIADNE_LOG(9,"eedx="<<eedx<<" eedyt="<<eedyt<<" eedz="<<eedz<<"\n");


    // Compute the differences
    IntervalVector rr=prod(mAE,ediv(esub(rx,emul(mx,rz)),mz))-ryt;
    IntervalVector dyt=prod(mSinv,rr);
    IntervalVector dz=-rz-feasibility_trmul(mA,dyt);
    IntervalVector dx=-ediv(IntervalVector(rx+emul(mx,dz)),mz);
    ARIADNE_LOG(9,"dx="<<dx<<" dyt="<<dyt<<" dz="<<dz<<"\n\n");

    IntervalVector nx,ny,nyt,nz; Float nt;
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
