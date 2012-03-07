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

// See Hande Y. Benson, David F. Shanno, And Robert J. Vanderbei,
// "Interior-point methods for nonconvex nonlinear programming: Jamming and comparative numerical testing"
// For some of the terminology used

#include <limits>

#include "config.h"

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
#include <include/multi_index-noaliasing.h>
#include <include/constraint_solver.h>

namespace Ariadne {

inline Sweeper default_sweeper() { return Sweeper(); }

static const double error =  1e-2;

typedef Vector<Float> FloatVector;
typedef Matrix<Float> FloatMatrix;
typedef VectorRange<FloatVector> FloatVectorRange;
typedef DiagonalMatrix<Float> FloatDiagonalMatrix;

typedef Vector<Interval> IntervalVector;
typedef Matrix<Interval> IntervalMatrix;
typedef VectorRange<IntervalVector> IntervalVectorRange;
typedef DiagonalMatrix<Interval> IntervalDiagonalMatrix;

template<class X> inline
DiagonalMatrix<X> const& diagonal_matrix(const Vector<X>& v) {
    return reinterpret_cast<DiagonalMatrix<X>const&>(v);
}

template<class X> inline
bool epos(const Vector<X>& x) {
    for(uint i=0; i!=x.size(); ++i) { if(x[i]<=0) { return false; } } return true;
}

template<class X> inline
bool eneg(const Vector<X>& x) {
    for(uint i=0; i!=x.size(); ++i) { if(x[i]>=0) { return false; } } return true;
}

template<class X, class XX> inline
bool egtr(const Vector<X>& x, const XX& s) {
    for(uint i=0; i!=x.size(); ++i) { if(x[i]<=s) { return false; } } return true;
}

template<class X, class XX> inline
bool elss(const Vector<X>& x, const XX& s) {
    for(uint i=0; i!=x.size(); ++i) { if(x[i]>=s) { return false; } } return true;
}

template<class X> inline
Vector<X> eadd(const Vector<X>& x, const Vector<X>& y) {
    Vector<X> r(x.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=x[i]+y[i]; } return r;
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

template<class X> inline
Vector<X> ediv(const X& s, const Vector<X>& z) {
    Vector<X> r(z.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=s/z[i]; } return r;
}

inline
Interval eivl(const FloatVector& x) {
    ARIADNE_ASSERT(x.size()>0); Interval r(x[0]); for(uint i=0; i!=x.size(); ++i) { r=hull(r,x[i]); } return r;
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

// Compute S+=A^TDA, where D is diagonal and S is symmetric.
template<class X>
void atda(Matrix<X>& S, const Matrix<X>& A, const DiagonalMatrix<X>& D)
{
    atda(S,A,D.diagonal());
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



enum ConstraintKind { EQUALITY, UPPER_BOUNDED, LOWER_BOUNDED, BOUNDED };

inline ConstraintKind constraint_kind(Interval C) {
    if(C.lower()==C.upper()) { return EQUALITY; }
    else if(C.lower()==-inf) { return UPPER_BOUNDED; }
    else if(C.upper()==+inf) { return LOWER_BOUNDED; }
    else { return BOUNDED; }
}





Bool OptimiserBase::
almost_feasible_point(IntervalVector D, IntervalVectorFunction g, IntervalVector C, FloatVector x, Float error) const
{
    if(!contains(D,x)) { return false; }
    IntervalVector gx=g(IntervalVector(x));
    return subset(gx,C+Vector<Interval>(C.size(),Interval(-error,+error)));
}


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
    ARIADNE_ASSERT(g.result_size()>0);
    IntervalVectorFunction ge(equality_constraints.size(),g.argument_size());
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
    FloatMatrix fltL = FloatDiagonalMatrix(fltD)*FloatMatrix(transpose(fltA));
    ARIADNE_LOG(7,"L="<<fltL<<"\n");

    IntervalMatrix ivlS = ivlA * make_exact(fltL);
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
    IntervalVector new_X = x - make_exact(fltL) * (ivlR * (ge(x)-ce) );
    ARIADNE_LOG(5,"old_X="<<X<<"\n");
    ARIADNE_LOG(5,"new_X="<<new_X<<"\n");
    IntervalVector reduced_X = intersection(X,new_X);
    ARIADNE_LOG(5,"reduced_X="<<reduced_X<<"\n");

    if(subset(new_X,X)) { return true; }
    else { return indeterminate; }
}




Bool OptimiserBase::
validate_feasibility(IntervalVector D, IntervalVectorFunction g, IntervalVector C,
                     FloatVector x0, FloatVector y0) const
{
    return this->validate_feasibility(D,g,C,x0);
}

Bool OptimiserBase::
validate_feasibility(IntervalVector D, IntervalVectorFunction g, IntervalVector C,
                     FloatVector x0) const
{
    ARIADNE_PRECONDITION(D.size()==g.argument_size());
    ARIADNE_PRECONDITION(C.size()==g.result_size());
    ARIADNE_PRECONDITION(x0.size()==D.size());
    ARIADNE_LOG(2,"validate_feasibility\n");
    ARIADNE_LOG(3,"D="<<D<<", g="<<g<<", C="<<C<<"\n");
    ARIADNE_LOG(3,"x0="<<x0<<"\n");

    IntervalVector x(x0);
    ARIADNE_LOG(3,"x="<<x<<"\n");

    IntervalVector gx=g(x);
    ARIADNE_LOG(4,"gx="<<gx<<"\n");

    List<Nat> equalities, inequalities;
    for(Nat i=0; i!=C.size(); ++i) {
        if(C[i].lower()==C[i].upper()) {
            equalities.append(i);
        } else {
            inequalities.append(i);
            if(!subset(gx[i],C[i])) {
                ARIADNE_LOG(3,"g["<<i<<"](x)="<<gx[i]<<", C["<<i<<"]="<<C[i]<<"\n");
                return false; }
        }
    }

    if(equalities.empty()) { ARIADNE_LOG(2,"feasible\n"); return true; }

    Nat k=equalities.size();
    Nat n=D.size();
    IntervalVectorFunction h(equalities.size(),g.argument_size());
    FloatVector c(equalities.size());
    for(uint i=0; i!=equalities.size(); ++i) {
        h[i] = g[equalities[i]];
        c[i] = C[equalities[i]].lower();
    }
    ARIADNE_LOG(5,"h="<<h<<" c="<<c<<" h(x)-c="<<(h(x0)-c)<<"\n");

    // Attempt to solve h(x0+AT*w)=0
    IntervalMatrix AT = transpose(h.jacobian(x0));
    ARIADNE_LOG(5,"A="<<transpose(AT)<<"\n");
    IntervalVector w0(k,Interval(0));

    bool found_solution=false;
    bool validated_solution=false;

    IntervalVector w(k), mw(k), nw(k);
    IntervalVector mx(n);

    for(uint ii=0; ii!=12; ++ii) {
        mw=midpoint(w);
        x=make_exact(x0)+AT*w;
        mx=make_exact(x0)+AT*mw;
        nw = mw - solve(h.jacobian(x)*AT,IntervalVector(h(mx)-make_exact(c)));
        ARIADNE_LOG(7,"w="<<w<<", h(x0+AT*w)="<<h(x)<<", nw="<<nw<<", subset="<<subset(nw,w)<<"\n");

        if(!found_solution) {
            if(subset(nw,w)) {
                found_solution=true;
                w=hull(w0,w);
            } else {
                w=hull(w,IntervalVector(2.0*nw)-w);
            }
        } else {
            if(subset(nw,w)) {
                validated_solution=true;
            } else if(validated_solution) {
                // Not a contraction, so presumably accurate enough
                break;
            }
            w=intersection(nw,w);
        }

    }
    ARIADNE_LOG(5,"w="<<w<<", validated="<<validated_solution<<"\n");

    if(!validated_solution) { return false; }

    // Compute x value
    x=make_exact(x0)+AT*w;
    ARIADNE_LOG(3,"x="<<x<<"\n");
    gx=g(x);
    ARIADNE_LOG(3,"g(x)="<<gx<<"\n");

    // Check that equality constraints are plausible
    ARIADNE_DEBUG_ASSERT(contains(h(x)-make_exact(c),FloatVector(k)));

    // Check inequality constraints once more
    for(uint i=0; i!=C.size(); ++i) {
        if(C[i].lower()==C[i].upper()) {
            ARIADNE_DEBUG_ASSERT(subset(C[i],gx[i]));
        } else {
            if(!subset(gx[i],C[i])) {
                return false;
            }
        }
    }
    return true;
}


Bool OptimiserBase::
validate_infeasibility(IntervalVector D, IntervalVectorFunction g, IntervalVector C,
                       FloatVector x, FloatVector y) const
{
    ARIADNE_PRECONDITION(D.size()==g.argument_size());
    ARIADNE_PRECONDITION(C.size()==g.result_size());
    ARIADNE_PRECONDITION(x.size()==D.size());
    ARIADNE_PRECONDITION(y.size()==C.size());
    ARIADNE_LOG(2,"validate_infeasibility\n");
    // Compute first-order approximation to g(D) centred at x.
    // For feasibilty, have yg(D) cap yC nonempty.
    // Estimate y g(X) = y g(x) + y Dg(X).(X-x)

    // Compute y.C
    Interval yC = dot(IntervalVector(y),C);

    // Compute Taylor estimate of y g(X)
    VectorTaylorFunction tg(D,g,default_sweeper());
    ScalarTaylorFunction tyg(D,default_sweeper());
    for(uint j=0; j!=y.size(); ++j) { tyg += y[j]*tg[j]; }
    Interval tygD = tyg(D);

    IntervalMatrix dgD = g.jacobian(D);
    IntervalVector ydgD = IntervalVector(y) * dgD;

    Interval ygx = dot(IntervalVector(y),g(IntervalVector(x)));

    Interval ygD = ygx;
    for(uint i=0; i!=x.size(); ++i) {
        ygD += ydgD[i] * (D[i]-x[i]);
    }

    ARIADNE_LOG(4,"yC="<<yC<<" tygD="<<tygD<<" ygD="<<ygD<<"\n");

    if(empty(intersection(yC,ygD))) { ARIADNE_LOG(3,"infeasible\n"); return true; }
    else { return false; }
}

// FIXME: Look at this code again, especially relating to generalised Lagrange multipliers
Bool OptimiserBase::
is_infeasibility_certificate(IntervalVector d, IntervalVectorFunction g, IntervalVector c, FloatVector y) const
{
    // Try to prove lambda.(g(y)-c) != 0
    const uint n=c.size();

    ScalarTaylorFunction tyg(d,default_sweeper());
    for(uint i=0; i!=n; ++i) {
        tyg+=Interval(y[i])*ScalarTaylorFunction(d,g[i],default_sweeper());
    }
    Interval iygx = tyg(d);

    Interval iyc = 0;
    for(uint i=0; i!=n; ++i) {
        iyc+=Interval(y[i])*c[i];
    }

    if(disjoint(iyc,iygx)) {
        return true;
    } else {
        return false;
    }
}






IntervalVector OptimiserBase::
minimise(IntervalScalarFunction f, IntervalVector D, IntervalVectorFunction g, IntervalVectorFunction h) const
{
    ARIADNE_LOG(2,"minimise::feasible(f,D,g,h)\n");
    IntervalVectorFunction gh=join(g,h);
    IntervalVector C(gh.result_size(),Interval(0.0));
    for(uint i=0; i!=g.result_size(); ++i) { C[i]=Interval(-inf,0.0); }
    return this->minimise(f,D,gh,C);
}



Tribool OptimiserBase::
feasible(IntervalVector D, IntervalVectorFunction g, IntervalVectorFunction h) const
{
    ARIADNE_LOG(2,"OptimiserBase::feasible(D,g,h)\n");
    IntervalVectorFunction gh=join(g,h);
    IntervalVector C(gh.result_size(),Interval(0.0));
    for(uint i=0; i!=g.result_size(); ++i) { C[i]=Interval(-inf,0.0); }
    return this->feasible(D,gh,C);
}


//------- NonlinearInfeasibleInteriorPointOptimiser -------------------------//

struct NonlinearInfeasibleInteriorPointOptimiser::PrimalDualData {
    FloatVector w,x,y;
};

struct NonlinearInfeasibleInteriorPointOptimiser::StepData : public PrimalDualData {
    FloatVector vl,wl,xl,zl,vu,wu,xu,zu; Float mu;
};

IntervalVector NonlinearInfeasibleInteriorPointOptimiser::
minimise(IntervalScalarFunction f, IntervalVector D, IntervalVectorFunction g, IntervalVector C) const
{
    ARIADNE_LOG(2,"NonlinearInfeasibleInteriorPointOptimiser::minimise(f,D,g,C)\n");
    ARIADNE_LOG(2,"  f="<<f<<", D="<<D<<", g="<<g<<", C="<<C<<"\n");

    static const double VALUE_TOLERANCE=1e-8;
    static const double STATE_TOLERANCE=1e-8;
    static const uint MAXIMUM_STEPS=24;

    ARIADNE_ASSERT(f.argument_size()==D.size());
    ARIADNE_ASSERT(g.argument_size()==D.size());
    ARIADNE_ASSERT(g.result_size()==C.size());
    StepData v;
    FloatVector& x=v.x;
    FloatVector& y=v.y;
    C=intersection(C,g(D)+IntervalVector(C.size(),Interval(-1,+1)));
    this->setup_feasibility(D,g,C,v);
    FloatVector oldx=x;

    static const float MU_MIN = 1e-12;

    // FIXME: Allow more steps
    for(uint i=0; i!=MAXIMUM_STEPS; ++i) {
        ARIADNE_LOG(4,"  f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x)<<"\n");
        oldx=x;
        Float oldfx=f(oldx);
        this->step(f,D,g,C,v);
        if(this->is_infeasibility_certificate(D,g,C,y)) {
            ARIADNE_LOG(2,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x)<<"\n");
            ARIADNE_LOG(2,"infeasible\n");
            std::cerr<<"EXCEPTION: "<<InfeasibleProblemException().what()<<"\n";
            throw InfeasibleProblemException();
        }
        Float fx=f(x);
        if(abs(fx-oldfx)<VALUE_TOLERANCE && norm(oldx-x)<STATE_TOLERANCE) {
            break;
        }
        if(v.mu<MU_MIN) {
            break;
        }
    }
    ARIADNE_LOG(2,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x)<<"\n");
    if(this->validate_feasibility(D,g,C,x)) {
        ARIADNE_LOG(2,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x)<<"\n");
        return x;
    }
    ARIADNE_LOG(2,"indeterminate_feasibility\n");
    throw IndeterminateFeasibilityException();
}

Tribool NonlinearInfeasibleInteriorPointOptimiser::
feasible(IntervalVector D, IntervalVectorFunction g, IntervalVector C) const
{
    ARIADNE_LOG(2,"NonlinearInfeasibleInteriorPointOptimiser::feasible(D,g,C)\n");
    ARIADNE_LOG(3,"D="<<D<<", g="<<g<<", C="<<C<<"\n");

    ARIADNE_ASSERT(g.argument_size()==D.size());
    ARIADNE_ASSERT(g.result_size()==C.size());

    StepData v;
    FloatVector& x=v.x;
    FloatVector& y=v.y;

    FloatScalarFunction f(D.size());
    IntervalVector R=intersection(g(D)+IntervalVector(C.size(),Interval(-1,+1)),C);
    this->setup_feasibility(D,g,R,v);

    static const float MU_MIN = 1e-12;

    // FIXME: Allow more steps
    for(uint i=0; i!=12; ++i) {
        ARIADNE_LOG(5,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x)<<"\n");
        this->step(f,D,g,R,v);
        if(this->validate_feasibility(D,g,C,x)) {
            ARIADNE_LOG(3,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x)<<"\n");
            ARIADNE_LOG(2,"feasible\n");
            return true;
        }
        if(this->is_infeasibility_certificate(D,g,C,y)) {
            ARIADNE_LOG(3,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x)<<"\n");
            ARIADNE_LOG(2,"infeasible\n");
            return false;
        }
        if(v.mu<MU_MIN) {
            break;
        }
    }
    ARIADNE_LOG(3,"f(x)="<<f(x)<<", x="<<x<<", y="<<y<<", g(x)="<<g(x)<<"\n");
    ARIADNE_LOG(2,"indeterminate\n");
    return indeterminate;
}

Void NonlinearInfeasibleInteriorPointOptimiser::
setup_feasibility(const IntervalVector& D, const FloatVectorFunctionInterface& g, const IntervalVector& C,
                  StepData& v) const
{
    Interval I(-1,+1);
    Nat m=C.size(); Nat n=D.size();

    v.x=midpoint(D);
    v.y=FloatVector(m,0.0);
    v.w=midpoint(C);

    //stp.xl=lower(D)-x;
    v.wl=FloatVector(m,-1.0);
    v.wu=FloatVector(m,+1.0);
    v.xl=lower_bounds(D)-v.x;
    v.xu=upper_bounds(D)-v.x;
    v.vl=FloatVector(m,-1.0);
    v.vu=FloatVector(m,+1.0);
    v.zl=FloatVector(n,-1.0);
    v.zu=FloatVector(n,+1.0);
    // FIXME: What should relaxation parameter be?
    v.mu=1.0;
}

Void
NonlinearInfeasibleInteriorPointOptimiser::step(
    const FloatScalarFunctionInterface& f, const IntervalVector& d, const FloatVectorFunctionInterface& g, const IntervalVector& c,
    StepData& v) const
{
    FloatVector& w=v.w; FloatVector& x=v.x; FloatVector& y=v.y; Float& mu=v.mu;
    FloatVector& wl=v.wl; FloatVector& wu=v.wu; FloatVector& xl=v.xl; FloatVector& xu=v.xu;
    FloatVector& vl=v.vl; FloatVector& vu=v.vu; FloatVector& zl=v.zl; FloatVector& zu=v.zu;
    FloatVector cl=lower_bounds(c); FloatVector cu=upper_bounds(c);
    FloatVector dl=lower_bounds(d); FloatVector du=upper_bounds(d);

    ARIADNE_LOG(4,"NonlinearInfeasibleInteriorPointOptimiser::step(f,D,g,C,...)\n");
    ARIADNE_LOG(5,"  f="<<f<<", D="<<d<<", g="<<g<<", C="<<c<<"\n");
    ARIADNE_LOG(5,"  w ="<<w<<",  x ="<<x<<", y ="<<y<<" mu="<<mu<<"\n");
    ARIADNE_LOG(5,"  wl="<<wl<<", wu="<<wu<<", xl="<<xl<<", xu="<<xu<<"\n");
    ARIADNE_LOG(5,"  vl="<<vl<<", vu="<<vu<<", zl="<<zl<<", zu="<<zu<<"\n");
    ARIADNE_LOG(9,"  cl-wl="<<cl-wl<<", dl-xl="<<dl-xl<<"\n");
    ARIADNE_LOG(9,"    w  ="<<w<<",   x  ="<<x<<"\n");
    ARIADNE_LOG(9,"  cu-wu="<<cu-wu<<", du-xu="<<du-xu<<"\n");
    static const double gamma=1.0/1024;
    static const double sigma=1.0/8;
    static const double scale=0.75;

    const uint n=d.size();
    const uint m=c.size();

    ARIADNE_ASSERT_MSG(f.argument_size()==d.size(),"f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_ASSERT_MSG(g.argument_size()==d.size(),"f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_ASSERT_MSG(g.result_size()==c.size(),  "f="<<f<<", D="<<d<<", g="<<g<<", C="<<c);
    ARIADNE_ASSERT(w.size()==m);
    ARIADNE_ASSERT(x.size()==n);
    ARIADNE_ASSERT(y.size()==m);

    mu = mu * sigma;

    FloatDifferential ddfx=f.evaluate(FloatDifferential::variables(2,x));
    ARIADNE_LOG(9,"ddfx="<<ddfx<<"\n");
    Vector<FloatDifferential> ddgx=g.evaluate(FloatDifferential::variables(2,x));
    ARIADNE_LOG(9,"ddgx="<<ddgx<<"\n");

    Float fx = ddfx.value();
    Vector<Float> gx = ddgx.value();
    ARIADNE_LOG(7,"f(x)="<<fx<<"\n");
    ARIADNE_LOG(7,"g(x)="<<gx<<"\n");
    Vector<Float> Jfx = ddfx.gradient();
    Matrix<Float> A = ddgx.jacobian();
    Matrix<Float>& Jgx = A;
    ARIADNE_LOG(7,"Df(x)="<<Jfx<<"\n");
    ARIADNE_LOG(7,"Dg(x)="<<Jgx<<"\n");

    // H is the Hessian matrix H of the Lagrangian $L(x,\lambda) = f(x) + \sum_k g_k(x) $
    Matrix<Float> YH = ddfx.hessian();
    for(uint i=0; i!=m; ++i) {
        YH+=y[i]*ddgx[i].hessian();
    }
    ARIADNE_LOG(7,"D2f(x)="<<ddfx.hessian()<<"\n");
    ARIADNE_LOG(7,"D2f(x)+Y.D2g(x)="<<YH<<"\n");

    // Set up the system of equations
    // (A^TDA + E - Y.H) dx = A^T(r_w-Dr_y)+r_x
    // dw = A \delta x + r_y
    // dy = r_w - D dw

    FloatDiagonalMatrix const& Vl=diagonal_matrix(vl);
    FloatDiagonalMatrix const& Vu=diagonal_matrix(vu);
    FloatDiagonalMatrix const& Wl=diagonal_matrix(wl);
    FloatDiagonalMatrix const& Wu=diagonal_matrix(wu);
    FloatDiagonalMatrix const& Xl=diagonal_matrix(xl);
    FloatDiagonalMatrix const& Xu=diagonal_matrix(xu);
    FloatDiagonalMatrix const& Zl=diagonal_matrix(zl);
    FloatDiagonalMatrix const& Zu=diagonal_matrix(zu);

    // Compute the diagonal matrices
    //   D=XL/ZL+XU/ZU  E=WL/VL+WU/VU
    FloatDiagonalMatrix Dl=Vl/Wl;
    FloatDiagonalMatrix Du=Vu/Wu;
    FloatDiagonalMatrix D=Dl+Du;
    ARIADNE_LOG(9,"D="<<D<<"\n");
    FloatDiagonalMatrix El=Zl/Xl;
    FloatDiagonalMatrix Eu=Zu/Xu;
    FloatDiagonalMatrix E=El+Eu;
    ARIADNE_LOG(9,"E="<<E<<"\n");

    // normal equation matrix
    FloatMatrix S=YH;
    atda(S,A,D);
    S+=E;

    //FloatMatrix EE(n,n); for(uint j=0; j!=n; ++j) { EE[j][j]=E[j]; }
    //FloatMatrix DD(m,m); for(uint i=0; i!=m; ++i) { DD[i][i]=E[i]; }

    ARIADNE_LOG(9,"S="<<S<<"\n");
    ARIADNE_DEBUG_ASSERT(norm(FloatMatrix(S-(YH+E+transpose(A)*(D*A))))/norm(S)<1e-8);
    FloatMatrix Sinv=inverse(S);
    ARIADNE_LOG(9,"Sinv="<<Sinv<<"\n");

    // Construct the residuals
    // The residual for the slack variable xl is given by the duality condition xl.zl=mu as mu/xl-zl
    // The residual for the dual variable zl is given by the slackness condition x-xl-cl
    // The residual for the auxiliary variable w is given by y-(vu-vl)
    // The residual for the dual variable y is given by g(x)-w
    FloatVector ew=(vl+vu)-y;
    FloatVector ex=Jfx+y*Jgx+(zl+zu);
    FloatVector ey=gx-w;
    FloatVector ewl=esub(vl,ediv(mu,wl));
    FloatVector ewu=esub(vu,ediv(mu,wu));
    FloatVector exl=esub(zl,ediv(mu,xl));
    FloatVector exu=esub(zu,ediv(mu,xu));
    FloatVector evl=w+wl-cl;
    FloatVector evu=w+wu-cu;
    FloatVector ezl=x+xl-dl;
    FloatVector ezu=x+xu-du;

    ARIADNE_LOG(9,"ew="<<ew<<", ex="<<ex<<", ey="<<ey<<"\n");
    ARIADNE_LOG(9,"ewl="<<ewl<<", ewu="<<ewu<<", exl="<<exl<<" exu="<<exu<<"\n");
    ARIADNE_LOG(9,"evl="<<evl<<", evu="<<evu<<", ezl="<<ezl<<" ezu="<<ezu<<"\n");

    FloatVector rw = ew - (ewl+ewu) + Dl*evl + Du*evu;
    FloatVector rx = ex - (exl+exu) + El*ezl + Eu*ezu;
    FloatVector& ry = ey;

    // Solve linear system
    // ( D   0  -I ) (dw)   (rw)
    // ( 0  H+E A^T) (dx) = (rx)
    // (-I   A   0 ) (dy) = (ry)

    FloatVector r = (rw+D*ry)*A+rx;
    ARIADNE_LOG(9,"rw="<<rw<<" rx="<<rx<<" ry="<<ry<<"\n");
    ARIADNE_LOG(9,"r="<<r<<"\n");

    // Compute the differences
    FloatVector dx = solve(S,r);
    ARIADNE_LOG(9,"S*dx="<<S*dx<<" r="<<r<<"\n");
    ARIADNE_LOG(9,"S*inverse(S)-I="<<S*inverse(S)-FloatMatrix::identity(n)<<"\n");
    ARIADNE_DEBUG_ASSERT(norm(S*dx - r)/max(1.0,norm(r))<1e-4);

    FloatVector dw = A*dx-ry;
    FloatVector dy = D*dw-rw;
    ARIADNE_LOG(9,"dw="<<dw<<" dx="<<dx<<" dy="<<dy<<"\n");

    ARIADNE_LOG(9,"YH*dx+E*dx+dy*A="<<(YH*dx+E*dx+dy*A)<<", rx="<<rx<<"\n");

    // Check solution of linear system for residuals
    ARIADNE_DEBUG_ASSERT(norm(D*dw-dy-rw)/max(1.0,norm(rw))<1e-4);
    ARIADNE_DEBUG_ASSERT(norm(YH*dx+E*dx+dy*A-rx)/max(1.0,norm(rx))<1e-2);
    ARIADNE_DEBUG_ASSERT(norm(-dw+A*dx-ry)/max(1.0,norm(ry))<1e-4);

    FloatVector dwl = evl-dw;
    FloatVector dwu = evu-dw;
    FloatVector dxl = ezl-dx;
    FloatVector dxu = ezu-dx;
    FloatVector dvl = ewl-Dl*dwl;
    FloatVector dvu = ewu-Du*dwu;
    FloatVector dzl = exl-El*dxl;
    FloatVector dzu = exu-Eu*dxu;

    ARIADNE_LOG(9,"dwl="<<dwl<<", dwu="<<dwu<<", dxl="<<dxl<<" dxu="<<dxu<<"\n");
    ARIADNE_LOG(9,"dvl="<<dvl<<", dvu="<<dvu<<", dzl="<<dzl<<" dzu="<<dzu<<"\n");

    ARIADNE_LOG(9,"YH*dx+dy*A+dzl+dzu="<<(YH*dx+dy*A+dzl+dzu)<<", ex="<<ex<<"\n");
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

    FloatVector nw; FloatVector nx; FloatVector ny;
    FloatVector nwl; FloatVector nwu; FloatVector nxl; FloatVector nxu;
    FloatVector nvl; FloatVector nvu; FloatVector nzl; FloatVector nzu;


    Float alpha=1.0;
    nx = x-alpha*dx;
    // Pick an update value which minimises the objective function
    Float fxmin=f(nx);
    Float alphamin=1.0;
    static const uint REDUCTION_STEPS=4;
    for(uint i=0; i!=REDUCTION_STEPS; ++i) {
        alpha*=scale;
        nx = x-alpha*dx;
        Float fnx=f(nx);
        if(fnx<fxmin*scale) {
            fxmin=fnx;
            alphamin=alpha;
        }
    }
    //alpha=alphamin;
    alpha=1.0;

    // Since we need to keep the point feasible, but the updates are linear
    // we need to validate feasibility directly.
    static const double MINIMUM_ALPHA=1e-16;
    bool allfeasible=false;
    while(alpha>MINIMUM_ALPHA && !allfeasible) {
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
    if(alpha<=MINIMUM_ALPHA) {
        ARIADNE_LOG(9," w="<<w<<"  x="<<x<<"  y="<<y<<"\n");
        ARIADNE_LOG(9," nw="<<nw<<"  nx="<<nx<<"  ny="<<ny<<"\n");
        throw NearBoundaryOfFeasibleDomainException(); }
    ARIADNE_LOG(5,"alpha="<<alpha<<"\n");
    ARIADNE_LOG(9,"nw="<<nw<<" nx="<<nx<<" ny="<<ny<<"\n");
    ARIADNE_LOG(9,"nwl="<<nwl<<", nwu="<<nwu<<", nxl="<<nxl<<" nxu="<<nxu<<"\n");
    ARIADNE_LOG(9,"nvl="<<nvl<<", nvu="<<nvu<<", nzl="<<nzl<<" nzu="<<nzu<<"\n");

    w=nw; x=nx; y=ny;
    wl=nwl; wu=nwu; xl=nxl; xu=nxu;
    vl=nvl; vu=nvu; zl=nzl; zu=nzu;

    Float nmu = 0.0;
    for(uint i=0; i!=m; ++i) {
        nmu = nmu + wl[i]*vl[i] + wu[i]*vu[i];
    }
    for(uint j=0; j!=n; ++j) {
        nmu = nmu + xl[j]*zl[j] + xu[j]*zu[j];
    }
    nmu /= (2*(m+n));
    mu = nmu;

    ARIADNE_LOG(9,"nmu="<<nmu<<"\n");

}





//------- NonlinearInteriorPointOptimiser -----------------------------------//

IntervalVector NonlinearInteriorPointOptimiser::
minimise(IntervalScalarFunction f, IntervalVector D, IntervalVectorFunction g, IntervalVector C) const
{
    ARIADNE_LOG(2,"NonlinearInteriorPointOptimiser::minimise(f,D,g,C)\n");
    ARIADNE_LOG(3,"f="<<f<<" D="<<D<<" g="<<g<<" C="<<C<<"\n");
    IntervalVectorFunction h(0,D.size());

    IntervalVector gD = g(D);
    if(disjoint(gD,C)) { throw InfeasibleProblemException(); }

    FloatVector x = midpoint(D);
    FloatVector w = midpoint(intersection(gD,C));

    FloatVector kappa(g.result_size(),0.0);
    FloatVector lambda(h.result_size(),0.0);
    Float mu = 1.0;


    for(uint i=0; i!=12; ++i) {
        this->minimisation_step(f,D,g,C,h, x,w, kappa,lambda, mu);
        if(i%3==0 && i<=10) { mu *= 0.25; }
    }

    return IntervalVector(x);
}



// See Hande Y. Benson, David F. Shanno, And Robert J. Vanderbei,
// "Interior-point methods for nonconvex nonlinear programming: Jamming and comparative numerical testing"
// For some of the terminology used


// min f(x) | x\in D & w\in C | g(x) = w & h(x) = 0
// Lagrange multipliers kappa d(g(x)-w); lambda dh(x)
Void NonlinearInteriorPointOptimiser::
minimisation_step(const FloatScalarFunction& f, const IntervalVector& d, const FloatVectorFunction& g, const IntervalVector& c, const FloatVectorFunction& h,
                  FloatVector& x, FloatVector& w, FloatVector& kappa, FloatVector& lambda, const Float& mu) const
{
    const uint n=x.size();
    const uint m=kappa.size();
    const uint l=lambda.size();

    ARIADNE_DEBUG_PRECONDITION(w.size()==kappa.size());
    ARIADNE_DEBUG_PRECONDITION(f.argument_size()==n);
    ARIADNE_DEBUG_PRECONDITION(g.argument_size()==n);
    ARIADNE_DEBUG_PRECONDITION(h.argument_size()==n);
    ARIADNE_DEBUG_PRECONDITION(g.result_size()==m);
    ARIADNE_DEBUG_PRECONDITION(h.result_size()==l);
    ARIADNE_DEBUG_PRECONDITION(contains(d,x));
    ARIADNE_DEBUG_PRECONDITION(contains(c,w));
    ARIADNE_DEBUG_PRECONDITION(mu>0);

    ARIADNE_LOG(4,"NonlinearInteriorPointOptimiser::minimisation_step(f,D,g,C,h, x,w, kappa,lambda, mu)\n");
    ARIADNE_LOG(5,"x="<<x<<"\n");
    ARIADNE_LOG(7,"w="<<w<<"\n");
    ARIADNE_LOG(7,"kappa="<<kappa<<"\n");
    ARIADNE_LOG(7,"lambda="<<lambda<<"\n");
    ARIADNE_LOG(7,"mu="<<mu<<"\n");

    FloatVector slack(2*n);
    FloatVectorRange slackl(slack,range(0,n));
    FloatVectorRange slacku(slack,range(n,2*n));

    FloatDifferential ddfx=f.evaluate(FloatDifferential::variables(2,x));
    Vector<FloatDifferential> ddgx=g.evaluate(FloatDifferential::variables(2,x));
    Vector<FloatDifferential> ddhx=h.evaluate(FloatDifferential::variables(2,x));

    // G is the constraint value vector
    Float fx = ddfx.value();
    FloatVector gx = ddgx.value();
    FloatVector hx = ddhx.value();
    ARIADNE_LOG(5,"f(x)="<<fx<<"\n");
    ARIADNE_LOG(5,"g(x)="<<gx<<"\n");
    ARIADNE_LOG(5,"h(x)="<<hx<<"\n");
    ARIADNE_LOG(9,"g(x)-w="<<(gx-w)<<"\n");

    // A, B are the derivative matrices aij=dgi/dxj
    // HACK: Need to explicitly set size of Jacobian if g or h have result_size of zero
    FloatVector df = ddfx.gradient();
    ARIADNE_LOG(9,"df(x)="<<df<<"\n");
    FloatMatrix A = ddgx.jacobian();
    if(m==0) { A=FloatMatrix(m,n); }
    ARIADNE_LOG(9,"A="<<A<<"\n");
    FloatMatrix B = ddhx.jacobian();
    if(l==0) { B=FloatMatrix(l,n); }
    ARIADNE_LOG(9,"B="<<B<<"\n");



    // H is the Hessian matrix H[i1,i2] = df/dx[i1]dx[i2] + Sum_[j]kappa[j]*dg[j]/dx[i1]dx[i2] + Sum[k]lambda[k]*dh[k]/dx[i1]dx[i2]
    FloatMatrix H = ddfx.hessian();
    for(uint j=0; j!=m; ++j) { H += kappa[j] * ddgx[j].hessian(); }
    for(uint k=0; k!=l; ++k) { H += lambda[k] * ddhx[k].hessian(); }
    ARIADNE_LOG(9,"H="<<H<<"\n");

    // Determines the weighting to give to the relaxation parameter mu
    // for equality constraints relative to other constraints
    static const double EQUALITY_RELAXATION_MULTIPLIER = 1.0;

    // Compute the residuals and contributions from slack in x and w
    //   rx = df/dx[i] + Sum[j] dg[j]/dx[i] * kappa[j] + Sum[k] dh[k]/dx[i] * lambda[j] + mu *( 1/(xu[i]-x[i]) - 1/(x[i]-xl[i]) )
    FloatVector rx = df + kappa * A + lambda * B;
    FloatDiagonalMatrix D(n);
    for(uint i=0; i!=n; ++i) {
        Float nuu = rec(d[i].upper()-x[i]);
        Float nul = rec(x[i]-d[i].lower());
        rx[i] += mu * ( nuu - nul );
        D[i] = mu * ( nuu*nuu + nul*nul );
    }

    //   rw = - kappa[j] + mu *( 1/(wu[i]-w[i]) - 1/(w[i]-wl[i]) )
    FloatVector rw = -kappa;
    FloatDiagonalMatrix C(m);
    for(uint j=0; j!=m; ++j) {
        Float nuu = rec(c[j].upper()-w[j]);
        Float nul = rec(w[j]-c[j].lower());
        rw[j] += (mu*EQUALITY_RELAXATION_MULTIPLIER) * ( nuu - nul );
        C[j] = (mu*EQUALITY_RELAXATION_MULTIPLIER) * ( nuu*nuu + nul*nul );
    }

    //   rkappa = g(x) - w
    FloatVector rkappa = gx - w;

    //   rlambda = h(x)
    FloatVector const& rlambda = hx;

    ARIADNE_LOG(9,"rx="<<rx<<"\n");
    ARIADNE_LOG(9,"rw="<<rw<<"\n");
    ARIADNE_LOG(9,"rkappa="<<rkappa<<"\n");
    ARIADNE_LOG(9,"rlambda="<<rlambda<<"\n");

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
    FloatMatrix& S=H;
    S+=D;
    S+=FloatMatrix(transpose(A))*C*A;
    ARIADNE_LOG(9,"S="<<S<<"\n");

    FloatMatrix Sinv=inverse(S);
    ARIADNE_LOG(9,"R=Sinv="<<Sinv<<"\n");

    FloatMatrix BSinvBT = (B*Sinv)*transpose(B);
    ARIADNE_LOG(9,"B*inverse(S)*BT="<<BSinvBT<<"\n");
    ARIADNE_LOG(9,"inverse(B*inverse(S)*BT)="<<inverse(BSinvBT)<<"\n");

    FloatVector rr = Sinv * (rx + (rkappa * C + rw) * A);
    FloatVector dlambda = inverse(BSinvBT) * (B * rr - rlambda);
    FloatVector dx = rr - dlambda * (B*Sinv);
    FloatVector dw = A * dx - rkappa;
    FloatVector dkappa = rw - C * dw;

    static const Float ALPHA_SCALE_FACTOR = 0.75;
    static const Float MINIMUM_ALPHA = 1e-16;

    // Compute distance to move variables preserving feasibility
    // FIXME: Current implementation might fail due to getting too close to boundary!
    FloatVector newx(n);
    FloatVector neww(m);
    Float alpha = 1.0;
    bool success = false;
    do {
        newx = x - alpha * dx;
        neww = w - alpha * dw;
        if (contains(d,newx) && contains(c,neww)) { success = true; }
        else { alpha *= ALPHA_SCALE_FACTOR; }
        if(alpha<MINIMUM_ALPHA) { throw NearBoundaryOfFeasibleDomainException(); }
    } while(!success);
    ARIADNE_LOG(9,"alpha="<<alpha<<"\n");

    FloatVector newlambda = lambda - alpha * dlambda;
    FloatVector newkappa = kappa - alpha * dkappa;

    ARIADNE_LOG(9,"newx="<<newx<<"\n");
    ARIADNE_LOG(9,"neww="<<neww<<"\n");
    ARIADNE_LOG(9,"newkappa="<<newkappa<<"\n");
    ARIADNE_LOG(9,"newlambda="<<newlambda<<"\n");

    x=newx; w=neww; kappa=newkappa; lambda=newlambda;

    if(verbosity>=6) { std::clog << "\n"; }
}



Tribool NonlinearInteriorPointOptimiser::
feasible(IntervalVector d, IntervalVectorFunction g, IntervalVector c) const
{
    ARIADNE_LOG(2,"NonlinearInteriorPointOptimiser::feasible(D,g,C,h)\n");
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
    static const double inf = std::numeric_limits<double>::infinity();

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

    Vector<FloatDifferential> ddgx=g.evaluate(FloatDifferential::variables(2,x));
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
        } else if(c[j].upper()==+inf) {
        } else if(c[j].lower()==-inf) {
        } else {
            ARIADNE_DEBUG_ASSERT(-inf<c[j].lower() && c[j].lower()<c[j].upper() && c[j].upper()<+inf);
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
    FloatVector rr=ryt + AE*ediv(FloatVector(rx-emul(x,rz)),z) - ryt;


    // Compute the differences
    FloatVector dyt=Sinv*rr;
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

    Vector<FloatDifferential> dg=g.evaluate(FloatDifferential::variables(1,y));
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



Float NonlinearInteriorPointOptimiser::
compute_mu(const IntervalVector& D, const FloatVectorFunction& g, const IntervalVector& C,
           const FloatVector& x, const FloatVector& lambda) const
{
    // Compute the relaxation parameter mu as the average of the product of the Lyapunov exponents and constraint satisfactions
    Float mu = 0.0;
    FloatVector gx = g(x);

    for(uint i=0; i!=C.size(); ++i) {
        if(C[i].lower()==C[i].upper()) { }
        else if(C[i].lower()==-inf) { mu += lambda[i] * (gx[i] - C[i].upper()); }
        else if(C[i].upper()==+inf) { mu += lambda[i] * (gx[i] - C[i].lower()); }
        else { // std::cerr<<"FIXME: Compute mu for bounded constraint\n";
            if (lambda[i] <=0.0) { mu += lambda[i] * (gx[i] - C[i].upper()); }
            else { mu += lambda[i] * (gx[i] - C[i].lower()); }
        }
    }
    mu /= C.size();
    return mu;
}


void NonlinearInteriorPointOptimiser::
setup_feasibility(const IntervalVector& d, const FloatVectorFunction& g, const IntervalVector& c,
                  FloatVector& x, FloatVector& y) const
{
    const uint l=2*(d.size()+c.size());
    y=midpoint(d);
    x=FloatVector(l,1.0/l);
    //compute_tz(d,g,c,y,t,z);
}




//------- PenaltyFunctionOptimiser ------------------------------------------//

PenaltyFunctionOptimiser* PenaltyFunctionOptimiser::
clone() const
{
    return new PenaltyFunctionOptimiser(*this);
}

IntervalVector PenaltyFunctionOptimiser::
minimise(IntervalScalarFunction f, IntervalVector D, IntervalVectorFunction g, IntervalVector C) const
{
    ARIADNE_NOT_IMPLEMENTED;
}

Tribool PenaltyFunctionOptimiser::
feasible(IntervalVector D, IntervalVectorFunction g, IntervalVector C) const
{
    ARIADNE_LOG(2,"PenaltyFunctionOptimiser::feasible(D,g,C)\n");
    ARIADNE_LOG(3,"D="<<D<<" g="<<g<<" C="<<C<<" \n");

    FloatVector x=midpoint(D);

    FloatVector w=midpoint(C);
    for(uint i=0; i!=C.size(); ++i) {
        if(C[i].upper()==+inf) { w[i]=C[i].lower()+1.0; }
        else if(C[i].lower()==-inf) { w[i]=C[i].upper()-1.0; }
    }

    FloatVector y(C.size(),0.0);

    ARIADNE_LOG(5,"x="<<x<<" w="<<w<<" y="<<y<<"\n");

    for(uint i=0; i!=10; ++i) {
        this->feasibility_step(D,g,C,x,y,w);
    }
    return this->check_feasibility(D,g,C,x,y);
}

Void PenaltyFunctionOptimiser::
feasibility_step(const IntervalVector& X, const FloatVectorFunction& g, const IntervalVector& W,
                 FloatVector& x, FloatVector& w, Float& mu) const
{
    FloatVectorFunction h(0u,X.size());
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
feasibility_step(const IntervalVector& X, const IntervalVectorFunction& g, const IntervalVector& W,
                 IntervalVector& x, IntervalVector& w) const
{
    ARIADNE_NOT_IMPLEMENTED;
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

    t=+inf;
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






Tribool ApproximateOptimiser::
feasible(IntervalVector D, IntervalVectorFunction h) const
{
    ARIADNE_LOG(2,"ApproximateOptimiser::feasible(D,h)\n");
    ARIADNE_LOG(3,"D="<<D<<", h="<<h<<"\n");
    FloatVector x=midpoint(D);
    FloatVector y(h.result_size(),0.0);

    for(uint i=0; i!=8; ++i) {
        this->feasibility_step(D,h,x,y);
    }

    if(norm(h(x))<1e-10) { return true; }

    if(!contains(dot(IntervalVector(y),h(D)),0.0)) { return false; }

    return indeterminate;
}

Void ApproximateOptimiser::
feasibility_step(const IntervalVector& D, const FloatVectorFunction& h,
                 FloatVector& x, FloatVector& y) const
{
    ARIADNE_LOG(4,"ApproximateOptimiser::feasibility_step(D,h,x,y)\n");
    ARIADNE_LOG(5,"x="<<x<<" y="<<y<<"\n");
    static const double SCALE_FACTOR = 0.75;
    const uint n=x.size();
    const uint m=y.size();
    // Solve equations y Dh(x) - 1/(x-xl) + 1/(xu-x) = 0; h(x) = 0
    Vector<FloatDifferential> ddhx=h.evaluate(FloatDifferential::variables(2,x));
    FloatMatrix A = ddhx.jacobian();
    ARIADNE_LOG(6,"A="<<A<<" b="<<ddhx.value()<<"\n");

    FloatMatrix H(n,n);
    for(uint i=0; i!=m; ++i) { H += y[i] * ddhx[i].hessian(); }
    for(uint j=0; j!=n; ++j) {
        H[j][j] += rec(sqr(x[j]-D[j].lower()));
        H[j][j] += rec(sqr(D[j].upper()-x[j]));
    }

    FloatVector rx = y * A;
    for(uint j=0; j!=n; ++j) {
        rx[j] -= rec(x[j]-D[j].lower());
        rx[j] += rec(D[j].upper()-x[j]);
    }
    FloatVector ry = ddhx.value();
    ARIADNE_LOG(5,"rx="<<rx<<" ry="<<ry<<"\n");

    // S = A Hinv AT
    // H dx + AT dy = rx; A dx = ry;
    //  dx = Hinv ( rx - AT dy )
    //  dy = Sinv ( A Hinv rx - ry )
    FloatMatrix Hinv=inverse(H);
    ARIADNE_LOG(6,"H="<<H<<" Hinv="<<Hinv<<"\n");
    FloatMatrix S=A*Hinv*transpose(A);
    FloatMatrix Sinv=inverse(S);
    ARIADNE_LOG(6,"S="<<S<<" Sinv="<<Sinv<<"\n");
    FloatVector dy = Sinv * ( A*(Hinv*rx) - ry );
    FloatVector dx = Hinv * ( rx - dy * A);
    ARIADNE_LOG(5,"dx="<<dx<<" dy="<<dy<<"\n");

    Float ax = 1.0;
    FloatVector nx = x-ax*dx;
    while(!contains(D,nx)) {
        ax*=SCALE_FACTOR;
        nx = x - ax * dx;
    }
    FloatVector ny = y-ax*dy;
    ARIADNE_LOG(5,"nx="<<nx<<" ax="<<ax<<" ny="<<ny<<"\n");
    ARIADNE_LOG(6,"h(x)="<<h(nx)<<"\n");

    x=nx; y=ny;
}


Tribool PenaltyFunctionOptimiser::
check_feasibility(IntervalVector D, IntervalVectorFunction g, IntervalVector C,
                     FloatVector fltx, FloatVector flty) const
{
    ARIADNE_PRECONDITION(D.size()==g.argument_size());
    ARIADNE_PRECONDITION(C.size()==g.result_size());
    ARIADNE_PRECONDITION(fltx.size()==D.size());
    ARIADNE_PRECONDITION(flty.size()==C.size());
    ARIADNE_LOG(2,"check_feasibility\n");
    ARIADNE_LOG(3,"D="<<D<<" C="<<C<<"\n");

    IntervalVector x(fltx);
    IntervalVector y(flty);
    IntervalVector gx=g(x);
    ARIADNE_LOG(3,"x="<<x<<" y="<<y<<" g(x)="<<gx<<"\n");

    tribool result = true;

    List<uint> equalities;
    for(uint j=0; j!=C.size(); ++j) {
        if(gx[j].upper()<C[j].lower() || gx[j].lower()>C[j].upper()) {
            return false;
        }
        if(C[j].lower()==C[j].upper()) {
            equalities.append(j);
        } else {
            if(!subset(gx[j],C[j])) { result = indeterminate; }
        }
    }

    if(definitely(result)) {
        if(equalities.empty()) { ARIADNE_LOG(2,"feasible\n"); return true; }

        IntervalVectorFunction h(equalities.size(),g.argument_size());
        FloatVector c(equalities.size());
        for(uint i=0; i!=equalities.size(); ++i) {
            h[i] = g[equalities[i]];
            c[i] = C[equalities[i]].lower();
        }
        ARIADNE_LOG(5,"g="<<g<<"\n");
        ARIADNE_LOG(5,"h="<<h<<" c="<<c<<" h(x)-c="<<FloatVector(h(fltx)-c)<<"\n");

        IntervalVector W(h.result_size(),Interval(-1e-8,1e-8));
        IntervalMatrix AT = transpose(midpoint(h.jacobian(fltx)));
        IntervalVector B = x+AT*W;
        IntervalMatrix IA = h.jacobian(B);
        ARIADNE_LOG(5,"AT="<<AT<<" IA="<<IA<<"\n");
        ARIADNE_LOG(5,"B="<<B<<"\n");

        // Perform an interval Newton step to try to attain feasibility
        IntervalVector nW = inverse(IA*AT) * IntervalVector(h(x)-make_exact(c));
        ARIADNE_LOG(4,"W="<<W<<"\nnew_W="<<nW<<"\n");
        if(subset(B,D) && subset(nW,W)) { ARIADNE_LOG(3,"feasible\n"); return true; }
        else { result=indeterminate; }
    }

    // Compute first-order approximation to g(D) centred at x.
    // For feasibilty, have yg(D) cap yC nonempty.
    // Estimate y g(X) = y g(x) + y Dg(X).(X-x)

    // Compute y.C
    Interval yC = dot(y,C);

    // Compute Taylor estimate of y g(X)
    VectorTaylorFunction tg(D,g,default_sweeper());
    ScalarTaylorFunction tyg(D,default_sweeper());
    for(uint j=0; j!=y.size(); ++j) { tyg += y[j]*tg[j]; }
    Interval tygD = tyg(D);

    IntervalMatrix dgD = g.jacobian(D);
    IntervalVector ydgD = y * dgD;

    Interval ygx = dot(y,gx);

    Interval ygD = ygx;
    for(uint i=0; i!=x.size(); ++i) {
        ygD += ydgD[i] * (D[i]-x[i]);
    }

    ARIADNE_LOG(4,"yC="<<yC<<" tygD="<<tygD<<" ygD="<<ygD<<"\n");

    if(empty(intersection(yC,ygD))) { ARIADNE_LOG(3,"infeasible\n"); return false; }
    else { return indeterminate; }
}







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
feasibility_step(const IntervalVector& D, const FloatVectorFunction& g, const IntervalVector& C,
                 FloatVector& x, FloatVector& y, FloatVector& z) const
{
    ARIADNE_LOG(2,"feasibility_step\n");
    FloatVector xl=lower_bounds(D); FloatVector xu=upper_bounds(D);
    FloatVector zl=lower_bounds(C); FloatVector zu=upper_bounds(C);

    const uint n=x.size();
    const uint m=y.size();

    ARIADNE_LOG(4,"x="<<x<<" y="<<y<<" z="<<z<<"\n");
    Vector<FloatDifferential> ddx = FloatDifferential::variables(2,x);
    Vector<FloatDifferential> ddgx = g.evaluate(ddx);

    FloatMatrix A=ddgx.jacobian();
    ARIADNE_LOG(6,"A="<<A<<"\n");
    FloatVector v = join(join(x,z),y);

    FloatVector r(n+2*m,n+2*m);
    project(r,range(0,n)) = y * A;
    for(uint i=0; i!=n; ++i) {
        r[i] += ( rec(x[i]-xl[i]) - rec(xu[i]-x[i]) );
    }
    for(uint j=0; j!=m; ++j) {
        if(zl[j]==zu[j]) { assert(zu[j]==zl[j]); r[n+j] = z[j]-zl[j]; }
        else { r[n+j] = ( rec(z[j]-zl[j]) - rec(zu[j]-z[j]) - y[j] ); }
    }
    project(r,range(n+m,n+2*m)) = ddgx.value() - z;
    r[n+2*m]=0.0;
    ARIADNE_LOG(5,"r="<<r<<"\n");

    FloatMatrix S(n+2*m+1,n+2*m+1);
    for(uint j=0; j!=m; ++j) {
        FloatMatrix H=ddgx[j].hessian();
        for(uint i1=0; i1!=n; ++i1) {
            for(uint i2=0; i2!=n; ++i2) {
                S[i1][i2]+=y[j]*H[i1][i2];
            }
        }
    }
    for(uint j=0; j!=m; ++j) {
        for(uint i=0; i!=n; ++i) {
            S[i][j+m+n]=A[j][i];
            S[j+m+n][i]=A[j][i];
        }
    }
    for(uint j=0; j!=m; ++j) {
        S[n+j][n+m+j] = -1.0;
        S[n+m+j][n+j] = -1.0;
        //if(zl[j]==zu[j]) { S[n+j][n+j] = -1.0; S[n+j][n+m+j] = 0.0; }
        if(zl[j]==zu[j]) { S[n+j][n+j] = +inf; }
        else { S[n+j][n+j] = - rec(sqr(z[j]-zl[j])) - rec(sqr(zu[j]-z[j])); }
    }
    for(uint i=0; i!=n; ++i) {
        S[i][i]-= rec(sqr(xu[i]-x[i]));
        S[i][i]-= rec(sqr(x[i]-xl[i]));
    }

    for(uint i=0; i!=n; ++i) {
        S[i][n+2*m] -= rec(xu[i]-x[i]);
        S[i][n+2*m] += rec(x[i]-xl[i]);
    }

    for(uint j=0; j!=n; ++j) {
        //S[n+m+j][n+m+j] = -1.0/1024;
    }

    ARIADNE_LOG(5,"S="<<S<<"\n");
    //ARIADNE_LOG(5,"S="<<std::fixed<<pretty(S)<<"\n");

    FloatMatrix Sinv = inverse(S);
    //ARIADNE_LOG(9,"Sinv="<<Sinv<<"\n);
    //ARIADNE_LOG(5,"Sinv="<<std::fixed<<pretty(Sinv)<<"\n");

    FloatVector dv = Sinv * r;
    ARIADNE_LOG(5,"dv="<<dv<<"\n");

    Float alpha = 1.0;
    FloatVector nv = v-dv;
    while(!contains(D,FloatVector(project(nv,range(0,n)))) || !contains(C,FloatVector(project(nv,range(n,n+m)))) ) {
        alpha *= 0.75;
        nv = v-alpha*dv;
    }

    ARIADNE_LOG(4,"nv="<<nv<<" a="<<alpha<<"\n");

    x=project(nv,range(0,n));
    z=project(nv,range(n,n+m));
    y=project(nv,range(n+m,n+2*m));
    ARIADNE_LOG(4,"g(x)-z="<<g(x)-z<<"\n");

}


// Solve equations y Dh(x) - 1/(x-xl) + 1/(xu-x) = 0; h(x) = 0
Tribool IntervalOptimiser::
feasible(IntervalVector D, IntervalVectorFunction h) const
{
    ARIADNE_LOG(2,"IntervalOptimiser::feasible(D,h)\n");
    ARIADNE_LOG(3,"D="<<D<<", h="<<h<<"\n");

    const uint n=D.size();

    IntervalVector zl(n), zu(n);
    FloatVector xl = Ariadne::lower_bounds(D);
    FloatVector xu = Ariadne::upper_bounds(D);

    IntervalVector x=D;
    IntervalVector y(h.result_size(),Interval(-1,+1));
    Interval mu(0,1);

    for(uint i=0; i!=8; ++i) {
        this->feasibility_step(xl,xu,h,x,y,zl,zu,mu);
    }

    return indeterminate;
}

Void IntervalOptimiser::
feasibility_step(const FloatVector& xl, const FloatVector& xu, const IntervalVectorFunction& h,
                 IntervalVector& x, IntervalVector& y, IntervalVector& zl, IntervalVector zu, Interval& mu) const
{
    ARIADNE_LOG(4,"IntervalOptimiser::feasibility_step(D,h,X,Lambda)\n");
    ARIADNE_LOG(5,"[x]="<<x<<" [lambda]="<<y<<", [zl]="<<zl<<", [zu]="<<zu<<" [mu]="<<mu<<"\n");

    const uint n=x.size();
    const uint m=y.size();

    IntervalVector mx=midpoint(x);
    IntervalVector my=midpoint(y);
    IntervalVector mzl=midpoint(zl);
    IntervalVector mzu=midpoint(zu);
    Interval mmu(midpoint(mu));
    ARIADNE_LOG(6,"x~"<<x<<" lambda~="<<y<<", mu~"<<mu<<"\n");

    // Solve equations y Dh(x) - zl + zu = 0; h(x) = 0; (x-xl).zl - mu = 0;  (xu-x).zu - mu = 0; Sum_j y_j^2 - mu = 0
    Vector<IntervalDifferential> ddhx=h.evaluate(IntervalDifferential::variables(2,x));
    Vector<IntervalDifferential> dhmx=h.evaluate(IntervalDifferential::variables(1,mx));
    IntervalMatrix A = ddhx.jacobian();
    IntervalMatrix mA = dhmx.jacobian();
    ARIADNE_LOG(6,"A="<<A<<" b="<<ddhx.value()<<"\n");

    IntervalVector rx = my * mA;
    for(uint j=0; j!=n; ++j) {
        rx[j] -= mmu*rec(mx[j]-xl[j]);
        rx[j] += mmu*rec(xu[j]-mx[j]);
    }
    IntervalVector ry = dhmx.value();
    IntervalVector rzl = esub(emul(IntervalVector(mx-make_exact(xl)),mzl),mmu);
    IntervalVector rzu = esub(emul(IntervalVector(make_exact(xu)-mx),mzu),mmu);
    ARIADNE_LOG(5,"rx="<<rx<<" ry="<<ry<<" rzl="<<rzl<<" rzu="<<rzu<<"\n");

    IntervalMatrix H(n,n);
    for(uint i=0; i!=m; ++i) { H += y[i] * ddhx[i].hessian(); }
    for(uint j=0; j!=n; ++j) {
        H[j][j] += mu*rec(sqr(x[j]-xl[j]));
        H[j][j] += mu*rec(sqr(xu[j]-x[j]));
    }

    // S = A Hinv AT
    // H dx + AT dy = rx; A dx = ry;
    //  dx = Hinv ( rx - AT dy )
    //  dy = Sinv ( A Hinv rx - ry )
    IntervalMatrix Hinv=inverse(H);
    ARIADNE_LOG(6,"H="<<H<<" Hinv="<<Hinv<<"\n");
    IntervalMatrix S=A*Hinv*transpose(A);
    IntervalMatrix Sinv=inverse(S);
    ARIADNE_LOG(6,"S="<<S<<" Sinv="<<Sinv<<"\n");
    IntervalVector dy = Sinv * ( A*(Hinv*rx) - ry );
    IntervalVector dx = Hinv * ( rx - dy * A);
    ARIADNE_LOG(5,"dx="<<dx<<" dy="<<dy<<"\n");

    IntervalVector nx = x-dx;
    IntervalVector ny = y-dy;
    ARIADNE_LOG(5,"nx="<<nx<<" ny="<<ny<<"\n");
    ARIADNE_LOG(6,"h(x)="<<h(nx)<<"\n");

    x = intersection(x,nx); y=intersection(y,ny);
    Interval nmu = 0;
    for(uint i=0; i!=m; ++i) { nmu += sqr(y[i]); }
    mu=intersection(mu,nmu);
}

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
        const X zero=arg[0].zero_element();
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
minimise(IntervalScalarFunction f, IntervalVector d, IntervalVectorFunction g, IntervalVector c) const
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
    t=Interval(+inf,+inf);
    for(uint j=0; j!=n; ++j) {
        t=min(t,c[j]-gy[j]);
        t=min(t,gy[j]-c[j]);
    }
    for(uint i=0; i!=m; ++i) {
        t=min(t,d[i]-y[i]);
        t=min(t,y[i]-d[i]);
    }

    // Find the range of possible values of the optimal t
    Float tmin=+inf;
    Float tmax=+inf;
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
minimisation_step(const IntervalScalarFunction& f, const IntervalVectorFunction& g,
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
