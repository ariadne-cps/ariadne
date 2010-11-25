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

#include "boost/multi_array.hpp"
#include "boost/array.hpp"

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

static int verbosity=0;

static const double error =  1e-2;

typedef bool Bool;
typedef tribool Tribool;
typedef tribool Tribool;

typedef Vector<Float> FloatVector;
typedef Matrix<Float> FloatMatrix;
typedef boost::numeric::ublas::vector_range<FloatVector> FloatVectorRange;

typedef Vector<Interval> IntervalVector;
typedef Matrix<Interval> IntervalMatrix;
typedef boost::numeric::ublas::vector_range<IntervalVector> IntervalVectorRange;

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

template<class X>
void simple_adat(Matrix<X>& S, const Matrix<X>& A, const Vector<X>& D)
{
    const uint m=A.row_size();
    const uint n=A.column_size();
    for(uint i1=0; i1!=m; ++i1) {
        for(uint j=0; j!=n; ++j) {
            for(uint i2=0; i2!=m; ++i2) {
                S[i1][i2]+=A[i1][j]*D[j]*A[i2][j];
            }
        }
    }
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

inline Vector<Interval> operator*(const Matrix<Float>& A, const Vector<Interval>& b) {
    return prod(A,b);
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


template<class X> Vector< Differential<X> > second_derivative(const RealVectorFunction& f, const Vector<X>& x) {
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
is_feasible_point(IntervalVector d, RealVectorFunction g, IntervalVector c, FloatVector y) const
{
    if(!contains(d,y)) { return false; }
    IntervalVector gy=g(IntervalVector(y));
    return subset(gy,c);
}


Bool OptimiserBase::
is_infeasibility_certificate(IntervalVector d, RealVectorFunction g, IntervalVector c, FloatVector x) const
{
    // Try to prove x.g(y) > 0
    const uint m=d.size();
    const uint n=c.size();
    VectorTaylorFunction tg(d,g);
    VectorTaylorFunction ti=VectorTaylorFunction::identity(d);
    ScalarTaylorFunction ts(d);
    for(uint i=0; i!=n; ++i) {
        ts+=x[i]*(tg[i]-c[i].upper())+x[i+n]*(c[i].lower()-tg[i]);
    }
    for(uint i=0; i!=m; ++i) {
        ts+=x[2*n+i]*(ti[i]-d[i].upper())+x[2*n+m+i]*(ti[i]-d[i].lower());
    }
    Interval xgd=ts(d);
    ARIADNE_LOG(2,"  x="<<x<<"  x.g(D)="<<xgd<<"\n");
    if(xgd.lower()>0.0) {
        return true;
    } else {
        return false;
    }
}



IntervalVector NonlinearInteriorPointOptimiser::
optimise(RealScalarFunction f, IntervalVector b, RealVectorFunction g, IntervalVector c) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


Tribool NonlinearInteriorPointOptimiser::
feasible(IntervalVector d, RealVectorFunction g, IntervalVector c) const
{
    ARIADNE_LOG(2,"NonlinearInteriorPointOptimiser::feasible(IntervalVector d, RealVectorFunction g, IntervalVector c)\n");
    ARIADNE_LOG(2,"  d="<<d<<", g="<<g<<", c="<<c<<"\n");

    ARIADNE_ASSERT(g.argument_size()==d.size());
    ARIADNE_ASSERT(g.result_size()==c.size());
    Float t;
    FloatVector x,y,z;

    setup_feasibility(d,g,c,x,y,z,t);

    // FIXME: Allow more steps
    for(uint i=0; i!=12; ++i) {
        ARIADNE_LOG(4,"  t="<<t<<", y="<<y<<", g(y)="<<g(y)<<", x="<<x<<", z="<<z<<"\n");
        this->feasibility_step(d,g,c,x,y,z,t);
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


void NonlinearInteriorPointOptimiser::feasibility_step (
        const FloatVectorFunction& g, FloatVector& x, FloatVector& y, FloatVector& z, Float& t) const
{
    const uint m=y.size();
    const uint n=x.size();

    Vector< Differential<Float> > ddg=g.evaluate(Differential<Float>::variables(2,y));

    // A is the transpose derivative matrix aij=dgj/dyi
    FloatMatrix A(m,n);
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    ARIADNE_LOG(9,"A="<<A<<"\n");

    // H is the Hessian matrix Hik = xj*dgj/dyidyk
    FloatMatrix H(m,m);
    for(uint j=0; j!=n; ++j) {
        add_hessian(H,x[j],ddg[j]);
    }
    ARIADNE_LOG(9," H="<<H<<"\n");

    FloatVector D(n);
    for(uint j=0; j!=n; ++j) { D[j]=x[j]/z[j]; }
    ARIADNE_LOG(9," D="<<D<<"\n");

    FloatMatrix& S=H;
    adat(S,A,D);
    ARIADNE_LOG(9,"S="<<S<<"\n");
    FloatMatrix Sinv=inverse(S);
    ARIADNE_LOG(9,"Sinv="<<Sinv<<"\n");

    FloatVector rx=prod(x,g.jacobian(y));
    FloatVector rz=g(y)+z; for(uint j=0; j!=n; ++j) { rz[j]+=t; }
    FloatVector rs(n); for(uint j=0; j!=n; ++j) { rs[j]=x[j]*z[j]; }
    Float rt=1; for(uint j=0; j!=n; ++j) { rt-=x[j]; }

    FloatVector tmp(n);
    for(uint j=0; j!=n; ++j) { tmp[j]=(rs[j]-rz[j]*x[j])/z[j]; }
    rx+=prod(A,tmp);

    for(uint j=0; j!=n; ++j) { rs[j]/=z[j]; }
    rx=prod(Sinv,rx);

    rz-=prod(rx,A);
    for(uint j=0; j!=n; ++j) { rs[j]-=rz[j]*x[j]/z[j]; }

    FloatVector& dx=rs;
    FloatVector& dy=rx;
    FloatVector& dz=rz;

    x-=dx;
    y-=dy;
    z-=dz;
}


void NonlinearInteriorPointOptimiser::
optimization_step(const FloatScalarFunction& f, const FloatVectorFunction& g,
                  FloatVector& x, FloatVector& y, FloatVector& z) const
{
    const uint m=y.size();
    const uint n=x.size();

    Differential<Float> ddf=f.evaluate(Differential<Float>::variables(2,y));
    Vector< Differential<Float> > ddg=g.evaluate(Differential<Float>::variables(2,y));

    // A is the transpose derivative matrix aij=dgj/dyi
    FloatMatrix A(m,n);
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    ARIADNE_LOG(9,"A="<<A<<"\n");

    // H is the Hessian matrix Hik = df/dyidyk - xj*dgj/dyidyk
    FloatMatrix H(m,m);
    set_hessian(H,ddf);
    for(uint j=0; j!=n; ++j) {
        add_hessian(H,-x[j],ddg[j]);
    }
    ARIADNE_LOG(9," H="<<H<<"\n");

    FloatVector D(n);
    for(uint j=0; j!=n; ++j) { D[j]=x[j]/z[j]; }
    ARIADNE_LOG(9," D="<<D<<"\n");

    FloatMatrix& S=H;
    adat(S,A,D);
    ARIADNE_LOG(9,"S="<<S<<"\n");
    FloatMatrix Sinv=inverse(S);
    ARIADNE_LOG(9,"Sinv="<<Sinv<<"\n");

    FloatVector rx=prod(x,g.jacobian(y));
    FloatVector rz=g(y)+z;
    FloatVector rs=emul(x,z);

    FloatVector tmp(n);
    for(uint j=0; j!=n; ++j) { tmp[j]=(rs[j]-rz[j]*x[j])/z[j]; }
    rx+=prod(A,tmp);

    for(uint j=0; j!=n; ++j) { rs[j]/=z[j]; }
    rx=prod(Sinv,rx);

    rz-=prod(rx,A);
    for(uint j=0; j!=n; ++j) { rs[j]-=rz[j]*x[j]/z[j]; }

    FloatVector& dx=rs;
    FloatVector& dy=rx;
    FloatVector& dz=rz;

    x-=dx;
    y-=dy;
    z-=dz;

}





void NonlinearInteriorPointOptimiser::feasibility_step(const IntervalVector& d, const FloatVectorFunction& g, const IntervalVector& c,
                                        FloatVector& x, FloatVector& y, FloatVector& z, Float& t) const
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
    ARIADNE_ASSERT(z.size()==o);

    Float mu=dot(x,z)/o;
    if(!egtr(emul(x,z),gamma*mu)) {
        if(verbosity>=1) { ARIADNE_WARN("Near-degeneracy in Lyapunov multipliers in interior-point solver:\n  x="<<x<<", y="<<y<<", z="<<z<<"\n"); }
        x=(1-sigma)*x+FloatVector(x.size(),sigma/x.size());
        mu=dot(x,z)/o;
    }

    FloatVector yt=join(y,t);
    ARIADNE_LOG(9,"m="<<m<<" n="<<n<<"\n");
    ARIADNE_LOG(9,"x="<<x<<" yt="<<yt<<" z="<<z<<"\n");

    Vector< Differential<Float> > ddg=g.evaluate(Differential<Float>::variables(2,y));
    ARIADNE_LOG(9,"  ddg="<<ddg<<"\n");

    // gy is the vector of values of g(y)
    FloatVector gy(n); for(uint j=0; j!=n; ++j) { gy[j]=ddg[j].value(); }
    ARIADNE_LOG(9,"  g(y)="<<gy<<" ");

    // A is the transpose derivative matrix aij=dgj/dyi, extended with a column of ones
    FloatMatrix A(m,n);
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=n; ++j) {
            A[i][j]=ddg[j][i];
        }
    }
    ARIADNE_LOG(9," A="<<A<<" ");

    // H is the Hessian matrix Hik = (xcuj-xclj)*dgj/dyidyk
    FloatMatrix H(m,m);
    for(uint j=0; j!=n; ++j) {
        add_hessian(H,x[j]-x[n+j],ddg[j]);
    }
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
        NonlinearInteriorPointOptimiser::compute_z(d,g,c,ny,nt,nz);
        allpositive = egtr(nx,0.0) && egtr(nz,0.0) && egtr(emul(nx,nz),gamma*mu);
    }
    ARIADNE_LOG(9,"alpha="<<alpha<<"\n");
    ARIADNE_LOG(9,"nx="<<nx<<" nyt="<<nyt<<" nz="<<nz<<" nxz="<<eivl(emul(nx,nz))<<"\n");

    x=nx; y=project(nyt,range(0,m)); z=nz; t=nyt[m];
}

void NonlinearInteriorPointOptimiser::linearised_feasibility_step(const IntervalVector& d, const RealVectorFunction& g, const IntervalVector& c,
                                                   FloatVector& x, FloatVector& y, FloatVector& z, Float& t) const
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
    ARIADNE_ASSERT(z.size()==o);

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
        NonlinearInteriorPointOptimiser::compute_z(d,g,c,ny,nt,nz);
        allpositive = egtr(nx,0.0) && egtr(nz,0.0) && egtr(emul(nx,nz),gamma*mu);
    }
    ARIADNE_LOG(9,"alpha="<<alpha<<"\n");
    ARIADNE_LOG(9,"nx="<<nx<<" nyt="<<nyt<<" nz="<<nz<<" nxz="<<eivl(emul(nx,nz))<<"\n");

    x=nx; y=project(nyt,range(0,m)); z=nz; t=nyt[m];

}



void NonlinearInteriorPointOptimiser::
setup_feasibility(const IntervalVector& d, const FloatVectorFunction& g, const IntervalVector& b,
                  FloatVector& x, FloatVector& y, FloatVector& z, Float& t) const
{
    const uint l=2*(d.size()+b.size());
    y=midpoint(d);
    x=FloatVector(l,1.0/l);
    z.resize(l);
    compute_tz(d,g,b,y,t,z);
}


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










struct KuhnTuckerFunctionBody : VectorFunctionMixin<KuhnTuckerFunctionBody,Real>
{
    RealScalarFunction f;
    Array<RealScalarFunction> g;
    Array<RealScalarFunction> df;
    Array<Array<RealScalarFunction> > dg;

    KuhnTuckerFunctionBody(RealScalarFunction _f, RealVectorFunction _g) {
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
    RealScalarFunction operator[](uint) const { ARIADNE_NOT_IMPLEMENTED; }
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

struct FeasibilityKuhnTuckerFunctionBody : VectorFunctionMixin<FeasibilityKuhnTuckerFunctionBody,Real>
{
    Array<RealScalarFunction> g;
    Array<Array<RealScalarFunction> > dg;

    FeasibilityKuhnTuckerFunctionBody(RealVectorFunction _g) {
        const uint m=_g.argument_size();
        const uint n=_g.result_size();
        g.resize(n); dg.resize(n); for(uint j=0; j!=n; ++j) { dg[j].resize(m); }
        for(uint j=0; j!=n; ++j) { g[j]=_g[j]; for(uint i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    uint result_size() const { return g.size()*2+g[0].argument_size()+1; }
    uint argument_size() const { return g.size()*2+g[0].argument_size()+1; }
    RealScalarFunction operator[](uint) const { ARIADNE_NOT_IMPLEMENTED; }
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



struct ConstrainedFeasibilityKuhnTuckerFunctionBody : VectorFunctionMixin<FeasibilityKuhnTuckerFunctionBody,Real>
{
    uint m;
    uint n;
    IntervalVector d;
    Array<RealScalarFunction> g;
    IntervalVector c;
    Array<Array<RealScalarFunction> > dg;

    ConstrainedFeasibilityKuhnTuckerFunctionBody(IntervalVector D, RealVectorFunction _g, IntervalVector C) {
        m=_g.argument_size();
        n=_g.result_size();
        d=D; c=C;
        g.resize(n); dg.resize(n); for(uint j=0; j!=n; ++j) { dg[j].resize(m); }
        for(uint j=0; j!=n; ++j) { g[j]=_g[j]; for(uint i=0; i!=m; ++i) { dg[j][i]=g[j].derivative(i); } }
    }

    uint result_size() const { return 5*m+4*n+1u; }
    uint argument_size() const { return 5*m+4*n+1u; }
    RealScalarFunction operator[](uint) const { ARIADNE_NOT_IMPLEMENTED; }
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
optimise(RealScalarFunction f, IntervalVector d, RealVectorFunction g, IntervalVector c) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


Tribool KrawczykOptimiser::
feasible(IntervalVector d, RealVectorFunction g, IntervalVector c) const
{
    ARIADNE_LOG(2,"KrawczykOptimiser::feasible(IntervalVector d, RealVectorFunction g, IntervalVector c)\n");
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



void KrawczykOptimiser::setup_feasibility(const IntervalVector& d, const RealVectorFunction& g, const IntervalVector& c,
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


void KrawczykOptimiser::compute_tz(const IntervalVector& d, const RealVectorFunction& g, const IntervalVector& c,
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
optimisation_step(const RealScalarFunction& f, const RealVectorFunction& g,
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



void KrawczykOptimiser::feasibility_step(const RealVectorFunction& g,
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
void KrawczykOptimiser::feasibility_step(const IntervalVector& d, const RealVectorFunction& g, const IntervalVector& c,
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




void KrawczykOptimiser::feasibility_step(const IntervalVector& d, const RealVectorFunction& g, const IntervalVector& c,
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




} // namespace Ariadne
