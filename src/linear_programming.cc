/***************************************************************************
 *            linear_programming.cc
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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
 
#include "config.h"

#include "tuple.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "linear_programming.h"

#include "macros.h"
#include "logging.h"

namespace Ariadne {

struct SingularLinearProgram : std::runtime_error { 
    SingularLinearProgram(const std::string& what)
        : std::runtime_error(what) { };
};

struct UnboundedLinearProgram : std::runtime_error { 
    UnboundedLinearProgram(const std::string& what)
        : std::runtime_error(what) { }
};

std::ostream& operator<<(std::ostream& os, VariableType t) {
    return os << (t==BASIS ? 'B' : t==LOWER ? 'L' : t==UPPER ? 'U' : '?');  
}

// Compare a matrix with a projection
// Used for checking whether a projection is the identity
template<class Mx, class X> bool operator==(const boost::numeric::ublas::matrix_range<Mx>& A, const boost::numeric::ublas::matrix<X>& B)
{
    if(A.size1()!=B.size1() || A.size2()!=B.size2()) {
        return false;
    }


    for(size_t i=0; i!=A.size1(); ++i) {
        for(size_t j=0; j!=A.size1(); ++j) {
            if(A(i,j)!=B(i,j)) { 
                return false; 
            }
        }
    }
    return true;
}


// Extend an array of size m to an array of size n 
// such that the first m elements are the same, 
// and the new array contains the elements [0,n)
array<size_t> extend_p(const array<size_t>& p, const size_t n)
{
    const size_t m=p.size();
    array<size_t> q(n);
    for(size_t j=0; j!=n; ++j) {
        q[j]=n;
    }
    for(size_t k=0; k!=m; ++k) {
        ARIADNE_ASSERT(p[k]<n);
        ARIADNE_ASSERT(q[p[k]]==n);
        q[p[k]]=k;
    }
    size_t k=m;
    for(size_t j=0; j!=n; ++j) {
        if(q[j]==n) { q[j]=k; ++k; }
    }
    array<size_t> r(n);
    for(size_t j=0; j!=n; ++j) {
        r[q[j]]=j;
    }
    for(size_t i=0; i!=m; ++i) { 
        ARIADNE_ASSERT(p[i]==r[i]);
    }
    return r;
}
      


// Compute a basis (p_1,\ldots,p_m) for the matrix A
// Assume that the matrix A has full row rank
template<class X>
pair< array<size_t>, Matrix<X> >
compute_basis(const Matrix<X>& A)
{
    ARIADNE_LOG(9,"compute_basis(A) with A="<<A<<"\n");
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    array<size_t> p(n);
    Matrix<X> I=Matrix<X>::identity(m);

 
    if(project(A,range(0,m),range(0,m))==I) {
        for(size_t i=0; i!=n; ++i) { p[i]=i; } 
        return make_pair(p,I); 
    }

    if(project(A,range(0,m),range(n-m,n))==I) {
        for(size_t i=0; i!=m; ++i) { p[i]=i+(n-m); } 
        for(size_t i=m; i!=n; ++i) { p[i]=i-m; } 
        return make_pair(p,I); 
    }

    Matrix<X> L=I;
    Matrix<X> U=A;
    
    for(size_t j=0; j!=n; ++j) { p[j]=j; }

    for(size_t k=0; k!=m; ++k) {
        // Find a good pivot row and swap entries
        if(U[k][k]==0) {
            size_t j=k+1;
            while(j<n && U[k][j]==0) {
                ++j;
            }
            if(j==n) { ARIADNE_THROW(SingularLinearProgram,"compute_basis"," matrix A="<<A<<" is singular"); }
            for(size_t i=0; i!=m; ++i) {
                std::swap(U[k][k],U[k][i]);
            }
            std::swap(p[k],p[j]);
        }
        
        // Update LU factorisation
        X r  = 1/U[k][k];
        for(size_t i=k+1; i!=m; ++i) {
            X s=U[i][k] * r;
            for(size_t j=0; j!=m; ++j) {
                L[i][j] -= s * L[k][j];
            }
            for(size_t j=k+1; j!=n; ++j) {
                U[i][j] -= s * U[k][j];
            }
            U[i][k] = 0;
        }
        for(size_t j=0; j!=m; ++j) {
            L[k][j] *= r;
        }
        for(size_t j=k+1; j!=n; ++j) {
            U[k][j] *= r;
        }
        U[k][k] = 1;
    }

    // Backsubstitute to find inverse of pivot columns

    for(size_t k=m; k!=0; ) {
        --k;
        for(size_t i=0; i!=k; ++i) {
            X s=U[i][k];
            for(size_t j=0; j!=m; ++j) {
                L[i][j] -= s * L[k][j];
            }
            U[i][k] = 0;
        }
    }

    return make_pair(p,L);

}
    

    
// Check that the matrix B is the inverse of the matrix A_B with columns p[0],...,p[m-1] of A.
template<class X>
void consistency_check(const Matrix<X>& A, const array<size_t>& p, const Matrix<X>& B) 
{
    static const X MAXIMUM_ERROR=1e-8;
    const size_t m=A.row_size();
    Matrix<X> A_B(m,m);
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        for(size_t i=0; i!=m; ++i) {
            A_B[i][k]=A[i][j];
        }
    }
    
    array<size_t> p_B(p.begin(),p.begin()+m);
    
    Matrix<X> Z=prod(B,A_B);
    ARIADNE_LOG(9,"        p_B="<<p_B<<" B="<<B<<" A_B="<<A_B<<" B*A_B-I="<<Z<<"\n");
    for(size_t i=0; i!=m; ++i) { Z[i][i]-=1; }
    ARIADNE_ASSERT(norm(Z)<MAXIMUM_ERROR);
}


// Check that the matrix B is the inverse of the matrix A_B with columns p[0],...,p[m-1] of A and that Ax=b.
template<class X>
void consistency_check(const Matrix<X>& A, const Vector<X>& b, const array<size_t>& p, const Matrix<X>& B, const Vector<X>& x) 
{
    static const X MAXIMUM_ERROR=1e-8;
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    Matrix<X> A_B(m,m);
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        for(size_t i=0; i!=m; ++i) {
            A_B[i][k]=A[i][j];
        }
    }
    
    Matrix<X> Z=prod(B,A_B);
    ARIADNE_LOG(9,"        B="<<B<<" A_B="<<A_B<<" B*A_B-I="<<Z<<"\n");
    for(size_t i=0; i!=m; ++i) { Z[i][i]-=1; }
    ARIADNE_ASSERT(norm(Z)<MAXIMUM_ERROR);

    Vector<X> z=prod(A,b);
    for(size_t j=0; j!=n; ++j) { z[j]-=x[j]; }
    ARIADNE_ASSERT(norm(z)<MAXIMUM_ERROR);
}


// Check that the matrix B is the inverse of the matrix A_B with columns p[0],...,p[m-1] of A, and that 
// the vector x is given by x_L=l_L, x_U=x_U and x_B=B^{-1} A_N x_N.
template<class X>
void consistency_check(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, const array<VariableType>& vt, const array<size_t>& p, const Matrix<X>& B, const Vector<X>& x) 
{
    ARIADNE_LOG(9,"        Checking consistency of B and x\n");
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Matrix<X> A_B(m,m);
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        for(size_t i=0; i!=m; ++i) {
            A_B[i][k]=A[i][j];
        }
    }
    
    array<size_t> p_B(p.begin(),p.begin()+m);
    
    Matrix<X> Z=prod(B,A_B);
    ARIADNE_LOG(9,"          p_B="<<p_B<<" B="<<B<<" A_B="<<A_B<<" B*A_B="<<Z<<"\n");
    for(size_t i=0; i!=m; ++i) { Z[i][i]-=1; }

    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        ARIADNE_ASSERT(vt[j]==LOWER || vt[j]==UPPER);
        X xj = (vt[j]==LOWER ? l[j] : u[j]);
        ARIADNE_ASSERT(x[j]==xj);
    }
    Vector<X> z=prod(A,x);
    ARIADNE_LOG(9,"          A="<<A<<" x="<<x<<" b="<<b<<" Ax="<<z<<"\n");
    z-=b;
    ARIADNE_ASSERT(norm(z)<1e-5);
}


// Compute the variable types from the permutation, taking m basic variables and all non-basic variables lower.
template<class X>
array<VariableType>
compute_vt(const Vector<X>& l, const Vector<X>& u, const array<size_t>& p, const size_t m)
{
    const size_t n=p.size();
    array<VariableType> vt(n);
    for(size_t k=0; k!=m; ++k) { 
        vt[p[k]]=BASIS; 
    }
    for(size_t k=m; k!=n; ++k) { 
        if(l[p[k]]==-inf<X>()) { 
            vt[p[k]] = UPPER;
        } else {
            vt[p[k]] = LOWER;
        }
    }
    return vt;
}

array<size_t>
compute_p(const array<VariableType>& tv)
{
    const size_t n=tv.size();
    array<size_t> p(n);
    size_t k=0;
    for(size_t j=0; j!=n; ++j) {
        if(tv[j]==BASIS) { p[k]=j; ++k; }
    }
    for(size_t j=0; j!=n; ++j) {
        if(tv[j]!=BASIS) { p[k]=j; ++k; }
    }
    return p;
}


template<class XX, class X>
Matrix<XX> compute_B(const Matrix<X>& A, const array<size_t>& p)
{
    const size_t m=A.row_size();
    Matrix<XX> A_B(m,m);
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        for(size_t i=0; i!=m; ++i) {
            A_B[i][k]=A[i][j];
        }
    }
    
    Matrix<XX> B=inverse(A_B);
    
    return B;
}

template<class X>
Vector<X> compute_c(const Matrix<X>& A, const array<size_t>& p, const Vector<X>& x)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    Vector<X> c(n);
    for(size_t k=m; k!=n; ++k) {
        if(x[p[k]]<0) { c[p[k]]=-1; }
    }
    return c;
}


template<class X, class XX>
Vector<XX> compute_x(const Matrix<X>& A, const Vector<X>& b, const array<size_t>& p, const Matrix<XX>& B)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Vector<XX> x(n);

    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        x[j]=0;
        for(size_t i=0; i!=m; ++i) {
            x[j]+=B[k][i]*b[i];
        }
    }

    return x;
}


template<class X, class XX>
Vector<XX>
compute_x(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, const array<VariableType>& vt, const array<size_t>& p, const Matrix<XX>& B)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    
    Vector<XX> w(m);
    Vector<XX> x(n);
    
    // Compute x_N
    for(size_t j=0; j!=n; ++j) {
        if(vt[j]==LOWER) { x[j]=l[j]; }
        else if(vt[j]==UPPER) { x[j]=u[j]; }
        else { x[j]=0; }
    }
    ARIADNE_LOG(9,"  x_N="<<x);

    // Compute w=b-A_N x_N
    for(size_t i=0; i!=m; ++i) {
        w[i]=b[i];
        for(size_t k=m; k!=n; ++k) {
            size_t j=p[k];
            w[i]-=A[i][j]*x[j];
        }
    }
    ARIADNE_LOG(9," w="<<w);

    // Compute x_B=B w
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        x[j]=0;
        for(size_t i=0; i!=m; ++i) {
            x[j]+=B[k][i]*w[i];
        }
    }

    ARIADNE_LOG(9," x="<<x<<"\n");

    Vector<XX> Axmb=prod(A,x)-b;
    ARIADNE_ASSERT(norm(Axmb)<0.00001);
    return x;
}


template<class X,class XX>
std::pair<Vector<XX>,Vector<XX> >
compute_wx(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, const array<size_t>& p, const Matrix<XX>& B)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    
    Vector<XX> w(m);
    Vector<XX> x(n);

    // Compute x_N
    for(size_t j=0; j!=n; ++j) {
        if(vt[j]==LOWER) { x[j]=l[j]; }
        else if(vt[j]==UPPER) { x[j]=u[j]; }
        else { x[j]=0; }
    }
    ARIADNE_LOG(9,"  x_N="<<x);

    // Compute w=b-A_N x_N
    for(size_t i=0; i!=m; ++i) {
        w[i]=b[i];
        for(size_t k=m; k!=n; ++k) {
            size_t j=p[k];
            w[i]-=A[i][j]*x[j];
        }
    }
    ARIADNE_LOG(9," w="<<w);

    // Compute x_B=B w
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        x[j]=0;
        for(size_t i=0; i!=m; ++i) {
            x[j]+=B[k][i]*w[i];
        }
    }

    ARIADNE_LOG(9," x="<<x<<"\n");

    Vector<X> Axmb=prod(A,x)-b;
    ARIADNE_ASSERT(norm(Axmb)<0.00001);
    return make_pair(w,x);
}


template<class X,class XX>
Vector<XX>
compute_y(const Vector<X>& c, const array<size_t>& p, const Matrix<XX>& B)
{
    const size_t m=B.row_size();
    Vector<XX> y(m);
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        for(size_t i=0; i!=m; ++i) {
            y[i]+=c[j]*B[k][i];
        }
    }
    return y;
}

template<class X, class XX>
Vector<XX>
compute_z(const Matrix<X>& A, const Vector<X>& c, const array<size_t>& p, const Vector<XX>& y)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    Vector<XX> z(n);
    for(size_t k=0; k!=m; ++k) {
        z[p[k]]=0;
    }
    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        z[j]=c[j];
        for(size_t i=0; i!=m; ++i) {
            z[j]-=y[i]*A[i][j];
        }
    }
    return z;
}

template<class X>
size_t
compute_s(const size_t m, const array<size_t>& p, const Vector<X>& z) 
{
    const size_t n=z.size();
    for(size_t k=0; k!=n; ++k) {
        if(z[p[k]]<0) { return k; }
    }
    return n;
}

template<class X>
size_t
compute_s(const size_t m, const array<VariableType>& vt, const array<size_t>& p, const Vector<X>& z) 
{
    const size_t n=z.size();
    for(size_t k=m; k!=n; ++k) {
        if( (vt[p[k]]==LOWER && z[p[k]]<0)
            || (vt[p[k]]==UPPER && z[p[k]]>0) ) { return k; }
    }
    return n;
}

// Compute the direction in which the basic variables move
// given by d=-B*A_s
template<class X>
Vector<X>
compute_d(const Matrix<X>& A, const array<size_t>& p, const Matrix<X>& B, const size_t ks)
{
    const size_t m=A.row_size();
    size_t js=p[ks];
    Vector<X> d(m);
    for(size_t k=0; k!=m; ++k) {
        for(size_t i=0; i!=m; ++i) {
            d[k]-=B[k][i]*A[i][js];
        }
    }
    return d;
}


template<class X>
pair<size_t,X>
compute_rt(const array<size_t>& p, const Vector<X>& x, const Vector<X>& d)
{
    const size_t m=d.size();
    X t=inf<X>();
    size_t r=m;
    for(size_t k=0; k!=m; ++k) {
        if(d[k]<0 && x[p[k]]>=0) { 
            X tk=(x[p[k]])/(-d[k]); 
            if(r==m || tk<t) { 
                t=tk; 
                r=k; 
            } 
        }
    }
    return make_pair(r,t);
}


template<class X>
std::pair<size_t,X>
compute_rt(const Vector<X>& l, const Vector<X>& u, const array<VariableType>& vt, const array<size_t>& p, const Vector<X>& x, const Vector<X>& d, const size_t s)
{
    const size_t m=d.size();
    size_t r=s;
    X ds=(vt[p[s]]==LOWER ? +1 : -1);
    X t=u[p[s]]-l[p[s]];
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if(d[k]<0 && x[j]>=l[j] && l[j] != -inf<X>()) { 
            X tk=(l[j]-x[j])/(ds*d[k]); 
            if(t==inf<X>() || tk<t) { t=tk; r=k; } 
        } else if( d[k]>0 && x[j]<=u[j] && u[j] != inf<X>() ) { 
            X tk=(u[j]-x[j])/(ds*d[k]); 
            if(t==inf<X>() || tk<t) { t=tk; r=k; } 
        }
    }
    t*=ds;
    ARIADNE_LOG(9,"       r="<<r<<" t="<<t<<"\n");

    return make_pair(r,t);
}


template<class X>
void
update_B(Matrix<X>& B, const Vector<X>& d, const size_t r)
{
    const size_t m=d.size();
    X dr=d[r];
    X drr=1/dr;
    Vector<X> e(m); e[r]=1;
    Vector<X> Br(m); for(uint j=0; j!=m; ++j) { Br[j]=B[r][j]; }
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=m; ++j) {
            B[i][j]-=(d[i]+e[i])*Br[j]*drr;
        }
    }
    return;
}


template<class X>
void
update_x(const array<size_t>& p, Vector<X>& x, const size_t s, const Vector<X>& d, const size_t r, const X& t)
{
    const size_t m=d.size();
    for(size_t i=0; i!=m; ++i) {
        x[p[i]]+=t*d[i];
    }
    x[p[r]] = 0.0;
    x[p[s]] = t;
    return;
}


template<class X>
void
update_x(const Vector<X>& l, const Vector<X>& u, const array<size_t>& p, Vector<X>& x, const size_t s, const Vector<X>& d, const size_t r, const X& t)
{
    ARIADNE_ASSERT(t>=0);
    const size_t m=d.size();
    for(size_t i=0; i!=m; ++i) {
        x[p[i]]+=t*d[i];
    }
    x[p[s]] += t;
    if(d[r]<0) { x[p[r]]=l[p[r]]; }
    else if(d[r]>0) { x[p[r]]=u[p[r]]; }
    return;
}

template<class X>
void
update_x(const Vector<X>& l, const Vector<X>& u, const array<size_t>& p, Vector<X>& x, const size_t ks, const Vector<X>& d, const X& t)
{
    ARIADNE_ASSERT(t>=0);
    const size_t m=d.size();
    for(size_t i=0; i!=m; ++i) {
        x[p[i]]+=t*d[i];
    }
    x[p[ks]]+=t;

    size_t js=p[ks];
    if(t>0) { x[js]=u[js]; } 
    else { x[js]=l[js]; }
}




template<class X>
bool 
lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, array<size_t>& p, Matrix<X>& B, Vector<X>& x)
{
    ARIADNE_LOG(6,"  lpstep(A,b,c,p,B,x)\n      A="<<A<<"\n      b="<<b<<"\n      c="<<c<<"\n      p="<<p<<"\n      B="<<B<<"\n      x="<<x<<"\n");

    const size_t m=A.row_size();
    const size_t n=A.column_size();
    Vector<X> y=compute_y(c,p,B);
    Vector<X> z=compute_z(A,c,p,y);
    ARIADNE_LOG(7,"      y="<<y<<"\n      z="<<z<<"\n");

    // Compute variable p[s] to enter the basis
    // If z is positive, we are done
    size_t s=compute_s(m,p,z);

    if(s==n) { 
        ARIADNE_LOG(5,"  No improvement possible.\n"); 
        return true; 
    }
    

    // Compute the vector d in which the current basic variables move
    // if the non-basic variable s increases by 1
    // d is given by d=-B*As where As is the s-th column of A 
    Vector<X> d=compute_d(A,p,B,s);
    
    // Compute distance t along d in which to move,
    // and the variable p[r] to leave the basis
    // If any variables are infeasible, we leave them infeasible 
    // and do not allow them to leave the basis
    size_t r; X t;
    make_lpair(r,t)=compute_rt(p,x,d);

    ARIADNE_LOG(7,"      r="<<r<<" s="<<s<<" t="<<t<<" d="<<d<<"\n");

    if(r==m) { 
        ARIADNE_LOG(5,"  Unbounded linear program.\n"); 
        ARIADNE_THROW(UnboundedLinearProgram,"feasible",""); 
    }

    
    ARIADNE_LOG(5,"  Swapping basic variable "<<p[r]<<" in position "<<r<<" with non-basic variable "<<p[s]<<" in position "<<s<<"\n");


    // Compute new value v of x[p[s]] and make sure that it's still feasible
    // Update matrix of inverses
    update_B(B,d,r);

    // Update state vector x
    update_x(p,x,s,d,r,t);

    // Update pivots
    std::swap(p[r],p[s]);

    ARIADNE_LOG(7,"      p ="<<p<<"\n");
    ARIADNE_LOG(7,"      B ="<<B<<"\n");
    ARIADNE_LOG(7,"      x ="<<x<<"\n");
    
    consistency_check(A,b,p,B,x);

    return false;
}


template<class X>
bool lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x)
{
    ARIADNE_LOG(9,"  lpstep(A,b,c,l,u,vt,p,V,x)\n    A="<<A<<" b="<<b<<" c="<<c<<"\n    p="<<p<<" B="<<B<<"\n    vt="<<vt<<" l="<<l<<" x="<<x<<" u="<<u<<"\n");

    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Vector<X> y=compute_y(c,p,B);
    Vector<X> z=compute_z(A,c,p,y);
    
    ARIADNE_LOG(9,"      y="<<y<<" z="<<z<<"\n");

    // Compute variable p[s] to enter the basis
    size_t s=compute_s(m,vt,p,z);
    if(s==n) { ARIADNE_LOG(9,"  No improvement possible.\n"); return true; }
    
    // Compute direction d in which to move the current basic variables
    // as the variable entering the basis changes by +1
    Vector<X> d=compute_d(A,p,B,s);
   

    // Compute distance t along d in which to move,
    // and the variable p[r] to leave the basis
    // The bounds on t are given by l <= x + t * d <= u 
    // Note that t is negative if an upper variable enters the basis
    size_t r; X t;
    make_lpair(r,t)=compute_rt(l,u,vt,p,x,d,s);
    ARIADNE_LOG(9,"      s="<<s<<" d="<<d<<" t="<<t);


    if(r==s) {
        VariableType nvts=(vt[p[s]]==LOWER ? UPPER : LOWER);
        ARIADNE_LOG(5,"   Changing non-basic variable "<<p[s]<<" in position "<<s<<" from type "<<vt[p[s]]<<" to type "<<nvts<<"\n");
    } else {
        ARIADNE_LOG(5,"   Swapping basic variable "<<p[r]<<" in position "<<r<<" with non-basic variable "<<p[s]<<" in position "<<s<<"\n");
    }


    if(r==s) {
        // Constraint is due to bounds on x_s
        // No change in basic variables or inverse basis matrix
        update_x(l,u,p,x,s,d,t);

        // Update variable type
        if(vt[p[s]]==LOWER) { vt[p[s]]=UPPER; }
        else { vt[p[s]]=LOWER; }
        
    } else {
        // Variable p[r] should leave basis, and variable p[s] enter
        ARIADNE_ASSERT(r<m);
        
        update_B(B,d,r);
        update_x(l,u,p,x,s,d,r,t);
   
    
        // Update pivots and variable types
        std::swap(p[r],p[s]);
        vt[p[r]] = BASIS; 
        vt[p[s]] = (d[r]>0 ? UPPER : LOWER);
    }

    ARIADNE_LOG(7,"      vt="<<vt<<"\n      p="<<p<<"\n");
    ARIADNE_LOG(7,"      B="<<B<<"\n      x="<<x<<"\n");
   
    consistency_check(A,b,l,u,vt,p,B,x);

    return false;
}





template<class X>
tribool _primal_feasible(const Matrix<X>& A, const Vector<X>& b, array<size_t>& p, Matrix<X>& B)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Vector<X> x=compute_x(A,b,p,B);
    Vector<X> c(n);
    bool infeasible=false;
    for(size_t j=0; j!=n; ++j) {
        if(x[j]<0) { c[j]=-1; infeasible=true; }
        else { c[j]=0; }
    }

    ARIADNE_LOG(5,"    x="<<x<<"\n    c="<<c<<"\n");

    while(infeasible) {
        infeasible=false;

        bool done=lpstep(A,b,c,p,B,x);
        ARIADNE_LOG(5,"    p="<<p<<"\n    B="<<B<<"\n");
        
        if(done) {
            return false;
        }

        for(size_t j=0; j!=n; ++j) {
            if(x[j]<0) { c[j]=-1; infeasible=true; }
            else { c[j]=0; }
        }
        ARIADNE_LOG(5,"    x="<<x<<"\n    c="<<c<<"\n");
    }
 
    // Check solution
    for(size_t i=0; i!=n; ++i) {
        ARIADNE_ASSERT(x[i]>=0.0);
    }
    Vector<X> Ax=prod(A,x);
    for(size_t i=0; i!=m; ++i) {
         ARIADNE_ASSERT(Ax[i]==b[i]);
    }

    ARIADNE_LOG(9,"\nFeasible point x="<<x<<"\n Ax="<<Vector<X>(prod(A,x))<<" b="<<b<<"\n");

    return true;
} 






template<class X>
tribool _constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Vector<X> c(n);
    Vector<X> ll(l);
    Vector<X> uu(u);

    bool infeasible=false;
    for(size_t j=0; j!=n; ++j) {
        if(x[j]<l[j]) { c[j]=-1; ll[j]=-inf<X>(); infeasible=true; }
        else if(x[j]>u[j]) { c[j]=+1; uu[j]=+inf<X>(); infeasible=true; }
        else { c[j]=0; }
    }
    ARIADNE_LOG(9,"    vt="<<vt<<" x="<<x<<" c="<<c<<"\n");

    while(infeasible) {
        infeasible=false;

        bool done=lpstep(A,b,c,ll,uu,vt,p,B,x);
        ARIADNE_LOG(9,"  Done changing basis\n");
        ARIADNE_LOG(9,"    p="<<p<<" B="<<B<<"\n");
        ARIADNE_LOG(9,"    vt="<<vt<<" x="<<x<<"\n");
        
        if(done) {
            ARIADNE_LOG(9,"  Cannot put infeasible variables into basis.");
            Vector<X> y=compute_y(c,p,B);
            Vector<X> yA=prod(y,A);
            ARIADNE_LOG(9,"\nCertificate of infeasibility:\n y="<<y<<"\n yA="<<yA<<"\n");
            return false;
        }

        for(size_t j=0; j!=n; ++j) {
            if(vt[j]==LOWER) { ARIADNE_ASSERT(x[j]==l[j]); }
            if(vt[j]==UPPER) { ARIADNE_ASSERT(x[j]==u[j]); }
            if(x[j]<l[j]) { c[j]=+1; ll[j]=-inf<X>(); infeasible=true; }
            else if(x[j]>u[j]) { c[j]=-1; uu[j]=+inf<X>(); infeasible=true; }
            else { c[j]=0; }
        }
        ARIADNE_LOG(9,"\n    vt="<<vt<<" x="<<x<<" c="<<c<<"\n");
    }
 
    ARIADNE_LOG(9,"  Checking solution\n");
    // Check solution
    for(size_t i=0; i!=n; ++i) {
        ARIADNE_ASSERT(x[i]>=l[i]);
        ARIADNE_ASSERT(x[i]<=u[i]);
    }
    Vector<X> Ax=prod(A,x);
    for(size_t i=0; i!=m; ++i) {
        ARIADNE_ASSERT(abs(Ax[i]-b[i])<0.0001);
    }

    ARIADNE_LOG(9,"\nFeasible point x="<<x<<"); l="<<l<<" u="<<u<<"\n Ax="<<Vector<X>(prod(A,x))<<" b="<<b);

    return true;
} 



// Check for feasibility of yA<=c. The solution is infeasible if there exists x with cx<0, Ax=0 and x>=0
template<class X>
tribool _dual_feasible(const Matrix<X>& A, const Vector<X>& c, array<size_t>& p, Matrix<X>& B) 
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    Vector<X> b(m);
    Vector<X> x(n);
    Vector<X> y(m);
    Vector<X> z(n);
    
    bool done;
    do {
        lpstep(A,b,c,p,B,x);
        y=compute_y(c,p,B);
        z=compute_z(A,c,p,y);
        
        done=true;
        for(size_t i=0; i!=n; ++i) {
            if(z[i]<0) { done=false; }
        }
        if(done) { return true; }
    } while(true);
    return false;
}


// Check for feasibility of Ax=b x>=0
template<class X>
tribool primal_feasible(const Matrix<X>& A, const Vector<X>& b)
{
    ARIADNE_LOG(2,"primal_feasible(A,b)\n");
    ARIADNE_LOG(3,"    A="<<A<<"\n    b="<<b<<"\n");

    ARIADNE_ASSERT(b.size()==A.row_size());
 
    array<size_t> p;
    Matrix<X> B;
    make_lpair(p,B)=compute_basis(A);

    return _primal_feasible(A,b,p,B);
}


// Check for feasibility of yA<=c. The solution is infeasible if there exists x with cx<0, Ax=0 and x>=0. 
template<class X>
tribool dual_feasible(const Matrix<X>& A, const Vector<X>& c) 
{
    ARIADNE_LOG(2,"dual_feasible(Matrix A, Vector c)\n");
    ARIADNE_LOG(3,"  A="<<A<<" c="<<c<<"\n");

    ARIADNE_ASSERT(c.size()==A.column_size());

    array<size_t> p;
    Matrix<X> B;
    make_lpair(p,B)=compute_basis(A);

    return _dual_feasible(A,c,p,B);
}



    

// Check for feasibility of Ax=b l<=b<=u
template<class X>
tribool constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u)
{
    ARIADNE_LOG(2,"constrained_feasible(A,b,l,u)\n");
    ARIADNE_LOG(3,"    A="<<A<<"\n    b="<<b<<"\n    l="<<l<<"\n    u="<<u<<"\n");
    ARIADNE_ASSERT(b.size()==A.row_size());
    ARIADNE_ASSERT(l.size()==A.column_size());
    ARIADNE_ASSERT(u.size()==A.column_size());

    array<size_t> p;
    Matrix<X> B;
    make_lpair(p,B)=compute_basis(A);

    array<VariableType> vt=compute_vt(l,u,p,A.row_size());

    ARIADNE_LOG(9,"    p="<<p<<" B="<<B<<"  (BA="<<Matrix<X>(prod(B,A))<<")\n");

    Vector<X> x=compute_x(A,b,l,u,vt,p,B);

    tribool fs = _constrained_feasible(A,b,l,u,vt,p,B,x);
    tribool vfs = verify_constrained_feasibility(A,b,l,u,vt);
    ARIADNE_ASSERT(indeterminate(vfs) || vfs==fs);

    return vfs;
}



template<class X> tribool
verify_primal_feasibility(const Matrix<X>& A, const Vector<X>& b, const array<VariableType>& vt)
{
    typedef Interval XX;
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    ARIADNE_ASSERT(b.size()==m);
    ARIADNE_ASSERT(vt.size()==n);

    const array<size_t> p=compute_p(vt);
    ARIADNE_LOG(9," p="<<p<<"\n");
    
    const Matrix<XX> B=compute_B<XX>(A,p);
    ARIADNE_LOG(9," B="<<midpoint(B)<<"\n   B*A="<<midpoint(prod(B,A))<<"\n");

    const Vector<XX> x=compute_x(A,b,p,B);
    ARIADNE_LOG(9," x="<<midpoint(x)<<"\n");
    
    tribool fs=true;
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if(!(0<x[j])) { 
            fs=indeterminate; 
            break; 
        }
    }
    if(fs==true) {
        return true; 
    }

    Vector<X> c(n);
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if(possibly(x[j]<=0)) { c[j]=-1; }
    }
    ARIADNE_LOG(9," c="<<midpoint(c)<<"\n");

    const Vector<XX> y=compute_y(c,p,B);
    ARIADNE_LOG(9," y="<<midpoint(y)<<"\n");

    const Vector<XX> z=compute_z(A,c,p,y);
    ARIADNE_LOG(9," z="<<midpoint(z)<<"\n");

    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        if(not(z[j]<0)) { return indeterminate; }
    }

    return false;
}
   



template<class X> tribool
verify_dual_feasibility(const Matrix<X>& A, const Vector<X>& c, const array<VariableType>& vt)
{
    ARIADNE_LOG(9,"verify_dual_feasibility(Matrix A, Vector c, VariableTypeArray vt):\n");
    typedef Interval XX;
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_ASSERT(c.size()==n);
    ARIADNE_ASSERT(vt.size()==n);
    array<size_t> p=compute_p(vt);
    ARIADNE_LOG(9," p="<<p<<"\n");
    Matrix<XX> B=compute_B<XX>(A,p);
    ARIADNE_LOG(9," B="<<midpoint(B)<<"\n");
    Vector<XX> y=compute_y(c,p,B);
    ARIADNE_LOG(9," y="<<midpoint(y)<<"\n");
    
    // Compute z=c-yA
    Vector<XX> z=compute_z(A,c,p,y);
    ARIADNE_LOG(9," z="<<midpoint(z)<<"\n");
    
    // Check feasibility of y
    tribool fs=true;
    for(size_t k=m; k!=n; ++k) {
        if(possibly(z[k]<=0)) { fs=false; break; }
    }
    if(fs==true) { return fs; }

    // Compute x_N
    Vector<XX> x(n);
    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        if(definitely(z[j]<0)) { x[j]=+1; }
    }
    // Compute w=-A_Nx_N
    Vector<XX> w(m);
    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        for(size_t i=0; i!=m; ++i) {
            w[i]-=A[i][j]*x[j];
        }
    }
    ARIADNE_LOG(9," w="<<midpoint(w)<<"\n");

    // Compute x_B=A_B^{-1}*w
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        for(size_t i=0; i!=m; ++i) {
            x[j]-=B[j][i]*w[i];
        }
    }
    
    ARIADNE_LOG(9," x="<<x<<"\n");
    ARIADNE_LOG(9," cx="<<dot(Vector<XX>(c),x)<<"\n");

    if(possibly(dot(Vector<XX>(c),x)>=0)) { return indeterminate; }
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if(possibly(x[j]<=0)) { return indeterminate; }
    }
    return false;
}

template<class X> tribool
verify_constrained_feasibility(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, const array<VariableType>& vt)
{
    ARIADNE_LOG(4,"verify_constrained_feasibility(Matrix A, Vector b, Vector l, Vector u, VariableTypeArray vt)\n");
    ARIADNE_LOG(5," A="<<A<<"\n b="<<b<<"\n l="<<l<<"\n u="<<u<<"\n t="<<vt<<"\n");

    typedef Interval XX;
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_ASSERT(b.size()==m);
    ARIADNE_ASSERT(l.size()==n);
    ARIADNE_ASSERT(u.size()==n);
    ARIADNE_ASSERT(vt.size()==n);

    const array<size_t> p=compute_p(vt);
    ARIADNE_LOG(9," p="<<p<<"\n");
    
    const Matrix<XX> B=compute_B<XX>(A,p);

    ARIADNE_LOG(9," B="<<midpoint(B)<<"\n   B*A="<<midpoint(prod(B,A))<<"\n");

    const Vector<XX> x=compute_x(A,b,l,u,vt,p,B);

    ARIADNE_LOG(9," x="<<midpoint(x)<<"\n");
    
    tribool fs=true;
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if(!(l[j]<x[j] && x[j]<u[j])) { 
            fs=indeterminate; 
            break; 
        }
    }
    if(fs==true) {
        return true; 
    }

    Vector<X> c(n);
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if(possibly(x[j]<=l[j])) { c[j]=-1; }
        if(possibly(x[j]>=u[j])) { c[j]=+1; }
    }
    ARIADNE_LOG(9," c="<<midpoint(c)<<"\n");

    const Vector<XX> y=compute_y(c,p,B);
    ARIADNE_LOG(9," y="<<midpoint(y)<<"\n");

    const Vector<XX> z=compute_z(A,c,p,y);
    ARIADNE_LOG(9," z="<<midpoint(z)<<"\n");

    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        if(vt[j]==LOWER && not(z[j]<0)) { return indeterminate; }
        if(vt[j]==UPPER && not(z[j]>0)) { return indeterminate; }
    }

    return false;
}
    
                

    
    



template<class X>
Vector<X> optimize(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_ASSERT(b.size()==m);
    ARIADNE_ASSERT(c.size()==n);
    ARIADNE_ASSERT(l.size()==n);
    ARIADNE_ASSERT(u.size()==n);

    array<VariableType> vt(n);
    array<size_t> p(n);
    Matrix<X> B(m,m);
    Vector<X> x(n);

    make_lpair(p,B)=compute_basis(A);
    for(size_t k=0; k!=m; ++k) { vt[p[k]]=BASIS; }
    for(size_t k=m; k!=n; ++k) { vt[p[k]]=LOWER; }
    x=compute_x(A,b,l,u,vt,p,B);
    feasible(A,b,l,u,vt,p,B,x);

    bool done=false;
    while(not done) {
        done=lpstep(A,b,c,l,u,vt,p,B,x);
    }
    
    return x;


}
    
template<class X>
tribool 
constrained_feasible_by_enumeration(const Matrix<X>& A, const Vector<X>& b, 
                                       const Vector<X>& l, const Vector<X>& u)
{
    bool fs=false;
    size_t m=A.row_size();
    size_t n=A.column_size();
    array<size_t> pb(m);
    for(size_t i=0; i!=m; ++i) { pb[i]=i; }
    do {
        array<size_t> p=extend_p(pb,n);  
        ARIADNE_LOG(9,"pb="<<pb<<" p="<<p<<"\n");
        Matrix<X> B=compute_B<X>(A,p);
        array<VariableType> vt(n);
        for(size_t k=0; k!=m; ++k) { vt[p[k]]=BASIS; } 
        for(size_t k=m; k!=n; ++k) { vt[p[k]]=LOWER; } 
        for(size_t cnt=0; cnt!=1u<<(n-m); ++cnt) {
            for(size_t k=m; k!=n; ++k) {
                if(vt[p[k]]==LOWER) { vt[p[k]]=UPPER; break; }
                else { vt[p[k]]=LOWER; }
            }
            ARIADNE_LOG(9," vt="<<vt);
            // Compute position x
            Vector<X> x(n);
            for(size_t k=m; k!=n; ++k) {
                x[p[k]]=(vt[p[k]]==LOWER ? l[p[k]] : u[p[k]]);
            }
            Vector<X> w(m);
            for(size_t k=m; k!=n; ++k) {
                for(size_t i=0; i!=m; ++i) {
                    w[i]+=A[i][p[k]]*x[p[k]];
                }
            }
            for(size_t k=0; k!=m; ++k) {
                for(size_t i=0; i!=m; ++i) {
                    x[p[k]]+=B[k][i]*w[i];
                }
            }
            ARIADNE_LOG(9," x="<<x);
            bool fsb=true;
            for(size_t i=0; i!=n; ++i) {
                if(l[i]>x[i] || x[i]>u[i]) {
                    fsb=false; 
                    break;
                }
            }
            fs = fs|fsb;
            ARIADNE_LOG(9," fsb="<<std::boolalpha<<fsb<<"\n");
        }
        
        // Increment basis elements 
        if(pb[m-1]+1<n) {
            ++pb[m-1];
        } else {
            for(size_t i=1; i!=m; ++i) {
                if(pb[m-i-1]+1<pb[m-i]) {
                    ++pb[m-i-1];
                    for(size_t j=m-i; j!=m; ++j) {
                        pb[m-j]=pb[m-j-1]+1;
                    }
                }
            }
        }
    } while(pb[0]!=(n-m));
    ARIADNE_LOG(9,"fs="<<std::boolalpha<<fs<<"\n");
    return fs;
}

#define ARIADNE_INSTANTIATE_LINEAR_PROGRAMMING(Class) \
    template tribool primal_feasible<Class>(const Matrix<Class>&, const Vector<Class>&); \
    template tribool dual_feasible<Class>(const Matrix<Class>&, const Vector<Class>&); \
    template tribool constrained_feasible<Class>(const Matrix<Class>&, const Vector<Class>&, const Vector<Class>&, const Vector<Class>&); \
                                                                        \
    template tribool constrained_feasible_by_enumeration<Class>(const Matrix<Class>&, const Vector<Class>&, const Vector<Class>&, const Vector<Class>&); \
                                                                        \
    template tribool verify_primal_feasibility<Class>(const Matrix<Class>&, const Vector<Class>&, const array<VariableType>&); \
    template tribool verify_dual_feasibility<Class>(const Matrix<Class>&, const Vector<Class>&, const array<VariableType>&); \
    template tribool verify_constrained_feasibility<Class>(const Matrix<Class>&, const Vector<Class>&, const Vector<Class>&, const Vector<Class>&, const array<VariableType>&); \
                                                                        \
    template std::pair< array<size_t>, Matrix<Class> > compute_basis(const Matrix<Class>&); \


ARIADNE_INSTANTIATE_LINEAR_PROGRAMMING(Float);

#ifdef HAVE_RATIONAL
ARIADNE_INSTANTIATE_LINEAR_PROGRAMMING(Rational);
#endif // HAVE_RATIONAL

#undef ARIADNE_INSTANTIATE_LINEAR_PROGRAMMING


} // namespace Ariadne

