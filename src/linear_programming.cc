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

namespace Ariadne {

std::ostream& operator<<(std::ostream& os, VariableType t) {
    return os << (t==INFEASIBLE ? 'I' : t==FEASIBLE ? 'F'  : t==BASIS ? 'B' : t==LOWER ? 'L' : 'U');  
}

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
array<size_t> extend(const array<size_t>& p, const size_t n)
{
    const size_t m=p.size();
    array<size_t> q(n);
    for(size_t j=0; j!=n; ++j) {
        q[j]=n;
    }
    for(size_t k=0; k!=m; ++k) {
        assert(p[k]<n);
        assert(q[p[k]]==n);
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
        assert(p[i]==r[i]);
    }
    return r;
}
      
// Compute a basis (p_1,\ldots,p_m) for the matrix A
// Assume 
template<class X>
pair< array<size_t>, Matrix<X> >
compute_basis(const Matrix<X>& A)
{
    //std::cerr<<"compute_basis(A) with A="<<A<<std::endl;
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
            if(j==n) { ARIADNE_THROW(std::runtime_error,"compute_basis"," matrix A="<<A<<" is singular"); }
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
    

    
template<class X>
void check(const Matrix<X>& A, const array<size_t>& p, const Matrix<X>& B) 
{
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
    std::cerr<<"        p_B="<<p_B<<" B="<<B<<" A_B="<<A_B<<" B*A_B-I="<<Z<<std::endl;
    for(size_t i=0; i!=m; ++i) { Z[i][i]-=1; }
    assert(norm(Z)<1e-5);
}


template<class X>
void check(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, const array<VariableType>& vt, const array<size_t>& p, const Matrix<X>& B, const Vector<X>& x) 
{
    std::cerr<<"        Checking consistency of B and x\n";
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
    std::cerr<<"          p_B="<<p_B<<" B="<<B<<" A_B="<<A_B<<" B*A_B="<<Z<<std::endl;
    for(size_t i=0; i!=m; ++i) { Z[i][i]-=1; }

    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        assert(vt[j]==LOWER || vt[j]==UPPER);
        X xj = (vt[j]==LOWER ? l[j] : u[j]);
        assert(x[j]==xj);
    }
    Vector<X> z=prod(A,x);
    std::cerr<<"          A="<<A<<" x="<<x<<" b="<<b<<" Ax="<<z<<std::endl;
    z-=b;
    assert(norm(z)<1e-5);
}


template<class X>
void verify_feasibility(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, const array<VariableType> vt)
{
    std::cerr<<"\nChecking feasibility with\n  A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<"\n  t="<<vt<<"\n"<<std::endl;
    typedef X XX;
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    // Compute basic variable indices
    array<uint> p(m); size_t j=0; for(size_t i=0; i!=m; ++i) { while(vt[j]!=BASIS) { ++j; } p[i]=j; } 
    
    // Compute inverse matrix
    Matrix<XX> A_B(m,m);
    for(size_t i=0; i!=m; ++i) { for(size_t j=0; j!=m; ++j) { A_B[i][j]=A[i][p[j]]; } }
    Matrix<XX> B=A_B.inverse();
    
    // Compute non-basic variables
    Vector<XX> x; 
    for(size_t i=0; i!=n; ++i) { if(vt[i]==LOWER) { x[i]=l[i]; } else if(vt[i]==UPPER) { x[i]=u[i]; } else { x[i]=0; } }

    Vector<XX> x_B = b-prod(B,prod(A,x));
    for(size_t i=0; i!=m; ++i) { assert(x_B[i]>l[p[i]]); assert(x_B[i]<u[p[i]]); }
}
    
    
   
    

template<class X>
void verify_infeasibility(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u,
                          const array<VariableType> t, const Vector<X>& y)
{
    // Check that (b-A_L l_L - A_U u_U).y > 0;  y A_L <= 0; y A_U >=0 
    std::cerr<<"\nChecking infeasibility with\n  A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<"\n  t="<<t<<" y="<<y<<"\n"<<std::endl;
    //const size_t m=A.row_size();
    const size_t n=A.column_size();
    Vector<X> x(n);
    for(uint i=0; i!=n; ++i) {
        x[i] = (t[i]==LOWER ? l[i] : t[i]==UPPER ? u[i] : X(0));
    }
    Vector<X> w=b-prod(A,x);
    X d=dot(w,y);
    Vector<X> yA=prod(y,A);
    
    std::cerr<<"  x="<<x<<" w="<<w<<" y="<<y<<" y^Tw="<<d<<" t="<<t<<" yA="<<yA<<"\n";
    assert(d>0);
    for(uint i=0; i!=n; ++i) {
        if(t[i]==LOWER) { if(yA[i]>=0) { std::cerr<<"  t["<<i<<"]=L, yA["<<i<<"]="<<yA[i]<<std::endl; } assert(yA[i]<0); }
        if(t[i]==UPPER) { if(yA[i]<=0) { std::cerr<<"  t["<<i<<"]=U, yA["<<i<<"]="<<yA[i]<<std::endl; } assert(yA[i]>0); }
    }
}



template<class X>
void verify_infeasibility(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& y)
{
    X z(0);
    Vector<X> yA=prod(y,A);
    X d=dot(y,b);

    if(d<=0) { std::cerr<< "dot(y,b)="<<d<<" is not strictly positive."; assert(false); }
    if(!(yA<=z)) { std::cerr<< "yA="<<prod(y,A)<<" is not negative."; assert(false); }
    
}



template<class X> inline X inf();
template<> inline Float inf<Float>() { return std::numeric_limits<double>::infinity(); }
//template<> inline Rational inf<Rational>() { return Rational(1000000); }
template<> inline Rational inf<Rational>() { return Rational(1,-0); }



template<class X>
tribool feasible(const Matrix<X>& A, const Vector<X>& b) {
    throw std::runtime_error("Not Implemented");
}


template<class X>
Matrix<X> compute_B(const Matrix<X>& A, const array<size_t>& p)
{
    const size_t m=A.row_size();
    Matrix<X> A_B(m,m);
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        for(size_t i=0; i!=m; ++i) {
            A_B[i][k]=A[i][j];
        }
    }
    
    Matrix<X> B=inverse(A_B);
    
    return B;
}



template<class X>
Vector<X> compute_x(const Matrix<X>& A, const Vector<X>& b, array<size_t>& p, Matrix<X>& B)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Vector<X> x(n);

    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        x[j]=0;
        for(size_t i=0; i!=m; ++i) {
            x[j]+=B[k][i]*b[i];
        }
    }

    return x;
}


template<class X>
Vector<X>
compute_x(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, const array<size_t>& p, const Matrix<X>& B)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    
    Vector<X> w(m);
    Vector<X> x(n);
    
    // Compute x_N
    for(size_t j=0; j!=n; ++j) {
        if(vt[j]==LOWER) { x[j]=l[j]; }
        else if(vt[j]==UPPER) { x[j]=u[j]; }
        else { x[j]=0; }
    }
    //std::cerr<<"  x_N="<<x<<std::flush;

    // Compute w=b-A_N x_N
    for(size_t i=0; i!=m; ++i) {
        w[i]=b[i];
        for(size_t k=m; k!=n; ++k) {
            size_t j=p[k];
            w[i]-=A[i][j]*x[j];
        }
    }
    //std::cerr<<" w="<<w<<std::flush;

    // Compute x_B=B w
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        x[j]=0;
        for(size_t i=0; i!=m; ++i) {
            x[j]+=B[k][i]*w[i];
        }
    }

    //std::cerr<<" x="<<x<<"\n"<<std::flush;

    Vector<X> Axmb=prod(A,x)-b;
    assert(norm(Axmb)<0.00001);
    return x;
}


template<class X>
std::pair<Vector<X>,Vector<X> >
compute_wx(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, const array<size_t>& p, const Matrix<X>& B)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    
    Vector<X> w(m);
    Vector<X> x(n);

    // Compute x_N
    for(size_t j=0; j!=n; ++j) {
        if(vt[j]==LOWER) { x[j]=l[j]; }
        else if(vt[j]==UPPER) { x[j]=u[j]; }
        else { x[j]=0; }
    }
    //std::cerr<<"  x_N="<<x<<std::flush;

    // Compute w=b-A_N x_N
    for(size_t i=0; i!=m; ++i) {
        w[i]=b[i];
        for(size_t k=m; k!=n; ++k) {
            size_t j=p[k];
            w[i]-=A[i][j]*x[j];
        }
    }
    //std::cerr<<" w="<<w<<std::flush;

    // Compute x_B=B w
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        x[j]=0;
        for(size_t i=0; i!=m; ++i) {
            x[j]+=B[k][i]*w[i];
        }
    }

    //std::cerr<<" x="<<x<<"\n"<<std::flush;

    Vector<X> Axmb=prod(A,x)-b;
    assert(norm(Axmb)<0.00001);
    return make_pair(w,x);
}



template<class X>
bool lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, array<size_t>& p, Matrix<X>& B, Vector<X>& x)
{
    std::cerr<<"lpstep(A,b,c,p,B,x)\n  A="<<A<<" b="<<b<<" c="<<c<<" p="<<p<<" B="<<B<<" x="<<x<<std::endl;
    assert(B==compute_B(A,p));
    assert(x==compute_x(A,b,p,B));
   
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    static const X plus_inf=inf<X>();
    Vector<X> y(m);
    Vector<X> z(n);
    
    // Compute y=c_B B
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        for(uint i=0; i!=m; ++i) {
            y[i]+=c[j]*B[k][i];
        }
    }

     
    // Compute reduced costs z_N=c_N-yA_N
    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        z[j]=c[j];
        for(uint i=0; i!=m; ++i) {
            z[j]-=y[i]*A[i][j];
        }
    }

    std::cerr<<"    x="<<x<<" y="<<y<<" z="<<z<<"\n";

    // Compute variable p[s] to enter the basis
    size_t s=n;
    for(size_t k=m; k!=n; ++k) {
        if(z[p[k]]<0) { s=k; }
    }
             
    if(s==n) { std::cerr<<"No improvement possible.\n"; return true; }
    
    // Compute direction d in which to move the current basic variables
    Vector<X> As(m); for(uint i=0; i!=m; ++i) { As[i]=A[i][p[s]]; }
    Vector<X> d=prod(B,As);
    
    // Compute distance t along d in which to move,
    // and the variable p[r] to leave the basis
    X t=inf<X>();
    size_t r=m;
    for(size_t i=0; i!=m; ++i) {
        if(d[i]>0) { X ti=(x[p[i]])/d[i]; if(r==m || ti<t) { t=ti; r=i; } }
    }

    if(r==m) { ARIADNE_THROW(std::runtime_error,"feasible","Unbounded linear program"); }
    
    // Compute new value v of x[p[s]] and make sure that it's still feasible
    
    
    std::cerr<<"  Swapping basic variable "<<p[r]<<" in position "<<r<<" with non-basic variable "<<p[s]<<" in position "<<s<<"\n";

    std::cerr<<"    r="<<r<<" s="<<s<<" t="<<t<<" d="<<d<<std::endl;

    // Update matrix of inverses
    X xr=x[p[r]];
    X dr=d[r];
    X drr=1/dr;
    Vector<X> e(m); e[r]=1;
    Vector<X> Br(m); for(uint j=0; j!=m; ++j) { Br[j]=B[r][j]; }
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=m; ++j) {
            B[i][j]-=(d[i]-e[i])*Br[j]*drr;
        }
    }

    // Update state vector x
    for(size_t i=0; i!=m; ++i) {
        x[p[i]]-=t*d[i];
    }
    x[p[r]] = 0.0;
    x[p[s]] = t;

    // Update pivots
    std::swap(p[r],p[s]);

    std::cerr<<"      B ="<<B<<" x ="<<x<<" p="<<p<<"\n";
    
    Matrix<X> BB=compute_B(A,p);
    Vector<X> xx=compute_x(A,b,p,B);
    
    std::cerr<<"      BB="<<BB<<" xx="<<xx<<"\n";
    check(A,p,B);

    return false;
}


template<class X>
bool lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x)
{
    std::cerr<<"  lpstep(A,b,c,l,u,vt,p,V,x)\n    A="<<A<<" b="<<b<<" c="<<c<<"\n    p="<<p<<" B="<<B<<"\n    vt="<<vt<<" l="<<l<<" x="<<x<<" u="<<u<<"\n";

    const X infinity=inf<X>();
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    static const X plus_inf=inf<X>();
    Vector<X> y(m);
    Vector<X> z(n);
    
    // Compute y=c_B B
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        for(uint i=0; i!=m; ++i) {
            y[i]+=c[j]*B[k][i];
        }
    }

     
    // Compute reduced costs z_N=c_N-yA_N
    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        z[j]=c[j];
        for(uint i=0; i!=m; ++i) {
            z[j]-=y[i]*A[i][j];
        }
    }

    std::cerr<<"      y="<<y<<" z="<<z<<"\n";

    // Compute variable p[s] to enter the basis
    size_t s=n;
    for(size_t k=m; k!=n; ++k) {
        if(z[p[k]]<0) { s=k; }
    }
             
    if(s==n) { std::cerr<<"No improvement possible.\n"; return true; }
    
    // Compute direction d in which to move the current basic variables
    X ds=(vt[p[s]]==LOWER) ? +1 : -1;
    Vector<X> As(m); for(uint i=0; i!=m; ++i) { As[i]=A[i][p[s]]; }
    Vector<X> d=(-ds)*prod(B,As);
   

    // Compute distance t along d in which to move,
    // and the variable p[r] to leave the basis
    // The bounds on t are given by l <= x +/- t * d <= u 
    X ts=u[p[s]]-l[p[s]];
    std::cerr<<"      s="<<s<<" d="<<d<<" ds="<<ds<<" ts="<<ts<<std::flush;

    size_t r=n;
    X tr=infinity;
    for(size_t i=0; i!=m; ++i) {
        size_t j=p[i];
        if(d[i]<0 && x[j]>=l[j]) { 
            X ti=(l[j]-x[j])/d[i]; if(tr==infinity || ti<tr) { tr=ti; r=i; } 
            //std::cerr<<"          i="<<i<<" j="<<j<<" xj="<<x[j]<<" lj="<<l[j]<<" di="<<d[i]<<" ti="<<ti<<" nxj="<<x[j]+d[i]*ti<<"\n"; 
        } else if( d[i]>0 && x[j]<=u[j] && u[j] != inf<X>() ) { 
            X ti=(u[j]-x[j])/d[i]; if(tr==infinity || ti<tr) { tr=ti; r=i; } 
            //std::cerr<<"          i="<<i<<" j="<<j<<" xj="<<x[j]<<" uj="<<u[j]<<" di="<<d[i]<<" ti="<<ti<<" nxj="<<x[j]+d[i]*ti<<"\n"; 
        }
    }

    std::cerr<<"       r="<<r<<" tr="<<tr<<" t="<<std::min(tr,ts)<<std::endl;

    if(ts==inf<X>() && tr==inf<X>()) { 
        ARIADNE_THROW(std::runtime_error,"feasible","Unbounded linear program"); 
    }

    X t; 
    if(ts<=tr) { r=s; t=ts; } else { t=tr; }

    if(r!=s) {
        std::cerr<<"   Swapping basic variable "<<p[r]<<" in position "<<r<<" with non-basic variable "<<p[s]<<" in position "<<s<<"\n";
    } else {
        std::cerr<<"   Changing non-basic variable "<<p[s]<<" in position "<<s<<" to type "<<vt[p[s]]<<"\n";
    }

    
    Vector<X> e(m); X dr,xr;
    if(r==s) {
        dr=1;
        // No change in basic variables. Constraint is due to bounds on x_s
        if(vt[p[s]]==LOWER) {
            vt[p[s]]=UPPER;
        } else {
            vt[p[s]]=LOWER;
        }
        xr=ds*(u[p[s]]-l[p[s]]);
    } else {
        e[r]=1;
        dr=d[r];
        xr=x[p[r]];
    }

    
    // Update matrix of inverses
    X drr=1/dr;
    Vector<X> Br(m); for(uint j=0; j!=m; ++j) { Br[j]=B[r][j]; }
    Vector<X> nd(m); for(uint i=0; i!=m; ++i) { nd[i]=(d[i]+e[i])*drr; }
    for(uint i=0; i!=m; ++i) { 
        for(uint j=0; j!=m; ++j) {
            B[i][j]-=nd[i]*Br[j];
        }
    }

    // Update state vector x to x-td
    for(size_t i=0; i!=m; ++i) {
        x[p[i]]+=d[i]*t;
    }
    x[p[s]] += ds*t;

    std::cerr<<"      B="<<B<<"  x="<<x<<"\n";
    
    // Update pivots and variable types
    std::swap(p[r],p[s]);
    for(size_t i=0; i!=m; ++i) { vt[p[i]]=BASIS; }
    vt[p[s]]=(dr>0 ? UPPER : LOWER);
    X xpr=x[p[r]];
    x[p[s]] = (dr>0 ? u[p[s]] : l[p[s]]);


    std::cerr<<"      vt="<<vt<<"  p="<<p<<"\n";
    
    Matrix<X> BB=compute_B(A,p);
    Vector<X> xx=compute_x(A,b,l,u,vt,p,BB);
    
    std::cerr<<"      BB="<<BB<<" xx="<<xx<<"\n";
    check(A,b,l,u,vt,p,B,x);

    return false;
    }



template<class X>
bool lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, array<size_t>& p, Matrix<X>& B)
{
    Vector<X> x=compute_x(A,b,l,u,vt,p,B);
    return lpstep(A,b,l,u,vt,p,B,x);
}



template<class X>
tribool _primal_feasible(const Matrix<X>& A, const Vector<X>& b, array<size_t>& p, Matrix<X>& B)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Vector<X> x=prod(B,b);
    Vector<X> c(n);

    bool infeasible=false;
    for(size_t j=0; j!=n; ++j) {
        if(x[j]<0) { c[j]=-1; infeasible=true; }
        else { c[j]=0; }
    }
    std::cerr << "    x="<<x<<" c="<<c<<"\n";

    while(infeasible) {
        infeasible=false;

        bool done=lpstep(A,b,c,p,B,x);
        std::cerr<<"  Done changing basis\n";
        std::cerr << "    p="<<p<<" B="<<B<<" x="<<x<<std::endl;
        
        if(done) {
            std::cerr<<"  Cannot put infeasible variables into basis.";
            return false;
        }

        for(size_t j=0; j!=n; ++j) {
            if(x[j]<0) { c[j]=-1; }
            else { c[j]=0; }
        }
        std::cerr << "\n    x="<<x<<" c="<<c<<"\n";
    }
 
    // Check solution
    for(size_t i=0; i!=n; ++i) {
        assert(x[i]>=0.0);
    }
    Vector<X> Ax=prod(A,x);
    for(size_t i=0; i!=m; ++i) {
         assert(Ax[i]==b[i]);
    }

    std::cerr<<"\nFeasible point x="<<x<<"\n Ax="<<Vector<X>(prod(A,x))<<" b="<<b<<std::endl;

    return true;
} 



template<class X>
tribool feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<VariableType>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Vector<X> c(n);

    bool infeasible=false;
    for(size_t j=0; j!=n; ++j) {
        if(x[j]<l[j]) { c[j]=-1; vt[j]=INFEASIBLE; infeasible=true; }
        else if(x[j]>u[j]) { c[j]=+1; vt[j]=INFEASIBLE; infeasible=true; }
        else { c[j]=0; if(vt[j]==BASIS) { vt[j]=FEASIBLE; } }
    }
    std::cerr << "    vt="<<vt<<" x="<<x<<" c="<<c<<"\n";

    while(infeasible) {
        infeasible=false;

        bool done=lpstep(A,b,c,l,u,vt,p,B,x);
        std::cerr<<"  Done changing basis\n";
        std::cerr << "    p="<<p<<" B="<<B<<std::endl;
        std::cerr << "    vt="<<vt<<" x="<<x<<std::endl;
        
        if(done) {
            std::cerr<<"  Cannot put infeasible variables into basis.";
            return false;
        }

        for(size_t j=0; j!=n; ++j) {
            if(vt[j]==LOWER) { assert(x[j]==l[j]); }
            if(vt[j]==UPPER) { assert(x[j]==u[j]); }
            if(x[j]<l[j]) { c[j]=+1; vt[j]=INFEASIBLE; infeasible=true; }
            else if(x[j]>u[j]) { c[j]=-1; vt[j]=INFEASIBLE; infeasible=true; }
            else { c[j]=0; if(vt[j]==INFEASIBLE || vt[j]==BASIS) { vt[j]=FEASIBLE; } }
        }
        std::cerr << "\n    vt="<<vt<<" x="<<x<<" c="<<c<<"\n";
    }
 
    // Check solution
    for(size_t i=0; i!=n; ++i) {
        assert(x[i]>=l[i]);
        assert(x[i]<=u[i]);
    }
    Vector<X> Ax=prod(A,x);
    for(size_t i=0; i!=m; ++i) {
         assert(Ax[i]==b[i]);
    }

    std::cerr<<"\nFeasible point x="<<x<<"; l="<<l<<" u="<<u<<"\n Ax="<<Vector<X>(prod(A,x))<<" b="<<b<<std::endl;

    return true;
} 



// Check for feasibility of Ax=b x>=0
template<class X>
tribool primal_feasible(const Matrix<X>& A, const Vector<X>& b)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_ASSERT(b.size()==n);
 
    array<size_t> p(n); // Array of basis indices followed by non-basis indices
    Matrix<X> B(m,m);
    Vector<X> x(n);
    

    std::cerr<<"\nprimal_feasible(A,b) with\nA="<<A<<" b="<<b<<" \n"<<std::endl;
    make_lpair(p,B)=compute_basis(A);

    std::cerr<<"    p="<<p<<" B="<<B<<"  (BA="<<Matrix<X>(prod(B,A))<<")\n";

    return _primal_feasible(A,b,p,B);
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
    Vector<X> c_B(m);
    
    bool done;
    do {
        lpstep(A,b,c,p,B,x);
        for(size_t k=0; k!=m; ++k) { c_B[k]=c[p[k]]; }
        y=prod(c_B,B);
        z=c-prod(y,A);
        
        done=true;
        for(size_t i=0; i!=n; ++i) {
            if(z[i]<0) { done=false; }
        }
        if(done) { return true; }
    } while(true);
    return false;
}


// Check for feasibility of yA<=c. The solution is infeasible if there exists x with cx<0, Ax=0 and x>=0. 
template<class X>
tribool dual_feasible(const Matrix<X>& A, const Vector<X>& c) 
{
    const size_t n=A.column_size();
    ARIADNE_ASSERT(c.size()==n);

    array<size_t> p;
    Matrix<X> B;
    make_lpair(p,B)=compute_basis(A);
    return _dual_feasible(A,c,p,B);
}


    
// Check for feasibility of yA<=c. The solution is infeasible if there exists x with cx<0, Ax=0 and x>=0. The values in array p are the initial pivot elements
template<class X>
tribool dual_feasible(const Matrix<X>& A, const Vector<X>& c, array<size_t>& p) 
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_ASSERT(c.size()==m);
    ARIADNE_ASSERT(p.size()==n || p.size()==m);
 
    if(p.size()==m) { p=extend(p,n); }
    Matrix<X> B=compute_b(A,p);
    
    return _dual_feasible(A,c,p,B);
}

// Check for feasibility of yA<=c. The solution is infeasible if there exists x with cx<0, Ax=0 and x>=0. The values in array p are the initial pivot elements, and the value of B is the initial A_B inverse.
template<class X>
tribool dual_feasible(const Matrix<X>& A, const Vector<X>& c, array<size_t>& p, Matrix<X>& B) 
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_ASSERT(p.size()==n || p.size()==m);
    ARIADNE_ASSERT(B.row_size()==m);
    ARIADNE_ASSERT(B.column_size()==m);
 
    if(p.size()==m) { p=extend(p,n); }

    return _dual_feasible(A,c,p,B);
}

    

// Check for feasibility of Ax=b l<=b<=u
template<class X>
tribool constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_ASSERT(b.size()==m);
    ARIADNE_ASSERT(l.size()==n);
    ARIADNE_ASSERT(u.size()==n);
    Vector<X> c(n);

    array<size_t> p(n); // Array of basis indices followed by non-basis indices
    Matrix<X> B(m,m);
    Vector<X> x(n);
    
    array<VariableType> vt(n);

    std::cerr<<"\nfeasible(A,b,l,u) with\nA="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<"\n"<<std::endl;
    make_lpair(p,B)=compute_basis(A);

    std::cerr<<"    p="<<p<<" B="<<B<<"  (BA="<<Matrix<X>(prod(B,A))<<")\n";

    for(size_t k=0; k!=m; ++k) { vt[p[k]]=BASIS; }
    for(size_t k=m; k!=n; ++k) { vt[p[k]]=LOWER; }

    x=compute_x(A,b,l,u,vt,p,B);

    feasible(A,b,l,u,vt,p,B,x);

    return true;
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
    

template tribool primal_feasible<Rational>(const Matrix<Rational>&, const Vector<Rational>&);
template tribool dual_feasible<Rational>(const Matrix<Rational>&, const Vector<Rational>&);
template tribool constrained_feasible<Rational>(const Matrix<Rational>&, const Vector<Rational>&, const Vector<Rational>&, const Vector<Rational>&);

//template void verify_infeasibility(const Matrix<Rational>& A, const Vector<Rational>& b, const Vector<Rational>& l, const Vector<Rational>& u, const array<size_t>& p, const Vector<Rational>& y);

template void verify_infeasibility(const Matrix<Rational>& A, const Vector<Rational>& b, const Vector<Rational>& y);

template std::pair< array<size_t>, Matrix<Rational> > compute_basis(const Matrix<Rational>&);

template bool lpstep(const Matrix<Rational>& A, const Vector<Rational>& b, const Vector<Rational>& c, array<size_t>& p, Matrix<Rational>& B, Vector<Rational>& x);

} // namespace Ariadne

