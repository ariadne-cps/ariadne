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
#include "affine.h"
#include "linear_programming.h"

#include "macros.h"
#include "logging.h"

#include <boost/numeric/ublas/matrix_sparse.hpp>

static const int verbosity=0;

namespace Ariadne {

class DegenerateFeasibilityProblemException : public std::exception { };

StandardLinearProgram::StandardLinearProgram(uint m, uint n)
    : A(m,n), b(m), c(n), x(n), y(m), z(n) { }


StandardLinearProgram
feasibility_problem(const Matrix<Float>& B, const Vector<Float>& bl, const Vector<Float>& bu,
                    const Matrix<Float>& C, const Vector<Float>& c,
                    const Vector<Float>& dl, const Vector<Float>& du)
{
    const uint nd=dl.size();
    const uint nc=C.column_size();
    const uint nb=B.column_size();
    const uint ns=2*nb+nc+2*nd;

    ARIADNE_ASSERT(B.row_size()==nd);
    ARIADNE_ASSERT(C.row_size()==nd);
    ARIADNE_ASSERT(bl.size()==nb);
    ARIADNE_ASSERT(bu.size()==nb);
    ARIADNE_ASSERT(c.size()==nc);
    ARIADNE_ASSERT(dl.size()==nd);
    ARIADNE_ASSERT(du.size()==nd);

    StandardLinearProgram lp(nd+1,ns);
    Matrix<Float>& AA=lp.A;
    Vector<Float>& bb=lp.b;
    Vector<Float>& cc=lp.c;
    Vector<Float>& xx=lp.x;
    Vector<Float>& yy=lp.y;
    Vector<Float>& zz=lp.z;

    Vector<Float> y(nd), yB(nb), yC(nc);
    y=(dl+du)/2; yB=prod(y,B); yC=prod(y,C);
    Vector<Float> e(nc,1.0), Ce(prod(C,e));

    for(uint i=0; i!=nd; ++i) {
        for(uint j=0; j!=nb; ++j) {
            AA[i][j]=B[i][j];
            AA[i][j+nb]=-B[i][j];
        }
        for(uint j=0; j!=nc; ++j) {
            AA[i][j+2*nb]=C[i][j];
        }
        AA[i][i+2*nb+nc]=1.0;
        AA[i][i+2*nb+nc+nd]=-1.0;
    }
    for(uint jj=0; jj!=2*nb+nc; ++jj) {
        AA[nd][jj]=1.0;
    }

    ARIADNE_LOG(2,"AA="<<AA<<"\n");

    // Set bb to be all zeros apart from last element
    bb[nd]=1.0;

    // Set cc to be upper constraints
    for(uint j=0; j!=nb; ++j) {
        cc[j]=bu[j];
        cc[j+nb]=-bl[j];
    }

    for(uint j=0; j!=nc; ++j) {
        cc[j+2*nb]=c[j];
    }

    for(uint j=0; j!=nd; ++j) {
        cc[j+2*nb+nc]=du[j];
        cc[j+2*nb+nc+nd]=-dl[j];
    }

    ARIADNE_LOG(2,"bb="<<bb<<" cc="<<cc<<"\n");

    // Set xx to be constant vector plus whatever is needed to make constraints satisfied
    for(uint jj=0; jj!=ns; ++jj) {
        xx[jj]=1.0/(2*nb+nc);
    }
    ARIADNE_LOG(2,"C="<<C<<" Ce="<<Ce<<" xx="<<xx<<" AA*xx="<<prod(AA,xx)<<"\n");
    for(uint i=0; i!=nd; ++i) {
        if(Ce[i]>0.0) { xx[i+2*nb+nc+nd]+=Ce[i]/(2*nb+nc); }
        else { xx[i+2*nb+nc]-=Ce[i]/(2*nb+nc); }
    }
    ARIADNE_LOG(2,"Ce="<<Ce<<" xx="<<xx<<" AA*xx="<<prod(AA,xx)<<"\n");

    // Set yy to be midpoint of domain
    for(uint i=0; i!=nd; ++i) {
        yy[i]=y[i];
    }

    // Set zz and t=yy[nd] to be midpoint of domain
    for(uint j=0; j!=nb; ++j) {
        zz[j]=bu[j]-yB[j];
        zz[j+nb]=yB[j]-bl[j];
    }
    for(uint j=0; j!=nc; ++j) {
        zz[j+2*nb]=c[j]-yC[j];
    }
    for(uint j=0; j!=nd; ++j) {
        zz[j+2*nb+nc]=du[j]-y[j];
        zz[j+2*nb+nc+nd]=y[j]-dl[j];
    }
    Float t=+inf<Float>();
    for(uint jj=0; jj!=2*nb+nc; ++jj) {
        t=min(t,zz[jj]);
    }
    t-=1.0/ns;
    yy[nd]=t;
    for(uint jj=0; jj!=2*nb+nc; ++jj) {
        zz[jj]-=t;
    }

    ARIADNE_LOG(2,"xx="<<xx<<"\nyy="<<yy<<" zz="<<zz<<"\n");

    ARIADNE_LOG(2,"AA*xx-bb="<<prod(AA,xx)-bb<<"\n");
    ARIADNE_LOG(2,"yy*AA+zz-cc="<<prod(yy,AA)+zz-cc<<"\n");

    return lp;
}



StandardLinearProgram
feasibility_problem(const Vector< Affine<Float> >& f, const Vector<Interval>& b,
                    const Vector< Affine<Float> >& g, const Vector<Float>& c,
                    const Vector<Interval>& d)
{
    ARIADNE_NOT_IMPLEMENTED;
}




struct SingularLinearProgram : std::runtime_error {
    SingularLinearProgram(const std::string& what)
        : std::runtime_error(what) { };
};

struct UnboundedLinearProgram : std::runtime_error {
    UnboundedLinearProgram(const std::string& what)
        : std::runtime_error(what) { }
};

std::ostream& operator<<(std::ostream& os, Slackness t) {
    return os << (t==BASIS ? 'B' : t==LOWER ? 'L' : t==UPPER ? 'U' : '?');
}


// Compute A A^T
template<class X> inline Matrix<X> sprod(const Matrix<X>& A) {
    Matrix<X> S(A.row_size(),A.row_size());
    // Compute lower triangle
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint j=0; j!=i+1u; ++j) {
            S[i][j]=0.0;
            for(uint k=0; k!=A.column_size(); ++k) {
                S[i][j]+=A[i][k]*A[j][k];
            }
        }
    }
    for(uint j=1; j!=A.row_size(); ++j) {
        for(uint i=0; i!=j; ++i) {
            S[i][j]=S[j][i];
        }
    }
    return S;
}


template<class X> inline Vector<X> operator*(const Matrix<X>& A, const Vector<X>& x) {
    return prod(A,x); }
template<class X> inline Vector<X> operator*(const Vector<X>& y, const Matrix<X>& A) {
    return prod(y,A); }

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
array<size_t>
extend_p(const array<size_t>& p, const size_t n)
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





// Check that the basic variable array p is consistent with the variable type array vt.
// There are two cases; p just lists the basic variables, or p lists all variables
// Returns the number of basic variables
size_t
consistency_check(const array<Slackness>& vt, const array<size_t>& p)
{
    if(p.size()!=vt.size()) {
        const size_t m=p.size();
        const size_t n=vt.size();
        ARIADNE_ASSERT(m<n);
        for(size_t i=0; i!=m; ++i) {
            ARIADNE_ASSERT_MSG(p[i]<n && vt[p[i]]==BASIS, "vt="<<vt<<" p="<<p);
        }
        return m;
    } else {
        const size_t n=vt.size();
        size_t m=n;
        for(size_t i=0; i!=m; ++i) {
            ARIADNE_ASSERT_MSG(p[i]<n, "vt="<<vt<<" p="<<p);
            if(vt[p[i]]!=BASIS) { m=n; break; }
        }
        for(size_t i=n; i!=n; ++i) {
            ARIADNE_ASSERT_MSG(p[i]<n && vt[p[i]]==BASIS, "vt="<<vt<<" p="<<p);
        }
        return m;
    }
}


// Check that the matrix B is the inverse of the matrix A_B with columns p[0],...,p[m-1] of A.
template<class X>
void
consistency_check(const Matrix<X>& A, const array<size_t>& p, const Matrix<X>& B)
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
void
consistency_check(const Matrix<X>& A, const Vector<X>& b, const array<size_t>& p, const Matrix<X>& B, const Vector<X>& x)
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
void
consistency_check(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, const array<Slackness>& vt, const array<size_t>& p, const Matrix<X>& B, const Vector<X>& x)
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

    Matrix<X> I=prod(B,A_B);
    ARIADNE_LOG(9,"          p_B="<<p_B<<" B="<<B<<" A_B="<<A_B<<" B*A_B="<<I<<"\n");
    Matrix<X> Z=I;
    for(size_t i=0; i!=m; ++i) { Z[i][i]-=1; }
    ARIADNE_ASSERT_MSG(norm(Z)<=1e-5,"vt="<<vt<<" p_B="<<p_B<<" B="<<B<<" A_B="<<A_B<<" B*A_B="<<I<<"\n");

    Vector<X> Ax=prod(A,x);
    ARIADNE_LOG(9,"          A="<<A<<" x="<<x<<" b="<<b<<" Ax="<<Ax<<"\n");

    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        ARIADNE_ASSERT(vt[j]==LOWER || vt[j]==UPPER);
        X xj = (vt[j]==LOWER ? l[j] : u[j]);
        ARIADNE_ASSERT_MSG(x[j]==xj,"A="<<A<<", b="<<b<<", l="<<l<<", u="<<u<<", vt="<<vt<<", p="<<p<<", x="<<x<<", Ax="<<Ax);
    }
    Vector<X> Axmb = Ax-b;
    ARIADNE_ASSERT_MSG(norm(Axmb)<1e-5,"A="<<A<<", b="<<b<<", l="<<l<<", u="<<u<<", vt="<<vt<<", p="<<p<<", x="<<x<<", Ax="<<Ax);
}




// Compute the cost function for a feasibility step given lower and upper bounds and the values of x.
template<class X>
Vector<X>
compute_c(const Vector<X>& l, const Vector<X>& u, const array<size_t>& p, const Vector<X>& x, size_t m) {
    const size_t n=x.size();
    Vector<X> c(n);
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if(possibly(x[j]<=l[j])) { c[j]=-1; }
        if(possibly(x[j]>=u[j])) { c[j]=+1; }
    }
    return c;
}

// Compute the variable types from the permutation, taking m basic variables and all non-basic variables lower.
template<class X>
array<Slackness>
compute_vt(const Vector<X>& l, const Vector<X>& u, const array<size_t>& p, const size_t m)
{
    const size_t n=p.size();
    array<Slackness> vt(n);
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
compute_p(const array<Slackness>& tv)
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


// Compute a basis (p_1,\ldots,p_m) for the matrix A
// Assume that the matrix A has full row rank
template<class X>
pair< array<size_t>, Matrix<X> >
SimplexSolver<X>::compute_basis(const Matrix<X>& A)
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


template<class X> template<class XX>
Matrix<XX>
SimplexSolver<X>::compute_B(const Matrix<X>& A, const array<size_t>& p)
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
Vector<X>
compute_c(const Matrix<X>& A, const array<size_t>& p, const Vector<X>& x)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    Vector<X> c(n);
    for(size_t k=m; k!=n; ++k) {
        if(x[p[k]]<0) { c[p[k]]=-1; }
    }
    return c;
}


template<class X> template<class XX>
Vector<XX>
SimplexSolver<X>::compute_x(const Matrix<X>& A, const Vector<X>& b, const array<size_t>& p, const Matrix<XX>& B)
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


template<class X> template<class XX>
Vector<XX>
SimplexSolver<X>::compute_x(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, const array<Slackness>& vt, const array<size_t>& p, const Matrix<XX>& B)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_ASSERT_MSG(p.size()==n, "vt="<<vt<<", p="<<p);

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

    Vector<XX> Ax=prod(A,x);
    Vector<XX> Axmb=Ax-b;
    ARIADNE_ASSERT_MSG(norm(Axmb)<0.00001,"A="<<A<<", b="<<b<<", l="<<l<<", u="<<u<<", vt="<<vt<<", p="<<p<<", x="<<x<<", Ax="<<Ax);
    return x;
}


template<class X,class XX>
std::pair<Vector<XX>,Vector<XX> >
compute_wx(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<Slackness>& vt, const array<size_t>& p, const Matrix<XX>& B)
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

template<class X,class XX>
Vector<XX>
compute_z(const Matrix<X>& A, const Vector<X>& c, const array<size_t>& p, const Vector<XX>& y)
{
    const double CUTOFF_THRESHOLD=1e-10;
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
        if(abs(z[j])<CUTOFF_THRESHOLD) {
            //z[j]=0;
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
compute_s(const size_t m, const array<Slackness>& vt, const array<size_t>& p, const Vector<X>& z)
{
    return compute_s_nocycling(m,vt,p,z);
}

// Compute the variable to enter the basis by finding the first which can increase the value function
template<class X>
size_t
compute_s_fast(const size_t m, const array<Slackness>& vt, const array<size_t>& p, const Vector<X>& z)
{
    const size_t n=z.size();
    for(size_t k=m; k!=n; ++k) {
        if( (vt[p[k]]==LOWER && z[p[k]]<0)
            || (vt[p[k]]==UPPER && z[p[k]]>0) ) { return k; }
    }
    return n;
}

// Compute the variable to enter the basis by giving the one with the highest rate of increase.
template<class X>
size_t
compute_s_best(const size_t m, const array<Slackness>& vt, const array<size_t>& p, const Vector<X>& z)
{
    const size_t n=z.size();
    size_t kmax=n;
    X abszmax=0;
    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        X posz=(vt[j]==LOWER ? -z[j] : z[j]);
        if(posz>abszmax) {
            kmax=k;
            abszmax=posz;
        }
    }
    return kmax;
}

// Compute the variable to enter the basis by using Bland's rule to avoid cycling.
template<class X>
size_t
compute_s_nocycling(const size_t m, const array<Slackness>& vt, const array<size_t>& p, const Vector<X>& z)
{
    const size_t n=z.size();
    for(size_t j=0; j!=n; ++j) {
        if( (vt[j]==LOWER && z[j]<0)
                || (vt[j]==UPPER && z[j]>0) ) {
            for(size_t k=m; k!=n; ++k) {
                if(p[k]==j) { return k; }
            }
            ARIADNE_ASSERT(false); // Should not reach here
        }
    }
    return n;
}

// Compute the direction in which the basic variables move if non-basic variable is increased.
// given by d=-B*A_s
template<class X>
Vector<X>
compute_d(const Matrix<X>& A, const array<size_t>& p, const Matrix<X>& B, const size_t ks)
{
    const double CUTOFF_THRESHOLD=1e-10;
    const size_t m=A.row_size();
    size_t js=p[ks];
    Vector<X> d(m);
    for(size_t k=0; k!=m; ++k) {
        for(size_t i=0; i!=m; ++i) {
            d[k]-=B[k][i]*A[i][js];
        }
        if(abs(d[k])<CUTOFF_THRESHOLD) { d[k]=0; }
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
compute_rt(const Vector<X>& l, const Vector<X>& u, const array<Slackness>& vt, const array<size_t>& p, const Vector<X>& x, const Vector<X>& d, const size_t s)
{
    // Choose variable to take out of basis
    // If the problem is degenerate, choose the variable with smallest index
    const size_t m=d.size();
    size_t r=s;
    X ds=(vt[p[s]]==LOWER ? +1 : -1);
    X t=u[p[s]]-l[p[s]];
    X tk=0;
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if(d[k]*ds<0 && x[j]>=l[j] && l[j] != -inf<X>()) {
            tk=(l[j]-x[j])/(ds*d[k]);
            if(t==inf<X>() || tk<t || (tk==t && p[k]<p[r]) ) { t=tk; r=k; }
        } else if( d[k]*ds>0 && x[j]<=u[j] && u[j] != inf<X>() ) {
            tk=(u[j]-x[j])/(ds*d[k]);
            if(t==inf<X>() || tk<t || (tk==t && p[k]<p[r])) { t=tk; r=k; }
        }
    }
    t*=ds;

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
    // Update x when variable p[s] becomes basic and variable p[r] becomes non-basic
    // The variable p[s] moves by t; the variables p[i] i<m by t*d[i]
    // The variable p[r] becomes an upper or lower variable depending on t*d[r]
    const size_t m=d.size();
    const size_t n=x.size();
    ARIADNE_ASSERT(r<m);
    ARIADNE_ASSERT(s>=m);
    ARIADNE_ASSERT(s<n);
    for(size_t i=0; i!=m; ++i) {
        x[p[i]]+=t*d[i];
    }
    x[p[s]] += t;
    if(t*d[r]<0) { x[p[r]]=l[p[r]]; }
    else if(t*d[r]>0) { x[p[r]]=u[p[r]]; }
    return;
}

template<class X>
void
update_x(const Vector<X>& l, const Vector<X>& u, const array<size_t>& p, Vector<X>& x, const size_t s, const Vector<X>& d, const X& t)
{
    // Update basis when a variable changes between lower and upper
    // The constant t determines how much the variable p[s] moves
    ARIADNE_ASSERT_MSG(s>=d.size(),"x="<<x<<" d="<<d<<" s="<<s<<"\n");
    ARIADNE_ASSERT_MSG(s<x.size(),"x="<<x<<" d="<<d<<" s="<<s<<"\n");
    const size_t m=d.size();
    for(size_t i=0; i!=m; ++i) {
        x[p[i]]+=t*d[i];
    }

    if(t>0) { x[p[s]]=u[p[s]]; }
    else if(t<0) { x[p[s]]=l[p[s]]; }
}




template<class X>
bool
SimplexSolver<X>::lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, array<size_t>& p, Matrix<X>& B, Vector<X>& x)
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
    ARIADNE_LOG(9,"  s="<<s<<" p[s]="<<p[s]<<"\n");


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
size_t lpenter(const Matrix<X>& A, const Vector<X>& c, const array<Slackness>& vt, const array<size_t>& p, const Matrix<X>& B)
{
    const size_t m=A.row_size();

    Vector<X> y=compute_y(c,p,B);
    Vector<X> z=compute_z(A,c,p,y);

    size_t s=compute_s(m,vt,p,z);
    ARIADNE_LOG(5,"  vt="<<vt<<" y="<<y<<" z="<<z<<" s="<<s<<" p[s]="<<p[s]<<"\n");
    return s;
}


template<class X>
void
SimplexSolver<X>::lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<Slackness>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x, size_t s)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    ARIADNE_ASSERT(s<=n);

    // Compute direction d in which to move the current basic variables
    // as the variable entering the basis changes by +1
    Vector<X> d=compute_d(A,p,B,s);


    // Compute distance t along d in which to move,
    // and the variable p[r] to leave the basis
    // The bounds on t are given by l <= x + t * d <= u
    // Note that t is negative if an upper variable enters the basis
    size_t r; X t;
    make_lpair(r,t)=compute_rt(l,u,vt,p,x,d,s);
    ARIADNE_LOG(7,"  s="<<s<<" p[s]="<<p[s]<<" r="<<r<<" p[r]="<<p[r]<<" d="<<d<<" t="<<t<<"\n");
    ARIADNE_LOG(5,"  s="<<s<<" p[s]="<<p[s]<<" r="<<r<<" p[r]="<<p[r]<<" d="<<d<<" t="<<t<<"\n");

    static const double DEGENERACY_THRESHOLD=1e-10;


    if(r==s) {
        Slackness nvts=(vt[p[s]]==LOWER ? UPPER : LOWER);
        ARIADNE_LOG(5,"   Changing non-basic variable "<<p[s]<<" in position "<<s<<" from type "<<vt[p[s]]<<" to type "<<nvts<<"\n");
    } else {
        ARIADNE_LOG(5,"   Swapping basic variable "<<p[r]<<" in position "<<r<<" with non-basic variable "<<p[s]<<" in position "<<s<<"\n");
    }

    if(abs(t)<DEGENERACY_THRESHOLD) {
        ARIADNE_LOG(3,"   Possible degeneracy changing non-basic variable "<<p[s]<<"=p["<<s<<"] with basic variable "<<p[r]<<"=p["<<r<<"] with d="<<d<<" and t="<<t<<"\n");
    }

    if(abs(t)<DEGENERACY_THRESHOLD && r!=s && false) {
        // Exact degeneracy, no need to update variables
        if(vt[p[s]]==LOWER) { vt[p[s]]=UPPER; }
        else { vt[p[s]]=LOWER; }
        vt[p[r]]=BASIS;
        std::swap(p[s],p[r]);
        B=compute_B<X>(A,p);
        x=compute_x(A,b,l,u,vt,p,B);
    } else if(r==s) {
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
        vt[p[s]] = BASIS;
        if(d[r]*t>0) {
            vt[p[r]] = UPPER;
        } else if(d[r]*t<0) {
            vt[p[r]] = LOWER;
        } else {
            size_t pr=p[r];
            ARIADNE_ASSERT(x[pr]==l[pr] || x[pr]==u[pr]);
            if(x[pr]==l[pr]) { vt[pr]=LOWER; } else { vt[pr]=UPPER; }
        }

        std::swap(p[r],p[s]);
    }

    ARIADNE_LOG(7,"      vt="<<vt<<"\n      p="<<p<<"\n");
    ARIADNE_LOG(7,"      B="<<B<<"\n      x="<<x<<"\n");

    consistency_check(A,b,l,u,vt,p,B,x);
}

template<class X>
bool
SimplexSolver<X>::lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u, array<Slackness>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x)
{
    ARIADNE_LOG(9,"  lpstep(A,b,c,l,u,vt,p,V,x)\n    A="<<A<<" b="<<b<<" c="<<c<<"\n    p="<<p<<" B="<<B<<"\n    vt="<<vt<<" l="<<l<<" x="<<x<<" u="<<u<<"\n");

    const size_t n=A.column_size();
    size_t s=lpenter(A,c,vt,p,B);
    if(s==n) { return true; }
    lpstep(A,b,l,u,vt,p,B,x,s);
    return false;
}





template<class X>
tribool
SimplexSolver<X>::_primal_feasible(const Matrix<X>& A, const Vector<X>& b, array<size_t>& p, Matrix<X>& B)
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
tribool
SimplexSolver<X>::_constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<Slackness>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x)
{
    ARIADNE_LOG(5,"\nInitial A="<<A<<" b="<<b<<"; l="<<l<<" u="<<u<<"\n  vt="<<vt<<"\n");
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Vector<X> cc(n);
    Vector<X> ll(l);
    Vector<X> uu(u);

    bool infeasible=false;
    for(size_t j=0; j!=n; ++j) {
        // If x[j] is infeasible by way of being to low, relax constraint x[j]>=l[j] to x[j]>=-inf.
        if(x[j]<l[j]) { cc[j]=-1; ll[j]=-inf<X>(); infeasible=true; }
        else if(x[j]>u[j]) { cc[j]=+1; uu[j]=+inf<X>(); infeasible=true; }
        else { cc[j]=0; }
    }
    ARIADNE_LOG(9,"    vt="<<vt<<" x="<<x<<" cc="<<cc<<"\n");

    static const int MAX_STEPS=1024;
    int steps=0;
    while(infeasible) {

        bool done=lpstep(A,b,cc,ll,uu,vt,p,B,x);
        ARIADNE_LOG(9,"  Done changing basis\n");
        ARIADNE_LOG(9,"    p="<<p<<" B="<<B<<"\n");
        ARIADNE_LOG(9,"    vt="<<vt<<" x="<<x<<"\n");

        if(done) {
            ARIADNE_LOG(9,"  Cannot put infeasible variables into basis.");
            Vector<X> y=compute_y(cc,p,B);
            Vector<X> yA=prod(y,A);
            X yb=dot(y,b);
            ARIADNE_LOG(5,"\nCertificate of infeasibility:\n y="<<y<<"\n yA="<<yA<<" yb="<<yb<<"\n");
            return false;
        }

        infeasible=false;
        for(size_t j=0; j!=n; ++j) {
            if(vt[j]==LOWER) { ARIADNE_ASSERT(x[j]==l[j]); }
            if(vt[j]==UPPER) { ARIADNE_ASSERT(x[j]==u[j]); }
            if(x[j]<l[j]) { cc[j]=-1; ll[j]=-inf<X>(); infeasible=true; }
            else if(x[j]>u[j]) { cc[j]=+1; uu[j]=+inf<X>(); infeasible=true; }
            else { cc[j]=0; ll[j]=l[j]; uu[j]=u[j]; }
        }
        ARIADNE_LOG(9,"\n    vt="<<vt<<" x="<<x<<" cc="<<cc<<"\n");

        ++steps;
        if(steps>=MAX_STEPS) {
            if(verbosity>0) {
                ARIADNE_WARN("WARNING: Maximum number of steps reached in constrained feasibility problem. "
                             <<"A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" cc="<<cc<<" ll="<<ll<<" uu="<<uu<<" vt="<<vt
                             <<" x="<<x<<" y="<<compute_y(cc,p,B)<<" Ay="<<Vector<X>(prod(compute_y(cc,p,B),A))<<"\n");
            }
            throw DegenerateFeasibilityProblemException();
        }
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

    ARIADNE_LOG(5,"\nFeasible point x="<<x<<"; l="<<l<<" u="<<u<<"\n Ax="<<Vector<X>(prod(A,x))<<" b="<<b<<"\n");

    return true;
}



// Check for feasibility of yA<=c. The solution is infeasible if there exists x with cx<0, Ax=0 and x>=0
template<class X>
tribool
SimplexSolver<X>::_dual_feasible(const Matrix<X>& A, const Vector<X>& c, array<size_t>& p, Matrix<X>& B)
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
tribool
SimplexSolver<X>::primal_feasible(const Matrix<X>& A, const Vector<X>& b)
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
tribool
SimplexSolver<X>::dual_feasible(const Matrix<X>& A, const Vector<X>& c)
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
tribool
SimplexSolver<X>::constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u)
{
    ARIADNE_LOG(2,"constrained_feasible(A,b,l,u)\n");
    ARIADNE_LOG(3,"    A="<<A<<"\n    b="<<b<<"\n    l="<<l<<"\n    u="<<u<<"\n");
    ARIADNE_ASSERT(b.size()==A.row_size());
    ARIADNE_ASSERT(l.size()==A.column_size());
    ARIADNE_ASSERT(u.size()==A.column_size());

    const size_t m=A.row_size();

    array<size_t> p;
    Matrix<X> B;
    make_lpair(p,B)=compute_basis(A);

    array<Slackness> vt=compute_vt(l,u,p,m);

    ARIADNE_LOG(9,"    p="<<p<<" B="<<B<<"  (BA="<<Matrix<X>(prod(B,A))<<")\n");

    Vector<X> x=compute_x(A,b,l,u,vt,p,B);

    tribool fs = _constrained_feasible(A,b,l,u,vt,p,B,x);
    tribool vfs = verify_constrained_feasibility(A,b,l,u,vt);
    ARIADNE_ASSERT(indeterminate(vfs) || vfs==fs);

    return vfs;
}



template<class X>
tribool
SimplexSolver<X>::constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<Slackness>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x, Vector<X>& y)
{
    const size_t m=A.row_size();
    //const size_t n=A.column_size();
    if(vt.size()==0) {
        make_lpair(p,B)=compute_basis(A);
        vt=compute_vt(l,u,p,m);
    }
    if(p.size()==0) {
        p=compute_p(vt);
        B=compute_B<X>(A,p);
    }

    x=compute_x(A,b,l,u,vt,p,B);
    ARIADNE_LOG(5,"A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<" x="<<x<<" Ax="<<Vector<X>(prod(A,x))<<"\n");

    tribool fs = _constrained_feasible(A,b,l,u,vt,p,B,x);

    Vector<X> c=compute_c(l,u,p,x,m);
    ARIADNE_LOG(9," c="<<midpoint(c)<<"\n");

    y=compute_y(c,p,B);
    ARIADNE_LOG(9," y="<<midpoint(y)<<"\n");

    return fs;
}

template<class X> struct RigorousNumericsTraits { typedef X Type; };
template<> struct RigorousNumericsTraits<Float> { typedef Interval Type; };

template<class X> tribool
SimplexSolver<X>::verify_primal_feasibility(const Matrix<X>& A, const Vector<X>& b, const array<Slackness>& vt)
{
    typedef typename RigorousNumericsTraits<X>::Type XX;
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
SimplexSolver<X>::verify_dual_feasibility(const Matrix<X>& A, const Vector<X>& c, const array<Slackness>& vt)
{
    ARIADNE_LOG(9,"verify_dual_feasibility(Matrix A, Vector c, VariableTypeArray vt):\n");
    typedef typename RigorousNumericsTraits<X>::Type XX;
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
SimplexSolver<X>::verify_constrained_feasibility(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, const array<Slackness>& vt)
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

    ARIADNE_LOG(9," x="<<x<<"\n   A*x="<<Vector<XX>(prod(A,x))<<"\n");

    tribool fs=true;
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        //if(possibly(x[j]<=l[j]) || possibly(x[j]>=u[j])) {
        if(possibly(x[j]<l[j]) || possibly(x[j]>u[j])) {
            ARIADNE_LOG(9," k="<<k<<" j="<<j<<" l[j]="<<l[j]<<" x[j]="<<x[j]<<" u[j]="<<u[j]<<"\n");
            fs=indeterminate;
            break;
        }
    }
    if(fs==true) {
        ARIADNE_LOG(9," true\n");
        return true;
    }

    Vector<X> c(n);
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        //if(possibly(x[j]<=l[j])) { c[j]=-1; }
        //if(possibly(x[j]>=u[j])) { c[j]=+1; }
        if(possibly(x[j]<l[j])) { c[j]=-1; }
        if(possibly(x[j]>u[j])) { c[j]=+1; }
    }
    ARIADNE_LOG(9," c="<<midpoint(c)<<"\n");

    const Vector<XX> y=compute_y(c,p,B);
    ARIADNE_LOG(9," y="<<y<<"\n");

    const Vector<XX> z=compute_z(A,c,p,y);
    ARIADNE_LOG(9," z="<<z<<"\n");

    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        if(vt[j]==LOWER && possibly(z[j]<=0)) { ARIADNE_LOG(9," indeterminate\n"); return indeterminate; }
        if(vt[j]==UPPER && possibly(z[j]>=0)) { ARIADNE_LOG(9," indeterminate\n"); return indeterminate; }
    }

    ARIADNE_LOG(9," false\n");
    return false;
}





template<class X>
Vector<X>
SimplexSolver<X>::optimize(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_ASSERT(b.size()==m);
    ARIADNE_ASSERT(c.size()==n);
    ARIADNE_ASSERT(l.size()==n);
    ARIADNE_ASSERT(u.size()==n);

    array<Slackness> vt(n);
    array<size_t> p(n);
    Matrix<X> B(m,m);
    Vector<X> x(n);

    make_lpair(p,B)=compute_basis(A);
    for(size_t k=0; k!=m; ++k) { vt[p[k]]=BASIS; }
    for(size_t k=m; k!=n; ++k) { vt[p[k]]=LOWER; }

    return optimize(A,b,c,l,u,vt,p,B);

}

template<class X>
Vector<X>
SimplexSolver<X>::optimize(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u, array<Slackness>& vt)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_ASSERT(b.size()==m);
    ARIADNE_ASSERT(c.size()==n);
    ARIADNE_ASSERT(l.size()==n);
    ARIADNE_ASSERT(u.size()==n);
    ARIADNE_ASSERT(vt.size()==n);
    ARIADNE_ASSERT(static_cast<size_t>(std::count(vt.begin(),vt.end(),BASIS))==m);

    array<size_t> p=compute_p(vt);
    Matrix<X> B=compute_B<X>(A,p);

    return optimize(A,b,c,l,u,vt,p,B);

}

template<class X>
Vector<X>
SimplexSolver<X>::optimize(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& c, const Vector<X>& l, const Vector<X>& u, array<Slackness>& vt, array<size_t>& p, Matrix<X>& B)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    if(p.size()==m) { p=extend_p(p,n); }

    ARIADNE_ASSERT(b.size()==m);
    ARIADNE_ASSERT(c.size()==n);
    ARIADNE_ASSERT(l.size()==n);
    ARIADNE_ASSERT(u.size()==n);
    ARIADNE_ASSERT(vt.size()==n);
    ARIADNE_ASSERT(p.size()==n);
    ARIADNE_ASSERT(B.row_size()==m);
    ARIADNE_ASSERT(B.column_size()==m);

    consistency_check(vt,p);
    consistency_check(A,p,B);

    Vector<X> x(n);
    x=compute_x(A,b,l,u,vt,p,B);
    ARIADNE_LOG(3,"Initial A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<" x="<<x<<" Ax="<<A*x<<"\n");

    _constrained_feasible(A,b,l,u,vt,p,B,x);
    ARIADNE_LOG(3,"Feasible A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<" x="<<x<<" Ax="<<A*x<<"\n");

    bool done=false;
    const int MAX_STEPS=1024;
    int steps=0;
    while(not done) {
        done=lpstep(A,b,c,l,u,vt,p,B,x);
        ++steps;
        ARIADNE_ASSERT_MSG(steps<MAX_STEPS,"Maximum number of steps reached for linear programming problem.");
    }
    ARIADNE_LOG(3,"Optimal A="<<A<<" b="<<b<<" c="<<c<<" l="<<l<<" u="<<u<<" vt="<<vt<<" p="<<p<<" x="<<x<<" Ax="<<A*x<<" cx="<<dot(x,x)<<"\n");

    return x;
}





template<class X>
tribool
SimplexSolver<X>::constrained_feasible_by_enumeration(const Matrix<X>& A, const Vector<X>& b,
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
        array<Slackness> vt(n);
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



template class SimplexSolver<Float>;

#ifdef HAVE_GMPXX_H
inline Interval operator+(const Interval& ivl, const Rational& q) { return ivl+Interval(q); }
inline Interval operator+(const Rational& q, const Interval& ivl) { return Interval(q)+ivl; }
inline Interval operator-(const Interval& ivl, const Rational& q) { return ivl-Interval(q); }
inline Interval operator-(const Rational& q, const Interval& ivl) { return Interval(q)-ivl; }
inline Interval operator*(const Interval& ivl, const Rational& q) { return ivl*Interval(q); }
inline Interval operator*(const Rational& q, const Interval& ivl) { return Interval(q)-ivl; }
template class SimplexSolver<Rational>;
#endif // HAVE_GMPXX_H



//--------------- Interior point solver -----------------------------------//


// Compute the product A:=XA, where X is diagonal
void ldmul(Matrix<Float>& A, const Vector<Float>& x) {
    ARIADNE_ASSERT(A.row_size()==x.size());
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint j=0; j!=A.row_size(); ++j) {
            A[i][j]*=x[i];
        }
    }
}

// Compute the product A:=AX, where X is diagonal
void rdmul(Matrix<Float>& A, const Vector<Float>& x) {
    ARIADNE_ASSERT(A.column_size()==x.size());
    for(uint j=0; j!=A.row_size(); ++j) {
        for(uint i=0; i!=A.row_size(); ++i) {
            A[i][j]*=x[j];
        }
    }
}

// Set r[i]=x[i]-c for i=1,...,n
Vector<Float> esub(const Vector<Float>& x, const Float& c) {
    Vector<Float> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]-c;
    }
    return r;
}

// Set r[i]=x[i]*y[i] for i=1,...,n
void _emul(Vector<Float>& r, const Vector<Float>& x, const Vector<Float>& y) {
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]*y[i];
    }
}

// Set r[i]=x[i]*y[i] for i=1,...,n
Vector<Float> emul(const Vector<Float>& x, const Vector<Float>& y) {
    Vector<Float> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]*y[i];
    }
    return r;
}

// Return r[i]=x[i]/z[i] for i=1,...,n
void _ediv(Vector<Float>& r, const Vector<Float>& x, const Vector<Float>& z) {
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]/z[i];
    }
}

// Return r[i]=x[i]*y[i] for i=1,...,n
Vector<Float> ediv(const Vector<Float>& x, const Vector<Float>& z) {
    Vector<Float> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]/z[i];
    }
    return r;
}

// Return r[i]=1/y[i] for i=1,...,n
Vector<Float> erev(const Vector<Float>& v) {
    Vector<Float> r(v.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=1.0/v[i];
    }
    return r;
}

// Compute R=ADA^T for diagonal D
Matrix<Float> mdtmul(const Matrix<Float>& A, const Vector<Float>& D)
{
    Matrix<Float> R(A.row_size(),A.row_size());
    ARIADNE_ASSERT(D.size()==A.column_size());
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint k=0; k!=A.column_size(); ++k) {
            Float ADik=A[i][k]*D[k];
            for(uint j=0; j!=A.row_size(); ++j) {
                R[i][j]+=ADik*A[j][k];
            }
        }
    }
    return R;
}

bool all_greater(const Vector<Float>& x, const Float& e) {
    for(uint i=0; i!=x.size(); ++i) {
        if(x[i]<=e) { return false; }
    }
    return true;
}



tuple<tribool,Vector<Float> >
InteriorPointSolver::constrained_dual_feasible(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
                                               const Vector<Float>& l, const Vector<Float>& u) const
{
    // Solve the constraint problem b<=Ax<=c; l<=x<=u.
    ARIADNE_NOT_IMPLEMENTED;
}

Interval dot(const Vector<Float>& c, const Vector<Interval>& x) {
    Interval r=0.0;
    for(uint i=0; i!=x.size(); ++i) {
        r+=c[i]*x[i];
    }
    return r;
}

Interval
InteriorPointSolver::validate(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c, const Vector<Float>& x, const Vector<Float>& y) const
{
    // Find a primal-dual feasible point (x,y) suitable as an input
    // to an optimization solver. The slack variables are A
    ARIADNE_ASSERT(A.row_size()==b.size());
    ARIADNE_ASSERT(A.column_size()==c.size());

    Matrix<Interval> IA=A;
    Matrix<Interval> IS=prod(transpose(IA),IA);
    Matrix<Interval> ISinv=inverse(IS);
    Vector<Interval> dx=IA*Vector<Interval>(x); dx=IS*dx; dx=dx*IA;
    Vector<Interval> Ix=x-dx;
    Vector<Interval> Iy(y);
    Vector<Interval> Iz=prod(y,IA)-c;
    for(uint i=0; i!=Ix.size(); ++i) {
        ARIADNE_ASSERT(Ix[i].lower()>=0.0);
    }
    for(uint i=0; i!=Iz.size(); ++i) {
        ARIADNE_ASSERT(Iz[i].lower()>=0.0);
    }
    return hull(dot(c,Ix),dot(b,Iy));
}

tribool
slow_feasible(const Matrix<Float>& A, const Vector<Float>& cl, const Vector<Float>& cu,
              const Vector<Float>& l, const Vector<Float>& u)
{
    const uint m=A.row_size();
    const uint n=A.column_size();

    ARIADNE_ASSERT(cl.size()==n && cu.size()==n && l.size()==m && u.size()==m);

    const uint mm=m+1;
    const uint nn=2*(m+n);

    Matrix<Float> AA(mm,nn);
    for(uint i=0; i!=m; ++i) {
        for(uint j=0; j!=n; ++j) {
            AA[i][j]=-A[i][j];
            AA[i][j+n]=A[i][j];
        }
        AA[i][2*n+i]=-1;
        AA[i][2*n+m+i]=+1;
    }
    for(uint jj=0; jj!=nn; ++jj) {
        AA[m][jj]=1.0;
    }

    ARIADNE_LOG(2,"AA="<<AA<<"\n");

    Vector<Float> bb(mm);
    bb[m]=1.0;

    ARIADNE_LOG(2,"bb="<<bb<<"\n");

    Vector<Float> cc(nn);
    for(uint j=0; j!=n; ++j) {
        cc[j]=-cl[j];
        cc[n+j]=cu[j];
    }
    for(uint i=0; i!=m; ++i) {
        cc[2*n+i]=-l[i];
        cc[2*n+m+i]=u[i];
    }

    ARIADNE_LOG(2,"cc="<<cc<<"\n");

    Vector<Float> xx(nn);
    Vector<Float> yy(mm);
    Vector<Float> zz(nn);

    for(uint jj=0; jj!=nn; ++jj) {
        xx[jj]=1.0/nn;
    }

    ARIADNE_LOG(2,"xx="<<xx<<"\n");

    Vector<Float> y(m);
    for(uint i=0; i!=m; ++i) {
        y[i]=(l[i]+u[i])/2;
        yy[i]=y[i];
        zz[2*n+i]=y[i]-l[i];
        zz[2*n+m+i]=u[i]-y[i];
    }

    Vector<Float> yA=prod(y,A);
    ARIADNE_LOG(2,"yA="<<yA<<"\n");

    Float t=+inf<Float>();
    for(uint j=0; j!=n; ++j) {
        zz[j]=yA[j]-cl[j];
        t=min(t,zz[j]);
        zz[n+j]=cu[j]-yA[j];
        t=min(t,zz[n+j]);
    }

    t-=1.0/nn;
    for(uint i=0; i!=nn; ++i) {
        zz[i]-=t;
    }
    yy[m]=t;

    ARIADNE_LOG(2,"yy="<<yy<<"\n");
    ARIADNE_LOG(2,"zz="<<zz<<"\n");

    InteriorPointSolver()._optimize(AA,bb,cc,xx,yy,zz);
    ARIADNE_LOG(2," xx="<<xx<<" yy="<<yy<<" zz="<<zz<<"\n");

    y=project(yy,range(0,m));
    t=yy[m];
    ARIADNE_LOG(2,"b="<<cl<<" yA="<<prod(y,A)<<" c="<<cu<<"\nl="<<l<<" y="<<y<<" u="<<u<<"\n");

    if(t>0.0) { return true; }
    else if (t<0.0) { return false; }
    else { return indeterminate; }

}


tribool
InteriorPointSolver::feasible(const Matrix<Float>& A, const Vector<Float>& cl, const Vector<Float>& cu,
                              const Vector<Float>& l, const Vector<Float>& u) const
{
    return slow_feasible(A,cl,cu,l,u);

    // Find a primal-dual feasible point (x,y) suitable as an input
    // to an optimization solver. The slack variables are A
    static const uint m=A.row_size();
    static const uint n=A.column_size();

    ARIADNE_ASSERT(cl.size()==n);
    ARIADNE_ASSERT(cu.size()==n);
    ARIADNE_ASSERT(l.size()==m);
    ARIADNE_ASSERT(u.size()==m);

    // Find vectors y,t,zcl,zcu,zl,zu such that A^Ty+et+zc=c etc
    Vector<Float> y(m);
    Vector<Float> z(2*(m+n));
    for(uint i=0; i!=m; ++i) {
        y[i]=(l[i]+u[i])/2;
        z[2*n+i]=y[i]-l[i];
        z[2*n+m+i]=u[i]-y[i];
    }

    Vector<Float> yA=prod(y,A);
    Float t=+inf<Float>();
    for(uint j=0; j!=n; ++j) {
        z[j]=yA[j]-cl[j];
        z[n+j]=cu[j]-yA[j];
        t=min(t,cu[j]-yA[j]);
        t=min(t,yA[j]-cl[j]);
    }
    if(t>0) {
        ARIADNE_LOG(2,"y="<<y<<" yA="<<yA<<" cl="<<cl<<" cu="<<cu<<"\n");
        return true;
    }
    t-=1.0/(2*(m+n));
    for(uint i=0; i!=2*(m+n); ++i) {
        z[i]-=t;
    }

    // Find a vector x=(xcl,xcu,xl,xu) such that A(xcu-xcl)+(xu-xl) = 0 summing to 1
    Vector<Float> x(2*(m+n));
    for(uint i=0; i!=2*(m+n); ++i) {
        x[i]=1.0/(2*(m+n));
    }

    Vector<Float> et(n,t);
    Vector<Float> zcu=project(z,range(n,2*n));
    Vector<Float> zcl=project(z,range(0,n));

    ARIADNE_LOG(2,"A="<<A<<" cl="<<cl<<" cu="<<cu<<" l="<<l<<" u="<<u<<"\n");
    ARIADNE_LOG(2,"x="<<x<<" t="<<t<<" y="<<y<<" z="<<z<<"\n");

    ARIADNE_ASSERT(norm(Vector<Float>(yA-et-zcl-cl))<1e-10);
    ARIADNE_ASSERT(norm(Vector<Float>(yA+et+zcu-cu))<1e-10);

    return _feasible(A,cl,cu,l,u,t,x,y,z);
}

tuple< Vector<Float>, Vector<Float>, Vector<Float> >
InteriorPointSolver::optimize(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c) const
{
    Vector<Float> x,y,z;
    //make_ltuple(x,y,z) = this->feasible(A,b,c);
    return this->_optimize(A,b,c,x,y,z);
}



void _set_es_matrix(Matrix<Float>& S, const Matrix<Float>& A, const Vector<Float>& D) {
    // D=diag(Db,Dc,Dl,Du)
    // S=[ A*(Dc+Db)*AT + (Du+Dl),     A*(Dc-Db)*e + (Du-Dl)*e ]
    //   [ eT*(Dc-Db)*AT + eT*(Du-Dl), eT*(Dc+Db+Du+Dl)*e      ]
    static const uint m=A.row_size();
    static const uint n=A.column_size();

    ARIADNE_ASSERT(S.row_size()==m+1);

    const Float* Db=&D[0];
    const Float* Dc=Db+n;
    const Float* Dl=Dc+n;
    const Float* Du=Dl+m;
    Float* Sptr=&S[0][0];

    //Clear matrix
    for(uint i=0; i!=(m+1u)*(m+1u); ++i) {
        Sptr[i]=0.0;
    }

    // Set upper-triangular values
    for(uint i=0; i!=m; ++i) {
        for(uint k=0; k!=n; ++k) {
            Float ADik=A[i][k]*(Db[k]+Dc[k]);
            for(uint j=i; j!=m; ++j) {
                S[i][j]+=ADik*A[j][k];
            }
            S[i][m]+=A[i][k]*(Dc[k]-Db[k]);
        }
        S[i][i]+=(Du[i]+Dl[i]);
        S[i][m]+=(Du[i]-Dl[i]);
    }
    for(uint h=0; h!=(2*(m+n)); ++h) {
        S[m][m]+=D[h];
    }

    // Set lower part
    for(uint i=1; i!=m+1; ++i) {
        for(uint j=0; j!=i; ++j) {
            S[i][j]=S[j][i];
        }
    }
}


void _set_s_matrix(Matrix<Float>& S, const Matrix<Float>& A, const Vector<Float>& Dc, const Vector<Float>& Dd) {
    // S=A*Dc*AT+Dd

    static const uint m=A.row_size();
    static const uint n=A.column_size();

    ARIADNE_ASSERT(S.row_size()==m);
    ARIADNE_ASSERT(S.column_size()==m);
    ARIADNE_ASSERT(Dc.size()==n);
    ARIADNE_ASSERT(Dd.size()==m);

    //Clear matrix
    Float* Sptr=&S[0][0];
    for(uint i=0; i!=m*m; ++i) {
        Sptr[i]=0.0;
    }

    // Set upper-triangular values
    for(uint i=0; i!=m; ++i) {
        S[i][i]=Dd[i];
        for(uint k=0; k!=n; ++k) {
            Float ADik=A[i][k]*Dc[k];
            for(uint j=i; j!=m; ++j) {
                S[i][j]+=ADik*A[j][k];
            }
        }
    }
    // Set lower part
    for(uint i=1; i!=m; ++i) {
        for(uint j=0; j!=i; ++j) {
            S[i][j]=S[j][i];
        }
    }
}



/*
// Consider feasibility of cl<=yA<=cu; l<=y<=u
// Augment system by considering b<=yA-te; yA+te<=c; max(t)
tribool
InteriorPointSolver::_feasible(const Matrix<Float>& A, const Vector<Float>& cl, const Vector<Float>& cu,
                               const Vector<Float>& yl, const Vector<Float>& yu,
                               Float& v, Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const
{
    const uint m=A.row_size();
    const uint n=A.column_size();

    Vector<Float> e(n,1.0);
    Vector<Float> yA(n);
    Vector<Float> Dc(n);
    Vector<Float> Db(n);

    Matrix<Float> S(m,m), Sinv(m,m);
    Vector<Float> yf(m);

    FloatVector dy(m); dx(2*(m+n)); dz(2*(m+n));
    FloatVector ny(m); nx(2*(m+n)); nz(2*(m+n));

    FloatVectorRange zcl=project(z,range(0,n));
    FloatVectorRange zcu=project(z,range(n,2*n));
    FloatVectorRange zbl=project(z,range(2*n,2*n+m));
    FloatVectorRange zbu=project(z,range(2*n+m,2*(n+m)));

    FloatVectorRange xcl=project(x,range(0,n));
    FloatVectorRange xcu=project(x,range(n,2*n));
    FloatVectorRange xbl=project(x,range(2*n,2*n+m));
    FloatVectorRange xbu=project(x,range(2*n+m,2*(n+m)));

    FloatVectorRange dzcl=project(dz,range(0,n));
    FloatVectorRange dzcu=project(dz,range(n,2*n));
    FloatVectorRange dzbl=project(dz,range(2*n,2*n+m));
    FloatVectorRange dzbu=project(dz,range(2*n+m,2*(n+m)));

    FloatVectorRange dxcl=project(dx,range(0,n));
    FloatVectorRange dxcu=project(dx,range(n,2*n));
    FloatVectorRange dxbl=project(dx,range(2*n,2*n+m));
    FloatVectorRange dxbu=project(dx,range(2*n+m,2*(n+m)));

    for(uint i=0; i!=10; ++i) {
        Float eps=0.01;

        Vector<Float> res(2*(m+n));
        for(uint i=0; i!=2*(m+n); ++i) {
            res[i]=x[i]*z[i]-eps;
        }

        Dc=ediv(xcu,zcu)-ediv(xcl,zcl);
        Db=ediv(xbu,zbu)-ediv(xbl,zbl);
        Drzc=erec(zcu)-erec(zcl);
        Drzb=erec(zbu)-erec(zbl);

        yf=prod(A,Drz)+Drzb;
        _set_s_matrix(S,A,Dc,Db);
        Sinv=inverse(S);
        dy=prod(Sinv,yf);
        dyA=prod(y,A);
        dv=0.0;

        dzcl=dyA-dv*e;
        dzcu=-dv*e-dyA;
        dzyl=dy;
        dzyu=-dy;

        dx=ediv(-res-emul(x,dz),z);

        bool allpositive=false;
        Float a=1/scale;
        while(!allpositive) {
            a=a*scale;
            nx=midpoint(x+a*dx);
            nz=midpoint(z+a*dz);
            nxz=emul(nx,nz);
            allpositive = all_greater(nx,0.0) && all_greater(nxz,eps);
            ARIADNE_LOG(2,"    a="<<a<<" e="<<neighbourhood<<" nx="<<nx<<" nt="<<nt<<" ny="<<ny<<" nz="<<nz<<"\n");
        }
        nt=midpoint(t+a*dt);
        ny=midpoint(y+a*dy);

        t=nt; y=ny; x=nx; z=nz;
    }

    if(t>0.0) { return true; }
    else if(cx<0.0) { return false; }
    else { return indeterminate; }
}
*/

// Consider feasibility of  bl<=yB<=bu; yC<=c; dl<=y<=du
// Augment system by considering b<=yA-ve; yA+ve<=c; max(t)
tribool
InteriorPointSolver::feasible(const Matrix<Float>& B, const Vector<Float>& bl, const Vector<Float>& bu,
                              const Matrix<Float>& C, const Vector<Float>& c,
                              const Vector<Float>& dl, const Vector<Float>& du) const
{
    StandardLinearProgram lp=feasibility_problem(B,bl,bu,C,c,dl,du);

    Matrix<Float>& AA=lp.A;
    Vector<Float>& bb=lp.b;
    Vector<Float>& cc=lp.c;
    Vector<Float>& xx=lp.x;
    Vector<Float>& yy=lp.y;
    Vector<Float>& zz=lp.z;

    Float yb, cx;

    const uint max_error=1e-8;
    yb=dot(yy,bb); cx=dot(cc,xx);

    while(yb<=0.0 && cx>=0.0 && (cx-yb)>max_error) {
        this->_optimization_step(AA,bb,cc,xx,yy,zz);
        yb=dot(yy,bb); cx=dot(cc,xx);
    }

    ARIADNE_LOG(2,"cx="<<cx<<"\nyb="<<yb<<"\n");
    Matrix<Interval> iAA(AA);
    Vector<Interval> iyy(yy);
    Vector<Interval> ixx(xx);
    Vector<Interval> icc(cc);
    Vector<Interval> ibb(bb);
    ARIADNE_LOG(2,"c-yA="<<icc-iyy*iAA<<"\n");

    Matrix<Interval> S=sprod(iAA);
    Matrix<Interval> Sinv=inverse(S);
    Vector<Interval> izz=ibb-iAA*ixx;
    ixx=transpose(iAA)*(Sinv*(izz))+ixx;
    ARIADNE_LOG(2,"Ax-b="<<prod(iAA,ixx)-bb<<"\n");
    ARIADNE_LOG(2,"x="<<ixx<<"\n");

    if(yb>0.0) { return true; }
    else if(cx<0.0) { return false; }
    else { return indeterminate; }
}


/*
// Consider feasibility of  bl<=yB<=bu; yC<=c; dl<=y<=du
// Augment system by considering b<=yA-ve; yA+ve<=c; max(t)
tribool
InteriorPointSolver::_feasible(const Matrix<Float>& B, const Vector<Float>& bl, const Vector<Float>& bu,
                               const Matrix<Float>& C, const Vector<Float>& c,
                               const Vector<Float>& dl, const Vector<Float>& du,
                               Float& t, Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const
{
    // For the linear system Ax+by=s; b^Tx+cy=t, a solution is given by
    // y=(t - b^T A^-1 s)/(c-b^T A^-1 b); x=A^-1 s - A^-1 b y

    typedef boost::numeric::ublas::vector_range< Vector<Float> > FloatVectorRange;
    typedef Vector<Float> FloatVector;

    const uint nd=y.size();
    const uint nc=C.column_size();
    const uint nb=B.column_size();
    const uint ns=2*nb+nc+2*nd;

    ARIADNE_ASSERT(B.row_size()==nd);
    ARIADNE_ASSERT(C.row_size()==nd);
    ARIADNE_ASSERT(bl.size()==nb);
    ARIADNE_ASSERT(bu.size()==nb);
    ARIADNE_ASSERT(c.size()==nc);
    ARIADNE_ASSERT(dl.size()==nd);
    ARIADNE_ASSERT(du.size()==nd);
    ARIADNE_ASSERT(x.size()==ns);
    ARIADNE_ASSERT(z.size()==ns);

    Float dt,nt;
    Vector<Float> dx(2*(m+n)),dy(m),dz(2*(m+n));
    Vector<Float> nx(2*(m+n)),ny(m),nz(2*(m+n)),nxz(2*(m+n));
    const Vector<Float> em(m,1u);
    const Vector<Float> en(n,1u);

    Matrix<Float> S(m+1,m+1),Sinv(m+1,m+1);
    Vector<Float> xdivz(2*(n+m));
    Vector<Float> xmulz(2*(n+m));
    Vector<Float> tau(2*(n+m));

    Vector<Float> c(2*(m+n));
    project(c,range(0,n))=cl;
    project(c,range(n,n+n))=cu;
    project(c,range(n+n,n+n+m))=l;
    project(c,range(n+n+m,n+n+m+m))=u;

    const Float scale=0.75;
    const Float maxerror=1e-3;
    const Float starterror=1.0/16;
    const Float errorscale=1.0/2;
    const Float startcentering=1.0/8;
    const Float gamma=1.0/1024;

    const uint maxsteps=10;
    uint steps=0;

    Float cx,yb;
    cx=dot(c,x);
    yb=t;
    ARIADNE_ASSERT(yb<=cx);

    ARIADNE_LOG(2,"A="<<A<<" cl="<<cl<<" cu="<<cu<<"\n");
    ARIADNE_LOG(2," t="<<t<<" x="<<x<<" y="<<y<<" z="<<z<<"\n");

    const Float centering=startcentering; // sigma
    Float duality; // mu = dot(x,z)/n
    Float neighbourhood; // gamma = sigma mu
    while(steps<maxsteps && cx>=0.0 && yb<=0.0) {
        duality=dot(x,z)/n; // mu
        neighbourhood=gamma*duality;
        _ediv(xdivz,x,z);
        _ediv(xmulz,x,z);
        _emul(tau,x,z);
        for(uint i=0; i!=2*(m+n); ++i) { tau[i]-=centering*duality; }
        for(uint i=0; i!=2*(m+n); ++i) { tau[i]/=z[i]; }
        ARIADNE_LOG(2,"  x*z="<<xmulz<<"  x/z="<<xdivz<<" tau="<<tau<<"\n");
        _set_es_matrix(S,A,xdivz);
        Sinv=inverse(S);
        ARIADNE_LOG(2,"  S="<<S<<"  inverse(S)="<<Sinv<<"\n");

        //dy=Sinv*(rb-A*(iZ*(X*rc))-A*(iZ*t))
        //dx=iZ*(X*(At*dy-rc)+t)
        //dz=iX*(t-Z*dx)

        // Compute the residuals on y
        // The general formula is A*Zinv*tau, where tau is the residual for the slack variables

        Vector<Float> rzx=tau;
        Vector<Float> ryt(m+1,0.0);
        project(ryt,range(0,m))=prod(A,(project(rzx,range(n,n+n))-project(rzx,range(0,n))))
            +(project(rzx,range(n+n+m,n+n+m+m))-project(rzx,range(n+n,n+n+m)));
        for(uint i=0; i!=2*(m+n); ++i) { ryt[m]+=rzx[i]; }
        Vector<Float> dyt=Sinv*ryt;
        dy=project(dyt,range(0,m));
        dt=dyt[m];

        FloatVector dyA(n);
        FloatVectorRange dzb=project(dz,range(0,n));
        FloatVectorRange dzc=project(dz,range(n,2*n));
        FloatVectorRange dzl=project(dz,range(2*n,2*n+m));
        FloatVectorRange dzu=project(dz,range(2*n+m,2*n+2*m));

        dyA=dy*A;

        dzb=dyA-dt*en;
        dzc=-dyA-dt*en;
        dzl=dy-dt*em;
        dzu=-dy-dt*em;

        dx=ediv(-tau-emul(x,dz),z);

        ARIADNE_LOG(2,"  dx="<<dx<<" dy="<<dy<<" dz="<<dz<<"\n");
        bool allpositive=false;
        Float a=1/scale;
        while(!allpositive) {
            a=a*scale;
            nt=midpoint(t+a*dt);
            nx=midpoint(x+a*dx);
            ny=midpoint(y+a*dy);
            nz=midpoint(z+a*dz);
            nxz=emul(nx,nz);
            allpositive = all_greater(nx,0.0) && all_greater(nxz,neighbourhood);
            ARIADNE_LOG(2,"    a="<<a<<" e="<<neighbourhood<<" nx="<<nx<<" nt="<<nt<<" ny="<<ny<<" nz="<<nz<<"\n");
           // ARIADNE_LOG(2,"    a="<<a<<" e="<<neighbourhood<<" nx="<<nx<<" ny="<<ny<<" nz="<<nz<<" nxz="<<nxz<<"\n");
        }
        t=nt;
        x=nx;
        y=ny;
        z=nz;
        ARIADNE_LOG(2,"  cx="<<dot(c,x)<<" yb="<<t<<"\n");
        ++steps;
    }


    if(yb>0.0) { return true; }
    else if(cx<0.0) { return false; }
    else { return indeterminate; }

}
*/


tuple< Vector<Float>, Vector<Float>, Vector<Float> >
InteriorPointSolver::_optimize(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
                               Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const
{
    ARIADNE_ASSERT(A.row_size()==b.size());
    ARIADNE_ASSERT(A.column_size()==c.size());
    ARIADNE_ASSERT(A.column_size()==x.size());
    ARIADNE_ASSERT(A.row_size()==y.size());
    ARIADNE_ASSERT(A.column_size()==z.size());

    const double maxerror=1e-3;
    const uint maxsteps=10;

    Float cx,yb;
    cx=dot(c,x);
    yb=dot(y,b);
    ARIADNE_ASSERT(yb<=cx);

    ARIADNE_LOG(2,"A="<<A<<" b="<<b<<" c="<<c<<"\n");
    ARIADNE_LOG(2,"x="<<x<<" y="<<y<<" z="<<z<<"\n");

    uint steps=0;
    while(steps<maxsteps && (cx-yb)>maxerror) {
        this->_optimization_step(A,b,c,x,y,z);
        ++steps;
    }

    return make_tuple(x,y,z);

}

void
InteriorPointSolver::_optimization_step(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
                                        Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const
{
    static const double gamma=1.0/1024;
    static const double sigma=1.0/8;
    static const double scale=0.75;

    const uint m=A.row_size();
    const uint n=A.column_size();

    Vector<Float> dx(n),dy(n),dz(n);
    Vector<Float> nx(n),ny(n),nz(n),nxz(n);
    Vector<Float> rb(m),rc(n),rs(n),ry;
    Matrix<Float> S(m,m),Sinv(m,m);

    Float mu=dot(x,z)/n; // duality constant

    rb=A*x-b;
    rc=y*A+z-c;
    rs=esub(emul(x,z),sigma*mu);

    ry=A*(ediv(rs-emul(x,rc),z))-rb;

    S=mdtmul(A,ediv(x,z));
    Sinv=inverse(S);

    ARIADNE_LOG(4,"  S="<<S<<"  inverse(S)="<<Sinv<<"\n");

    dy=Sinv*ry;
    dz=-rc-dy*A;
    dx=-ediv((rs+emul(x,dz)),z);

    ARIADNE_LOG(4,"  dx="<<dx<<" dy="<<dy<<" dz="<<dz<<"\n");
    bool allpositive=false;
    Float a=1/scale;
    while(!allpositive) {
        a=a*scale;
        nx=midpoint(x+a*dx);
        nz=midpoint(z+a*dz);
        allpositive = all_greater(nx,0.0) && all_greater(nz,0.0) && all_greater(emul(nx,nz),gamma*mu);
        ARIADNE_LOG(6,"    a="<<a<<" e="<<eps<<" nx="<<nx<<" nz="<<nz<<"\n");
    }
    ny=midpoint(y+a*dy);

    x=nx; y=ny; z=nz;

    ARIADNE_LOG(2,"  cx="<<dot(c,x)<<" yb="<<dot(y,b)<<"\n");
    ARIADNE_LOG(2,"  x="<<x<<" y="<<y<<" z="<<z<<"\n");
    ARIADNE_LOG(2,"  Ax-b="<<prod(A,x)-b<<" c-yA="<<c-prod(y,A)<<"\n");
}


} // namespace Ariadne

