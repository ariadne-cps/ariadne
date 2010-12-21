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

class DegenerateFeasibilityProblemException : public std::runtime_error {
  public:
    DegenerateFeasibilityProblemException() : std::runtime_error("") { }
    DegenerateFeasibilityProblemException(const std::string& what) : std::runtime_error(what) { }
};

template<class X> inline Vector<X> operator*(const Matrix<X>& A, const Vector<X>& x) {
    return prod(A,x); }
template<class X> inline Vector<X> operator*(const Vector<X>& y, const Matrix<X>& A) {
    return prod(y,A); }
template<class X> inline Matrix<X> operator*(const Matrix<X>& A, const Matrix<X>& B) {
    return prod(A,B); }

StandardLinearProgram::StandardLinearProgram(uint m, uint n)
    : A(m,n), b(m), c(n), x(n), y(m), z(n) { }

// Threshold for matrix diagonal elements below which it may be considered singular
static const double SINGULARITY_THRESHOLD = std::numeric_limits<double>::epsilon();

// Threshold for slack variables z_N = (c_N - c_B A_B^-1 A_N) below which basis
// may be considered optimal
static const double PROGRESS_THRESHOLD = std::numeric_limits<double>::epsilon() * 1024;

// Threshold for primal variable direction of change below which we do not
// need to see if it would violate its constraints
const double CUTOFF_THRESHOLD=std::numeric_limits<double>::epsilon() * 16;

// Threshold for primal variables to exceed their bounds
const double BOUNDS_TOLERANCE=std::numeric_limits<double>::epsilon() * 4098;

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
    return os << (t==BASIS ? 'B' : t==LOWER ? 'L' : t==UPPER ? 'U' : t==FIXED ? 'E' : '?');
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


template<class X>
size_t
SimplexSolver<X>::consistency_check(const array<Slackness>& vt, const array<size_t>& p) {
    return Ariadne::consistency_check(vt,p);
}

// Check that the matrix B is the inverse of the matrix A_B with columns p[0],...,p[m-1] of A.
template<class X>
void
SimplexSolver<X>::consistency_check(const Matrix<X>& A, const array<size_t>& p, const Matrix<X>& B)
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
    ARIADNE_ASSERT_MSG(norm(Z)<MAXIMUM_ERROR, "A="<<A<<"\np="<<p<<"\nB="<<B<<"\nZ=B*A_B-I="<<Z<<"\nnorm(Z)="<<norm(Z));
}


// Check that the matrix B is the inverse of the matrix A_B with columns p[0],...,p[m-1] of A and that Ax=b.
template<class X>
void
SimplexSolver<X>::consistency_check(const Matrix<X>& A, const Vector<X>& b,
                                    const array<size_t>& p, const Matrix<X>& B, const Vector<X>& x)
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
    ARIADNE_ASSERT_MSG(norm(Z)<MAXIMUM_ERROR, "A="<<A<<"p="<<p<<"\nB*A_B-I="<<Z<<"\nnorm(Z)="<<norm(Z));

    Vector<X> z=prod(A,b);
    for(size_t j=0; j!=n; ++j) { z[j]-=x[j]; }
    ARIADNE_ASSERT(norm(z)<MAXIMUM_ERROR);
}


// Check that the matrix B is the inverse of the matrix A_B with columns p[0],...,p[m-1] of A, and that
// the vector x is given by x_L=l_L, x_U=x_U and x_B=B^{-1} A_N x_N.
template<class X>
void
SimplexSolver<X>::consistency_check(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u,
                                    const array<Slackness>& vt, const array<size_t>& p, const Matrix<X>& B, const Vector<X>& x)
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
        ARIADNE_ASSERT_MSG(vt[j]==LOWER || vt[j]==UPPER,
                           "vt["<<j<<"]="<<vt[j]<<"\n  A="<<A<<", b="<<b<<", l="<<l<<", u="<<u<<", vt="<<vt<<", p="<<p<<", x="<<x<<", Ax="<<Ax);
        X xj = (vt[j]==LOWER ? l[j] : u[j]);
        ARIADNE_ASSERT_MSG(x[j]==xj,"x["<<j<<"]="<<x[j]<<" xj="<<xj<<"\n  A="<<A<<", b="<<b<<", l="<<l<<", u="<<u<<", vt="<<vt<<", p="<<p<<", x="<<x<<", Ax="<<Ax);
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
// Throws an error if the matrix A has full row rank
template<class X>
pair< array<size_t>, Matrix<X> >
SimplexSolver<X>::compute_basis(const Matrix<X>& A)
{
    ARIADNE_LOG(9,"compute_basis(A) with A="<<A<<"\n");
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_DEBUG_ASSERT(n>=m);

    array<size_t> p(n);
    for(size_t j=0; j!=n; ++j) { p[j]=j; }

    // Factorise into lower and upper triangular matrices L and U
    Matrix<X> L=Matrix<X>::identity(m);
    Matrix<X> U=A;

    for(size_t k=0; k!=m; ++k) {
        // Find a good pivot column j and swap entries of jth and kth columns

        // Look for a column which is the unit vector ek below k
        bool found_unit=false;
        size_t j;
        for(j=k; j!=n; ++j) {
            if(U[k][j]==+1) {
                found_unit=true;
                for(uint i=k+1; i!=m; ++i) {
                    if(U[i][j]!=0) { found_unit=false; break; }
                }
            }
            if(found_unit) { break; }
        }

        // Look for a column with largest U[k][j]
        if(!found_unit) {
            X Ukjmax = abs(U[k][k]);
            size_t jmax=k;
            for(j=k+1; j!=n; ++j) {
                X absUkj=abs(U[k][j]);
                if(absUkj>Ukjmax) {
                    Ukjmax=absUkj;
                    jmax=j;
                }
            }
            j=jmax;
        }

        if (abs(U[k][j]) < SINGULARITY_THRESHOLD) { ARIADNE_THROW(SingularLinearProgram,"compute_basis"," matrix A="<<A<<" is singular or very nearly singular"); }

        if(j!=k) {
            std::swap(p[k],p[j]);
            for(size_t i=0; i!=m; ++i) {
                std::swap(U[i][k],U[i][j]);
            }
        }

        ARIADNE_DEBUG_ASSERT(U[k][k]!=0);

        if(!found_unit) {
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

    } // end loop on diagonal k

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

template<class X,class XX,class XXX>
Vector<XX>
compute_z(const Matrix<X>& A, const Vector<XXX>& c, const array<size_t>& p, const Vector<XX>& y)
{
    const double CUTOFF_THRESHOLD=1e-10;
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    Vector<XX> z(n);

    // Shortcut if we can assume c_B - y A_B = 0
    for(size_t k=0; k!=m; ++k) {
        z[p[k]]=0;
    }
    for(size_t k=0; k!=n; ++k) { // Change to k=m,...,n if c_B - y A_B = 0
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
        if(z[p[k]]< -PROGRESS_THRESHOLD) { return k; }
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
        if( (vt[p[k]]==LOWER && z[p[k]]< -PROGRESS_THRESHOLD)
            || (vt[p[k]]==UPPER && z[p[k]]> +PROGRESS_THRESHOLD) ) { return k; }
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
        if( (vt[j]==LOWER && z[j]< -PROGRESS_THRESHOLD)
                || (vt[j]==UPPER && z[j]> +PROGRESS_THRESHOLD) ) {
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
template<class X, class XX>
Vector<XX>
compute_d(const Matrix<X>& A, const array<size_t>& p, const Matrix<XX>& B, const size_t ks)
{
    const size_t m=A.row_size();
    size_t js=p[ks];
    Vector<XX> d(m);
    for(size_t k=0; k!=m; ++k) {
        for(size_t i=0; i!=m; ++i) {
            d[k]-=B[k][i]*A[i][js];
        }
    }
    return d;
}

template<> Interval inf<Interval>() { return Interval(inf<Float>()); }

template<class X>
pair<size_t,X>
compute_rt(const array<size_t>& p, const Vector<X>& x, const Vector<X>& d)
{
    const size_t m=d.size();
    X t=inf<X>();
    size_t r=m;
    for(size_t k=0; k!=m; ++k) {
        if(d[k] < -CUTOFF_THRESHOLD && x[p[k]] >= CUTOFF_THRESHOLD) {
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
    const X inf=Ariadne::inf<X>();

    // Choose variable to take out of basis
    // If the problem is degenerate, choose the variable with smallest index
    const size_t m=d.size();
    const size_t n=x.size();
    size_t r=n;
    X ds=(vt[p[s]]==LOWER ? +1 : -1);
    X t=u[p[s]]-l[p[s]];
    if(t<inf) { r=s; }
    X tk=0.0;
    ARIADNE_LOG(7,"   l="<<l<<" x="<<x<<" u="<<u<<"\n");
    ARIADNE_LOG(7,"   vt="<<vt<<" p="<<p<<" d="<<d<<"\n");
    ARIADNE_LOG(7,"   s="<<s<<" p[s]="<<p[s]<<" vt[p[s]]="<<vt[p[s]]<<" ds="<<ds<<" l[p[s]]="<<l[p[s]]<<" u[p[s]]="<<u[p[s]]<<" r="<<r<<" t[r]="<<t<<"\n");
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if( d[k]*ds<-CUTOFF_THRESHOLD && x[j]>=l[j] && l[j] != -inf) {
            tk=(l[j]-x[j])/(ds*d[k]);
            //if( r==n || tk<t || (tk==t && p[k]<p[r]) ) { t=tk; r=k; }
            if( tk<t || (tk==t && p[k]<p[r]) ) { t=tk; r=k; }
        } else if( d[k]*ds>CUTOFF_THRESHOLD && x[j]<=u[j] && u[j] != inf ) {
            tk=(u[j]-x[j])/(ds*d[k]);
            //if( r==n || tk<t || (tk==t && p[k]<p[r])) { t=tk; r=k; }
            if( tk<t || (tk==t && p[k]<p[r])) { t=tk; r=k; }
        } else {
            tk=inf;
        }
        ARIADNE_LOG(7,"    k="<<k<<" j=p[k]="<<j<<" l[j]="<<l[j]<<" x[j]="<<x[j]<<" u[j]="<<u[j]<<" d[k]="<<d[k]<<" t[k]="<<tk<<" r="<<r<<" t[r]="<<t<<"\n");
    }
    t*=ds;

    if(r==n) {
        // Problem is either highly degenerate or optimal do nothing.
        ARIADNE_WARN("SimplexSolver<X>::compute_rt(...): "<<
                     "Cannot find compute variable to exit basis\n"<<
                     "  l="<<l<<" x="<<x<<" u="<<u<<" vt="<<vt<<" p="<<p<<" d="<<d<<"\n");
    }
    return make_pair(r,t);
}

std::pair<size_t,Interval>
compute_rt(const Vector<Float>& l, const Vector<Float>& u, const array<Slackness>& vt, const array<size_t>& p, const Vector<Interval>& x, const Vector<Interval>& d, const size_t s)
{
    typedef Float X;
    typedef Interval XX;
    const X inf=Ariadne::inf<X>();

    // Choose variable to take out of basis
    // If the problem is degenerate, choose the variable with smallest index
    const size_t m=d.size();
    const size_t n=x.size();
    size_t r=n;
    X ds=(vt[p[s]]==LOWER ? +1 : -1);
    XX t=XX(u[p[s]])-XX(l[p[s]]);
    if(definitely(t<inf)) { r=s; }
    XX tk=0.0;
    ARIADNE_LOG(7,"   l="<<l<<" x="<<x<<" u="<<u<<"\n");
    ARIADNE_LOG(7,"   vt="<<vt<<" p="<<p<<" d="<<d<<"\n");
    ARIADNE_LOG(7,"   s="<<s<<" p[s]="<<p[s]<<" vt[p[s]]="<<vt[p[s]]<<" ds="<<ds<<" l[p[s]]="<<l[p[s]]<<" u[p[s]]="<<u[p[s]]<<" r="<<r<<" t[r]="<<t<<"\n");
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if( definitely(d[k]*ds<0.0) && definitely(x[j]>=l[j]) && l[j] != -inf) {
            tk=(l[j]-x[j])/(ds*d[k]);
            //if( r==n || tk<t || (tk==t && p[k]<p[r]) ) { t=tk; r=k; }
            if( tk<t || (tk==t && p[k]<p[r]) ) { t=tk; r=k; }
        } else if( definitely(d[k]*ds>0.0) && definitely(x[j]<=u[j]) && u[j] != inf ) {
            tk=(u[j]-x[j])/(ds*d[k]);
            //if( r==n || tk<t || (tk==t && p[k]<p[r])) { t=tk; r=k; }
            if( tk<t || (tk==t && p[k]<p[r])) { t=tk; r=k; }
        } else {
            tk=inf;
        }
        ARIADNE_LOG(7,"    k="<<k<<" j=p[k]="<<j<<" l[j]="<<l[j]<<" x[j]="<<x[j]<<" u[j]="<<u[j]<<" d[k]="<<d[k]<<" t[k]="<<tk<<" r="<<r<<" t[r]="<<t<<"\n");
    }
    t*=ds;

    if(r==n) {
        // Problem is either highly degenerate or optimal do nothing.
        ARIADNE_WARN("SimplexSolver<X>::compute_rt(...): "<<
                     "Cannot find compute variable to exit basis\n"<<
                     "  l="<<l<<" x="<<x<<" u="<<u<<" vt="<<vt<<" p="<<p<<" d="<<d<<"\n");
    }
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
void
update_y(const Vector<X>& l, const Vector<X>& u, const array<size_t>& p, Vector<X>& y, const size_t s, const Vector<X>& d, const X& t)
{
    ARIADNE_NOT_IMPLEMENTED;
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
tribool
SimplexSolver<X>::validated_constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u)
{
    ARIADNE_LOG(4,"A="<<A<<" b="<<b<<"\n");
    ARIADNE_LOG(4,"l="<<l<<" u="<<u<<"\n");

    array<size_t> p(A.column_size());
    array<Slackness> vt(A.column_size());
    Matrix<X> B(A.row_size(),A.row_size());
    make_lpair(p,B)=this->compute_basis(A);
    vt=compute_vt(l,u,p,A.row_size());

    bool done = false;
    while(!done) {
        done=this->validated_feasibility_step(A,b,l,u,vt,p);
    }
    return this->verify_constrained_feasibility(A,b,l,u,vt);
}

template<class X>
bool
SimplexSolver<X>::validated_feasibility_step(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<Slackness>& vt, array<size_t>& p)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    static const X inf = Ariadne::inf<X>();

    typedef Interval XX;

    ARIADNE_LOG(9,"vt="<<vt<<" p="<<p<<"\n");
    Matrix<XX> B=compute_B<XX>(A,p);
    ARIADNE_LOG(9," B="<<B<<"\n");
    Vector<XX> x=compute_x(A,b,l,u,vt,p,B);
    ARIADNE_LOG(9," x="<<x<<"\n");

    tribool feasible=true;

    Vector<X> c(n);
    Vector<X> ll(l);
    Vector<X> uu(u);
    for(uint i=0; i!=m; ++i) {
        size_t j=p[i];
        if(possibly(x[p[i]]<=l[p[i]])) { c[j]=-1; ll[j]=-inf; feasible=indeterminate; }
        if(possibly(x[p[i]]>=u[p[i]])) { c[j]=+1; uu[j]=+inf; feasible=indeterminate; }
    }
    ARIADNE_LOG(9," c="<<c<<"\n");
    if(definitely(feasible)) { return true; }

    const Vector<XX> y=compute_y(c,p,B);
    ARIADNE_LOG(9," y="<<y<<"\n");

    const Vector<XX> z=compute_z(A,c,p,y);
    ARIADNE_LOG(9," z="<<z<<"\n");

    size_t s = n;
    feasible=false;
    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        if(vt[j]==LOWER) { if(possibly(z[j]<=0)) { feasible=indeterminate; if(definitely(z[j]<0)) { s=k; break; } } }
        if(vt[j]==UPPER) { if(possibly(z[j]>=0)) { feasible=indeterminate; if(definitely(z[j]>0)) { s=k; break; } } }
    }
    ARIADNE_LOG(9," s="<<s<<"\n");
    if(definitely(!feasible)) { return true; }
    if(s==n) { ARIADNE_LOG(9," Cannot find variable to exit basis; no improvement can be made\n"); return true; }

    Vector<XX> d=compute_d(A,p,B,s);
    ARIADNE_LOG(9," d="<<d<<"\n");

    // Compute distance t along d in which to move,
    // and the variable p[r] to leave the basis
    // The bounds on t are given by l <= x + t * d <= u
    // Note that t is negative if an upper variable enters the basis
    size_t r; XX t;
    make_lpair(r,t)=compute_rt(l,u,vt,p,x,d,s);
    if(r==n) {
        ARIADNE_LOG(3,"   Cannot find variable to enter basis; no improvement can be made\n");
        return true;
    }

    ARIADNE_LOG(5,"  s="<<s<<" p[s]="<<p[s]<<" r="<<r<<" p[r]="<<p[r]<<" d="<<d<<" t="<<t<<"\n");

    if(r==s) {
        // Update variable type
        if(vt[p[s]]==LOWER) { vt[p[s]]=UPPER; }
        else { vt[p[s]]=LOWER; }
    } else {
        // Variable p[r] should leave basis, and variable p[s] enter
        ARIADNE_ASSERT(r<m);

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

    return false;

}



template<class X>
size_t
SimplexSolver<X>::lpstep(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<Slackness>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x, size_t s)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    ARIADNE_ASSERT(s<=n);
    ARIADNE_ASSERT(vt[p[s]]!=BASIS);

    // Compute direction d in which to move the current basic variables
    // as the variable entering the basis changes by +1
    Vector<X> d=compute_d(A,p,B,s);

    // Compute distance t along d in which to move,
    // and the variable p[r] to leave the basis
    // The bounds on t are given by l <= x + t * d <= u
    // Note that t is negative if an upper variable enters the basis
    size_t r; X t;
    make_lpair(r,t)=compute_rt(l,u,vt,p,x,d,s);
    if(r==n) {
        ARIADNE_LOG(3,"   Cannot find variable to enter basis; no improvement can be made\n");
        return r;
    }

    ARIADNE_LOG(5,"  s="<<s<<" p[s]="<<p[s]<<" r="<<r<<" p[r]="<<p[r]<<" d="<<d<<" t="<<t<<"\n");

    if(r==s) {
        Slackness nvts=(vt[p[s]]==LOWER ? UPPER : LOWER);
        ARIADNE_LOG(5,"   Changing non-basic variable x["<<p[s]<<"]=x[p["<<s<<"]] from type "<<vt[p[s]]<<" to type "<<nvts<<"\n");
    } else {
        ARIADNE_LOG(5,"   Swapping non-basic variable x["<<p[s]<<"]=x[p["<<s<<"]] with basic variable x["<<p[r]<<"]=x[p["<<r<<"]]\n");
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

    const double ERROR_TOLERANCE = std::numeric_limits<float>::epsilon();

    // Recompute B and x if it appears that there are problems with numerical degeneracy
    bool possible_degeneracy=false;
    for(uint i=0; i!=m; ++i) {
        if(l[p[i]]>x[p[i]] || x[p[i]]>u[p[i]]) {
            possible_degeneracy=true;
            break;
        }
    }
    B=compute_B<X>(A,p);
    x=compute_x<X>(A,b,l,u, vt,p,B);
    for(uint i=0; i!=m; ++i) {
        if(x[p[i]]<l[p[i]]) {
            ARIADNE_ASSERT(x[p[i]]>l[p[i]]-ERROR_TOLERANCE);
            x[p[i]]=l[p[i]];
        } else if(x[p[i]]>u[p[i]]) {
            ARIADNE_ASSERT(x[p[i]]<u[p[i]]+ERROR_TOLERANCE);
            x[p[i]]=u[p[i]];
        }
    }

    ARIADNE_LOG(7,"      vt="<<vt<<"\n      p="<<p<<"\n");
    ARIADNE_LOG(7,"      B="<<B<<"\n      x="<<x<<"\n");

    consistency_check(A,b,l,u,vt,p,B,x);

    return r;
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
    static const X inf = Ariadne::inf<X>();

    Vector<X> cc(n);
    Vector<X> ll(l);
    Vector<X> uu(u);

    // It seems that using this threshold does not work...
    static const double ROBUST_FEASIBILITY_THRESHOLD = std::numeric_limits<double>::epsilon() * 0;

    bool infeasible=false;
    for(size_t j=0; j!=n; ++j) {
        // If x[j] is (almost) infeasible by way of being to low, relax constraint x[j]>=l[j] to x[j]>=-inf.
        if(x[j]<l[j]) { cc[j]=-1; ll[j]=-inf; infeasible=true; }
        else if(x[j]>u[j]) { cc[j]=+1; uu[j]=+inf; infeasible=true; }
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
            if(x[j]<l[j]+ROBUST_FEASIBILITY_THRESHOLD) { cc[j]=-1; ll[j]=-inf; infeasible=true; }
            else if(x[j]>u[j]-ROBUST_FEASIBILITY_THRESHOLD) { cc[j]=+1; uu[j]=+inf; infeasible=true; }
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
        ARIADNE_ASSERT_MSG(x[i]>=l[i]-BOUNDS_TOLERANCE, "A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" B="<<B<<" x="<<x );
        ARIADNE_ASSERT_MSG(x[i]<=u[i]+BOUNDS_TOLERANCE, "A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<" B="<<B<<" x="<<x );
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
    Vector<X> y(m);

    return this->hotstarted_constrained_feasible(A,b,l,u,vt,p,B,x,y);
}



template<class X>
tribool
SimplexSolver<X>::hotstarted_constrained_feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, array<Slackness>& vt, array<size_t>& p, Matrix<X>& B, Vector<X>& x, Vector<X>& y)
{
    ARIADNE_LOG(5,"A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<"\n");
    ARIADNE_LOG(5,"vt="<<vt<<" p="<<p<<"\n");

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
    consistency_check(A,p,B);

    x=compute_x(A,b,l,u,vt,p,B);
    consistency_check(A,b,l,u,vt,p,B,x);

    tribool fs = _constrained_feasible(A,b,l,u,vt,p,B,x);

    ARIADNE_LOG(7,"vt="<<vt<<" p="<<p<<" fs="<<fs<<"\n");
    Vector<X> c=compute_c(l,u,p,x,m);
    y=compute_y(c,p,B);
    Vector<X> z=compute_z(A,c,p,y);
    ARIADNE_LOG(7,"x="<<x<<" c="<<c<<" y="<<y<<" z="<<z<<"\n");

    tribool vfs = verify_constrained_feasibility(A,b,l,u,vt);
    if(!indeterminate(vfs) && vfs!=fs) {
        if(verbosity>0) {
            ARIADNE_WARN("Approximate feasibility algorithm for\n  A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<"\nyielded basic variables "<<vt<<
                         " and result "<<fs<<", but validation code gave "<<vfs<<".\n");
        }
    }
    //FIXME: Return correct value
    // return vfs;
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

template<class X, class XX> Vector<X> compute_c(const size_t m, const Vector<X>& l, const Vector<X>& u, const array<size_t>& p, const Vector<XX>& x) {
    const size_t n=x.size();
    Vector<X> c(n);
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if((x[j]<=l[j])) { c[j]=-1; }
        else if((x[j]>=u[j])) { c[j]=+1; }
    }
    return c;
}

template<> Vector<Float> compute_c(const size_t m, const Vector<Float>& l, const Vector<Float>& u, const array<size_t>& p, const Vector<Interval>& x) {
    const size_t n=x.size();
    Vector<Float> c(n);
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if(possibly(x[j]<=l[j])) { c[j]=-1; }
        else if(possibly(x[j]>=u[j])) { c[j]=+1; }
        if(possibly(x[j]<l[j]) && possibly(x[j]>u[j])) {
            ARIADNE_FAIL_MSG("Unhandled case in checking feasibility of linear program. Basic variable x["<<j<<"]="<<x[j]<<" may violate both lower bound "<<l[j]<<" and upper bound "<<u[j]<<".");
        }
    }
    return c;
}



// A point x is strictly feasible for the basis B with lower variables L and upper variables U if
//   x_B = A_B^{-1} (b - A_L x_L - A_U x_U) is definitely within (x_B), the open interval (x_BL,x_BU).
// To prove infeasibility, choose a vector c (typically c_L=c_U=0, (c_B)_i = +1 if x_i>u_i and -1 if x_i<u_i)
// and set dual vector y = c_B A_B^{-1} .
//  y (b - A_L x_L - A_U x_U) > 0
//  z = c - y A  satisfies z_U < 0 and z_L > 0; by construction z_B = 0.
template<class X> tribool
SimplexSolver<X>::verify_constrained_feasibility(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u, const array<Slackness>& vt)
{
    ARIADNE_LOG(4,"verify_constrained_feasibility(Matrix A, Vector b, Vector l, Vector u, VariableTypeArray vt)\n");
    ARIADNE_LOG(5,"A="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<" vt="<<vt<<"\n");
    const array<size_t> p=compute_p(vt);

    typedef Interval XX;
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_ASSERT(b.size()==m);
    ARIADNE_ASSERT(l.size()==n);
    ARIADNE_ASSERT(u.size()==n);
    ARIADNE_ASSERT(vt.size()==n);

    // Ensure singleton constraints for x are non-basic
    for(size_t j=0; j!=n; ++j) {
        if(l[j]==u[j]) { ARIADNE_ASSERT(!vt[j]==BASIS);}
    }

    {
        const Matrix<X> B=compute_B<X>(A,p);
        const Vector<X> x=compute_x(A,b,l,u,vt,p,B);
        const Vector<X> c=compute_c(m,l,u,p,x);
        const Vector<X> y=compute_y(c,p,B);
        const Vector<X> z=compute_z(A,c,p,y);
        ARIADNE_LOG(7,"x="<<x<<" c="<<c<<" y="<<y<<" z="<<z<<"\n");
    }


    const Matrix<XX> B=compute_B<XX>(A,p);
    ARIADNE_LOG(9," B="<<B<<"; B*A="<<midpoint(prod(B,A))<<"\n");

    const Vector<XX> x=compute_x(A,b,l,u,vt,p,B);
    ARIADNE_LOG(9," x="<<x<<"; A*x="<<Vector<XX>(prod(A,x))<<"\n");


    const Vector<X> c=compute_c(m,l,u,p,x);
    ARIADNE_LOG(9," c="<<c<<"\n");

    const Vector<XX> y=compute_y(c,p,B);
    ARIADNE_LOG(9," y="<<y<<"\n");

    const Vector<XX> z=compute_z(A,c,p,y);
    ARIADNE_LOG(9," z="<<z<<"\n");

    ARIADNE_LOG(5,"x="<<x<<" c="<<c<<" y="<<y<<" z="<<z<<"\n");

    tribool fs=true;
    for(size_t k=0; k!=m; ++k) {
        size_t j=p[k];
        if(possibly(x[j]<=l[j]) || possibly(x[j]>=u[j])) {
            ARIADNE_LOG(9," k="<<k<<" j="<<j<<" l[j]="<<l[j]<<" x[j]="<<x[j]<<" u[j]="<<u[j]<<"\n");
            fs=indeterminate;
            if(definitely(x[j]<l[j]) || definitely(x[j]>u[j])) {
                fs=false;
                break;
            }
        }
    }

    if(fs==true) {
        return fs;
    }

    // The basis is optimal for min c x if z_L >=0  and z_U <= 0.
    // We have definite infeasibility if z_L > 0  and z_U < 0
    // If z_j < 0 for j lower, or z_j > 0 for j upper, then the simplex algorithm has not terminated correctly.

    for(size_t k=m; k!=n; ++k) {
        size_t j=p[k];
        if(vt[j]==LOWER && possibly(z[j]<=0)) { fs=indeterminate; break; }
        if(vt[j]==UPPER && possibly(z[j]>=0)) { fs=indeterminate; break; }
    }

    if(fs==false) {
        return fs;
    }

    // TODO: Code below is too enthusiastic about declaring an error.
    // For some degenerate, there may be no strictly feasible basic solution,
    // but it is easy to find a strictly feasible non-basic solution.
    // Should make feasibility testing more robust by testing the possibility of
    // obtaining strict feasibility by a partial step.

    return fs;
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
    this->consistency_check(A,p,B);

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

// Set r[i]=x[i]*y[i]+z for i=1,...,n
Vector<Float> efma(const Vector<Float>& x, const Vector<Float>& y, const Float& z) {
    Vector<Float> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=x[i]*y[i]+z;
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
Vector<Float> erec(const Vector<Float>& v) {
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

bool all_less(const Vector<Float>& x, const Float& e) {
    for(uint i=0; i!=x.size(); ++i) {
        if(x[i]>=e) { return false; }
    }
    return true;
}

bool all_greater(const Vector<Float>& x1, const Vector<Float>& x2) {
    for(uint i=0; i!=x1.size(); ++i) {
        if(x1[i]<=x2[i]) { return false; }
    }
    return true;
}



Interval dot(const Vector<Float>& c, const Vector<Interval>& x) {
    Interval r=0.0;
    for(uint i=0; i!=x.size(); ++i) {
        r+=c[i]*x[i];
    }
    return r;
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



Float compute_mu(const Vector<Float>& xl, const Vector<Float>& xu,
                 const Vector<Float>& x, const Vector<Float>& zl, const Vector<Float>& zu)
{
    const uint n=x.size();
    Float mu = 0.0;
    for(uint i=0; i!=n; ++i) {
        if(xl[i]!=-infty) { mu += ((x[i]-xl[i])*zl[i]); }
        if(xu[i]!=+infty) { mu += ((xu[i]-x[i])*zu[i]); }
    }
    mu /= (2.0*n);
    return mu;
}

Float compute_mu(const Vector<Float>& x, const Vector<Float>& z)
{
    const uint n=x.size();
    Float mu = 0.0;
    for(uint i=0; i!=n; ++i) {
        mu += x[i]*z[i];
    }
    mu /= n;
    return mu;
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
InteriorPointSolver::validate_feasibility(const Matrix<Float>& A, const Vector<Float>& b,
                                          const Vector<Float>& x, const Vector<Float>& y) const
{
    const uint m=A.row_size();
    const uint n=A.column_size();

    ARIADNE_LOG(2,"InteriorPointSolver::validate_feasibility(A,b,x,y)\n");
    ARIADNE_LOG(3,"A="<<A<<" b="<<b<<"\n");
    ARIADNE_LOG(3,"x="<<x<<" y="<<y<<"\n");
    ARIADNE_LOG(4,"Ax-b="<<Vector<Float>(A*x-b)<<" yA="<<y*A<<" yb="<<dot(y,b)<<"\n");

    // If yA < 0 and yb > 0, then problem is infeasible
    tribool result = false;

    set_rounding_upward();
    for(uint j=0; j!=x.size(); ++j) {
        Float yAj=0.0;
        for(uint i=0; i!=m; ++i) {
            yAj += y[i]*A[i][j];
        }
        if(yAj>=0.0) { result=indeterminate; break; }
    }
    set_rounding_downward();
    if(result==false) {
        Float yb=0.0;
        for(uint i=0; i!=y.size(); ++i) {
            yb += b[i]*y[i];
        }
        if(yb <= 0.0) { result=indeterminate; }
    }
    set_rounding_to_nearest();
    if(definitely(!result)) { return result; }

    // x should be an approximate solution to Ax=b
    // Use the fact that for any x, x'=(x + A^T (AA^T)^{-1}(b-Ax)) satisfies Ax'=0
    Vector<Interval> ivlx = x;
    Vector<Interval> ivle = b-A*ivlx;

    ARIADNE_LOG(4,"[e] = "<<ivle << "\n");
    Matrix<Interval> ivlS(m,m);
    for(uint i1=0; i1!=m; ++i1) {
        for(uint i2=i1; i2!=m; ++i2) {
            for(uint j=0; j!=n; ++j) {
                ivlS[i1][i2] += mul_ivl(A[i1][j],A[i2][j]);
            }
        }
    }
    for(uint i1=0; i1!=m; ++i1) {
        for(uint i2=0; i2!=i1; ++i2) {
            ivlS[i1][i2]=ivlS[i2][i1];
        }
    }

    Vector<Interval> ivld = solve(ivlS,ivle) * A;
    ARIADNE_LOG(4,"[d] = "<<ivld << "\n");

    ivlx = x + ivld;

    ARIADNE_LOG(3,"[x] = "<<ivlx << " A[x]-b="<<Vector<Interval>(A*ivlx-b)<<"\n");
    result=true;
    for(uint i=0; i!=n; ++i) {
        if(ivlx[i].lower()<=0.0) {
            result = indeterminate; break;
        }
    }

    return result;

}


tribool
InteriorPointSolver::validate_constrained_feasibility(const Matrix<Float>& A, const Vector<Float>& b,
                                                      const Vector<Float>& xl, const Vector<Float>& xu,
                                                      const Vector<Float>& x, const Vector<Float>& y) const
{
    const uint m=A.row_size();
    const uint n=A.column_size();

    // If yb - max(yA,0) xu + min(yA,0) xl > 0, then problem is infeasible
    // Evaluate lower bound for yb - max(z,0) xu + min(z,0) xl
    Vector<Interval> z=Vector<Interval>(y)*A;
    Float mx = 0.0;
    set_rounding_downward();
    for(uint i=0; i!=y.size(); ++i) {
        mx += b[i]*y[i];
    }
    for(uint i=0; i!=x.size(); ++i) {
        Float neg_xui = -xu[i];
        if(z[i].upper()>0.0) { mx += z[i].upper() * neg_xui; }
        if(z[i].lower()<0.0) { mx += z[i].lower() * xl[i]; }
    }
    set_rounding_to_nearest();
    if(mx>0.0) { return false; }

    // x should be an approximate solution to Ax=b
    // Use the fact that for any x, x'=(x + A^T (AA^T)^{-1}(b-Ax)) satisfies Ax'=0
    Vector<Interval> ivlx = x;
    Vector<Interval> ivle = b-A*ivlx;

    Matrix<Interval> ivlS(m,m);
    for(uint i1=0; i1!=m; ++i1) {
        for(uint i2=i1; i2!=m; ++i2) {
            for(uint j=0; j!=n; ++j) {
                ivlS[i1][i2] += mul_ivl(A[i1][j],A[i2][j]);
            }
        }
    }
    for(uint i1=0; i1!=m; ++i1) {
        for(uint i2=0; i2!=i1; ++i2) {
            ivlS[i1][i2]=ivlS[i2][i1];
        }
    }

    Vector<Interval> ivld = solve(ivlS,ivle) * A;

    ivlx += ivld;

    ARIADNE_LOG(2,"[x] = "<<ivlx);
    tribool result=true;
    for(uint i=0; i!=n; ++i) {
        if(ivlx[i].lower()<=xl[i] || ivlx[i].upper()>=xu[i]) {
            result = indeterminate; break;
        }
    }

    return result;

}


tuple< Vector<Float>, Vector<Float>, Vector<Float> >
InteriorPointSolver::optimize(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c) const
{
    Vector<Float> x,y,z;
    return this->hotstarted_optimize(A,b,c,x,y,z);
}


tuple< Vector<Float>, Vector<Float>, Vector<Float> >
InteriorPointSolver::hotstarted_optimize(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
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


tuple<Float, Vector<Float>, Vector<Float> >
InteriorPointSolver::constrained_optimize(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
                                          const Vector<Float>& xl, const Vector<Float>& xu) const
{
    ARIADNE_LOG(2,"InteriorPointSolver::constrained_optimize(A,b,c,xl,xu)\n");
    ARIADNE_LOG(2,"A="<<A<<", b="<<b<<", c="<<c<<"\n");
    ARIADNE_LOG(2,"xl="<<xl<<", xu="<<xu<<"\n");

    const uint m = b.size();
    const uint n = c.size();
    Vector<Float> y(m, 0.0);
    Vector<Float> x(n);
    Vector<Float> zl(n);
    Vector<Float> zu(n);
    for(uint i=0; i!=n; ++i) {
        if(xl[i]==-infty) {
            if(xu[i]==+infty) { x[i]=0.0; } else { x[i] = xu[i]-1.0; }
        } else {
            if(xu[i]==+infty) { x[i]=xl[i]+1.0; } else { ARIADNE_ASSERT(xl[i]<xu[i]); x[i] = (xl[i]+xu[i])/2; }
        }
        if(xl[i]==-infty) { zl[i] = 0.0; } else { zl[i] = 1.0; }
        if(xu[i]==+infty) { zu[i] = 0.0; } else { zu[i] = 1.0; }
    }

    bool done;
    do {
        done = this->_constrained_optimization_step(A,b,c,xl,xu, x,y,zl,zu);
    } while(!done);

    do {
        this->_constrained_optimization_step(A,b,c,xl,xu, x,y,zl,zu);
    } while(dot(c,x)-dot(y,b)>1e-4);

    // Todo: check for optimality
    return make_tuple(dot(c,x),x,y);
}




// Consider feasibility of  bl<=yB<=bu; yC<=c; dl<=y<=du
// Augment system by considering b<=yA-ve; yA+ve<=c; max(t)
tribool
InteriorPointSolver::primal_feasible(const Matrix<Float>& A, const Vector<Float>& b) const
{
    ARIADNE_LOG(2,"InteriorPointSolver::primal_feasible(A,b)\n");
    ARIADNE_LOG(3,"A="<<A<<" b="<<b<<"\n");
    Vector<Float> x(A.column_size(),1.0);
    Vector<Float> y(A.row_size(),0.0);
    Vector<Float> z(A.column_size(),1.0);

    const double THRESHOLD = 1e-6;

    while(true) {
        tribool result=this->_primal_feasibility_step(A,b,x,y,z);
        if(!indeterminate(result)) { return result; }
        if(compute_mu(x,z)<THRESHOLD) { return indeterminate; }
    }

}



tribool
InteriorPointSolver::constrained_feasible(const Matrix<Float>& A, const Vector<Float>& b,
                                          const Vector<Float>& xl, const Vector<Float>& xu) const
{
    ARIADNE_LOG(2,"InteriorPointSolver::constrained_feasible(A,b,xl,xu)\n");
    ARIADNE_LOG(2,"A="<<A<<", b="<<b<<"\n");
    ARIADNE_LOG(2,"xl="<<xl<<", xu="<<xu<<"\n");

    const uint m = A.row_size();
    const uint n = A.column_size();
    Vector<Float> c(n, 0.0);
    Vector<Float> y(m, 0.0);
    Vector<Float> x(n);
    Vector<Float> zl(n);
    Vector<Float> zu(n);
    for(uint i=0; i!=n; ++i) {
        if(xl[i]==-infty) {
            if(xu[i]==+infty) { x[i]=0.0; } else { x[i] = xu[i]-1.0; }
        } else {
            if(xu[i]==+infty) { x[i]=xl[i]+1.0; } else { ARIADNE_ASSERT(xl[i]<xu[i]); x[i] = (xl[i]+xu[i])/2; }
        }
        if(xl[i]==-infty) { zl[i] = 0.0; } else { zl[i] = 1.0; }
        if(xu[i]==+infty) { zu[i] = 0.0; } else { zu[i] = 1.0; }
    }

    const double THRESHOLD = 1e-6;
    while(true) {
        tribool result=this->_constrained_feasibility_step(A,b,xl,xu, x,y,zl,zu);
        if(!indeterminate(result)) { return result; }
        if(compute_mu(xl,xu, x,zl,zu)<THRESHOLD ) { return indeterminate; }
    }

}

tribool
InteriorPointSolver::constrained_dual_feasible(const Matrix<Float>& A, const Vector<Float>& cl, const Vector<Float>& cu,
                                               const Vector<Float>& l, const Vector<Float>& u) const
{
    ARIADNE_NOT_IMPLEMENTED;
}









tribool
InteriorPointSolver::_optimization_step(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
                                        Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const
{
    ARIADNE_LOG(4,"InteriorPointSolver::_optimization_step(A,b,c, x,y,z)\n");
    ARIADNE_LOG(5,"x="<<x<<" y="<<y<<" z="<<z<<"\n");
    static const double gamma=1.0/1024;
    static const double sigma=1.0/8;
    static const double scale=0.75;

    const uint m=A.row_size();
    const uint n=A.column_size();

    Vector<Float> dx(n),dy(n),dz(n);
    Vector<Float> nx(n),ny(n),nz(n),nxz(n);
    Vector<Float> rx(m),ry(n),rz(n),r(n);
    Matrix<Float> S(m,m),Sinv(m,m);

    Float mu=sigma*dot(x,z)/n; // duality constant
    ARIADNE_LOG(5,"mu="<<mu<<"\n");

    // rx = Ax-b; ry=yA+z-c; rz=x.z-mu
    rx=A*x-b;
    ry=y*A+z-c;
    rz=esub(emul(x,z),mu);
    ARIADNE_LOG(5,"rx="<<rx<<" ry="<<ry<<" rz="<<rz<<"\n");

    r=rx-A*(ediv(rz-emul(x,ry),z));

    S=mdtmul(A,ediv(x,z));
    Sinv=inverse(S);
    ARIADNE_LOG(5,"S="<<S<<"  inverse(S)="<<Sinv<<"\n");

    // Adx = rx; ATdy + dz = ry; Zdx + X dz = rz
    // dx = (rz - X dz)/Z
    // dz = ry - AT dy
    // A(X/Z)AT dy = rx - A Zinv (rz - X ry)
    dy=Sinv*r;
    dz=ry-dy*A;
    dx=ediv(rz-emul(x,dz),z);
    ARIADNE_LOG(5,"dx="<<dx<<" dy="<<dy<<" dz="<<dz<<"\n");

    // Try to enforce feasibility or dual feasibility
    Float alphax=1.0;
    nx=x-dx;
    while ( !all_greater(emul(nx,z),gamma*mu) ) {
        alphax=alphax*scale;
        nx=(x-alphax*dx);
    }

    Float alphaz=1.0;
    nz=z-dz;
    while ( !all_greater(emul(x,nz),gamma*mu) ) {
        alphaz=alphaz*scale;
        nz=(z-alphaz*dz);
    }
    ny=(y-alphaz*dy);
    ARIADNE_LOG(6,"alphax="<<alphax<<" nx="<<nx<<" alphaz="<<alphaz<<" ny="<<ny<<" nz="<<nz<<"\n");

    x=nx; y=ny; z=nz;
    ARIADNE_LOG(5,"cx="<<dot(c,x)<<" yb="<<dot(y,b)<<" Ax-b="<<Vector<Float>(prod(A,x)-b)<<" yA+z-c="<<Vector<Float>(y*A+z-c)<<"\n");

    if(alphax==1.0) { return true; }
    else if(alphaz==1.0 && dot(y,b) > 0.0) { return false; }
    else { return indeterminate; }
}


tribool
InteriorPointSolver::_constrained_optimization_step(const Matrix<Float>& A, const Vector<Float>& b, const Vector<Float>& c,
                                                    const Vector<Float>& xl, const Vector<Float>& xu,
                                                    Vector<Float>& x, Vector<Float>& y, Vector<Float>& zl, Vector<Float>& zu) const
{
    ARIADNE_LOG(4,"x="<<x<<", y="<<y<<", zl="<<zl<<", zu="<<zu<<"\n");

    static const double gamma=1.0/256;
    static const double sigma=1.0/8;
    static const double scale=0.75;


    const uint m=A.row_size();
    const uint n=A.column_size();

    Vector<Float> dx(n),dy(n),dzl(n),dzu(n);
    Vector<Float> nx,ny,nzl,nzu;
    Vector<Float> rx(m),ry(n),rzl(n),rzu(n);
    Matrix<Float> S(m,m),Sinv(m,m);
    DiagonalMatrix<Float> Xl(xl), Xu(xu), X(x), Zl(zl), Zu(zu);

    Float mu = compute_mu(xl,xu, x,zl,zu) * sigma;
    ARIADNE_LOG(4,"mu="<<mu<<"\n");

    // rx = Ax-b; ry=yA+zl-zu-c; rzl=(x-xl).zl-mu; rzu=(xu-x).zu-mu.
    rx=A*x-b;
    ry=y*A+(zl-zu)-c;
    for(uint i=0; i!=n; ++i) {
        if(xl[i]!=-infty) { rzl[i] = (x[i]-xl[i])*zl[i] - mu; }
        if(xu[i]!=+infty) { rzu[i] = (xu[i]-x[i])*zu[i] - mu; }
    }
    ARIADNE_LOG(5,"rx="<<rx<<", ry="<<ry<<", rzl="<<rzl<<", rzu="<<rzu<<"\n");

    // A dx = rx;  AT dy + dzl - dzu = ry;  Zl dx + (X-Xl) dzl = rzl; -Zu dx + (Xu-X) dzu = rzu;

    //   ryz = ( ry - rzl/(X-Xl) + rzu/(Xu-X) )
    //   D = 1/(Zu/(Xu-X) + Zl/(X-Xl))
    DiagonalMatrix<Float> D(erec(ediv(zu,xu-x)+ediv(zl,x-xl)));
    Vector<Float> ryz = ry - ediv(rzl,Vector<Float>(x-xl)) + ediv(rzu,Vector<Float>(xu-x));
    S=mdtmul(A,D.diagonal());
    ARIADNE_LOG(5,"S="<<S<<"  inverse(S)="<<Sinv<<"\n");

    // dzl = (rzl - Zl dx) / (X-Xl)
    // dzu = (rzu + Zu dx) / (Xu-X)
    // dx = D (AT dy - ryz)
    // (A D AT) dy = rx - A D ryz

    dy = solve(S, Vector<Float>( rx - A * (D * ryz) ) );
    dx = D * Vector<Float>(dy*A - ryz);
    dzl = Vector<Float>(rzl-Zl*dx)/(X-Xl);
    dzu = Vector<Float>(rzu+Zu*dx)/(Xu-X);
    ARIADNE_LOG(5,"dx="<<dx<<" dy="<<dy<<" dzl="<<dzl<<" dzu="<<dzu<<"\n");

    // Try to enforce feasibility or dual feasibility
    Float alphax=1.0;
    nx=x-dx;
    while ( !all_greater(emul(nx-xl,zl),gamma*mu) || !all_greater(emul(xu-nx,zu),gamma*mu) ) {
        alphax=alphax*scale;
        nx=(x-alphax*dx);
    }

    Float alphaz=1.0;
    nzl=zl-dzl; nzu=zu-dzu;
    while ( !all_greater(emul(nx-xl,nzl),gamma*mu) || !all_greater(emul(xu-nx,nzu),gamma*mu) ) {
        alphaz=alphaz*scale;
        nzl=(zl-alphaz*dzl);
        nzu=(zu-alphaz*dzu);
    }
    ny=(y-alphaz*dy);
    ARIADNE_LOG(7,"alphax="<<alphax<<" nx="<<nx<<" alphaz="<<alphaz<<" ny="<<ny<<" nzl="<<nzl<<" nzu="<<nzu<<"\n");

    x=nx; y=ny; zl=nzl; zu=nzu;
    ARIADNE_LOG(5,"cx="<<dot(c,x)<<" yb="<<dot(y,b)<<" Ax-b="<<Vector<Float>(prod(A,x)-b)<<" yA+(zl-zu)-c="<<Vector<Float>(y*A+(zl-zu)-c)<<"\n");

    if(alphax==1.0) { return true; }
    else if(alphaz==1.0 && dot(y,b)>dot(xl,zl)+dot(xu,zu)) { return false; }
    else { return indeterminate; }
}



tribool
InteriorPointSolver::_primal_feasibility_step(const Matrix<Float>& A, const Vector<Float>& b,
                                              Vector<Float>& x, Vector<Float>& y, Vector<Float>& z) const
{
    Vector<Float> c(x.size(),0.0);
    return this->_optimization_step(A,b,c, x,y,z);

    ARIADNE_LOG(4,"primal_feasibility_step(A,b,x,y)\n");
    ARIADNE_LOG(5,"x="<<x<<" y="<<y<<" z="<<z<<"\n");
    static const double gamma=1.0/1024;
    static const double sigma=1.0/8;
    static const double scale=0.75;

    const uint m=A.row_size();
    const uint n=A.column_size();

    Vector<Float> dx(n),dy(n),dz(n);
    Vector<Float> nx(n),ny(n),nz(n);
    Vector<Float> rx(m),ry(n),rz(n);
    Matrix<Float> S(m,m),Sinv(m,m);

    DiagonalMatrix<Float> Z(z);
    DiagonalMatrix<Float> X(x);

    Float mu = dot(x,z)/n * sigma;
    ARIADNE_LOG(5,"mu="<<mu<<"\n");

    rx = y*A+z;
    ry = A*x-b;
    rz = efma(x,z,-mu);
    ARIADNE_LOG(5,"rx=yA+z="<<rx<<" ry=Ax-b="<<ry<<" rz=x.z-mu="<<rz<<"\n");

    S=mdtmul(A,ediv(x,z));
    Sinv=inverse(S);

    ARIADNE_LOG(5,"S="<<S<<"  inverse(S)="<<Sinv<<"\n");

    // S dy = ry - A Zinv (rz - X rx)
    // dz = rx - AT dy
    // dx = Zinv * (rz - X dz)
    dy = Sinv * Vector<Float>( ry - A * Vector<Float>(rz-X*rx)/Z );
    dz = rx - dy * A;
    dx = Vector<Float>(rz-X*dz)/Z;
    ARIADNE_LOG(4,"  dx="<<dx<<" dy="<<dy<<" dz="<<dz<<"\n");

    // Check for feasibility
    nx=x-dx;
    if(all_greater(nx,0.0)) {
        ARIADNE_LOG(4,"nx="<<nx<<" Ax-b="<<Vector<Float>(A*nx-b)<<"\n");
        x=nx; return true;
    }

    // Check for infeasibility
    nz=z-dz; ny=y-dy;
    if(all_greater(nz,0.0) && dot(y,b)>0.0) {
        ARIADNE_LOG(4,"ny="<<ny<<" nz="<<nz<<" yA="<<Vector<Float>(ny*A)<<" yb="<<dot(ny,b)<<" yA+z="<<Vector<Float>(ny*A+nz)<<"\n");
        x=nx; y=ny; z=nz; return false;
    }

    Float alpha=1.0;
    while ( !all_greater(nx,0.0) || !all_greater(emul(nx,nz),gamma*mu) ) {
        alpha *= scale;
        nx=(x-alpha*dx);
        nz=(z-alpha*dz);
        ARIADNE_LOG(6,"    alpha="<<alpha<<" nx="<<nx<<" nz="<<nz<<"\n");
    }
    ny=(y-alpha*dy);

    x=nx; y=ny; z=nz;

    ARIADNE_LOG(4,"nx="<<x<<" ny="<<y<<" nz="<<z<<"\n");
    ARIADNE_LOG(5,"yb="<<dot(y,b)<<" yA="<<y*A<<" Ax-b="<<Vector<Float>(A*x-b)<<" yA+z="<<Vector<Float>(y*A+z)<<" x.z="<<emul(x,z)<<"\n");

    return indeterminate;
}



tribool
InteriorPointSolver::_constrained_feasibility_step(const Matrix<Float>& A, const Vector<Float>& b,
                                                   const Vector<Float>& xl, const Vector<Float>& xu,
                                                   Vector<Float>& x, Vector<Float>& y, Vector<Float>& zl, Vector<Float>& zu) const
{
    Vector<Float> c(x.size(),0.0);
    return this->_constrained_optimization_step(A,b,c,xl,xu,x,y,zl,zu);


    ARIADNE_LOG(4,"InteriorPointSolver::constrained_feasibility_step(A,b,xl,xu, x,y,zl,zu)\n");
    ARIADNE_LOG(5,"x="<<x<<" y="<<y<<" zl="<<zl<<" zu="<<zu<<"\n");
    static const double gamma=1.0/1024;
    static const double sigma=1.0/8;
    static const double scale=0.75;

    const uint m=A.row_size();
    const uint n=A.column_size();

    Vector<Float> wl(n),wu(n);
    Vector<Float> dx(n),dy(n),dz(n),dzl(n),dzu(n);
    Vector<Float> nx(n),ny(n),nzl(n),nzu(n);
    Vector<Float> rx(m),ry(n),rzl(n),rzu(n);
    Matrix<Float> S(m,m),Sinv(m,m);

    DiagonalMatrix<Float> Xl(xl), Xu(xu), X(x), Zl(zl), Zu(zu);

    // Compute central path parameter mu
    Float mu = 0.0;
    for(uint i=0; i!=n; ++i) {
        if(xl[i]!=-infty) { mu += (x[i]-xl[i])*zl[i]; }
        if(xu[i]!=+infty) { mu += (xu[i]-x[i])*zu[i]; }
    }
    mu *= ( sigma / (2*n) );
    ARIADNE_LOG(5,"mu="<<mu<<"\n");

    rx = y*A+(zl-zu);
    ry = A*x-b;
    rzl = efma(x-xl,zl,-mu);
    rzu = efma(xu-x,zu,-mu);
    //Vector<Float> ryz = ry - ediv(rzl,x-xl) + ediv(rzu,xu-x);
    ARIADNE_LOG(5,"rx=yA+z="<<rx<<" ry=Ax-b="<<ry<<" rzl=(x-xl).zl-mu="<<rzl<<" rzu=(xu-x).zu-mu="<<rzu<<"\n");

    DiagonalMatrix<Float> D(erec(ediv(wl,zl)+ediv(wu,zu)));
    S=mdtmul(A,D.diagonal());
    Sinv=inverse(S);

    ARIADNE_LOG(5,"S="<<S<<"  inverse(S)="<<Sinv<<"\n");

    // Ax-b=rx=0; yA+zl-zu=ry=0; (x-xl).zl-mu=rzl=0; (xu-x).zu-mu=rzu=0
    //
    // dzl = (rzl - Zl dx) / (X-Xl)
    // dzu = (rzu + Zu dx) / (Xu-X)
    // (Zu/(Xu-X) + Zl/(X-Xl)) dx = AT dy - ( ry - rzl/(X-Xl) + rzu/(Xu-X) )
    // (A D AT) dy = rx - A D ( ry - rzl/(X-Xl) + rzu/(Xu-X) )

    //dy = Sinv * Vector<Float>( rx - A * (D * ryz) );
    //dx = D * (dy*A - ryz);
    //dzl = Vector<Float>(rzl-Zl*dx)/(X-Xl);
    //dzu = Vector<Float>(rzu-Zu*dx)/(Xu-X);
    ARIADNE_LOG(4,"  dx="<<dx<<" dy="<<dy<<" dzl="<<dzl<<" dzu="<<dzu<<"\n");

    // Check for feasibility
    nx=x-dx;
    if(all_greater(nx,0.0)) {
        ARIADNE_LOG(4,"nx="<<nx<<" Ax-b="<<Vector<Float>(A*nx-b)<<"\n");
        x=nx; return true;
    }

    // Check for infeasibility
    nzl=zl-dzl; nzu=zu-dzu; ny=y-dy;
    if(all_greater(nzl,0.0) && all_greater(nzu,0.0) && dot(y,b)>0.0) {
        ARIADNE_LOG(4,"ny="<<ny<<" nzl="<<nzl<<" nzu="<<nzu<<" yA="<<Vector<Float>(ny*A)<<" yb="<<dot(ny,b)<<" yA+zl-zu="<<Vector<Float>(ny*A+nzl-nzu)<<"\n");
        x=nx; y=ny; zl=nzl; zu=nzu; return false;
    }

    Float alpha=1.0;
    while ( !all_greater(nx,0.0) || !all_greater(emul(nx-xl,nzl),gamma*mu) || !all_greater(emul(xu-nx,nzu),gamma*mu) ) {
        alpha *= scale;
        nx=(x-alpha*dx);
        nzl=(zl-alpha*dzl);
        nzu=(zu-alpha*dzu);
        ARIADNE_LOG(6,"    alpha="<<alpha<<" nx="<<nx<<" nzl="<<nzl<<" nzu="<<nzu<<"\n");
    }
    ny=(y-alpha*dy);

    x=nx; y=ny; zl=nzl; zu=nzu;

    ARIADNE_LOG(4,"nx="<<x<<" ny="<<y<<" nzl="<<zl<<" nzu="<<zu<<"\n");
    ARIADNE_LOG(5,"yb="<<dot(y,b)<<" yA="<<y*A<<" Ax-b="<<Vector<Float>(A*x-b)<<" yA+zl-zu="<<Vector<Float>(y*A+zl-zu)<<
                  " (x-xl).zl="<<emul(x-xl,zl)<<" (xu-x).zu"<<emul(xu-x,zu)<<"\n");

    return indeterminate;
}



} // namespace Ariadne

