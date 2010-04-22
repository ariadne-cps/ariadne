/***************************************************************************
 *            matrix.cc
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

#include "rounding.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"

template class boost::numeric::ublas::matrix<Ariadne::Float>;
template class boost::numeric::ublas::matrix<Ariadne::Interval>;

#ifdef HAVE_RATIONAL
template class boost::numeric::ublas::matrix<Ariadne::Rational>;
#endif // HAVE_RATIONAL


namespace Ariadne {

Matrix<Float>
midpoint(const Matrix<Interval>& A) {
    Matrix<Float> R(A.row_size(),A.column_size());
    for(size_t i=0; i!=A.row_size(); ++i) {
        for(size_t j=0; j!=A.column_size(); ++j) {
            R[i][j]=A[i][j].midpoint();
        }
    }
    return R;
}

template<class X>
Matrix<X>
lu_inverse(const Matrix<X>& M)
{
    typedef X RealType;
    ARIADNE_ASSERT_MSG(M.row_size()==M.column_size(),"A="<<M);

    size_t m=M.row_size();
    size_t n=M.column_size();
    Matrix<RealType> A=M;
    Matrix<RealType> B=Matrix<RealType>::identity(n);

    // Array of row pivots. The value p[i] gives the row
    // swapped with the ith row in the ith stage.
    array<size_t> p(m);
    for(size_t k=0; k!=m; ++k) { p[k]=k; }


    for(size_t k=0; k!=std::min(m,n); ++k) {
        // Choose a pivot row
        size_t iamax=k;
        X amax=0;
        for(size_t i=k; i!=m; ++i) {
            if(abs(A[i][k])>amax) {
                iamax=i;
                amax=abs(A[i][k]);
            }
        }

        // Set pivot row
        size_t i=iamax;
        p[k]=i;

        // Swap rows of both A and B
        for(size_t j=k; j!=n; ++j) {
            X tmp=A[k][j];
            A[k][j]=A[i][j];
            A[i][j]=tmp;
        }

        for(size_t j=0; j!=n; ++j) {
            X tmp=B[k][j];
            B[k][j]=B[i][j];
            B[i][j]=tmp;
        }

        RealType r  = 1/A[k][k];
        for(size_t i=k+1; i!=n; ++i) {
            RealType s=A[i][k] * r;
            for(size_t j=0; j!=n; ++j) {
                B[i][j] -= s * B[k][j];
            }
            for(size_t j=k+1; j!=n; ++j) {
                A[i][j] -= s * A[k][j];
            }
            A[i][k] = 0;
        }
        for(size_t j=0; j!=n; ++j) {
            B[k][j] *= r;
        }
        for(size_t j=k+1; j!=n; ++j) {
            A[k][j] *= r;
        }
        A[k][k] = 1;
    }

    // Backsubstitute to find inverse
    for(size_t k=n; k!=0; ) {
        --k;
        for(size_t i=0; i!=k; ++i) {
            RealType s=A[i][k];
            for(size_t j=0; j!=n; ++j) {
                B[i][j] -= s * B[k][j];
            }
            A[i][k] = 0;
        }
    }

    // No need to repivot!

    return B;
}

template<class X> Matrix<X> lu_solve(const Matrix<X>& A, const Matrix<X>& B);
template<class X> Matrix<X> gs_solve(const Matrix<X>& A, const Matrix<X>& B);
template<class X> Matrix<X> dd_solve(const Matrix<X>& A, const Matrix<X>& B);

template<class X> Matrix<X> lu_solve(const Matrix<X>& A, const Matrix<X>& B) {
    return prod(lu_inverse(A),B);
}

// Find a starting solution for a diagonally dominant system
Matrix<Interval> dd_solve(const Matrix<Interval>& A, const Matrix<Interval>& B)
{
    const size_t n=B.row_size();
    const size_t m=B.column_size();

    //Compute 1/(|aii|-sum|aij|) using outward rounding
    rounding_mode_t rounding_mode=get_rounding_mode();
    set_rounding_mode(upward);
    Vector<Float> c(n,0.0);
    for(size_t i=0; i!=n; ++i) {
        volatile double& ci=c[i];
        for(size_t j=0; j!=n; ++j) {
            if(j!=i) {
                ci+=mag(A[i][j]);
            }
        }
        ci=ci-mig(A[i][i]);
        ci=-ci;
        if(ci<=0.0) {
            ARIADNE_THROW(std::runtime_error,"dd_solve(Matrix<Interval> A, Matrix<Interval> B)",
                          "Interval matrix A="<<A<<" is not diagonally-dominant.");
        }
        ci=1.0/ci;
    }
    set_rounding_mode(rounding_mode);

    // Compute initial solution
    Matrix<Interval> R(n,m);
    for(size_t i=0; i!=n; ++i) {
        Interval ci(-c[i],+c[i]);
        for(size_t j=0; j!=m; ++j) {
            R[i][j]=B[i][j]*ci;
        }
    }

    return R;
}


template<> Matrix<Interval> gs_solve(const Matrix<Interval>& A, const Matrix<Interval>& B)
{
    ARIADNE_ASSERT(A.row_size()==A.column_size());
    ARIADNE_ASSERT(B.row_size()==A.column_size());
    const size_t n=B.row_size();
    const size_t m=B.column_size();

    // Precondition A and B
    Matrix<Float> J;
    try {
        J=inverse(midpoint(A));
    } catch(const SingularMatrixException& e) {
        ARIADNE_FAIL_MSG("SingularMatrixException catche");
    }
    Matrix<Interval> JA=prod(J,A);
    Matrix<Interval> JB=prod(J,B);
    //std::cerr<<"J="<<J<<"\nJA="<<JA<<"\nJB="<<JB<<"\n";

    //Matrix<Interval> R=dd_solve(JA,JB);
    Matrix<Interval> R=lu_inverse(A);
    //std::cerr<<"R="<<R<<"\n";

    // Perform Gauss-Seidel iteration
    size_t step=0;
    while(step<1) {
        for(size_t i=0; i!=n; ++i) {
            for(size_t j=0; j!=m; ++j) {
                // compute R'[i][j] := (JB[i][j] - Sum{k!=j}JA[i][k]*R[k][j]) / JA[i][i]
                Interval Rij=JB[i][j];
                for(size_t k=0; k!=n; ++k) {
                    if(k!=i) {
                        Rij-=JA[i][k]*R[k][j];
                    }
                }
                Rij/=JA[i][i];
                R[i][j]=intersection(R[i][j],Rij);
            }
        }
        ++step;
        //std::cerr<<"R="<<R<<"\n";
    }

    return R;
}

template<class X> Matrix<X> gs_inverse(const Matrix<X>& A) {
    return gs_solve(A,Matrix<X>::identity(A.row_size()));
}


template<class X> Matrix<X> inverse(const Matrix<X>& A) {
    return lu_inverse(A);
}

template<> Matrix<Interval> inverse<Interval>(const Matrix<Interval>& A) {
    try {
        return lu_inverse(A);
    } catch(const DivideByZeroException& e) {
        throw SingularMatrixException();
    }
}


template<class X> Matrix<X> solve(const Matrix<X>& A, const Matrix<X>& B) {
    try {
        return prod(inverse(A),B);
    } catch(const SingularMatrixException& e) {
        ARIADNE_FAIL_MSG("SingularMatrixException");
    }
}

template<> Matrix<Interval> solve(const Matrix<Interval>& A, const Matrix<Interval>& B) {
    return gs_solve(A,B);
}

template<class X> Vector<X> solve(const Matrix<X>& A, const Vector<X>& b) {
    return prod(lu_inverse(A),b);
}


template<class X>
Matrix<X>
solve(const PLUMatrix<X>& A, const Matrix<X>& B)
{
    const PivotMatrix& P=A.P;
    const Matrix<X>& L=A.L;
    const Matrix<X>& U=A.U;

    // Only work with square matrices L, U
    ARIADNE_ASSERT(L.row_size()==L.column_size());
    ARIADNE_ASSERT(U.row_size()==U.column_size());

    // Check sizes for consistency
    ARIADNE_ASSERT(L.column_size()==U.row_size());
    ARIADNE_ASSERT(P.size()==L.row_size());
    ARIADNE_ASSERT(B.row_size()==P.size());

    // Set constants for sizes
    const size_t m=B.column_size();
    const size_t n=B.row_size();

    Matrix<X> R=B;

    // PivotMatrix rows of B
    for(size_t k=0; k!=n; ++k) {
        size_t i=P[k];
        for(size_t j=0; j!=m; ++j) {
            X tmp=R[i][j];
            R[i][k]=R[k][j];
            R[k][j]=tmp;
        }
    }

    // Backsubstitute on L
    for(size_t k=0; k!=n; ++k) {
        for(size_t i=k+1; i!=n; ++i) {
            X s=L[i][k];
            for(size_t j=0; j!=n; ++j) {
                R[i][j] -= s * R[k][j];
            }
        }
    }

    // Backsubstitute on U with row scaling
    for(size_t k=n; k!=0; ) {
        --k;
        X s=1/U[k][k];
        for(size_t j=0; j!=m; ++j) {
            R[k][j] *= s;
        }
        for(size_t i=0; i!=k; ++i) {
            X s=U[i][k];
            for(size_t j=0; j!=m; ++j) {
                R[i][j] -= s * R[k][j];
            }
        }
    }

}


// Compute the vector of row norms (sums of absolute values) of A
Vector<Float>
row_norms(const Matrix<Float>& A)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    Vector<Float> r(m);

    rounding_mode_t prev_rounding_mode=get_rounding_mode();
    set_rounding_mode(upward);
    for(size_t i=0; i!=m; ++i) {
        r[i]=0.0;
        for(size_t j=0; j!=n; ++j) {
            r[i]+=abs(A[i][j]);
        }
    }
    set_rounding_mode(prev_rounding_mode);

    return r;

}


// Scale the rows of A to each have some of absolute values equal to one
Matrix<Float>
normalise_rows(const Matrix<Float>& A)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Matrix<Float> R=A;

    array<Float> row_asums(m);
    rounding_mode_t prev_rounding_mode=get_rounding_mode();
    set_rounding_mode(upward);
    for(size_t i=0; i!=m; ++i) {
        row_asums[i]=0.0;
        for(size_t j=0; j!=n; ++j) {
            row_asums[i]+=abs(A[i][j]);
        }
    }
    set_rounding_mode(toward_zero);
    for(size_t i=0; i!=m; ++i) {
        for(size_t j=0; j!=n; ++j) {
            R[i][j]/=row_asums[i];
        }
    }
    set_rounding_mode(prev_rounding_mode);

    return R;
}


// Returns a pivot P and matrices L and U such that L is unit lower-triangular,
// U is upper-trianguler and A=PLU.
tuple< PivotMatrix, Matrix<Float>, Matrix<Float> >
triangular_decomposition(const Matrix<Float>& A)
{
    typedef Float RealType;
    ARIADNE_ASSERT(A.row_size()==A.column_size());

    size_t m=A.row_size();
    size_t n=A.column_size();
    Matrix<RealType> L=Matrix<RealType>::identity(m);
    Matrix<RealType> U=A;

    // Array of row pivots. The value P[i] gives the row
    // swapped with the ith row in the ith stage.
    PivotMatrix P(m);
    for(size_t k=0; k!=m; ++k) { P[k]=k; }

    for(size_t k=0; k!=std::min(m,n); ++k) {
        // Choose a pivot row
        size_t iamax=k;
        RealType amax=0;
        for(size_t i=k; i!=m; ++i) {
            if(abs(A[i][k])>amax) {
                iamax=i;
                amax=abs(A[i][k]);
            }
        }

        // Set pivot row
        size_t l=iamax;
        P[k]=l;

        // Swap rows of L and U
        for(size_t j=0; j!=k; ++j) {
            RealType tmp=L[k][j];
            L[k][j]=L[l][j];
            L[l][j]=tmp;
        }

        for(size_t j=k; j!=n; ++j) {
            RealType tmp=U[k][j];
            U[k][j]=U[l][j];
            U[l][j]=tmp;
        }

        RealType r  = 1/U[k][k];
        for(size_t i=k+1; i!=m; ++i) {
            RealType s=U[i][k] * r;
            for(size_t j=k+1; j!=n; ++j) {
                U[i][j] -= s * U[k][j];
            }
            U[i][k] = 0;
            L[i][k] = s;
        }

    }

    return make_tuple(P,L,U);

}

// Returns a matrix R such that A=OR where O is a square matrix with
// orthogonal rows and the sum of the absolute values of the rows of R are
// at most one.
//
// Use Householder transformation H=I-vv' where v=u/|u|
// and u=a +/- |a|e with a and e the working column of A
// and corresponding unit vector.
tuple< Matrix<Float>, PivotMatrix>
triangular_factor(const Matrix<Float>& A)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Matrix<Float> R=A;
    PivotMatrix P(n); for(uint i=0; i!=n; ++i) { P[i]=i; }

    // Array of column norm squares
    array<Float> cns(n);

    Vector<Float> u(m);

    for(size_t k=0; k!=min(m,n); ++k) {
        //std::cerr<<"k="<<k<<" R="<<R<<std::endl;

        bool pivoting=true;
        if(pivoting) {
            // Compute column norms
            for(size_t j=k; j!=n; ++j) {
                cns[j]=0.0;
                for(size_t i=k; i!=m; ++i) {
                    cns[j]+=sqr(R[i][j]);
                }
            }

            // Find largest column norm
            size_t l=k; Float cnsmax=cns[k];
            for(size_t j=k+1; j!=n; ++j) {
                if(cns[j]>cnsmax) {
                    l=j; cnsmax=cns[j];
                }
            }

            // Swap columns l and k
            for(size_t i=0; i!=m; ++i) {
                Float tmp=R[i][k];
                R[i][k]=R[i][l];
                R[i][l]=tmp;
            }

            // Set pivot element to pivot column
            P[k]=l;
        }

        // Compute |a| where a is the working column
        Float nrmas=0.0;
        for(size_t i=k; i!=m; ++i) {
            nrmas+=R[i][k]*R[i][k];
        }
        Float nrma=sqrt(nrmas);

        // Compute u=a +/- |a|e
        for(size_t i=k; i!=m; ++i) {
            u[i]=R[i][k];
        }
        if(u[k]>=0) { u[k]+=nrma; }
        else { u[k]-=nrma; }

        // Compute -2/u.u
        Float nrmus=0.0;
        for(size_t i=k; i!=m; ++i) {
            nrmus+=sqr(u[i]);
        }
        Float mtdnu=(-2)/nrmus;

        // For each column b, compute b-=2u(u.b)/u.u
        for(size_t j=k; j!=n; ++j) {
            Float udtb=0.0;
            for(size_t i=k; i!=m; ++i) {
                udtb+=u[i]*R[i][j];
            }
            Float scl=udtb*mtdnu;
            for(size_t i=k; i!=m; ++i) {
                R[i][j]+=scl*u[i];
            }
        }

        // For the kth column, set R[k][k]=-/+ |a|
        // and R[i][k]=0 for i>k
        for(size_t i=k+1; i!=m; ++i) {
            R[i][k]=0.0;
        }

    } // end of loop on working column k

    // Scale the rows of R to have sum of absolute values equal to 1

    rounding_mode_t prev_rounding_mode=get_rounding_mode();
    for(size_t i=0; i!=m; ++i) {
        set_rounding_mode(upward);
        Float rsum=0.0;
        for(size_t j=i; j!=n; ++j) {
            rsum+=abs(R[i][j]);
        }
        set_rounding_mode(toward_zero);
        for(size_t j=i; j!=n; ++j) {
            R[i][j]/=rsum;
        }
    }
    set_rounding_mode(prev_rounding_mode);

    return make_tuple(R,P);

}

// Returns a matrix T such that the matrix O=AT is a square matrix with
// orthogonal rows, and such that the row absolute value sums of the inverse of
// T are at most 1. Then the image of the unit box under the matrix T is
// an over-approximation of the unit box.
Matrix<Float>
triangular_multiplier(const Matrix<Float>& A)
{
    ARIADNE_ASSERT(A.row_size()<=A.column_size());
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Matrix<Float> R; PivotMatrix P;
    make_ltuple(R,P)=triangular_factor(A);

    Matrix<Float> T(n,m); for(size_t i=0; i!=m; ++i) { T[i][i]=1.0; }

    for(size_t i=0; i!=m; ++i) { assert(R[i][i]!=0.0); }

    for(size_t k=m-1; k!=size_t(-1); --k) {
        Float r=1.0/R[k][k];
        for(size_t i=0; i!=k; ++i) {
            Float s=R[i][k]*r;
            for(size_t j=k; j!=m; ++j) {
                T[i][j]-=s*T[k][j];
            }
        }
        for(size_t j=k; j!=m; ++j) {
            T[k][j]*=r;
        }
    }

    for(size_t k=m-1; k!=size_t(-1); --k) {
        size_t p=P[k];
        for(size_t i=0; i!=m; ++i) {
            Float tmp=T[p][i]; T[p][i]=T[k][i]; T[k][i]=tmp;
        }
    }

    return T;
}



Matrix<Float> pivot_matrix(const array<size_t>& pv)
{
    const size_t n=pv.size();
    array<size_t> perm(n); for(uint i=0; i!=n; ++i) { perm[i]=i; }
    for(size_t i=0; i!=n; ++i) {
        std::swap(perm[i],perm[pv[i]]);
    }
    Matrix<Float> P(n,n);
    for(size_t i=0; i!=n; ++i) {
        P[i][perm[i]]=1;
    }
    return P;
}

// Compute the orthogonal decomposition A=QR without column pivoting. The
// matrix Q is built up as a composition of elementary Householder
// transformations H=I-vv' with |v|=1. Note that inv(H)=H'=H. The vector v is
// chosen to be a multiple of the first working column of A.
tuple< Matrix<Float>, Matrix<Float>, PivotMatrix >
orthogonal_decomposition(const Matrix<Float>& A)
{
    typedef Float X;

    size_t m=A.row_size();
    size_t n=A.column_size();
    Matrix<X> Q(m,m);
    Matrix<X> R(A);
    PivotMatrix P(n);

    array<X> p(n);
    Vector<X> u(m);

    for(size_t i=0; i!=m; ++i) {
        for(size_t j=0; j!=m; ++j) {
            Q[i][j]=0.0;
        }
        Q[i][i]=1.0;
    }

    for(size_t k=0; k!=min(m,n); ++k) {
        //std::cerr<<"k="<<k<<" Q="<<Q<<" R="<<R<<std::flush;

        bool pivoting=true;
        if(pivoting) {
            // Find a pivot column
            size_t pivot_column=k;
            X max_column_norm=0.0;
            for(size_t j=k; j!=n; ++j) {
                X column_norm=0.0;
                for(size_t i=k; i!=m; ++i) {
                    column_norm+=abs(R[i][j]);
                }
                if(column_norm>max_column_norm) {
                    pivot_column=j;
                    max_column_norm=column_norm;
                }
            }
            size_t l=pivot_column;

            // Swap working column and pivot column
            for(size_t i=0; i!=m; ++i) {
                std::swap(R[i][l],R[i][k]);
            }

            // Set pivot column in result
            P[k]=l;
        }

        // Compute |a| where a is the working column
        Float nrmas=0.0;
        for(size_t i=k; i!=m; ++i) {
            nrmas+=R[i][k]*R[i][k];
        }
        Float nrma=sqrt(nrmas);

        // Compute u=a +/- |a|e
        for(size_t i=0; i!=k; ++i) {
            u[i]=0.0;
        }
        for(size_t i=k; i!=m; ++i) {
            u[i]=R[i][k];
        }
        if(u[k]>=0) { u[k]+=nrma; }
        else { u[k]-=nrma; }

        // Compute -2/u.u
        Float nrmus=0.0;
        for(size_t i=k; i!=m; ++i) {
            nrmus+=sqr(u[i]);
        }
        Float mtdnu=(-2)/nrmus;

        // Compute H=(1-2uu'/u'u)
        // Matrix<Float> H(n,n); for(size_t i=0; i!=n; ++i) {
        // H[i][i]=1.0; for(size_t j=0; j!=n; ++j) { H[i][j]+=u[i]*u[j]*mtdnu; } }

        // For each column b of R, compute b-=2u(u.b)/u.u
        for(size_t j=k; j!=n; ++j) {
            Float udtb=0.0;
            for(size_t i=k; i!=m; ++i) {
                udtb+=u[i]*R[i][j];
            }
            Float scl=udtb*mtdnu;
            for(size_t i=k; i!=m; ++i) {
                R[i][j]+=scl*u[i];
            }
        }

        // For the kth column, set R[k][k]=-/+ |a|
        // and R[i][k]=0 for i>k
        for(size_t i=k+1; i!=m; ++i) {
            R[i][k]=0.0;
        }
        if(u[k]>=0) { R[k][k]=-nrma; } else { R[k][k]=nrma; }

        // Update Q'=QH = Q(I-2uu'/u'u)
        // For each row q, compute q-=2u(u.q)/(u.u)
        for(size_t i=0; i!=m; ++i) {
            Float qdtu=0.0;
            for(size_t j=k; j!=m; ++j) {
                qdtu+=Q[i][j]*u[j];
            }
            Float scl=qdtu*mtdnu;
            for(size_t j=k; j!=m; ++j) {
                Q[i][j]+=scl*u[j];
            }
        }

    }

    return make_tuple(Q,R,P);
}

/*
tuple< Matrix<Float>, Matrix<Float> >
orthogonal_decomposition(const Matrix<Float>& A)
{
    typedef Float X;

    size_t m=A.row_size();
    size_t n=A.column_size();
    Matrix<X> O(m,m,0.0);
    Matrix<X> R(A);

    array<X> p(n);

    for(size_t c=0; c!=min(m,n); ++c) {

        // Find a pivot column
        size_t pivot_column=c;
        X max_column_norm=0.0;
        for(size_t j=c; j!=n; ++j) {
            X column_norm=0.0;
            for(size_t i=c; i!=m; ++i) {
                column_norm+=abs(R[i][j]);
            }
            if(column_norm>max_column_norm) {
                pivot_column=j;
                max_column_norm=column_norm;
            }
        }
        size_t j=pivot_column;

        // Swap first column and pivot column
        // FIXME: We do not keep track of column pivoting here.
        for(size_t i=c; i!=m; ++i) {
            std::swap(R[i][j],R[i][c]);
        }

        // Compute inner product of pivot column with remaining columns
        X pivot_norm_square=0.0;
        for(size_t i=c; i!=m; ++i) {
            pivot_norm_square += R[i][c]*R[i][c];
        }

        X inner_product_sum=0.0;
        for(size_t j=c; j!=n; ++j) {
            X inner_product=0.0;
            for(size_t i=c; i!=m; ++i) {
                inner_product += R[i][c]*R[i][j];
            }
            p[j]=inner_product / pivot_norm_square;
            inner_product_sum += inner_product;
        }

        X scale_factor = inner_product_sum / pivot_norm_square;
        for(size_t i=c; i!=m; ++i) {
            O[i][c]=scale_factor * R[i][c];
        }

        for(size_t j=c+1; j!=n; ++j) {
            for(size_t i=c; i!=m; ++i) {
                R[i][j] -= R[c][j] * p[j];
            }
        }
        for(size_t i=c; i!=m; ++i) {
            R[i][c] /= pivot_norm_square;
        }

    }

    return make_tuple(O,R);
}
*/


template class Matrix<Float>;
template Matrix<Float> inverse(const Matrix<Float>&);
template Matrix<Float> solve(const Matrix<Float>&, const Matrix<Float>&);
template Vector<Float> solve(const Matrix<Float>&, const Vector<Float>&);
template class Matrix<Interval>;
template Matrix<Interval> inverse(const Matrix<Interval>&);
template Matrix<Interval> lu_inverse(const Matrix<Interval>&);
template Matrix<Interval> gs_inverse(const Matrix<Interval>&);
template Matrix<Interval> solve(const Matrix<Interval>&, const Matrix<Interval>&);
template Matrix<Interval> lu_solve(const Matrix<Interval>&, const Matrix<Interval>&);
template Matrix<Interval> gs_solve(const Matrix<Interval>&, const Matrix<Interval>&);
template Vector<Interval> solve(const Matrix<Interval>&, const Vector<Interval>&);
#ifdef HAVE_RATIONAL
template class Matrix<Rational>;
template Matrix<Rational> inverse(const Matrix<Rational>&);
template Matrix<Rational> solve(const Matrix<Rational>&, const Matrix<Rational>&);
template Vector<Rational> solve(const Matrix<Rational>&, const Vector<Rational>&);
#endif // HAVE_RATIONAL

} // namespace Ariadne

