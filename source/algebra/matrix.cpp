/***************************************************************************
 *            algebra/matrix.cpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
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

#include "../numeric/module.hpp"

#include "../config.hpp"

#include "matrix.hpp"
#include "symmetric_matrix.hpp"

#include "../utility/exceptions.hpp"
#include "../numeric/float.hpp"
#include "../numeric/dyadic.hpp"
#include "../numeric/rational.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/covector.hpp"

#include "matrix.tpl.hpp"

namespace Ariadne {



template<class X>
Matrix<X>
lu_inverse(const Matrix<X>& M)
{
    typedef X RealType;
    ARIADNE_ASSERT_MSG(M.row_size()==M.column_size(),"A="<<M);

    SizeType m=M.row_size();
    SizeType n=M.column_size();
    Matrix<RealType> A=M;
    Matrix<RealType> B=Matrix<RealType>::identity(n,M.zero_element());

    // Array of row pivots. The value p[i] gives the row
    // swapped with the ith row in the ith stage.
    Array<SizeType> p(m);
    for(SizeType k=0; k!=m; ++k) { p[k]=k; }


    for(SizeType k=0; k!=std::min(m,n); ++k) {
        // Choose a pivot row
        SizeType iamax=k;
        X amax(0);
        for(SizeType i=k; i!=m; ++i) {
            if(decide(abs(A[i][k])>amax)) {
                iamax=i;
                amax=abs(A[i][k]);
            }
        }

        // Set pivot row
        SizeType q=iamax;
        p[k]=q;

        // Swap rows of both A and B
        for(SizeType j=k; j!=n; ++j) {
            X tmp=A[k][j];
            A[k][j]=A[q][j];
            A[q][j]=tmp;
        }

        for(SizeType j=0; j!=n; ++j) {
            X tmp=B[k][j];
            B[k][j]=B[q][j];
            B[q][j]=tmp;
        }

        RealType r  = rec(A[k][k]);
        for(SizeType i=k+1; i!=n; ++i) {
            RealType s=A[i][k] * r;
            for(SizeType j=0; j!=n; ++j) {
                B[i][j] -= s * B[k][j];
            }
            for(SizeType j=k+1; j!=n; ++j) {
                A[i][j] -= s * A[k][j];
            }
            A[i][k] = 0;
        }
        for(SizeType j=0; j!=n; ++j) {
            B[k][j] *= r;
        }
        for(SizeType j=k+1; j!=n; ++j) {
            A[k][j] *= r;
        }
        A[k][k] = 1;
    }

    // Backsubstitute to find inverse
    for(SizeType k=n; k!=0; ) {
        --k;
        for(SizeType i=0; i!=k; ++i) {
            RealType s=A[i][k];
            for(SizeType j=0; j!=n; ++j) {
                B[i][j] -= s * B[k][j];
            }
            A[i][k] = 0;
        }
    }

    // No need to repivot!

    return B;
}

template<class X> Matrix<Bounds<X>> dd_solve_bounds(const Matrix<Bounds<X>>& A, const Matrix<Bounds<X>>& B);
template<class X> Matrix<Bounds<X>> gs_solve_bounds(const Matrix<Bounds<X>>& A, const Matrix<Bounds<X>>& B);

template<class X> Matrix<X> lu_solve(const Matrix<X>& A, const Matrix<X>& B);
template<class X> Matrix<X> gs_solve(const Matrix<X>& A, const Matrix<X>& B) {
    return gs_solve_bounds(A,B); }
template<class X> Matrix<X> dd_solve(const Matrix<X>& A, const Matrix<X>& B) {
    return gs_solve_bounds(A,B); }

template<class X> Matrix<X> lu_solve(const Matrix<X>& A, const Matrix<X>& B) {
    return lu_inverse(A)*B;
}

template<class X> Vector<X> lu_solve(const Matrix<X>& A, const Vector<X>& b) {
    return lu_inverse(A)*b;
}

// Find a starting solution for a diagonally dominant system
template<class X> Matrix<Bounds<X>> dd_solve_bounds(const Matrix<Bounds<X>>& A, const Matrix<Bounds<X>>& B)
{
    const SizeType n=B.row_size();
    const SizeType m=B.column_size();

    //Compute an upper bound for 1/(|aii|-sum|aij|) using outward rounding
    Vector<UpperBound<X>> c(n,UpperBound<X>(0u));
    for(SizeType i=0; i!=n; ++i) {
        LowerBound<X> rci=mig(A[i][i]);
        for(SizeType j=0; j!=n; ++j) {
            if(j!=i) {
                rci=rci-mag(A[i][j] );
            }
        }
        if(possibly(rci<=0)) {
            ARIADNE_THROW(std::runtime_error,"dd_solve(Matrix<Bounds<X>> A, Matrix<Bounds<X>> B)",
                          "Matrix A="<<A<<" is not diagonally-dominant.");
        }
        c[i]=rec(cast_positive(rci));
    }

    // Compute initial solution
    Matrix<Bounds<X>> R(n,m);
    for(SizeType i=0; i!=n; ++i) {
        Bounds<X> ci(-c[i],+c[i]);
        for(SizeType j=0; j!=m; ++j) {
            R[i][j]=B[i][j]*ci;
        }
    }

    return R;
}


template<class X> Void gs_step(const Matrix<Bounds<X>>& A, const Vector<Bounds<X>>& b, Vector<Bounds<X>>& x)
{
    // Perform Gauss-Seidel iteration
    const SizeType n=x.size();
    for(SizeType i=0; i!=n; ++i) {
        // compute R'[i][j] := (JB[i][j] - Sum{k!=j}JA[i][k]*R[k][j]) / JA[i][i]
        Bounds<X> ri=b[i];
        for(SizeType k=0; k!=n; ++k) {
            if(k!=i) {
                ri-=A[i][k]*x[k];
            }
        }
        if(definitely(A[i][i].lower()>0 || A[i][i].upper()<0)) {
            ri/=A[i][i];
            x[i]=refinement(x[i],ri);
        }
    }
}

template<class X> Matrix<Bounds<X>> gs_solve_bounds(const Matrix<Bounds<X>>& A, const Matrix<Bounds<X>>& B)
{
    ARIADNE_ASSERT(A.row_size()==A.column_size());
    ARIADNE_ASSERT(B.row_size()==A.column_size());
    const SizeType n=B.row_size();
    const SizeType m=B.column_size();

    // Precondition A and B
    Matrix<Approximation<X>> mA(A);

    Matrix<Value<X>> J=cast_exact(inverse(mA));
    Matrix<Bounds<X>> JA=J*A;
    Matrix<Bounds<X>> JB=J*B;

    //Matrix<Bounds<X>> R=dd_solve(JA,JB);
    Matrix<Bounds<X>> R=lu_solve(A,B);
    //std::cerr<<"R="<<R<<"\n";

    // Perform Gauss-Seidel iteration
    SizeType step=0;
    while(step<1) {
        for(SizeType i=0; i!=n; ++i) {
            for(SizeType j=0; j!=m; ++j) {
                // compute R'[i][j] := (JB[i][j] - Sum{k!=j}JA[i][k]*R[k][j]) / JA[i][i]
                Bounds<X> Rij=JB[i][j];
                for(SizeType k=0; k!=n; ++k) {
                    if(k!=i) {
                        Rij-=JA[i][k]*R[k][j];
                    }
                }
                if(definitely(JA[i][i].lower()>0 || JA[i][i].upper()<0)) {
                    Rij/=JA[i][i];
                    R[i][j]=refinement(R[i][j],Rij);
                }
                // FIXME: Use FloatBounds or Float here?
                //R[i][j]=FloatBounds(max(R[i][j].lower(),Rij.lower()),min(R[i][j].upper(),Rij.upper()));
            }
        }
        ++step;
        //std::cerr<<"R="<<R<<"\n";
    }

    return R;
}

template<class X> Vector<X> gs_solve(const Matrix<X>& A, const Vector<X>& b) {
    Matrix<X> B(b.size(),1u); for(SizeType i=0; i!=b.size(); ++i) { B[i][0]=b[i]; }
    Matrix<X> R=gs_solve(A,B);
    Vector<X> r(R.row_size()); for(SizeType i=0; i!=r.size(); ++i) { r[i]=R[i][0]; }
    return r;
}


template<class X> Matrix<X> gs_inverse(const Matrix<X>& A) {
    return gs_solve(A,Matrix<X>::identity(A.row_size()));
}


template<class X> Matrix<ArithmeticType<X>> inverse(const Matrix<X>& A) {
    return lu_inverse(A);
}

template<> Matrix<FloatDPBounds> inverse<>(const Matrix<FloatDPValue>& A) {
    return lu_inverse(Matrix<FloatDPBounds>(A));
}

template<> Matrix<FloatDPBounds> inverse<FloatDPBounds>(const Matrix<FloatDPBounds>& A) {
    try {
        return lu_inverse(A);
    } catch(const DivideByZeroException& e) {
        throw SingularMatrixException();
    }
}


template<class X1, class X2> Matrix<ArithmeticType<X1,X2>> solve(const Matrix<X1>& A, const Matrix<X2>& B) {
    //return inverse(A)*B;
    auto Ainv=inverse(A); return Ainv*B;
}

template<> Matrix<FloatDPBounds> solve(const Matrix<FloatDPBounds>& A, const Matrix<FloatDPBounds>& B) {
    return gs_solve(A,B);
}

template<class X1, class X2> Vector<ArithmeticType<X1,X2>> solve(const Matrix<X1>& A, const Vector<X2>& b) {
    return lu_inverse(A)*b;
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
    const SizeType m=B.column_size();
    const SizeType n=B.row_size();

    Matrix<X> R=B;

    // PivotMatrix rows of B
    for(SizeType k=0; k!=n; ++k) {
        SizeType i=P[k];
        for(SizeType j=0; j!=m; ++j) {
            X tmp=R[i][j];
            R[i][k]=R[k][j];
            R[k][j]=tmp;
        }
    }

    // Backsubstitute on L
    for(SizeType k=0; k!=n; ++k) {
        for(SizeType i=k+1; i!=n; ++i) {
            X s=L[i][k];
            for(SizeType j=0; j!=n; ++j) {
                R[i][j] -= s * R[k][j];
            }
        }
    }

    // Backsubstitute on U with row scaling
    for(SizeType k=n; k!=0; ) {
        --k;
        X s=1/U[k][k];
        for(SizeType j=0; j!=m; ++j) {
            R[k][j] *= s;
        }
        for(SizeType i=0; i!=k; ++i) {
            X t=U[i][k];
            for(SizeType j=0; j!=m; ++j) {
                R[i][j] -= t * R[k][j];
            }
        }
    }

}


// Compute the vector of row norms (sums of absolute values) of A
template<class X>
Vector<RowNormType<X>>
row_norms(const Matrix<X>& A)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    Vector<RowNormType<X>> r(m);

    for(SizeType i=0; i!=m; ++i) {
        r[i]=A.zero_element();
        for(SizeType j=0; j!=n; ++j) {
            r[i]+=abs(A[i][j]);
        }
    }

    return r;

}


// Scale the rows of A to each have sum of absolute values less than or equal to one
// FIXME: This code is inherantly approximate and relies on changing the rounding mode
template<class X>
Void
normalise_rows(Matrix<X>& A)
{

    const SizeType m=A.row_size();
    const SizeType n=A.column_size();

    auto prev_rounding_mode=X::get_rounding_mode();
    X::set_rounding_upward();
    Array<X> row_asums(m);
    for(SizeType i=0; i!=m; ++i) {
        row_asums[i]=0.0;
        for(SizeType j=0; j!=n; ++j) {
            row_asums[i]+=abs(A[i][j]);
        }
    }
    X::set_rounding_toward_zero();
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            A[i][j]/=row_asums[i];
        }
    }
    X::set_rounding_mode(prev_rounding_mode);
}


// Returns a pivot P and matrices L and U such that L is unit lower-triangular,
// U is upper-trianguler and A=PLU.
template<class X>
Tuple< PivotMatrix, Matrix<X>, Matrix<X> >
triangular_decomposition(const Matrix<X>& A)
{
    ARIADNE_ASSERT(A.row_size()==A.column_size());

    SizeType m=A.row_size();
    SizeType n=A.column_size();
    Matrix<X> L=Matrix<X>::identity(m);
    Matrix<X> U=A;

    // Array of row pivots. The value P[i] gives the row
    // swapped with the ith row in the ith stage.
    PivotMatrix P(m);
    for(SizeType k=0; k!=m; ++k) { P[k]=k; }

    for(SizeType k=0; k!=std::min(m,n); ++k) {
        // Choose a pivot row
        SizeType iamax=k;
        X amax=A.zero_element();
        for(SizeType i=k; i!=m; ++i) {
            if(decide(abs(A[i][k])>amax)) {
                iamax=i;
                amax=abs(A[i][k]);
            }
        }

        // Set pivot row
        SizeType l=iamax;
        P[k]=l;

        // Swap rows of L and U
        for(SizeType j=0; j!=k; ++j) {
            X tmp=L[k][j];
            L[k][j]=L[l][j];
            L[l][j]=tmp;
        }

        for(SizeType j=k; j!=n; ++j) {
            X tmp=U[k][j];
            U[k][j]=U[l][j];
            U[l][j]=tmp;
        }

        X r  = 1/U[k][k];
        for(SizeType i=k+1; i!=m; ++i) {
            X s=U[i][k] * r;
            for(SizeType j=k+1; j!=n; ++j) {
                U[i][j] -= s * U[k][j];
            }
            U[i][k] = 0;
            L[i][k] = s;
        }

    }

    return Tuple<PivotMatrix,Matrix<X>,Matrix<X>>{P,L,U};

}

// Returns a matrix R such that A=OR where O is a square matrix with
// orthogonal rows and the sum of the absolute values of the rows of R are
// at most one.
//
// Use Householder transformation H=I-vv' where v=u/|u|
// and u=a +/- |a|e with a and e the working column of A
// and corresponding unit vector.
template<class X> Tuple< Matrix<X>, PivotMatrix>
triangular_factor(const Matrix<X>& A)
{
    X zero = A.zero_element();

    const SizeType m=A.row_size();
    const SizeType n=A.column_size();

    Matrix<X> R=A;
    PivotMatrix P(n); for(Nat i=0; i!=n; ++i) { P[i]=i; }

    // Array of column norm squares
    Array<X> cns(n);

    Vector<X> u(m);

    for(SizeType k=0; k!=std::min(m,n); ++k) {
        //std::cerr<<"k="<<k<<" R="<<R<<std::endl;

        Bool pivoting=true;
        if(pivoting) {
            // Compute column norms
            for(SizeType j=k; j!=n; ++j) {
                cns[j]=zero;
                for(SizeType i=k; i!=m; ++i) {
                    cns[j]+=sqr(R[i][j]);
                }
            }

            // Find largest column norm
            SizeType l=k; X cnsmax=cns[k];
            for(SizeType j=k+1; j!=n; ++j) {
                if(cns[j]>cnsmax) {
                    l=j; cnsmax=cns[j];
                }
            }

            // Swap columns l and k
            for(SizeType i=0; i!=m; ++i) {
                X tmp=R[i][k];
                R[i][k]=R[i][l];
                R[i][l]=tmp;
            }

            // Set pivot element to pivot column
            P[k]=l;
        }

        // Compute |a| where a is the working column
        X nrmas=zero;
        for(SizeType i=k; i!=m; ++i) {
            nrmas+=R[i][k]*R[i][k];
        }
        X nrma=sqrt(nrmas);

        // Compute u=a +/- |a|e
        for(SizeType i=k; i!=m; ++i) {
            u[i]=R[i][k];
        }
        if(u[k]>=0) { u[k]+=nrma; }
        else { u[k]-=nrma; }

        // Compute -2/u.u
        X nrmus=zero;
        for(SizeType i=k; i!=m; ++i) {
            nrmus+=sqr(u[i]);
        }
        X mtdnu=(-2)/nrmus;

        // For each column b, compute b-=2u(u.b)/u.u
        for(SizeType j=k; j!=n; ++j) {
            X udtb=zero;
            for(SizeType i=k; i!=m; ++i) {
                udtb+=u[i]*R[i][j];
            }
            X scl=udtb*mtdnu;
            for(SizeType i=k; i!=m; ++i) {
                R[i][j]+=scl*u[i];
            }
        }

        // For the kth column, set R[k][k]=-/+ |a|
        // and R[i][k]=0 for i>k
        for(SizeType i=k+1; i!=m; ++i) {
            R[i][k]=zero;
        }

    } // end of loop on working column k

    // Scale the rows of R to have sum of absolute values equal to 1
    normalise_rows(R);

    return Tuple<Matrix<X>,PivotMatrix>{R,P};

}

// Returns a matrix T such that the matrix O=AT is a square matrix with
// orthogonal rows, and such that the row absolute value sums of the inverse of
// T are at most 1. Then the image of the unit box under the matrix T is
// an over-approximation of the unit box.
template<class X> Matrix<X>
triangular_multiplier(const Matrix<X>& A)
{
    ARIADNE_ASSERT(A.row_size()<=A.column_size());
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();

    Matrix<X> R; PivotMatrix P;
    make_ltuple(R,P)=triangular_factor(A);

    Matrix<X> T(n,m); for(SizeType i=0; i!=m; ++i) { T[i][i]=1.0; }

    for(SizeType i=0; i!=m; ++i) { assert(R[i][i]!=0.0); }

    for(SizeType k=m-1; k!=SizeType(-1); --k) {
        X r=1/R[k][k];
        for(SizeType i=0; i!=k; ++i) {
            X s=R[i][k]*r;
            for(SizeType j=k; j!=m; ++j) {
                T[i][j]-=s*T[k][j];
            }
        }
        for(SizeType j=k; j!=m; ++j) {
            T[k][j]*=r;
        }
    }

    for(SizeType k=m-1; k!=SizeType(-1); --k) {
        SizeType p=P[k];
        for(SizeType i=0; i!=m; ++i) {
            X tmp=T[p][i]; T[p][i]=T[k][i]; T[k][i]=tmp;
        }
    }

    return T;
}


template<class X> Matrix<X> pivot_matrix(const Array<SizeType>& pv)
{
    const SizeType n=pv.size();
    Array<SizeType> perm(n); for(Nat i=0; i!=n; ++i) { perm[i]=i; }
    for(SizeType i=0; i!=n; ++i) {
        std::swap(perm[i],perm[pv[i]]);
    }
    Matrix<X> P(n,n);
    for(SizeType i=0; i!=n; ++i) {
        P[i][perm[i]]=1;
    }
    return P;
}

template<class X> PivotMatrix::operator Matrix<X> () const {
    return pivot_matrix<X>(this->_ary);
}

OutputStream& operator<<(OutputStream& os, const PivotMatrix& pv) {
    return os << "PivotMatrix(" << static_cast< Matrix<Int> >(pv) << ")";
}

// Compute the orthogonal decomposition A=QR with or without column pivoting. The
// matrix Q is built up as a composition of elementary Householder
// transformations H=I-vv' with |v|=1. Note that inv(H)=H'=H. The vector v is
// chosen to be a multiple of the first working column of A.
template<class X> Tuple< Matrix<X>, Matrix<X>, PivotMatrix >
orthogonal_decomposition(const Matrix<X>& A, Bool allow_pivoting)
{
    X zero=A.zero_element();
    X one=zero+1;

    SizeType m=A.row_size();
    SizeType n=A.column_size();
    Matrix<X> Q(m,m);
    Matrix<X> R(A);
    PivotMatrix P(n);

    Array<X> p(n);
    Vector<X> u(m);

    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=m; ++j) {
            Q[i][j]=zero;
        }
        Q[i][i]=one;
    }

    for(SizeType k=0; k!=std::min(m,n); ++k) {
        //std::cerr<<"k="<<k<<" Q="<<Q<<" R="<<R<<std::flush;

        if(allow_pivoting) {
            // Find a pivot column
            SizeType pivot_column=k;
            X max_column_norm=zero;
            for(SizeType j=k; j!=n; ++j) {
                X column_norm=zero;
                for(SizeType i=k; i!=m; ++i) {
                    column_norm+=abs(R[i][j]);
                }
                if(decide(column_norm>max_column_norm)) {
                    pivot_column=j;
                    max_column_norm=column_norm;
                }
            }
            SizeType l=pivot_column;

            // Swap working column and pivot column
            for(SizeType i=0; i!=m; ++i) {
                std::swap(R[i][l],R[i][k]);
            }

            // Set pivot column in result
            P[k]=l;
        }

        // Compute |a| where a is the working column
        X nrmas=zero;
        for(SizeType i=k; i!=m; ++i) {
            nrmas+=R[i][k]*R[i][k];
        }
        X nrma=sqrt(nrmas);

        // Compute u=a +/- |a|e
        for(SizeType i=0; i!=k; ++i) {
            u[i]=zero;
        }
        for(SizeType i=k; i!=m; ++i) {
            u[i]=R[i][k];
        }
        if(decide(u[k]>=0)) { u[k]+=nrma; }
        else { u[k]-=nrma; }

        // Compute -2/u.u
        X nrmus=zero;
        for(SizeType i=k; i!=m; ++i) {
            nrmus+=sqr(u[i]);
        }
        X mtdnu=(-2)/nrmus;

        // Compute H=(1-2uu'/u'u)
        // Matrix<X> H(n,n); for(SizeType i=0; i!=n; ++i) {
        // H[i][i]=1.0; for(SizeType j=0; j!=n; ++j) { H[i][j]+=u[i]*u[j]*mtdnu; } }

        // For each column b of R, compute b-=2u(u.b)/u.u
        for(SizeType j=k; j!=n; ++j) {
            X udtb=zero;
            for(SizeType i=k; i!=m; ++i) {
                udtb+=u[i]*R[i][j];
            }
            X scl=udtb*mtdnu;
            for(SizeType i=k; i!=m; ++i) {
                R[i][j]+=scl*u[i];
            }
        }

        // For the kth column, set R[k][k]=-/+ |a|
        // and R[i][k]=0 for i>k
        for(SizeType i=k+1; i!=m; ++i) {
            R[i][k]=zero;
        }
        if(decide(u[k]>=0)) { R[k][k]=-nrma; } else { R[k][k]=nrma; }

        // Update Q'=QH = Q(I-2uu'/u'u)
        // For each row q, compute q-=2u(u.q)/(u.u)
        for(SizeType i=0; i!=m; ++i) {
            X qdtu=zero;
            for(SizeType j=k; j!=m; ++j) {
                qdtu+=Q[i][j]*u[j];
            }
            X scl=qdtu*mtdnu;
            for(SizeType j=k; j!=m; ++j) {
                Q[i][j]+=scl*u[j];
            }
        }

    }

    return std::make_tuple(Q,R,P);
}

template<class X>
Tuple< Matrix<X>, Matrix<X> >
orthogonal_decomposition(const Matrix<X>& A)
{
    SizeType m=A.row_size();
    SizeType n=A.column_size();
    X z=A.zero_element();
    Matrix<X> O(m,m,z);
    Matrix<X> R(A);

    Array<X> p(n);

    for(SizeType c=0; c!=std::min(m,n); ++c) {

        // Find a pivot column
        SizeType pivot_column=c;
        X max_column_norm=z;
        for(SizeType j=c; j!=n; ++j) {
            X column_norm=z;
            for(SizeType i=c; i!=m; ++i) {
                column_norm+=abs(R[i][j]);
            }
            if(decide(column_norm>max_column_norm)) {
                pivot_column=j;
                max_column_norm=column_norm;
            }
        }
        SizeType j=pivot_column;

        // Swap first column and pivot column
        // FIXME: We do not keep track of column pivoting here.
        for(SizeType i=c; i!=m; ++i) {
            std::swap(R[i][j],R[i][c]);
        }

        // Compute inner product of pivot column with remaining columns
        X pivot_norm_square=z;
        for(SizeType i=c; i!=m; ++i) {
            pivot_norm_square += R[i][c]*R[i][c];
        }

        X inner_product_sum=z;
        for(SizeType k=c; k!=n; ++k) {
            X inner_product=z;
            for(SizeType i=c; i!=m; ++i) {
                inner_product += R[i][c]*R[i][k];
            }
            p[k]=inner_product / pivot_norm_square;
            inner_product_sum += inner_product;
        }

        X scale_factor = inner_product_sum / pivot_norm_square;
        for(SizeType i=c; i!=m; ++i) {
            O[i][c]=scale_factor * R[i][c];
        }

        for(SizeType k=c+1; k!=n; ++k) {
            for(SizeType i=c; i!=m; ++i) {
                R[i][k] -= R[c][k] * p[k];
            }
        }
        for(SizeType i=c; i!=m; ++i) {
            R[i][c] /= pivot_norm_square;
        }

    }

    return make_tuple(O,R);
}


template<class X> Matrix<MidpointType<X>> midpoint(Matrix<X> const& A) {
    Matrix<MidpointType<X>> R(A.row_size(),A.column_size(),midpoint(A.zero_element()));
    for(SizeType i=0; i!=A.row_size(); ++i) {
        for(SizeType j=0; j!=A.column_size(); ++j) {
            R.at(i,j)=midpoint(A.at(i,j));
        }
    }
    return R;
}

template<class X> Matrix<SingletonType<X>> cast_singleton(Matrix<X> const& A) {
    Matrix<SingletonType<X>> R(A.row_size(),A.column_size(),cast_singleton(A.zero_element()));
    for(SizeType i=0; i!=A.row_size(); ++i) {
        for(SizeType j=0; j!=A.column_size(); ++j) {
            R.at(i,j)=cast_singleton(A.at(i,j));
        }
    }
    return R;
}


template<class X> Matrix<Value<X>> cast_exact(Matrix<Approximation<X>> const& A) {
    return reinterpret_cast<Matrix<Value<X>> const&>(A);
}

template<class AX> Matrix<decltype(cast_exact(declval<AX>()))> cast_exact(Matrix<AX> const& A) {
    typedef decltype(cast_exact(declval<AX>())) EX;
    return reinterpret_cast<Matrix<EX> const&>(A);
}


template class Matrix<FloatDP>;
template Matrix<FloatDP> inverse(const Matrix<FloatDP>&);
template Vector<FloatDP> solve(const Matrix<FloatDP>&, const Vector<FloatDP>&);
template Void normalise_rows(Matrix<FloatDP>&);

template class Matrix<FloatDPApproximation>;
template Matrix<FloatDPApproximation> inverse(const Matrix<FloatDPApproximation>&);
template Matrix<FloatDPApproximation> solve(const Matrix<FloatDPApproximation>&, const Matrix<FloatDPApproximation>&);
template Vector<FloatDPApproximation> solve(const Matrix<FloatDPApproximation>&, const Vector<FloatDPApproximation>&);
//template Matrix<FloatDPApproximation> lu_solve(const Matrix<FloatDPApproximation>&, const Matrix<FloatDPApproximation>&);
//template Matrix<FloatDPApproximation> lu_inverse(const Matrix<FloatDPApproximation>&, const Matrix<FloatDPApproximation>&);
template Tuple<PivotMatrix,Matrix<FloatDPApproximation>,Matrix<FloatDPApproximation>> triangular_decomposition(Matrix<FloatDPApproximation> const&);
template Tuple<Matrix<FloatDPApproximation>,Matrix<FloatDPApproximation>,PivotMatrix> orthogonal_decomposition(Matrix<FloatDPApproximation> const&, Bool);
template Tuple<Matrix<FloatDPApproximation>,Matrix<FloatDPApproximation>> orthogonal_decomposition(Matrix<FloatDPApproximation> const&);

template Vector<FloatDPApproximation> row_norms(Matrix<FloatDPApproximation> const&);

template class Matrix<FloatDPBounds>;
template Matrix<FloatDPBounds> lu_inverse(const Matrix<FloatDPBounds>&);
template Matrix<FloatDPBounds> gs_inverse(const Matrix<FloatDPBounds>&);
template Matrix<FloatDPBounds> lu_solve(const Matrix<FloatDPBounds>&, const Matrix<FloatDPBounds>&);
template Vector<FloatDPBounds> solve(const Matrix<FloatDPBounds>&, const Vector<FloatDPBounds>&);
template Vector<FloatDPBounds> lu_solve(const Matrix<FloatDPBounds>&, const Vector<FloatDPBounds>&);
template Vector<FloatDPBounds> gs_solve(const Matrix<FloatDPBounds>&, const Vector<FloatDPBounds>&);
template Matrix<FloatDPBounds> gs_solve(const Matrix<FloatDPBounds>&, const Matrix<FloatDPBounds>&);
template Matrix<MidpointType<FloatDPBounds>> midpoint(Matrix<FloatDPBounds> const&);
template Tuple<PivotMatrix,Matrix<FloatDPBounds>,Matrix<FloatDPBounds>> triangular_decomposition(Matrix<FloatDPBounds> const&);
template Tuple<Matrix<FloatDPBounds>,Matrix<FloatDPBounds>> orthogonal_decomposition(Matrix<FloatDPBounds> const&);

template class Matrix<FloatDPValue>;

template class Matrix<FloatMPApproximation>;
template class Matrix<FloatMPBounds>;
template class Matrix<FloatMPValue>;

template Matrix<FloatMPApproximation> inverse(const Matrix<FloatMPApproximation>&);
template Vector<FloatMPApproximation> solve(const Matrix<FloatMPApproximation>&, const Vector<FloatMPApproximation>&);
template Matrix<FloatMPApproximation> solve(const Matrix<FloatMPApproximation>&, const Matrix<FloatMPApproximation>&);
template Tuple<PivotMatrix,Matrix<FloatMPBounds>,Matrix<FloatMPBounds>> triangular_decomposition(Matrix<FloatMPBounds> const&);
template Tuple<Matrix<FloatMPApproximation>,Matrix<FloatMPApproximation>> orthogonal_decomposition(Matrix<FloatMPApproximation> const&);
template Vector<FloatMPApproximation> row_norms(Matrix<FloatMPApproximation> const&);

template Matrix<FloatMPBounds> inverse(const Matrix<FloatMPBounds>&);
template Vector<FloatMPBounds> solve(const Matrix<FloatMPBounds>&, const Vector<FloatMPBounds>&);
template Matrix<FloatMPBounds> solve(const Matrix<FloatMPBounds>&, const Matrix<FloatMPBounds>&);
template Matrix<FloatMPBounds> lu_inverse(const Matrix<FloatMPBounds>&);
template Matrix<FloatMPBounds> gs_inverse(const Matrix<FloatMPBounds>&);
template Matrix<FloatMPBounds> lu_solve(const Matrix<FloatMPBounds>&, const Matrix<FloatMPBounds>&);
template Vector<FloatMPBounds> gs_solve(const Matrix<FloatMPBounds>&, const Vector<FloatMPBounds>&);
template Matrix<FloatMPBounds> gs_solve(const Matrix<FloatMPBounds>&, const Matrix<FloatMPBounds>&);
template Tuple<PivotMatrix,Matrix<FloatMPApproximation>,Matrix<FloatMPApproximation>> triangular_decomposition(Matrix<FloatMPApproximation> const&);
template Tuple<Matrix<FloatMPBounds>,Matrix<FloatMPBounds>> orthogonal_decomposition(Matrix<FloatMPBounds> const&);

template class Matrix<Real>;

template PositiveFloatDPUpperBound sup_norm(const Matrix<FloatDPBounds>& A);
template FloatDPUpperBound log_norm(const Matrix<FloatDPBounds>& A);

template class Matrix<Dyadic>;

template class Matrix<Rational>;
template Matrix<Rational> inverse(const Matrix<Rational>&);
template Matrix<Rational> solve(const Matrix<Rational>&, const Matrix<Rational>&);
template Vector<Rational> solve(const Matrix<Rational>&, const Vector<Rational>&);
Rational midpoint(Rational);
template<> Matrix<Rational> midpoint(Matrix<Rational> const& A) { return A; }

template class SymmetricMatrix<FloatDPApproximation>;

} // namespace Ariadne

#include "../geometry/interval.hpp"

namespace Ariadne {
template class Matrix<FloatDPUpperInterval>;
template Matrix<SingletonType<FloatDPUpperInterval>> cast_singleton(Matrix<FloatDPUpperInterval> const&);
template Matrix<MidpointType<FloatDPUpperInterval>> midpoint(Matrix<FloatDPUpperInterval> const&);
template Matrix<FloatDPUpperInterval> inverse(const Matrix<FloatDPUpperInterval>&);
template Vector<FloatDPUpperInterval> solve(const Matrix<FloatDPUpperInterval>&, const Vector<FloatDPUpperInterval>&);

template class Matrix<FloatMPUpperInterval>;
} // namespace Ariadne

