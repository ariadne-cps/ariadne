/***************************************************************************
 *            matrix.cc
 *
 *  Copyright 2008-14  Alberto Casagrande, Pieter Collins
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

#include "numeric/module.h"

#include "config.h"

#include "algebra/matrix.h"

#include "numeric/float.h"
#include "numeric/rational.h"
#include "algebra/vector.h"
#include "algebra/covector.h"

#include "matrix.tcc"

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
    Matrix<RealType> B=Matrix<RealType>::identity(n);

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
        SizeType i=iamax;
        p[k]=i;

        // Swap rows of both A and B
        for(SizeType j=k; j!=n; ++j) {
            X tmp=A[k][j];
            A[k][j]=A[i][j];
            A[i][j]=tmp;
        }

        for(SizeType j=0; j!=n; ++j) {
            X tmp=B[k][j];
            B[k][j]=B[i][j];
            B[i][j]=tmp;
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

template<class X> Matrix<X> lu_solve(const Matrix<X>& A, const Matrix<X>& B);
template<class X> Matrix<X> gs_solve(const Matrix<X>& A, const Matrix<X>& B);
template<class X> Matrix<X> dd_solve(const Matrix<X>& A, const Matrix<X>& B);

template<class X> Matrix<X> lu_solve(const Matrix<X>& A, const Matrix<X>& B) {
    return lu_inverse(A)*B;
}

template<class X> Vector<X> lu_solve(const Matrix<X>& A, const Vector<X>& b) {
    return lu_inverse(A)*b;

}

// Find a starting solution for a diagonally dominant system
template<> Matrix<BoundFloat> dd_solve(const Matrix<BoundFloat>& A, const Matrix<BoundFloat>& B)
{
    const SizeType n=B.row_size();
    const SizeType m=B.column_size();

    //Compute an upper bound for 1/(|aii|-sum|aij|) using outward rounding
    Vector<UpperFloat> c(n,0.0_x);
    for(SizeType i=0; i!=n; ++i) {
        LowerFloat rci=mig(A[i][i]);
        for(SizeType j=0; j!=n; ++j) {
            if(j!=i) {
                rci=rci-mag(A[i][j] );
            }
        }
        if(rci<=0) {
            ARIADNE_THROW(std::runtime_error,"dd_solve(Matrix<BoundFloat> A, Matrix<BoundFloat> B)",
                          "Matrix A="<<A<<" is not diagonally-dominant.");
        }
        c[i]=rec(rci);
    }

    // Compute initial solution
    Matrix<BoundFloat> R(n,m);
    for(SizeType i=0; i!=n; ++i) {
        BoundFloat ci(-c[i],+c[i]);
        for(SizeType j=0; j!=m; ++j) {
            R[i][j]=B[i][j]*ci;
        }
    }

    return R;
}


template<> Void gs_step(const Matrix<BoundFloat>& A, const Vector<BoundFloat>& b, Vector<BoundFloat>& x)
{
    // Perform Gauss-Seidel iteration
    const SizeType n=x.size();
    for(SizeType i=0; i!=n; ++i) {
        // compute R'[i][j] := (JB[i][j] - Sum{k!=j}JA[i][k]*R[k][j]) / JA[i][i]
        BoundFloat ri=b[i];
        for(SizeType k=0; k!=n; ++k) {
            if(k!=i) {
                ri-=A[i][k]*x[k];
            }
        }
        if(A[i][i].lower()>0 || A[i][i].upper()<0) {
            ri/=A[i][i];
            x[i]=refinement(x[i],ri);
        }
    }
}

template<> Matrix<BoundFloat> gs_solve(const Matrix<BoundFloat>& A, const Matrix<BoundFloat>& B)
{
    ARIADNE_ASSERT(A.row_size()==A.column_size());
    ARIADNE_ASSERT(B.row_size()==A.column_size());
    const SizeType n=B.row_size();
    const SizeType m=B.column_size();

    // Precondition A and B
    Matrix<ApproximateFloat> mA(A);

    Matrix<ExactFloat> J=make_exact(inverse(mA));
    Matrix<BoundFloat> JA=J*A;
    Matrix<BoundFloat> JB=J*B;

    //Matrix<BoundFloat> R=dd_solve(JA,JB);
    Matrix<BoundFloat> R=lu_solve(A,B);
    //std::cerr<<"R="<<R<<"\n";

    // Perform Gauss-Seidel iteration
    SizeType step=0;
    while(step<1) {
        for(SizeType i=0; i!=n; ++i) {
            for(SizeType j=0; j!=m; ++j) {
                // compute R'[i][j] := (JB[i][j] - Sum{k!=j}JA[i][k]*R[k][j]) / JA[i][i]
                BoundFloat Rij=JB[i][j];
                for(SizeType k=0; k!=n; ++k) {
                    if(k!=i) {
                        Rij-=JA[i][k]*R[k][j];
                    }
                }
                if(JA[i][i].lower()>0 || JA[i][i].upper()<0) {
                    Rij/=JA[i][i];
                    R[i][j]=refinement(R[i][j],Rij);
                }
                // FIXME: Use FloatInterval or Flt here?
                //R[i][j]=FloatInterval(max(R[i][j].lower(),Rij.lower()),min(R[i][j].upper(),Rij.upper()));
            }
        }
        ++step;
        //std::cerr<<"R="<<R<<"\n";
    }

    return R;
}

template<class X> Vector<X> gs_solve(const Matrix<X>& A, const Vector<X>& b) {
    Matrix<X> B(b.size(),1u); for(SizeType i=0; i!=b.size(); ++i) { B[i][1]=b[i]; }
    Matrix<X> R=gs_solve(A,B);
    Vector<X> r(R.row_size()); for(SizeType i=0; i!=r.size(); ++i) { r[i]=R[i][1]; }
    return r;
}


template<class X> Matrix<X> gs_inverse(const Matrix<X>& A) {
    return gs_solve(A,Matrix<X>::identity(A.row_size()));
}


Matrix<ValidatedFloat> inverse(const Matrix<ExactFloat>& A) {
    return lu_inverse(Matrix<ValidatedFloat>(A));
}

template<class X> Matrix<X> inverse(const Matrix<X>& A) {
    return lu_inverse(A);
}

struct DivideByZeroException { };
struct SingularMatrixException { };

template<> Matrix<BoundFloat> inverse<BoundFloat>(const Matrix<BoundFloat>& A) {
    try {
        return lu_inverse(A);
    } catch(const DivideByZeroException& e) {
        throw SingularMatrixException();
    }
}


template<class X> Matrix<X> solve(const Matrix<X>& A, const Matrix<X>& B) {
    //return inverse(A)*B;
    Matrix<X> Ainv=inverse(A); return Ainv*B;
}

template<> Matrix<BoundFloat> solve(const Matrix<BoundFloat>& A, const Matrix<BoundFloat>& B) {
    return gs_solve(A,B);
}

template<class X> Vector<X> solve(const Matrix<X>& A, const Vector<X>& b) {
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
            X s=U[i][k];
            for(SizeType j=0; j!=m; ++j) {
                R[i][j] -= s * R[k][j];
            }
        }
    }

}


// Compute the vector of row norms (sums of absolute values) of A
Vector<Flt>
row_norms(const Matrix<Flt>& A)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();
    Vector<Flt> r(m);

    for(SizeType i=0; i!=m; ++i) {
        r[i]=0.0;
        for(SizeType j=0; j!=n; ++j) {
            r[i]+=abs(A[i][j]);
        }
    }

    return r;

}


// Scale the rows of A to each have some of absolute values equal to one
Matrix<Flt>
normalise_rows(const Matrix<Flt>& A)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();

    Matrix<Flt> R=A;

    Array<Flt> row_asums(m);
    auto prev_rounding_mode=get_rounding_mode();
    set_rounding_upward();
    for(SizeType i=0; i!=m; ++i) {
        row_asums[i]=0.0;
        for(SizeType j=0; j!=n; ++j) {
            row_asums[i]+=abs(A[i][j]);
        }
    }
    set_rounding_toward_zero();
    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=n; ++j) {
            R[i][j]/=row_asums[i];
        }
    }
    set_rounding_mode(prev_rounding_mode);

    return R;
}


// Returns a pivot P and matrices L and U such that L is unit lower-triangular,
// U is upper-trianguler and A=PLU.
Tuple< PivotMatrix, Matrix<Flt>, Matrix<Flt> >
triangular_decomposition(const Matrix<Flt>& A)
{
    typedef Flt RealType;
    ARIADNE_ASSERT(A.row_size()==A.column_size());

    SizeType m=A.row_size();
    SizeType n=A.column_size();
    Matrix<RealType> L=Matrix<RealType>::identity(m);
    Matrix<RealType> U=A;

    // Array of row pivots. The value P[i] gives the row
    // swapped with the ith row in the ith stage.
    PivotMatrix P(m);
    for(SizeType k=0; k!=m; ++k) { P[k]=k; }

    for(SizeType k=0; k!=std::min(m,n); ++k) {
        // Choose a pivot row
        SizeType iamax=k;
        RealType amax=0;
        for(SizeType i=k; i!=m; ++i) {
            if(abs(A[i][k])>amax) {
                iamax=i;
                amax=abs(A[i][k]);
            }
        }

        // Set pivot row
        SizeType l=iamax;
        P[k]=l;

        // Swap rows of L and U
        for(SizeType j=0; j!=k; ++j) {
            RealType tmp=L[k][j];
            L[k][j]=L[l][j];
            L[l][j]=tmp;
        }

        for(SizeType j=k; j!=n; ++j) {
            RealType tmp=U[k][j];
            U[k][j]=U[l][j];
            U[l][j]=tmp;
        }

        RealType r  = 1/U[k][k];
        for(SizeType i=k+1; i!=m; ++i) {
            RealType s=U[i][k] * r;
            for(SizeType j=k+1; j!=n; ++j) {
                U[i][j] -= s * U[k][j];
            }
            U[i][k] = 0;
            L[i][k] = s;
        }

    }

    return Tuple<PivotMatrix,Matrix<Flt>,Matrix<Flt>>{P,L,U};

}

using std::min;
// Returns a matrix R such that A=OR where O is a square matrix with
// orthogonal rows and the sum of the absolute values of the rows of R are
// at most one.
//
// Use Householder transformation H=I-vv' where v=u/|u|
// and u=a +/- |a|e with a and e the working column of A
// and corresponding unit vector.
Tuple< Matrix<Flt>, PivotMatrix>
triangular_factor(const Matrix<Flt>& A)
{
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();

    Matrix<Flt> R=A;
    PivotMatrix P(n); for(uint i=0; i!=n; ++i) { P[i]=i; }

    // Array of column norm squares
    Array<Flt> cns(n);

    Vector<Flt> u(m);

    for(SizeType k=0; k!=min(m,n); ++k) {
        //std::cerr<<"k="<<k<<" R="<<R<<std::endl;

        bool pivoting=true;
        if(pivoting) {
            // Compute column norms
            for(SizeType j=k; j!=n; ++j) {
                cns[j]=0.0;
                for(SizeType i=k; i!=m; ++i) {
                    cns[j]+=sqr(R[i][j]);
                }
            }

            // Find largest column norm
            SizeType l=k; Flt cnsmax=cns[k];
            for(SizeType j=k+1; j!=n; ++j) {
                if(cns[j]>cnsmax) {
                    l=j; cnsmax=cns[j];
                }
            }

            // Swap columns l and k
            for(SizeType i=0; i!=m; ++i) {
                Flt tmp=R[i][k];
                R[i][k]=R[i][l];
                R[i][l]=tmp;
            }

            // Set pivot element to pivot column
            P[k]=l;
        }

        // Compute |a| where a is the working column
        Flt nrmas=0.0;
        for(SizeType i=k; i!=m; ++i) {
            nrmas+=R[i][k]*R[i][k];
        }
        Flt nrma=sqrt(nrmas);

        // Compute u=a +/- |a|e
        for(SizeType i=k; i!=m; ++i) {
            u[i]=R[i][k];
        }
        if(u[k]>=0) { u[k]+=nrma; }
        else { u[k]-=nrma; }

        // Compute -2/u.u
        Flt nrmus=0.0;
        for(SizeType i=k; i!=m; ++i) {
            nrmus+=sqr(u[i]);
        }
        Flt mtdnu=(-2)/nrmus;

        // For each column b, compute b-=2u(u.b)/u.u
        for(SizeType j=k; j!=n; ++j) {
            Flt udtb=0.0;
            for(SizeType i=k; i!=m; ++i) {
                udtb+=u[i]*R[i][j];
            }
            Flt scl=udtb*mtdnu;
            for(SizeType i=k; i!=m; ++i) {
                R[i][j]+=scl*u[i];
            }
        }

        // For the kth column, set R[k][k]=-/+ |a|
        // and R[i][k]=0 for i>k
        for(SizeType i=k+1; i!=m; ++i) {
            R[i][k]=0.0;
        }

    } // end of loop on working column k

    // Scale the rows of R to have sum of absolute values equal to 1

    RoundingModeType prev_rounding_mode=get_rounding_mode();
    for(SizeType i=0; i!=m; ++i) {
        set_rounding_upward();
        Flt rsum=0.0;
        for(SizeType j=i; j!=n; ++j) {
            rsum+=abs(R[i][j]);
        }
        set_rounding_toward_zero();
        for(SizeType j=i; j!=n; ++j) {
            R[i][j]/=rsum;
        }
    }
    set_rounding_mode(prev_rounding_mode);

    return Tuple<Matrix<Flt>,PivotMatrix>{R,P};

}

// Returns a matrix T such that the matrix O=AT is a square matrix with
// orthogonal rows, and such that the row absolute value sums of the inverse of
// T are at most 1. Then the image of the unit box under the matrix T is
// an over-approximation of the unit box.
Matrix<Flt>
triangular_multiplier(const Matrix<Flt>& A)
{
    ARIADNE_ASSERT(A.row_size()<=A.column_size());
    const SizeType m=A.row_size();
    const SizeType n=A.column_size();

    Matrix<Flt> R; PivotMatrix P;
    std::make_tuple(R,P)=triangular_factor(A);

    Matrix<Flt> T(n,m); for(SizeType i=0; i!=m; ++i) { T[i][i]=1.0; }

    for(SizeType i=0; i!=m; ++i) { assert(R[i][i]!=0.0); }

    for(SizeType k=m-1; k!=SizeType(-1); --k) {
        Flt r=1.0/R[k][k];
        for(SizeType i=0; i!=k; ++i) {
            Flt s=R[i][k]*r;
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
            Flt tmp=T[p][i]; T[p][i]=T[k][i]; T[k][i]=tmp;
        }
    }

    return T;
}


template<class X> Matrix<X> pivot_matrix(const Array<SizeType>& pv)
{
    const SizeType n=pv.size();
    Array<SizeType> perm(n); for(uint i=0; i!=n; ++i) { perm[i]=i; }
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
    return os << static_cast< Matrix<ExactFloat> >(pv);
}

// Compute the orthogonal decomposition A=QR with or without column pivoting. The
// matrix Q is built up as a composition of elementary Householder
// transformations H=I-vv' with |v|=1. Note that inv(H)=H'=H. The vector v is
// chosen to be a multiple of the first working column of A.
Tuple< Matrix<Flt>, Matrix<Flt>, PivotMatrix >
orthogonal_decomposition(const Matrix<Flt>& A, bool allow_pivoting)
{
    typedef Flt X;

    SizeType m=A.row_size();
    SizeType n=A.column_size();
    Matrix<X> Q(m,m);
    Matrix<X> R(A);
    PivotMatrix P(n);

    Array<X> p(n);
    Vector<X> u(m);

    for(SizeType i=0; i!=m; ++i) {
        for(SizeType j=0; j!=m; ++j) {
            Q[i][j]=0.0;
        }
        Q[i][i]=1.0;
    }

    for(SizeType k=0; k!=min(m,n); ++k) {
        //std::cerr<<"k="<<k<<" Q="<<Q<<" R="<<R<<std::flush;

        if(allow_pivoting) {
            // Find a pivot column
            SizeType pivot_column=k;
            X max_column_norm=0.0;
            for(SizeType j=k; j!=n; ++j) {
                X column_norm=0.0;
                for(SizeType i=k; i!=m; ++i) {
                    column_norm+=abs(R[i][j]);
                }
                if(column_norm>max_column_norm) {
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
        Flt nrmas=0.0;
        for(SizeType i=k; i!=m; ++i) {
            nrmas+=R[i][k]*R[i][k];
        }
        Flt nrma=sqrt(nrmas);

        // Compute u=a +/- |a|e
        for(SizeType i=0; i!=k; ++i) {
            u[i]=0.0;
        }
        for(SizeType i=k; i!=m; ++i) {
            u[i]=R[i][k];
        }
        if(u[k]>=0) { u[k]+=nrma; }
        else { u[k]-=nrma; }

        // Compute -2/u.u
        Flt nrmus=0.0;
        for(SizeType i=k; i!=m; ++i) {
            nrmus+=sqr(u[i]);
        }
        Flt mtdnu=(-2)/nrmus;

        // Compute H=(1-2uu'/u'u)
        // Matrix<Flt> H(n,n); for(SizeType i=0; i!=n; ++i) {
        // H[i][i]=1.0; for(SizeType j=0; j!=n; ++j) { H[i][j]+=u[i]*u[j]*mtdnu; } }

        // For each column b of R, compute b-=2u(u.b)/u.u
        for(SizeType j=k; j!=n; ++j) {
            Flt udtb=0.0;
            for(SizeType i=k; i!=m; ++i) {
                udtb+=u[i]*R[i][j];
            }
            Flt scl=udtb*mtdnu;
            for(SizeType i=k; i!=m; ++i) {
                R[i][j]+=scl*u[i];
            }
        }

        // For the kth column, set R[k][k]=-/+ |a|
        // and R[i][k]=0 for i>k
        for(SizeType i=k+1; i!=m; ++i) {
            R[i][k]=0.0;
        }
        if(u[k]>=0) { R[k][k]=-nrma; } else { R[k][k]=nrma; }

        // Update Q'=QH = Q(I-2uu'/u'u)
        // For each row q, compute q-=2u(u.q)/(u.u)
        for(SizeType i=0; i!=m; ++i) {
            Flt qdtu=0.0;
            for(SizeType j=k; j!=m; ++j) {
                qdtu+=Q[i][j]*u[j];
            }
            Flt scl=qdtu*mtdnu;
            for(SizeType j=k; j!=m; ++j) {
                Q[i][j]+=scl*u[j];
            }
        }

    }

    return std::make_tuple(Q,R,P);
}

/*
Tuple< Matrix<Flt>, Matrix<Flt> >
orthogonal_decomposition(const Matrix<Flt>& A)
{
    typedef Flt X;

    SizeType m=A.row_size();
    SizeType n=A.column_size();
    Matrix<X> O(m,m,0.0);
    Matrix<X> R(A);

    Array<X> p(n);

    for(SizeType c=0; c!=min(m,n); ++c) {

        // Find a pivot column
        SizeType pivot_column=c;
        X max_column_norm=0.0;
        for(SizeType j=c; j!=n; ++j) {
            X column_norm=0.0;
            for(SizeType i=c; i!=m; ++i) {
                column_norm+=abs(R[i][j]);
            }
            if(column_norm>max_column_norm) {
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
        X pivot_norm_square=0.0;
        for(SizeType i=c; i!=m; ++i) {
            pivot_norm_square += R[i][c]*R[i][c];
        }

        X inner_product_sum=0.0;
        for(SizeType j=c; j!=n; ++j) {
            X inner_product=0.0;
            for(SizeType i=c; i!=m; ++i) {
                inner_product += R[i][c]*R[i][j];
            }
            p[j]=inner_product / pivot_norm_square;
            inner_product_sum += inner_product;
        }

        X scale_factor = inner_product_sum / pivot_norm_square;
        for(SizeType i=c; i!=m; ++i) {
            O[i][c]=scale_factor * R[i][c];
        }

        for(SizeType j=c+1; j!=n; ++j) {
            for(SizeType i=c; i!=m; ++i) {
                R[i][j] -= R[c][j] * p[j];
            }
        }
        for(SizeType i=c; i!=m; ++i) {
            R[i][c] /= pivot_norm_square;
        }

    }

    return make_tuple(O,R);
}
*/

template<class AX> Matrix<decltype(make_exact(declval<AX>()))> make_exact(Matrix<AX> const& A) {
    typedef decltype(make_exact(declval<AX>())) EX;
    return reinterpret_cast<Matrix<EX> const&>(A);
}


template class Matrix<Flt>;

template class Matrix<ApproximateFloat>;
template Matrix<ApproximateFloat> inverse(const Matrix<ApproximateFloat>&);
template Matrix<ApproximateFloat> solve(const Matrix<ApproximateFloat>&, const Matrix<ApproximateFloat>&);
template Vector<ApproximateFloat> solve(const Matrix<ApproximateFloat>&, const Vector<ApproximateFloat>&);
//template Matrix<ApproximateFloat> lu_solve(const Matrix<ApproximateFloat>&, const Matrix<ApproximateFloat>&);
//template Matrix<ApproximateFloat> lu_inverse(const Matrix<ApproximateFloat>&, const Matrix<ApproximateFloat>&);
template class Matrix<MetricFloat>;
template Vector<MetricFloat> solve(const Matrix<MetricFloat>&, const Vector<MetricFloat>&);

template class Matrix<BoundFloat>;
template Matrix<BoundFloat> inverse(const Matrix<BoundFloat>&);
template Matrix<BoundFloat> lu_inverse(const Matrix<BoundFloat>&);
template Matrix<BoundFloat> gs_inverse(const Matrix<BoundFloat>&);
template Matrix<BoundFloat> solve(const Matrix<BoundFloat>&, const Matrix<BoundFloat>&);
template Matrix<BoundFloat> lu_solve(const Matrix<BoundFloat>&, const Matrix<BoundFloat>&);
template Matrix<BoundFloat> gs_solve(const Matrix<BoundFloat>&, const Matrix<BoundFloat>&);
template Vector<BoundFloat> solve(const Matrix<BoundFloat>&, const Vector<BoundFloat>&);
template Vector<BoundFloat> lu_solve(const Matrix<BoundFloat>&, const Vector<BoundFloat>&);
template Vector<BoundFloat> gs_solve(const Matrix<BoundFloat>&, const Vector<BoundFloat>&);

template class Matrix<ExactFloat>;

template class Matrix<Real>;

template ErrorFloat sup_norm(const Matrix<BoundFloat>& A);
template UpperFloat log_norm(const Matrix<BoundFloat>& A);

template Matrix<ExactFloat64> make_exact(const Matrix<ApprxFloat64>& mx);

#ifdef HAVE_GMPXX_H
template class Matrix<Rational>;
template Matrix<Rational> inverse(const Matrix<Rational>&);
template Matrix<Rational> solve(const Matrix<Rational>&, const Matrix<Rational>&);
template Vector<Rational> solve(const Matrix<Rational>&, const Vector<Rational>&);
#endif

} // namespace Ariadne

