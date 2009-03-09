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

template class boost::numeric::ublas::matrix<Ariadne::Float>;
template class boost::numeric::ublas::matrix<Ariadne::Interval>;

#ifdef HAVE_GMPXX_H
template class boost::numeric::ublas::matrix<Ariadne::Rational>;
#endif // HAVE_GMPXX_H


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
    assert(M.row_size()==M.column_size());
    
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


Matrix<Float>
inverse(const Matrix<Float>& A)
{
    return lu_inverse(A);
}

Matrix<Interval>
inverse(const Matrix<Interval>& A)
{
    return lu_inverse(A);
}

#ifdef HAVE_GMPXX_H
Matrix<Rational>
inverse(const Matrix<Rational>& A)
{
    return lu_inverse(A);
}
#endif // HAVE_GMPXX_H

// Returns a matrix R such that A=OR where O is a square matrix with
// orthogonal rows and the sum of the absolute values of the rows of R are
// at most one.
//
// Use Householder transformation H=I-vv' where v=u/|u|
// and u=a +/- |a|e with a and e the working column of A
// and corresponding unit vector.
Matrix<Float>
triangular_factor(const Matrix<Float>& A)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();

    Matrix<Float> R=A;

    // Array of column norm squares
    array<Float> cns(n);

    Vector<Float> u(m);

    for(size_t k=0; k!=min(m,n); ++k) {
        //std::cerr<<"k="<<k<<" R="<<R<<std::endl;

        bool pivoting=false;
        if(pivoting) {
            // Compute column norms
            for(size_t j=k; j!=n; ++j) {
                cns[j]=0.0;
                for(size_t i=k; i!=m; ++i) {
                    cns[j]+=sqr(R[i][j]);
                }
            }
            
            // Find largest column norm
            size_t p=k; Float cnsmax=cns[k];
            for(size_t j=k+1; j!=n; ++j) {
                if(cns[j]>cnsmax) {
                    p=j; cnsmax=cns[j];
                }
            }
            
            // Swap columns p and k
            for(size_t i=0; i!=m; ++i) {
                Float tmp=R[i][k];
                R[i][k]=R[i][p];
                R[i][p]=tmp;
            }
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
    return R;
        
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

    Matrix<Float> R=triangular_factor(A);
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

    return T;
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


// Compute the orthogonal decomposition A=QR without column pivoting. The
// matrix Q is built up as a composition of elementary Householder 
// transformations H=I-vv' with |v|=1. Note that inv(H)=H'=H. The vector v is 
// chosen to be a multiple of the first working column of A. 
tuple< Matrix<Float>, Matrix<Float> > 
orthogonal_decomposition(const Matrix<Float>& A)
{
    typedef Float X;

    size_t m=A.row_size();
    size_t n=A.column_size();
    Matrix<X> Q(m,m);
    Matrix<X> R(A);
    
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

        // No pivoting
/* 
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
        for(size_t i=c; i!=m; ++i) {
            std::swap(R[i][j],R[i][c]);
        }
*/        

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

        // Update Q'=QH = Q(I-2vv')
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

    return make_tuple(Q,R);
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
template class Matrix<Interval>;
#ifdef HAVE_GMPXX_H
template class Matrix<Rational>;
#endif

} // namespace Ariadne

