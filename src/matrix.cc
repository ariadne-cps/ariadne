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
        for(size_t j=0; j!=A.row_size(); ++j) {
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
    
    size_t n=M.row_size();
    Matrix<RealType> A=M;
    Matrix<RealType> B=Matrix<RealType>::identity(n);
    

    for(size_t k=0; k!=n; ++k) {
        // FIXME: Use column pivoting
        assert(A[k][k]!=0);
        RealType r  = 1/A[k][k];
        for(size_t i=k+1; i!=n; ++i) {
            RealType s=A[i][k] * r;
            for(size_t j=0; j<=k; ++j) {
                B[i][j] -= s * B[k][j];
            }
            for(size_t j=k+1; j!=n; ++j) {
                A[i][j] -= s * A[k][j];
            }
            A[i][k] = 0;
        }
        for(size_t j=0; j<=k; ++j) {
            B[k][j] *= r;
        }
        for(size_t j=k+1; j!=n; ++j) {
            A[k][j] *= r;
        }
        A[k][k] = 1;
    }

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

Matrix<Float> triangular_multiplier(const Matrix<Float>& A, size_t b);

// Returns a square matrix R such that the first columns of AR are
// orthogonal, and such that the row absolute value sums are at
// most 1. Then the image of the unit box under the matrix AR is
// an over-approximation of the image of the unit box under A.
Matrix<Float>
triangular_multiplier(const Matrix<Float>& A)
{
    return triangular_multiplier(A,0u);
    typedef Float X;
    size_t m=A.row_size();
    size_t n=A.column_size();
    Matrix<X> B=A;
    Matrix<X> M=Matrix<X>::identity(n);
    for(size_t c=0; c!=std::min(m,n); ++c) {
        Matrix<X> T=triangular_multiplier(B,c);
        Matrix<X> Tinv=inverse(T);
        std::cerr<<"    B="<<B<<"  T["<<c<<"]="<<T<<"  Tinv="<<inverse(T)<<"\n";
        B=prod(B,T);
        M=prod(M,T);
    }
    return M;
}

Matrix<Float>
triangular_multiplier(const Matrix<Float>& A, size_t b)
{
    typedef Float X;

    size_t m=A.row_size();
    size_t n=A.column_size();
    
    Matrix<X> T=Matrix<X>::identity(n);
    
    array<X> p(n);
    
    //for(size_t c=0; c!=std::min(m,n); ++c) {
    for(size_t c=b; c!=b+1u; ++c) {

        for(size_t j=c; j!=n; ++j) {
            X ip=0.0;
            for(size_t i=c; i!=m; ++i) {
                ip+=A[i][c]*A[i][j];
            }
            p[j]=ip;
        }
        
        //T[c][c]=1/p[c];
        for(size_t j=c+1u; j!=n; ++j) {
            T[c][j]=-p[j]/p[c];
            T[c][c]+=abs(T[c][j]);
        }
        
    }

    return T;

}    
     

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

template class Matrix<Float>;
template class Matrix<Interval>;
#ifdef HAVE_GMPXX_H
template class Matrix<Rational>;
#endif

} // namespace Ariadne

