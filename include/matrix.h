/***************************************************************************
 *            matrix.h
 *
 *  Copyright 2005-8  Alberto Casagrande, Pieter Collins
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
 
/*! \file matrix.h
 *  \brief Matrices.
 */

#ifndef ARIADNE_MATRIX_H
#define ARIADNE_MATRIX_H 

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <cstdarg>

#include "macros.h"
#include "exceptions.h"
#include "numeric.h"
#include "tuple.h"

using namespace boost::numeric;

namespace Ariadne {

/// A matrix over a field. See also \link Ariadne::Vector \c Vector<X> \endlink.
template<class X>
class Matrix
    : public ublas::matrix<X>
{
  public:
    //@{
    //! \name Constructors

    //! Default constructor makes a \f$0\times0\f$ matrix.
    Matrix()
        : ublas::matrix<X>() { }

    //! Construct a matrix with \a r rows and \a c columns with values initialised to zero.
    Matrix(size_t r, size_t c)
        : ublas::matrix<X>(r,c) { for(size_t i=0; i!=r; ++i) { for(size_t j=0; j!=c; ++j) { (*this)(i,j)=0; } } }
    
    //! Construct a matrix with \a r rows and \a c columns, with values initialised from the C-style array beginning at \a ptr in row-major format. The value in the \a i<sup>th</sup> row and \a j<sup>th</sup> column of the resulting matrix is \a ptr[i*c+j].
    template<class XX> Matrix(size_t r, size_t c, const XX* ptr)
        : ublas::matrix<X>(r,c) { for(size_t i=0; i!=r; ++i) { for(size_t j=0; j!=c; ++j) { (*this)(i,j)=ptr[i*c+j]; } } }
    
    //! Construct a matrix with \a r rows and \a c columns, with values initialised from the C-style array beginning at \a ptr. The C-style array stores a matrix with row increment \a ri and column increment \a ci.
    template<class XX> Matrix(size_t r, size_t c, const XX* ptr, int ri, int ci)
        : ublas::matrix<X>(r,c) { for(size_t i=0; i!=r; ++i) { for(size_t j=0; j!=c; ++j) { (*this)(i,j)=ptr[i*ri+j*ci]; } } }
    //! Construct a matrix with \a r rows and \a c columns, with values initialised from a variadic argument list. WARNING: The values in the list must all be double-precision type; in particular, constants must be floating-point values \c 2.0 rather integer values \c 2 .
    Matrix(size_t r, size_t c, const double& x0, ... );
    //! Construct a matrix from a string in Matlab format, with entries in a row separated with commas, and rows separated with semicolons. e.g. <tt>"[a00, a01, a02; a10, a11, a12]"</tt>.
    Matrix(const std::string& str)
        : ublas::matrix<X>() { std::stringstream ss(str); ss >> *this; }

    template<class AE> Matrix(const ublas::matrix_expression<AE> &ae)
        : ublas::matrix<X>(ae) { }
    template<class AE> Matrix<X>& operator=(const ublas::matrix_expression<AE> &ae) { this->ublas::matrix<X>::operator=(ae); return *this; }
    //@}

    //@{
    //! \name Comparison operators
    
    //! \brief Equality operator. Allows comparison with a matrix using another real type.
    template<class XX> bool operator==(const Matrix<XX>& A) const;
    //@}

    //@{
    //! \name Data access
    
  
    //! \brief The number of rows of the matrix.
    size_t row_size() const { return this->size1(); }
    //! \brief The number of columns of the matrix.
    size_t column_size() const { return this->size2(); }
    //! \brief Get the value stored in the \a i<sup>th</sup> row and \a j<sup>th</sup> column.
    const X& get(size_t i, size_t j) const { return (*this)[i][j]; }
    //! \brief Set the value stored in the \a i<sup>th</sup> row and \a j<sup>th</sup> column to \a x.
    template<class T> void set(size_t i, size_t j, const T& x) { (*this)[i][j] = x; }

#ifdef DOXYGEN
    //! \brief C-style subscripting operator.
    X& operator[][](size_t i, size_t j);
    //! \brief C-style constant subscripting operator.
    const X& operator[][](size_t i, size_t j) const;
#else
    const X* operator[](size_t r) const { return &this->operator()(r,0); }
    X* operator[](size_t r) { return &this->operator()(r,0); }
#endif
    //@}

    //@{
    //! \name Static constructors

    //! \brief The zero matrix with \a r rows and \a c columns.
    static Matrix<X> zero(size_t r, size_t c) { return Matrix<X>(r,c); }
    //! \brief The itentity matrix with \a n rows and \a n columns.
    static Matrix<X> identity(size_t n) { Matrix<X> I(n,n); for(size_t i=0; i!=n; ++i) { I[i][i]=1; } return I; }
    //@}
    

#ifdef DOXYGEN
    //! \brief The supremum norm of the matrix \a A, defined as \f$||A||_\infty=\max_i \sum_j |A_{ij}|\f$.
    friend template<class X> X norm(const Matrix<X>& A);

     //! \brief %Matrix negation.
    friend template<class X> Matrix<X> operator-(const Matrix<X>& A);
    //! \brief %Matrix addition.
    friend template<class X> Matrix<X> operator+(const Matrix<X>& A1, const Matrix<X>& A2);
    //! \brief %Matrix subtraction.
    friend template<class X> Matrix<X> operator-(const Matrix<X>& A1, const Matrix<X>& A2);
    //! \brief %Scalar multiplication.
    friend template<class X> Matrix<X> operator*(const X& s, const Matrix<X>& A);
    //! \brief %Scalar multiplication.
    friend template<class X> Matrix<X> operator*(const Matrix<X>& A, const X& s);
    //! \brief %Scalar division.
    friend template<class X> Matrix<X> operator/(const Matrix<X>& A, const X& s);

    //! \brief %Matrix-vector multiplication.
    friend template<class X> Vector<X> operator*(const Matrix<X>& A, const Vector<X>& v);
    //! \brief %Matrix-matrix multiplication.
    friend template<class X> Vector<X> operator*(const Matrix<X>& A1, const Matrix<X>& A2);

    //! \brief Write to an output stream.
    friend template<class X> std::ostream& operator<<(std::ostream& os, const Matrix<X>& A);
    //! \brief Read from an output stream.
    friend template<class X> std::istream& operator>>(std::istream& is, Matrix<X>& A);
#endif

    };

template<class X> X norm(const Matrix<X>& A);

template<class X> std::ostream& operator<<(std::ostream& os, const Matrix<X>& A);
template<class X> std::istream& operator>>(std::istream& is, Matrix<X>& A);

Matrix<Float> inverse(const Matrix<Float>& A);
Matrix<Interval> inverse(const Matrix<Interval>& A);

#ifdef HAVE_GMPXX_H
Matrix<Rational> inverse(const Matrix<Rational>& A);
#endif // HAVE_GMPXX_H

Matrix<Float> triangular_multiplier(const Matrix<Float>& A);
tuple< Matrix<Float>, Matrix<Float> > orthogonal_decomposition(const Matrix<Float>&);

Matrix<Float> midpoint(const Matrix<Interval>&);


/*
template<class X> inline Vector<X> operator*(const Matrix<X>& A, const Vector<X>& v) {
    return ublas::prod(A,v);
}
 
template<class X> inline Matrix<X> operator*(const Matrix<X>& A, const Matrix<X>& B) {
    return ublas::prod(A,B);
}
 

template<class E1, class E2>
inline
typename ublas::matrix_vector_binary1_traits<typename E1::value_type, E1,
                                             typename E2::value_type, E2>::result_type
operator* (const ublas::matrix_expression<E1> &e1,
           const ublas::vector_expression<E2> &e2) 
{
    typedef typename ublas::matrix_vector_binary1_traits<typename E1::value_type, E1,
    typename E2::value_type, E2>::expression_type expression_type;
    return expression_type (e1 (), e2 ());
}


template<class E1, class E2>
inline
typename ublas::matrix_matrix_binary_traits<typename E1::value_type, E1,
                                            typename E2::value_type, E2>::result_type
operator* (const ublas::matrix_expression<E1> &e1,
           const ublas::matrix_expression<E2> &e2) 
{
    typedef typename ublas::matrix_matrix_binary_traits<typename E1::value_type, E1,
        typename E2::value_type, E2>::expression_type expression_type;
    return expression_type (e1 (), e2 ());
}

template<class E1, class T2>
inline
typename ublas::matrix_binary_scalar2_traits<E1, const T2, ublas::scalar_multiplies<typename E1::value_type, T2> >::result_type
operator* (const ublas::matrix_expression<E1> &e1,
           const T2 &e2) {
    typedef typename ublas::matrix_binary_scalar2_traits<E1, const T2, 
        ublas::scalar_multiplies<typename E1::value_type, T2> >::expression_type expression_type;
    return expression_type (e1 (), e2);
}

*/


template<class X> Matrix<X>::Matrix(size_t r, size_t c, const double& x0, ...) 
    : ublas::matrix<X>(r,c) 
{
    assert(r>=1 && c>=1); va_list args; va_start(args,x0);
    (*this)[0][0]=x0; 
    for(size_t j=1; j!=c; ++j) { (*this)[0][j]=va_arg(args,double); } 
    for(size_t i=1; i!=r; ++i) { for(size_t j=0; j!=c; ++j) { (*this)[i][j]=va_arg(args,double); } }
    va_end(args);
}

template<class X> template<class XX> bool Matrix<X>::operator==(const Matrix<XX>& A2) const 
{
    const Matrix<X>& A1=*this;
    if(A1.row_size()!=A2.row_size() || A1.column_size() != A2.column_size()) {
        return false;
    }
    for(size_t i=0; i!=A1.row_size(); ++i) {
        for(size_t j=0; j!=A1.column_size(); ++j) {
            if(A1[i][j]!=A2[i][j]) {
                return false;
            }
        }
    }
    return true;
}


template<class X> X norm(const Matrix<X>& A) 
{
    X result=0;
    for(size_t i=0; i!=A.row_size(); ++i) {
        X row_sum=0;
        for(size_t j=0; j!=A.column_size(); ++j) {
            row_sum+=abs(A[i][j]);
        }
        if(row_sum>result) { result=row_sum; }
    }
    return result;
}



template<class X> std::ostream& operator<<(std::ostream& os, const Matrix<X>& A) {
    if(A.row_size()==0 || A.column_size()==0) { os << "["; }
    for(size_t i=0; i!=A.row_size(); ++i) { 
        for(size_t j=0; j!=A.column_size(); ++j) { 
            os << (j==0 ? (i==0 ? "[" : "; ") : ",") << A(i,j); } }
    return os << "]";
}


template<class X> std::istream& operator>>(std::istream& is, Matrix<X>& A) {
    char c;
    is >> c;
    is.putback(c);
    if(c=='[') {
        is >> c;
        /* Representation as a literal [a11,a12,...,a1n; a21,a22,...a2n; ... ; am1,am2,...,amn] */
        std::vector< std::vector<X> > v;
        X x;
        c=';';
        while(is && c==';') {
            v.push_back(std::vector<X>());
            c=',';
            while(is && c==',') {
                is >> x;
                v.back().push_back(x);
                is >> c;
            }
        }
        if(is) {
            A=Matrix<X>(v.size(),v.front().size());
            for(size_t i=0; i!=A.row_size(); ++i) {
                if(v[i].size()!=A.column_size()) {
                    ARIADNE_THROW(InvalidInput,"Matrix::read(istream&)","row[0].size()="<<v[0].size()<<", row["<<i<<"].size()="<<v[i].size());
                }
                for(size_t j=0; j!=A.column_size(); ++j) {
                    A(i,j)=v[i][j];
                }
            }
        }
    }
    else {
        ARIADNE_THROW(InvalidInput,"Matrix::read(istream&)"," separator c="<<c);
    }
    return is;
}


} // namespace Ariadne

#endif
