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

namespace Ariadne {

/// A matrix over a field.
template<class X>
class Matrix
    : public boost::numeric::ublas::matrix<X>
{
  public:
    Matrix()
        : boost::numeric::ublas::matrix<X>() { }
    template<class AE> Matrix(const boost::numeric::ublas::matrix_expression<AE> &ae)
        : boost::numeric::ublas::matrix<X>(ae) { }
    Matrix(size_t r, size_t c)
        : boost::numeric::ublas::matrix<X>(r,c) { for(size_t i=0; i!=r; ++i) { for(size_t j=0; j!=c; ++j) { (*this)(i,j)=0; } } }
    template<class XX> Matrix(size_t r, size_t c, const XX* ptr)
        : boost::numeric::ublas::matrix<X>(r,c) { for(size_t i=0; i!=r; ++i) { for(size_t j=0; j!=c; ++j) { (*this)(i,j)=ptr[i*c+j]; } } }
    template<class XX> Matrix(size_t r, size_t c, const XX* ptr, int ri, int ci)
        : boost::numeric::ublas::matrix<X>(r,c) { for(size_t i=0; i!=r; ++i) { for(size_t j=0; j!=c; ++j) { (*this)(i,j)=ptr[i*ri+j*ci]; } } }
    Matrix(size_t r, size_t c, const Float& x0, ... );
    size_t row_size() const { return this->size1(); }
    size_t column_size() const { return this->size2(); }
    template<class XX> bool operator==(const Matrix<XX>& mx) const;
    const X* operator[](size_t r) const { return &this->operator()(r,0); }
    X* operator[](size_t r) { return &this->operator()(r,0); }
    const X& get(size_t i, size_t j) const { return (*this)[i][j]; }
    template<class T> void set(size_t i, size_t j, const T& x) { (*this)[i][j] = x; }
  
    static Matrix<X> zero(size_t r, size_t c) { return Matrix<X>(r,c); }
    static Matrix<X> identity(size_t n) { Matrix<X> I(n,n); for(size_t i=0; i!=n; ++i) { I[i][i]=1; } return I; }

};

template<class X> Matrix<X>::Matrix(size_t r, size_t c, const Float& x0, ...) 
    : ublas::matrix<X>(r,c) 
{
    assert(r>=1 && c>=1); va_list args; va_start(args,x0);
    (*this)[0][0]=x0; 
    for(size_t j=1; j!=c; ++j) { (*this)[0][j]=va_arg(args,Float); } 
    for(size_t i=1; i!=r; ++i) { for(size_t j=0; j!=c; ++j) { (*this)[i][j]=va_arg(args,Float); } }
    va_end(args);
}

template<class X> Matrix<X> make_matrix(const std::string& str) {
    Matrix<X> res;
    std::stringstream ss(str);
    ss >> res;
    return res;
}
  

template<class X> Matrix<X> inverse(const Matrix<X>& A);
template<class X> X norm(const Matrix<X>& A);

/*
template<class X> inline Vector<X> operator*(const Matrix<X>& A, const Vector<X>& v) {
    return boost::numeric::ublas::prod(A,v);
}
 
template<class X> inline Matrix<X> operator*(const Matrix<X>& A, const Matrix<X>& B) {
    return boost::numeric::ublas::prod(A,B);
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



template<class X> Matrix<X> inverse(const Matrix<X>& A) 
{
    assert(A.row_size()==A.column_size());
    if(A.row_size()==1) {
        Matrix<X> I(1,1);
        I[0][0]=1/A[0][0];
        return I;
    }
    if(A.row_size()==2) {
        Matrix<X> I(2,2);
        X rdet = 1/(A[0][0]*A[1][1] - A[0][1] * A[1][0]);
        I[0][0]=A[1][1]*rdet;
        I[0][1]=-A[0][1]*rdet;
        I[1][0]=-A[1][0]*rdet;
        I[1][1]=A[0][0]*rdet;
        return I;
    }
    assert(false);
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


Matrix<Float> midpoint(const Matrix<Interval>& A);



template<class X> std::ostream& operator<<(std::ostream& os, const Matrix<X>& A) {
    if(A.row_size()==0 || A.column_size()==0) { os << '['; }
    for(size_t i=0; i!=A.row_size(); ++i) { 
        for(size_t j=0; j!=A.column_size(); ++j) { 
            os << (j==0 ? (i==0 ? '[' : ';') : ',') << A(i,j); } }
    return os << ']';
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
