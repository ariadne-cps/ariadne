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

#include <iostream>
#include <iomanip>
#include <initializer_list>
#include <limits>

#include "macros.h"
#include "exceptions.h"
#include "numeric.h"
#include "tuple.h"

namespace Ariadne {

template<class T> class Pretty {
    const T& _t;
  public:
    Pretty(const T& t) : _t(t) { }
    operator const T& () const { return this->_t; }
};
template<class T> Pretty<T> pretty(const T& t) { return Pretty<T>(t); }
template<class T> std::ostream& operator<<(std::ostream& os, const Pretty<T>& p);

class SingularMatrixException {
};

template<class X> class Matrix;
typedef Matrix<Float> FloatMatrix;
typedef Matrix<Interval> IntervalMatrix;
typedef Matrix<Dyadic> DyadicMatrix;

struct Uninitialised { };

template<class M> struct MatrixExpression {
    const M& operator()() const { return static_cast<const M&>(*this); }
};

template<class M> struct MatrixContainer : public MatrixExpression<M> { };

template<class X> struct MatrixRow {
    Matrix<X>& _A; size_t _i; MatrixRow(Matrix<X>& A,size_t i) : _A(A), _i(i) { }
    const X& operator[](size_t j) const { return _A.get(_i,j); }
    X& operator[](size_t j) { return _A.at(_i,j); }
};
template<class X> struct ConstMatrixRow {
    const Matrix<X>& _A; size_t _i; ConstMatrixRow(const Matrix<X>& A,size_t i) : _A(A), _i(i) { }
    const X& operator[](size_t j) const { return _A.get(_i,j); }
};


//! \ingroup LinearAlgebraModule
//! \brief A matrix over a field. See also \link Ariadne::Vector \c Vector<X> \endlink.
//!
//! \par Python interface
//!
//! In the Python interface, %Ariadne matrices can be constructed from Python literals of the form <br>
//! \c   [[a11,a12,...,a1n],[a21,a22,...,a2n],...,[am1,am2,...,amn]].
//!
template<class X>
class Matrix
    : public MatrixContainer< Matrix<X> >
{
    size_t _nr; size_t _nc; X* _ptr;
  public:
    //@{
    //! \name Type definitions.

    typedef X ValueType;

    //@}

    //@{
    //! \name Constructors

    //! Destructor
    ~Matrix() { delete[] _ptr; }

    //! Default constructor makes a \f$0\times0\f$ matrix.
    Matrix() : _nr(0), _nc(0), _ptr(0) { }

    //! Construct a matrix with \a r rows and \a c columns with values uninitialised.
    Matrix(size_t r, size_t c, Uninitialised) : _nr(r), _nc(c), _ptr(new X[_nr*_nc]) { }

    //! Construct a matrix with \a r rows and \a c columns with values initialised to zero.
    Matrix(size_t r, size_t c) : _nr(r), _nc(c), _ptr(new X[_nr*_nc]) {
        for(size_t k=0; k!=r*c; ++k) { this->_ptr[k]=0; } }

    //! Construct a matrix with \a r rows and \a c columns with values initialised to \a x.
    Matrix(size_t r, size_t c, const X& x) : _nr(r), _nc(c), _ptr(new X[_nr*_nc]) {
        for(size_t k=0; k!=r*c; ++k) { this->_ptr[k]=x; } }

    //! Construct a matrix with \a r rows and \a c columns, with values initialised from the C-style Array beginning at \a ptr in row-major format. The value in the \a i<sup>th</sup> row and \a j<sup>th</sup> column of the resulting matrix is \a ptr[i*c+j].
    template<class XX> Matrix(size_t r, size_t c, const XX* ptr) : _nr(r), _nc(c), _ptr(new X[_nr*_nc]) {
        for(size_t k=0; k!=r*c; ++k) { this->_ptr[k] = ptr[k]; } }

    //! Construct a matrix with \a r rows and \a c columns, with values initialised from the C-style Array beginning at \a ptr. The C-style Array stores a matrix with row increment \a ri and column increment \a ci.
    template<class XX> Matrix(size_t r, size_t c, const XX* ptr, int ri, int ci) : _nr(r), _nc(c), _ptr(new X[_nr*_nc]) {
        for(size_t i=0; i!=r; ++i) { for(size_t j=0; j!=c; ++j) { this->set(i,j,ptr[i*ri+j*ci]); } } }
    //! Construct a matrix using initializer lists.
    Matrix(std::initializer_list< std::initializer_list<X> > const lst);
    //! Construct a matrix from a string in Matlab format, with entries in a row separated with commas, and rows separated with semicolons. e.g. <tt>"[a00, a01, a02; a10, a11, a12]"</tt>.
    explicit Matrix(const std::string& str)
        : _nr(0), _nc(0), _ptr(0) { std::stringstream ss(str); ss >> *this; }

    //! \brief Copy constructor.
    Matrix(const Matrix<X>& other) : _nr(other._nr), _nc(other._nc), _ptr(new X[_nr*_nc]) {
        for(size_t k=0; k!=_nr*_nc; ++k) { this->_ptr[k] = other._ptr[k]; } }

    //! \brief Copy assignment.
    Matrix<X>& operator=(const Matrix<X>& other) {
        if(this!=&other) {
            this->resize(other.row_size(),other.column_size());
            for(size_t k=0; k!=_nr*_nc; ++k) { this->_ptr[k] = other._ptr[k]; }
        }
        return *this;
    }

    template<class ME> Matrix(const MatrixExpression<ME>& Ae) : _nr(Ae().row_size()), _nc(Ae().column_size()), _ptr(new X[_nr*_nc]) {
        for(size_t i=0; i!=this->_nr; ++i) { for(size_t j=0; j!=this->_nc; ++j) { this->set(i,j,static_cast<X>(Ae().get(i,j))); } } }
    template<class ME> Matrix<X>& operator=(const MatrixExpression<ME>& Ae) {
        this->resize(Ae().row_size(),Ae().column_size());
        for(size_t i=0; i!=this->_nr; ++i) { for(size_t j=0; j!=this->_nc; ++j) { this->set(i,j,Ae().get(i,j)); } } return *this; }
    //@}

    //@{
    //! \name Data access

    void _check_data_access(size_t i, size_t j) const {
        ARIADNE_PRECONDITION_MSG(i<this->row_size()&&j<this->column_size(),"A="<<*this<<" i="<<i<<" j="<<j); }

    void resize(size_t r, size_t c) {
        if(this->_nr*this->_nc!=r*c) { X* new_ptr=new X[r*c]; delete[] this->_ptr; this->_ptr=new_ptr; } this->_nr=r; this->_nc=c; }

    //! \brief The number of rows of the matrix.
    size_t row_size() const { return this->_nr; }
    //! \brief The number of columns of the matrix.
    size_t column_size() const { return this->_nc; }
    //! \brief Get the value stored in the \a i<sup>th</sup> row and \a j<sup>th</sup> column.
    X& at(size_t i, size_t j) const { this->_check_data_access(i,j); return this->_ptr[i*_nc+j]; }
    //! \brief Get the value stored in the \a i<sup>th</sup> row and \a j<sup>th</sup> column.
    const X& get(size_t i, size_t j) const { this->_check_data_access(i,j); return this->_ptr[i*_nc+j]; }
    //! \brief Set the value stored in the \a i<sup>th</sup> row and \a j<sup>th</sup> column to \a x.
    template<class T> void set(size_t i, size_t j, const T& x) { this->_check_data_access(i,j); this->_ptr[i*_nc+j] = x; }
    //! \brief A pointer to the first element of the data storage.
    const X* begin() const { return _ptr; }

#ifdef DOXYGEN
    //! \brief C-style subscripting operator.
    X& operator[][](size_t i, size_t j);
    //! \brief C-style constant subscripting operator.
    const X& operator[][](size_t i, size_t j) const;
#else
    ConstMatrixRow<X> operator[](size_t i) const { return ConstMatrixRow<X>(*this,i); }
    MatrixRow<X> operator[](size_t i) { return MatrixRow<X>(*this,i); }
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

     //! \brief %Matrix unary plus.
    friend template<class X> Matrix<X> operator+(const Matrix<X>& A);
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

    //! \brief The operator norm using the supremum norm on Euclidean space.
    //! Given explicitly by \f$\max_i \sum_j |a_{ij}|\f$.
    friend template<class X> X norm(const Matrix<X>& A);

    //! \brief The transpose matrix.
    friend template<class X> Matrix<X> transpose(const Matrix<X>& A);
    //! \brief The inverse matrix.
    friend template<class X> Matrix<X> inverse(const Matrix<X>& A);

    //! \brief Solve a system of linear equations.
    //!
    //! Uses Gaussian elimination with column row pivoting for floating-point matrices.
    //!
    //! Uses Gauss-Seidel iteration on a preconditioned system for interval matrices.
    //! First, an approximate inverse \f$S\f$ to \f$m([A])\f$ is computed using floating-point arithmetic.
    //! Next, the system is reformulated as \f$S[A][x]=S[b]\f$, or \f$[J][x]=[c]\f$
    //! The interval matrix \f$S[A]\f$ should be close to the identity
    //! (if the interval matrix \f$[A]\f$ is well-conditioned
    //! and the radii of the coefficients are sufficiently small).
    //!
    //! If (the preconditioned) \f$[A]\f$ is diagonally-dominant, then an initial bound for the solutions to
    //! to \f$[A][x]=[b]\f$ can be given as \f$|x_i|\leq |b_i|/(|a_{ii}|-\sum_{j\neq i}|a_{ij}|)\f$.
    //!
    //! Any bound \f$[x]\f$ for the solutions can then be iteratively updated using the Gauss-Seidel scheme
    //! \f$x_i' = x_i\cap (b_i-\sum_{j\neq i} a_{ij}x_j)/a_{ii}\f$.
    friend template<class X> Vector<X> solve(const Matrix<X>& A, const Vector<X>& b);
    //! \brief Simultaneously solve multiple linear equations.
    friend template<class X> Matrix<X> solve(const Matrix<X>& A, const Matrix<X>& b);

    //! \brief The midpoint of an interval matrix.
    friend Matrix<Float> midpoint(const Matrix<Interval>& A);

    //! \brief Write to an output stream.
    friend template<class X> std::ostream& operator<<(std::ostream& os, const Matrix<X>& A);
    //! \brief Read from an output stream.
    friend template<class X> std::istream& operator>>(std::istream& is, Matrix<X>& A);


#endif

};

template<class X> std::ostream& operator<<(std::ostream& os, const Matrix<X>& A);
template<class X> std::istream& operator>>(std::istream& is, Matrix<X>& A);

template<class M> std::ostream& operator<<(std::ostream& os, const MatrixExpression<M>& Ae);


template<class X> Matrix<X>::Matrix(std::initializer_list< std::initializer_list<X> > lst)
    : _nr(lst.size()), _nc(lst.size()==0?0u:lst.begin()->size()), _ptr(new X[_nr*_nc])
{
    X* pointer = this->_ptr;
    for(typename std::initializer_list< std::initializer_list<X> >::const_iterator row_iter=lst.begin();
        row_iter!=lst.end(); ++row_iter)
    {
        assert(row_iter->size()==this->_nc);
        for(typename std::initializer_list<X>::const_iterator iter=row_iter->begin();
            iter!=row_iter->end(); ++iter)
        {
            *pointer = *iter;
            ++pointer;
        }
    }
}

template<class X1, class X2> bool operator==(const Matrix<X1>& A1, const Matrix<X2>& A2)
{
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

template<class M1, class M2> bool operator==(const MatrixExpression<M1>& A1e, const MatrixExpression<M2>& A2e) {
    return Matrix<typename M1::ValueType>(A1e) == Matrix<typename M2::ValueType>(A2e);
}



template<class M> struct MatrixContainerRange
    : public MatrixContainer< MatrixContainerRange<M> >
{
    typedef typename M::ValueType ValueType;
    M& _A; Range _rng1; Range _rng2;
    MatrixContainerRange(M& A, Range rng1, Range rng2) : _A(A), _rng1(rng1), _rng2(rng2) { }
    size_t row_size() const { return _rng1.size(); }
    size_t column_size() const { return _rng2.size(); }
    ValueType get(size_t i, size_t j) const { return _A.get(i+_rng1.start(),j+_rng2.start()); }
    void set(size_t i, size_t j, const ValueType& x) { _A.set(i+_rng1.start(),j+_rng2.start(),x); }
    template<class ME> MatrixContainerRange<M>& operator=(const MatrixExpression<ME>& Ae) {
        ARIADNE_PRECONDITION(this->row_size()==Ae().row_size() && this->column_size()==Ae().column_size());
        for(size_t i=0; i!=this->row_size(); ++i) { for(size_t j=0; j!=this->column_size(); ++j) { this->set(i,j,Ae().get(i,j)); } }
        return *this; }
};

template<class M> struct MatrixRange
    : public MatrixExpression< MatrixRange<M> >
{
    typedef typename M::ValueType ValueType;
    const M& _A; Range _rng1; Range _rng2;
    MatrixRange(const M& A, Range rng1, Range rng2) : _A(A), _rng1(rng1), _rng2(rng2) { }
    size_t row_size() const { return _rng1.size(); }
    size_t column_size() const { return _rng2.size(); }
    ValueType get(size_t i, size_t j) const { return _A.get(i+_rng1.start(),j+_rng2.start()); }
};

template<class M> struct MatrixNegation
    : public MatrixExpression< MatrixNegation<M> >
{
    typedef typename M::ValueType ValueType;
    const M& _A;
    MatrixNegation(const M& A) : _A(A) { }
    size_t row_size() const { return _A.row_size(); }
    size_t column_size() const { return _A.column_size(); }
    ValueType get(size_t i,size_t j) const { return -_A.get(i,j); }
};

template<class M1, class M2> struct MatrixSum
    : public MatrixExpression< MatrixSum<M1,M2> >
{
    typedef typename Arithmetic<typename M1::ValueType, typename M2::ValueType>::ResultType ValueType;
    const M1& _A1; const M2& _A2;
    MatrixSum(const M1& A1, const M2& A2) : _A1(A1), _A2(A2) { }
    size_t row_size() const { return _A1.row_size(); }
    size_t column_size() const { return _A1.column_size(); }
    ValueType get(size_t i,size_t j) const { return _A1.get(i,j)+_A2.get(i,j); }
};

template<class M1, class M2> struct MatrixDifference
    : public MatrixExpression< MatrixDifference<M1,M2> >
{
    typedef typename Arithmetic<typename M1::ValueType, typename M2::ValueType>::ResultType ValueType;
    const M1& _A1; const M2& _A2;
    MatrixDifference(const M1& A1, const M2& A2) : _A1(A1), _A2(A2) { }
    size_t row_size() const { return _A1.row_size(); }
    size_t column_size() const { return _A1.column_size(); }
    ValueType get(size_t i,size_t j) const { return _A1.get(i,j)-_A2.get(i,j); }
};

template<class M1, class X2> struct MatrixScalarProduct
    : public MatrixExpression< MatrixScalarProduct<M1,X2> >
{
    typedef typename Arithmetic<typename M1::ValueType, X2>::ResultType ValueType;
    const M1& _A1; const X2& _x2;
    MatrixScalarProduct(const M1& A1, const X2& x2) : _A1(A1), _x2(x2) { }
    size_t row_size() const { return _A1.row_size(); }
    size_t column_size() const { return _A1.column_size(); }
    ValueType get(size_t i,size_t j) const { return _A1.get(i,j)*_x2; }
};

template<class M1, class X2> struct MatrixScalarQuotient
    : public MatrixExpression< MatrixScalarProduct<M1,X2> >
{
    typedef typename Arithmetic<typename M1::ValueType, X2>::ResultType ValueType;
    const M1& _A1; const X2& _x2;
    MatrixScalarQuotient(const M1& A1, const X2& x2) : _A1(A1), _x2(x2) { }
    size_t row_size() const { return _A1.row_size(); }
    size_t column_size() const { return _A1.column_size(); }
    ValueType get(size_t i,size_t j) const { return _A1.get(i,j)/_x2; }
};

template<class M> struct MatrixTranspose
    : public MatrixExpression< MatrixTranspose<M> >
{
    const M& _A;
    typedef typename M::ValueType ValueType;
    MatrixTranspose(const M& A) : _A(A) { }
    size_t row_size() const { return _A.column_size(); }
    size_t column_size() const { return _A.row_size(); }
    ValueType get(size_t i, size_t j) const { return _A.get(j,i); }
};



template<class M> inline
const M&
operator+(const MatrixExpression<M>& Ae) { return Ae(); }

template<class M> inline
MatrixNegation<M>
operator-(const MatrixExpression<M>& Ae) { return MatrixNegation<M>(Ae()); }

template<class M1, class M2> inline
typename EnableIf< IsDefined<typename Arithmetic<typename M1::ValueType,typename M2::ValueType>::ResultType>, MatrixSum< M1, M2 > >::Type
operator+(const MatrixExpression<M1>& A1e, const MatrixExpression<M2>& A2e) { return MatrixSum<M1,M2>(A1e(),A2e()); }

template<class M1, class M2> inline
typename EnableIf< IsDefined<typename Arithmetic<typename M1::ValueType,typename M2::ValueType>::ResultType>, MatrixDifference< M1, M2 > >::Type
operator-(const MatrixExpression<M1>& A1e, const MatrixExpression<M2>& A2e) { return MatrixDifference<M1,M2>(A1e(),A2e()); }

template<class X1, class M2> inline
typename EnableIf< IsDefined<typename Arithmetic<X1,typename M2::ValueType>::ResultType>, MatrixScalarProduct< M2, X1 > >::Type
operator*(const X1& x1, const MatrixExpression<M2>& A2e) { return MatrixScalarProduct<M2,X1>(A2e(),x1); }

template<class M1, class X2> inline
typename EnableIf< IsDefined<typename Arithmetic<typename M1::ValueType,X2>::ResultType>, MatrixScalarProduct< M1, X2 > >::Type
operator*(const MatrixExpression<M1>& A1e, const X2& x2) { return MatrixScalarProduct<M1,X2>(A1e(),x2); }

template<class M1, class X2> inline
typename EnableIf< IsDefined<typename Arithmetic<typename M1::ValueType,X2>::ResultType>, MatrixScalarQuotient< M1, X2 > >::Type
operator/(const MatrixExpression<M1>& A1e, const X2& x2) { return MatrixScalarQuotient<M1,X2>(A1e(),x2); }


template<class M, class V>
Vector<typename Arithmetic<typename M::ValueType,typename V::ValueType>::ResultType>
operator*(const MatrixExpression<M>& Ae, const VectorExpression<V>& ve) {
    const M& A=Ae();
    const V& v=ve();
    ARIADNE_PRECONDITION(A.column_size()==v.size());
    Vector<typename Arithmetic<typename M::ValueType,typename V::ValueType>::ResultType> r(A.row_size());
    for(size_t j=0; j!=v.size(); ++j) {
        typename V::ValueType vj = v[j];
        for(size_t i=0; i!=r.size(); ++i) {
            r[i] += A.get(i,j) * vj;
        }
    }
    return r;
}

template<class V, class M>
Vector<typename Arithmetic<typename V::ValueType,typename M::ValueType>::ResultType>
operator*(const VectorExpression<V>& ve, const MatrixExpression<M>& Ae) {
    const V& v=ve();
    const M& A=Ae();
    ARIADNE_PRECONDITION(A.row_size()==v.size());
    Vector<typename Arithmetic<typename V::ValueType,typename M::ValueType>::ResultType> r(A.column_size());
    for(size_t i=0; i!=v.size(); ++i) {
        typename V::ValueType vi = v[i];
        for(size_t j=0; j!=r.size(); ++j) {
            r[j] += vi * A.get(i,j);
        }
    }
    return r;
}

template<class X1, class X2>
Matrix<typename Arithmetic<X1,X2>::ResultType>
operator*(const Matrix<X1>& A1, const Matrix<X2>& A2) {
    ARIADNE_PRECONDITION(A1.column_size()==A2.row_size());
    Matrix<typename Arithmetic<X1,X2>::ResultType> r(A1.row_size(),A2.column_size());
    for(size_t k=0; k!=A2.row_size(); ++k) {
        for(size_t j=0; j!=A2.column_size(); ++j) {
            for(size_t i=0; i!=A1.row_size(); ++i) {
                r.at(i,j) += A1.get(i,k) * A2.get(k,j);
            }
        }
    }
    return r;
}

template<class M1, class M2> inline
Matrix<typename Arithmetic<typename M1::ValueType,typename M2::ValueType>::ResultType>
operator*(const MatrixExpression<M1>& A1e, const MatrixExpression<M2>& A2e) {
    return Matrix<typename M1::ValueType>(A1e) * Matrix<typename M2::ValueType>(A2e);
}

template<class X, class M> inline Matrix<X>& operator+=(Matrix<X>& A, const MatrixExpression<M>& Be) {
    const M& B=Be();
    ARIADNE_PRECONDITION(A.row_size()==B.row_size() && A.column_size()==B.column_size());
    for(size_t i=0; i!=A.row_size(); ++i) { for(size_t j=0; j!=A.column_size(); ++j) { A.at(i,j)+=B.get(i,j); } }
    return A;
}

template<class X, class M> inline Matrix<X>& operator-=(Matrix<X>& A, const MatrixExpression<M>& Be) {
    const M& B=Be();
    ARIADNE_PRECONDITION(A.row_size()==B.row_size() && A.column_size()==B.column_size());
    for(size_t i=0; i!=A.row_size(); ++i) { for(size_t j=0; j!=A.column_size(); ++j) { A.at(i,j)-=B.get(i,j); } }
    return A;
}

template<class X1, class X2> inline typename EnableIfNumeric<X2,Matrix<X1>&>::Type operator*=(Matrix<X1>& A, const X2& x) {
    for(size_t i=0; i!=A.row_size(); ++i) { for(size_t j=0; j!=A.column_size(); ++j) { A.at(i,j)*=x; } }
    return A;
}

template<class X1, class X2> inline typename EnableIfNumeric<X2,Matrix<X1>&>::Type operator/=(Matrix<X1>& A, const X2& x) {
    for(size_t i=0; i!=A.row_size(); ++i) { for(size_t j=0; j!=A.column_size(); ++j) { A.at(i,j)/=x; } }
    return A;
}

template<class M> inline MatrixRange<M> project(const MatrixExpression<M>& Ae, Range rw_rng, Range cl_rng) {
    return MatrixRange<M>(Ae(),rw_rng,cl_rng); }

template<class X> inline MatrixContainerRange< Matrix<X> > project(Matrix<X>& A, Range rw_rng, Range cl_rng) {
    return MatrixContainerRange< Matrix<X> >(A,rw_rng,cl_rng); }

template<class M> inline MatrixTranspose<M> transpose(const MatrixExpression<M>& Ae) {
    return MatrixTranspose<M>(Ae());
}



template<class X> X norm(const Matrix<X>& A)
{
    X result=0;
    for(size_t i=0; i!=A.row_size(); ++i) {
        X row_sum=0;
        for(size_t j=0; j!=A.column_size(); ++j) {
            row_sum+=abs(A[i][j]);
        }
        // NOTE: The arguments must be this way round to propagate a nan row_sum
        result=max(row_sum,result);
    }
    return result;
}

template<class X> Vector<X> row(const Matrix<X>& A, size_t i)
{
    Vector<X> r(A.column_size());
    for(size_t j=0; j!=r.size(); ++j) { r[j]=A.get(i,j); }
    return r;
}

template<class X> Vector<X> column(const Matrix<X>& A, size_t j)
{
    Vector<X> r(A.row_size());
    for(size_t i=0; i!=r.size(); ++i) { r[i]=A.get(i,j); }
    return r;
}



template<class M> std::ostream& operator<<(std::ostream& os, const MatrixExpression<M>& Ae) {
    return os << Matrix<typename M::ValueType>(Ae);
}

template<class X> std::ostream& operator<<(std::ostream& os, const Matrix<X>& A) {
    if(A.row_size()==0 || A.column_size()==0) { os << "["; }
    for(size_t i=0; i!=A.row_size(); ++i) {
        for(size_t j=0; j!=A.column_size(); ++j) {
            os << (j==0 ? (i==0 ? "[" : "; ") : ",") << A.get(i,j); } }
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
                    A[i][j]=v[i][j];
                }
            }
        }
    }
    else {
        ARIADNE_THROW(InvalidInput,"Matrix::read(istream&)"," separator c="<<c);
    }
    return is;
}



class PivotMatrix;
template<class X> class PLUMatrix;
template<class X> class QRMatrix;

template<class X> Matrix<X> inverse(const Matrix<X>& A);
template<class X> Matrix<X> solve(const Matrix<X>& A, const Matrix<X>& B);
template<class X> Vector<X> solve(const Matrix<X>& A, const Vector<X>& B);

// Compute the inverse using lower/upper triangular factorization
template<class X> Matrix<X> lu_inverse(const Matrix<X>& A);
template<class X> Matrix<X> lu_solve(const Matrix<X>& A, const Matrix<X>& B);
template<class X> Vector<X> lu_solve(const Matrix<X>& A, const Vector<X>& b);
// Compute the inverse using Gauss-Seidel iteration
template<class X> Matrix<X> gs_inverse(const Matrix<X>& A);
template<class X> Matrix<X> gs_solve(const Matrix<X>& A, const Matrix<X>& B);
template<class X> Vector<X> gs_solve(const Matrix<X>& A, const Vector<X>& b);

template<class X> Matrix<X> inverse(const PLUMatrix<X>& A);
template<class X> Matrix<X> solve(const PLUMatrix<X>& A, const Matrix<X>& B);
template<class X> Vector<X> solve(const PLUMatrix<X>& A, const Vector<X>& B);

template<class X> PLUMatrix<X> plu(const Matrix<X>& A);
template<class X> QRMatrix<X> qr(const Matrix<X>& A);

template<class X>
class DiagonalMatrix {
    Vector<X> _x;
  public:
    DiagonalMatrix(size_t n) : _x(n) { }
    DiagonalMatrix(const Vector<X>& x) : _x(x) { }
    size_t size() const { return _x.size(); }
    const X& operator[](size_t i) const { return _x[i]; }
    X& operator[](size_t i) { return _x[i]; }
    const X& get(size_t i) const { return _x[i]; }
    void set(size_t i, const X& x) { _x[i]=x; }
    const Vector<X>& diagonal() const { return _x; }
    template<class XX> Vector<XX> solve(const Vector<XX>& v) const {
        Vector<XX> result(_x.size()); for(uint i=0; i!=_x.size(); ++i) { result[i]=v[i]/_x[i]; } return result; }
};
template<class X> DiagonalMatrix<X> operator+(const DiagonalMatrix<X>& D1, const DiagonalMatrix<X>& D2) {
    return DiagonalMatrix<X>(D1.diagonal()+D2.diagonal());
}
template<class X> DiagonalMatrix<X> operator-(const DiagonalMatrix<X>& D1, const DiagonalMatrix<X>& D2) {
    return DiagonalMatrix<X>(D1.diagonal()-D2.diagonal());
}
template<class X> DiagonalMatrix<X> operator*(const DiagonalMatrix<X>& D1, const DiagonalMatrix<X>& D2) {
    DiagonalMatrix<X> r(D1.size()); for(size_t i=0; i!=r.size(); ++i) { r[i] = D1[i]*D2[i]; } return r;
}
template<class X> DiagonalMatrix<X> operator/(const DiagonalMatrix<X>& D1, const DiagonalMatrix<X>& D2) {
    DiagonalMatrix<X> r(D1.size()); for(size_t i=0; i!=r.size(); ++i) { r[i] = D1[i]/D2[i]; } return r;
}
template<class X, class XX> Vector<XX> operator*(const DiagonalMatrix<X>& D, const Vector<XX>& v) {
    Vector<XX> r(v.size()); for(size_t i=0; i!=r.size(); ++i) { r[i] = D[i]*v[i]; } return r;
}
template<class X, class XX> inline Vector<XX> operator*(const Vector<XX>& v1, const DiagonalMatrix<X>& D2) {
    Vector<XX> r(v1.size()); for(uint i=0; i!=r.size(); ++i) { r[i] = v1[i] * D2[i]; } return r;
}
template<class X, class XX> Vector<XX> operator/(const Vector<XX>& v, const DiagonalMatrix<X>& D) {
    Vector<XX> r(v.size()); for(size_t i=0; i!=r.size(); ++i) { r[i] = v[i]/D[i]; } return r;
}
template<class X, class XX> Matrix<XX> operator*(const DiagonalMatrix<X>& D, const Matrix<XX>& A) {
    ARIADNE_ASSERT_MSG(D.size()==A.row_size(),"D*A: D="<<D<<" A="<<A);
    Matrix<XX> R(A.row_size(),A.column_size());
    for(size_t i=0; i!=R.row_size(); ++i) { for(size_t j=0; j!=R.column_size(); ++j) { R[i][j] = D[i]*A[i][j]; } }
    return R;
}
template<class X, class XX> inline Matrix<XX> operator*(const Matrix<XX>& A, const DiagonalMatrix<X>& B) {
    Matrix<XX> R(A.row_size(),A.column_size());
    for(uint i=0; i!=A.row_size(); ++i) { for(uint j=0; j!=A.column_size(); ++j) { R[i][j]=A[i][j]*B.diagonal()[j]; } }
    return R;
}
template<class X> Matrix<X>& operator+=(Matrix<X>& A, const DiagonalMatrix<X>& D) {
    ARIADNE_ASSERT_MSG(D.size()==A.row_size() && D.size()==A.column_size(),"D*A: D="<<D<<" A="<<A);
    for(size_t i=0; i!=D.size(); ++i) { A[i][i] += D[i]; }
    return A;
}
template<class X> Matrix<X> operator+(const Matrix<X>& A, const DiagonalMatrix<X>& D) {
    Matrix<X> R(A); R+=D; return R;
}
template<class X> DiagonalMatrix<X> inverse(const DiagonalMatrix<X>& D) {
    DiagonalMatrix<X> r(D.size()); for(size_t i=0; i!=r.size(); ++i) { r[i] = 1/D[i]; } return r;
}
template<class X> std::ostream& operator<<(std::ostream& os, const DiagonalMatrix<X>& D) {
    return os << D.diagonal();
}

struct PivotMatrix {
    Array<size_t> _ary;
    PivotMatrix(size_t n=0u) : _ary(n) {
        for(uint i=0; i!=n; ++i) { _ary[i]=i; } }
    size_t size() const { return _ary.size(); }
    size_t const& operator[](size_t i) const { return _ary[i]; }
    size_t& operator[](size_t i) { return _ary[i]; }
    operator Matrix<Float> () const;
};
std::ostream& operator<<(std::ostream& os, const PivotMatrix& pv);

template<class X> struct PLUMatrix {
    PivotMatrix P; Matrix<X> LU;
};

template<class X> struct QRMatrix {
    Matrix<X> Q; Matrix<X> R;
};


Tuple< PivotMatrix, Matrix<Float>, Matrix<Float> > triangular_decomposition(const Matrix<Float>& A);

Vector<Float> row_norms(const Matrix<Float>& A);
Tuple< Matrix<Float>, PivotMatrix> triangular_factor(const Matrix<Float>& A);
Matrix<Float> triangular_multiplier(const Matrix<Float>& A);
Tuple< Matrix<Float>, Matrix<Float>, PivotMatrix > orthogonal_decomposition(const Matrix<Float>&, bool allow_pivoting=true);
Matrix<Float> normalise_rows(const Matrix<Float>& A);

Matrix<Float> midpoint(const Matrix<Interval>&);

inline std::ostream& operator<<(std::ostream& os, const Pretty< Interval >& pI) {
    Interval const& I(pI); return os << std::setprecision(5) << std::fixed << "{" << std::setw(8) << I.lower() << ":" << std::setw(8) << I.upper() << "}";
}

inline std::ostream& operator<<(std::ostream& os, const Pretty< Float >& px) {
    Float const& x(px); return os << std::setprecision(5) << std::fixed << std::setw(8) << x;
}

template<class X> inline std::ostream& operator<<(std::ostream& os, const Pretty< Vector<X> >& pv) {
    const Vector<X>& v(pv);
    os << "  ( ";
    for(size_t j=0; j!=v.size(); ++j) {
        if(j!=0) { os << ", "; }
        os << pretty(v[j]);
    }
    return os;
}

template<class X> inline std::ostream& operator<<(std::ostream& os, const Pretty< Matrix<X> >& pA) {
    const Matrix<X>& A(pA);
    for(size_t i=0; i!=A.row_size(); ++i) {
        os << "  ( ";
        for(size_t j=0; j!=A.column_size(); ++j) {
            if(j!=0) { os << ", "; }
            os << pretty(A[i][j]);
        }
        os << " )\n";
    }
    return os;
}

} // namespace Ariadne

#endif
