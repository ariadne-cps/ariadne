/***************************************************************************
 *            algebra/matrix.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file algebra/matrix.hpp
 *  \brief
 */

#ifndef ARIADNE_MATRIX_HPP
#define ARIADNE_MATRIX_HPP

#include <initializer_list>

#include "vector.hpp"
#include "covector.hpp"

namespace Ariadne {

template<class T> using InitializerList = InitializerList<T>;

/************ Matrix *********************************************************/

template<class X> class Vector;
template<class X> class Covector;

template<class X> class Matrix;
template<class M> struct MatrixRow;
template<class M> struct MatrixRows;
template<class M> struct MatrixColumn;
template<class M> struct MatrixTranspose;
template<class M1, class M2> struct MatrixMatrixProduct;
template<class M1, class V2> struct MatrixVectorProduct;

template<class M> using Row = MatrixRow<M>;
template<class M> using Column = MatrixColumn<M>;
template<class M> using Transpose = MatrixTranspose<M>;

class PivotMatrix;
template<class X> struct PLUMatrix;
template<class X> struct QRMatrix;

class SingularMatrixException { };


class DeclareMatrixOperations {
    template<class X1, class X2> friend Matrix<ProductType<Scalar<X1>,X2>> operator*(X1 const& s, Matrix<X2> const& A);
    template<class X1, class X2> friend Matrix<ProductType<X1,Scalar<X2>>> operator*(Matrix<X1> const& A, X2 const& s);
    template<class X1, class X2> friend Matrix<QuotientType<X1,Scalar<X2>>> operator/(Matrix<X1> const& A, X2 const& s);
    template<class X1, class X2> friend Matrix<InplaceProductType<X1,X2>>& operator*=(Matrix<X1>& A, X2 const& s);
    template<class X1, class X2> friend Matrix<InplaceQuotientType<X1,X2>>& operator/=(Matrix<X1>& A, X2 const& s);

    template<class X> friend Matrix<X> operator+(Matrix<X> A);
    template<class X> friend Matrix<X> operator-(Matrix<X> A);
    template<class X1, class X2> friend Matrix<SumType<X1,X2>> operator+(Matrix<X1> const& A1, Matrix<X2> const& A2);
    template<class X1, class X2> friend Matrix<DifferenceType<X1,X2>> operator-(Matrix<X1> const& A1, Matrix<X2> const& A2);
    template<class X1, class X2> friend Matrix<ArithmeticType<X1,X2>> operator*(Matrix<X1> const& A1, Matrix<X2> const& A2);
    template<class X1, class X2> friend Vector<ArithmeticType<X1,X2>> operator*(Matrix<X1> const& A1, Vector<X2> const& v2);
    template<class X1, class X2> friend Covector<ArithmeticType<X1,X2>> operator*(Covector<X1> const& u1, Matrix<X2> const& A2);
    template<class X1, class X2> friend Matrix<ProductType<X1,X2>> operator*(Vector<X1> const& v1, Covector<X2> const& u2);

    template<class X1, class X2> friend Matrix<ArithmeticType<X1,X2>> operator*(Matrix<X1> const& A1, Transpose<Matrix<X2>> const& A2);
    template<class X1, class X2> friend Matrix<ArithmeticType<X1,X2>> operator*(Transpose<Matrix<X1>> const& A1, Matrix<X2> const& A2);
    template<class X1, class X2> friend Vector<ArithmeticType<X1,X2>> operator*(Transpose<Matrix<X1>> const& A1, Vector<X2> const& v2);

    template<class X> friend decltype(abs(declval<X>()+declval<X>())) norm(Matrix<X> const& A);

    template<class X1, class X2> friend decltype(declval<X1>()==declval<X2>()) operator==(Matrix<X1> const& A1, Matrix<X2> const& A2);
};

template<class M> struct MatrixExpression : DeclareMatrixOperations { const M& operator()() const { return static_cast<const M&>(*this); } };
template<class M> struct MatrixContainer : public MatrixExpression<M> { };

/*

template<class X1, class X2=Void> struct MatrixOperations;

class DispatchMatrixOperations {
    template<class X1, class X2> friend Matrix<ProductType<Scalar<X1>,X2>> operator*(X1 const& s, Matrix<X2> const& A) { return MatrixOperations<X1,X2>::_mul(s,A); }
    template<class X1, class X2> friend Matrix<ProductType<X1,Scalar<X2>>> operator*(Matrix<X1> const& A, X2 const& s) { return MatrixOperations<X1,X2>::_mul(A,s); }
    template<class X1, class X2> friend Matrix<QuotientType<X1,Scalar<X2>>> operator/(Matrix<X1> const& A, X2 const& s) { return MatrixOperations<X1,X2>::_div(A,s); }

    template<class X> friend Matrix<X> operator+(Matrix<X> A) { return MatrixOperations<X>::_pos(A); }
    template<class X> friend Matrix<X> operator-(Matrix<X> A) { return MatrixOperations<X>::_neg(A); }
    template<class X1, class X2> friend Matrix<SumType<X1,X2>> operator+(Matrix<X1> const& A1, Matrix<X2> const& A2) { return MatrixOperations<X1,X2>::_add(A1,A2); }
    template<class X1, class X2> friend Matrix<DifferenceType<X1,X2>> operator-(Matrix<X1> const& A1, Matrix<X2> const& A2) { return MatrixOperations<X1,X2>::_sub(A1,A2); }
    template<class X1, class X2> friend Matrix<ArithmeticType<X1,X2>> operator*(Matrix<X1> const& A1, Matrix<X2> const& A2) { return MatrixOperations<X1,X2>::_mul(A1,A2); }
    template<class X1, class X2> friend Vector<ArithmeticType<X1,X2>> operator*(Matrix<X1> const& A1, Vector<X2> const& v2) { return MatrixOperations<X1,X2>::_mul(A1,v2); }
    template<class X1, class X2> friend Covector<ArithmeticType<X1,X2>> operator*(Covector<X1> const& u1, Matrix<X2> const& A2) { return MatrixOperations<X1,X2>::_mul(u1,A2); }
    template<class X1, class X2> friend Matrix<ProductType<X1,X2>> operator*(Vector<X1> const& v1, Covector<X2> const& u2) { return MatrixOperations<X1,X2>::_mul(v1,u2); }

    template<class X1, class X2> friend decltype(declval<X1>()==declval<X2>()) operator==(Matrix<X1> const& A1, Matrix<X2> const& A2);
};
*/

//! \ingroup LinearAlgebraModule
//! \brief Matrices over some type \a X.
template<class X> class Matrix
    : public MatrixContainer<Matrix<X>>
{
    X _zero;
    SizeType _rs;
    SizeType _cs;
    Array<X> _ary;
  public:
    typedef X ScalarType;
    typedef X ValueType;
  public:

    //@{
    //! \name Constructors

    //! Destructor
    ~Matrix();

    //! Default constructor makes a \f$0\times0\f$ matrix.
    Matrix();

    //! Construct a matrix with \a m rows and \a n columns with values uninitialised.
    //! The values should be initialised using placement new.
    Matrix(SizeType m, SizeType n, Uninitialised);

    //! Construct a matrix with \a r rows and \a c columns with values initialised to zero.
    Matrix(SizeType m, SizeType n);

    //! Construct a matrix with \a r rows and \a c columns with values initialised to \a x.
    Matrix(SizeType m, SizeType n, const X& z);

    //! Construct a matrix from parameters of \a X.
    template<class... PRS, EnableIf<IsConstructible<X,PRS...>> =dummy> explicit Matrix(SizeType m, SizeType n, PRS... prs) : Matrix(m,n,X(prs...)) { }

    //! Construct a matrix with \a r rows and \a c columns, with values initialised from the C-style Array beginning at \a ptr in row-major format. The value in the \a i<sup>th</sup> row and \a j<sup>th</sup> column of the resulting matrix is \a ptr[i*c+j].
    Matrix(SizeType m, SizeType n, const X* p);

    //! Construct a matrix using initializer lists.
    Matrix(InitializerList<InitializerList<X>> lst);

    //! Construct a matrix using initializer lists.
    template<class... PRS, EnableIf<IsConstructible<X,ExactDouble,PRS...>> =dummy> Matrix(InitializerList<InitializerList<double>> lst, PRS... pr);

    //! Construct a matrix as a linear map from functionals.
    Matrix(Vector<Covector<X>>);

    //@{
    //! \name Static constructors

    //! \brief The zero matrix with \a r rows and \a c columns.
    static Matrix<X> zero(SizeType m, SizeType n);
    //! \brief The itentity matrix with \a n rows and \a n columns.
    static Matrix<X> identity(SizeType n);
    //! Construct the identity matrix from parameters of \a X.
    template<class... PRS, EnableIf<IsConstructible<X,PRS...>> =dummy> static Matrix<X> identity(SizeType n, PRS... prs);
    //@}

    template<class M, EnableIf<And<IsMatrixExpression<M>,IsConvertible<typename M::ScalarType,X>>> =dummy>
        Matrix(const M& A);
    template<class M, EnableIf<And<IsMatrixExpression<M>,IsConvertible<typename M::ScalarType,X>>> =dummy>
        Matrix<X>& operator=(const M& A);
    template<class M, EnableIf<And<IsMatrixExpression<M>,IsConstructible<X,typename M::ScalarType>,Not<IsConvertible<typename M::ScalarType,X>>>> =dummy>
        explicit Matrix(const M& A);

    template<class M, class... PRS, EnableIf<And<IsMatrixExpression<M>,IsConstructible<X,typename M::ScalarType,PRS...>>> =dummy>
        explicit Matrix(const M& A, PRS... prs);


    template<class M> Matrix<X>& operator+=(const M& A);
    template<class M> Matrix<X>& operator-=(const M& A);

    //! \brief The number of rows of the matrix.
    SizeType row_size() const;
    //! \brief The number of columns of the matrix.
    SizeType column_size() const;
    //! \brief Resize to an \a m by \a n matrix.
    Void resize(SizeType m, SizeType n);
    //! \brief Get the value stored in the \a i<sup>th</sup> row and \a j<sup>th</sup> column.
    X& at(SizeType i, SizeType j);
    const X& at(SizeType i, SizeType j) const;
    //! \brief Get the value stored in the \a i<sup>th</sup> row and \a j<sup>th</sup> column.
    const X& get(SizeType i, SizeType j) const ;
    //! \brief Set the value stored in the \a i<sup>th</sup> row and \a j<sup>th</sup> column to \a x.
    Void set(SizeType i, SizeType j, const X& x);
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy> Void set(SizeType i, SizeType j, const Y& y);
    //! \brief A pointer to the first element of the data storage.
    const X* begin() const;

#ifdef DOXYGEN
    //! \brief C-style subscripting operator.
    X& operator[][](SizeType i, SizeType j);
    //! \brief C-style constant subscripting operator.
    const X& operator[][](SizeType i, SizeType j) const;
    //! \brief C-style range subscripting operator.
    X& operator[][](Range is, Range js);
    //! \brief C-style constant subscripting operator.
    const X& operator[][](Range i, Range js) const;
#else
    MatrixRow<const Matrix<X>> operator[](SizeType i) const;
    MatrixRow<Matrix<X>> operator[](SizeType i);
    MatrixRows<const Matrix<X>> operator[](Range is) const;
    MatrixRows<Matrix<X>> operator[](Range is);
#endif
    //! \brief The zero element of the field/algebra of the matrix.
    X zero_element() const;
  public:
    static Matrix<X> _mul(const Matrix<X>& A1, const Matrix<X>& A2);
  private:
    Void _check_data_access(SizeType i, SizeType j) const;
    OutputStream& write(OutputStream& os) const;
    InputStream& read(InputStream& is);

    template<class T> friend OutputStream& operator<<(OutputStream& os, Matrix<T>const& A);
    template<class T> friend InputStream& operator>>(InputStream& is, Matrix<T>& A);
};

template<class X> struct IsMatrix<Matrix<X>> : True { };
template<class X> struct IsMatrixExpression<Matrix<X>> : True { };

template<class M1, class M2, EnableIf<And<IsMatrixExpression<M1>,IsMatrixExpression<M2>>> =dummy>
auto operator==(M1 const& A1, M2 const& A2) -> decltype(declval<ScalarType<M1>>()==declval<ScalarType<M2>>());

template<class X> OutputStream& operator<<(OutputStream& os, Matrix<X> const& A);
template<class M, EnableIf<IsMatrixExpression<M>> =dummy> OutputStream& operator<<(OutputStream& os, M const& A);

/************ Combinatorial Matrices *********************************************************/

class PivotMatrix {
    Array<SizeType> _ary;
  public:
    PivotMatrix(SizeType n=0u) : _ary(n) {
        for(SizeType i=0; i!=n; ++i) { _ary[i]=i; } }
    SizeType size() const { return _ary.size(); }
    SizeType const& operator[](SizeType i) const { return _ary[i]; }
    SizeType& operator[](SizeType i) { return _ary[i]; }
    template<class X> operator Matrix<X> () const;
};
OutputStream& operator<<(OutputStream& os, const PivotMatrix& pv);
template<class X> Vector<X> operator*(PivotMatrix, Vector<X>);
template<class X> Matrix<X> operator*(PivotMatrix, Matrix<X>);

template<class X> struct PLUMatrix {
    PivotMatrix P; Matrix<X> LU;
};

template<class X> struct QRMatrix {
    Matrix<X> Q; Matrix<X> R;
};


/************ Matrix expressions *********************************************************/


// NOTE: Can't just use project(MatrixRow,Range) as the MatrixRow lifetime is too short
template<class M> struct MatrixRowColumns
    : CovectorExpression<MatrixRowColumns<M>>
{
  public:
    M& _A; SizeType _i; Range _js;
  public:
    typedef typename M::ScalarType ScalarType;
    MatrixRowColumns(M& A, SizeType i, Range js) : _A(A), _i(i), _js(js) { }
    SizeType size() const { return _js.size(); }
    auto zero_element() const -> decltype(_A.zero_element()) { return _A.zero_element(); }
    decltype(auto) operator[](SizeType j) { return _A[_i][_js[j]]; }
    decltype(auto) operator[](SizeType j) const { return _A[_i][_js[j]]; }
    MatrixRowColumns<M>& operator=(Covector<ScalarType> const& u) {
        for(SizeType j=0; j!=u.size(); ++j) { _A.set(_i,_js[j],u[j]); } return *this; }
};
template<class M> struct IsCovectorExpression<MatrixRowColumns<M>> : True { };

template<class M> struct MatrixRow
    : public CovectorExpression< MatrixRow<M> >
{
    M& _A; SizeType _i;
  public:
    typedef typename M::ScalarType ScalarType;
    MatrixRow(M& A, SizeType i) : _A(A), _i(i) { }
    SizeType size() const { return _A.column_size(); }
    auto zero_element() const -> decltype(_A.zero_element()) { return _A.zero_element(); }
    auto operator[](SizeType j) const -> decltype(_A.at(_i,j)) { return _A.at(_i,j); }
    auto operator[](SizeType j) -> decltype(_A.at(_i,j)) { return _A.at(_i,j); }
    MatrixRowColumns<M> operator[](Range js) { return MatrixRowColumns<M>(_A,_i,js); }
    MatrixRow<M>& operator=(Covector<ScalarType> const& u) { for(SizeType j=0; j!=u.size(); ++j) { _A.set(_i,j,u[j]); } return *this; }
};
template<class M> struct IsCovectorExpression<MatrixRow<M>> : True { };

template<class M> struct MatrixColumn
    : public VectorExpression< MatrixColumn<M> >
{
     M& _A; SizeType _j;
  public:
    typedef typename M::ScalarType ScalarType;
    MatrixColumn(M& A, SizeType j) : _A(A), _j(j) { }
    SizeType size() const { return _A.row_size(); }
    auto zero_element() const -> decltype(_A.zero_element()) { return _A.zero_element(); }
    auto operator[](SizeType i) const -> decltype(_A.at(i,_j)) { return _A[i][_j]; }
    auto operator[](SizeType i) -> decltype(_A.at(i,_j)) { return _A.at(i,_j); }
    operator Vector<ScalarType> () const {
        Vector<ScalarType> r(size(),zero_element()); for(SizeType i=0; i!=size(); ++i) { r[i]=_A.at(i,_j); } return r; }
};
template<class M> struct IsVectorExpression<MatrixColumn<M>> : True { };

template<class M> struct MatrixTranspose {
    M const& _AT;
  public:
    typedef typename M::ScalarType ScalarType;
    MatrixTranspose(M const& AT) : _AT(AT) { }
    SizeType row_size() const { return _AT.column_size(); }
    SizeType column_size() const { return _AT.row_size(); }
    ScalarType zero_element() const { return _AT.zero_element(); }
    ScalarType const& at(SizeType i, SizeType j) const { return _AT.at(j,i); }
    ScalarType const& get(SizeType i, SizeType j) const { return _AT.at(j,i); }
};
template<class M> struct IsMatrixExpression<MatrixTranspose<M>> : True { };

template<class M1, class M2> struct MatrixMatrixProduct {
    M1 const& _A1; M2 const& _A2;
  public:
    typedef ArithmeticType<typename M1::ScalarType,typename M2::ScalarType> ScalarType;
    MatrixMatrixProduct(M1 const& A1, M2 const& A2) : _A1(A1), _A2(A2) { }
    SizeType row_size() const { return _A1.row_size(); }
    SizeType column_size() const { return _A2.column_size(); }
    ScalarType zero_element() const { return _A1.zero_element()*_A2.zero_element(); }
    ScalarType at(SizeType i, SizeType j) const { ScalarType r=this->zero_element();
        for(SizeType k=0; k!=_A1.row_size(); ++k) { r+=_A1.at(i,k)*_A2.at(k,j); } return std::move(r); }
};
template<class M1,class M2> struct IsMatrixExpression<MatrixMatrixProduct<M1,M2>> : True { };

template<class M1, class V2> struct MatrixVectorProduct {
   typedef ArithmeticType<typename M1::ScalarType,typename V2::ScalarType> ScalarType;
  public:
    M1 const& _A1; V2 const& _v2;
    MatrixVectorProduct(M1 const& A1, V2 const& v2) : _A1(A1), _v2(v2) { }
    SizeType size() const { return _A1.row_size(); }
    ScalarType zero_element() const { return _A1.zero_element()*_v2.zero_element(); }
    ScalarType operator[](SizeType i) const { ScalarType r=this->zero_element();
        for(SizeType j=0; j!=_v2.size(); ++j) { r+=_A1.at(i,j)*_v2.at(j); } return std::move(r); }
};
template<class M1,class V2> struct IsVectorExpression<MatrixVectorProduct<M1,V2>> : True { };

template<class M1, class X2> struct MatrixScalarQuotient {
    const M1& _a1; const X2& _x2;
  public:
    typedef QuotientType<typename M1::ScalarType,X2> ScalarType;
    MatrixScalarQuotient(const M1& a1, const X2& x2) : _a1(a1), _x2(x2) { }
    SizeType row_size() const { return _a1.row_size(); }
    SizeType column_size() const { return _a1.column_size(); }
    auto at(SizeType i, SizeType j) const -> decltype(_a1.at(i,j)/_x2) { return _a1.at(i,j)/_x2; }
};
template<class M1, class X2> struct IsMatrixExpression<MatrixScalarQuotient<M1,X2>> : True { };

template<class M> struct MatrixRange
    : public MatrixContainer< MatrixRange<M> >
{
    M& _A; Range _rng1; Range _rng2;
  public:
    typedef typename M::ScalarType ScalarType;
    MatrixRange(M& A, Range rng1, Range rng2) : _A(A), _rng1(rng1), _rng2(rng2) { }
    SizeType row_size() const { return _rng1.size(); }
    SizeType column_size() const { return _rng2.size(); }
    ScalarType zero_element() const { return _A.zero_element(); }
    ScalarType get(SizeType i, SizeType j) const { return _A.get(i+_rng1.start(),j+_rng2.start()); }
    decltype(auto) operator[](SizeType i) { return _A[_rng1[i]][_rng2]; }
    Void set(SizeType i, SizeType j, ScalarType const& x) const { _A.set(i+_rng1.start(),j+_rng2.start(),x); }
    template<class ME> MatrixRange<M>& operator=(const MatrixExpression<ME>& Ae) {
        ARIADNE_PRECONDITION(this->row_size()==Ae().row_size() && this->column_size()==Ae().column_size());
        for(SizeType i=0; i!=this->row_size(); ++i) { for(SizeType j=0; j!=this->column_size(); ++j) { this->set(i,j,Ae().get(i,j)); } }
        return *this; }
};
template<class M> struct IsMatrixExpression<MatrixRange<M>> : True { };

template<class M> struct MatrixRowsColumn
    : VectorExpression<MatrixRowsColumn<M>>
{
  public:
    M& _A; Range _is; SizeType _j;
  public:
    typedef typename M::ScalarType ScalarType;
    MatrixRowsColumn(M& A, Range is, SizeType j) : _A(A), _is(is), _j(j) { }
    SizeType size() const { return _is.size(); }
    auto zero_element() const -> decltype(_A.zero_element()) { return _A.zero_element(); }
    decltype(auto) operator[](SizeType i) { return _A[_is[i]][_j]; }
    decltype(auto) operator[](SizeType i) const { return _A[_is[i]][_j]; }
    MatrixRowsColumn<M>& operator=(Vector<ScalarType> const& v) {
        for(SizeType i=0; i!=v.size(); ++i) { _A.set(_is[i],_j,v[i]); } return *this; }
};
template<class M> struct IsVectorExpression<MatrixRowsColumn<M>> : True { };

template<class M> struct MatrixRows {
  public:
    M& _A; Range _is;
  public:
    typedef typename M::ScalarType ScalarType;
    MatrixRows(M& A, Range is) : _A(A), _is(is) { }
    SizeType row_size() const { return _is.size(); }
    SizeType column_size() const { return _A.column_size(); }
    auto zero_element() const -> decltype(_A.zero_element()) { return _A.zero_element(); }
    // Note: Can't project to VectorExpression below since j in MatrixColumn(_A,j) would be destroyed before use
    MatrixRowsColumn<M> operator[](SizeType j) { return MatrixRowsColumn<M>(_A,_is,j); }
    MatrixRange<M> operator[](Range js) { return MatrixRange<M>(_A,_is,js); }
    MatrixRows<M>& operator=(Matrix<ScalarType> const& B) {
        for(SizeType i=0; i!=B.row_size(); ++i) { SizeType p_i=_is[i]; for(SizeType j=0; j!=B.column_size(); ++j) { _A.set(p_i,j,B.get(i,j)); } } return *this; }
};
template<class M> struct IsMatrixExpression<MatrixRows<M>> : True { };

template<class X> inline MatrixRows<const Matrix<X>> Matrix<X>::operator[](Range is) const {
    return MatrixRows<const Matrix<X>>(*this,is);
}
template<class X> inline MatrixRows<Matrix<X>> Matrix<X>::operator[](Range is) {
    return MatrixRows<Matrix<X>>(*this,is);
}

template<class M1, class M2, EnableIf<And<IsMatrixExpression<M1>,IsMatrixExpression<M2>>>> inline
auto operator==(M1 const& A1, M2 const& A2) -> decltype(declval<ScalarType<M1>>()==declval<ScalarType<M2>>()) {
    typedef ScalarType<M1> X1;typedef ScalarType<M2> X2; return Matrix<X1>(A1)==Matrix<X2>(A2); }

/* Dispatching Matrix expression template operators


template<class M, EnableIf<IsMatrix<M>> =dummy> inline
MatrixScalarQuotient<M,typename M::ScalarType> operator/(const M& A1, typename M::ScalarType const& x2) {
    return MatrixScalarQuotient<M,typename M::ScalarType>(A1,x2); }

template<class X1, class X2> inline MatrixVectorProduct<X1,X2> operator*(Matrix<X1> const& A1, Vector<X2> const& v2) {
    return MatrixVectorProduct<X1,X2>(A1,v2);
}

template<class X1, class V2, EnableIf<IsVectorExpression<V2>> =dummy> inline MatrixVectorProduct<X1,ScalarType<V2>> operator*(Matrix<X1> const& A1, V2 const& v2) {
    typedef typename V2::ScalarType X2; return MatrixVectorProduct<X1,X2>(A1,Vector<X2>(v2));
}

*/

/************ Matrix inline *********************************************************/

template<class X> inline Matrix<X>::~Matrix()
{
}

template<class X> inline Matrix<X>::Matrix()
    : _zero(), _rs(0), _cs(0), _ary() {
}

template<class X> Matrix<X>::Matrix(SizeType m, SizeType n)
    : _zero(), _rs(m), _cs(n), _ary(m*n) {
}

template<class X> Matrix<X>::Matrix(SizeType m, SizeType n, const X& x)
    : _zero(create_zero(x)), _rs(m), _cs(n), _ary(m*n,x) {
}

template<class X> inline Void Matrix<X>::_check_data_access(SizeType i, SizeType j) const {
    ARIADNE_PRECONDITION_MSG(i<this->row_size()&&j<this->column_size(),"A="<<*this<<" i="<<i<<" j="<<j);
}

template<class X> inline MatrixRow<const Matrix<X>> Matrix<X>::operator[](SizeType i) const {
    return MatrixRow<const Matrix<X>>(*this,i);
}

template<class X> inline MatrixRow<Matrix<X>> Matrix<X>::operator[](SizeType i) {
    return MatrixRow<Matrix<X>>(*this,i);
}

template<class X> inline SizeType Matrix<X>::row_size() const {
    return this->_rs;
}

template<class X> inline SizeType Matrix<X>::column_size() const {
    return this->_cs;
}

template<class X> inline const X& Matrix<X>::at(SizeType i, SizeType j) const {
    this->_check_data_access(i,j); return this->_ary[i*this->_cs+j];
}

template<class X> inline X& Matrix<X>::at(SizeType i, SizeType j) {
    this->_check_data_access(i,j); return this->_ary[i*this->_cs+j];
}

template<class X> inline Void Matrix<X>::set(SizeType i, SizeType j, const X& c) {
    this->_ary[i*this->_cs+j]=c;
}

template<class X> inline const X& Matrix<X>::get(SizeType i, SizeType j) const {
    return this->_ary[i*this->_cs+j];
}

template<class X> inline const X* Matrix<X>::begin() const {
    return this->_ary.begin();
}

template<class X> inline OutputStream& operator<<(OutputStream& os, Matrix<X> const& A) {
    A.write(os); return os;
}

template<class X> inline InputStream& operator>>(InputStream& is, Matrix<X>& A) {
    A.read(is); return is;
}

template<class X> OutputStream& Matrix<X>::write(OutputStream& os) const {
    const Matrix<X>& A=*this;
    if(A.row_size()==0 || A.column_size()==0) { os << "["; }
    for(SizeType i=0; i!=A.row_size(); ++i) {
        for(SizeType j=0; j!=A.column_size(); ++j) {
            os << (j==0 ? (i==0 ? "[" : "; ") : ",") << A.at(i,j); } }
    return os << "]";
}



template<class M> inline MatrixRange<const M> project(const MatrixExpression<M>& Ae, Range rw_rng, Range cl_rng) {
    return MatrixRange<const M>(Ae(),rw_rng,cl_rng);
}

template<class X> inline MatrixRange<Matrix<X>> project(Matrix<X>& A, Range rw_rng, Range cl_rng) {
    return MatrixRange<Matrix<X>>(A,rw_rng,cl_rng);
}

template<class M> inline MatrixColumn<const M> column(const M& A, SizeType j) {
    return MatrixColumn<const M>(A,j);
}

template<class M> inline MatrixRow<const M> row(const M& A, SizeType j) {
    return MatrixRow<const M>(A,j);
}



#ifdef ARIADNE_UNDEF
template<class X> Matrix<X>& Matrix<X>::operator=(const MatrixMatrixProduct<X,X>& A1mulA2) {
    Matrix<X> const& A1=A1mulA2._a1;
    Matrix<X> const& A2=A1mulA2._a2;
    if(this==&A1 || this==&A2) {
        Matrix<X> A0(A1.row_size(),A2.column_size());
        _mul_assign(A0,A1,A2);
        *this = std::move(A0);
    } else {
        Matrix<X>& A0=*this;
        A0.resize(A1.row_size(),A2.column_size());
        _mul_assign(A0,A1,A2);
    }
    return *this;
}
#endif

template<class X> template<class M, EnableIf<And<IsMatrixExpression<M>,IsConvertible<typename M::ScalarType,X>>>> Matrix<X>::Matrix(const M& A)
    : Matrix(A.row_size(),A.column_size(),A.zero_element())
{
    this->operator=(A);
}

template<class X> template<class M, EnableIf<And<IsMatrixExpression<M>,IsConstructible<X,typename M::ScalarType>,Not<IsConvertible<typename M::ScalarType,X>>>>>
Matrix<X>::Matrix(const M& A) : Matrix(A.row_size(),A.column_size(),X(A.zero_element()))
{
    for(SizeType i=0; i!=this->row_size(); ++i) {
        for(SizeType j=0; j!=this->column_size(); ++j) {
            this->at(i,j)=X(A.get(i,j));
        }
    }
}

template<class X> template<class M, class... PRS, EnableIf<And<IsMatrixExpression<M>,IsConstructible<X,typename M::ScalarType,PRS...>>>>
Matrix<X>::Matrix(const M& A, PRS... prs) : Matrix(A.row_size(),A.column_size(),X(A.zero_element(),prs...)) {
    for(SizeType i=0; i!=this->row_size(); ++i) {
        for(SizeType j=0; j!=this->column_size(); ++j) {
            this->at(i,j)=X(A.get(i,j),prs...);
        }
    }
}


template<class X> template<class M, EnableIf<And<IsMatrixExpression<M>,IsConvertible<typename M::ScalarType,X>>>> Matrix<X>& Matrix<X>::operator=(const M& A) {
    this->resize(A.row_size(),A.column_size());
    for(SizeType i=0; i!=this->row_size(); ++i) {
        for(SizeType j=0; j!=this->column_size(); ++j) {
            this->at(i,j)=A.get(i,j);
        }
    }
    return *this;
}

template<class M, EnableIf<IsMatrixExpression<M>>> inline OutputStream& operator<<(OutputStream& os, const M& A) {
    return os << Matrix<ScalarType<M>>(A); }


template<class X> template<class M> Matrix<X>& Matrix<X>::operator+=(const M& A) {
    for(SizeType i=0; i!=this->row_size(); ++i) {
        for(SizeType j=0; j!=this->column_size(); ++j) {
            this->at(i,j)+=A.get(i,j);
        }
    }
    return *this;
}

template<class X> template<class M> Matrix<X>& Matrix<X>::operator-=(const M& A) {
    for(SizeType i=0; i!=this->row_size(); ++i) {
        for(SizeType j=0; j!=this->column_size(); ++j) {
            this->at(i,j)-=A.get(i,j);
        }
    }
    return *this;
}

template<class X0, class X1, class X2> Void _mul_assign(Matrix<X0>& A0, Matrix<X1> const& A1, Matrix<X2> const& A2) {
    for(SizeType i=0; i!=A0.row_size(); ++i) {
        for(SizeType j=0; j!=A0.column_size(); ++j) {
            A0.at(i,j)=0;
            for(SizeType k=0; k!=A1.column_size(); ++k) {
                A0.at(i,j)+=A1.at(i,k)*A2.at(k,j);
            }
        }
    }
}

struct ProvideMatrixOperations {

    template<class X> friend inline Matrix<X> operator+(Matrix<X> A) {
        return std::move(A);
    }

    template<class X> friend inline Matrix<X> operator-(Matrix<X> A) {
        for(SizeType i=0; i!=A.row_size(); ++i) {
            for(SizeType j=0; j!=A.column_size(); ++j) {
                A[i][j]=-A[i][j];
            }
        }
        return std::move(A);
    }

    template<class X1, class X2> friend inline Matrix<SumType<X1,X2>> operator+(Matrix<X1> const& A1, Matrix<X2> const& A2) {
        ARIADNE_PRECONDITION(A1.row_size()==A2.row_size());
        ARIADNE_PRECONDITION(A1.column_size()==A2.column_size());
        Matrix<SumType<X1,X2>> R(A1.row_size(),A1.column_size(),A1.zero_element()+A2.zero_element());
        for(SizeType i=0; i!=A1.row_size(); ++i) {
            for(SizeType j=0; j!=A1.column_size(); ++j) {
                R[i][j]=A1[i][j]+A2[i][j];
            }
        }
        return std::move(R);
    }

    template<class X1, class X2> friend inline Matrix<DifferenceType<X1,X2>> operator-(Matrix<X1> const& A1, Matrix<X2> const& A2) {
        ARIADNE_PRECONDITION(A1.row_size()==A2.row_size());
        ARIADNE_PRECONDITION(A1.column_size()==A2.column_size());
        Matrix<DifferenceType<X1,X2>> R(A1.row_size(),A1.column_size(),A1.zero_element()-A2.zero_element());
        for(SizeType i=0; i!=A1.row_size(); ++i) {
            for(SizeType j=0; j!=A1.column_size(); ++j) {
                R[i][j]=A1[i][j]-A2[i][j];
            }
        }
        return std::move(R);
    }

    template<class X1, class X2> friend inline Matrix<ProductType<Scalar<X1>,X2>> operator*(X1 const& s1, Matrix<X2> const& A2) {
        Matrix<ProductType<X1,X2>> R(A2.row_size(),A2.row_size(),s1*A2.zero_element());
        for(SizeType i=0; i!=A2.row_size(); ++i) {
            for(SizeType j=0; j!=A2.column_size(); ++j) {
                R[i][j]=s1*A2[i][j];
            }
        }
        return std::move(R);
    }

    template<class X1, class X2> friend inline Matrix<ProductType<X1,Scalar<X2>>> operator*(Matrix<X1> const& A1, X2 const& s2) {
        Matrix<ProductType<X1,X2>> R(A1.row_size(),A1.row_size(),A1.zero_element()*s2);
        for(SizeType i=0; i!=A1.row_size(); ++i) {
            for(SizeType j=0; j!=A1.column_size(); ++j) {
                R[i][j]=A1[i][j]*s2;
            }
        }
        return std::move(R);
    }

    template<class X1, class X2> friend inline Matrix<QuotientType<X1,Scalar<X2>>> operator/(Matrix<X1> const& A1, X2 const& s2) {
        Matrix<QuotientType<X1,X2>> R(A1.row_size(),A1.row_size(),A1.zero_element()/s2);
        for(SizeType i=0; i!=A1.row_size(); ++i) {
            for(SizeType j=0; j!=A1.column_size(); ++j) {
                R[i][j]=A1[i][j]/s2;
            }
        }
        return std::move(R);
    }

    template<class X1, class X2> friend inline Matrix<InplaceProductType<X1,X2>>& operator*=(Matrix<X1>& A1, X2 const& s2) {
        for(SizeType i=0; i!=A1.row_size(); ++i) {
            for(SizeType j=0; j!=A1.column_size(); ++j) {
                A1[i][j]*=s2;
            }
        }
        return A1;
    }

    template<class X1, class X2> friend inline Matrix<InplaceQuotientType<X1,X2>>& operator/=(Matrix<X1>& A1, X2 const& s2) {
        for(SizeType i=0; i!=A1.row_size(); ++i) {
            for(SizeType j=0; j!=A1.column_size(); ++j) {
                A1[i][j]/=s2;
            }
        }
        return A1;
    }

    //template<class X1, class X2> friend inline MatrixMatrixProduct<X1,X2> operator*(Matrix<X1> const& A1, Matrix<X2> const& A2) {
    //    return MatrixMatrixProduct<X1,X2>(A1,A2); }

    template<class X1, class X2> friend Matrix<ArithmeticType<X1,X2>> operator*(Matrix<X1> const& A1, Matrix<X2> const& A2) {
        typedef ArithmeticType<X1,X2> X0;
        ARIADNE_PRECONDITION(A1.column_size()==A2.row_size());
        Matrix<X0> A0(A1.row_size(), A2.column_size(),A1.zero_element()*A2.zero_element());
        for(SizeType i=0; i!=A0.row_size(); ++i) {
            for(SizeType j=0; j!=A0.column_size(); ++j) {
                for(SizeType k=0; k!=A1.column_size(); ++k) {
                    A0.at(i,j)+=A1.at(i,k)*A2.at(k,j);
                }
            }
        }
        return std::move(A0);
    }

    template<class X1, class X2> friend Vector<ArithmeticType<X1,X2>> operator*(Matrix<X1> const& A1, Vector<X2> const& v2) {
        typedef ArithmeticType<X1,X2> X0;
        ARIADNE_PRECONDITION(A1.column_size()==v2.size());
        Vector<X0> v0(A1.row_size(),A1.zero_element()*v2.zero_element());
        for(SizeType i=0; i!=v0.size(); ++i) {
            for(SizeType j=0; j!=v2.size(); ++j) {
                v0.at(i)+=A1.at(i,j)*v2.at(j);
            }
        }
        return std::move(v0);
    }

    template<class X1, class X2> friend Covector<ArithmeticType<X1,X2>> operator*(Covector<X1> const& u1, Matrix<X2> const& A2) {
        typedef ArithmeticType<X1,X2> X0;
        ARIADNE_PRECONDITION(u1.size()==A2.row_size());
        Covector<X0> u0(A2.column_size(),u1.zero_element()*A2.zero_element());
        for(SizeType j=0; j!=u0.size(); ++j) {
            for(SizeType i=0; i!=u1.size(); ++i) {
                u0.at(j)+=u1.at(i)*A2.at(i,j);
            }
        }
        return std::move(u0);
    }

    template<class X1, class X2> friend Matrix<ArithmeticType<X1,X2>> operator*(Matrix<X1> const& A1, Transpose<Matrix<X2>> const& A2) {
        typedef ArithmeticType<X1,X2> X0;
        ARIADNE_PRECONDITION(A1.column_size()==A2.row_size());
        Matrix<X0> A0(A1.row_size(), A2.column_size(),A1.zero_element()*A2.zero_element());
        for(SizeType i=0; i!=A0.row_size(); ++i) {
            for(SizeType j=0; j!=A0.column_size(); ++j) {
                for(SizeType k=0; k!=A1.column_size(); ++k) {
                    A0.at(i,j)+=A1.at(i,k)*A2.at(k,j);
                }
            }
        }
        return std::move(A0);
    }

    template<class X1, class X2> friend Matrix<ArithmeticType<X1,X2>> operator*(Transpose<Matrix<X1>> const& A1, Matrix<X2> const& A2) {
        Matrix<X1> const& A1T=A1._AT;
        ARIADNE_PRECONDITION(A1T.row_size()==A2.row_size());
        Matrix<ArithmeticType<X1,X2>> A0(A1T.column_size(),A2.column_size(),A1.zero_element()*A2.zero_element());
        for(SizeType i=0; i!=A1T.column_size(); ++i) {
            for(SizeType j=0; j!=A2.column_size(); ++j) {
                for(SizeType k=0; k!=A1T.row_size(); ++k) {
                    A0.at(i,j)+=A1T.at(k,i)*A2.at(k,j);
                }
            }
        }
        return std::move(A0);
    }

    template<class X1, class X2> friend Vector<ArithmeticType<X1,X2>> operator*(MatrixTranspose<Matrix<X1>> const& A1, Vector<X2> const& v2) {
        Matrix<X1> const& A1T=A1._AT;
        ARIADNE_PRECONDITION(A1T.row_size()==v2.size());
        Vector<ArithmeticType<X1,X2>> v0(A1T.column_size(),A1.zero_element()*v2.zero_element());
        for(SizeType i=0; i!=v0.size(); ++i) {
            for(SizeType j=0; j!=v2.size(); ++j) {
                v0.at(i)+=A1T.at(j,i)*v2.at(j);
            }
        }
        return std::move(v0);
    }

    template<class X> friend inline Matrix<X>& operator*=(Matrix<X>& A, typename Matrix<X>::ScalarType const& s) {
        for(SizeType i=0; i!=A.row_size(); ++i) {
            for(SizeType j=0; j!=A.column_size(); ++j) {
                A.at(i,j)*=s;
            }
        }
        return A;
    }

    template<class X1, class X2> friend inline
    auto operator==(Matrix<X1> const& A1, Matrix<X2> const& A2) -> decltype(declval<X1>()==declval<X2>()) {
        if(A1.row_size()!=A2.row_size() || A1.column_size()!=A2.column_size()) { return false; }
        decltype(declval<X1>()==declval<X2>()) r=true;
        for(SizeType i=0; i!=A1.row_size(); ++i) {
            for(SizeType j=0; j!=A1.column_size(); ++j) {
                r=r and (A1.at(i,j)==A2.at(i,j));
            }
        }
        return r;
    }

    template<class X> friend decltype(abs(declval<X>()+declval<X>())) norm(Matrix<X> const& A) {
        typedef decltype(abs(declval<X>()+declval<X>())) R;
        R r=abs(A.zero_element());
        for(SizeType i=0; i!=A.row_size(); ++i) {
            R s=abs(A.zero_element());
            for(SizeType j=0; j!=A.column_size(); ++j) {
                s+=abs(A.at(i,j));
            }
            r=max(r,s);
        }
        return r;
    }

};

template<class X> template<class... PRS, EnableIf<IsConstructible<X,ExactDouble,PRS...>>>
Matrix<X>::Matrix(InitializerList<InitializerList<double>> lst, PRS... prs) : _rs(lst.size()), _cs(lst.begin()->size()), _ary(_rs*_cs,X(prs...)) {
    typename InitializerList<InitializerList<double>>::const_iterator row_iter=lst.begin();
    for(SizeType i=0; i!=this->row_size(); ++i, ++row_iter) {
        ARIADNE_PRECONDITION(row_iter->size()==this->column_size());
        typename InitializerList<double>::const_iterator col_iter=row_iter->begin();
        for(SizeType j=0; j!=this->column_size(); ++j, ++col_iter) {
            this->at(i,j)=X(ExactDouble(*col_iter),prs...);
        }
    }
}

template<class X> template<class... PRS, EnableIf<IsConstructible<X,PRS...>>> auto
Matrix<X>::identity(SizeType n, PRS... prs) -> Matrix<X> {
    Matrix<X> I(n,n,X(prs...));
    for(SizeType i=0; i!=n; ++i) {
        I.at(i,i)=1u;
    }
    return I;
}

template<class X> template<class Y, EnableIf<IsAssignable<X,Y>>>
Void Matrix<X>::set(SizeType i, SizeType j, Y const& y) {
    this->at(i,j)=y;
}

template<class X> Matrix<MidpointType<X>> midpoint(const Matrix<X>&);
template<class X> Matrix<SingletonType<X>> cast_singleton(const Matrix<X>&);


// Construct transpose
template<class X> inline Transpose<Matrix<X>> transpose(const Matrix<X>& A) { return Transpose<Matrix<X>>(A); }

// Invert matrices and solve linear systems
template<class X> Matrix<ArithmeticType<X>> inverse(Matrix<X> const& A);
template<class X1, class X2> Vector<ArithmeticType<X1,X2>> solve(Matrix<X1> const& A, Vector<X2> const& b);
template<class X1, class X2> Matrix<ArithmeticType<X1,X2>> solve(Matrix<X1> const& A, Matrix<X2> const& B);

// Compute the inverse using lower/upper triangular factorization
template<class X> Matrix<X> lu_inverse(const Matrix<X>& A);
template<class X> Matrix<X> lu_solve(const Matrix<X>& A, const Matrix<X>& B);
template<class X> Vector<X> lu_solve(const Matrix<X>& A, const Vector<X>& b);
// Compute the inverse using Gauss-Seidel iteration
template<class X> Matrix<X> gs_inverse(const Matrix<X>& A);
template<class X> Matrix<X> gs_solve(const Matrix<X>& A, const Matrix<X>& B);
template<class X> Vector<X> gs_solve(const Matrix<X>& A, const Vector<X>& b);

// Use Gauss-Seidel iteration hotstarted by iX
template<class X> Vector<X> gs_solve(Matrix<X> const& A, Vector<X> const& b, Vector<X> iX);
template<class X> Void gs_step(Matrix<X> const& A, Vector<X> const& b, Vector<X>& iX);

// Compute an LU decomposition
template<class X> Tuple< PivotMatrix, Matrix<X>, Matrix<X> > triangular_decomposition(const Matrix<X>& A);
// Compute an QR decomposition
template<class X> Tuple< Matrix<X>, Matrix<X> > orthogonal_decomposition(const Matrix<X>& A);

template<class X> using RowNormType = decltype(abs(declval<X>())+abs(declval<X>()));

template<class X> Vector<RowNormType<X>> row_norms(const Matrix<X>& A);
template<class X> Matrix<X> normalise_rows(const Matrix<X>& A);

Tuple< Matrix<FloatDPApproximation>, PivotMatrix> triangular_factor(const Matrix<FloatDPApproximation>& A);
Matrix<FloatDPApproximation> triangular_multiplier(const Matrix<FloatDPApproximation>& A);


} // namespace Ariadne

#endif
