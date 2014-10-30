/***************************************************************************
 *            algebra/matrix.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file algebra/matrix.h
 *  \brief
 */



#ifndef ARIADNE_MATRIX_H
#define ARIADNE_MATRIX_H

#include <initializer_list>

#include "vector.h"

namespace Ariadne {

template<class T> using InitializerList = std::initializer_list<T>;

/************ Matrix *********************************************************/

template<class X> class Vector;
template<class X> class Covector;

template<class X> class Matrix;
template<class M> class MatrixRow;
template<class M> class MatrixColumn;
template<class X1, class X2> struct MatrixMatrixProduct;
template<class X1, class X2> struct MatrixVectorProduct;

class PivotMatrix;
template<class X> class PLUMatrix;

//! \ingroup LinearAlgebraSubModule
//! \brief Matrices over some type \a X.
template<class X> class Matrix {
    X _zero;
    SizeType _rs;
    SizeType _cs;
    Array<X> _ary;
  public:
    typedef X ScalarType;
    ~Matrix();
    Matrix();
    Matrix(SizeType m, SizeType n);
    Matrix(SizeType m, SizeType n, const X& x);
    Matrix(InitializerList<InitializerList<X>> lst);
    Matrix(Vector<Covector<X>>);
    static Matrix<X> zero(SizeType m, SizeType n);
    static Matrix<X> identity(SizeType n);

    template<class M, EnableIf<And<IsMatrixExpression<M>,IsConvertible<typename M::ScalarType,X>>> =dummy> Matrix(const M& A);
    template<class M, EnableIf<And<IsMatrixExpression<M>,IsConvertible<typename M::ScalarType,X>>> =dummy> Matrix<X>& operator=(const M& A);
    template<class M> Matrix<X>& operator+=(const M& A);
    SizeType row_size() const;
    SizeType column_size() const;
    Void resize(SizeType m, SizeType n);
    MatrixRow<const Matrix<X>> operator[](SizeType i) const;
    MatrixRow<Matrix<X>> operator[](SizeType i);
    const X& at(SizeType i, SizeType j) const;
    X& at(SizeType i, SizeType j);
    Void set(SizeType i, SizeType j, const X& c);
    X zero_element() const;
    OutputStream& write(OutputStream& os) const;
    InputStream& read(InputStream& is);

    MatrixRow<const Matrix<X>> row(SizeType i) const;
    MatrixColumn<const Matrix<X>> column(SizeType j) const;
  public:
    static Matrix<X> _mul(const Matrix<X>& A1, const Matrix<X>& A2);
};


template<class X> struct IsMatrix<Matrix<X>> : True { };
template<class X> struct IsMatrixExpression<Matrix<X>> : True { };

template<class X1, class X2, EnableIf<IsScalar<X1>> =dummy> Matrix<ProductType<X1,X2>> operator*(X1 const& s, Matrix<X2> const& A);
template<class X1, class X2, EnableIf<IsScalar<X2>> =dummy> Matrix<ProductType<X1,X2>> operator*(Matrix<X1> const& A, X2 const& s);
template<class X1, class X2, EnableIf<IsScalar<X2>> =dummy> Matrix<QuotientType<X1,X2>> operator/(Matrix<X1> const& A, X2 const& s);

template<class X> Matrix<X> operator+(Matrix<X> A);
template<class X> Matrix<X> operator-(Matrix<X> A);
template<class X1, class X2> Matrix<SumType<X1,X2>> operator+(Matrix<X1> const& A1, Matrix<X2> const& A2);
template<class X1, class X2> Matrix<DifferenceType<X1,X2>> operator-(Matrix<X1> const& A1, Matrix<X2> const& A2);
template<class X1, class X2> Matrix<ArithmeticType<X1,X2>> operator*(Matrix<X1> const& A1, Matrix<X2> const& A2);
template<class X1, class X2> Vector<ArithmeticType<X1,X2>> operator*(Matrix<X1> const& A1, Vector<X2> const& v2);
template<class X1, class X2> Covector<ArithmeticType<X1,X2>> operator*(Covector<X1> const& u1, Matrix<X2> const& A2);
template<class X1, class X2> Matrix<ProductType<X1,X2>> operator*(Vector<X1> const& v1, Covector<X2> const& u2);

template<class X1, class X2> decltype(declval<X1>()==declval<X2>()) operator==(Matrix<X1> const& A1, Matrix<X2> const& A2);

template<class M1, class M2, EnableIf<And<IsMatrixExpression<M1>,IsMatrixExpression<M2>>> =dummy>
auto operator==(M1 const& A1, M2 const& A2) -> decltype(declval<ScalarType<M1>>()==declval<ScalarType<M2>>());

template<class X> OutputStream& operator<<(OutputStream& os, Matrix<X> const& A);
template<class M, EnableIf<IsMatrixExpression<M>> =dummy> OutputStream& operator<<(OutputStream& os, M const& A);

/************ Combinatorial Matrices *********************************************************/

struct PivotMatrix {
    Array<SizeType> _ary;
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


/************ Matrix inlines *********************************************************/

/************ Matrix expressions *********************************************************/

template<class M> class MatrixRow {
    M* _Ap; SizeType _i;
  public:
    typedef typename M::ScalarType ScalarType;
    MatrixRow(M* Ap, SizeType i) : _Ap(Ap), _i(i) { }
    SizeType size() const { return _Ap->column_size(); }
    auto operator[](SizeType j) -> decltype(_Ap->at(_i,j)) { return _Ap->at(_i,j); }
};
template<class M> struct IsCovectorExpression<MatrixRow<M>> : True { };

template<class M> class MatrixColumn {
    M* _Ap; SizeType _j;
  public:
    typedef typename M::ScalarType ScalarType;
    MatrixColumn(M* Ap, SizeType j) : _Ap(Ap), _j(j) { }
    SizeType size() const { return _Ap->row_size(); }
    auto operator[](SizeType i) -> decltype(_Ap->at(i,_j)) { return _Ap->at(i,_j); }
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
};
template<class M> struct IsMatrixExpression<MatrixTranspose<M>> : True { };

template<class M1, class M2> struct MatrixMatrixProduct {
    typedef ArithmeticType<typename M1::ScalarType,typename M2::ScalarType> ScalarType;
    M1 const& _A1; M2 const& _A2;
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
     M1 const& _A1; V2 const& _v2;
    MatrixVectorProduct(M1 const& A1, V2 const& v2) : _A1(A1), _v2(v2) { }
    SizeType size() const { return _A1.row_size(); }
    ScalarType zero_element() const { return _A1.zero_element()*_v2.zero_element(); }
    ScalarType operator[](SizeType i) const { ScalarType r=this->zero_element();
        for(SizeType j=0; j!=_v2.size(); ++j) { r+=_A1.at(i,j)*_v2.at(j); } return std::move(r); }
};
template<class M1,class V2> struct IsVectorExpression<MatrixVectorProduct<M1,V2>> : True { };

template<class M1, class X2> struct MatrixScalarQuotient {
    typedef QuotientType<typename M1::ScalarType,X2> ScalarType;
    const M1& _a1; const X2& _x2;
    MatrixScalarQuotient(const M1& a1, const X2& x2) : _a1(a1), _x2(x2) { }
    SizeType row_size() const { return _a1.row_size(); }
    SizeType column_size() const { return _a1.column_size(); }
    auto at(SizeType i, SizeType j) const -> decltype(_a1.at(i,j)/_x2) { return _a1.at(i,j)/_x2; }
};
template<class M1, class X2> struct IsMatrixExpression<MatrixScalarQuotient<M1,X2>> : True { };

/* Dispatching Matrix expression template operators

template<class M1, class M2, EnableIf<And<IsMatrixExpression<M1>,IsMatrixExpression<M2>>>> inline
auto operator==(M1 const& A1, M2 const& A2) -> decltype(declval<ScalarType<M1>>()==declval<ScalarType<M2>>()) {
    typedef ScalarType<M1> X1;typedef ScalarType<M2> X2; return Matrix<X1>(A1)==Matrix<X2>(A2); }


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

template<class X> inline MatrixRow<const Matrix<X>> Matrix<X>::operator[](SizeType i) const {
    return MatrixRow<const Matrix<X>>{this,i};
}

template<class X> inline MatrixRow<Matrix<X>> Matrix<X>::operator[](SizeType i) {
    return MatrixRow<Matrix<X>>{this,i};
}

template<class X> inline const X& Matrix<X>::at(SizeType i, SizeType j) const {
    return this->_ary[i*this->_cs+j];
}

template<class X> inline X& Matrix<X>::at(SizeType i, SizeType j) {
    return this->_ary[i*this->_cs+j];
}

template<class X> inline Void Matrix<X>::set(SizeType i, SizeType j, const X& c) {
    this->_ary[i*this->_cs+j]=c;
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

template<class X> template<class M, EnableIf<And<IsMatrixExpression<M>,IsConvertible<typename M::ScalarType,X>>>> Matrix<X>& Matrix<X>::operator=(const M& A) {
    this->resize(A.row_size(),A.column_size());
    for(SizeType i=0; i!=this->row_size(); ++i) {
        for(SizeType j=0; j!=this->column_size(); ++j) {
            this->at(i,j)=A.at(i,j);
        }
    }
    return *this;
}

template<class X> OutputStream& operator<<(OutputStream& os, Matrix<X> const& A) {
    return A.write(os);
}

template<class M, EnableIf<IsMatrixExpression<M>>> inline OutputStream& operator<<(OutputStream& os, const M& A) {
    return os << Matrix<ScalarType<M>>(A); }


template<class X> template<class M> Matrix<X>& Matrix<X>::operator+=(const M& A) {
    for(SizeType i=0; i!=this->row_size(); ++i) {
        for(SizeType j=0; j!=this->column_size(); ++j) {
            this->at(i,j)+=A.at(i,j);
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

template<class X> inline Matrix<X> operator+(Matrix<X> A) {
    return std::move(A);
}

template<class X> inline Matrix<X> operator-(Matrix<X> A) {
    for(SizeType i=0; i!=A.row_size(); ++i) {
        for(SizeType j=0; j!=A.column_size(); ++j) {
            A[i][j]=-A[i][j];
        }
    }
    return std::move(A);
}

template<class X1, class X2> inline Matrix<SumType<X1,X2>> operator+(Matrix<X1> const& A1, Matrix<X2> const& A2) {
    ARIADNE_PRECONDITION(A1.row_size()==A2.row_size());
    ARIADNE_PRECONDITION(A1.column_size()==A2.column_size());
    Matrix<SumType<X1,X2>> R(A1.row_size(),A1.column_size());
    for(SizeType i=0; i!=A1.row_size(); ++i) {
        for(SizeType j=0; j!=A1.column_size(); ++j) {
            R[i][j]=A1[i][j]+A2[i][j];
        }
    }
    return std::move(R);
}

template<class X1, class X2> inline Matrix<DifferenceType<X1,X2>> operator-(Matrix<X1> const& A1, Matrix<X2> const& A2) {
    ARIADNE_PRECONDITION(A1.row_size()==A2.row_size());
    ARIADNE_PRECONDITION(A1.column_size()==A2.column_size());
    Matrix<DifferenceType<X1,X2>> R(A1.row_size(),A1.column_size());
    for(SizeType i=0; i!=A1.row_size(); ++i) {
        for(SizeType j=0; j!=A1.column_size(); ++j) {
            R[i][j]=A1[i][j]-A2[i][j];
        }
    }
    return std::move(R);
}

template<class X1, class X2, EnableIf<IsScalar<X1>>> inline Matrix<ProductType<X1,X2>> operator*(X1 const& s1, Matrix<X2> const& A2) {
    Matrix<ProductType<X1,X2>> R(A2.row_size(),A2.row_size());
    for(SizeType i=0; i!=A2.row_size(); ++i) {
        for(SizeType j=0; j!=A2.column_size(); ++j) {
            R[i][j]=s1*A2[i][j];
        }
    }
    return std::move(R);
}

template<class X1, class X2, EnableIf<IsScalar<X2>>> inline Matrix<ProductType<X1,X2>> operator*(Matrix<X1> const& A1, X2 const& s2) {
    Matrix<ProductType<X1,X2>> R(A1.row_size(),A1.row_size());
    for(SizeType i=0; i!=A1.row_size(); ++i) {
        for(SizeType j=0; j!=A1.column_size(); ++j) {
            R[i][j]=A1[i][j]*s2;
        }
    }
    return std::move(R);
}

template<class X1, class X2,EnableIf<IsScalar<X2>>> inline Matrix<QuotientType<X1,X2>> operator/(Matrix<X1> const& A1, X2 const& s2) {
    Matrix<QuotientType<X1,X2>> R(A1.row_size(),A1.row_size());
    for(SizeType i=0; i!=A1.row_size(); ++i) {
        for(SizeType j=0; j!=A1.column_size(); ++j) {
            R[i][j]=A1[i][j]/s2;
        }
    }
    return std::move(R);
}

//template<class X1, class X2> inline MatrixMatrixProduct<X1,X2> operator*(Matrix<X1> const& A1, Matrix<X2> const& A2) {
//    return MatrixMatrixProduct<X1,X2>(A1,A2); }

template<class X1, class X2> Matrix<ArithmeticType<X1,X2>> operator*(Matrix<X1> const& A1, Matrix<X2> const& A2) {
    typedef ArithmeticType<X1,X2> X0;
    ARIADNE_PRECONDITION(A1.column_size()==A2.row_size());
    Matrix<X0> A0(A1.row_size(), A2.column_size());
    for(SizeType i=0; i!=A0.row_size(); ++i) {
        for(SizeType j=0; j!=A0.column_size(); ++j) {
            for(SizeType k=0; k!=A1.column_size(); ++k) {
                A0.at(i,j)+=A1.at(i,k)*A2.at(k,j);
            }
        }
    }
    return std::move(A0);
}

template<class X1, class X2> Vector<ArithmeticType<X1,X2>> operator*(Matrix<X1> const& A1, Vector<X2> const& v2) {
    typedef ArithmeticType<X1,X2> X0;
    ARIADNE_PRECONDITION(A1.column_size()==v2.size());
    Vector<X0> v0(A1.row_size());
    for(SizeType i=0; i!=v0.size(); ++i) {
        for(SizeType j=0; j!=v2.size(); ++j) {
            v0.at(i)+=A1.at(i,j)*v2.at(j);
        }
    }
    return std::move(v0);
}

template<class X1, class X2> Covector<ArithmeticType<X1,X2>> operator*(Covector<X1> const& u1, Matrix<X2> const& A2) {
    typedef ArithmeticType<X1,X2> X0;
    ARIADNE_PRECONDITION(u1.size()==A2.row_size());
    Covector<X0> u0(A2.column_size());
    for(SizeType j=0; j!=u0.size(); ++j) {
        for(SizeType i=0; i!=u1.size(); ++i) {
            u0.at(j)+=u1.at(i)*A2.at(i,j);
        }
    }
    return std::move(u0);
}

template<class X1, class X2> inline
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

template<class X> auto norm(Matrix<X> const& A) -> decltype(abs(declval<X>())+abs(declval<X>()))  {
    typedef decltype(abs(declval<X>())+abs(declval<X>())) R;
    R r=0;
    for(SizeType i=0; i!=A.row_size(); ++i) {
        R s=0;
        for(SizeType j=0; j!=A.column_size(); ++j) {
            s+=abs(A.at(i,j));
        }
        r=max(r,s);
    }
    return r;
}

template<class X> inline Matrix<X>& operator*=(Matrix<X>& A, typename Matrix<X>::ScalarType const& s) {
    for(SizeType i=0; i!=A.row_size(); ++i) {
        for(SizeType j=0; j!=A.column_size(); ++j) {
            A.at(i,j)*=s;
        }
    }
    return A;
}

template<class X> MatrixTranspose<Matrix<X>> transpose(const Matrix<X>& A) {
    return MatrixTranspose<Matrix<X>>(A);
}

template<class X1,class X2> Matrix<ArithmeticType<X1,X2>> operator*(MatrixTranspose<Matrix<X1>> const& A1, Matrix<X2> A2) {
    Matrix<X1> const& A1T=A1._AT;
    ARIADNE_PRECONDITION(A1T.row_size()==A2.row_size());
    Matrix<ArithmeticType<X1,X2>> A0(A1T.column_size(),A2.column_size());
    for(SizeType i=0; i!=A1T.column_size(); ++i) {
        for(SizeType j=0; j!=A2.column_size(); ++j) {
            for(SizeType k=0; k!=A1.row_size(); ++k) {
                A0.at(i,j)+=A1T.at(k,i)*A2.at(k,j);
            }
        }
    }
    return std::move(A0);
}

template<class X1, class X2> Vector<ArithmeticType<X1,X2>> operator*(MatrixTranspose<Matrix<X1>> const& A1, Vector<X2> v2) {
    Matrix<X1> const& A1T=A1._AT;
    ARIADNE_PRECONDITION(A1T.row_size()==v2.size());
    Vector<ArithmeticType<X1,X2>> v0(A1T.column_size());
    for(SizeType i=0; i!=v0.size(); ++i) {
        for(SizeType j=0; j!=v2.size(); ++j) {
            v0.at(i)+=A1T.at(j,i)*v2.at(j);
        }
    }
    return std::move(v0);
}


template<class X> Matrix<X> inverse(Matrix<X> const& A);
template<class X> Vector<X> solve(Matrix<X> const& A, Vector<X> const& b);

template<class X> Vector<X> gs_solve(Matrix<X> const& A, Vector<X> const& b);
template<class X> Vector<X> gs_solve(Matrix<X> const& A, Vector<X> const& b, Vector<X> iX);
template<class X> Void gs_step(Matrix<X> const& A, Vector<X> const& b, Vector<X>& iX);

template<class AX> Matrix<decltype(make_exact(declval<AX>()))> make_exact(Matrix<AX> const&);

inline Matrix<ExactFloat64>& make_exact(Matrix<ApprxFloat64>& mx) {
    return reinterpret_cast<Matrix<ExactFloat64>&>(mx); }


} // namespace Ariadne

#endif
