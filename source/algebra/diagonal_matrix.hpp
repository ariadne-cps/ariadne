/***************************************************************************
 *            algebra/diagonal_matrix.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file algebra/diagonal_matrix.hpp
 *  \brief
 */



#ifndef ARIADNE_DIAGONAL_MATRIX_HPP
#define ARIADNE_DIAGONAL_MATRIX_HPP

#include <initializer_list>

#include "vector.hpp"
#include "matrix.hpp"

namespace Ariadne {

template<class T> using InitializerList = InitializerList<T>;

template<class X> class Matrix;
template<class X> class SymmetricMatrix;


template<class X> class DiagonalMatrix;

template<class X, class G> requires InvocableReturning<X,G,SizeType,SizeType>
Matrix<X>& assign(Matrix<X>& A, G const& g) {
    for (SizeType i=0; i!=A.row_size(); ++i) {
        for (SizeType j=0; j!=A.column_size(); ++j) {
            A.at(i,j)=g(i,j);
        }
    }
    return A;
}

struct DiagonalMatrixOperations {

    //! \brief <p/>
    template<class X> friend OutputStream& operator<<(OutputStream& os, DiagonalMatrix<X> const& D) {
        return D._write(os);
    }

    //! \brief <p/>
    template<class X> friend DiagonalMatrix<NegationType<X>> operator-(DiagonalMatrix<X> D) {
        if constexpr (Same<NegationType<X>,X>) {
            for(SizeType i=0; i!=D.size(); ++i) { D._at(i)=-D._at(i); } return D;
        } else {
            return DiagonalMatrix<NegationType<X>>(D.size(),[&D](SizeType i){return -D._at(i);});
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend DiagonalMatrix<SumType<X1,X2>> operator+(DiagonalMatrix<X1> D1, DiagonalMatrix<X2> const& D2) {
        ARIADNE_PRECONDITION(D1.size()==D2.size());
        if constexpr (Same<SumType<X1,X2>,X1>) {
            for(SizeType i=0; i!=D1.size(); ++i) { D1._at(i)+=D2._at(i); } return D1;
        } else {
            return DiagonalMatrix<SumType<X1,X2>>(D1.size(),[&D1,&D2](SizeType i){return D1._at(i)+D2._at(i);});
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend DiagonalMatrix<DifferenceType<X1,X2>> operator-(DiagonalMatrix<X1> D1, DiagonalMatrix<X2> const& D2) {
        ARIADNE_PRECONDITION(D1.size()==D2.size());
        if constexpr (Same<DifferenceType<X1,X2>,X1>) {
            for(SizeType i=0; i!=D1.size(); ++i) { D1._at(i)-=D2._at(i); } return D1;
        } else {
            return DiagonalMatrix<DifferenceType<X1,X2>>(D1.size(),[&D1,&D2](SizeType i){return D1._at(i)-D2._at(i);});
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend DiagonalMatrix<ProductType<X1,X2>> operator*(X1 const& s1, DiagonalMatrix<X2> D2) {
        if constexpr (Same<ProductType<X1,X2>,X2>) {
            for(SizeType i=0; i!=D2.size(); ++i) { D2._at(i)*=s1; } return D2;
        } else {
            return DiagonalMatrix<ProductType<X1,X2>>(D2.size(),[&s1,&D2](SizeType i){return s1*D2._at(i);});
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend DiagonalMatrix<ProductType<X1,Scalar<X2>>> operator*(DiagonalMatrix<X1> D1, X2 const& s2) {
        if constexpr (Same<ProductType<X1,X2>,X1>) {
            for(SizeType i=0; i!=D1.size(); ++i) { D1._at(i)*=s2; } return D1;
        } else {
            return DiagonalMatrix<ProductType<X1,X2>>(D1.size(),[&D1,&s2](SizeType i){return D1._at(i)*s2;});
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend DiagonalMatrix<QuotientType<X1,X2>> operator/(DiagonalMatrix<X1> D1, X2 const& s2) {
        if constexpr (Same<QuotientType<X1,X2>,X1>) {
            for(SizeType i=0; i!=D1.size(); ++i) { D1._at(i)/=s2; } return D1;
        } else {
            return DiagonalMatrix<QuotientType<X1,X2>>(D1.size(),[&D1,&s2](SizeType i){return D1._at(i)/s2;});
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend DiagonalMatrix<ProductType<X1,X2>> operator*(DiagonalMatrix<X1> D1, DiagonalMatrix<X2> const& D2) {
        ARIADNE_PRECONDITION(D1.size()==D2.size());
        if constexpr (Same<ProductType<X1,X2>,X1>) {
            for(SizeType i=0; i!=D1.size(); ++i) { D1._at(i)*=D2._at(i); } return D1;
        } else {
            return DiagonalMatrix<ProductType<X1,X2>>(D1.size(),[&D1,&D2](SizeType i){return D1._at(i)*D2._at(i);});
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend DiagonalMatrix<QuotientType<X1,X2>> operator/(DiagonalMatrix<X1> D1, DiagonalMatrix<X2> const& D2) {
        ARIADNE_PRECONDITION(D1.size()==D2.size());
        if constexpr (Same<QuotientType<X1,X2>,X1>) {
            for(SizeType i=0; i!=D1.size(); ++i) { D1._at(i)/=D2._at(i); } return D1;
        } else {
            return DiagonalMatrix<QuotientType<X1,X2>>(D1.size(),[&D1,&D2](SizeType i){return D1._at(i)/D2._at(i);});
        }
    }

    //! \brief <p/>
    template<class X> friend DiagonalMatrix<ReciprocalType<X>> inverse(DiagonalMatrix<X> D) {
        if constexpr (Same<ReciprocalType<X>,X>) {
            for(SizeType i=0; i!=D.size(); ++i) { D._at(i)=rec(D._at(i)); } return D;
        } else {
            return DiagonalMatrix<ReciprocalType<X>>(D.size(),[&D](SizeType i){return rec(D._at(i));});
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend Matrix<SumType<X1,X2>> operator+(Matrix<X1> A, DiagonalMatrix<X2> const& D) {
        ARIADNE_PRECONDITION(A.row_size()==D.size() && A.column_size()==D.size());
        if constexpr (Same<SumType<X1,X2>,X1>) {
            for(SizeType i=0; i!=D.size(); ++i) { A.at(i,i)+=D._at(i); } return A;
        } else {
            Matrix<SumType<X1,X2>> R=A; for(SizeType i=0; i!=D.size(); ++i) { R.at(i,i)+=D._at(i); } return R;
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend Matrix<SumType<X1,X2>> operator+(DiagonalMatrix<X1> const& D, Matrix<X2> A) {
        ARIADNE_PRECONDITION(A.row_size()==D.size() && A.column_size()==D.size());
        if constexpr (Same<SumType<X1,X2>,X2>) {
            for(SizeType i=0; i!=D.size(); ++i) { A.at(i,i)+=D._at(i); } return A;
        } else {
            Matrix<SumType<X1,X2>> R=A; for(SizeType i=0; i!=D.size(); ++i) { R.at(i,i)+=D._at(i); } return R;
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend Matrix<DifferenceType<X1,X2>> operator-(Matrix<X1> A, DiagonalMatrix<X2> const& D) {
        ARIADNE_PRECONDITION(A.row_size()==D.size() && A.column_size()==D.size());
        if constexpr (Same<DifferenceType<X1,X2>,X1>) {
            for(SizeType i=0; i!=D.size(); ++i) { A.at(i,i)-=D._at(i); } return A;
        } else {
            Matrix<SumType<X1,X2>> R=A; for(SizeType i=0; i!=D.size(); ++i) { R.at(i,i)-=D._at(i); } return R;
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend Matrix<DifferenceType<X1,X2>> operator-(DiagonalMatrix<X1> const& D, Matrix<X2> A) {
        ARIADNE_PRECONDITION(A.row_size()==D.size() && A.column_size()==D.size());
        if constexpr (Same<DifferenceType<X1,X2>,X2>) {
            for(SizeType i=0; i!=A.row_size(); ++i) {
                for(SizeType j=0; j!=A.column_size(); ++j) { A.at(i,j)=-A.at(i,j); }
                A.at(i,i)+=D._at(i);
            }
            return A;
        } else {
            Matrix<DifferenceType<X1,X2>> R(A.row_size(),A.column_size(),[&A](SizeType i, SizeType j){return -A.at(i,j);});
            for(SizeType i=0; i!=D.size(); ++i) { R.at(i,i)+=D._at(i); } return R;
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend Matrix<ProductType<X1,X2>> operator*(Matrix<X1> A, DiagonalMatrix<X2> const& D) {
        ARIADNE_PRECONDITION(A.column_size()==D.size());
        auto g = [&A,&D](SizeType i, SizeType j){return A.at(i,j)*D._at(j);};
        if constexpr (Same<ProductType<X1,X2>,X1>) { assign(A,g); return A; }
        else { return Matrix<ProductType<X1,X2>>(A.row_size(),A.column_size(),g); }
    }

    //! \brief <p/>
    template<class X1, class X2> friend Matrix<ProductType<X1,X2>> operator*(DiagonalMatrix<X1> const& D, Matrix<X2> A) {
        ARIADNE_PRECONDITION(D.size()==A.row_size())
        auto g = [&D,&A](SizeType i, SizeType j){return D._at(i)*A.at(i,j);};
        if constexpr (Same<ProductType<X1,X2>,X2>) { assign(A,g); return A; }
        else { return Matrix<ProductType<X1,X2>>(A.row_size(),A.column_size(),g); }
    }

    //! \brief <p/>
    template<class X1, class X2> friend Matrix<ProductType<X1,X2>> operator*(DiagonalMatrix<X1> const& D, Transpose<Matrix<X2>> const& A) {
        ARIADNE_PRECONDITION(D.size()==A.row_size())
        auto g = [&D,&A](SizeType i, SizeType j){return D._at(i)*A.at(i,j);};
        return Matrix<ProductType<X1,X2>>(A.row_size(),A.column_size(),g);
    }

    //! \brief <p/>
    template<class X1, class X2> friend Vector<ProductType<X1,X2>> operator*(DiagonalMatrix<X1> const& D, Vector<X2> v) {
        ARIADNE_PRECONDITION(D.size()==v.size())
        if constexpr (Same<ProductType<X1,X2>,X2>) {
            for(SizeType i=0; i!=v.size(); ++i) {
                v.at(i)*=D._at(i);
            }
            return v;
        } else {
            return Vector<ProductType<X1,X2>>(v.size(),[&D,&v](SizeType i){return D._at(i)*v.at(i);});
        }
    }

    //! \brief <p/>
    template<class X1, class X2> friend Vector<ProductType<X1,X2>> operator*(Vector<X1> v, DiagonalMatrix<X2> const& D) {
        ARIADNE_PRECONDITION(D.size()==v.size())
        if constexpr (Same<ProductType<X1,X2>,X1>) {
            for(SizeType i=0; i!=v.size(); ++i) {
                v.at(i)*=D._at(i);
            }
            return v;
        } else {
            return Vector<ProductType<X1,X2>>(v.size(),[&v,&D](SizeType i){return v.at(i)*D._at(i);});
        }
    }

    //! \brief <p/>
    template<class X> friend Vector<X> operator/(Vector<X> v, DiagonalMatrix<X> const& D) {
        ARIADNE_PRECONDITION(D.size()==v.size())
        for(SizeType i=0; i!=v.size(); ++i) {
            v.at(i)/=D._at(i);
        }
        return v;
    }

    //! \brief The product \f$A^TDA\f$.
    template<class X1, class X2> friend SymmetricMatrix<ArithmeticType<X1,X2>> outer(Matrix<X1> const& A, DiagonalMatrix<X2> D);

};


//! \ingroup LinearAlgebraModule
//! \brief Diagonal matrices over some type \a X.
template<class X> class DiagonalMatrix
    : DiagonalMatrixOperations
{
    X _zero;
    Array<X> _ary;
  public:
    template<class... PRS> requires Constructible<X,Nat,PRS...> explicit DiagonalMatrix(SizeType n, PRS...); //!< <p/>
    explicit DiagonalMatrix(SizeType n, X const& z); //!< <p/>
    explicit DiagonalMatrix(Array<X>); //!< <p/>
    explicit DiagonalMatrix(Vector<X>);
    template<class Y, class... PRS> requires Constructible<X,Y,PRS...> explicit DiagonalMatrix(DiagonalMatrix<Y> const&, PRS...); //!< <p/>
    template<class G> requires InvocableReturning<X,G,SizeType> DiagonalMatrix(SizeType n, G const& g); //!< <p/>
    SizeType size() const; //!< <p/>
    SizeType row_size() const; //!< <p/>
    SizeType column_size() const; //!< <p/>
    X const& zero_element() const; //!< <p/>
    X const& operator[](SizeType i) const; //!< <p/>
    X const& at(SizeType i, SizeType j) const; //!< <p/>
    X const& get(SizeType i, SizeType j) const; //!< <p/>
    X& operator[](SizeType i); //!< <p/>
    Void set(SizeType i, SizeType j, X const& x); //!< <p/>
    Vector<X> diagonal() const; //!< <p/>
    operator Matrix<X>() const; //!< <p/>
    operator SymmetricMatrix<X>() const; //!< <p/>
  private: public:
    X& _at(SizeType i);
    X const& _at(SizeType i) const;
    OutputStream& _write(OutputStream&) const;
};

template<class X> template<class... PRS> requires Constructible<X,Nat,PRS...>
DiagonalMatrix<X>::DiagonalMatrix(SizeType n, PRS... prs)
    : _zero(0u,prs...), _ary(n,_zero)
{ }

template<class X> template<class G> requires InvocableReturning<X,G,SizeType>
DiagonalMatrix<X>::DiagonalMatrix(SizeType n, G const& g)
    : _zero(nul(g(0))), _ary(n,g)
{
}

template<class X> DiagonalMatrix<X>::DiagonalMatrix(Array<X> ary)
    : _zero(create_zero(ary[0])), _ary(ary)
{ }

template<class X> DiagonalMatrix<X>::DiagonalMatrix(SizeType n, X const& z)
    : _zero(nul(z)), _ary(n,z)
{ }

template<class X> DiagonalMatrix<X>::DiagonalMatrix(Vector<X> vec)
    : _zero(vec.zero_element()), _ary(vec.array())
{ }

template<class X> template<class Y, class... PRS> requires Constructible<X,Y,PRS...>
DiagonalMatrix<X>::DiagonalMatrix(DiagonalMatrix<Y> const& D, PRS... prs)
    : _zero(D.zero_element(),prs...), _ary(D.diagonal().array(),prs...)
{ }

template<class X> SizeType DiagonalMatrix<X>::size() const {
    return this->_ary.size();
}

template<class X> SizeType DiagonalMatrix<X>::row_size() const {
    return this->_ary.size();
}

template<class X> SizeType DiagonalMatrix<X>::column_size() const {
    return this->_ary.size();
}

template<class X> X const& DiagonalMatrix<X>::zero_element() const {
    return this->_zero;
}

template<class X> X const& DiagonalMatrix<X>::operator[](SizeType i) const {
    return this->_ary[i];
}

template<class X> X& DiagonalMatrix<X>::operator[](SizeType i) {
    return this->_ary[i];
}

template<class X> X const& DiagonalMatrix<X>::at(SizeType i, SizeType j) const {
    if (i==j) { return _ary[i]; } else { return _zero; }
}

template<class X> X const& DiagonalMatrix<X>::get(SizeType i, SizeType j) const {
    if (i==j) { return _ary[i]; } else { return _zero; }
}

template<class X> Void DiagonalMatrix<X>::set(SizeType i, SizeType j, X const& x) {
    ARIADNE_PRECONDITION(i==j);
    _ary[i]=x;
}

template<class X> Vector<X> DiagonalMatrix<X>::diagonal () const {
    return Vector<X>(this->_ary);
}

template<class X> DiagonalMatrix<X>::operator Matrix<X> () const {
    Matrix<X> A(this->row_size(),this->column_size(),this->zero_element());
    for(SizeType i=0; i!=this->size(); ++i) { A[i][i]=this->_ary[i]; }
    return A;
}

template<class X> OutputStream& DiagonalMatrix<X>::_write(OutputStream& os) const {
    return os << "diag(" << this->_ary << ")";
}

template<class X> X& DiagonalMatrix<X>::_at(SizeType i) {
    return _ary[i];
}

template<class X> X const& DiagonalMatrix<X>::_at(SizeType i) const {
    return _ary[i];
}

} // namespace Ariadne

#endif
