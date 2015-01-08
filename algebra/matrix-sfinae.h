/***************************************************************************
 *            algebra/matrix-sfinae.h
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

/*! \file algebra/matrix-sfinae.h
 *  \brief 
 */



#ifndef ARIADNE_MATRIX_H
#define ARIADNE_MATRIX_H

#include "vector.h"

namespace Ariadne {

/************ Matrix *********************************************************/


template<class M> class MatrixRow {
    M* _Ap; SizeType _i;
  public:
    MatrixRow(M* Ap, SizeType i) : _Ap(Ap), _i(i) { }
    auto operator[](SizeType j) -> decltype(_Ap->at(_i,j)) { return _Ap->at(_i,j); }
};

template<class X> class Matrix {
    X _zero;
    SizeType _rs;
    SizeType _cs;
    Array<X> _ary;
  public:
    Matrix(SizeType m, SizeType n) : _zero(), _rs(m), _cs(n), _ary(m*n) { }
    Matrix(SizeType m, SizeType n, const X& x) : _zero(x*0), _rs(m), _cs(n), _ary(m*n,x) { }
//    Matrix(InitializerList<X> lst) : _zero(*lst.begin()*0), _ary(lst.begin(),lst.end()) { }
    SizeType row_size() const { return _rs; }
    SizeType column_size() const { return _cs; }
    MatrixRow<const Matrix<X>> operator[](SizeType i) const { return MatrixRow<const Matrix<X>>{this,i}; }
    MatrixRow<Matrix<X>> operator[](SizeType i) { return MatrixRow<Matrix<X>>{this,i}; }
    const X& at(SizeType i, SizeType j) const { return this->_ary[i*this->_rs+this->_cs]; }
    X& at(SizeType i, SizeType j) { return this->_ary[i*this->_rs+this->_cs]; }
    Void set(SizeType i, SizeType j, const X& c) const { this->_ary[i*this->_rs+this->_cs]=c; }
    X zero_element() const { return _zero; }
    OutputStream& write(OutputStream& os) const;
};
template<class X> struct IsMatrix<Matrix<X>> : True { };

template<class X> class SymmetricMatrix : public Matrix<X>
{
  public:
    SymmetricMatrix(SizeType n) : Matrix<X>(n,n) { }
};
template<class X> struct IsMatrix<SymmetricMatrix<X>> : True { };


} // namespace Ariadne

#endif
