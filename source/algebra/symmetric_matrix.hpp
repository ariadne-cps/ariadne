/***************************************************************************
 *            algebra/symmetric_matrix.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file algebra/symmetric_matrix.hpp
 *  \brief
 */



#ifndef ARIADNE_SYMMETRIC_MATRIX_HPP
#define ARIADNE_SYMMETRIC_MATRIX_HPP

#include <initializer_list>

#include "vector.hpp"
#include "matrix.hpp"

namespace Ariadne {

/************ SymmetricMatrix *********************************************************/

template<class X> class SymmetricMatrix : public Matrix<X>
{
  public:
    SymmetricMatrix(SizeType n) : Matrix<X>(n,n) { }
    template<class X1, class X2> friend Covector<ArithmeticType<X1,X2>> operator*(SymmetricMatrix<X1> const& S, Vector<X2> const& v);
};


} // namespace Ariadne

#endif
