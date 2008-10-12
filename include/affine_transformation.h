/***************************************************************************
 *            affine_transformation.h
 *
 *  Copyright 2008  Pieter Collins
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
 
/*! \file affine_transformation.h
 *  \brief Affine transformations of Euclidean space.
 */

#ifndef ARIADNE_AFFINE_TRANSFORMATION_H
#define ARIADNE_AFFINE_TRANSFORMATION_H

#include "vector.h"
#include "matrix.h"

namespace Ariadne {

template<class X> class AffineTransformation;

//! The affine transformation \f$x\mapsto A(x-c)+b \f$.
template<class X>
class AffineTransformation
{
  Vector<Float> _centre;
  Vector<X> _b;
  Matrix<X> _A;
 public:
  AffineTransformation(const Vector<Float>& c, const Vector<X>& b, const Matrix<X>& A)
    : _centre(c), _b(b), _A(A) { assert(A.column_size()==c.size()); assert(A.row_size()==b.size()); }
  uint result_size() const { return _b.size(); }
  uint argument_size() const { return _centre.size(); }
  const Vector<Float>& centre() const { return _centre; }
  const Vector<Float>& c() const { return _centre; }
  const Vector<X>& b() const { return _b; }
  const Matrix<X>& A() const { return _A; }
};

} // namespace Ariadne

#endif /* ARIADNE_AFFINE_TRANSFORMATION_H */

