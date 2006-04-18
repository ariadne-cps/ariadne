/***************************************************************************
 *            Tensor.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
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

#include "real_typedef.h"

#include "linear_algebra/tensor.h"
#include "linear_algebra/tensor.tpl"

namespace Ariadne {
  namespace LinearAlgebra {
    
    template class Tensor<Float64>;
    template class Tensor<Real>;
    template class Tensor<Field>;

    template Vector<Float64> product(const Tensor<Float64>&, const Vector<Float64>&, const Vector<Float64>&);
    template Vector<Real> product(const Tensor<Real>&, const Vector<Real>&, const Vector<Real>&);
    template Vector<Field> product(const Tensor<Field>&, const Vector<Field>&, const Vector<Field>&);
    
    template Matrix<Float64> product(const Tensor<Float64>&, const Vector<Float64>&);
    template Matrix<Real> product(const Tensor<Real>&, const Vector<Real>&);
    template Matrix<Field> product(const Tensor<Field>&, const Vector<Field>&);
    
    template Tensor<Float64> product(const Tensor<Float64>&, const Matrix<Float64>&);
    template Tensor<Real> product(const Tensor<Real>&, const Matrix<Real>&);
    template Tensor<Field> product(const Tensor<Field>&, const Matrix<Field>&);
    
    template Tensor<Float64> product(const Tensor<Float64>&, const Matrix<Float64>&, const Matrix<Float64>&);
    template Tensor<Real> product(const Tensor<Real>&, const Matrix<Real>&, const Matrix<Real>&);
    template Tensor<Field> product(const Tensor<Field>&, const Matrix<Field>&, const Matrix<Field>&);
    
    template std::ostream& operator<<(std::ostream&, const Tensor<Float64>&);
    template std::ostream& operator<<(std::ostream&, const Tensor<Real>&);
    template std::ostream& operator<<(std::ostream&, const Tensor<Field>&);
    
  }
}
