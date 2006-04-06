/***************************************************************************
 *            tensor.cc
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
    
    template class tensor<Float64>;
    template class tensor<Real>;
    template class tensor<Field>;

    template vector<Float64> product(const tensor<Float64>&, const vector<Float64>&, const vector<Float64>&);
    template vector<Real> product(const tensor<Real>&, const vector<Real>&, const vector<Real>&);
    template vector<Field> product(const tensor<Field>&, const vector<Field>&, const vector<Field>&);
    
    template matrix<Float64> product(const tensor<Float64>&, const vector<Float64>&);
    template matrix<Real> product(const tensor<Real>&, const vector<Real>&);
    template matrix<Field> product(const tensor<Field>&, const vector<Field>&);
    
    template tensor<Float64> product(const tensor<Float64>&, const matrix<Float64>&);
    template tensor<Real> product(const tensor<Real>&, const matrix<Real>&);
    template tensor<Field> product(const tensor<Field>&, const matrix<Field>&);
    
    template tensor<Float64> product(const tensor<Float64>&, const matrix<Float64>&, const matrix<Float64>&);
    template tensor<Real> product(const tensor<Real>&, const matrix<Real>&, const matrix<Real>&);
    template tensor<Field> product(const tensor<Field>&, const matrix<Field>&, const matrix<Field>&);
    
    template std::ostream& operator<<(std::ostream&, const tensor<Float64>&);
    template std::ostream& operator<<(std::ostream&, const tensor<Real>&);
    template std::ostream& operator<<(std::ostream&, const tensor<Field>&);
    
  }
}
