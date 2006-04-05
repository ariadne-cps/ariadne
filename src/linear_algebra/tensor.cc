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
 
#include "linear_algebra/tensor.h"
#include "linear_algebra/tensor.tpl"


namespace Ariadne {
  namespace LinearAlgebra {
    
    template class tensor<Float64>;
    template class tensor<Dyadic>;
    template class tensor<Rational>;

    template vector<Float64> product(const tensor<Float64>&, const vector<Float64>&, const vector<Float64>&);
    template vector<Dyadic> product(const tensor<Dyadic>&, const vector<Dyadic>&, const vector<Dyadic>&);
    template vector<Rational> product(const tensor<Rational>&, const vector<Rational>&, const vector<Rational>&);
    
    template matrix<Float64> product(const tensor<Float64>&, const vector<Float64>&);
    template matrix<Dyadic> product(const tensor<Dyadic>&, const vector<Dyadic>&);
    template matrix<Rational> product(const tensor<Rational>&, const vector<Rational>&);
    
    template tensor<Float64> product(const tensor<Float64>&, const matrix<Float64>&);
    template tensor<Dyadic> product(const tensor<Dyadic>&, const matrix<Dyadic>&);
    template tensor<Rational> product(const tensor<Rational>&, const matrix<Rational>&);
    
    template tensor<Float64> product(const tensor<Float64>&, const matrix<Float64>&, const matrix<Float64>&);
    template tensor<Dyadic> product(const tensor<Dyadic>&, const matrix<Dyadic>&, const matrix<Dyadic>&);
    template tensor<Rational> product(const tensor<Rational>&, const matrix<Rational>&, const matrix<Rational>&);
    
    template std::ostream& operator<<(std::ostream&, const tensor<Float64>&);
    template std::ostream& operator<<(std::ostream&, const tensor<Dyadic>&);
    template std::ostream& operator<<(std::ostream&, const tensor<Rational>&);
    
  }
}
