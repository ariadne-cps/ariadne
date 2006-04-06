/***************************************************************************
 *            matrix.cc
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
 
#include "linear_algebra/interval_matrix.h"
#include "linear_algebra/interval_matrix.tpl"

namespace Ariadne {
  namespace LinearAlgebra {
    template class interval_matrix<Dyadic>;
    template class interval_matrix<Rational>;
    template class interval_matrix<Float64>;
      
    template interval_vector<Dyadic> prod(const matrix<Dyadic>&, const interval_vector<Dyadic>&);
    template interval_vector<Rational> prod(const matrix<Rational>&, const interval_vector<Rational>&);
    template interval_vector<Float64> prod(const matrix<Float64>&, const interval_vector<Float64>&);
    
    template interval_vector<Dyadic> prod(const interval_matrix<Dyadic>&, const vector<Dyadic>&);
    template interval_vector<Rational> prod(const interval_matrix<Rational>&, const vector<Rational>&);
    template interval_vector<Float64> prod(const interval_matrix<Float64>&, const vector<Float64>&);
    
    template interval_vector<Dyadic> prod(const interval_matrix<Dyadic>&, const interval_vector<Dyadic>&);
    template interval_vector<Rational> prod(const interval_matrix<Rational>&, const interval_vector<Rational>&);
    template interval_vector<Float64> prod(const interval_matrix<Float64>&, const interval_vector<Float64>&);
    
    template interval_matrix<Dyadic> prod(const matrix<Dyadic>&, const interval_matrix<Dyadic>&);
    template interval_matrix<Rational> prod(const matrix<Rational>&, const interval_matrix<Rational>&);
    template interval_matrix<Float64> prod(const matrix<Float64>&, const interval_matrix<Float64>&);
    
    template interval_matrix<Dyadic> prod(const interval_matrix<Dyadic>&, const matrix<Dyadic>&);
    template interval_matrix<Rational> prod(const interval_matrix<Rational>&, const matrix<Rational>&);
    template interval_matrix<Float64> prod(const interval_matrix<Float64>&, const matrix<Float64>&);
    
    template interval_matrix<Dyadic> prod(const interval_matrix<Dyadic>&, const interval_matrix<Dyadic>&);
    template interval_matrix<Rational> prod(const interval_matrix<Rational>&, const interval_matrix<Rational>&);
    template interval_matrix<Float64> prod(const interval_matrix<Float64>&, const interval_matrix<Float64>&);
    
    template interval_matrix<Dyadic> fprod(const matrix<Rational>&, const interval_matrix<Dyadic>&);
     
    template std::ostream& operator<<(std::ostream&, const interval_matrix<Dyadic>&);
    template std::ostream& operator<<(std::ostream&, const interval_matrix<Rational>&);
    template std::ostream& operator<<(std::ostream&, const interval_matrix<Float64>&);
    
    template matrix<Dyadic> centre(const interval_matrix<Dyadic>&);
    template matrix<Rational> centre(const interval_matrix<Rational>&);
    template matrix<Float64> centre(const interval_matrix<Float64>&);

    template Dyadic radius(const interval_matrix<Dyadic>&);
    template Rational radius(const interval_matrix<Rational>&);
    template Float64 radius(const interval_matrix<Float64>&);

    template Interval<Dyadic> norm(const interval_matrix<Dyadic>&);
    template Interval<Rational> norm(const interval_matrix<Rational>&);
    template Interval<Float64> norm(const interval_matrix<Float64>&);

    template Dyadic upper_log_norm(const interval_matrix<Dyadic>&);
    template Rational upper_log_norm(const interval_matrix<Rational>&);
    template Float64 upper_log_norm(const interval_matrix<Float64>&);
      
    template matrix<Dyadic> over_approximation(const interval_matrix<Dyadic>&);
    template matrix<Rational> over_approximation(const interval_matrix<Rational>&);
    template matrix<Float64> over_approximation(const interval_matrix<Float64>&);


  }
}
