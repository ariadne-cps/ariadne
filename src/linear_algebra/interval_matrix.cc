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
 
#include "real_typedef.h"

#include "linear_algebra/interval_matrix.h"
#include "linear_algebra/interval_matrix.tpl"

namespace Ariadne {
  namespace LinearAlgebra {
    template class interval_matrix<Real>;
    template class interval_matrix<Field>;
      
    template interval_vector<Real> prod(const matrix<Real>&, const interval_vector<Real>&);
    template interval_vector<Field> prod(const matrix<Field>&, const interval_vector<Field>&);
    
    template interval_vector<Real> prod(const interval_matrix<Real>&, const vector<Real>&);
    template interval_vector<Field> prod(const interval_matrix<Field>&, const vector<Field>&);
    
    template interval_vector<Real> prod(const interval_matrix<Real>&, const interval_vector<Real>&);
    template interval_vector<Field> prod(const interval_matrix<Field>&, const interval_vector<Field>&);
    
    template interval_matrix<Real> prod(const matrix<Real>&, const interval_matrix<Real>&);
    template interval_matrix<Field> prod(const matrix<Field>&, const interval_matrix<Field>&);
    
    template interval_matrix<Real> prod(const interval_matrix<Real>&, const matrix<Real>&);
    template interval_matrix<Field> prod(const interval_matrix<Field>&, const matrix<Field>&);
    
    template interval_matrix<Real> prod(const interval_matrix<Real>&, const interval_matrix<Real>&);
    template interval_matrix<Field> prod(const interval_matrix<Field>&, const interval_matrix<Field>&);
    
    template interval_matrix<Real> fprod(const matrix<Field>&, const interval_matrix<Real>&);
     
    template std::ostream& operator<<(std::ostream&, const interval_matrix<Real>&);
    template std::ostream& operator<<(std::ostream&, const interval_matrix<Field>&);
    
    template matrix<Real> centre(const interval_matrix<Real>&);
    template matrix<Field> centre(const interval_matrix<Field>&);

    template Real radius(const interval_matrix<Real>&);
    template Field radius(const interval_matrix<Field>&);

    template Interval<Real> norm(const interval_matrix<Real>&);
    template Interval<Field> norm(const interval_matrix<Field>&);

    template Real upper_log_norm(const interval_matrix<Real>&);
    template Field upper_log_norm(const interval_matrix<Field>&);
       
    template matrix<Real> over_approximation(const interval_matrix<Real>&);
    template matrix<Field> over_approximation(const interval_matrix<Field>&);

    template interval_matrix<Real> approximate(const matrix<Field>& A,const Real& e); 

  }
}
