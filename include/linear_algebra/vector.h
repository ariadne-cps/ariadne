/***************************************************************************
 *            vector.h
 *
 *  Mon May  3 12:31:15 2004
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file vector.h
 *  \brief Vectors and vector operations.
  */

#ifndef _ARIADNE_VECTOR_H
#define _ARIADNE_VECTOR_H 

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <base/basic_type.h>
#include <base/numerical_type.h>
#include <base/interval.h>

namespace Ariadne {
  namespace LinearAlgebra {

    using boost::numeric::ublas::vector;
        
    template <typename Real>
    inline vector<Real> zero_vector(dimension_type dim) {
      vector<Real> v(dim);
      for (dimension_type j=0; j< dim; j++) {
        v(j)=0.0;
      }
      return v;
    }
    
    template <typename Real>
    inline Integer common_denominator(const vector<Real>& b) 
    {
      Integer denom=1;
      for (dimension_type i=0; i< b.size(); ++i) {
        denom=lcm( denom, denominator(b(i)) );
      }
      return denom;
    }

  }
}


#endif /* _ARIADNE_VECTOR_H */
