/***************************************************************************
 *            affine_integrator.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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


#include "numeric/float.h"

#include "evaluation/affine_integrator.h"
#include "evaluation/affine_integrator.code.h"

namespace Ariadne {
  namespace Evaluation {
    using namespace Numeric;

#ifdef ENABLE_FLOAT64
   template Float64 gexp_up(const Float64& x, uint k);

    template LinearAlgebra::Vector<Interval64> 
    gexp(const LinearAlgebra::Matrix<Interval64>& A, 
         const LinearAlgebra::Vector<Interval64>& b, 
         const Interval64& t, 
         const uint& k);
    
    template LinearAlgebra::Matrix<Interval64>
    gexp(const LinearAlgebra::Matrix<Interval64>& A, const Interval64& t, const uint& k);
    
  template class AffineIntegrator< Zonotope<Float64> >;
#endif
  
#ifdef ENABLE_FLOATMP
   template FloatMP gexp_up(const FloatMP& x, uint k);

    template LinearAlgebra::Vector<IntervalMP> 
    gexp(const LinearAlgebra::Matrix<IntervalMP>& A, 
         const LinearAlgebra::Vector<IntervalMP>& b, 
         const IntervalMP& t, const uint& k);
    
    template LinearAlgebra::Matrix<IntervalMP> 
    gexp(const LinearAlgebra::Matrix<IntervalMP>& A, 
         const IntervalMP& t, const uint& k);
    
    template class AffineIntegrator< Zonotope<FloatMP> >;
#endif

  }
}
