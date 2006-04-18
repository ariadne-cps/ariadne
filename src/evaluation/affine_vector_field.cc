/***************************************************************************
 *            affine_vector_field.cc
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

#include "evaluation/affine_vector_field.h"
#include "evaluation/affine_vector_field.tpl"

#include "real_typedef.h"

namespace Ariadne {
  namespace Evaluation {

    template class AffineVectorField<Real>;

  }
}


namespace Ariadne {
  namespace LinearAlgebra {
    
    template Matrix<Real> exp_Ah_approx(const Matrix<Real>& A, 
                                          const Real& h, 
                                          const Real& e);
    
    template Matrix<Real> exp_Ah_sub_id_div_A_approx(const Matrix<Real>& A, 
                                                       const Real& h, 
                                                       const Real& e);

  }
}
