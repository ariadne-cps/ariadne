/***************************************************************************
 *            affine_vector_field.tpl
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
#include "../system/affine_vector_field.h"

namespace Ariadne {
  namespace System {

    template<typename R>
    AffineVectorField<R>::~AffineVectorField() 
    {
    }
    
    template<typename R>
    LinearAlgebra::Vector<R>
    AffineVectorField<R>::apply(const Geometry::Point<R>& s) const 
    { 
      return this->_A*s.position_vector()+this->_b; 
    }
    
    template<typename R>
    LinearAlgebra::IntervalVector<R> 
    AffineVectorField<R>::apply(const Geometry::Rectangle<R>& r) const 
    {
      LinearAlgebra::IntervalVector<R> iv(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        iv(i)=r.interval(i);
      }
      return LinearAlgebra::IntervalVector<R>(this->_A*iv)+(this->_b);
    }
  
    template<typename R>
    LinearAlgebra::Matrix<R>
    AffineVectorField<R>::derivative(const state_type& x) const 
    { 
      return this->_A; 
    }
    
    template<typename R>
    LinearAlgebra::IntervalMatrix<R> 
    AffineVectorField<R>::derivative(const Rectangle& r) const { 
    return LinearAlgebra::IntervalMatrix<R>(this->_A);
    }
    
  }
}


namespace Ariadne {
  namespace LinearAlgebra {

    template <typename R>
    Matrix<R>
    exp_Ah_approx(const Matrix<R>& A, 
                  const R& h, 
                  const R& e) 
    {
      Matrix<R> result=identity_matrix<R>(A.size1());
      
      R norm_Ah=h*norm(A);
      Matrix<R> AhpowNdivfN=result;
      uint n=0;
      while(norm(AhpowNdivfN)*n >= e*(n-norm_Ah)) {
        ++n;
        AhpowNdivfN=(h/n)*(AhpowNdivfN*A);
        result=result+AhpowNdivfN;
      }
      return result;
    }
    
    template <typename R> 
    Matrix<R> 
    exp_Ah_sub_id_div_A_approx(const Matrix<R>& A, 
                               const R& h, 
                               const R& e)
    {
      Matrix<R> result=h*identity_matrix<R>(A.size1());
      
      R norm_Ah=h*norm(A);
      Matrix<R> AhpowNdivfN=result;
      uint n=0;
      while(norm(AhpowNdivfN)*n >= e*(n-norm_Ah)) {
        ++n;
        AhpowNdivfN=(h/(n+1))*(AhpowNdivfN*A);
        result=result+AhpowNdivfN;
      }
      return result;
      return A;
    }      

  }
}
