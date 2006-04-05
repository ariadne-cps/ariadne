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
 
#include "../evaluation/affine_vector_field.h"

namespace Ariadne {
  namespace Evaluation {

    template<typename R>
    AffineVectorField<R>::~AffineVectorField() 
    {
    }
    
    template<typename R>
    LinearAlgebra::vector<R>
    AffineVectorField<R>::apply(const Geometry::Point<R>& s) const 
    { 
      return this->_A*s.position_vector()+this->_b; 
    }
    
    template<typename R>
    LinearAlgebra::interval_vector<R> 
    AffineVectorField<R>::apply(const Geometry::Rectangle<R>& r) const 
    {
      LinearAlgebra::interval_vector<R> iv(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        iv(i)=r.interval(i);
      }
      return LinearAlgebra::interval_vector<R>(this->_A*iv)+(this->_b);
    }
  
    template<typename R>
    LinearAlgebra::matrix<R>
    AffineVectorField<R>::derivative(const State& x) const 
    { 
      return this->_A; 
    }
    
    template<typename R>
    LinearAlgebra::interval_matrix<R> 
    AffineVectorField<R>::derivative(const Rectangle& r) const { 
      LinearAlgebra::interval_matrix<R> result(this->dimension(),this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        for(dimension_type j=0; j!=this->dimension(); ++j) {
          result(i,j)=this->_A(i,j);
        }
      }
      return result;
    }
    
  }
}


namespace Ariadne {
  namespace LinearAlgebra {

    template <typename R>
    matrix<R>
    exp_Ah_approx(const matrix<R>& A, 
                  const R& h, 
                  const R& e) 
    {
      matrix<R> result=identity_matrix<R>(A.size1());
      
      R norm_Ah=h*norm(A);
      matrix<R> AhpowNdivfN=result;
      uint n=0;
      while(norm(AhpowNdivfN)*n >= e*(n-norm_Ah)) {
        ++n;
        AhpowNdivfN=(h/n)*(AhpowNdivfN*A);
        result=result+AhpowNdivfN;
      }
      return result;
    }
    
    template <typename R> 
    matrix<R> 
    exp_Ah_sub_id_div_A_approx(const matrix<R>& A, 
                               const R& h, 
                               const R& e)
    {
      matrix<R> result=h*identity_matrix<R>(A.size1());
      
      R norm_Ah=h*norm(A);
      matrix<R> AhpowNdivfN=result;
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
