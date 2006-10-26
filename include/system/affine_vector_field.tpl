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
 
#include "affine_vector_field.h"

#include "../numeric/numerical_traits.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/zonotope.h"
#include "../geometry/parallelotope.h"
#include "../geometry/simplex.h"
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace System {

    template<class R>
    AffineVectorField<R>::~AffineVectorField() 
    {
    }
    
    template<class R>
    LinearAlgebra::Vector<typename AffineVectorField<R>::F>
    AffineVectorField<R>::image(const Geometry::Point<R>& s) const 
    { 
      return LinearAlgebra::Vector<F>(this->_A*LinearAlgebra::Vector<F>(s.position_vector())+this->_b); 
    }
    
    template<class R>
    LinearAlgebra::Vector< Interval<R> > 
    AffineVectorField<R>::image(const Geometry::Rectangle<R>& r) const 
    {
      LinearAlgebra::Vector< Interval<R> > iv=r.position_vectors();
      return LinearAlgebra::Vector< Interval<R> >(this->_A*iv)+(this->_b);
    }
  
    template<class R>
    LinearAlgebra::Matrix<typename AffineVectorField<R>::F>
    AffineVectorField<R>::derivative(const Geometry::Point<R>& x) const 
    { 
      return this->_A; 
    }
    
    template<class R>
    LinearAlgebra::Matrix< Interval<R> > 
    AffineVectorField<R>::derivative(const Geometry::Rectangle<R>& r) const 
    { 
      return LinearAlgebra::Matrix< Interval<R> >(this->_A);
    }
    
    template<class R> 
    std::ostream& 
    operator<<(std::ostream& os, const AffineVectorField<R>& vf)
    {
      return os << "AffineVectorField(\n  matrix=" << vf.A() << ",\n"
                << "  vector=" << vf.b() << "\n)\n";
    }

    
    
  }
}


namespace Ariadne {
  namespace LinearAlgebra {

    template<class R>
    Matrix<typename Numeric::traits<R>::arithmetic_type>
    exp_Ah_approx(const Matrix<R>& A, 
                  const R& h, 
                  const R& e) 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      Matrix<F> result=identity_matrix<R>(A.number_of_rows());
      
      F norm_Ah=F(h)*norm(A);
      Matrix<F> AhpowNdivfN=result;
      uint n=0;
      while( (norm(AhpowNdivfN)*static_cast<F>(n)) >= (e*(static_cast<R>(n)-norm_Ah)) ) {
        ++n;
        AhpowNdivfN=(F(h)/R(n))*(AhpowNdivfN*A);
        result=result+AhpowNdivfN;
      }
      return result;
    }
    
    template<class R> 
    Matrix<typename Numeric::traits<R>::arithmetic_type> 
    exp_Ah_sub_id_div_A_approx(const Matrix<R>& A, 
                               const R& h, 
                               const R& e)
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      Matrix<F> result=static_cast<F>(h)*identity_matrix<R>(A.number_of_rows());
      
      F norm_Ah=F(h)*norm(A);
      Matrix<F> AhpowNdivfN=result;
      uint n=0;
      while(norm(AhpowNdivfN)*static_cast<R>(n) >= e*(static_cast<R>(n)-norm_Ah)) {
        ++n;
        AhpowNdivfN=(F(h)/R(n+1))*(AhpowNdivfN*A);
        result=result+AhpowNdivfN;
      }
      return result;
      return A;
    }      

  }
}
