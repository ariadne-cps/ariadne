/***************************************************************************
 *            affine_map.code.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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
 

#include "affine_map.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/polytope.h"




namespace Ariadne {
  namespace System {

    template<class R>
    Geometry::Point<typename AffineMap<R>::F>
    AffineMap<R>::image(const Geometry::Point<F>& pt) const
    {
      check_argument_dimension(*this,pt,__PRETTY_FUNCTION__);
      LinearAlgebra::Vector<F> image(this->A()*LinearAlgebra::Vector<F>(pt.position_vector())+this->b());
      return Geometry::Point<F>(image);
    }
    
    
    template<class R>
    Geometry::Zonotope<typename AffineMap<R>::F>
    AffineMap<R>::image(const Geometry::Zonotope<F>& z) const
    {
      //std::cerr << __PRETTY_FUNCTION__ << std::endl;

      check_argument_dimension(*this,z,__PRETTY_FUNCTION__);
      Geometry::Point<F> c=z.centre();
      LinearAlgebra::Matrix<F> G=z.generators();
      const LinearAlgebra::Matrix<R>& A=this->A();
      const LinearAlgebra::Vector<R>& b=this->b();
      
      return Geometry::Zonotope<F>(
        Geometry::Point<F>(A*c.position_vector()+b),A*G
      );
    }   
    
    template<class R>
    Geometry::Polytope<typename AffineMap<R>::F>
    AffineMap<R>::image(const Geometry::Polytope<F>& p) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }   
    
    template<class R>
    LinearAlgebra::Matrix<typename AffineMap<R>::F>
    AffineMap<R>::jacobian(const Geometry::Point<F>& pt) const
    {
      return this->_A;
    }
    
    template<class R>
    std::ostream& 
    AffineMap<R>::write(std::ostream& os) const
    {
      return os << "AffineMap( A=" << this->A()
                << ", b=" << this->b() << " )";
    }
    


  }
}
