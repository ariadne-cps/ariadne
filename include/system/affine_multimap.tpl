/***************************************************************************
 *            affine_multimap.tpl
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
 

#include "../system/affine_map.h"
#include "../system/affine_multimap.h"

namespace Ariadne {
  namespace System {

    template <typename R, template<typename> class BS>
    BS<R>
    AffineMultiMap<R,BS>::operator() (const Geometry::Point<R>& s) const
    {
      using namespace Ariadne::LinearAlgebra;

      const Matrix_type& A=this->A();
      const Vector_type& b=this->b();
      if (A.size2()!=s.dimension()) {
        throw std::domain_error("AffineMultiMap<R,BS>::operator() (const Point& s): the map does not have the same dimension of the point.");
      }
      const Vector_type& pv=s.position_vector();
      Vector_type v=prod(A,pv);
      
      AffineMap<R> map(identity_matrix<R>(s.dimension()),v+b); 
      return map(this->_S);
    }
    
    /* TO IMPROVE */
    template <typename R, template<typename> class BS>
    Geometry::ListSet<R,Geometry::Zonotope>
    AffineMultiMap<R,BS>::operator() (const Geometry::GridMaskSet<R>& cms) const
    {
      using namespace Ariadne::Geometry;

      ListSet<R,Rectangle> lrs=ListSet<R,Rectangle>(cms);

      ListSet<R,Geometry::Zonotope> output(cms.dimension());

      for (size_t i=0; i< lrs.size(); i++) {
	Zonotope<R> p=(*this)(Zonotope<R>(lrs[i]));

	output.push_back(p);
      }
      
      return output;
    }  
  }
}
