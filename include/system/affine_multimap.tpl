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
#include "affine_multimap.h"

namespace Ariadne {
  namespace System {

    template <typename R, template<typename> class BS>
    BS<R>
    AffineMultiMap<R,BS>::operator() (const Geometry::Point<R>& pt) const
    {
      using namespace Ariadne::LinearAlgebra;

      if (this->argument_dimension()!=pt.dimension()) {
        throw std::domain_error("AffineMultiMap<R,BS>::operator() (const Point&): the map does not have the same dimension of the point.");
      }
      IntervalVector<R> iv=this->A()*IntervalVector<R>(pt.position_vector());
      return minkowski_sum(this->S(),BS<R>(Geometry::Rectangle<R>(iv)));
    }
    
    template <typename R, template<typename> class BS>
    BS<R>
    AffineMultiMap<R,BS>::operator() (const BS<R>& bs) const
    {
      using namespace Ariadne::LinearAlgebra;
      using namespace Ariadne::Geometry;

      if (this->argument_dimension()!=bs.dimension()) {
        throw std::domain_error("AffineMultiMap<R,BS>::operator() (const Point&): the map does not have the same dimension of the point.");
      }
      return minkowski_sum(AffineMap<R>(this->A())(bs),this->S());
    }
    
    /* TO IMPROVE */
    template <typename R, template<typename> class BS>
    Geometry::ListSet<R,BS>
    AffineMultiMap<R,BS>::operator() (const Geometry::GridMaskSet<R>& gms) const
    {
      using namespace Ariadne::Geometry;

      ListSet<R,Rectangle> lrs=ListSet<R,Rectangle>(gms);
      ListSet<R,BS> output(gms.dimension());

      for (size_type i=0; i< lrs.size(); ++i) {
        BS<R> rz(lrs[i]);
        BS<R> p=this->operator()(rz);
        output.push_back(p);
      }
      
      return output;
    }  
  }
}
