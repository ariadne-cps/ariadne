/***************************************************************************
 *            vector_field.tpl
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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
 
#include "vector_field.h"

#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_vector.h"
#include "../linear_algebra/interval_matrix.h"
#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"



namespace Ariadne {
  namespace System {

    template<typename R>
    VectorField<R>::~VectorField() 
    {
    }
     
    template<typename R>
    LinearAlgebra::Vector<R> 
    VectorField<R>::operator() (const Geometry::Point<R>& x) const
    {
      throw std::invalid_argument(this->name()+"::operator() (Point) not implemented."); 
    }
    
    template<typename R>
    LinearAlgebra::IntervalVector<R>
    VectorField<R>::operator() (const Geometry::Rectangle<R>& x) const
    {
      throw std::invalid_argument(this->name()+"::operator() (Rectangle) not implemented."); 
    }
    
    template<typename R>
    LinearAlgebra::Matrix<R> 
    VectorField<R>::derivative(const Geometry::Point<R>& x) const 
    {
      throw std::invalid_argument(this->name()+"::derivative(Point) not implemented."); 
    }
    
    template<typename R>
    LinearAlgebra::IntervalMatrix<R> 
    VectorField<R>::derivative(const Geometry::Rectangle<R>& A) const 
    {
      throw std::invalid_argument(this->name()+"::derivative(Rectangle) not implemented."); 
    }
    
  }
}
