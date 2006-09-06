/***************************************************************************
 *            rectangle.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include "geometry/rectangle.h"
#include "geometry/rectangle.tpl"

#include "numeric/float64.h"
#include "numeric/rational.h"
#include "real_typedef.h"


namespace Ariadne {
  namespace Geometry {
    template class Rectangle<Real>; 
    template bool subset(const Rectangle<Real>&, 
                         const ListSet<Real,Geometry::Rectangle>&);

    #ifndef RATIONAL_REAL
    template class Rectangle<Rational>; 
    #endif
      
    #ifndef DOUBLE_REAL  
    template class Rectangle<Float64>; 
    #endif
    
  }
}
