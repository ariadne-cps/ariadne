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

#include "real_typedef.h"

namespace Ariadne {
  namespace Geometry {

    template class Rectangle<Real>; 

    template Rectangle<Real> regular_intersection(const Rectangle<Real>&, 
                                         const Rectangle<Real>&);

    template Rectangle<Real> minkowski_sum(const Rectangle<Real>& A, 
                                  const Rectangle<Real>& B);

    template Rectangle<Real> minkowski_difference(const Rectangle<Real>& A, 
                                         const Rectangle<Real>& B);


    template Rectangle<Real> scale(const Rectangle<Real>& s, const Real& factor);

    template std::ostream& operator<<(std::ostream&, const Rectangle<Real>&);
    template std::istream& operator>>(std::istream&, Rectangle<Real>&);

    template std::ostream& operator<<(std::ostream&, const Rectangle<double>&);
  }
}
