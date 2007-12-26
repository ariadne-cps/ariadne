/***************************************************************************
 *            zonotope.inline.h
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
 
#include "numeric/arithmetic.h"

#include "geometry/rectangle_expression.h"
#include "geometry/box.h"


namespace Ariadne {

template<class R> template<class RR> inline
Geometry::Zonotope<R,Geometry::ExactTag>::Zonotope(const Zonotope<RR,ExactTag>& z) 
  : _centre(z.centre()), _generators(z.generators())
{
}


template<class R, class Tag> inline      
tribool
Geometry::empty(const Zonotope<R,Tag>& z) 
{
  return false;
}

template<class R, class Tag> inline      
tribool
Geometry::bounded(const Zonotope<R,Tag>& z) 
{
  return true;
}

template<class R, class Tag> inline      
tribool
Geometry::subset(const Box<R>& r, const Zonotope<R,Tag>& z) 
{
  return superset(z,r);
}


template<class R, class Tag> inline      
R 
Geometry::radius(const Zonotope<R,Tag>& z) 
{
  return radius(z.centre()+z.generators()*z.domain().position_vectors());
}


template<class R, class Tag> inline      
tribool
Geometry::disjoint(const Box<R>& r, Zonotope<R,Tag>& z) 
{
  return disjoint(z,r);
}

template<class R, class Tag> inline      
tribool
Geometry::subset(const Box<R>& r, Zonotope<R,Tag>& z) 
{
  return superset(z,r);
}








} // namespace Ariadne
