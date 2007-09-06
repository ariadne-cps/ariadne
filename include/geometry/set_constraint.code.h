/***************************************************************************
 *            set_constraint.code.h
 *
 *  Copyright  2007  Pieter Collins
 *  pieter.collins@cwi.nl
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
 
#include "set_constraint.h"

namespace Ariadne {

template<class R>
Geometry::SetConstraint<R>::SetConstraint(const Geometry::SetInterface<R>& s, bool i)
  : _set_ptr(s.clone()), _inside(i)
{
}

template<class R>
Geometry::SetConstraint<R>::~SetConstraint() 
{
}

template<class R>
Geometry::SetConstraint<R>* 
Geometry::SetConstraint<R>::clone() const 
{
  return new SetConstraint<R>(*this->_set_ptr,this->_inside);
}


template<class R>
dimension_type 
Geometry::SetConstraint<R>::dimension() const 
{
  return this->_set_ptr->dimension();
}

template<class R>
smoothness_type 
Geometry::SetConstraint<R>::smoothness() const
{
  return 0u;
}

template<class R>
Geometry::Comparison
Geometry::SetConstraint<R>::comparison() const
{
  return this->_inside==true ? less : greater;
}

template<class R>
std::ostream& 
Geometry::SetConstraint<R>::write(std::ostream& os) const
{
  return os << "SetConstraint( set=" << *this->_set_ptr << " )";
}

template<class R>
typename Geometry::SetConstraint<R>::A 
Geometry::SetConstraint<R>::value(const Point<A>& pt) const
{
  Rectangle<R> r(pt);
  tribool inside=this->_set_ptr->superset(r);
  tribool outside=this->_set_ptr->disjoint(r);
  if(inside==true) {
    return -1; 
  } else if(outside==true) {
    return +1;
  } else { 
    return 0;
  }
}

template<class R>
LinearAlgebra::Vector<typename Geometry::SetConstraint<R>::A> 
Geometry::SetConstraint<R>::gradient(const Point<A>& pt) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
const Geometry::SetInterface<R>& 
Geometry::SetConstraint<R>::set() const
{
  return *this->_set_ptr;
}

template<class R>
void 
instantiate()
{
}



} // namespace Ariadne
