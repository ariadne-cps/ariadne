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
#include "geometry/box.h"
#include "geometry/euclidean_space.h"

namespace Ariadne {

template<class R>
SetConstraint<R>::SetConstraint(const SetInterface< Box<R> >& s, bool i)
  : _set_ptr(s.clone()), _inside(i)
{
}

template<class R>
SetConstraint<R>::~SetConstraint() 
{
}

template<class R>
SetConstraint<R>* 
SetConstraint<R>::clone() const 
{
  return new SetConstraint<R>(*this->_set_ptr,this->_inside);
}


template<class R>
EuclideanSpace
SetConstraint<R>::space() const 
{
  return this->_set_ptr->space();
}


template<class R>
dimension_type
SetConstraint<R>::dimension() const 
{
  return this->_set_ptr->space().dimension();
}

template<class R>
smoothness_type 
SetConstraint<R>::smoothness() const
{
  return 0u;
}

template<class R>
Comparison
SetConstraint<R>::comparison() const
{
  return this->_inside==true ? less : greater;
}

template<class R>
std::ostream& 
SetConstraint<R>::write(std::ostream& os) const
{
  return os << "SetConstraint( set=" << *this->_set_ptr << " )";
}

template<class R>
typename SetConstraint<R>::A 
SetConstraint<R>::value(const Point<A>& pt) const
{
  Box<R> r(pt);
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
Vector<typename SetConstraint<R>::A> 
SetConstraint<R>::gradient(const Point<A>& pt) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
const SetInterface< Box<R> >& 
SetConstraint<R>::set() const
{
  return *this->_set_ptr;
}

template<class R>
bool
SetConstraint<R>::inside() const
{
  return this->_inside;
}

template<class R>
void 
instantiate()
{
}



} // namespace Ariadne
