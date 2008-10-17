/***************************************************************************
 *            set.cc
 *
 *  Copyright 2008  Pieter Collins
 * 
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
 
#include "set.h"
#include "macros.h"

namespace Ariadne {

 
ImageSet::ImageSet(const BoxType& dom, const FunctionInterface& fn)
  : _domain(dom), _function(fn.clone()) 
{ 
  ARIADNE_ASSERT(dom.size()==fn.argument_size()); 
}


ImageSet*
ImageSet::clone() const 
{ 
  return new ImageSet(*this);
}


uint
ImageSet::dimension() const 
{ 
  return this->_function->result_size();
}


tribool
ImageSet::disjoint(const Vector<Interval>& bx) const 
{ 
  ARIADNE_NOT_IMPLEMENTED;
}


tribool
ImageSet::intersects(const Vector<Interval>& bx) const 
{ 
  ARIADNE_NOT_IMPLEMENTED;
}


tribool
ImageSet::subset(const Vector<Interval>& bx) const 
{ 
  ARIADNE_NOT_IMPLEMENTED;
}


Vector<Interval>
ImageSet::bounding_box() const 
{ 
  ARIADNE_NOT_IMPLEMENTED;
}


std::ostream&
ImageSet::write(std::ostream& os) const 
{ 
  ARIADNE_NOT_IMPLEMENTED;
}





ConstraintSet::ConstraintSet(const BoxType& codom, const FunctionInterface& fn) 
  : _codomain(codom), _function(fn.clone()) 
{ 
  ARIADNE_ASSERT(codom.size()==fn.result_size());
}

ConstraintSet*
ConstraintSet::clone() const 
{ 
  return new ConstraintSet(*this);
}


uint
ConstraintSet::dimension() const 
{ 
  return this->_function->argument_size();
}


tribool
ConstraintSet::disjoint(const Vector<Interval>& bx) const 
{ 
  ARIADNE_NOT_IMPLEMENTED;
}


tribool
ConstraintSet::intersects(const Vector<Interval>& bx) const 
{ 
  ARIADNE_NOT_IMPLEMENTED;
}


tribool
ConstraintSet::superset(const Vector<Interval>& bx) const 
{ 
  ARIADNE_NOT_IMPLEMENTED;
}


Vector<Interval>
ConstraintSet::bounding_box() const 
{ 
  ARIADNE_NOT_IMPLEMENTED;
}


std::ostream&
ConstraintSet::write(std::ostream& os) const 
{ 
  ARIADNE_NOT_IMPLEMENTED;
}



} // namespace Ariadne
