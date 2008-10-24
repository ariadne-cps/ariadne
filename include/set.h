/***************************************************************************
 *            set.h
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
 
/*! \file set.h
 *  \brief Image and preimages of boxes in Euclidean space.
 */

#ifndef ARIADNE_SET_H
#define ARIADNE_SET_H

#include <iosfwd>

#include "boost/shared_ptr.hpp"
#include "macros.h"
#include "numeric.h"
#include "vector.h"
#include "set_interface.h"
#include "function_interface.h"


namespace Ariadne {

class ImageSet
  : public LocatedSetInterface
{
  Vector<Interval> _domain;
  boost::shared_ptr<FunctionInterface> _function_ptr;
 public:
  typedef Vector<Interval> BoxType;
  ImageSet();
  ImageSet(const BoxType& dom);
  ImageSet(const BoxType& dom, const FunctionInterface& fn);
  const BoxType& domain() { return this->_domain; }
  const FunctionInterface& function() { return *this->_function_ptr; }
  bool operator==(const ImageSet& ims) const { return this->_domain==ims._domain && this->_function_ptr==ims._function_ptr; }
  ImageSet* clone() const;
  uint dimension() const;
  tribool disjoint(const Vector<Interval>&) const;
  tribool intersects(const Vector<Interval>&) const;
  tribool subset(const Vector<Interval>&) const;
  Vector<Interval> bounding_box() const;
  std::ostream& write(std::ostream&) const;
};

class ConstraintSet
  : public RegularSetInterface
{
  Vector<Interval> _codomain;
  boost::shared_ptr<FunctionInterface> _function_ptr;
 public:
  typedef Vector<Interval> BoxType;
  ConstraintSet(const BoxType& codom, const FunctionInterface& fn);
  const BoxType& codomain() { return this->_codomain; }
  const FunctionInterface& function() { return *this->_function_ptr; };

  ConstraintSet* clone() const;
  uint dimension() const;
  tribool disjoint(const Vector<Interval>&) const;
  tribool intersects(const Vector<Interval>&) const;
  tribool superset(const Vector<Interval>&) const;
  Vector<Interval> bounding_box() const;
  std::ostream& write(std::ostream&) const;
};

} // namespace Ariadne 

#endif
