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
 *  \brief Images and preimages of boxes in Euclidean space.
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
  const BoxType& domain() const { return this->_domain; }
  const FunctionInterface& function() const { return *this->_function_ptr; }
  bool operator==(const ImageSet& ims) const { return this->_domain==ims._domain && this->_function_ptr==ims._function_ptr; }
  ImageSet* clone() const;
  uint dimension() const;
  tribool disjoint(const Vector<Interval>&) const;
  tribool intersects(const Vector<Interval>&) const;
  tribool subset(const Vector<Interval>&) const;
  Vector<Interval> bounding_box() const;
  std::ostream& write(std::ostream&) const;
};

template<class Mdl>
class ModelSet
  : public CompactSetInterface
{
  Mdl _model;
 public:
  typedef Vector<Interval> BoxType;
  ModelSet(const Mdl& mdl) : _model(mdl) { }
  const BoxType& domain() const { return this->_model.domain(); }
  ModelSet* clone() const { return new ModelSet<Mdl>(*this); }
  uint dimension() const { return this->_model.result_size(); }
  tribool disjoint(const Vector<Interval>& bx) const { 
    return Ariadne::disjoint(this->_model,bx); }
  tribool subset(const Vector<Interval>& bx) const { 
    return Ariadne::subset(this->_model.range(),bx) || indeterminate; }
  Vector<Interval> bounding_box() const { 
    return this->_model.range(); }
  std::ostream& write(std::ostream& os) const {
    return os << "ModelSet( " << this->_model << ")"; }
};
 


class ConstraintSet
  : public RegularSetInterface
{
  Vector<Interval> _codomain;
  boost::shared_ptr<FunctionInterface> _function_ptr;
 public:
  typedef Vector<Interval> BoxType;
  ConstraintSet(const BoxType& codom, const FunctionInterface& fn);
  const BoxType& codomain() const { return this->_codomain; }
  const FunctionInterface& function() const { return *this->_function_ptr; };

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
