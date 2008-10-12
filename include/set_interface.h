/***************************************************************************
 *            set_interface.h
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
 
/*! \file set_interface.h
 *  \brief Interfaces for open, closed, overt and compact subsets of Euclidean space.
 */

#ifndef ARIADNE_SET_INTERFACE_H
#define ARIADNE_SET_INTERFACE_H

#include <iosfwd>

#include "tribool.h"

namespace Ariadne {

class Interval;
template<class X> class Vector;

class SetInterfaceBase 
{
 public:
  virtual ~SetInterfaceBase() { };
  virtual SetInterfaceBase* clone() const = 0;
  virtual uint dimension() const = 0;
  virtual std::ostream& write(std::ostream& os) const = 0;
};

class OvertSetInterface 
  : public virtual SetInterfaceBase 
{
 public:
  virtual OvertSetInterface* clone() const = 0;
  virtual tribool intersects(const Vector<Interval>& bx) const = 0;
};

class OpenSetInterface 
  : public virtual OvertSetInterface 
{
 public:
  virtual OpenSetInterface* clone() const = 0;
  virtual tribool superset(const Vector<Interval>& bx) const = 0;
};

class ClosedSetInterface
  : public virtual SetInterfaceBase {
 public:
  virtual ClosedSetInterface* clone() const = 0;
  virtual tribool disjoint(const Vector<Interval>& bx) const = 0;
};

class CompactSetInterface
  : public virtual ClosedSetInterface {
 public:
  virtual CompactSetInterface* clone() const = 0;
  virtual tribool subset(const Vector<Interval>& bx) const = 0;
  virtual Vector<Interval> bounding_box() const = 0;
};

class RegularSetInterface 
  : public virtual OpenSetInterface,
    public virtual ClosedSetInterface
{
  virtual RegularSetInterface* clone() const = 0;
};

class LocatedSetInterface 
  : public virtual OvertSetInterface,
    public virtual CompactSetInterface
{
  virtual LocatedSetInterface* clone() const = 0;
};

class SetInterface 
  : public virtual RegularSetInterface,
    public virtual LocatedSetInterface
{
  virtual SetInterface* clone() const = 0;
};
    

inline std::ostream& operator<<(std::ostream& os, const SetInterfaceBase& s) {
  return s.write(os); 
}


} // namespace Ariadne


#endif // ARIADNE_SET_INTERFACE
