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

class BoundedSetInterface;
class OpenSetInterface;
class ClosedSetInterface;
class OvertSetInterface;
class CompactSetInterface;
class RegularSetInterface;
class LocatedSetInterface;

//! \brief Base class for sets described by predicates involving boxes.
class SetInterfaceBase 
{
 public:
  //! \brief Virtual destructor.
  virtual ~SetInterfaceBase() { };
  //! \brief Construct a dynamically-allocated copy.
  virtual SetInterfaceBase* clone() const = 0;
  //! \brief The dimension of the set.
  virtual uint dimension() const = 0;
  //! \brief Write to an output stream.
  virtual std::ostream& write(std::ostream& os) const = 0;
};

//! \brief Interface for bounded sets. 
class BoundedSetInterface
  : public virtual SetInterfaceBase {
 public:
  virtual BoundedSetInterface* clone() const = 0;
  //! \brief Tests if the set is a subset of \a bx.
  virtual tribool subset(const Vector<Interval>& bx) const = 0;
  //! \brief Returns a bounding box for the set.
  //! If the set is empty, then the first component of the result should be empty. 
  virtual Vector<Interval> bounding_box() const = 0;
};


//! \brief Interface for overt sets. 
class OvertSetInterface 
  : public virtual SetInterfaceBase 
{
 public:
  virtual OvertSetInterface* clone() const = 0;
  //! \brief Tests if the set intersects \a bx.
  virtual tribool intersects(const Vector<Interval>& bx) const = 0;
  //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
  friend tribool overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops, const Float& eps);
};

//! \brief Interface for open sets. 
class OpenSetInterface 
  : public virtual OvertSetInterface 
{
 public:
  virtual OpenSetInterface* clone() const = 0;
  //! \brief Tests if the set is a superset of \a bx.
  virtual tribool superset(const Vector<Interval>& bx) const = 0;
  //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
  friend tribool overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops, const Float& eps);
  //! \brief Tests if \a ls is a subset of \a rs, to a tolerance of \a eps.
  friend tribool subset(const CompactSetInterface& ls, const OpenSetInterface& rs, const Float& eps);
};

//! \brief Interface for closed sets. 
class ClosedSetInterface
  : public virtual SetInterfaceBase 
{
 public:
  virtual ClosedSetInterface* clone() const = 0;
  //! \brief Tests if the set is disjoint from \a bx.
  virtual tribool disjoint(const Vector<Interval>& bx) const = 0;
  //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
  friend tribool disjoint(const CompactSetInterface& cps, const ClosedSetInterface& cls, const Float& eps);
};

//! \brief Interface for compact (closed and bounded) sets. 
class CompactSetInterface
  : public virtual BoundedSetInterface,
    public virtual ClosedSetInterface 
{
 public:
  virtual CompactSetInterface* clone() const = 0;
  //virtual tribool empty() const = 0;
  //! \brief Tests if \a ls is a subset of \a rs, to a tolerance of \a eps.
  friend tribool subset(const CompactSetInterface& ls, const OpenSetInterface& rs, const Float& eps);
  //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
  friend tribool disjoint(const CompactSetInterface& cps, const ClosedSetInterface& cls, const Float& eps);

};

//! \brief Interface for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure. 
class RegularSetInterface 
  : public virtual OpenSetInterface,
    public virtual ClosedSetInterface
{
  virtual RegularSetInterface* clone() const = 0;
  //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
  friend tribool overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);

  //! \brief Tests if \a ls is a subset of \a rs, to a tolerance of \a eps.
  friend tribool subset(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);
  //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
  friend tribool disjoint(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);
};


//! \brief Interface for located (overt and compact) sets. 
class LocatedSetInterface 
  : public virtual OvertSetInterface,
    public virtual CompactSetInterface
{
  virtual LocatedSetInterface* clone() const = 0;
  //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
  friend tribool overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);

  //! \brief Tests if \a ls is a subset of \a rs, to a tolerance of \a eps.
  friend tribool subset(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);
  //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
  friend tribool disjoint(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);
};

//! \brief Complete set interface for bounded regular sets.
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
