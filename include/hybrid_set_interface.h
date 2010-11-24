/***************************************************************************
 *      hybrid_set_interface.h
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

/*! \file hybrid_set_interface.h
 *  \brief Interfaces for open, closed, overt and compact subsets of hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_SET_INTERFACE_H
#define ARIADNE_HYBRID_SET_INTERFACE_H

#include <iosfwd>
#include <utility>
#include <map>


#include "tribool.h"
#include "set_interface.h"
#include "discrete_location.h"

namespace Ariadne {

class Point;
class Box;

class HybridSpace;
template<class BS> class HybridBasicSet;

typedef HybridBasicSet<Point> HybridPoint;
typedef HybridBasicSet<Box> HybridBox;
typedef std::map<DiscreteLocation,Box> HybridBoxes;

//! \brief Base class for sets in a hybrid space.
class HybridSetInterfaceBase
{
  public:
    virtual ~HybridSetInterfaceBase() { };
    virtual HybridSetInterfaceBase* clone() const = 0;
    virtual HybridSpace space() const = 0;
    virtual SetInterfaceBase const& operator[](DiscreteLocation) const = 0;
    virtual std::ostream& write(std::ostream& os) const = 0;
    friend std::ostream& operator<<(std::ostream& os, const HybridSetInterfaceBase& hs);
};

//! \brief Interface for bounded sets in a hybrid space.
class HybridBoundedSetInterface
    : public virtual HybridSetInterfaceBase
{
  public:
    virtual HybridBoundedSetInterface* clone() const = 0;
    virtual tribool inside(const HybridBoxes& bx) const = 0;
    virtual HybridBoxes bounding_box() const = 0;
    virtual BoundedSetInterface const& operator[](DiscreteLocation) const = 0;
};

//! \brief Interface for overt sets in a hybrid space.
class HybridOvertSetInterface
    : public virtual HybridSetInterfaceBase
{
  public:
    virtual HybridOvertSetInterface* clone() const = 0;
    virtual tribool overlaps(const HybridBox& bx) const = 0;
    virtual OvertSetInterface const& operator[](DiscreteLocation) const = 0;
};

//! \brief Interface for open sets in a hybrid space.
class HybridOpenSetInterface
    : public virtual HybridOvertSetInterface
{
  public:
    virtual HybridOpenSetInterface* clone() const = 0;
    virtual tribool covers(const HybridBox& bx) const = 0;
    virtual OpenSetInterface const& operator[](DiscreteLocation) const = 0;
};

//! \brief Interface for closed sets in a hybrid space.
class HybridClosedSetInterface
    : public virtual HybridSetInterfaceBase
{
  public:
    virtual HybridClosedSetInterface* clone() const = 0;
    virtual tribool disjoint(const HybridBox& bx) const = 0;
    virtual ClosedSetInterface const& operator[](DiscreteLocation) const = 0;
};

//! \brief Interface for compact (closed and bounded) sets in a hybrid space.
class HybridCompactSetInterface
    : public virtual HybridBoundedSetInterface
    , public virtual HybridClosedSetInterface
{
  public:
    virtual HybridCompactSetInterface* clone() const = 0;
    virtual CompactSetInterface const& operator[](DiscreteLocation) const = 0;
};

//! \brief Interface for regular (open and closed) sets in a hybrid space.
class HybridRegularSetInterface
    : public virtual HybridOpenSetInterface,
      public virtual HybridClosedSetInterface
{
    virtual HybridRegularSetInterface* clone() const = 0;
    virtual RegularSetInterface const& operator[](DiscreteLocation) const = 0;
};

//! \brief Interface for located (overt and compact) sets in a hybrid space.
class HybridLocatedSetInterface
    : public virtual HybridOvertSetInterface,
      public virtual HybridCompactSetInterface
{
    virtual HybridLocatedSetInterface* clone() const = 0;
    virtual LocatedSetInterface const& operator[](DiscreteLocation) const = 0;
};

//! \brief Complete set interface for bounded regular sets in a hybrid space.
class HybridSetInterface
    : public virtual HybridRegularSetInterface,
      public virtual HybridLocatedSetInterface
{
    virtual HybridSetInterface* clone() const = 0;
    virtual SetInterface const& operator[](DiscreteLocation) const = 0;
};


inline std::ostream& operator<<(std::ostream& os, const HybridSetInterfaceBase& s) {
    return s.write(os);
}


} // namespace Ariadne


#endif // ARIADNE_SET_INTERFACE
