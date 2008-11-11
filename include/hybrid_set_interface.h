/***************************************************************************
 *            hybrid_set_interface.h
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

namespace Ariadne {

class Interval;
class Box;

class HybridSpace;
typedef int DiscreteState;

typedef std::pair<DiscreteState,Box> HybridBox;
typedef std::map<DiscreteState,Box> HybridBoxes;

class HybridSetInterfaceBase 
{
  public:
    virtual ~HybridSetInterfaceBase() { };
    virtual HybridSetInterfaceBase* clone() const = 0;
    virtual HybridSpace space() const = 0;
    virtual std::ostream& write(std::ostream& os) const = 0;
};

class HybridOvertSetInterface 
    : public virtual HybridSetInterfaceBase 
{
  public:
    virtual HybridOvertSetInterface* clone() const = 0;
    virtual tribool intersects(const HybridBox& bx) const = 0;
    virtual OvertSetInterface const& operator[](DiscreteState) const = 0;
};

class HybridOpenSetInterface 
    : public virtual HybridOvertSetInterface 
{
  public:
    virtual HybridOpenSetInterface* clone() const = 0;
    virtual tribool superset(const HybridBox& bx) const = 0;
    virtual OpenSetInterface const& operator[](DiscreteState) const = 0;
};

class HybridClosedSetInterface
    : public virtual HybridSetInterfaceBase {
  public:
    virtual HybridClosedSetInterface* clone() const = 0;
    virtual tribool disjoint(const HybridBox& bx) const = 0;
    virtual ClosedSetInterface const& operator[](DiscreteState) const = 0;
};

class HybridCompactSetInterface
    : public virtual HybridClosedSetInterface {
  public:
    virtual HybridCompactSetInterface* clone() const = 0;
    virtual tribool subset(const HybridBoxes& bx) const = 0;
    virtual HybridBoxes bounding_box() const = 0;
    virtual CompactSetInterface const& operator[](DiscreteState) const = 0;
};

class HybridRegularSetInterface 
    : public virtual HybridOpenSetInterface,
      public virtual HybridClosedSetInterface
{
    virtual HybridRegularSetInterface* clone() const = 0;
    virtual RegularSetInterface const& operator[](DiscreteState) const = 0;
};

class HybridLocatedSetInterface 
    : public virtual HybridOvertSetInterface,
      public virtual HybridCompactSetInterface
{
    virtual HybridLocatedSetInterface* clone() const = 0;
    virtual LocatedSetInterface const& operator[](DiscreteState) const = 0;
};

class HybridSetInterface 
    : public virtual HybridRegularSetInterface,
      public virtual HybridLocatedSetInterface
{
    virtual HybridSetInterface* clone() const = 0;
    virtual SetInterface const& operator[](DiscreteState) const = 0;
};
    

inline std::ostream& operator<<(std::ostream& os, const HybridSetInterfaceBase& s) {
    return s.write(os); 
}


} // namespace Ariadne


#endif // ARIADNE_SET_INTERFACE
