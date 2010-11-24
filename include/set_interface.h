/***************************************************************************
 *      set_interface.h
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
#include "numeric.h"

namespace Ariadne {

class Box;

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

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for bounded sets.
class BoundedSetInterface
    : public virtual SetInterfaceBase {
  public:
    virtual BoundedSetInterface* clone() const = 0;
    //! \brief Tests if the set is a inside of \a bx.
    virtual tribool inside(const Box& bx) const = 0;
    //! \brief Returns a bounding box for the set.
    //! If the set is empty, then the first component of the result should be empty.
    virtual Box bounding_box() const = 0;
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for overt sets, for which intersection with an open box is verifiable.
class OvertSetInterface
    : public virtual SetInterfaceBase
{
  public:
    virtual OvertSetInterface* clone() const = 0;
    //! \brief Tests if the set overlaps \a bx.
    virtual tribool overlaps(const Box& bx) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend tribool overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops, const Float& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for open sets.
class OpenSetInterface
    : public virtual OvertSetInterface
{
  public:
    virtual OpenSetInterface* clone() const = 0;
    //! \brief Tests if the set is a covers of \a bx.
    virtual tribool covers(const Box& bx) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend tribool overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops, const Float& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend tribool inside(const CompactSetInterface& ls, const OpenSetInterface& rs, const Float& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for closed sets.
class ClosedSetInterface
    : public virtual SetInterfaceBase
{
  public:
    virtual ClosedSetInterface* clone() const = 0;
    //! \brief Tests if the set is disjoint from \a bx.
    virtual tribool disjoint(const Box& bx) const = 0;
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend tribool disjoint(const CompactSetInterface& cps, const ClosedSetInterface& cls, const Float& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for compact (closed and bounded) sets.
class CompactSetInterface
    : public virtual BoundedSetInterface,
      public virtual ClosedSetInterface
{
  public:
    virtual CompactSetInterface* clone() const = 0;
    //virtual tribool empty() const = 0;
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend tribool inside(const CompactSetInterface& ls, const OpenSetInterface& rs, const Float& eps);
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend tribool disjoint(const CompactSetInterface& cps, const ClosedSetInterface& cls, const Float& eps);

};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
class RegularSetInterface
    : public virtual OpenSetInterface,
      public virtual ClosedSetInterface
{
    virtual RegularSetInterface* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend tribool overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);

    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend tribool inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend tribool disjoint(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);
};


class GridTreeSet;


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for located (overt and compact) sets.
class LocatedSetInterface
    : public virtual OvertSetInterface,
      public virtual CompactSetInterface
{
    virtual LocatedSetInterface* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend tribool overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);

    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend tribool inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend tribool disjoint(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
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


//! \brief A Euclidean space \f$\R^d\f$ of dimension \a d.
class EuclideanSpace
{
  public:
    //! \brief The interface satisified by bounded sets in the space.
    typedef BoundedSetInterface BoundedSetInterfaceType;
    //! \brief The interface satisified by overt sets in the space.
    typedef OvertSetInterface OvertSetInterfaceType;
    //! \brief The interface satisified by over sets in the space.
    typedef OpenSetInterface OpenSetInterfaceType;
    //! \brief The interface satisified by closed sets in the space.
    typedef ClosedSetInterface ClosedSetInterfaceType;
    //! \brief The interface satisified by compact sets in the space.
    typedef CompactSetInterface CompactSetInterfaceType;
    //! \brief The interface satisified by regular sets in the space.
    typedef RegularSetInterface RegularSetInterfaceType;
    //! \brief The interface satisified by located sets in the space.
    typedef LocatedSetInterface LocatedSetInterfaceType;
    //! \brief The type of approximations to sets in the space.
    typedef GridTreeSet SetApproximationType;
  public:
    EuclideanSpace(const size_t& d) : _dimension(d) { }
    const size_t& dimension() const { return this->_dimension; }
  private:
    size_t _dimension;
};


} // namespace Ariadne


#endif // ARIADNE_SET_INTERFACE
