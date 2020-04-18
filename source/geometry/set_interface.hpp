/***************************************************************************
 *            geometry/set_interface.hpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file geometry/set_interface.hpp
 *  \brief Interfaces for open, closed, overt and compact subsets of Euclidean space.
 */

#ifndef ARIADNE_SET_INTERFACE_HPP
#define ARIADNE_SET_INTERFACE_HPP

#include <iosfwd>

#include "../utility/declarations.hpp"
#include "../utility/tribool.hpp"
#include "../utility/writable.hpp"
#include "../numeric/numeric.hpp"

#include "box.decl.hpp"

namespace Ariadne {

template<class X> class Vector;

template<class UB> class Interval;
typedef FloatDPExactInterval ExactIntervalType;

template<class IVL> class Box;
typedef FloatDPExactBox ExactBoxType;
typedef FloatDPUpperBox UpperBoxType;

class BoundedSetInterface;
class OpenSetInterface;
class ClosedSetInterface;
class OvertSetInterface;
class CompactSetInterface;
class RegularSetInterface;
class LocatedSetInterface;
class RegularLocatedSetInterface;

struct Overlap {
    constexpr const char* code() const { return "overlap"; }
    template<class S> decltype(auto) operator()(S const& s, ExactBoxType const& bx, Effort eff) const { return s.overlaps(bx,eff); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2, Effort eff) const { return overlap(s1,s2,eff); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2) const { return overlap(s1,s2); }
};

struct Separated {
    constexpr const char* code() const { return "separated"; }
    template<class S> decltype(auto) operator()(S const& s, ExactBoxType const& bx, Effort eff) const { return s.separated(bx,eff); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2, Effort eff) const { return separated(s1,s2,eff); }
};

struct Inside {
    constexpr const char* code() const { return "inside"; }
    template<class S> decltype(auto) operator()(S const& s, ExactBoxType const& bx, Effort eff) const { return s.inside(bx,eff); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2, Effort eff) const { return inside(s1,s2,eff); }
};

struct Covers {
    constexpr const char* code() const { return "covers"; }
    template<class S> decltype(auto) operator()(S const& s, ExactBoxType const& bx, Effort eff) const { return s.covers(bx,eff); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2, Effort eff) const { return inside(s2,s1,eff); }
};

//! \brief Base class for sets described by predicates involving boxes.
class SetInterfaceBase : public virtual WritableInterface
{
  public:
    //! \brief Virtual destructor.
    virtual ~SetInterfaceBase() = default;
    //! \brief Construct a dynamically-allocated copy.
    virtual SetInterfaceBase* clone() const = 0;
    //! \brief The dimension of the set.
    virtual DimensionType dimension() const = 0;
    //! \brief Write to an output stream.
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for bounded sets.
class BoundedSetInterface
    : public virtual SetInterfaceBase {
  public:
    virtual BoundedSetInterface* clone() const = 0;
    //! \brief Tests if the set is a inside of \a bx.
    //! A set \a A is \em inside \a B if the closure of \a A is a subset of the interior of \a B.
    //! A set \f$A\f$ is \em inside \f$B\f$ if \f$\,\overline{\!A} \subset B^\circ\f$.
    virtual LowerKleenean inside(const ExactBoxType& bx) const = 0;
    virtual ValidatedLowerKleenean inside(const ExactBoxType& bx, Effort eff) const = 0;
    //! \brief Returns a bounding box for the set.
    //! If the set is empty, then the first component of the result should be empty.
    virtual UpperBoxType bounding_box() const = 0;
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for overt sets, for which intersection with an open box is verifiable.
class OvertSetInterface
    : public virtual SetInterfaceBase
{
  public:
    virtual OvertSetInterface* clone() const = 0;
    //! \brief Tests if the set overlaps \a bx.
    //! Sets \a A and \a B \em overlap if the interiors of \a A and \a B intersect.
    //! Sets \f$A\f$ and \f$B\f$ \em overlap if \f$A^\circ \cap B^\circ \neq \emptyset\f$.
    virtual LowerKleenean overlaps(const ExactBoxType& bx) const = 0;
    virtual ValidatedLowerKleenean overlaps(const ExactBoxType& bx, Effort eff) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend LowerKleenean overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops);
    friend ValidatedLowerKleenean overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for open sets.
class OpenSetInterface
    : public virtual OvertSetInterface
{
  public:
    virtual OpenSetInterface* clone() const = 0;
    //! \brief Tests if the set covers of \a bx.
    //! A set \a A \em covers \a B if the interiors of \a A is a superset of the closure of \a B.
    //! A set \f$A\f$ \em covers \f$B\f$ if \f$A^\circ \supset \overline{B}\f$.
    virtual LowerKleenean covers(const ExactBoxType& bx) const = 0;
    virtual ValidatedLowerKleenean covers(const ExactBoxType& bx, Effort eff) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend LowerKleenean overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops);
    friend ValidatedLowerKleenean overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops, const RawFloatDP& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend LowerKleenean inside(const CompactSetInterface& ls, const OpenSetInterface& rs);
    friend ValidatedLowerKleenean inside(const CompactSetInterface& ls, const OpenSetInterface& rs, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for closed sets.
class ClosedSetInterface
    : public virtual SetInterfaceBase
{
  public:
    virtual ClosedSetInterface* clone() const = 0;
    //! \brief Tests if the set is separated from \a bx.
    //! A set \a A is \em separated from \a B if the closures of \a A and \a B are disjoint.
    //! A set \f$A\f$ is \em separated from \f$B\f$ if \f$\,\overline{\!A} \cap \overline{B} = \emptyset\f$.
    virtual LowerKleenean separated(const ExactBoxType& bx) const = 0;
    virtual ValidatedLowerKleenean separated(const ExactBoxType& bx, Effort eff) const = 0;
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend LowerKleenean separated(const CompactSetInterface& cps, const ClosedSetInterface& cls);
    friend ValidatedLowerKleenean separated(const CompactSetInterface& cps, const ClosedSetInterface& cls, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for compact (closed and bounded) sets.
class CompactSetInterface
    : public virtual BoundedSetInterface,
      public virtual ClosedSetInterface
{
  public:
    virtual CompactSetInterface* clone() const = 0;
    //virtual ValidatedSierpinskian empty() const = 0;
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend LowerKleenean inside(const CompactSetInterface& ls, const OpenSetInterface& rs);
    friend ValidatedLowerKleenean inside(const CompactSetInterface& ls, const OpenSetInterface& rs, const RawFloatDP& eps);
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend LowerKleenean separated(const CompactSetInterface& cps, const ClosedSetInterface& cls);
    friend ValidatedLowerKleenean separated(const CompactSetInterface& cps, const ClosedSetInterface& cls, const RawFloatDP& eps);

};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
class RegularSetInterface
    : public virtual OpenSetInterface,
      public virtual ClosedSetInterface
{
    virtual RegularSetInterface* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend Kleenean overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs);
    friend ValidatedKleenean overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend Kleenean inside(const LocatedSetInterface& ls, const RegularSetInterface& rs);
    friend ValidatedKleenean inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend Kleenean separated(const LocatedSetInterface& ls, const RegularSetInterface& rs);
    friend ValidatedKleenean separated(const LocatedSetInterface& ls, const RegularSetInterface& rs, const RawFloatDP& eps);
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for located (overt and compact) sets.
class LocatedSetInterface
    : public virtual OvertSetInterface,
      public virtual CompactSetInterface
{
    virtual LocatedSetInterface* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend Kleenean overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs);
    friend ValidatedKleenean overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend Kleenean inside(const LocatedSetInterface& ls, const RegularSetInterface& rs);
    friend ValidatedKleenean inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend Kleenean separated(const LocatedSetInterface& ls, const RegularSetInterface& rs);
    friend ValidatedKleenean separated(const LocatedSetInterface& ls, const RegularSetInterface& rs, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Complete set interface for bounded regular sets.
class RegularLocatedSetInterface
    : public virtual RegularSetInterface,
      public virtual LocatedSetInterface
{
  public:
    virtual RegularLocatedSetInterface* clone() const = 0;
};

using SetInterface = RegularLocatedSetInterface;

inline OutputStream& operator<<(OutputStream& os, const SetInterfaceBase& s) {
    return s._write(os);
}



class ValidatedBoundedSetInterface;
class ValidatedOpenSetInterface;
class ValidatedClosedSetInterface;
class ValidatedOvertSetInterface;
class ValidatedCompactSetInterface;
class ValidatedRegularSetInterface;
class ValidatedLocatedSetInterface;
class ValidatedRegularLocatedSetInterface;

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for bounded sets.
class ValidatedBoundedSetInterface
    : public virtual SetInterfaceBase {
  public:
    virtual ValidatedBoundedSetInterface* clone() const = 0;
    //! \brief Tests if the set is a inside of \a bx.
    virtual ValidatedLowerKleenean inside(const ExactBoxType& bx) const = 0;
    //! \brief Returns a bounding box for the set.
    virtual UpperBoxType bounding_box() const = 0;
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for overt sets, for which intersection with an open box is verifiable.
class ValidatedOvertSetInterface
    : public virtual SetInterfaceBase
{
  public:
    virtual ValidatedOvertSetInterface* clone() const = 0;
    //! \brief Tests if the set overlaps \a bx.
    virtual ValidatedLowerKleenean overlaps(const ExactBoxType& bx) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend ValidatedLowerKleenean overlap(const ValidatedOvertSetInterface& ovs, const ValidatedOpenSetInterface& ops);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for open sets.
class ValidatedOpenSetInterface
    : public virtual ValidatedOvertSetInterface
{
  public:
    virtual ValidatedOpenSetInterface* clone() const = 0;
    //! \brief Tests if the set covers of \a bx.
    virtual ValidatedLowerKleenean covers(const ExactBoxType& bx) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend ValidatedLowerKleenean overlap(const ValidatedOvertSetInterface& ovs, const ValidatedOpenSetInterface& ops);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedLowerKleenean inside(const ValidatedCompactSetInterface& ls, const ValidatedOpenSetInterface& rs);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for closed sets.
class ValidatedClosedSetInterface
    : public virtual SetInterfaceBase
{
  public:
    virtual ValidatedClosedSetInterface* clone() const = 0;
    //! \brief Tests if the set is separated from \a bx.
    virtual ValidatedLowerKleenean separated(const ExactBoxType& bx) const = 0;
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend ValidatedLowerKleenean separated(const ValidatedCompactSetInterface& cps, const ValidatedClosedSetInterface& cls);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for compact (closed and bounded) sets.
class ValidatedCompactSetInterface
    : public virtual ValidatedBoundedSetInterface,
      public virtual ValidatedClosedSetInterface
{
  public:
    virtual ValidatedCompactSetInterface* clone() const = 0;
    //virtual ValidatedSierpinskian empty() const = 0;
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedLowerKleenean inside(const ValidatedCompactSetInterface& ls, const ValidatedOpenSetInterface& rs);
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend ValidatedLowerKleenean separated(const ValidatedCompactSetInterface& cps, const ValidatedClosedSetInterface& cls);

};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
class ValidatedRegularSetInterface
    : public virtual ValidatedOpenSetInterface,
      public virtual ValidatedClosedSetInterface
{
    virtual ValidatedRegularSetInterface* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean overlap(const ValidatedLocatedSetInterface& ls, const ValidatedRegularSetInterface& rs);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean inside(const ValidatedLocatedSetInterface& ls, const ValidatedRegularSetInterface& rs);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean separated(const ValidatedLocatedSetInterface& ls, const ValidatedRegularSetInterface& rs);
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for located (overt and compact) sets.
class ValidatedLocatedSetInterface
    : public virtual ValidatedOvertSetInterface,
      public virtual ValidatedCompactSetInterface
{
    virtual ValidatedLocatedSetInterface* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean overlap(const ValidatedLocatedSetInterface& ls, const ValidatedRegularSetInterface& rs);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean inside(const ValidatedLocatedSetInterface& ls, const ValidatedRegularSetInterface& rs);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean separated(const ValidatedLocatedSetInterface& ls, const ValidatedRegularSetInterface& rs);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Complete set interface for bounded regular sets.
class ValidatedRegularLocatedSetInterface
    : public virtual ValidatedRegularSetInterface,
      public virtual ValidatedLocatedSetInterface
{
  public:
    virtual ValidatedRegularLocatedSetInterface* clone() const = 0;
};



class GridTreePaving;

//! \brief A Euclidean space \f$\R^d\f$ of dimension \a d.
class EuclideanSpace
{
  public:
    //! \brief The canonical type used for bounding sets in the space.
    typedef ExactBoxType BoundingDomainType;
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
    //! \brief The interface satisified by bounded regular sets.
    typedef RegularLocatedSetInterface RegularLocatedSetInterfaceType;
    //! \brief The type of approximations to sets in the space.
    typedef GridTreePaving SetApproximationType;
  public:
    EuclideanSpace(const SizeType& d) : _dimension(d) { }
    const SizeType& dimension() const { return this->_dimension; }
  private:
    SizeType _dimension;
};


} // namespace Ariadne


#endif // ARIADNE_SET_INTERFACE
