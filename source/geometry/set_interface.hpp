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
 *  \brief Interface<T>s for open, closed, overt and compact subsets of Euclidean space.
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

using DimensionOne = SizeOne;

template<class T> struct SetTraits;
template<> struct SetTraits<Real> {
    typedef DimensionOne DimensionType;
    typedef FloatDPExactInterval BasicSetType;
    typedef FloatDPUpperInterval BoundingSetType;
};
template<> struct SetTraits<RealVector> {
    typedef Ariadne::DimensionType DimensionType;
    typedef FloatDPExactBox BasicSetType;
    typedef FloatDPUpperBox BoundingSetType;
};
template<class T> using DimensionOfType = typename SetTraits<T>::DimensionType;
template<class T> using BasicSetType = typename SetTraits<T>::BasicSetType;
template<class T> using BoundingSetType = typename SetTraits<T>::BoundingSetType;

using EuclideanSetTraits = SetTraits<RealVector>;

template<class UB> class Interval;
typedef FloatDPExactInterval ExactIntervalType;

template<class IVL> class Box;
typedef FloatDPExactBox ExactBoxType;
typedef FloatDPUpperBox UpperBoxType;

template<class T> class SetInterfaceBase;
template<class T> class BoundedSetInterface;
template<class T> class OpenSetInterface;
template<class T> class ClosedSetInterface;
template<class T> class OvertSetInterface;
template<class T> class CompactSetInterface;
template<class T> class RegularSetInterface;
template<class T> class LocatedSetInterface;
template<class T> class RegularLocatedSetInterface;
template<class T> using SetInterface = RegularLocatedSetInterface<T>;

using EuclideanSetInterfaceBase = SetInterfaceBase<RealVector>;

using EuclideanBoundedSetInterface = BoundedSetInterface<RealVector>;
using EuclideanOpenSetInterface = OpenSetInterface<RealVector>;
using EuclideanClosedSetInterface = ClosedSetInterface<RealVector>;
using EuclideanOvertSetInterface = OvertSetInterface<RealVector>;
using EuclideanCompactSetInterface = CompactSetInterface<RealVector>;
using EuclideanRegularSetInterface = RegularSetInterface<RealVector>;
using EuclideanLocatedSetInterface = LocatedSetInterface<RealVector>;
using EuclideanRegularLocatedSetInterface = RegularLocatedSetInterface<RealVector>;
using EuclideanSetInterface = SetInterface<RealVector>;

template<class T> class ValidatedBoundedSetInterface;
template<class T> class ValidatedOpenSetInterface;
template<class T> class ValidatedClosedSetInterface;
template<class T> class ValidatedOvertSetInterface;
template<class T> class ValidatedCompactSetInterface;
template<class T> class ValidatedRegularSetInterface;
template<class T> class ValidatedLocatedSetInterface;
template<class T> class ValidatedRegularLocatedSetInterface;
using ValidatedEuclideanBoundedSetInterface = ValidatedBoundedSetInterface<RealVector>;
using ValidatedEuclideanOpenSetInterface = ValidatedOpenSetInterface<RealVector>;
using ValidatedEuclideanClosedSetInterface = ValidatedClosedSetInterface<RealVector>;
using ValidatedEuclideanOvertSetInterface = ValidatedOvertSetInterface<RealVector>;
using ValidatedEuclideanCompactSetInterface = ValidatedCompactSetInterface<RealVector>;
using ValidatedEuclideanRegularSetInterface = ValidatedRegularSetInterface<RealVector>;
using ValidatedEuclideanLocatedSetInterface = ValidatedLocatedSetInterface<RealVector>;
using ValidatedEuclideanRegularLocatedSetInterface = ValidatedRegularLocatedSetInterface<RealVector>;

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
template<class T> class SetInterfaceBase
    : public virtual WritableInterface
{
  public:
    //! \brief The type of element of the set.
    typedef T ElementType;
    //! \brief The type representing the dimension of the set.
    typedef Ariadne::DimensionOfType<T> DimensionType;
    //! \brief The type of basic set in the space.
    typedef Ariadne::BasicSetType<T> BasicSetType;
    //! \brief The type of basic set in the space.
    typedef Ariadne::BoundingSetType<T> BoundingSetType;

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
template<class T> class BoundedSetInterface
    : public virtual SetInterfaceBase<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual BoundedSetInterface<T>* clone() const = 0;
    //! \brief Tests if the set is a inside of \a bx.
    //! A set \a A is \em inside \a B if the closure of \a A is a subset of the interior of \a B.
    //! A set \f$A\f$ is \em inside \f$B\f$ if \f$\,\overline{\!A} \subset B^\circ\f$.
    virtual LowerKleenean inside(const BasicSetType& bx) const = 0;
    virtual ValidatedLowerKleenean inside(const BasicSetType& bx, Effort eff) const = 0;
    //! \brief Returns a bounding box for the set.
    //! If the set is empty, then the first component of the result should be empty.
    virtual UpperBoxType bounding_box() const = 0;
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for overt sets, for which intersection with an open box is verifiable.
template<class T> class OvertSetInterface
    : public virtual SetInterfaceBase<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual OvertSetInterface<T>* clone() const = 0;
    //! \brief Tests if the set overlaps \a bx.
    //! Sets \a A and \a B \em overlap if the interiors of \a A and \a B intersect.
    //! Sets \f$A\f$ and \f$B\f$ \em overlap if \f$A^\circ \cap B^\circ \neq \emptyset\f$.
    virtual LowerKleenean overlaps(const BasicSetType& bx) const = 0;
    virtual ValidatedLowerKleenean overlaps(const BasicSetType& bx, Effort eff) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend LowerKleenean overlap(const OvertSetInterface<T>& ovs, const OpenSetInterface<T>& ops);
    friend ValidatedLowerKleenean overlap(const OvertSetInterface<T>& ovs, const OpenSetInterface<T>& ops, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for open sets.
template<class T> class OpenSetInterface
    : public virtual OvertSetInterface<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual OpenSetInterface<T>* clone() const = 0;
    //! \brief Tests if the set covers of \a bx.
    //! A set \a A \em covers \a B if the interiors of \a A is a superset of the closure of \a B.
    //! A set \f$A\f$ \em covers \f$B\f$ if \f$A^\circ \supset \overline{B}\f$.
    virtual LowerKleenean covers(const BasicSetType& bx) const = 0;
    virtual ValidatedLowerKleenean covers(const BasicSetType& bx, Effort eff) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend LowerKleenean overlap(const OvertSetInterface<T>& ovs, const OpenSetInterface<T>& ops);
    friend ValidatedLowerKleenean overlap(const OvertSetInterface<T>& ovs, const OpenSetInterface<T>& ops, const RawFloatDP& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend LowerKleenean inside(const CompactSetInterface<T>& ls, const OpenSetInterface<T>& rs);
    friend ValidatedLowerKleenean inside(const CompactSetInterface<T>& ls, const OpenSetInterface<T>& rs, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for closed sets.
template<class T> class ClosedSetInterface
    : public virtual SetInterfaceBase<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual ClosedSetInterface<T>* clone() const = 0;
    //! \brief Tests if the set is separated from \a bx.
    //! A set \a A is \em separated from \a B if the closures of \a A and \a B are disjoint.
    //! A set \f$A\f$ is \em separated from \f$B\f$ if \f$\,\overline{\!A} \cap \overline{B} = \emptyset\f$.
    virtual LowerKleenean separated(const BasicSetType& bx) const = 0;
    virtual ValidatedLowerKleenean separated(const BasicSetType& bx, Effort eff) const = 0;
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend LowerKleenean separated(const CompactSetInterface<T>& cps, const ClosedSetInterface<T>& cls);
    friend ValidatedLowerKleenean separated(const CompactSetInterface<T>& cps, const ClosedSetInterface<T>& cls, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for compact (closed and bounded) sets.
template<class T> class CompactSetInterface
    : public virtual BoundedSetInterface<T>,
      public virtual ClosedSetInterface<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual CompactSetInterface<T>* clone() const = 0;
    //virtual ValidatedSierpinskian empty() const = 0;
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend LowerKleenean inside(const CompactSetInterface<T>& ls, const OpenSetInterface<T>& rs);
    friend ValidatedLowerKleenean inside(const CompactSetInterface<T>& ls, const OpenSetInterface<T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend LowerKleenean separated(const CompactSetInterface<T>& cps, const ClosedSetInterface<T>& cls);
    friend ValidatedLowerKleenean separated(const CompactSetInterface<T>& cps, const ClosedSetInterface<T>& cls, const RawFloatDP& eps);

};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class RegularSetInterface
    : public virtual OpenSetInterface<T>,
      public virtual ClosedSetInterface<T>
{
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual RegularSetInterface<T>* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend Kleenean overlap(const LocatedSetInterface<T>& ls, const RegularSetInterface<T>& rs);
    friend ValidatedKleenean overlap(const LocatedSetInterface<T>& ls, const RegularSetInterface<T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend Kleenean inside(const LocatedSetInterface<T>& ls, const RegularSetInterface<T>& rs);
    friend ValidatedKleenean inside(const LocatedSetInterface<T>& ls, const RegularSetInterface<T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend Kleenean separated(const LocatedSetInterface<T>& ls, const RegularSetInterface<T>& rs);
    friend ValidatedKleenean separated(const LocatedSetInterface<T>& ls, const RegularSetInterface<T>& rs, const RawFloatDP& eps);
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for located (overt and compact) sets.
template<class T> class LocatedSetInterface
    : public virtual OvertSetInterface<T>,
      public virtual CompactSetInterface<T>
{
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual LocatedSetInterface<T>* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend Kleenean overlap(const LocatedSetInterface<T>& ls, const RegularSetInterface<T>& rs);
    friend ValidatedKleenean overlap(const LocatedSetInterface<T>& ls, const RegularSetInterface<T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend Kleenean inside(const LocatedSetInterface<T>& ls, const RegularSetInterface<T>& rs);
    friend ValidatedKleenean inside(const LocatedSetInterface<T>& ls, const RegularSetInterface<T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend Kleenean separated(const LocatedSetInterface<T>& ls, const RegularSetInterface<T>& rs);
    friend ValidatedKleenean separated(const LocatedSetInterface<T>& ls, const RegularSetInterface<T>& rs, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Complete set interface for bounded regular sets.
template<class T> class RegularLocatedSetInterface
    : public virtual RegularSetInterface<T>,
      public virtual LocatedSetInterface<T>
{
  public:
    virtual RegularLocatedSetInterface<T>* clone() const = 0;
};


inline OutputStream& operator<<(OutputStream& os, const WritableInterface& w);



template<class T> class ValidatedBoundedSetInterface;
template<class T> class ValidatedOpenSetInterface;
template<class T> class ValidatedClosedSetInterface;
template<class T> class ValidatedOvertSetInterface;
template<class T> class ValidatedCompactSetInterface;
template<class T> class ValidatedRegularSetInterface;
template<class T> class ValidatedLocatedSetInterface;
template<class T> class ValidatedRegularLocatedSetInterface;

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for bounded sets.
template<class T> class ValidatedBoundedSetInterface
    : public virtual SetInterfaceBase<T> {
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual ValidatedBoundedSetInterface<T>* clone() const = 0;
    //! \brief Tests if the set is a inside of \a bx.
    virtual ValidatedLowerKleenean inside(const BasicSetType& bx) const = 0;
    //! \brief Returns a bounding box for the set.
    virtual UpperBoxType bounding_box() const = 0;
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for overt sets, for which intersection with an open box is verifiable.
template<class T> class ValidatedOvertSetInterface
    : public virtual SetInterfaceBase<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual ValidatedOvertSetInterface<T>* clone() const = 0;
    //! \brief Tests if the set overlaps \a bx.
    virtual ValidatedLowerKleenean overlaps(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend ValidatedLowerKleenean overlap(const ValidatedOvertSetInterface<T>& ovs, const ValidatedOpenSetInterface<T>& ops);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for open sets.
template<class T> class ValidatedOpenSetInterface
    : public virtual ValidatedOvertSetInterface<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual ValidatedOpenSetInterface<T>* clone() const = 0;
    //! \brief Tests if the set covers of \a bx.
    virtual ValidatedLowerKleenean covers(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend ValidatedLowerKleenean overlap(const ValidatedOvertSetInterface<T>& ovs, const ValidatedOpenSetInterface<T>& ops);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedLowerKleenean inside(const ValidatedCompactSetInterface<T>& ls, const ValidatedOpenSetInterface<T>& rs);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for closed sets.
template<class T> class ValidatedClosedSetInterface
    : public virtual SetInterfaceBase<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual ValidatedClosedSetInterface<T>* clone() const = 0;
    //! \brief Tests if the set is separated from \a bx.
    virtual ValidatedLowerKleenean separated(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend ValidatedLowerKleenean separated(const ValidatedCompactSetInterface<T>& cps, const ValidatedClosedSetInterface<T>& cls);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for compact (closed and bounded) sets.
template<class T> class ValidatedCompactSetInterface
    : public virtual ValidatedBoundedSetInterface<T>,
      public virtual ValidatedClosedSetInterface<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual ValidatedCompactSetInterface<T>* clone() const = 0;
    //virtual ValidatedSierpinskian empty() const = 0;
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedLowerKleenean inside(const ValidatedCompactSetInterface<T>& ls, const ValidatedOpenSetInterface<T>& rs);
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend ValidatedLowerKleenean separated(const ValidatedCompactSetInterface<T>& cps, const ValidatedClosedSetInterface<T>& cls);

};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class ValidatedRegularSetInterface
    : public virtual ValidatedOpenSetInterface<T>,
      public virtual ValidatedClosedSetInterface<T>
{
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual ValidatedRegularSetInterface<T>* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean overlap(const ValidatedLocatedSetInterface<T>& ls, const ValidatedRegularSetInterface<T>& rs);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean inside(const ValidatedLocatedSetInterface<T>& ls, const ValidatedRegularSetInterface<T>& rs);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean separated(const ValidatedLocatedSetInterface<T>& ls, const ValidatedRegularSetInterface<T>& rs);
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for located (overt and compact) sets.
template<class T> class ValidatedLocatedSetInterface
    : public virtual ValidatedOvertSetInterface<T>,
      public virtual ValidatedCompactSetInterface<T>
{
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual ValidatedLocatedSetInterface<T>* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean overlap(const ValidatedLocatedSetInterface<T>& ls, const ValidatedRegularSetInterface<T>& rs);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean inside(const ValidatedLocatedSetInterface<T>& ls, const ValidatedRegularSetInterface<T>& rs);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean separated(const ValidatedLocatedSetInterface<T>& ls, const ValidatedRegularSetInterface<T>& rs);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Complete set interface for bounded regular sets.
template<class T> class ValidatedRegularLocatedSetInterface
    : public virtual ValidatedRegularSetInterface<T>,
      public virtual ValidatedLocatedSetInterface<T>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    virtual ValidatedRegularLocatedSetInterface<T>* clone() const = 0;
};



class GridTreePaving;

//! \brief A Euclidean space \f$\R^d\f$ of dimension \a d.
class EuclideanSpace
{
    using T=RealVector;
  public:
    //! \brief The canonical type used for bounding sets in the space.
    typedef ExactBoxType BoundingDomainType;
    //! \brief The interface satisified by bounded sets in the space.
    typedef BoundedSetInterface<T> BoundedSetInterfaceType;
    //! \brief The interface satisified by overt sets in the space.
    typedef OvertSetInterface<T> OvertSetInterfaceType;
    //! \brief The interface satisified by over sets in the space.
    typedef OpenSetInterface<T> OpenSetInterfaceType;
    //! \brief The interface satisified by closed sets in the space.
    typedef ClosedSetInterface<T> ClosedSetInterfaceType;
    //! \brief The interface satisified by compact sets in the space.
    typedef CompactSetInterface<T> CompactSetInterfaceType;
    //! \brief The interface satisified by regular sets in the space.
    typedef RegularSetInterface<T> RegularSetInterfaceType;
    //! \brief The interface satisified by located sets in the space.
    typedef LocatedSetInterface<T> LocatedSetInterfaceType;
    //! \brief The interface satisified by bounded regular sets.
    typedef RegularLocatedSetInterface<T> RegularLocatedSetInterfaceType;
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
