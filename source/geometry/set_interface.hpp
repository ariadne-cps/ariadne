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

#include "utility/declarations.hpp"
#include "utility/tribool.hpp"
#include "utility/writable.hpp"
#include "numeric/numeric.hpp"

#include "set.decl.hpp"
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
template<class T> class BoundedSetInterface<EffectiveTag,T>
    : public virtual SetInterfaceBase<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;
    using typename SetInterfaceBase<T>::BoundingSetType;

    virtual BoundedSetInterface<EffectiveTag,T>* clone() const = 0;
    //! \brief Tests if the set is a inside of \a bx.
    //! A set \a A is \em inside \a B if the closure of \a A is a subset of the interior of \a B.
    //! A set \f$A\f$ is \em inside \f$B\f$ if \f$\,\overline{\!A} \subset B^\circ\f$.
    virtual LowerKleenean inside(const BasicSetType& bx) const = 0;
    //! \brief Returns a bounding box for the set.
    //! If the set is empty, then the first component of the result should be empty.
    virtual BoundingSetType bounding_box() const = 0;
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for overt sets, for which intersection with an open box is verifiable.
template<class T> class OvertSetInterface<EffectiveTag,T>
    : public virtual SetInterfaceBase<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual OvertSetInterface<EffectiveTag,T>* clone() const = 0;
    //! \brief Tests if the set overlaps \a bx.
    //! Sets \a A and \a B \em overlap if the interiors of \a A and \a B intersect.
    //! Sets \f$A\f$ and \f$B\f$ \em overlap if \f$A^\circ \cap B^\circ \neq \emptyset\f$.
    virtual LowerKleenean overlaps(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend LowerKleenean overlap(const OvertSetInterface<EffectiveTag,T>& ovs, const OpenSetInterface<EffectiveTag,T>& ops);
    friend ValidatedLowerKleenean overlap(const OvertSetInterface<EffectiveTag,T>& ovs, const OpenSetInterface<EffectiveTag,T>& ops, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for open sets.
template<class T> class OpenSetInterface<EffectiveTag,T>
    : public virtual OvertSetInterface<EffectiveTag,T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual OpenSetInterface<EffectiveTag,T>* clone() const = 0;
    //! \brief Tests if the set covers of \a bx.
    //! A set \a A \em covers \a B if the interiors of \a A is a superset of the closure of \a B.
    //! A set \f$A\f$ \em covers \f$B\f$ if \f$A^\circ \supset \overline{B}\f$.
    virtual LowerKleenean covers(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend LowerKleenean overlap(const OvertSetInterface<EffectiveTag,T>& ovs, const OpenSetInterface<EffectiveTag,T>& ops);
    friend ValidatedLowerKleenean overlap(const OvertSetInterface<EffectiveTag,T>& ovs, const OpenSetInterface<EffectiveTag,T>& ops, const RawFloatDP& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend LowerKleenean inside(const CompactSetInterface<EffectiveTag,T>& ls, const OpenSetInterface<EffectiveTag,T>& rs);
    friend ValidatedLowerKleenean inside(const CompactSetInterface<EffectiveTag,T>& ls, const OpenSetInterface<EffectiveTag,T>& rs, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for closed sets.
template<class T> class ClosedSetInterface<EffectiveTag,T>
    : public virtual SetInterfaceBase<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual ClosedSetInterface<EffectiveTag,T>* clone() const = 0;
    //! \brief Tests if the set is separated from \a bx.
    //! A set \a A is \em separated from \a B if the closures of \a A and \a B are disjoint.
    //! A set \f$A\f$ is \em separated from \f$B\f$ if \f$\,\overline{\!A} \cap \overline{B} = \emptyset\f$.
    virtual LowerKleenean separated(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend LowerKleenean separated(const CompactSetInterface<EffectiveTag,T>& cps, const ClosedSetInterface<EffectiveTag,T>& cls);
    friend ValidatedLowerKleenean separated(const CompactSetInterface<EffectiveTag,T>& cps, const ClosedSetInterface<EffectiveTag,T>& cls, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for compact (closed and bounded) sets.
template<class T> class CompactSetInterface<EffectiveTag,T>
    : public virtual BoundedSetInterface<EffectiveTag,T>,
      public virtual ClosedSetInterface<EffectiveTag,T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual CompactSetInterface<EffectiveTag,T>* clone() const = 0;
    //virtual ValidatedSierpinskian empty() const = 0;
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend LowerKleenean inside(const CompactSetInterface<EffectiveTag,T>& ls, const OpenSetInterface<EffectiveTag,T>& rs);
    friend ValidatedLowerKleenean inside(const CompactSetInterface<EffectiveTag,T>& ls, const OpenSetInterface<EffectiveTag,T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend LowerKleenean separated(const CompactSetInterface<EffectiveTag,T>& cps, const ClosedSetInterface<EffectiveTag,T>& cls);
    friend ValidatedLowerKleenean separated(const CompactSetInterface<EffectiveTag,T>& cps, const ClosedSetInterface<EffectiveTag,T>& cls, const RawFloatDP& eps);

};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class RegularSetInterface<EffectiveTag,T>
    : public virtual OpenSetInterface<EffectiveTag,T>,
      public virtual ClosedSetInterface<EffectiveTag,T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual RegularSetInterface<EffectiveTag,T>* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend Kleenean overlap(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs);
    friend ValidatedKleenean overlap(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend Kleenean inside(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs);
    friend ValidatedKleenean inside(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend Kleenean separated(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs);
    friend ValidatedKleenean separated(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const RawFloatDP& eps);
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for located (overt and compact) sets.
template<class T> class LocatedSetInterface<EffectiveTag,T>
    : public virtual OvertSetInterface<EffectiveTag,T>,
      public virtual CompactSetInterface<EffectiveTag,T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual LocatedSetInterface<EffectiveTag,T>* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend Kleenean overlap(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs);
    friend ValidatedKleenean overlap(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend Kleenean inside(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs);
    friend ValidatedKleenean inside(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend Kleenean separated(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs);
    friend ValidatedKleenean separated(const LocatedSetInterface<EffectiveTag,T>& ls, const RegularSetInterface<EffectiveTag,T>& rs, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Complete set interface for bounded regular sets.
template<class T> class RegularLocatedSetInterface<EffectiveTag,T>
    : public virtual RegularSetInterface<EffectiveTag,T>,
      public virtual LocatedSetInterface<EffectiveTag,T>
{
  public:
    virtual RegularLocatedSetInterface<EffectiveTag,T>* clone() const = 0;
};



template<class T> class BoundedSetInterface<ValidatedTag,T>;
template<class T> class OpenSetInterface<ValidatedTag,T>;
template<class T> class ClosedSetInterface<ValidatedTag,T>;
template<class T> class OvertSetInterface<ValidatedTag,T>;
template<class T> class CompactSetInterface<ValidatedTag,T>;
template<class T> class RegularSetInterface<ValidatedTag,T>;
template<class T> class LocatedSetInterface<ValidatedTag,T>;
template<class T> class RegularLocatedSetInterface<ValidatedTag,T>;

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for bounded sets.
template<class T> class BoundedSetInterface<ValidatedTag,T>
    : public virtual SetInterfaceBase<T> {
  public:
    using typename SetInterfaceBase<T>::BasicSetType;
    using typename SetInterfaceBase<T>::BoundingSetType;

    virtual BoundedSetInterface<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set is a inside of \a bx.
    virtual ValidatedLowerKleenean inside(const BasicSetType& bx) const = 0;
    //! \brief Returns a bounding box for the set.
    virtual BoundingSetType bounding_box() const = 0;
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for overt sets, for which intersection with an open box is verifiable.
template<class T> class OvertSetInterface<ValidatedTag,T>
    : public virtual SetInterfaceBase<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual OvertSetInterface<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set overlaps \a bx.
    virtual ValidatedLowerKleenean overlaps(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend ValidatedLowerKleenean overlap(const OvertSetInterface<ValidatedTag,T>& ovs, const OpenSetInterface<ValidatedTag,T>& ops);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for open sets.
template<class T> class OpenSetInterface<ValidatedTag,T>
    : public virtual OvertSetInterface<ValidatedTag,T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual OpenSetInterface<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set covers of \a bx.
    virtual ValidatedLowerKleenean covers(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend ValidatedLowerKleenean overlap(const OvertSetInterface<ValidatedTag,T>& ovs, const OpenSetInterface<ValidatedTag,T>& ops);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedLowerKleenean inside(const CompactSetInterface<ValidatedTag,T>& ls, const OpenSetInterface<ValidatedTag,T>& rs);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for closed sets.
template<class T> class ClosedSetInterface<ValidatedTag,T>
    : public virtual SetInterfaceBase<T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual ClosedSetInterface<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set is separated from \a bx.
    virtual ValidatedLowerKleenean separated(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend ValidatedLowerKleenean separated(const CompactSetInterface<ValidatedTag,T>& cps, const ClosedSetInterface<ValidatedTag,T>& cls);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for compact (closed and bounded) sets.
template<class T> class CompactSetInterface<ValidatedTag,T>
    : public virtual BoundedSetInterface<ValidatedTag,T>,
      public virtual ClosedSetInterface<ValidatedTag,T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual CompactSetInterface<ValidatedTag,T>* clone() const = 0;
    //virtual ValidatedSierpinskian empty() const = 0;
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedLowerKleenean inside(const CompactSetInterface<ValidatedTag,T>& ls, const OpenSetInterface<ValidatedTag,T>& rs);
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend ValidatedLowerKleenean separated(const CompactSetInterface<ValidatedTag,T>& cps, const ClosedSetInterface<ValidatedTag,T>& cls);

};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class RegularSetInterface<ValidatedTag,T>
    : public virtual OpenSetInterface<ValidatedTag,T>,
      public virtual ClosedSetInterface<ValidatedTag,T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual RegularSetInterface<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean overlap(const LocatedSetInterface<ValidatedTag,T>& ls, const RegularSetInterface<ValidatedTag,T>& rs);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean inside(const LocatedSetInterface<ValidatedTag,T>& ls, const RegularSetInterface<ValidatedTag,T>& rs);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean separated(const LocatedSetInterface<ValidatedTag,T>& ls, const RegularSetInterface<ValidatedTag,T>& rs);
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for located (overt and compact) sets.
template<class T> class LocatedSetInterface<ValidatedTag,T>
    : public virtual OvertSetInterface<ValidatedTag,T>,
      public virtual CompactSetInterface<ValidatedTag,T>
{
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual LocatedSetInterface<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean overlap(const LocatedSetInterface<ValidatedTag,T>& ls, const RegularSetInterface<ValidatedTag,T>& rs);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean inside(const LocatedSetInterface<ValidatedTag,T>& ls, const RegularSetInterface<ValidatedTag,T>& rs);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean separated(const LocatedSetInterface<ValidatedTag,T>& ls, const RegularSetInterface<ValidatedTag,T>& rs);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Complete set interface for bounded regular sets.
template<class T> class RegularLocatedSetInterface<ValidatedTag,T>
    : public virtual RegularSetInterface<ValidatedTag,T>,
      public virtual LocatedSetInterface<ValidatedTag,T>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    virtual RegularLocatedSetInterface<ValidatedTag,T>* clone() const = 0;
};



class GridTreePaving;

//! \brief A Euclidean space \f$\R^d\f$ of dimension \a d.
class EuclideanSpace
{
    using P=EffectiveTag; using T=RealVector;
  public:
    //! \brief The canonical type used for bounding sets in the space.
    typedef ExactBoxType BoundingDomainType;
    //! \brief The interface satisified by bounded sets in the space.
    typedef BoundedSetInterface<P,T> BoundedSetInterfaceType;
    //! \brief The interface satisified by overt sets in the space.
    typedef OvertSetInterface<P,T> OvertSetInterfaceType;
    //! \brief The interface satisified by over sets in the space.
    typedef OpenSetInterface<P,T> OpenSetInterfaceType;
    //! \brief The interface satisified by closed sets in the space.
    typedef ClosedSetInterface<P,T> ClosedSetInterfaceType;
    //! \brief The interface satisified by compact sets in the space.
    typedef CompactSetInterface<P,T> CompactSetInterfaceType;
    //! \brief The interface satisified by regular sets in the space.
    typedef RegularSetInterface<P,T> RegularSetInterfaceType;
    //! \brief The interface satisified by located sets in the space.
    typedef LocatedSetInterface<P,T> LocatedSetInterfaceType;
    //! \brief The interface satisified by bounded regular sets.
    typedef RegularLocatedSetInterface<P,T> RegularLocatedSetInterfaceType;
    //! \brief The type of approximations to sets in the space.
    typedef GridTreePaving SetApproximationType;
  public:
    EuclideanSpace(const SizeType& d) : _dimension(d) { }
    const SizeType& dimension() const { return this->_dimension; }
  private:
    SizeType _dimension;
};


} // namespace Ariadne


#endif // ARIADNE_SET_INTERFACE_HPP
