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
    typedef EffectiveNumber ElementType;
    typedef ValidatedNumber ValidatedElementType;
};
template<> struct SetTraits<RealVector> {
    typedef Ariadne::DimensionType DimensionType;
    typedef FloatDPExactBox BasicSetType;
    typedef FloatDPUpperBox BoundingSetType;
    typedef Vector<EffectiveNumber> ElementType;
    typedef Vector<ValidatedNumber> ValidatedElementType;
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
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2, Accuracy acc) const { return overlap(s1,s2,acc); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2) const { return overlap(s1,s2); }
};

struct Separated {
    constexpr const char* code() const { return "separated"; }
    template<class S> decltype(auto) operator()(S const& s, ExactBoxType const& bx, Effort eff) const { return s.separated(bx,eff); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2, Effort eff) const { return separated(s1,s2,eff); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2, Accuracy acc) const { return separated(s1,s2,acc); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2) const { return separated(s1,s2); }
};

struct Inside {
    constexpr const char* code() const { return "inside"; }
    template<class S> decltype(auto) operator()(S const& s, ExactBoxType const& bx, Effort eff) const { return s.inside(bx,eff); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2, Effort eff) const { return inside(s1,s2,eff); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2, Accuracy acc) const { return inside(s1,s2,acc); }
};

struct Covers {
    constexpr const char* code() const { return "covers"; }
    template<class S> decltype(auto) operator()(S const& s, ExactBoxType const& bx, Effort eff) const { return s.covers(bx,eff); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2, Effort eff) const { return inside(s2,s1,eff); }
    template<class S1, class S2> decltype(auto) operator()(S1 const& s1, S2 const& s2, Accuracy acc) const { return covers(s1,s2,acc); }
};

template<class T> class SetOperations {

  public:
    typedef typename SetTraits<T>::ElementType ElementType;
    typedef typename SetTraits<T>::BasicSetType BasicSetType;

    static ValidatedLowerKleenean contains(EffectiveOpenSetInterface<T> const& ops, ElementType const& pt, Effort eff);
    static ValidatedUpperKleenean contains(EffectiveClosedSetInterface<T> const& cls, ElementType const& pt, Effort eff);
    static ValidatedKleenean contains(EffectiveRegularSetInterface<T> const& rs, ElementType const& pt, Effort eff);
    static ValidatedKleenean contains(EffectiveRegularSetInterface<T> const& rs, ElementType const& pt, Accuracy acc);

    static ValidatedLowerKleenean overlap(EffectiveOvertSetInterface<T> const& ovs, EffectiveOpenSetInterface<T> const& ops, Effort eff);
    static ValidatedLowerKleenean inside(EffectiveCompactSetInterface<T> const& cps, EffectiveOpenSetInterface<T> const& ops, Effort eff);
    static ValidatedLowerKleenean separated(EffectiveCompactSetInterface<T> const& cps, EffectiveClosedSetInterface<T> const& cls, Effort eff);
    static ValidatedKleenean overlap(EffectiveLocatedSetInterface<T> const& ls, EffectiveRegularSetInterface<T> const& rs, Effort eff);
    static ValidatedKleenean inside(EffectiveLocatedSetInterface<T> const& ls, EffectiveRegularSetInterface<T> const& rs, Effort eff);
    static ValidatedKleenean separated(EffectiveLocatedSetInterface<T> const& ls, EffectiveRegularSetInterface<T> const& rs, Effort eff);
    static ValidatedKleenean overlap(EffectiveLocatedSetInterface<T> const& ls, EffectiveRegularSetInterface<T> const& rs, Accuracy acc);
    static ValidatedKleenean inside(EffectiveLocatedSetInterface<T> const& ls, EffectiveRegularSetInterface<T> const& rs, Accuracy acc);
    static ValidatedKleenean separated(EffectiveLocatedSetInterface<T> const& ls, EffectiveRegularSetInterface<T> const& rs, Accuracy acc);

    static ValidatedLowerKleenean overlap(ValidatedOvertSetInterface<T> const& ovs, ValidatedOpenSetInterface<T> const& ops);
    static ValidatedLowerKleenean inside(ValidatedCompactSetInterface<T> const& cps, ValidatedOpenSetInterface<T> const& ops);
    static ValidatedLowerKleenean separated(ValidatedCompactSetInterface<T> const& cps, ValidatedClosedSetInterface<T> const& cls);
    static ValidatedKleenean overlap(ValidatedLocatedSetInterface<T> const& ls, ValidatedRegularSetInterface<T> const& rs);
    static ValidatedKleenean inside(ValidatedLocatedSetInterface<T> const& ls, ValidatedRegularSetInterface<T> const& rs);
    static ValidatedKleenean separated(ValidatedLocatedSetInterface<T> const& ls, ValidatedRegularSetInterface<T> const& rs);
};


//! \brief Base class for sets described by predicates involving boxes.
template<class T> class SetInterfaceBase
    : public virtual WritableInterface
//    , public std::enable_shared_from_this<SetInterfaceBase<T>>
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
    //! \brief The type of element of the set.
    typedef typename SetTraits<T>::ValidatedElementType ValidatedElementType;

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
    using P=EffectiveTag;
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual BoundedSetInterface<EffectiveTag,T>* clone() const = 0;
    //! \brief Tests if the set is a inside of \a bx.
    //! A set \a A is \em inside \a B if the closure of \a A is a subset of the interior of \a B.
    //! A set \f$A\f$ is \em inside \f$B\f$ if \f$\,\overline{\!A} \subset B^\circ\f$.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const { return this->_inside(bx); }
    inline ValidatedLowerKleenean inside(const BasicSetType& bx, Effort eff) const { return this->_inside(bx); }
        //! \brief Returns a bounding box for the set.
    //! If the set is empty, then the first component of the result should be empty.
    virtual UpperBoxType bounding_box() const = 0;
  private: public:
    virtual ValidatedLowerKleenean _inside(const BasicSetType& bx) const = 0;
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for overt sets, for which intersection with an open box is verifiable.
template<class T> class OvertSetInterface<EffectiveTag,T>
    : public virtual SetInterfaceBase<T>
{
    using P=EffectiveTag;
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    inline OvertSetInterface<P,T>* _copy() const { return this->clone(); }
    virtual OvertSetInterface<P,T>* clone() const = 0;
    //! \brief Tests if the set overlaps \a ops.
    LowerKleenean overlaps(const OpenSetInterface<P,T>& ops) const;
    //! \brief Tests if the set overlaps \a ops.
    ValidatedLowerKleenean overlaps(const OpenSetInterface<P,T>& ops, Effort eff) const {
        return SetOperations<T>::overlap(*this,ops,eff); }
    //! \brief Tests if the set overlaps \a bx.
    //! Sets \a A and \a B \em overlap if the interiors of \a A and \a B intersect.
    //! Sets \f$A\f$ and \f$B\f$ \em overlap if \f$A^\circ \cap B^\circ \neq \emptyset\f$.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const {
        return this->_overlaps(bx); }
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx, Effort eff) const {
        return this->_overlaps(bx,eff); }
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a acc.
    friend LowerKleenean overlap(const OvertSetInterface<P,T>& ovs, const OpenSetInterface<P,T>& ops) {
        return ovs.overlaps(ops); }
    friend ValidatedLowerKleenean overlap(const OvertSetInterface<P,T>& ovs, const OpenSetInterface<P,T>& ops, Effort eff) {
        return ovs.overlaps(ops,eff); }
  private: public:
    virtual ValidatedLowerKleenean _overlaps(const BasicSetType& bx) const = 0;
    virtual ValidatedLowerKleenean _overlaps(const BasicSetType& bx, Effort eff) const;
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for open sets.
template<class T> class OpenSetInterface<EffectiveTag,T>
    : public virtual OvertSetInterface<EffectiveTag,T>
{
    using P=EffectiveTag;
  public:
    using typename SetInterfaceBase<T>::BasicSetType;
    using typename SetInterfaceBase<T>::ElementType;

    inline OpenSetInterface<P,T>* _copy() const { return this->clone(); }
    virtual OpenSetInterface<P,T>* clone() const = 0;
    //! \brief Tests if the set contains \a pt.
    LowerKleenean contains(const ElementType& pt) const;
    ValidatedLowerKleenean contains(const ElementType& pt, Effort eff) const {
        return SetOperations<T>::contains(*this,pt,eff); }
    //! \brief Tests if the set covers of \a bx.
    //! A set \a A \em covers \a B if the interiors of \a A is a superset of the closure of \a B.
    //! A set \f$A\f$ \em covers \f$B\f$ if \f$A^\circ \supset \overline{B}\f$.
    inline ValidatedLowerKleenean covers(const BasicSetType& bx) const {
        return this->_covers(bx); }
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a acc.
    friend LowerKleenean overlap(const OvertSetInterface<P,T>& ovs, const OpenSetInterface<P,T>& ops);
    friend ValidatedLowerKleenean overlap(const OvertSetInterface<P,T>& ovs, const OpenSetInterface<P,T>& ops, Effort eff);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a acc.
    friend LowerKleenean inside(const CompactSetInterface<P,T>& ls, const OpenSetInterface<P,T>& rs);
    friend ValidatedLowerKleenean inside(const CompactSetInterface<P,T>& ls, const OpenSetInterface<P,T>& rs, Effort eff);
  private: public:
    virtual ValidatedLowerKleenean _covers(const BasicSetType& bx) const = 0;
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for closed sets.
template<class T> class ClosedSetInterface<EffectiveTag,T>
    : public virtual SetInterfaceBase<T>
{
    using P=EffectiveTag;
  public:
    using typename SetInterfaceBase<T>::ElementType;
    using typename SetInterfaceBase<T>::BasicSetType;

    inline ClosedSetInterface<P,T>* _copy() const { return this->clone(); }
    virtual ClosedSetInterface<P,T>* clone() const = 0;
    //! \brief Tests if the set (does not) contain \a pt.
    UpperKleenean contains(const ElementType& pt) const;
    ValidatedUpperKleenean contains(const ElementType& pt, Effort eff) const {
        return SetOperations<T>::contains(*this,pt,eff); }
    //! \brief Tests if the set is separated from \a bx.
    //! A set \a A is \em separated from \a B if the closures of \a A and \a B are disjoint.
    //! A set \f$A\f$ is \em separated from \f$B\f$ if \f$\,\overline{\!A} \cap \overline{B} = \emptyset\f$.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const {
        return this->_separated(bx); }
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a acc.
    friend LowerKleenean separated(const CompactSetInterface<P,T>& cps, const ClosedSetInterface<P,T>& cls);
    friend ValidatedLowerKleenean separated(const CompactSetInterface<P,T>& cps, const ClosedSetInterface<P,T>& cls, Effort eff);
  private: public:
    virtual ValidatedLowerKleenean _separated(const BasicSetType& bx) const = 0;
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for compact (closed and bounded) sets.
template<class T> class CompactSetInterface<EffectiveTag,T>
    : public virtual BoundedSetInterface<EffectiveTag,T>,
      public virtual ClosedSetInterface<EffectiveTag,T>
{
    using P=EffectiveTag;
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    inline CompactSetInterface<P,T>* _copy() const { return this->clone(); }
    virtual CompactSetInterface<P,T>* clone() const = 0;

    //! \brief Tests if the set is inside \a ops.
    LowerKleenean inside(const OpenSetInterface<P,T>& ops) const;
        ValidatedLowerKleenean inside(const OpenSetInterface<P,T>& ops, Effort eff) const {
            return SetOperations<T>::inside(*this,ops,eff); }

    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a acc.
    LowerKleenean separated(const ClosedSetInterface<P,T>& cls) const;
    ValidatedLowerKleenean separated(const ClosedSetInterface<P,T>& cls, Effort eff) const {
        return SetOperations<T>::separated(*this,cls,eff); }

    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const {
        return this->_inside(bx); }
    inline ValidatedLowerKleenean inside(const BasicSetType& bx, Effort eff) const {
        return this->_inside(bx); }
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const {
        return this->_separated(bx); }

    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a acc.
    friend LowerKleenean inside(const CompactSetInterface<P,T>& cps, const OpenSetInterface<P,T>& ops) {
        return cps.inside(ops); }
    friend ValidatedLowerKleenean inside(const CompactSetInterface<P,T>& cps, const OpenSetInterface<P,T>& ops, Effort eff) {
        return cps.inside(ops,eff); }
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a acc.
    friend LowerKleenean separated(const CompactSetInterface<P,T>& cps, const ClosedSetInterface<P,T>& cls) {
        return cps.separated(cls); }
    friend ValidatedLowerKleenean separated(const CompactSetInterface<P,T>& cps, const ClosedSetInterface<P,T>& cls, Effort eff) {
        return cps.separated(cls,eff); }
  private: public:
    virtual ValidatedLowerKleenean _inside(const BasicSetType& bx) const = 0;
    virtual ValidatedLowerKleenean _separated(const BasicSetType& bx) const = 0;
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class RegularSetInterface<EffectiveTag,T>
    : public virtual OpenSetInterface<EffectiveTag,T>,
      public virtual ClosedSetInterface<EffectiveTag,T>
{
    using P=EffectiveTag;
  public:
    using typename SetInterfaceBase<T>::ElementType;
    using typename SetInterfaceBase<T>::BasicSetType;

    inline RegularSetInterface<P,T>* _copy() const { return this->clone(); }
    virtual RegularSetInterface<P,T>* clone() const = 0;
    //! \brief Tests if the set contains \a pt.
    Kleenean contains(const ElementType& pt) const;
    ValidatedKleenean contains(const ElementType& pt, Effort eff) const {
        return SetOperations<T>::contains(*this,pt,eff); }
    ValidatedKleenean contains(const ElementType& pt, Accuracy acc) const {
        return SetOperations<T>::contains(*this,pt,acc); }
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a acc.
    friend Kleenean overlap(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs);
    friend ValidatedKleenean overlap(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs, Effort eff);
    friend ValidatedKleenean overlap(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs, Accuracy acc);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a acc.
    friend Kleenean inside(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs);
    friend ValidatedKleenean inside(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs, Effort eff);
    friend ValidatedKleenean inside(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs, Accuracy acc);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a acc.
    friend Kleenean separated(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs);
    friend ValidatedKleenean separated(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs, Effort eff);
    friend ValidatedKleenean separated(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs, Accuracy acc);
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for located (overt and compact) sets.
template<class T> class LocatedSetInterface<EffectiveTag,T>
    : public virtual OvertSetInterface<EffectiveTag,T>,
      public virtual CompactSetInterface<EffectiveTag,T>
{
    using P=EffectiveTag;
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    inline LocatedSetInterface<P,T>* _copy() const { return this->clone(); }
    virtual LocatedSetInterface<P,T>* clone() const = 0;

    using OvertSetInterface<P,T>::overlaps;
    using CompactSetInterface<P,T>::inside;
    using CompactSetInterface<P,T>::separated;

    //! \brief Tests if the set overlaps (intersects) the set \a rs.
    Kleenean overlaps(const RegularSetInterface<P,T>& rs) const;
    ValidatedKleenean overlaps(const RegularSetInterface<P,T>& rs, Effort eff) const {
        return SetOperations<T>::overlap(*this,rs,eff); }
    ValidatedKleenean overlaps(const RegularSetInterface<P,T>& rs, Accuracy acc) const {
        return SetOperations<T>::overlap(*this,rs,acc); }
    virtual ValidatedLowerKleenean _overlaps(const BasicSetType& bx) const = 0;
    //! \brief Tests if the set is inside (a subset of) the set \a rs.
    Kleenean inside(const RegularSetInterface<P,T>& rs) const;
    ValidatedKleenean inside(const RegularSetInterface<P,T>& rs, Effort eff) const {
        return SetOperations<T>::inside(*this,rs,eff); }
    ValidatedKleenean inside(const RegularSetInterface<P,T>& rs, Accuracy acc) const {
        return SetOperations<T>::inside(*this,rs,acc); }
    virtual ValidatedLowerKleenean _inside(const BasicSetType& bx) const = 0;
    //! \brief Tests if the set is separated (disjoint) from the set \a rs.
    Kleenean separated(const RegularSetInterface<P,T>& rs) const;
    ValidatedKleenean separated(const RegularSetInterface<P,T>& rs, Effort eff) const {
        return SetOperations<T>::separated(*this,rs,eff); }
    ValidatedKleenean separated(const RegularSetInterface<P,T>& rs, Accuracy acc) const{
        return SetOperations<T>::separated(*this,rs,acc); }
    virtual ValidatedLowerKleenean _separated(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a acc.
    friend Kleenean overlap(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs) {
        return ls.overlaps(rs); }
    friend ValidatedKleenean overlap(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs, Effort eff) {
        return ls.overlaps(rs,eff); }
    friend ValidatedKleenean overlap(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs, Accuracy acc) {
        return ls.overlaps(rs,acc); }
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a acc.
    friend Kleenean inside(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs) {
        return ls.inside(rs); }
    friend ValidatedKleenean inside(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs, Effort eff) {
        return ls.inside(rs,eff); }
    friend ValidatedKleenean inside(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs, Accuracy acc) {
        return ls.inside(rs,acc); }
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a acc.
    friend Kleenean separated(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs) {
        return ls.separated(rs); }
    friend ValidatedKleenean separated(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs, Effort eff) {
        return ls.separated(rs,eff); }
    friend ValidatedKleenean separated(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs, Accuracy acc) {
        return ls.separated(rs,acc); }
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Complete set interface for bounded regular sets.
template<class T> class RegularLocatedSetInterface<EffectiveTag,T>
    : public virtual RegularSetInterface<EffectiveTag,T>,
      public virtual LocatedSetInterface<EffectiveTag,T>
{
    using P=EffectiveTag;
  public:
    inline RegularLocatedSetInterface<P,T>* _copy() const { return this->clone(); }
    virtual RegularLocatedSetInterface<P,T>* clone() const = 0;
};


inline OutputStream& operator<<(OutputStream& os, const WritableInterface& w);

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

    virtual BoundedSetInterface<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set is a inside of \a bx.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const {
        return this->_inside(bx); }
    //! \brief Returns a bounding box for the set.
    virtual UpperBoxType bounding_box() const = 0;
  private: public:
    virtual ValidatedLowerKleenean _inside(const BasicSetType& bx) const = 0;
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for overt sets, for which intersection with an open box is verifiable.
template<class T> class OvertSetInterface<ValidatedTag,T>
    : public virtual SetInterfaceBase<T>
{
    using P=ValidatedTag;
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual OvertSetInterface<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set overlaps \a ops.
    ValidatedLowerKleenean overlaps(const OpenSetInterface<ValidatedTag,T>& ops) const;
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const {
        return this->_overlaps(bx); }
    //! \brief Tests if \a ovs overlaps \a ops.
    friend ValidatedLowerKleenean overlap(const OvertSetInterface<ValidatedTag,T>& ovs, const OpenSetInterface<ValidatedTag,T>& ops);
  private: public:
    virtual ValidatedLowerKleenean _overlaps(const BasicSetType& bx) const = 0;
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for open sets.
template<class T> class OpenSetInterface<ValidatedTag,T>
    : public virtual OvertSetInterface<ValidatedTag,T>
{
    using P=ValidatedTag;
  public:
    using typename SetInterfaceBase<T>::ValidatedElementType;
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual OpenSetInterface<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set contains \a pt.
    ValidatedLowerKleenean contains(const ValidatedElementType& pt) const;
    //! \brief Tests if the set covers \a bx.
    inline ValidatedLowerKleenean covers(const BasicSetType& bx) const {
        return this->_covers(bx); }
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a acc.
    friend ValidatedLowerKleenean overlap(const OvertSetInterface<ValidatedTag,T>& ovs, const OpenSetInterface<ValidatedTag,T>& ops);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a acc.
    friend ValidatedLowerKleenean inside(const CompactSetInterface<ValidatedTag,T>& ls, const OpenSetInterface<ValidatedTag,T>& rs);
  private: public:
    virtual ValidatedLowerKleenean _covers(const BasicSetType& bx) const = 0;
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for closed sets.
template<class T> class ClosedSetInterface<ValidatedTag,T>
    : public virtual SetInterfaceBase<T>
{
    using P=ValidatedTag;
  public:
    using typename SetInterfaceBase<T>::ValidatedElementType;
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual ClosedSetInterface<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set (does not) contain \a pt.
    ValidatedUpperKleenean contains(const ValidatedElementType& pt) const;
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const {
        return this->_separated(bx); }
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a acc.
    friend ValidatedLowerKleenean separated(const CompactSetInterface<ValidatedTag,T>& cps, const ClosedSetInterface<ValidatedTag,T>& cls);
  private: public:
    virtual ValidatedLowerKleenean _separated(const BasicSetType& bx) const = 0;
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for compact (closed and bounded) sets.
template<class T> class CompactSetInterface<ValidatedTag,T>
    : public virtual BoundedSetInterface<ValidatedTag,T>,
      public virtual ClosedSetInterface<ValidatedTag,T>
{
    using P=ValidatedTag;
  public:
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual CompactSetInterface<ValidatedTag,T>* clone() const = 0;
    using BoundedSetInterface<ValidatedTag,T>::inside;
    using ClosedSetInterface<ValidatedTag,T>::separated;
    //! \brief Tests if the set is inside (a subset of) \a ops.
    ValidatedLowerKleenean inside(const OpenSetInterface<ValidatedTag,T>& ops) const;
    //! \brief Tests if the set is separated from \a cls.
    ValidatedLowerKleenean separated(const ClosedSetInterface<ValidatedTag,T>& cls) const;
    //virtual ValidatedSierpinskian empty() const = 0;
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a acc.
    friend ValidatedLowerKleenean inside(const CompactSetInterface<ValidatedTag,T>& ls, const OpenSetInterface<ValidatedTag,T>& rs);
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a acc.
    friend ValidatedLowerKleenean separated(const CompactSetInterface<ValidatedTag,T>& cps, const ClosedSetInterface<ValidatedTag,T>& cls);

};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class RegularSetInterface<ValidatedTag,T>
    : public virtual OpenSetInterface<ValidatedTag,T>,
      public virtual ClosedSetInterface<ValidatedTag,T>
{
    using P=ValidatedTag;
  public:
    using typename SetInterfaceBase<T>::ValidatedElementType;
    using typename SetInterfaceBase<T>::BasicSetType;

    virtual RegularSetInterface<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set contains \a pt.
    ValidatedKleenean contains(const ValidatedElementType& pt) const;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a acc.
    friend ValidatedKleenean overlap(const LocatedSetInterface<ValidatedTag,T>& ls, const RegularSetInterface<ValidatedTag,T>& rs);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a acc.
    friend ValidatedKleenean inside(const LocatedSetInterface<ValidatedTag,T>& ls, const RegularSetInterface<ValidatedTag,T>& rs);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a acc.
    friend ValidatedKleenean separated(const LocatedSetInterface<ValidatedTag,T>& ls, const RegularSetInterface<ValidatedTag,T>& rs);
};


//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Interface for located (overt and compact) sets.
template<class T> class LocatedSetInterface<ValidatedTag,T>
    : public virtual OvertSetInterface<ValidatedTag,T>,
      public virtual CompactSetInterface<ValidatedTag,T>
{
    using P=ValidatedTag;
  public:
    using typename CompactSetInterface<P,T>::BasicSetType;

    virtual LocatedSetInterface<P,T>* clone() const = 0;
    using OvertSetInterface<P,T>::overlaps;
    using CompactSetInterface<P,T>::inside;
    using CompactSetInterface<P,T>::separated;
    //! \brief Tests if the set is inside (a subset of) \a rs.
    ValidatedKleenean inside(const RegularSetInterface<P,T>& rs) const;
    //! \brief Tests if the set overlaps (intersects) the set \a rs.
    ValidatedKleenean overlaps(const RegularSetInterface<P,T>& rs) const;
    //! \brief Tests if the set is separated (disjoint) from the set \a rs.
    ValidatedKleenean separated(const RegularSetInterface<P,T>& rs) const;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a acc.
    friend ValidatedKleenean overlap(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a acc.
    friend ValidatedKleenean inside(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a acc.
    friend ValidatedKleenean separated(const LocatedSetInterface<P,T>& ls, const RegularSetInterface<P,T>& rs);
};

//! \ingroup GeometryModule SetInterfaceSubModule
//! \brief Complete set interface for bounded regular sets.
template<class T> class RegularLocatedSetInterface<ValidatedTag,T>
    : public virtual RegularSetInterface<ValidatedTag,T>,
      public virtual LocatedSetInterface<ValidatedTag,T>
{
    using P=ValidatedTag;
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
