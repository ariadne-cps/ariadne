/***************************************************************************
 *            geometry/set_mixin.hpp
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

/*! \file geometry/set_mixin.hpp
 *  \brief Mixins for open, closed, overt and compact subsets of Euclidean space.
 */

#ifndef ARIADNE_SET_MIXIN_HPP
#define ARIADNE_SET_MIXIN_HPP

#include <iosfwd>

#include "../utility/declarations.hpp"
#include "../utility/tribool.hpp"
#include "../utility/writable.hpp"
#include "../numeric/numeric.hpp"

#include "set_interface.hpp"

namespace Ariadne {

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

template<class S, class P, class T> class OvertSetMixin;
template<class S, class P, class T> class OpenSetMixin;
template<class S, class P, class T> class ClosedSetMixin;
template<class S, class P, class T> class CompactSetMixin;

template<class T> class OvertSetMixin<Void,EffectiveTag,T>
    : public OvertSetInterface<EffectiveTag,T>
{
    using P=EffectiveTag;
  public:
    LowerKleenean overlaps(const OpenSet<P,T>& ops) const;
    ValidatedLowerKleenean overlaps(const OpenSet<P,T>& ops, Effort eff) const;
    virtual ValidatedLowerKleenean overlaps(const BasicSetType& bx, Effort eff) const;
    virtual LowerKleenean overlaps(const BasicSetType& bx) const = 0;
};

template<class T> class OpenSetMixin<Void,EffectiveTag,T>
    : public OpenSetInterface<EffectiveTag,T>
{
    using P=EffectiveTag;
  public:
    LowerKleenean overlaps(const OpenSet<P,T>& ops) const;
    ValidatedLowerKleenean overlaps(const OpenSet<P,T>& ops, Effort eff) const;
    virtual ValidatedLowerKleenean overlaps(const BasicSetType& bx, Effort eff) const = 0;
};

//! \ingroup GeometryModule SetMixinSubModule
//! \brief Mixin for open sets.
template<class T> class OpenSetMixin<EffectiveTag,T>
    : public virtual OvertSetMixin<EffectiveTag,T>
{
    using P=EffectiveTag;
  public:
    using typename SetMixinBase<T>::BasicSetType;
    using typename SetMixinBase<T>::ElementType;

    virtual OpenSetMixin<EffectiveTag,T>* clone() const = 0;
    //! \brief Tests if the set contains \a pt.
    LowerKleenean contains(const ElementType& pt) const;
    ValidatedLowerKleenean contains(const ElementType& pt, Effort eff) const;
    //! \brief Tests if the set covers of \a bx.
    //! A set \a A \em covers \a B if the interiors of \a A is a superset of the closure of \a B.
    //! A set \f$A\f$ \em covers \f$B\f$ if \f$A^\circ \supset \overline{B}\f$.
    virtual LowerKleenean covers(const BasicSetType& bx) const = 0;
    virtual ValidatedLowerKleenean covers(const BasicSetType& bx, Effort eff) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend LowerKleenean overlap(const OvertSetMixin<EffectiveTag,T>& ovs, const OpenSetMixin<EffectiveTag,T>& ops);
    friend ValidatedLowerKleenean overlap(const OvertSetMixin<EffectiveTag,T>& ovs, const OpenSetMixin<EffectiveTag,T>& ops, const RawFloatDP& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend LowerKleenean inside(const CompactSetMixin<EffectiveTag,T>& ls, const OpenSetMixin<EffectiveTag,T>& rs);
    friend ValidatedLowerKleenean inside(const CompactSetMixin<EffectiveTag,T>& ls, const OpenSetMixin<EffectiveTag,T>& rs, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetMixinSubModule
//! \brief Mixin for closed sets.
template<class T> class ClosedSetMixin<EffectiveTag,T>
    : public virtual SetMixinBase<T>
{
    using P=EffectiveTag;
  public:
    using typename SetMixinBase<T>::ElementType;
    using typename SetMixinBase<T>::BasicSetType;

    virtual ClosedSetMixin<EffectiveTag,T>* clone() const = 0;
    //! \brief Tests if the set (does not) contain \a pt.
    UpperKleenean contains(const ElementType& pt) const;
    ValidatedUpperKleenean contains(const ElementType& pt, Effort eff) const;
    //! \brief Tests if the set is separated from \a bx.
    //! A set \a A is \em separated from \a B if the closures of \a A and \a B are disjoint.
    //! A set \f$A\f$ is \em separated from \f$B\f$ if \f$\,\overline{\!A} \cap \overline{B} = \emptyset\f$.
    virtual LowerKleenean separated(const BasicSetType& bx) const = 0;
    virtual ValidatedLowerKleenean separated(const BasicSetType& bx, Effort eff) const = 0;
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend LowerKleenean separated(const CompactSetMixin<EffectiveTag,T>& cps, const ClosedSetMixin<EffectiveTag,T>& cls);
    friend ValidatedLowerKleenean separated(const CompactSetMixin<EffectiveTag,T>& cps, const ClosedSetMixin<EffectiveTag,T>& cls, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetMixinSubModule
//! \brief Mixin for compact (closed and bounded) sets.
template<class T> class CompactSetMixin<EffectiveTag,T>
    : public virtual BoundedSetMixin<EffectiveTag,T>,
      public virtual ClosedSetMixin<EffectiveTag,T>
{
    using P=EffectiveTag;
  public:
    using typename SetMixinBase<T>::BasicSetType;

    virtual CompactSetMixin<P,T>* clone() const = 0;

    using BoundedSetMixin<P,T>::inside;
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eff.
    LowerKleenean inside(const OpenSet<P,T>& rs) const;
    ValidatedLowerKleenean inside(const OpenSet<P,T>& rs, Effort eff) const;

    using ClosedSetMixin<P,T>::separated;
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    LowerKleenean separated(const ClosedSet<P,T>& cls) const;
    ValidatedLowerKleenean separated(const ClosedSetMixin<P,T>& cls, Effort eff) const;

    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend LowerKleenean inside(const CompactSetMixin<P,T>& ls, const OpenSetMixin<P,T>& rs);
    friend ValidatedLowerKleenean inside(const CompactSetMixin<P,T>& ls, const OpenSetMixin<P,T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend LowerKleenean separated(const CompactSetMixin<P,T>& cps, const ClosedSetMixin<P,T>& cls);
    friend ValidatedLowerKleenean separated(const CompactSetMixin<P,T>& cps, const ClosedSetMixin<P,T>& cls, const RawFloatDP& eps);

};

//! \ingroup GeometryModule SetMixinSubModule
//! \brief Mixin for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class RegularSetMixin<EffectiveTag,T>
    : public virtual OpenSetMixin<EffectiveTag,T>,
      public virtual ClosedSetMixin<EffectiveTag,T>
{
    using P=EffectiveTag;
  public:
    using typename SetMixinBase<T>::ElementType;
    using typename SetMixinBase<T>::BasicSetType;

    virtual RegularSetMixin<P,T>* clone() const = 0;
    //! \brief Tests if the set contains \a pt.
//    virtual Kleenean contains(const ElementType& pt) const;
    ValidatedKleenean contains(const ElementType& pt, Effort eff) const;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend Kleenean overlap(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs);
    friend ValidatedKleenean overlap(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend Kleenean inside(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs);
    friend ValidatedKleenean inside(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend Kleenean separated(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs);
    friend ValidatedKleenean separated(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs, const RawFloatDP& eps);
};


//! \ingroup GeometryModule SetMixinSubModule
//! \brief Mixin for located (overt and compact) sets.
template<class T> class LocatedSetMixin<EffectiveTag,T>
    : public virtual OvertSetMixin<EffectiveTag,T>,
      public virtual CompactSetMixin<EffectiveTag,T>
{
    using P=EffectiveTag;
  public:
    using typename SetMixinBase<T>::BasicSetType;

    virtual LocatedSetMixin<P,T>* clone() const = 0;

    using OvertSetMixin<P,T>::overlaps;
    using CompactSetMixin<P,T>::inside;
    using CompactSetMixin<P,T>::separated;

    //! \brief Tests if the set overlaps (intersects) the set \a rs.
    Kleenean overlaps(const RegularSetMixin<P,T>& rs) const;
    ValidatedKleenean overlaps(const RegularSetMixin<P,T>& rs, Effort eff) const;
    //! \brief Tests if the set is inside (a subset of) the set \a rs.
    Kleenean inside(const RegularSetMixin<P,T>& rs) const;
    ValidatedKleenean inside(const RegularSetMixin<P,T>& rs, Effort eff) const;
    //! \brief Tests if the set is separated (disjoint) from the set \a rs.
    Kleenean separated(const RegularSetMixin<P,T>& rs) const;
    ValidatedKleenean separated(const RegularSetMixin<P,T>& rs, Effort eff) const;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend Kleenean overlap(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs);
    friend ValidatedKleenean overlap(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend Kleenean inside(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs);
    friend ValidatedKleenean inside(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs, const RawFloatDP& eps);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend Kleenean separated(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs);
    friend ValidatedKleenean separated(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs, const RawFloatDP& eps);
};

//! \ingroup GeometryModule SetMixinSubModule
//! \brief Complete set interface for bounded regular sets.
template<class T> class RegularLocatedSetMixin<EffectiveTag,T>
    : public virtual RegularSetMixin<EffectiveTag,T>,
      public virtual LocatedSetMixin<EffectiveTag,T>
{
  public:
    virtual RegularLocatedSetMixin<EffectiveTag,T>* clone() const = 0;
};


inline OutputStream& operator<<(OutputStream& os, const WritableMixin& w);



template<class T> class BoundedSetMixin<ValidatedTag,T>;
template<class T> class OpenSetMixin<ValidatedTag,T>;
template<class T> class ClosedSetMixin<ValidatedTag,T>;
template<class T> class OvertSetMixin<ValidatedTag,T>;
template<class T> class CompactSetMixin<ValidatedTag,T>;
template<class T> class RegularSetMixin<ValidatedTag,T>;
template<class T> class LocatedSetMixin<ValidatedTag,T>;
template<class T> class RegularLocatedSetMixin<ValidatedTag,T>;

//! \ingroup GeometryModule SetMixinSubModule
//! \brief Mixin for bounded sets.
template<class T> class BoundedSetMixin<ValidatedTag,T>
    : public virtual SetMixinBase<T> {
  public:
    using typename SetMixinBase<T>::BasicSetType;

    virtual BoundedSetMixin<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set is a inside of \a bx.
    virtual ValidatedLowerKleenean inside(const BasicSetType& bx) const = 0;
    //! \brief Returns a bounding box for the set.
    virtual UpperBoxType bounding_box() const = 0;
};

//! \ingroup GeometryModule SetMixinSubModule
//! \brief Mixin for overt sets, for which intersection with an open box is verifiable.
template<class T> class OvertSetMixin<ValidatedTag,T>
    : public virtual SetMixinBase<T>
{
    using P=ValidatedTag;
  public:
    using typename SetMixinBase<T>::BasicSetType;

    virtual OvertSetMixin<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set overlaps \a ops.
    ValidatedLowerKleenean overlaps(const OpenSetMixin<ValidatedTag,T>& ops) const;
    //! \brief Tests if the set overlaps \a bx.
    virtual ValidatedLowerKleenean overlaps(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops.
    friend ValidatedLowerKleenean overlap(const OvertSetMixin<ValidatedTag,T>& ovs, const OpenSetMixin<ValidatedTag,T>& ops);
};

//! \ingroup GeometryModule SetMixinSubModule
//! \brief Mixin for open sets.
template<class T> class OpenSetMixin<ValidatedTag,T>
    : public virtual OvertSetMixin<ValidatedTag,T>
{
    using P=ValidatedTag;
  public:
    using typename SetMixinBase<T>::ValidatedElementType;
    using typename SetMixinBase<T>::BasicSetType;

    virtual OpenSetMixin<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set contains \a pt.
    ValidatedLowerKleenean contains(const ValidatedElementType& pt) const;
    //! \brief Tests if the set covers \a bx.
    virtual ValidatedLowerKleenean covers(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
    friend ValidatedLowerKleenean overlap(const OvertSetMixin<ValidatedTag,T>& ovs, const OpenSetMixin<ValidatedTag,T>& ops);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedLowerKleenean inside(const CompactSetMixin<ValidatedTag,T>& ls, const OpenSetMixin<ValidatedTag,T>& rs);
};

//! \ingroup GeometryModule SetMixinSubModule
//! \brief Mixin for closed sets.
template<class T> class ClosedSetMixin<ValidatedTag,T>
    : public virtual SetMixinBase<T>
{
    using P=ValidatedTag;
  public:
    using typename SetMixinBase<T>::ValidatedElementType;
    using typename SetMixinBase<T>::BasicSetType;

    virtual ClosedSetMixin<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set (does not) contain \a pt.
    ValidatedUpperKleenean contains(const ValidatedElementType& pt) const;
    //! \brief Tests if the set is separated from \a bx.
    virtual ValidatedLowerKleenean separated(const BasicSetType& bx) const = 0;
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend ValidatedLowerKleenean separated(const CompactSetMixin<ValidatedTag,T>& cps, const ClosedSetMixin<ValidatedTag,T>& cls);
};

//! \ingroup GeometryModule SetMixinSubModule
//! \brief Mixin for compact (closed and bounded) sets.
template<class T> class CompactSetMixin<ValidatedTag,T>
    : public virtual BoundedSetMixin<ValidatedTag,T>,
      public virtual ClosedSetMixin<ValidatedTag,T>
{
    using P=ValidatedTag;
  public:
    using typename SetMixinBase<T>::BasicSetType;

    virtual CompactSetMixin<ValidatedTag,T>* clone() const = 0;
    using BoundedSetMixin<ValidatedTag,T>::inside;
    using ClosedSetMixin<ValidatedTag,T>::separated;
    //! \brief Tests if the set is inside (a subset of) \a ops.
    ValidatedLowerKleenean inside(const OpenSetMixin<ValidatedTag,T>& ops) const;
    //! \brief Tests if the set is separated from \a cls.
    ValidatedLowerKleenean separated(const ClosedSetMixin<ValidatedTag,T>& cls) const;
    //virtual ValidatedSierpinskian empty() const = 0;
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedLowerKleenean inside(const CompactSetMixin<ValidatedTag,T>& ls, const OpenSetMixin<ValidatedTag,T>& rs);
    //! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
    friend ValidatedLowerKleenean separated(const CompactSetMixin<ValidatedTag,T>& cps, const ClosedSetMixin<ValidatedTag,T>& cls);

};

//! \ingroup GeometryModule SetMixinSubModule
//! \brief Mixin for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class RegularSetMixin<ValidatedTag,T>
    : public virtual OpenSetMixin<ValidatedTag,T>,
      public virtual ClosedSetMixin<ValidatedTag,T>
{
    using P=ValidatedTag;
  public:
    using typename SetMixinBase<T>::ValidatedElementType;
    using typename SetMixinBase<T>::BasicSetType;

    virtual RegularSetMixin<ValidatedTag,T>* clone() const = 0;
    //! \brief Tests if the set contains \a pt.
    ValidatedKleenean contains(const ValidatedElementType& pt) const;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean overlap(const LocatedSetMixin<ValidatedTag,T>& ls, const RegularSetMixin<ValidatedTag,T>& rs);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean inside(const LocatedSetMixin<ValidatedTag,T>& ls, const RegularSetMixin<ValidatedTag,T>& rs);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean separated(const LocatedSetMixin<ValidatedTag,T>& ls, const RegularSetMixin<ValidatedTag,T>& rs);
};


//! \ingroup GeometryModule SetMixinSubModule
//! \brief Mixin for located (overt and compact) sets.
template<class T> class LocatedSetMixin<ValidatedTag,T>
    : public virtual OvertSetMixin<ValidatedTag,T>,
      public virtual CompactSetMixin<ValidatedTag,T>
{
    using P=ValidatedTag;
  public:
    using typename CompactSetMixin<P,T>::BasicSetType;

    virtual LocatedSetMixin<P,T>* clone() const = 0;
    using OvertSetMixin<P,T>::overlaps;
    using CompactSetMixin<P,T>::inside;
    using CompactSetMixin<P,T>::separated;
    //! \brief Tests if the set is inside (a subset of) \a rs.
    ValidatedKleenean inside(const RegularSetMixin<P,T>& rs) const;
    //! \brief Tests if the set overlaps (intersects) the set \a rs.
    ValidatedKleenean overlaps(const RegularSetMixin<P,T>& rs) const;
    //! \brief Tests if the set is separated (disjoint) from the set \a rs.
    ValidatedKleenean separated(const RegularSetMixin<P,T>& rs) const;
    //! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean overlap(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs);
    //! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean inside(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs);
    //! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
    friend ValidatedKleenean separated(const LocatedSetMixin<P,T>& ls, const RegularSetMixin<P,T>& rs);
};

//! \ingroup GeometryModule SetMixinSubModule
//! \brief Complete set interface for bounded regular sets.
template<class T> class RegularLocatedSetMixin<ValidatedTag,T>
    : public virtual RegularSetMixin<ValidatedTag,T>,
      public virtual LocatedSetMixin<ValidatedTag,T>
{
    using P=ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    virtual RegularLocatedSetMixin<ValidatedTag,T>* clone() const = 0;
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
    typedef BoundedSetMixin<P,T> BoundedSetMixinType;
    //! \brief The interface satisified by overt sets in the space.
    typedef OvertSetMixin<P,T> OvertSetMixinType;
    //! \brief The interface satisified by over sets in the space.
    typedef OpenSetMixin<P,T> OpenSetMixinType;
    //! \brief The interface satisified by closed sets in the space.
    typedef ClosedSetMixin<P,T> ClosedSetMixinType;
    //! \brief The interface satisified by compact sets in the space.
    typedef CompactSetMixin<P,T> CompactSetMixinType;
    //! \brief The interface satisified by regular sets in the space.
    typedef RegularSetMixin<P,T> RegularSetMixinType;
    //! \brief The interface satisified by located sets in the space.
    typedef LocatedSetMixin<P,T> LocatedSetMixinType;
    //! \brief The interface satisified by bounded regular sets.
    typedef RegularLocatedSetMixin<P,T> RegularLocatedSetMixinType;
    //! \brief The type of approximations to sets in the space.
    typedef GridTreePaving SetApproximationType;
  public:
    EuclideanSpace(const SizeType& d) : _dimension(d) { }
    const SizeType& dimension() const { return this->_dimension; }
  private:
    SizeType _dimension;
};


} // namespace Ariadne


#endif // ARIADNE_SET_MIXIN_HPP
