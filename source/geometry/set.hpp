/***************************************************************************
 *            geometry/set.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file geometry/set.hpp
 *  \brief Handles for open, closed, overt and compact subsets of Euclidean space.
 */

#ifndef ARIADNE_SET_HPP
#define ARIADNE_SET_HPP

#include <iosfwd>

#include "../utility/tribool.hpp"
#include "../geometry/box.hpp"
#include "../utility/handle.hpp"
#include "../geometry/set.decl.hpp"
#include "../geometry/set_interface.hpp"

namespace Ariadne {

//! \ingroup GeometryModule SetSubModule
//! \brief Base handle class for sets.
template<class T> class SetBase
    : public Handle<SetInterfaceBase<T>>
{
  public:
    template<class ...Args> SetBase(Args&&... args) : Handle<SetInterfaceBase<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
};


//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded sets.
template<class T> class BoundedSet<EffectiveTag,T>
    : public Handle<BoundedSetInterface<EffectiveTag,T>>
{
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
  public:
    template<class ...Args> BoundedSet(Args&&... args) : Handle<BoundedSetInterface<EffectiveTag,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is a inside of \a bx.
    //! A set \a A is \em inside \a B if the closure of \a A is a subset of the interior of \a B.
    //! A set \f$A\f$ is \em inside \f$B\f$ if \f$\,\overline{\!A} \subset B^\circ\f$.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    //! If the set is empty, then the first component of the result should be empty.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};


//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for overt sets, for which intersection with an open box is verifiable.
template<class T> class OvertSet<EffectiveTag,T>
    : public Handle<OvertSetInterface<EffectiveTag,T>>
{
    using P=EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
  public:
    template<typename ...Args> OvertSet(Args&&... args) : Handle<OvertSetInterface<EffectiveTag,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a ops.
    inline LowerKleenean overlaps(const OpenSet<P,T>& ops) const;

    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const {
        return this->reference().overlaps(bx); }
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx, Effort eff) const {
        return this->reference().overlaps(bx,eff); }
/*
    friend ValidatedLowerKleenean overlap(const OvertSet<P,T>& ovs, const OpenSet<P,T>& ops, Effort eff) {
        return ovs.reference().overlaps(ops.reference(),eff); }
*/
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for open sets.
template<class T> class OpenSet<EffectiveTag,T>
    : public Handle<OpenSetInterface<EffectiveTag,T>>
{
    using P=EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::ElementType ElementType;
    typedef OpenSet<P,T> OpenSetType;
  public:
    template<typename ...Args> OpenSet(Args&&... args) : Handle<OpenSetInterface<EffectiveTag,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set contains \a pt.
    inline LowerKleenean contains(const ElementType& pt) const;
    //! \brief Tests if the set overlaps \a ops.
    inline LowerKleenean overlaps(const OpenSetType& ops) const;

    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx, Effort eff) const {
        return this->reference().overlaps(bx,eff); }
    //! \brief Tests if the set covers of \a bx.
    //! A set \a A \em covers \a B if the interiors of \a A is a superset of the closure of \a B.
    //! A set \f$A\f$ \em covers \f$B\f$ if \f$A^\circ \supset \overline{B}\f$.
    inline ValidatedLowerKleenean covers(const BasicSetType& bx) const { return this->reference().covers(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for closed sets.
template<class T> class ClosedSet<EffectiveTag,T>
    : public Handle<ClosedSetInterface<EffectiveTag,T>>
{
    using P=EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::ElementType ElementType;
  public:
    template<typename ...Args> ClosedSet(Args&&... args) : Handle<ClosedSetInterface<EffectiveTag,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set (does not) contain \a pt.
    UpperKleenean contains(const ElementType& pt) const;
    //! \brief Tests if the set is separated from \a bx.
    //! A set \a A is \em separated from \a B if the closures of \a A and \a B are disjoint.
    //! A set \f$A\f$ is \em separated from \f$B\f$ if \f$\,\overline{\!A} \cap \overline{B} = \emptyset\f$.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for compact (closed and bounded) sets.
template<class T> class CompactSet<EffectiveTag,T>
    : public Handle<CompactSetInterface<EffectiveTag,T>>
{
    using P=EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::ElementType ElementType;
    typedef OpenSet<P,T> OpenSetType;
    typedef ClosedSet<P,T> ClosedSetType;
  public:
    template<typename ...Args> CompactSet(Args&&... args) : Handle<CompactSetInterface<EffectiveTag,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }

    //! \brief Tests if the set is inside (a subset of) \a ops.
    LowerKleenean inside(const OpenSetType& ops) const;
    //! \brief Tests if the set is separated (disjoint) from \a cls.
    LowerKleenean separated(const ClosedSetType& cls) const;

    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    inline ValidatedLowerKleenean inside(const BasicSetType& bx, Effort eff) const { return this->reference().inside(bx,eff); }
    //! \brief Returns a bounding box for the set.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
/*
    friend ValidatedLowerKleenean inside(const CompactSet<P,T>& cps, const ClosedSet<P,T>& ops, Effort eff) {
        return cps.reference().inside(ops.reference(),eff); }
    friend ValidatedLowerKleenean separated(const CompactSet<P,T>& cps, const ClosedSet<P,T>& cls, Effort eff) {
        return cps.reference().separated(cls.reference(),eff); }
*/
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class RegularSet<EffectiveTag,T>
    : public Handle<RegularSetInterface<EffectiveTag,T>>
{
    using P=EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::ElementType ElementType;
  public:
    template<typename ...Args> RegularSet(Args&&... args) : Handle<RegularSetInterface<EffectiveTag,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set (does not) contain \a pt.
    Kleenean contains(const ElementType& pt) const;
    ValidatedKleenean contains(const ElementType& pt, Effort eff) const {
        return SetOperations<T>::contains(*this,pt,eff); }
    ValidatedKleenean contains(const ElementType& pt, Accuracy acc) const {
        return SetOperations<T>::contains(*this,pt,acc); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx, Effort eff) const { return this->reference().overlaps(bx,eff); }
    //! \brief Tests if the set covers of \a bx.
    inline ValidatedLowerKleenean covers(const BasicSetType& bx) const { return this->reference().covers(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for located (overt and compact) sets.
template<class T> class LocatedSet<EffectiveTag,T>
    : public Handle<LocatedSetInterface<EffectiveTag,T>>
{
    using P=EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::ElementType ElementType;
  public:
    template<typename ...Args> LocatedSet(Args&&... args) : Handle<LocatedSetInterface<EffectiveTag,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a ops.
    LowerKleenean overlaps(const OpenSet<P,T>& ops) const;
    //! \brief Tests if the set is a inside of \a ops.
    LowerKleenean inside(const OpenSet<P,T>& ops) const;
    //! \brief Tests if the set is separated from \a cls.
    LowerKleenean separated(const ClosedSet<P,T>& cls) const;
    //! \brief Tests if the set overlaps \a rs.
    Kleenean overlaps(const RegularSet<P,T>& rs) const;
    //! \brief Tests if the set is a inside of \a rs.
    Kleenean inside(const RegularSet<P,T>& rs) const;
    //! \brief Tests if the set is separated from \a rs.
    Kleenean separated(const RegularSet<P,T>& rs) const;
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx, Effort eff) const { return this->reference().overlaps(bx,eff); }
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    inline ValidatedLowerKleenean inside(const BasicSetType& bx, Effort eff) const { return this->reference().inside(bx,eff); }
    //! \brief Returns a bounding box for the set.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
/*
    friend ValidatedKleenean overlap(const LocatedSet<P,T>& ls, const RegularSet<P,T>& rs, Effort eff) {
        return ls.reference().overlaps(rs.reference(),eff); }
    friend ValidatedKleenean inside(const LocatedSet<P,T>& ls, const RegularSet<P,T>& rs, Effort eff) {
        return ls.reference().inside(rs.reference(),eff); }
    friend ValidatedKleenean separated(const LocatedSet<P,T>& ls, const RegularSet<P,T>& cls, Effort eff) {
        return ls.reference().separated(cls.reference(),eff); }
*/
    friend ValidatedKleenean overlap(const LocatedSet<P,T>& ls, const RegularSet<P,T>& rs, Accuracy acc) {
        return ls.reference().overlaps(rs.reference(),acc); }
    friend ValidatedKleenean inside(const LocatedSet<P,T>& ls, const RegularSet<P,T>& rs, Accuracy acc) {
        return ls.reference().inside(rs.reference(),acc); }
    friend ValidatedKleenean separated(const LocatedSet<P,T>& ls, const RegularSet<P,T>& cls, Accuracy acc) {
        return ls.reference().separated(cls.reference(),acc); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded regular sets.
template<class T> class RegularLocatedSet<EffectiveTag,T>
    : public Handle<RegularLocatedSetInterface<EffectiveTag,T>>
{
    using P=EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::ElementType ElementType;
  public:
    template<typename ...Args> RegularLocatedSet(Args&&... args) : Handle<RegularLocatedSetInterface<P,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx, Effort eff) const { return this->reference().overlaps(bx,eff); }
    //! \brief Tests if the set covers of \a bx.
    inline ValidatedLowerKleenean covers(const BasicSetType& bx) const { return this->reference().covers(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    inline ValidatedLowerKleenean inside(const BasicSetType& bx, Effort eff) const { return this->reference().inside(bx,eff); }
    //! \brief Returns a bounding box for the set.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};






//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded sets.
template<class T> class BoundedSet<ValidatedTag,T>
    : public Handle<BoundedSetInterface<ValidatedTag,T>>
{
    using P=ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
  public:
    template<class ...Args> BoundedSet(Args&&... args) : Handle<BoundedSetInterface<P,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is a inside of \a bx.
    //! A set \a A is \em inside \a B if the closure of \a A is a subset of the interior of \a B.
    //! A set \f$A\f$ is \em inside \f$B\f$ if \f$\,\overline{\!A} \subset B^\circ\f$.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    //! If the set is empty, then the first component of the result should be empty.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};


//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for overt sets, for which intersection with an open box is verifiable.
template<class T> class OvertSet<ValidatedTag,T>
    : public Handle<OvertSetInterface<ValidatedTag,T>>
{
    using P=ValidatedTag;
  public:
    typedef OpenSet<P,T> OpenSetType;
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
  public:
    template<typename ...Args> OvertSet(Args&&... args) : Handle<OvertSetInterface<P,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a ops.
    inline ValidatedLowerKleenean overlaps(const OpenSetType& ops) const { return this->reference().overlaps(ops); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for open sets.
template<class T> class OpenSet<ValidatedTag,T>
    : public Handle<OpenSetInterface<ValidatedTag,T>>
{
    using P=ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::ValidatedElementType ValidatedElementType;
  public:
    template<typename ...Args> OpenSet(Args&&... args) : Handle<OpenSetInterface<ValidatedTag,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set contains \a pt.
    ValidatedLowerKleenean contains(const ValidatedElementType& pt) const;
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set covers of \a bx.
    //! A set \a A \em covers \a B if the interiors of \a A is a superset of the closure of \a B.
    //! A set \f$A\f$ \em covers \f$B\f$ if \f$A^\circ \supset \overline{B}\f$.
    inline ValidatedLowerKleenean covers(const BasicSetType& bx) const { return this->reference().covers(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for closed sets.
template<class T> class ClosedSet<ValidatedTag,T>
    : public Handle<ClosedSetInterface<ValidatedTag,T>>
{
    using P=ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::ElementType ElementType;
    typedef typename SetTraits<T>::ValidatedElementType ValidatedElementType;
  public:
    template<typename ...Args> ClosedSet(Args&&... args) : Handle<ClosedSetInterface<ValidatedTag,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set contains \a pt.
    ValidatedUpperKleenean contains(const ValidatedElementType& pt) const;
    //! \brief Tests if the set is separated from \a bx.
    //! A set \a A is \em separated from \a B if the closures of \a A and \a B are disjoint.
    //! A set \f$A\f$ is \em separated from \f$B\f$ if \f$\,\overline{\!A} \cap \overline{B} = \emptyset\f$.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for compact (closed and bounded) sets.
template<class T> class CompactSet<ValidatedTag,T>
    : public Handle<CompactSetInterface<ValidatedTag,T>>
{
    using P=ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::ElementType ElementType;
    typedef OpenSet<P,T> OpenSetType;
    typedef ClosedSet<P,T> ClosedSetType;
  public:
    template<typename ...Args> CompactSet(Args&&... args) : Handle<CompactSetInterface<P,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is separated from \a bx.
    ValidatedLowerKleenean separated(const ClosedSetType& cls) const;
    //! \brief Tests if the set is a inside of \a bx.
    ValidatedLowerKleenean inside(const OpenSetType& ops) const;
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class RegularSet<ValidatedTag,T>
    : public Handle<RegularSetInterface<ValidatedTag,T>>
{
    using P=ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::ElementType ElementType;
    typedef typename SetTraits<T>::ValidatedElementType ValidatedElementType;
  public:
    template<typename ...Args> RegularSet(Args&&... args) : Handle<RegularSetInterface<P,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set contains \a pt.
    ValidatedKleenean contains(const ValidatedElementType& pt) const;
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set covers of \a bx.
    inline ValidatedLowerKleenean covers(const BasicSetType& bx) const { return this->reference().covers(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for located (overt and compact) sets.
template<class T> class LocatedSet<ValidatedTag,T>
    : public Handle<LocatedSetInterface<ValidatedTag,T>>
{
    using P=ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::ElementType ElementType;
  public:
    template<typename ...Args> LocatedSet(Args&&... args) : Handle<LocatedSetInterface<P,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded regular sets.
template<class T> class RegularLocatedSet<ValidatedTag,T>
    : public Handle<RegularLocatedSetInterface<ValidatedTag,T>>
{
    using P=ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::ElementType ElementType;
  public:
    template<typename ...Args> RegularLocatedSet(Args&&... args) : Handle<RegularLocatedSetInterface<P,T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set covers of \a bx.
    inline ValidatedLowerKleenean covers(const BasicSetType& bx) const { return this->reference().covers(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};

} // namespace Ariadne


#endif // ARIADNE_SET_HPP
