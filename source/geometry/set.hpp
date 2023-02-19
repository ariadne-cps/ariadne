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

#include "utility/tribool.hpp"
#include "geometry/box.hpp"
#include "utility/handle.hpp"
#include "geometry/set.decl.hpp"
#include "geometry/set_interface.hpp"

#include "geometry/geometry_concepts.hpp"

namespace Ariadne {

//! \ingroup GeometryModule SetSubModule
//! \brief Base handle class for sets.
template<class T> class SetBase
    : public Handle<SetInterfaceBase<T>>
{
  public:
    using Handle<SetInterfaceBase<T>>::Handle;
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
};


//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded sets.
template<class T> class BoundedSet<EffectiveTag,T>
    : public Handle<BoundedSetInterface<EffectiveTag,T>>
{
    using P = EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;

    using Handle<BoundedSetInterface<P,T>>::Handle;
    template<class S> requires ABoundedSet<S,P,T> BoundedSet(S const& s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is a inside of \a bx.
    //! A set \a A is \em inside \a B if the closure of \a A is a subset of the interior of \a B.
    //! A set \f$A\f$ is \em inside \f$B\f$ if \f$\,\overline{\!A} \subset B^\circ\f$.
    inline LowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    //! If the set is empty, then the first component of the result should be empty.
    inline BoundingSetType bounding_box() const { return this->reference().bounding_box(); }
};


//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for overt sets, for which intersection with an open box is verifiable.
template<class T> class OvertSet<EffectiveTag,T>
    : public Handle<OvertSetInterface<EffectiveTag,T>>
{
    using P = EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;

    using Handle<OvertSetInterface<P,T>>::Handle;
    template<class S> requires AnOvertSet<S,P,T> OvertSet(S const& s);

    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline LowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for open sets.
template<class T> class OpenSet<EffectiveTag,T>
    : public Handle<OpenSetInterface<EffectiveTag,T>>
{
    using P = EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;

    using Handle<OpenSetInterface<P,T>>::Handle;
    template<class S> requires AnOpenSet<S,P,T> OpenSet(S const& s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline LowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set covers of \a bx.
    //! A set \a A \em covers \a B if the interiors of \a A is a superset of the closure of \a B.
    //! A set \f$A\f$ \em covers \f$B\f$ if \f$A^\circ \supset \overline{B}\f$.
    inline LowerKleenean covers(const BasicSetType& bx) const { return this->reference().covers(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for closed sets.
template<class T> class ClosedSet<EffectiveTag,T>
    : public Handle<ClosedSetInterface<EffectiveTag,T>>
{
    using P = EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;

    using Handle<ClosedSetInterface<P,T>>::Handle;
    template<class S> requires AClosedSet<S,P,T> ClosedSet(S const& s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is separated from \a bx.
    //! A set \a A is \em separated from \a B if the closures of \a A and \a B are disjoint.
    //! A set \f$A\f$ is \em separated from \f$B\f$ if \f$\,\overline{\!A} \cap \overline{B} = \emptyset\f$.
    inline LowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for compact (closed and bounded) sets.
template<class T> class CompactSet<EffectiveTag,T>
    : public Handle<CompactSetInterface<EffectiveTag,T>>
{
    using P = EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;

    using Handle<CompactSetInterface<P,T>>::Handle;
    template<class S> requires ACompactSet<S,P,T> CompactSet(S const& s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is separated from \a bx.
    inline LowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline LowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline BoundingSetType bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class RegularSet<EffectiveTag,T>
    : public Handle<RegularSetInterface<EffectiveTag,T>>
{
    using P = EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;

    using Handle<RegularSetInterface<P,T>>::Handle;
    template<class S> requires ARegularSet<S,P,T> RegularSet(S const& s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline LowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set covers of \a bx.
    inline LowerKleenean covers(const BasicSetType& bx) const { return this->reference().covers(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline LowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for located (overt and compact) sets.
template<class T> class LocatedSet<EffectiveTag,T>
    : public Handle<LocatedSetInterface<EffectiveTag,T>>
{
    using P = EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;

    using Handle<LocatedSetInterface<P,T>>::Handle;
    template<class S> requires ALocatedSet<S,P,T> LocatedSet(S const& s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline LowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline LowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline LowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline BoundingSetType bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded regular sets.
template<class T> class RegularLocatedSet<EffectiveTag,T>
    : public Handle<RegularLocatedSetInterface<EffectiveTag,T>>
{
    using P = EffectiveTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;

    using Handle<RegularLocatedSetInterface<P,T>>::Handle;
    template<class S> requires ARegularLocatedSet<S,P,T> RegularLocatedSet(S const& s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline LowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set covers of \a bx.
    inline LowerKleenean covers(const BasicSetType& bx) const { return this->reference().covers(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline LowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline LowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline BoundingSetType bounding_box() const { return this->reference().bounding_box(); }
};






//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded sets.
template<class T> class BoundedSet<ValidatedTag,T>
    : public Handle<BoundedSetInterface<ValidatedTag,T>>
{
    using P = ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;

    using Handle<BoundedSetInterface<P,T>>::Handle;
    template<class S> requires ABoundedSet<S,P,T> BoundedSet(S s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is a inside of \a bx.
    //! A set \a A is \em inside \a B if the closure of \a A is a subset of the interior of \a B.
    //! A set \f$A\f$ is \em inside \f$B\f$ if \f$\,\overline{\!A} \subset B^\circ\f$.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    //! If the set is empty, then the first component of the result should be empty.
    inline BoundingSetType bounding_box() const { return this->reference().bounding_box(); }
};


//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for overt sets, for which intersection with an open box is verifiable.
template<class T> class OvertSet<ValidatedTag,T>
    : public Handle<OvertSetInterface<ValidatedTag,T>>
{
    using P = ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;

    using Handle<OvertSetInterface<P,T>>::Handle;
    template<class S> requires AnOvertSet<S,P,T> OvertSet(S s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for open sets.
template<class T> class OpenSet<ValidatedTag,T>
    : public Handle<OpenSetInterface<ValidatedTag,T>>
{
    using P = ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;

    using Handle<OpenSetInterface<P,T>>::Handle;
    template<class S> requires AnOpenSet<S,P,T> OpenSet(S s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
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
    using P = ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;

    using Handle<ClosedSetInterface<P,T>>::Handle;
    template<class S> requires AClosedSet<S,P,T> ClosedSet(S s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
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
    using P = ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;

    using Handle<CompactSetInterface<P,T>>::Handle;
    template<class S> requires ACompactSet<S,P,T> CompactSet(S s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline BoundingSetType bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class RegularSet<ValidatedTag,T>
    : public Handle<RegularSetInterface<ValidatedTag,T>>
{
    using P = ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;

    using Handle<RegularSetInterface<P,T>>::Handle;
    template<class S> requires ARegularSet<S,P,T> RegularSet(S s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
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
    using P = ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;

    using Handle<LocatedSetInterface<P,T>>::Handle;
    template<class S> requires ALocatedSet<S,P,T> LocatedSet(S s);
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const BasicSetType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const BasicSetType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline ValidatedLowerKleenean inside(const BasicSetType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline BoundingSetType bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded regular sets.
template<class T> class RegularLocatedSet<ValidatedTag,T>
    : public Handle<RegularLocatedSetInterface<ValidatedTag,T>>
{
    using P = ValidatedTag;
  public:
    typedef typename SetTraits<T>::BasicSetType BasicSetType;
    typedef typename SetTraits<T>::BoundingSetType BoundingSetType;

    using Handle<RegularLocatedSetInterface<P,T>>::Handle;
    template<class S> requires ARegularLocatedSet<S,P,T> RegularLocatedSet(S s);

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
    inline BoundingSetType bounding_box() const { return this->reference().bounding_box(); }
};

} // namespace Ariadne


#endif // ARIADNE_SET_HPP
