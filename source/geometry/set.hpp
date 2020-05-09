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
#include "../geometry/set_interface.hpp"

namespace Ariadne {

template<class T> class SetBase;
using EuclideanSetBase = SetBase<RealVector>;
template<class T> class BoundedSet;
template<class T> class OpenSet;
template<class T> class ClosedSet;
template<class T> class OvertSet;
template<class T> class CompactSet;
template<class T> class RegularSet;
template<class T> class LocatedSet;
template<class T> class RegularLocatedSet;

using EuclideanSetBase = SetBase<RealVector>;
using EuclideanBoundedSet = BoundedSet<RealVector>;
using EuclideanOpenSet = OpenSet<RealVector>;
using EuclideanClosedSet = ClosedSet<RealVector>;
using EuclideanOvertSet = OvertSet<RealVector>;
using EuclideanCompactSet = CompactSet<RealVector>;
using EuclideanRegularSet = RegularSet<RealVector>;
using EuclideanLocatedSet = LocatedSet<RealVector>;
using EuclideanRegularLocatedSet = RegularLocatedSet<RealVector>;
//using EuclideanSet = EuclideanRegularLocatedSet;

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
template<class T> class BoundedSet
    : public Handle<BoundedSetInterface<T>>
{
  public:
    template<class ...Args> BoundedSet(Args&&... args) : Handle<BoundedSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is a inside of \a bx.
    //! A set \a A is \em inside \a B if the closure of \a A is a subset of the interior of \a B.
    //! A set \f$A\f$ is \em inside \f$B\f$ if \f$\,\overline{\!A} \subset B^\circ\f$.
    inline LowerKleenean inside(const ExactBoxType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    //! If the set is empty, then the first component of the result should be empty.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};


//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for overt sets, for which intersection with an open box is verifiable.
template<class T> class OvertSet
    : public Handle<OvertSetInterface<T>>
{
  public:
    template<typename ...Args> OvertSet(Args&&... args) : Handle<OvertSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline LowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for open sets.
template<class T> class OpenSet
    : public Handle<OpenSetInterface<T>>
{
  public:
    template<typename ...Args> OpenSet(Args&&... args) : Handle<OpenSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline LowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set covers of \a bx.
    //! A set \a A \em covers \a B if the interiors of \a A is a superset of the closure of \a B.
    //! A set \f$A\f$ \em covers \f$B\f$ if \f$A^\circ \supset \overline{B}\f$.
    inline LowerKleenean covers(const ExactBoxType& bx) const { return this->reference().covers(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for closed sets.
template<class T> class ClosedSet
    : public Handle<ClosedSetInterface<T>>
{
  public:
    template<typename ...Args> ClosedSet(Args&&... args) : Handle<ClosedSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is separated from \a bx.
    //! A set \a A is \em separated from \a B if the closures of \a A and \a B are disjoint.
    //! A set \f$A\f$ is \em separated from \f$B\f$ if \f$\,\overline{\!A} \cap \overline{B} = \emptyset\f$.
    inline LowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for compact (closed and bounded) sets.
template<class T> class CompactSet
    : public Handle<CompactSetInterface<T>>
{
  public:
    template<typename ...Args> CompactSet(Args&&... args) : Handle<CompactSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is separated from \a bx.
    inline LowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline LowerKleenean inside(const ExactBoxType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class RegularSet
    : public Handle<RegularSetInterface<T>>
{
  public:
    template<typename ...Args> RegularSet(Args&&... args) : Handle<RegularSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline LowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set covers of \a bx.
    inline LowerKleenean covers(const ExactBoxType& bx) const { return this->reference().covers(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline LowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for located (overt and compact) sets.
template<class T> class LocatedSet
    : public Handle<LocatedSetInterface<T>>
{
  public:
    template<typename ...Args> LocatedSet(Args&&... args) : Handle<LocatedSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline LowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline LowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline LowerKleenean inside(const ExactBoxType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded regular sets.
template<class T> class RegularLocatedSet
    : public Handle<RegularLocatedSetInterface<T>>
{
  public:
    template<typename ...Args> RegularLocatedSet(Args&&... args) : Handle<RegularLocatedSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline LowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set covers of \a bx.
    inline LowerKleenean covers(const ExactBoxType& bx) const { return this->reference().covers(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline LowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline LowerKleenean inside(const ExactBoxType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};






//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded sets.
template<class T> class ValidatedBoundedSet
    : public Handle<ValidatedBoundedSetInterface<T>>
{
  public:
    template<class ...Args> ValidatedBoundedSet(Args&&... args) : Handle<ValidatedBoundedSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is a inside of \a bx.
    //! A set \a A is \em inside \a B if the closure of \a A is a subset of the interior of \a B.
    //! A set \f$A\f$ is \em inside \f$B\f$ if \f$\,\overline{\!A} \subset B^\circ\f$.
    inline ValidatedLowerKleenean inside(const ExactBoxType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    //! If the set is empty, then the first component of the result should be empty.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};


//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for overt sets, for which intersection with an open box is verifiable.
template<class T> class ValidatedOvertSet
    : public Handle<ValidatedOvertSetInterface<T>>
{
  public:
    template<typename ...Args> ValidatedOvertSet(Args&&... args) : Handle<ValidatedOvertSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for open sets.
template<class T> class ValidatedOpenSet
    : public Handle<ValidatedOpenSetInterface<T>>
{
  public:
    template<typename ...Args> ValidatedOpenSet(Args&&... args) : Handle<ValidatedOpenSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set covers of \a bx.
    //! A set \a A \em covers \a B if the interiors of \a A is a superset of the closure of \a B.
    //! A set \f$A\f$ \em covers \f$B\f$ if \f$A^\circ \supset \overline{B}\f$.
    inline ValidatedLowerKleenean covers(const ExactBoxType& bx) const { return this->reference().covers(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for closed sets.
template<class T> class ValidatedClosedSet
    : public Handle<ValidatedClosedSetInterface<T>>
{
  public:
    template<typename ...Args> ValidatedClosedSet(Args&&... args) : Handle<ValidatedClosedSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is separated from \a bx.
    //! A set \a A is \em separated from \a B if the closures of \a A and \a B are disjoint.
    //! A set \f$A\f$ is \em separated from \f$B\f$ if \f$\,\overline{\!A} \cap \overline{B} = \emptyset\f$.
    inline ValidatedLowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for compact (closed and bounded) sets.
template<class T> class ValidatedCompactSet
    : public Handle<ValidatedCompactSetInterface<T>>
{
  public:
    template<typename ...Args> ValidatedCompactSet(Args&&... args) : Handle<ValidatedCompactSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline ValidatedLowerKleenean inside(const ExactBoxType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
template<class T> class ValidatedRegularSet
    : public Handle<ValidatedRegularSetInterface<T>>
{
  public:
    template<typename ...Args> ValidatedRegularSet(Args&&... args) : Handle<ValidatedRegularSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set covers of \a bx.
    inline ValidatedLowerKleenean covers(const ExactBoxType& bx) const { return this->reference().covers(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for located (overt and compact) sets.
template<class T> class ValidatedLocatedSet
    : public Handle<ValidatedLocatedSetInterface<T>>
{
  public:
    template<typename ...Args> ValidatedLocatedSet(Args&&... args) : Handle<ValidatedLocatedSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline ValidatedLowerKleenean inside(const ExactBoxType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded regular sets.
template<class T> class ValidatedRegularLocatedSet
    : public Handle<ValidatedRegularLocatedSetInterface<T>>
{
  public:
    template<typename ...Args> ValidatedRegularLocatedSet(Args&&... args) : Handle<ValidatedRegularLocatedSetInterface<T>>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
    //! \brief Tests if the set covers of \a bx.
    inline ValidatedLowerKleenean covers(const ExactBoxType& bx) const { return this->reference().covers(bx); }
    //! \brief Tests if the set is separated from \a bx.
    inline ValidatedLowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
    //! \brief Tests if the set is a inside of \a bx.
    inline ValidatedLowerKleenean inside(const ExactBoxType& bx) const { return this->reference().inside(bx); }
    //! \brief Returns a bounding box for the set.
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};

} // namespace Ariadne


#endif // ARIADNE_SET_INTERFACE
