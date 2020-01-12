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

//! \ingroup GeometryModule SetSubModule
//! \brief Base handle class for sets.
class SetBase
    : public Handle<SetInterfaceBase>
{
  public:
    template<class ...Args> SetBase(Args&&... args) : Handle<SetInterfaceBase>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded sets.
class BoundedSet
    : public Handle<BoundedSetInterface>
{
  public:
    template<class ...Args> BoundedSet(Args&&... args) : Handle<BoundedSetInterface>(std::forward<Args>(args)...) { }
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
class OvertSet
    : public Handle<OvertSetInterface>
{
  public:
    template<typename ...Args> OvertSet(Args&&... args) : Handle<OvertSetInterface>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline LowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for open sets.
class OpenSet
    : public Handle<OpenSetInterface>
{
  public:
    template<typename ...Args> OpenSet(Args&&... args) : Handle<OpenSetInterface>(std::forward<Args>(args)...) { }
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
class ClosedSet
    : public Handle<ClosedSetInterface>
{
  public:
    template<typename ...Args> ClosedSet(Args&&... args) : Handle<ClosedSetInterface>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is separated from \a bx.
    //! A set \a A is \em separated from \a B if the closures of \a A and \a B are disjoint.
    //! A set \f$A\f$ is \em separated from \f$B\f$ if \f$\,\overline{\!A} \cap \overline{B} = \emptyset\f$.
    inline LowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for compact (closed and bounded) sets.
class CompactSet
    : public Handle<CompactSetInterface>
{
  public:
    template<typename ...Args> CompactSet(Args&&... args) : Handle<CompactSetInterface>(std::forward<Args>(args)...) { }
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
class RegularSet
    : public Handle<RegularSetInterface>
{
  public:
    template<typename ...Args> RegularSet(Args&&... args) : Handle<RegularSetInterface>(std::forward<Args>(args)...) { }
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
class LocatedSet
    : public Handle<LocatedSetInterface>
{
  public:
    template<typename ...Args> LocatedSet(Args&&... args) : Handle<LocatedSetInterface>(std::forward<Args>(args)...) { }
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
class RegularLocatedSet
    : public Handle<RegularLocatedSetInterface>
{
  public:
    template<typename ...Args> RegularLocatedSet(Args&&... args) : Handle<RegularLocatedSetInterface>(std::forward<Args>(args)...) { }
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
class ValidatedBoundedSet
    : public Handle<ValidatedBoundedSetInterface>
{
  public:
    template<class ...Args> ValidatedBoundedSet(Args&&... args) : Handle<ValidatedBoundedSetInterface>(std::forward<Args>(args)...) { }
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
class ValidatedOvertSet
    : public Handle<ValidatedOvertSetInterface>
{
  public:
    template<typename ...Args> ValidatedOvertSet(Args&&... args) : Handle<ValidatedOvertSetInterface>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set overlaps \a bx.
    inline ValidatedLowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for open sets.
class ValidatedOpenSet
    : public Handle<ValidatedOpenSetInterface>
{
  public:
    template<typename ...Args> ValidatedOpenSet(Args&&... args) : Handle<ValidatedOpenSetInterface>(std::forward<Args>(args)...) { }
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
class ValidatedClosedSet
    : public Handle<ValidatedClosedSetInterface>
{
  public:
    template<typename ...Args> ValidatedClosedSet(Args&&... args) : Handle<ValidatedClosedSetInterface>(std::forward<Args>(args)...) { }
    //! \brief The dimension of the set.
    inline DimensionType dimension() const { return this->reference().dimension(); }
    //! \brief Tests if the set is separated from \a bx.
    //! A set \a A is \em separated from \a B if the closures of \a A and \a B are disjoint.
    //! A set \f$A\f$ is \em separated from \f$B\f$ if \f$\,\overline{\!A} \cap \overline{B} = \emptyset\f$.
    inline ValidatedLowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for compact (closed and bounded) sets.
class ValidatedCompactSet
    : public Handle<ValidatedCompactSetInterface>
{
  public:
    template<typename ...Args> ValidatedCompactSet(Args&&... args) : Handle<ValidatedCompactSetInterface>(std::forward<Args>(args)...) { }
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
class ValidatedRegularSet
    : public Handle<ValidatedRegularSetInterface>
{
  public:
    template<typename ...Args> ValidatedRegularSet(Args&&... args) : Handle<ValidatedRegularSetInterface>(std::forward<Args>(args)...) { }
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
class ValidatedLocatedSet
    : public Handle<ValidatedLocatedSetInterface>
{
  public:
    template<typename ...Args> ValidatedLocatedSet(Args&&... args) : Handle<ValidatedLocatedSetInterface>(std::forward<Args>(args)...) { }
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
class ValidatedRegularLocatedSet
    : public Handle<ValidatedRegularLocatedSetInterface>
{
  public:
    template<typename ...Args> ValidatedRegularLocatedSet(Args&&... args) : Handle<ValidatedRegularLocatedSetInterface>(std::forward<Args>(args)...) { }
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
