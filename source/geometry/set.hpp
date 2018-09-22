/***************************************************************************
 *            set.hpp
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file set.hpp
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
    inline DimensionType dimension() const { return this->reference().dimension(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for singleton sets.
class BoundedSet
    : public Handle<BoundedSetInterface>
{
  public:
    template<class ...Args> BoundedSet(Args&&... args) : Handle<BoundedSetInterface>(std::forward<Args>(args)...) { }
    inline DimensionType dimension() const { return this->reference().dimension(); }
    inline LowerKleenean inside(const ExactBoxType& bx) const { return this->reference().inside(bx); }
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};


//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for overt sets, for which intersection with an open box is verifiable.
class OvertSet
    : public Handle<OvertSetInterface>
{
  public:
    template<typename ...Args> OvertSet(Args&&... args) : Handle<OvertSetInterface>(std::forward<Args>(args)...) { }
    inline DimensionType dimension() const { return this->reference().dimension(); }
    inline LowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for open sets.
class OpenSet
    : public Handle<OpenSetInterface>
{
  public:
    template<typename ...Args> OpenSet(Args&&... args) : Handle<OpenSetInterface>(std::forward<Args>(args)...) { }
    inline DimensionType dimension() const { return this->reference().dimension(); }
    inline LowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
    inline LowerKleenean covers(const ExactBoxType& bx) const { return this->reference().covers(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for closed sets.
class ClosedSet
    : public Handle<ClosedSetInterface>
{
  public:
    template<typename ...Args> ClosedSet(Args&&... args) : Handle<ClosedSetInterface>(std::forward<Args>(args)...) { }
    inline DimensionType dimension() const { return this->reference().dimension(); }
    inline LowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for compact (closed and singleton) sets.
class CompactSet
    : public Handle<CompactSetInterface>
{
  public:
    template<typename ...Args> CompactSet(Args&&... args) : Handle<CompactSetInterface>(std::forward<Args>(args)...) { }
    inline DimensionType dimension() const { return this->reference().dimension(); }
    inline LowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
    inline LowerKleenean inside(const ExactBoxType& bx) const { return this->reference().inside(bx); }
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
class RegularSet
    : public Handle<RegularSetInterface>
{
  public:
    template<typename ...Args> RegularSet(Args&&... args) : Handle<RegularSetInterface>(std::forward<Args>(args)...) { }
    inline DimensionType dimension() const { return this->reference().dimension(); }
    inline LowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
    inline LowerKleenean covers(const ExactBoxType& bx) const { return this->reference().covers(bx); }
    inline LowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
};



//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for located (overt and compact) sets.
class LocatedSet
    : public Handle<LocatedSetInterface>
{
  public:
    template<typename ...Args> LocatedSet(Args&&... args) : Handle<LocatedSetInterface>(std::forward<Args>(args)...) { }
    inline DimensionType dimension() const { return this->reference().dimension(); }
    inline LowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
    inline LowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
    inline LowerKleenean inside(const ExactBoxType& bx) const { return this->reference().inside(bx); }
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for singleton regular sets.
class RegularLocatedSet
    : public Handle<RegularLocatedSetInterface>
{
  public:
    template<typename ...Args> RegularLocatedSet(Args&&... args) : Handle<RegularLocatedSetInterface>(std::forward<Args>(args)...) { }
    inline DimensionType dimension() const { return this->reference().dimension(); }
    inline LowerKleenean overlaps(const ExactBoxType& bx) const { return this->reference().overlaps(bx); }
    inline LowerKleenean covers(const ExactBoxType& bx) const { return this->reference().covers(bx); }
    inline LowerKleenean separated(const ExactBoxType& bx) const { return this->reference().separated(bx); }
    inline LowerKleenean inside(const ExactBoxType& bx) const { return this->reference().inside(bx); }
    inline UpperBoxType bounding_box() const { return this->reference().bounding_box(); }
};



} // namespace Ariadne


#endif // ARIADNE_SET_INTERFACE
