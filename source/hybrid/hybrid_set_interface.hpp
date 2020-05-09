/***************************************************************************
 *            hybrid/hybrid_set_interface.hpp
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

/*! \file hybrid/hybrid_set_interface.hpp
 *  \brief Interfaces for open, closed, overt and compact subsets of hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_SET_INTERFACE_HPP
#define ARIADNE_HYBRID_SET_INTERFACE_HPP

#include <iosfwd>
#include <utility>
#include <map>


#include "../utility/tribool.hpp"
#include "../geometry/set_interface.hpp"
#include "../hybrid/discrete_location.hpp"
#include "../symbolic/space.hpp"
#include "../geometry/set.hpp"
#include "../symbolic/expression_set.hpp"
#include "../hybrid/hybrid_set.decl.hpp"

namespace Ariadne {

//! \ingroup HybridSetSubmodule
//! \ingroup EuclideanExpressionSet
//! \brief Base class for sets in a hybrid space.
class HybridSetInterfaceBase
{
  public:
    virtual ~HybridSetInterfaceBase() = default;
    virtual HybridSetInterfaceBase* clone() const = 0;
    virtual Set<RealVariable> variables(DiscreteLocation) const = 0;
    inline EuclideanSetBase euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
    virtual OutputStream& _write(OutputStream& os) const = 0;
    friend OutputStream& operator<<(OutputStream& os, const HybridSetInterfaceBase& hs);
  protected:
    virtual EuclideanSetInterfaceBase* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \brief Interface for bounded sets in a hybrid space.
class HybridBoundedSetInterface
    : public virtual HybridSetInterfaceBase
{
  public:
    virtual HybridBoundedSetInterface* clone() const = 0;
    virtual Set<DiscreteLocation> locations() const = 0;
    inline EuclideanBoundedSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
    virtual LowerKleenean inside(const HybridExactBoxes& bx) const = 0;
    virtual HybridUpperBoxes bounding_box() const = 0;
  protected:
    virtual EuclideanBoundedSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \brief Interface for overt sets in a hybrid space.
class HybridOvertSetInterface
    : public virtual HybridSetInterfaceBase
{
  public:
    virtual HybridOvertSetInterface* clone() const = 0;
    inline EuclideanOvertSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
    virtual LowerKleenean overlaps(const HybridExactBox& bx) const = 0;
  protected:
    virtual EuclideanOvertSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \brief Interface for open sets in a hybrid space.
class HybridOpenSetInterface
    : public virtual HybridOvertSetInterface
{
  public:
    virtual HybridOpenSetInterface* clone() const = 0;
    inline EuclideanOpenSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
    virtual LowerKleenean covers(const HybridExactBox& bx) const = 0;
  protected:
    virtual EuclideanOpenSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \brief Interface for closed sets in a hybrid space.
class HybridClosedSetInterface
    : public virtual HybridSetInterfaceBase
{
  public:
    virtual HybridClosedSetInterface* clone() const = 0;
    inline EuclideanClosedSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
    virtual LowerKleenean separated(const HybridExactBox& bx) const = 0;
  protected:
    virtual EuclideanClosedSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \brief Interface for compact (closed and bounded) sets in a hybrid space.
class HybridCompactSetInterface
    : public virtual HybridBoundedSetInterface
    , public virtual HybridClosedSetInterface
{
  public:
    virtual HybridCompactSetInterface* clone() const = 0;
    inline EuclideanCompactSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
  protected:
    virtual EuclideanCompactSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \brief Interface for regular (open and closed) sets in a hybrid space.
class HybridRegularSetInterface
    : public virtual HybridOpenSetInterface,
      public virtual HybridClosedSetInterface
{
  public:
    virtual HybridRegularSetInterface* clone() const = 0;
    inline EuclideanRegularSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
  protected:
    virtual EuclideanRegularSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \brief Interface for located (overt and compact) sets in a hybrid space.
class HybridLocatedSetInterface
    : public virtual HybridOvertSetInterface,
      public virtual HybridCompactSetInterface
{
  public:
    virtual HybridLocatedSetInterface* clone() const = 0;
    inline EuclideanLocatedSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
  protected:
    virtual EuclideanLocatedSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \brief Complete set interface for bounded regular sets in a hybrid space.
class HybridRegularLocatedSetInterface
    : public virtual HybridRegularSetInterface,
      public virtual HybridLocatedSetInterface
{
  public:
    virtual HybridRegularLocatedSetInterface* clone() const = 0;
    inline EuclideanRegularLocatedSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
  protected:
    virtual EuclideanRegularLocatedSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

using HybridSetInterface = HybridRegularLocatedSetInterface;

inline OutputStream& operator<<(OutputStream& os, const HybridSetInterfaceBase& s) {
    return s._write(os);
}


//! \ingroup GeometryModule EuclideanSetInterfaceSubModule
//! \brief Interface for bounded sets.
class HybridValidatedBoundedSetInterface
    : public virtual HybridSetInterfaceBase {
  public:
    virtual HybridValidatedBoundedSetInterface* clone() const = 0;
    virtual Set<DiscreteLocation> locations() const = 0;
//    inline ValidatedEuclideanBoundedSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
    virtual ValidatedLowerKleenean inside(const HybridExactBoxes& bx) const = 0;
    virtual HybridUpperBoxes bounding_box() const = 0;
  protected:
    virtual ValidatedEuclideanBoundedSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \brief Interface for overt sets in a hybrid space.
class HybridValidatedOvertSetInterface
    : public virtual HybridSetInterfaceBase
{
  public:
    virtual HybridValidatedOvertSetInterface* clone() const = 0;
//    inline ValidatedEuclideanOvertSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
    virtual ValidatedLowerKleenean overlaps(const HybridExactBox& bx) const = 0;
  protected:
    virtual ValidatedEuclideanOvertSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \ingroup GeometryModule EuclideanSetInterfaceSubModule
//! \brief Interface for open sets.
class HybridValidatedOpenSetInterface
    : public virtual HybridValidatedOvertSetInterface
{
  public:
    virtual HybridValidatedOpenSetInterface* clone() const = 0;
//    inline EuclideanOpenSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
    virtual ValidatedLowerKleenean covers(const ExactBoxType& bx) const = 0;
  protected:
    virtual ValidatedEuclideanOpenSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \ingroup GeometryModule EuclideanSetInterfaceSubModule
//! \brief Interface for closed sets.
class HybridValidatedClosedSetInterface
    : public virtual HybridSetInterfaceBase
{
  public:
    virtual HybridValidatedClosedSetInterface* clone() const = 0;
//    inline ValidatedEuclideanClosedSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
    virtual ValidatedLowerKleenean separated(const HybridExactBox& bx) const = 0;
  protected:
    virtual ValidatedEuclideanClosedSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \ingroup GeometryModule EuclideanSetInterfaceSubModule
//! \brief Interface for compact (closed and bounded) sets.
class HybridValidatedCompactSetInterface
    : public virtual HybridValidatedBoundedSetInterface,
      public virtual HybridValidatedClosedSetInterface
{
  public:
    virtual HybridValidatedCompactSetInterface* clone() const = 0;
//    inline ValidatedEuclideanCompactSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
  protected:
    virtual ValidatedEuclideanCompactSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \ingroup GeometryModule EuclideanSetInterfaceSubModule
//! \brief Interface for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
class HybridValidatedRegularSetInterface
    : public virtual HybridValidatedOpenSetInterface,
      public virtual HybridValidatedClosedSetInterface
{
    virtual HybridValidatedRegularSetInterface* clone() const = 0;
//    inline ValidatedEuclideanRegularSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
  protected:
    virtual ValidatedEuclideanRegularSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \ingroup GeometryModule EuclideanSetInterfaceSubModule
//! \brief Interface for located (overt and compact) sets.
class HybridValidatedLocatedSetInterface
    : public virtual HybridValidatedOvertSetInterface,
      public virtual HybridValidatedCompactSetInterface
{
    virtual HybridValidatedLocatedSetInterface* clone() const = 0;
//    inline ValidatedEuclideanLocatedSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
  protected:
    virtual ValidatedEuclideanLocatedSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};

//! \ingroup GeometryModule EuclideanSetInterfaceSubModule
//! \brief Complete set interface for bounded regular sets.
class HybridValidatedRegularLocatedSetInterface
    : public virtual HybridValidatedRegularSetInterface,
      public virtual HybridValidatedLocatedSetInterface
{
  public:
    virtual HybridValidatedRegularLocatedSetInterface* clone() const = 0;
//    inline ValidatedEuclideanRegularLocatedSet euclidean_set(DiscreteLocation loc, RealSpace spc) const { return this->_euclidean_set(loc,spc); }
  protected:
    virtual ValidatedEuclideanRegularLocatedSetInterface* _euclidean_set(DiscreteLocation,RealSpace) const = 0;
};



} // namespace Ariadne


#endif // ARIADNE_SET_INTERFACE
