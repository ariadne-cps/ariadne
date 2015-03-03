/***************************************************************************
 *            hybrid_orbit.h
 *
 *  Copyright 2007-10  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file hybrid_orbit.h
 *  \brief Orbits of hybrid dynamical systems
 */

#ifndef ARIADNE_HYBRID_ORBIT_H
#define ARIADNE_HYBRID_ORBIT_H

#include <utility>
#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <boost/scoped_ptr.hpp>

#include "numeric/numeric.h"
#include "output/graphics_interface.h"
#include "geometry/function_set.h"
#include "hybrid/hybrid_enclosure.h"

#include "dynamics/orbit.h"

namespace Ariadne {

typedef Float64 Time;

template<class ES> class Orbit;

template<class BS> class ListSet;
template<class X> class Point;
typedef Point<ExactNumericType> ExactPoint;
class Grid;
class GridCell;
class GridTreeSet;
class HybridGrid;
class HybridGridCell;
class HybridGridTreeSet;
class HybridTime;

class DiscreteLocation;
template<class BS> class HybridBasicSet;

class HybridPoint;
class HybridBoxType;
typedef HybridBasicSet<InterpolatedCurve> HybridInterpolatedCurve;

template<class ES> OutputStream& operator<<(OutputStream&, const Orbit<ES>&);

template<>
class Orbit<HybridPoint>
    : public HybridDrawableInterface
{
  public:
    Orbit(const HybridPoint& hpt);
    Void insert(HybridTime ht, const HybridPoint& hpt);
    Nat size() const;
    const InterpolatedCurve& curve(Nat m) const;
    const std::vector<HybridInterpolatedCurve>& curves() const { return *this->_curves_ptr; }
    Void draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d& v) const;
  private:
    std::shared_ptr<std::vector<HybridInterpolatedCurve> > _curves_ptr;
};

template<>
OutputStream&
operator<<(OutputStream& os, const Orbit< HybridPoint >& orb);

template<>
class Orbit<HybridGridCell>
{
    class Data;
  public:
    typedef HybridGridCell EnclosureType;
    typedef HybridGridTreeSet EnclosureListType;

    Orbit(const HybridGrid&, const HybridGridCell&);
    Orbit(const HybridGridTreeSet&);
    Orbit(const HybridGridTreeSet&, const HybridGridTreeSet&,
          const HybridGridTreeSet&, const HybridGridTreeSet&);
    HybridGrid const& grid() const;
    HybridGridTreeSet const& initial() const;
    HybridGridTreeSet const& reach() const;
    HybridGridTreeSet const& intermediate() const;
    HybridGridTreeSet const& final() const;
  private:
    std::shared_ptr<Data> _data;
};

template<>
class Orbit<HybridEnclosure>
    : public HybridDrawableInterface
{
  public:
    typedef HybridEnclosure EnclosureType;
    typedef ListSet<HybridEnclosure> EnclosureListType;

    Orbit() { }
    Orbit(const EnclosureType& set) : _initial(set) { }
    Void adjoin_reach(const EnclosureType& set) { this->_reach.adjoin(set); }
    Void adjoin_intermediate(const EnclosureType& set) { this->_intermediate.adjoin(set); }
    Void adjoin_final(const EnclosureType& set) { this->_final.adjoin(set); }

    Void adjoin_reach(const EnclosureListType& set) { this->_reach.adjoin(set); }
    Void adjoin_intermediate(const EnclosureListType& set) { this->_intermediate.adjoin(set); }
    Void adjoin_final(const EnclosureListType& set) { this->_final.adjoin(set); }

    EnclosureType const& initial() const { return this->_initial; }
    EnclosureListType const& reach() const { return this->_reach; }
    EnclosureListType const& intermediate() const { return this->_intermediate; }
    EnclosureListType const& final() const { return this->_final; }

    Void draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d& v) const;
  private:
    EnclosureType _initial;
    EnclosureListType _reach;
    EnclosureListType _intermediate;
    EnclosureListType _final;
};

template<>
OutputStream&
operator<<(OutputStream& os, const Orbit< HybridEnclosure >& orb);


} // namespace Ariadne

#endif // ARIADNE_HYBRID_ORBIT_H
