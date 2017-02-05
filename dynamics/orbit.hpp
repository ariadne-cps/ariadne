/***************************************************************************
 *            orbit.h
 *
 *  Copyright 2007-17  Pieter Collins
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

/*! \file orbit.h
 *  \brief Orbits of dynamic systems
 */

#ifndef ARIADNE_ORBIT_H
#define ARIADNE_ORBIT_H

#include <utility>
#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <boost/scoped_ptr.hpp>

#include "numeric/numeric.h"
#include "output/graphics_interface.h"
#include "geometry/function_set.h"



namespace Ariadne {

#ifdef DOXYGEN
//! \brief Class for storing evolution data.
template<class E> class Orbit {
  public:
    //! The type used to store a single enclosure set.
    typedef E EnclosureType;
    //! The type used to store a list of enclosure sets.
    typedef EL EnclosureListType;
    //! The initial set of the orbit.
    EnclosureType const& initial() const;
    //! The set of points reached by evolution from the given initial set over the evolution time.
    EnclosureListType const& reach() const;
    //! The set of points reached by evolution from the given initial set at the final evolution time.
    EnclosureListType const& final() const;
};
#endif

typedef Real TimeType;

template<class ES> class Orbit;
template<class ES> OutputStream& operator<<(OutputStream&, const Orbit<ES>&);

template<class BS> class ListSet;

template<class X> class Point;
typedef Point<ExactNumericType> ExactPoint;

class InterpolatedCurve;
class Grid;
class GridCell;
class GridTreeSet;


template<>
class Orbit<ExactPoint>
{
  public:
    Orbit(const ExactPoint& pt);
    Void insert(Float64Value t, const ExactPoint& hpt);
    const InterpolatedCurve& curve() const { return *this->_curve; }
  private:
    std::shared_ptr< InterpolatedCurve > _curve;
};

template<>
class Orbit<GridCell>
{
    struct Data;
  public:
    typedef GridCell EnclosureType;
    typedef GridTreeSet EnclosureListType;

    Orbit(const Grid&, const GridCell&);
    Orbit(const GridTreeSet&);
    Orbit(const GridTreeSet&, const GridTreeSet&,
          const GridTreeSet&, const GridTreeSet&);
    Grid const& grid() const;
    GridTreeSet const& initial() const;
    GridTreeSet const& reach() const;
    GridTreeSet const& intermediate() const;
    GridTreeSet const& final() const;
  private:
    std::shared_ptr<Data> _data;
};


template<class ES>
class Orbit
{
    typedef ListSet<ES> ESL;
  public:
    typedef ES EnclosureType;
    typedef ListSet<ES> EnclosureListType;

    Orbit(const ES& set) : _initial(set) { }
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
  private:
    ES _initial;
    ESL _reach;
    ESL _intermediate;
    ESL _final;
};

template<class ES> OutputStream& operator<<(OutputStream& os, const Orbit< ES >& orb);

template<class ES>
OutputStream&
operator<<(OutputStream& os, const Orbit< ES >& orb)
{
    os << "Orbit(\n  initial=" << orb.initial()
       << "\n  intermediate=" << orb.intermediate()
       << "\n  reach=" << orb.reach()
       << "\n  final=" << orb.final()
       << ")\n";
    return os;
}

template<class ES> Void draw(FigureInterface& figure, const Orbit<ES>& orbit) {
    draw(figure,orbit.reach());
    draw(figure,orbit.initial());
    draw(figure,orbit.final());
}

template<class ES> FigureInterface& operator<<(FigureInterface& figure, const Orbit<ES>& orbit) {
    draw(figure,orbit); return figure;
}

} // namespace Ariadne

#endif // ARIADNE_ORBIT_H
