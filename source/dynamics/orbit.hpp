/***************************************************************************
 *            orbit.hpp
 *
 *  Copyright 2007-17  Pieter Collins
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

/*! \file orbit.hpp
 *  \brief Orbits of dynamic systems
 */

#ifndef ARIADNE_ORBIT_HPP
#define ARIADNE_ORBIT_HPP

#include <utility>
#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "../numeric/numeric.hpp"
#include "../output/graphics_interface.hpp"
#include "../geometry/function_set.hpp"
#include "../geometry/list_set.hpp"
#include "../geometry/enclosure.hpp"


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
class GridTreePaving;


template<>
class Orbit<ExactPoint>
{
  public:
    Orbit(const ExactPoint& pt);
    Void insert(FloatDPValue t, const ExactPoint& hpt);
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
    typedef GridTreePaving EnclosureListType;

    Orbit(const Grid&, const GridCell&);
    Orbit(const GridTreePaving&);
    Orbit(const GridTreePaving&, const GridTreePaving&,
          const GridTreePaving&, const GridTreePaving&);
    Grid const& grid() const;
    GridTreePaving const& initial() const;
    GridTreePaving const& reach() const;
    GridTreePaving const& intermediate() const;
    GridTreePaving const& final() const;
  private:
    std::shared_ptr<Data> _data;
};


template<>
class Orbit<Enclosure>
{
    typedef Enclosure ES;
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

#endif // ARIADNE_ORBIT_HPP
