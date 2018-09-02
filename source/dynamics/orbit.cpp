/***************************************************************************
 *            orbit.cpp
 *
 *  Copyright 2007--17 Pieter Collins
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

#include "../function/functional.hpp"

#include "../config.hpp"

#include <utility>

#include "../dynamics/orbit.hpp"

#include "../geometry/box.hpp"
#include "../geometry/point.hpp"
#include "../geometry/curve.hpp"
#include "../geometry/function_set.hpp"
#include "../geometry/list_set.hpp"
#include "../geometry/grid_paving.hpp"

namespace Ariadne {


Orbit<ExactPoint>::Orbit(const ExactPoint& pt)
    : _curve(new InterpolatedCurve(0,pt))
{ }

Void
Orbit<ExactPoint>::insert(FloatDPValue t, const ExactPoint& pt)
{
    this->_curve->insert(t,pt);
}




struct Orbit<GridCell>::Data {
    Data(const Grid& grid) : initial(grid), reach(grid), intermediate(grid), final(grid) { }
    GridTreePaving initial;
    GridTreePaving reach;
    GridTreePaving intermediate;
    GridTreePaving final;
};

Orbit<GridCell>::
Orbit(const GridTreePaving& initial_set)
    : _data(new Data(initial_set.grid()))
{
    this->_data->initial=initial_set;
}


Orbit<GridCell>::
Orbit(const GridTreePaving& initial_set,
      const GridTreePaving& reach_set,
      const GridTreePaving& intermediate_set,
      const GridTreePaving& final_set)
    : _data(new Data(initial_set.grid()))
{
    this->_data->initial=initial_set;
    this->_data->reach=reach_set;
    this->_data->intermediate=intermediate_set;
    this->_data->final=final_set;
}


GridTreePaving const&
Orbit<GridCell>::
initial() const
{
    return this->_data->initial;
}

GridTreePaving const&
Orbit<GridCell>::
reach() const
{
    return this->_data->reach;
}

GridTreePaving const&
Orbit<GridCell>::
intermediate() const
{
    return this->_data->intermediate;
}

GridTreePaving const&
Orbit<GridCell>::
final() const
{
    return this->_data->final;
}





} // namespace Ariadne
