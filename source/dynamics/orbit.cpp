/***************************************************************************
 *            orbit.cpp
 *
 *  Copyright 2007--17 Pieter Collins
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

#include "../function/functional.hpp"

#include "../config.hpp"

#include <utility>

#include "../dynamics/orbit.hpp"

#include "../geometry/box.hpp"
#include "../geometry/point.hpp"
#include "../geometry/curve.hpp"
#include "../geometry/function_set.hpp"
#include "../geometry/list_set.hpp"
#include "../geometry/grid_set.hpp"

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
    GridTreeSet initial;
    GridTreeSet reach;
    GridTreeSet intermediate;
    GridTreeSet final;
};

Orbit<GridCell>::
Orbit(const GridTreeSet& initial_set)
    : _data(new Data(initial_set.grid()))
{
    this->_data->initial=initial_set;
}


Orbit<GridCell>::
Orbit(const GridTreeSet& initial_set,
      const GridTreeSet& reach_set,
      const GridTreeSet& intermediate_set,
      const GridTreeSet& final_set)
    : _data(new Data(initial_set.grid()))
{
    this->_data->initial=initial_set;
    this->_data->reach=reach_set;
    this->_data->intermediate=intermediate_set;
    this->_data->final=final_set;
}


GridTreeSet const&
Orbit<GridCell>::
initial() const
{
    return this->_data->initial;
}

GridTreeSet const&
Orbit<GridCell>::
reach() const
{
    return this->_data->reach;
}

GridTreeSet const&
Orbit<GridCell>::
intermediate() const
{
    return this->_data->intermediate;
}

GridTreeSet const&
Orbit<GridCell>::
final() const
{
    return this->_data->final;
}





} // namespace Ariadne
