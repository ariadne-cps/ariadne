/***************************************************************************
 *            dynamics/orbit.cpp
 *
 *  Copyright  2007-20  Pieter Collins
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

#include "function/functional.hpp"

#include "config.hpp"

#include <utility>

#include "dynamics/orbit.hpp"

#include "geometry/box.hpp"
#include "geometry/point.hpp"
#include "geometry/curve.hpp"
#include "geometry/function_set.hpp"
#include "geometry/list_set.hpp"
#include "io/graphics.hpp"

namespace Ariadne {


template<class F> Orbit<ExactPoint<F>>::Orbit(const ExactPoint<F>& pt)
    : _curve(new InterpolatedCurve(0,pt))
{ }

template<class F> Void
Orbit<ExactPoint<F>>::insert(Value<F> t, const ExactPoint<F>& pt)
{
    this->_curve->insert(t,pt);
}


struct Orbit<Storage>::Data {
    Data(const Grid& grid, const EffectiveVectorMultivariateFunction& auxiliary_mapping)
        : initial(grid,auxiliary_mapping), reach(initial), intermediate(initial), final(initial) { }
    Storage initial;
    Storage reach;
    Storage intermediate;
    Storage final;
};

Orbit<Storage>::
Orbit(const Storage& initial_set)
    : _data(new Data(initial_set.grid(),initial_set.auxiliary_mapping()))
{
    this->_data->initial=initial_set;
}


Orbit<Storage>::
Orbit(const Storage& initial_set,
      const Storage& reach_set,
      const Storage& intermediate_set,
      const Storage& final_set)
    : _data(new Data(initial_set.grid(),initial_set.auxiliary_mapping()))
{
    this->_data->initial=initial_set;
    this->_data->reach=reach_set;
    this->_data->intermediate=intermediate_set;
    this->_data->final=final_set;
}


Storage const&
Orbit<Storage>::
initial() const
{
    return this->_data->initial;
}

Storage const&
Orbit<Storage>::
reach() const
{
    return this->_data->reach;
}

Storage const&
Orbit<Storage>::
intermediate() const
{
    return this->_data->intermediate;
}

Storage const&
Orbit<Storage>::
final() const
{
    return this->_data->final;
}





} // namespace Ariadne
