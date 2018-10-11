/***************************************************************************
 *            map.hpp
 *
 *  Copyright  2004-8  Alberto Casagrande, Pieter Collins
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

/*! \file map.hpp
 *  \brief Main continuous dynamics system class.
 */

#ifndef ARIADNE_MAP_HPP
#define ARIADNE_MAP_HPP

#include <memory>

#include "../function/function.hpp"
#include "../geometry/set_interface.hpp"
#include "../geometry/grid.hpp"

namespace Ariadne {

/*! \brief An iterated function system in Euclidean space.
 */
class IteratedMap
{
  public:
    //! \brief The type used to represent time.
    typedef Integer TimeType;
    //! \brief The type used to represent real numbers.
    typedef Real RealType ;
    //! \brief The type used to describe the state space.
    typedef EuclideanSpace StateSpaceType;
  public:
    IteratedMap(const EffectiveVectorMultivariateFunction& f) : _function(f) { }
    virtual IteratedMap* clone() const { return new IteratedMap(*this); }
    virtual ~IteratedMap() = default;
    const EffectiveVectorMultivariateFunction& function() const { return _function; }
    Grid grid() const { return Grid(_function.argument_size()); }
  private:
    EffectiveVectorMultivariateFunction _function;
};

inline OutputStream& operator<<(OutputStream& os, const IteratedMap& vf) {
    return os << "IteratedMap( " << vf.function() << " )";
}


} // namespace Ariadne

#endif // ARIADNE_MAP_HPP
