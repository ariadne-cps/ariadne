/***************************************************************************
 *            dynamics/scalings.hpp
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

/*! \file dynamics/scalings.hpp
 *  \brief Scalings for real variables.
 */

#ifndef ARIADNE_SCALINGS_HPP
#define ARIADNE_SCALINGS_HPP

#include <iostream>
#include <map>

#include "../utility/container.hpp"
#include "../utility/stlio.hpp"

#include "../numeric/floatdp.hpp"
#include "../numeric/float_value.hpp"
#include "../geometry/grid.hpp"
#include "../symbolic/variables.hpp"

namespace Ariadne {

typedef void Void;
typedef std::ostream OutputStream;

typedef Value<FloatDP> FloatDPValue;

//! \ingroup DynamicsModule
//! \brief A class which defines the state space grid to use given the continuous state variables \a spc.
class Scalings {
    FloatDPValue _default_scaling;
    Map<Identifier,FloatDPValue> _scalings;
  public:
    Scalings(FloatDPValue default_scaling)
        : _default_scaling(default_scaling), _scalings() { }
    Scalings(FloatDPValue default_scaling, Map<RealVariable,FloatDPValue> const& scalings);
    Void set_scaling(RealVariable const& var, FloatDPValue scal) {
        ARIADNE_ASSERT(decide(scal>0)); _scalings[var.name()]=scal; }
    FloatDPValue scaling(const RealVariable& var) const {
        return (this->_scalings.has_key(var.name())) ? this->_scalings[var.name()] : this->_default_scaling; }
    Grid grid(RealSpace const& spc) const;
    friend OutputStream& operator<<(OutputStream& os, Scalings const& s) {
        return os << "Scalings( " << s._scalings << " )"; }
};

inline Scalings::Scalings(FloatDPValue default_scaling, Map<RealVariable,FloatDPValue> const& scalings)
    : _default_scaling(default_scaling), _scalings()
{
    ARIADNE_ASSERT(default_scaling>0);
    for(auto var_scal : scalings) {
        this->_scalings.insert(var_scal.first.name(),var_scal.second);
    }
}

inline Grid Scalings::grid(RealSpace const& spc) const {
    return Grid(Vector<FloatDP>(spc.dimension(),[this,&spc](SizeType i){return this->scaling(spc[i]).raw();}));
}

} // namespace Ariadne

#endif // ARIADNE_SCALINGS_HPP

