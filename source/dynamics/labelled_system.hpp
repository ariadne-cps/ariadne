/***************************************************************************
 *            dynamics/labelled_system.hpp
 *
 *  Copyright  2020  Pieter Collins
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

/*! \file dynamics/labelled_system.hpp
 *  \brief Labelling systems with named variables.
 */

#ifndef ARIADNE_LABELLED_SYSTEM_HPP
#define ARIADNE_LABELLED_SYSTEM_HPP

#include "../function/function.hpp"
#include "../symbolic/space.hpp"

namespace Ariadne {

class RealParameter;

/*! \brief A base class for systems with named state and auxiliary variables. */
class ExtendedSystemMixin {
    EffectiveVectorMultivariateFunction _auxiliary_map;
    List<Identifier> _state_variables;
    List<Identifier> _auxiliary_variables;
  protected:
    ~ExtendedSystemMixin() = default;
    ExtendedSystemMixin(SizeType n);
    ExtendedSystemMixin(List<RealAssignment> const& auxiliary, RealSpace const& state_space);
  public:
    const EffectiveVectorMultivariateFunction& auxiliary_map() const { return this->_auxiliary_map; }

    RealSpace state_space() const;
    RealSpace auxiliary_space() const;
    RealSpace state_auxiliary_space() const;
  public:
    OutputStream& write(OutputStream& os) const {
        return os << "state_variables=" << this->_state_variables<<", "
                  << "auxiliary_map=" << this->_auxiliary_map << ", "
                  << "auxiliary_variables=" << this->_auxiliary_variables;

    }
};

} // namespace Ariadne

#endif // ARIADNE_LABELLED_SYSTEM_HPP
