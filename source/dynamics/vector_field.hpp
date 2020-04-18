/***************************************************************************
 *            dynamics/vector_field.hpp
 *
 *  Copyright  2004-20  Alberto Casagrande, Pieter Collins
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

/*! \file dynamics/vector_field.hpp
 *  \brief Main continuous dynamics system class.
 */

#ifndef ARIADNE_VECTOR_FIELD_HPP
#define ARIADNE_VECTOR_FIELD_HPP

#include <memory>

#include "../function/function.hpp"
#include "../geometry/set_interface.hpp"
#include "../geometry/grid.hpp"
#include "../symbolic/expression.decl.hpp"

namespace Ariadne {

class Enclosure;
class Storage;
class VectorFieldEvolver;

//! \brief A vector field in Euclidean space.
class VectorField
{
  public:
    //! \brief The type used to represent time.
    typedef Real TimeType;
    //! \brief The type used to represent real numbers.
    typedef Real RealType;
    //! \brief The type used to describe the state space.
    typedef EuclideanSpace StateSpaceType;
    //! \brief The generic type used to compute the system evolution.
    typedef VectorFieldEvolver EvolverType;
    typedef Enclosure EnclosureType;
    //! \brief The type used to define global pavings of reach and evolve sets.
    typedef Storage StorageType;
  public:
    VectorField(List<DottedRealAssignment> const& dynamics);
    VectorField(EffectiveVectorMultivariateFunction const& function);
    virtual ~VectorField() = default;
    virtual VectorField* clone() const { return new VectorField(*this); }
    SizeType dimension() const { return _function.result_size(); }
    RealSpace state_space() const;
    const EffectiveVectorMultivariateFunction& function() const { return _function; }
    Grid grid() const { return Grid(_function.result_size()); }
    friend OutputStream& operator<<(OutputStream& os, const VectorField& vf) {
        return os << "VectorField( " << vf.function() << " )"; }
  private:
    List<Identifier> _variable_names;
    EffectiveVectorMultivariateFunction _function;
};

} // namespace Ariadne

#endif // ARIADNE_VECTOR_FIELD_HPP
