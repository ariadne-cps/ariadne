/***************************************************************************
 *            vector_field.hpp
 *
 *  Copyright  2004-8  Alberto Casagrande, Pieter Collins
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

/*! \file vector_field.hpp
 *  \brief Main continuous dynamics system class.
 */

#ifndef ARIADNE_VECTOR_FIELD_HPP
#define ARIADNE_VECTOR_FIELD_HPP

#include <memory>

#include "function/function.hpp"
#include "geometry/set_interface.hpp"
#include "geometry/grid.hpp"
#include "expression/expression.decl.hpp"

namespace Ariadne {

//! \brief A vector field in Euclidean space.
class VectorField
{
  public:
    //! \brief The type used to represent time.
    typedef Real TimeType;
    //! \brief The type used to represent real numbers.
    typedef Real RealType ;
    //! \brief The type used to describe the state space.
    typedef EuclideanSpace StateSpaceType;
  public:
    VectorField(List<DottedRealAssignment> const& dynamics);
    VectorField(EffectiveVectorFunction const& function) : _function(function) { }
    virtual ~VectorField() { }
    virtual VectorField* clone() const { return new VectorField(*this); }
    SizeType dimension() const { return _function.result_size(); }
    const EffectiveVectorFunction& function() const { return _function; }
    Grid grid() const { return Grid(_function.argument_size()); }
    friend OutputStream& operator<<(OutputStream& os, const VectorField& vf) {
        return os << "VectorField( " << vf.function() << " )"; }
  private:
    List<std::string> _variable_names;
    EffectiveVectorFunction _function;
};

} // namespace Ariadne

#endif // ARIADNE_VECTOR_FIELD_HPP
