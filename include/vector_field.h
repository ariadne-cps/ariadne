/***************************************************************************
 *            vector_field.h
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

/*! \file vector_field.h
 *  \brief Main continuous dynamics system class.
 */

#ifndef ARIADNE_VECTOR_FIELD_H
#define ARIADNE_VECTOR_FIELD_H

#include <memory>

#include "function.h"
#include "set_interface.h"
#include "grid.h"

namespace Ariadne {

/*! \brief A vector field in Euclidean space.
 */
class VectorField
{
  public:
    //! \brief The type used to represent time.
    typedef Float TimeType;
    //! \brief The type used to represent real numbers.
    typedef Float RealType ;
    //! \brief The type used to describe the state space.
    typedef EuclideanSpace StateSpaceType;
  public:
    VectorField(const RealVectorFunction& f) : _function(f) { }
    virtual ~VectorField() { }
    virtual VectorField* clone() const { return new VectorField(*this); }
    const RealVectorFunction& function() const { return _function; }
    Grid grid() const { return Grid(_function.argument_size()); }
  private:
    RealVectorFunction _function;
};

inline std::ostream& operator<<(std::ostream& os, const VectorField& vf) {
    return os << "VectorField( " << vf.function() << " )";
}


} // namespace Ariadne

#endif // ARIADNE_VECTOR_FIELD_H
