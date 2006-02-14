/***************************************************************************
 *            vector_field.h
 *
 *  Thu Feb  3 21:06:54 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
#ifndef _ARIADNE_VECTOR_FIELD_H
#define _ARIADNE_VECTOR_FIELD_H

#include <base/exception.h>

namespace Ariadne {
namespace VectorField{

enum VectorFieldKind {
  LINEAR,
  AFFINE,
  MULTIVALUE,
  GENERAL
};

}
}

template<typename R> class Ariadne::Geometry::Rectangle;
template<typename R> class Ariadne::Geometry::Polyhedron;
template<typename R> class Ariadne::Geometry::Point;
template<typename R, template<typename> class BS> class Ariadne::Geometry::ListSet;
template<typename R, template<typename> class S> class Ariadne::Map::Map;


namespace Ariadne {
namespace VectorField{

   
template <typename R, template<typename> class S>
class VectorField {
 public:
  typedef R Real;
  typedef S<R> State;
  
  virtual State operator() (const State& x) const = 0;
  virtual Geometry::Rectangle<R> operator() (const Geometry::Rectangle<R>& A) const {
    throw std::invalid_argument("Not implemented."); }
  virtual Geometry::Polyhedron<R> operator() (const Geometry::Polyhedron<R>& A) const {
    throw std::invalid_argument("Not implemented."); };

  size_t dimension() const = 0;
};

    
}
}

#endif /* _ARIADNE_VECTOR_FIELD_H */
