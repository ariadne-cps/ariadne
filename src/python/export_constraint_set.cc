/***************************************************************************
 *            python/export_constraint_set.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "python/python_float.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/constraint_set.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_constraint_set() 
{
  typedef typename Numeric::traits<R>::arithmetic_type A;

  class_< ConstraintSet<R>, bases< SetInterface<R> > >("ConstraintSet",init< const FunctionInterface<R>& >())
    .def("function", (Point<A>(ConstraintSet<R>::*)(const Point<A>&)const)&ConstraintSet<R>::function)
    .def("contains", (tribool(ConstraintSet<R>::*)(const Point<A>&)const)&ConstraintSet<R>::contains)
    .def("dimension", &ConstraintSet<R>::dimension)
    .def("contains", (tribool(ConstraintSet<R>::*)(const Point<R>&)const)&ConstraintSet<R>::contains)
    .def("superset", &ConstraintSet<R>::superset)
    .def("intersects", &ConstraintSet<R>::intersects)
    .def("disjoint", &ConstraintSet<R>::disjoint)
    .def("subset", &ConstraintSet<R>::subset)
    .def("bounded", &ConstraintSet<R>::bounded)
    .def("bounding_box", &ConstraintSet<R>::bounding_box)
    .def(self_ns::str(self))
  ;
}


template void export_constraint_set<Float>();
