/***************************************************************************
 *            python/export_level_set.cc
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

#include "python/float.h"

#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/level_set.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_level_set() 
{
  typedef typename Numeric::traits<R>::arithmetic_type A;

  class_< LevelSet<R>, bases< SetInterface<R> > >("LevelSet",init< const FunctionInterface<R>& >())
    .def("function", (Point<A>(LevelSet<R>::*)(const Point<A>&)const) &LevelSet<R>::function)
    .def("separates", (tribool(LevelSet<R>::*)(const Point<A>&,const Point<A>&)const) &LevelSet<R>::separates)
    .def("dimension", &LevelSet<R>::dimension)
    .def("contains", &LevelSet<R>::contains)
    .def("superset", &LevelSet<R>::superset)
    .def("intersects", &LevelSet<R>::intersects)
    .def("disjoint", &LevelSet<R>::disjoint)
    .def("subset", &LevelSet<R>::subset)
    .def("bounded", &LevelSet<R>::bounded)
    .def("bounding_box", &LevelSet<R>::bounding_box)
    .def(self_ns::str(self))
  ;
}


template void export_level_set<FloatPy>();
