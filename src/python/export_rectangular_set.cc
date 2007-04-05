/***************************************************************************
 *            python/export_rectangular_set.cc
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

#include "geometry/rectangular_set.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_rectangular_set() 
{
  class_< RectangularSet<R>, bases< SetInterface<R>, Rectangle<R> > >("RectangularSet",init< Rectangle<R> >())
    .def(init< std::string >())
    .def(init< Rectangle<R> >())
    .def("dimension", &RectangularSet<R>::dimension)
    .def("contains", &RectangularSet<R>::contains)
    .def("superset", &RectangularSet<R>::superset)
    .def("intersects", &RectangularSet<R>::intersects)
    .def("disjoint", &RectangularSet<R>::disjoint)
    .def("subset", &RectangularSet<R>::subset)
    .def("bounding_box", &RectangularSet<R>::bounding_box)
    .def(self_ns::str(self))
  ;
}


template void export_rectangular_set<Float>();
