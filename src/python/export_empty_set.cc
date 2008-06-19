/***************************************************************************
 *            python/export_empty_set.cc
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

#include "geometry/empty_set.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_empty_set() 
{
  typedef Box<R> BS;
  class_< EmptySet<R>, bases< SetInterface<BS> > >("EmptySet",init< dimension_type >())
    .def("space", &EmptySet<R>::space)
    .def("contains", &EmptySet<R>::contains)
    .def("superset", &EmptySet<R>::superset)
    .def("intersects", &EmptySet<R>::intersects)
    .def("disjoint", &EmptySet<R>::disjoint)
    .def("subset", &EmptySet<R>::subset)
    .def("bounding_box", &EmptySet<R>::bounding_box)
    .def(self_ns::str(self))
  ;
}


template void export_empty_set<FloatPy>();
