/***************************************************************************
 *            python/export_polyhedral_set.cc
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

#include "geometry/polyhedral_set.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_polyhedral_set() 
{
  class_< PolyhedralSet<R>, bases< SetInterface<R> > >("PolyhedralSet",init< Polyhedron<R> >())
    .def(init< Matrix<R>,Vector<R> >())
    .def(init< Box<R> >())
    .def(init< Polyhedron<R> >())
    .def("polyhedron", &PolyhedralSet<R>::operator const Polyhedron<R>&,return_value_policy<copy_const_reference>())
    .def("dimension", &PolyhedralSet<R>::dimension)
    .def("contains", &PolyhedralSet<R>::contains)
    .def("superset", &PolyhedralSet<R>::superset)
    .def("intersects", &PolyhedralSet<R>::intersects)
    .def("disjoint", &PolyhedralSet<R>::disjoint)
    .def("subset", &PolyhedralSet<R>::subset)
    .def("bounded", &PolyhedralSet<R>::bounded)
    .def("bounding_box", &PolyhedralSet<R>::bounding_box)
    .def(self_ns::str(self))
  ;
}


template void export_polyhedral_set<FloatPy>();
