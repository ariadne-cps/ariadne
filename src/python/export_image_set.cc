/***************************************************************************
 *            python/export_image_set.cc
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
#include "geometry/image_set.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_image_set() 
{
  typedef typename traits<R>::arithmetic_type A;

  class_< ImageSet<R>, bases< SetInterface< Box<R> > > >("ImageSet",init< const Box<R>&, const FunctionInterface<R>& >())
    .def(init< const Point<R>& >())
    .def(init< const Box<R>& >())
    .def(init< const ImageSet<R>& >())
    .def("domain", (Box<R>(ImageSet<R>::*)()const)&ImageSet<R>::domain)
    .def("function", (FunctionInterface<R>*(ImageSet<R>::*)()const)&ImageSet<R>::function,
         return_value_policy<manage_new_object>())
    .def("dimension", &ImageSet<R>::dimension)
    .def("contains", (tribool(ImageSet<R>::*)(const Point<R>&)const)&ImageSet<R>::contains)
    .def("superset", &ImageSet<R>::superset)
    .def("intersects", &ImageSet<R>::intersects)
    .def("disjoint", &ImageSet<R>::disjoint)
    .def("subset", &ImageSet<R>::subset)
    .def("bounded", &ImageSet<R>::bounded)
    .def("bounding_box", &ImageSet<R>::bounding_box)
    .def(self_ns::str(self))
  ;
}


template void export_image_set<FloatPy>();
