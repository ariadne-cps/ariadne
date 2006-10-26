/***************************************************************************
 *            python/export_henon_map.cc
 *
 *  8 February 2006
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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

#include "real_typedef.h"

#include "system/henon_map.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;

#include <boost/python.hpp>
using namespace boost::python;
  
template<class R>
void export_henon_map() {
  typedef typename Numeric::traits<R>::arithmetic_type F;
  
  class_< HenonMap<R>, bases< Map<R> > >("HenonMap",init<R,R>())
    .def("argument_dimension", &HenonMap<R>::argument_dimension)
    .def("result_dimension", &HenonMap<R>::result_dimension)
    .def("__call__", (Point<F>(HenonMap<R>::*)(const Point<R>&)const)(&HenonMap<R>::image))
    .def("__call__", (Rectangle<R>(HenonMap<R>::*)(const Rectangle<R>&)const)(&HenonMap<R>::image))
    .def("jacobian", (Matrix<F>(HenonMap<R>::*)(const Point<R>&)const)(&HenonMap<R>::jacobian))
    .def("jacobian", (Matrix< Interval<R> >(HenonMap<R>::*)(const Rectangle<R>&)const)(&HenonMap<R>::jacobian))
    .def(self_ns::str(self))
  ;
}

template void export_henon_map<Real>();
