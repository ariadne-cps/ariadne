/***************************************************************************
 *            python/export_orthotope.cc
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

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/orthotope.h"
#include "geometry/list_set.h"

#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;


template<class R> inline
Orthotope<R,R> over_approximation_of_minkowski_sum(const Orthotope<R,R>& z1, const Orthotope<R,R>& z2) {
  return over_approximation(minkowski_sum(z1,z2));
}

template<class R> inline
Orthotope<R,R> over_approximation_of_minkowski_difference(const Orthotope<R,R>& z1, const Orthotope<R,R>& z2) {
  return over_approximation(minkowski_difference(z1,z2));
}


template<class R>
void export_orthotope() 
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  
  typedef Interval<R> I;
  typedef Matrix<R> RMatrix;
  typedef Matrix<I> IMatrix;
  typedef Point<R> RPoint;
  typedef Point<F> FPoint;
  typedef Rectangle<R> RRectangle;
  typedef Orthotope<R,R> ROrthotope;
  typedef Orthotope<F,R> EOrthotope;
  typedef Orthotope<F,F> IOrthotope;
    
  class_<EOrthotope>("Orthotope",init<int>())
    .def(init<EOrthotope>())
    .def(init<ROrthotope>())
    .def(init<RRectangle>())
    .def("centre",(const FPoint&(EOrthotope::*)()const)&EOrthotope::centre,return_value_policy<copy_const_reference>())
    .def("generators",(const RMatrix&(EOrthotope::*)()const)&EOrthotope::generators,return_value_policy<copy_const_reference>())
    .def("dimension", &EOrthotope::dimension)
    .def("empty", &EOrthotope::empty)
    .def("contains", (tribool(EOrthotope::*)(const RPoint&)const)&EOrthotope::contains)
    .def("bounding_box", &EOrthotope::bounding_box)
    .def("subdivide", &EOrthotope::subdivide)
    .def("divide", &EOrthotope::divide)
    .def(self_ns::str(self))
  ;
}  

template void export_orthotope<Float>();
