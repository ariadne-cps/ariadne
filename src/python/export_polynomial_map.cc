/***************************************************************************
 *            python/export_rectangle.cc
 *
 *  21 October 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include <base/numerical_type.h>

#include <geometry/point.h>

#include <evaluation/polynomial_map.h>

#include <boost/python.hpp>

#include "python/real_typedef.h"

typedef Ariadne::Geometry::Point<Real> RPoint;
typedef Ariadne::Evaluation::PolynomialMap<Real> RPolynomialMap;

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

void export_polynomial_map() {

  class_<RPolynomialMap>("PolynomialMap",init<int,int>())
    .def("argument_dimension", &RPolynomialMap::argument_dimension)
    .def("result_dimension", &RPolynomialMap::result_dimension)
    .def(str(self))    // __str__
  ;
}
