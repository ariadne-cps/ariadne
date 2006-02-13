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
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "henon_map.h"

#include <boost/python.hpp>

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::bases;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

#include "real_typedef.h"
  
typedef Ariadne::Interval<Real> RInterval;
typedef Ariadne::Geometry::Point<Real> RPoint;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::Parallelopiped<Real> RParallelopiped;
typedef Ariadne::LinearAlgebra::matrix<Real> RMatrix;
typedef Ariadne::LinearAlgebra::matrix<RInterval> IMatrix;
typedef Ariadne::Evaluation::Map<Real> RMap;
typedef Ariadne::Evaluation::HenonMap<Real> RHenonMap;

typedef RPoint (RHenonMap::* PointMap) (const RPoint&) const;
typedef RRectangle (RHenonMap::* RectangleMap) (const RRectangle&) const;
typedef RParallelopiped (RHenonMap::* ParallelopipedMap) (const RParallelopiped&) const;
typedef RMatrix (RHenonMap::* PointDerivative) (const RPoint&) const;
typedef IMatrix (RHenonMap::* RectangleDerivative) (const RRectangle&) const;

void export_henon_map() {
  class_<RHenonMap, bases<RMap> >("HenonMap",init<Real,Real>())
    .def("argument_dimension", &RHenonMap::argument_dimension)
    .def("result_dimension", &RHenonMap::result_dimension)
    .def("__call__", PointMap(&RHenonMap::apply))
    .def("__call__", RectangleMap(&RHenonMap::apply))
    .def("__call__", ParallelopipedMap(&RHenonMap::apply))
    .def("apply", PointMap(&RHenonMap::apply))
    .def("apply", RectangleMap(&RHenonMap::apply))
    .def("apply", ParallelopipedMap(&RHenonMap::apply))
    .def("derivative", PointDerivative(&RHenonMap::derivative))
    .def("derivative", RectangleDerivative(&RHenonMap::derivative))
    .def(str(self))    // __str__
  ;
}
