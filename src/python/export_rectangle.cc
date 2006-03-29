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

#include "base/numerical_type.h"

#include "geometry/rectangle.h"
#include "geometry/list_set.h"

#include <boost/python.hpp>

#include "python/real_typedef.h"

typedef Ariadne::Geometry::Point<Real> RPoint;
typedef Ariadne::Interval<Real> RInterval;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Rectangle> RRectangleListSet;

using Ariadne::Geometry::regular_intersection;
using Ariadne::Geometry::intersection;

using Ariadne::Geometry::interiors_intersect;
using Ariadne::Geometry::disjoint;
using Ariadne::Geometry::inner_subset;
using Ariadne::Geometry::subset;

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

void export_rectangle() {
  typedef bool (*RectBinPred) (const RRectangle&, const RRectangle&);
  typedef RRectangle (*RectBinFunc) (const RRectangle&, const RRectangle&);
  
  def("regular_intersection", RectBinFunc(&regular_intersection));
  def("interiors_intersect", RectBinPred(&interiors_intersect));
  def("disjoint", RectBinPred(&disjoint));
  def("inner_subset", RectBinPred(&inner_subset));

  def("intersection", RectBinFunc(&intersection));
  def("subset", RectBinPred(&subset));
  
  class_<RRectangle>("Rectangle",init<int>())
    .def(init<RPoint,RPoint>())
    .def(init<RRectangle>())
    .def(init<std::string>())
    .def("empty", &RRectangle::empty)
    .def("empty_interior", &RRectangle::empty_interior)
    .def("dimension", &RRectangle::dimension)
    .def("contains", &RRectangle::contains)
    .def("interior_contains", &RRectangle::interior_contains)
    .def("__getitem__", &RRectangle::interval)
    .def("__setitem__", &RRectangle::set_interval)
    .def("set_lower_bound", &RRectangle::set_lower_bound)
    .def("set_upper_bound", &RRectangle::set_upper_bound)
    .def("bounding_box", &RRectangle::bounding_box)
    .def("lower_corner", &RRectangle::lower_corner)
    .def("upper_corner", &RRectangle::upper_corner)
    .def("lower_bound", &RRectangle::lower_bound, return_value_policy<copy_const_reference>())
    .def("upper_bound", &RRectangle::upper_bound, return_value_policy<copy_const_reference>())
    .def(str(self))    // __str__
  ;
}
