/***************************************************************************
 *            python/geometry.cc
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

#include <iostream>
#include "numerical_type.h"
#include "rectangle.h"
#include "list_set.h"

#include <boost/python.hpp>

#include "real_typedef.h"
#include "container_utilities.h"

typedef Ariadne::Geometry::State<Real> RState;
typedef Ariadne::Interval<Real> RInterval;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Rectangle> RRectangleListSet;

using Ariadne::Geometry::regular_intersection;
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
  typedef RRectangle (*RectBinFun) (const RRectangle&, const RRectangle&);
  RectBinFun rect_regular_intersection=&regular_intersection<Real>;
  RectBinPred rect_interiors_intersect=&interiors_intersect<Real>;
  RectBinPred rect_disjoint=&disjoint<Real>;
  RectBinPred rect_inner_subset=&inner_subset<Real>;
  RectBinPred rect_subset=&subset<Real>;

  def("regular_intersection", rect_regular_intersection);
  def("interiors_intersect", rect_interiors_intersect);
  def("disjoint", rect_disjoint);
  def("inner_subset", rect_inner_subset);
  def("subset", rect_subset);

  class_<RRectangle>("Rectangle",init<int>())
    .def(init<RState,RState>())
    .def(init<RRectangle>())
    .def("dimension", &RRectangle::dimension)
    .def("__getitem__", &RRectangle::get)
    .def("__setitem__", &RRectangle::set)
    .def("set_lower", &RRectangle::set_lower)
    .def("set_upper", &RRectangle::set_upper)
    .def("lower", &RRectangle::lower)
    .def("upper", &RRectangle::upper)
    .def(str(self))    // __str__
  ;
}
