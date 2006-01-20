/***************************************************************************
 *            python/export_list_set.cc
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

#include "numerical_type.h"
#include "rectangle.h"
#include "list_set.h"

#include <boost/python.hpp>

#include "real_typedef.h"
#include "python_utilities.h"

typedef Ariadne::Geometry::Point<Real> RPoint;
typedef Ariadne::Interval<Real> RInterval;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Rectangle> RRectangleListSet;

using Ariadne::Geometry::Rectangle;
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

void export_list_set() {
  typedef bool (*RectLSBinPred) (const RRectangleListSet&, const RRectangleListSet&);
  typedef RRectangleListSet (*RectLSBinFun) (const RRectangleListSet&, const RRectangleListSet&);
  RectLSBinFun rectls_regular_intersection=&regular_intersection<Real,Rectangle>;
  RectLSBinPred rectls_interiors_intersect=&interiors_intersect<Real,Rectangle>;
  RectLSBinPred rectls_disjoint=&disjoint<Real,Rectangle>;
  RectLSBinPred rectls_inner_subset=&inner_subset<Real,Rectangle>;
  RectLSBinPred rectls_subset=&subset<Real,Rectangle>;

  def("regular_intersection", rectls_regular_intersection);
  def("interiors_intersect", rectls_interiors_intersect);
  def("disjoint", rectls_disjoint);
  def("inner_subset", rectls_inner_subset);
  def("subset", rectls_subset);

  class_<RRectangleListSet>("RectangleListSet",init<int>())
    .def(init<RRectangle>())
    .def(init<RRectangleListSet>())
    .def("dimension", &RRectangleListSet::dimension)
    .def("push_back", &RRectangleListSet::push_back)
    .def("__len__", &RRectangleListSet::size)
    .def("__getitem__", &RRectangleListSet::get, return_value_policy<copy_const_reference>())
    .def("__setitem__", &RRectangleListSet::set)
    .def(str(self))    // __str__
  ;
}
