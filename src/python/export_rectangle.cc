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
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */



#include "geometry/rectangle.h"

#include "geometry/zonotope.h"
#include "geometry/list_set.h"

#include "python/typedefs.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

void export_rectangle() {
  typedef bool (*RectBinPred) (const RRectangle&, const RRectangle&);
  typedef RRectangle (*RectBinFunc) (const RRectangle&, const RRectangle&);
  typedef RZonotope (*RectZntpBinFunc) (const RRectangle&, const RZonotope&);

  typedef bool (RRectangle::*RectPred)(const RRectangle &) const;
  typedef bool (RRectangle::*PointPred) (const RPoint&) const;

  def("convex_hull", RectBinFunc(&rectangular_hull));
  def("regular_intersection", RectBinFunc(&regular_intersection));
  def("touching_intersection", RectBinFunc(&regular_intersection));
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
    .def("contains", RectPred(&RRectangle::contains) )
    .def("contains", PointPred(&RRectangle::contains))
    .def("interior_contains", &RRectangle::interior_contains)
    .def("centre", &RRectangle::centre)
    .def("radius", &RRectangle::radius)
    .def("__getitem__", &RRectangle::interval)
    .def("__setitem__", &RRectangle::set_interval)
    .def("set_lower_bound", &RRectangle::set_lower_bound)
    .def("set_upper_bound", &RRectangle::set_upper_bound)
    .def("bounding_box", &RRectangle::bounding_box)
    .def("lower_corner", &RRectangle::lower_corner)
    .def("upper_corner", &RRectangle::upper_corner)
    .def("lower_bound", &RRectangle::lower_bound, return_value_policy<copy_const_reference>())
    .def("upper_bound", &RRectangle::upper_bound, return_value_policy<copy_const_reference>())
    .def("__add__", (RectBinFunc)(&minkowski_sum))
    .def("__add__", (RectZntpBinFunc)(&minkowski_sum))
    .def("__sub__", (RectBinFunc)(&minkowski_difference))
    .def("__sub__", (RectZntpBinFunc)(&minkowski_difference))
    .def(self_ns::str(self))    // __self_ns::str__
  ;
}
