/***************************************************************************
 *            python/ariadne.cc
 *
 *  22 June 2005
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
#include "interval.h"
#include "state.h"
#include "rectangle.h"
#include "list_set.h"

#include <boost/python.hpp>

using namespace Ariadne;
using namespace Ariadne::Geometry;

typedef Ariadne::Rational Real;

typedef Ariadne::Geometry::State<Real> RState;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::ListSet<Real,Rectangle> RRectangleListSet;

BOOST_PYTHON_MODULE(geometry)
{
  using boost::python::class_;
  using boost::python::init;
  using boost::python::self;
  using boost::python::return_value_policy;
  using boost::python::copy_const_reference;
  using boost::python::def;
  
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

  typedef bool (*RectLSBinPred) (const RRectangleListSet&, const RRectangleListSet&);
  typedef RRectangleListSet (*RectLSBinFun) (const RRectangleListSet&, const RRectangleListSet&);
  //  RectLSBinPred qrect_interiors_intersect=&interiors_intersect<Real>;
  //RectLSBinPred rectls_disjoint=&disjoint<Real,RRectangle>;
  //RectLSBinPred rectls_inner_subset=&inner_subset<Real,RRectangle>;
  //RectLSBinPred rectls_subset=&subset<Real,RRectangle>;

  //  def("regular_intersection", qrect_regular_intersection);
  //  def("interiors_intersect", qrect_interiors_intersect);
  //def("disjoint", rectls_disjoint);
  //def("inner_subset", rectls_inner_subset);
  //def("subset", rectls_subset);

  class_<RState>("State")
    .def(init<int>())
    .def(init<int,Real>())
    .def(init<RState>())
    .def("dimension", &RState::dimension)
    .def("__len__", &RState::dimension)
    .def("__getitem__", &RState::get)
    .def("__setitem__", &RState::set)
    .def("__eq__", &RState::operator==)
    .def("__ne__", &RState::operator!=)
    .def(boost::python::self_ns::str(self))    // __str__
    ;

  class_< RRectangle >("Rectangle")
    .def(init<int>())
    .def(init<RState,RState>())
    .def(init<RRectangle>())
    .def("dimension", &RRectangle::dimension)
    .def("lower_corner", &RRectangle::lower_corner)
    .def("upper_corner", &RRectangle::upper_corner)
    .def("set_lower", &RRectangle::set_lower)
    .def("set_upper", &RRectangle::set_upper)
    .def("__len__", &RRectangle::dimension)
    .def("__getitem__", &RRectangle::get)
    .def("__setitem__", &RRectangle::set)
    .def("__eq__", &RRectangle::operator==)
    .def("__ne__", &RRectangle::operator!=)
    //.def(boost::python::self_ns::str(self))    // __str__
    ;

  class_<RRectangleListSet>("RectangleListSet")
    .def(init<RRectangle>())
    .def(init<RRectangleListSet>())
    .def("dimension", &RRectangleListSet::dimension)
    .def("push_back", &RRectangleListSet::push_back)
    .def("__len__", &RRectangleListSet::size)
    .def("__getitem__", &RRectangleListSet::get, return_value_policy<copy_const_reference>())
    .def("__setitem__", &RRectangleListSet::set)
    //.def(boost::python::self_ns::str(self))    // __str__
    ;

}
