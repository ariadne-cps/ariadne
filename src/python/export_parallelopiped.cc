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

#include "numerical_type.h"
#include "rectangle.h"
#include "parallelopiped.h"
#include "list_set.h"

#include <boost/python.hpp>

#include "real_typedef.h"
#include "python_utilities.h"

typedef Ariadne::LinearAlgebra::matrix<Real> RMatrix;
typedef Ariadne::Geometry::Point<Real> RPoint;
typedef Ariadne::Interval<Real> RInterval;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::Parallelopiped<Real> RParallelopiped;
typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Parallelopiped> RParallelopipedListSet;

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

void export_parallelopiped() {
  typedef bool (*PlpdPlpdBinPred) (const RParallelopiped&, const RParallelopiped&);
  typedef bool (*PlpdRectBinPred) (const RParallelopiped&, const RRectangle&);
  typedef bool (*RectPlpdBinPred) (const RRectangle&, const RParallelopiped&);
  PlpdPlpdBinPred plpd_plpd_interiors_intersect=&interiors_intersect<Real>;
  PlpdRectBinPred plpd_rect_interiors_intersect=&interiors_intersect<Real>;
  RectPlpdBinPred rect_plpd_interiors_intersect=&interiors_intersect<Real>;
  PlpdPlpdBinPred plpd_plpd_disjoint=&disjoint<Real>;
  PlpdRectBinPred plpd_rect_disjoint=&disjoint<Real>;
  RectPlpdBinPred rect_plpd_disjoint=&disjoint<Real>;
  PlpdPlpdBinPred plpd_plpd_inner_subset=&inner_subset<Real>;
  PlpdRectBinPred plpd_rect_inner_subset=&inner_subset<Real>;
  RectPlpdBinPred rect_plpd_inner_subset=&inner_subset<Real>;
  PlpdPlpdBinPred plpd_plpd_subset=&subset<Real>;

  def("interiors_intersect", plpd_plpd_interiors_intersect);
  def("interiors_intersect", plpd_rect_interiors_intersect);
  def("interiors_intersect", rect_plpd_interiors_intersect);
  def("disjoint", plpd_plpd_disjoint);
  def("disjoint", plpd_rect_disjoint);
  def("disjoint", rect_plpd_disjoint);
  def("inner_subset", plpd_plpd_inner_subset);
  def("inner_subset", plpd_rect_inner_subset);
  def("inner_subset", rect_plpd_inner_subset);

  def("subset", plpd_plpd_subset);

  class_<RParallelopiped>("Parallelopiped",init<int>())
    .def(init<RPoint,RMatrix>())
    .def(init<RParallelopiped>())
    .def(init<RRectangle>())
    .def(init<std::string>())
    .def("bounding_box", &RParallelopiped::bounding_box)
    .def("subdivide", &RParallelopiped::subdivide)
    .def("empty", &RParallelopiped::empty)
    .def("empty_interior", &RParallelopiped::empty_interior)
    .def("dimension", &RParallelopiped::dimension)
    .def("contains", &RParallelopiped::contains)
    .def("interior_contains", &RParallelopiped::interior_contains)
    .def(str(self))    // __str__
  ;
}
