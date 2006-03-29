/***************************************************************************
 *            python/export_parallelotope.cc
 *
 *  22 March 2006
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"

#include <boost/python.hpp>

#include "python/real_typedef.h"
#include "python/python_utilities.h"

typedef Ariadne::LinearAlgebra::matrix<Real> RMatrix;
typedef Ariadne::Geometry::Point<Real> RPoint;
typedef Ariadne::Interval<Real> RInterval;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::Parallelotope<Real> RParallelotope;
typedef Ariadne::Geometry::Zonotope<Real> RZonotope;
typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Parallelotope> RParallelotopeListSet;

using Ariadne::Geometry::interiors_intersect;
using Ariadne::Geometry::disjoint;
using Ariadne::Geometry::inner_subset;
using Ariadne::Geometry::subset;
using Ariadne::Geometry::minkowski_sum;

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

void export_zonotope() {
  typedef RZonotope (*ZntpZntpBinFunc) (const RZonotope&, const RZonotope&);
  typedef bool (*ZntpZntpBinPred) (const RZonotope&, const RZonotope&);
  typedef bool (*ZntpRectBinPred) (const RZonotope&, const RRectangle&);
  typedef bool (*RectZntpBinPred) (const RRectangle&, const RZonotope&);
  typedef bool (*PltpZntpBinPred) (const RParallelotope&, const RZonotope&);
  typedef bool (*ZntpPltpBinPred) (const RZonotope&, const RParallelotope&);
 
  def("interiors_intersect", ZntpZntpBinPred(&interiors_intersect));
  def("interiors_intersect", ZntpRectBinPred(&interiors_intersect));
  def("interiors_intersect", RectZntpBinPred(&interiors_intersect));
  def("interiors_intersect", PltpZntpBinPred(&interiors_intersect));
  def("interiors_intersect", ZntpPltpBinPred(&interiors_intersect));
  def("disjoint", ZntpZntpBinPred(&disjoint));
  def("disjoint", ZntpRectBinPred(&disjoint));
  def("disjoint", RectZntpBinPred(&disjoint));
  def("disjoint", PltpZntpBinPred(&disjoint));
  def("disjoint", ZntpPltpBinPred(&disjoint));
  def("inner_subset", ZntpZntpBinPred(&inner_subset));
  def("inner_subset", ZntpRectBinPred(&inner_subset));
  def("inner_subset", RectZntpBinPred(&inner_subset));
  def("inner_subset", PltpZntpBinPred(&inner_subset));
  def("inner_subset", ZntpPltpBinPred(&inner_subset));
  def("subset", ZntpZntpBinPred(&subset));
  def("subset", ZntpRectBinPred(&subset));
  def("subset", RectZntpBinPred(&subset));
  def("subset", PltpZntpBinPred(&subset));
  def("subset", ZntpPltpBinPred(&subset));
  def("minkowski_sum", ZntpZntpBinFunc(&minkowski_sum));

  class_<RZonotope>("Zonotope",init<int>())
    .def(init<RPoint,RMatrix>())
    .def(init<RZonotope>())
    .def(init<RParallelotope>())
    .def(init<RRectangle>())
    .def(init<std::string>())
    .def("bounding_box", &RZonotope::bounding_box)
    //.def("subdivide", &RZonotope::subdivide)
    .def("empty", &RZonotope::empty)
    .def("empty_interior", &RZonotope::empty_interior)
    .def("dimension", &RZonotope::dimension)
    .def("contains", &RZonotope::contains)
    //.def("==", &RZonotope::operator==)
    //.def("!=", &RZonotope::operator!=)
    .def("interior_contains", &RZonotope::interior_contains)
    .def(str(self))    // __str__
  ;
}
