/***************************************************************************
 *            python/export_parallelotope.cc
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
#include "geometry/parallelotope.h"

#include <boost/python.hpp>

#include "python/real_typedef.h"
#include "python/python_utilities.h"

typedef Ariadne::LinearAlgebra::matrix<Real> RMatrix;
typedef Ariadne::Geometry::Point<Real> RPoint;
typedef Ariadne::Interval<Real> RInterval;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::Parallelotope<Real> RParallelotope;
typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Parallelotope> RParallelotopeListSet;

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

void export_parallelotope() {
  typedef bool (*PltpPltpBinPred) (const RParallelotope&, const RParallelotope&);
  typedef bool (*PltpRectBinPred) (const RParallelotope&, const RRectangle&);
  typedef bool (*RectPltpBinPred) (const RRectangle&, const RParallelotope&);
  PltpPltpBinPred pltp_pltp_interiors_intersect=&interiors_intersect<Real>;
  PltpRectBinPred pltp_rect_interiors_intersect=&interiors_intersect<Real>;
  RectPltpBinPred rect_pltp_interiors_intersect=&interiors_intersect<Real>;
  PltpPltpBinPred pltp_pltp_disjoint=&disjoint<Real>;
  PltpRectBinPred pltp_rect_disjoint=&disjoint<Real>;
  RectPltpBinPred rect_pltp_disjoint=&disjoint<Real>;
  PltpPltpBinPred pltp_pltp_inner_subset=&inner_subset<Real>;
  PltpRectBinPred pltp_rect_inner_subset=&inner_subset<Real>;
  RectPltpBinPred rect_pltp_inner_subset=&inner_subset<Real>;
  PltpPltpBinPred pltp_pltp_subset=&subset<Real>;

  def("interiors_intersect", pltp_pltp_interiors_intersect);
  def("interiors_intersect", pltp_rect_interiors_intersect);
  def("interiors_intersect", rect_pltp_interiors_intersect);
  def("disjoint", pltp_pltp_disjoint);
  def("disjoint", pltp_rect_disjoint);
  def("disjoint", rect_pltp_disjoint);
  def("inner_subset", pltp_pltp_inner_subset);
  def("inner_subset", pltp_rect_inner_subset);
  def("inner_subset", rect_pltp_inner_subset);

  def("subset", pltp_pltp_subset);

  class_<RParallelotope>("Parallelotope",init<int>())
    .def(init<RPoint,RMatrix>())
    .def(init<RParallelotope>())
    .def(init<RRectangle>())
    .def(init<std::string>())
    .def("bounding_box", &RParallelotope::bounding_box)
    .def("subdivide", &RParallelotope::subdivide)
    .def("empty", &RParallelotope::empty)
    .def("empty_interior", &RParallelotope::empty_interior)
    .def("dimension", &RParallelotope::dimension)
    .def("contains", &RParallelotope::contains)
    .def("interior_contains", &RParallelotope::interior_contains)
    .def(str(self))    // __str__
  ;
}
