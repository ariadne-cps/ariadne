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
 
  ZntpZntpBinFunc zntp_zntp_minkowski_sum=&minkowski_sum<Real>;
  
  ZntpZntpBinPred zntp_zntp_interiors_intersect=&interiors_intersect<Real>;
  ZntpRectBinPred zntp_rect_interiors_intersect=&interiors_intersect<Real>;
  RectZntpBinPred rect_zntp_interiors_intersect=&interiors_intersect<Real>;
  ZntpPltpBinPred zntp_pltp_interiors_intersect=&interiors_intersect<Real>;
  PltpZntpBinPred pltp_zntp_interiors_intersect=&interiors_intersect<Real>;
  
  ZntpZntpBinPred zntp_zntp_disjoint=&disjoint<Real>;
  ZntpRectBinPred zntp_rect_disjoint=&disjoint<Real>;
  RectZntpBinPred rect_zntp_disjoint=&disjoint<Real>;
  ZntpPltpBinPred zntp_pltp_disjoint=&disjoint<Real>;
  PltpZntpBinPred pltp_zntp_disjoint=&disjoint<Real>;

  ZntpZntpBinPred zntp_zntp_inner_subset=&inner_subset<Real>;
  ZntpRectBinPred zntp_rect_inner_subset=&inner_subset<Real>;
  RectZntpBinPred rect_zntp_inner_subset=&inner_subset<Real>;
  ZntpPltpBinPred zntp_pltp_inner_subset=&inner_subset<Real>;
  PltpZntpBinPred pltp_zntp_inner_subset=&inner_subset<Real>;
  
  ZntpZntpBinPred zntp_zntp_subset=&subset<Real>;
  ZntpRectBinPred zntp_rect_subset=&subset<Real>;
  RectZntpBinPred rect_zntp_subset=&subset<Real>;
  ZntpPltpBinPred zntp_pltp_subset=&subset<Real>;
  PltpZntpBinPred pltp_zntp_subset=&subset<Real>;
  
  def("interiors_intersect", zntp_zntp_interiors_intersect);
  def("interiors_intersect", zntp_rect_interiors_intersect);
  def("interiors_intersect", rect_zntp_interiors_intersect);
  def("interiors_intersect", pltp_zntp_interiors_intersect);
  def("interiors_intersect", zntp_pltp_interiors_intersect);
  def("disjoint", zntp_zntp_disjoint);
  def("disjoint", zntp_rect_disjoint);
  def("disjoint", rect_zntp_disjoint);
  def("disjoint", pltp_zntp_disjoint);
  def("disjoint", zntp_pltp_disjoint);
  def("inner_subset", zntp_zntp_inner_subset);
  def("inner_subset", zntp_rect_inner_subset);
  def("inner_subset", rect_zntp_inner_subset);
  def("inner_subset", pltp_zntp_inner_subset);
  def("inner_subset", zntp_pltp_inner_subset);
  def("subset", zntp_zntp_subset);
  def("subset", zntp_rect_subset);
  def("subset", rect_zntp_subset);
  def("subset", pltp_zntp_subset);
  def("subset", zntp_pltp_subset);
  def("minkowski_sum", zntp_zntp_minkowski_sum);

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
