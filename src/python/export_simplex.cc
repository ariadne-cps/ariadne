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
#include "geometry/simplex.h"
#include "geometry/list_set.h"

#include <boost/python.hpp>

#include "python/real_typedef.h"

typedef Ariadne::LinearAlgebra::matrix<Real> RMatrix;
typedef Ariadne::Geometry::Point<Real> RPoint;
typedef Ariadne::Interval<Real> RInterval;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::Simplex<Real> RSimplex;
typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Simplex> RSimplexListSet;

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

void export_simplex() {
  typedef bool (*SmplxBinPred) (const RSimplex&, const RSimplex&);
  typedef bool (*SmplxRectBinPred) (const RSimplex&, const RRectangle&);
  SmplxRectBinPred smplx_rect_interiors_intersect=&interiors_intersect<Real>;
  SmplxBinPred smplx_interiors_intersect=&interiors_intersect<Real>;
  SmplxBinPred smplx_disjoint=&disjoint<Real>;
  SmplxBinPred smplx_inner_subset=&inner_subset<Real>;
  SmplxBinPred smplx_subset=&subset<Real>;

  def("interiors_intersect", smplx_interiors_intersect);
  def("interiors_intersect", smplx_rect_interiors_intersect);
  def("disjoint", smplx_disjoint);
  def("inner_subset", smplx_inner_subset);

  def("subset", smplx_subset);

  class_<RSimplex>("Simplex",init<int>())
    .def(init< std::vector<RPoint> >())
    .def(init< Ariadne::array<RPoint> >())
    .def(init<RSimplex>())
    .def(init<std::string>())
    .def("empty", &RSimplex::empty)
    .def("empty_interior", &RSimplex::empty_interior)
    .def("dimension", &RSimplex::dimension)
    .def("contains", &RSimplex::contains)
    .def("interior_contains", &RSimplex::interior_contains)
    .def(str(self))    // __str__
  ;
}
