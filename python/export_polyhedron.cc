/***************************************************************************
 *            python/export_polyhedron.cc
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
#include "polyhedron.h"

#include <boost/python.hpp>

#include "real_typedef.h"
#include "python_utilities.h"

typedef Ariadne::Geometry::State<Real> RState;
typedef Ariadne::Geometry::Polyhedron<Real> RPolyhedron;

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

void export_polyhedron() {
  typedef bool (*PolyBinPred) (const RPolyhedron&, const RPolyhedron&);
  typedef RPolyhedron (*PolyBinFunc) (const RPolyhedron&, const RPolyhedron&);
  PolyBinFunc poly_regular_intersection=&regular_intersection<Real>;
  PolyBinPred poly_interiors_intersect=&interiors_intersect<Real>;
  PolyBinPred poly_disjoint=&disjoint<Real>;
  PolyBinPred poly_inner_subset=&inner_subset<Real>;
  PolyBinPred poly_subset=&subset<Real>;

  def("regular_intersection", poly_regular_intersection);
  def("interiors_intersect", poly_interiors_intersect);
  def("disjoint", poly_disjoint);
  def("inner_subset", poly_inner_subset);
  def("subset", poly_subset);

  class_<RPolyhedron>("Polyhedron",init<int>())
    .def(init<RPolyhedron>())
    .def("dimension", &RPolyhedron::dimension)
    .def(str(self))    // __str__
  ;
}
