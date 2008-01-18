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
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "python/float.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/box.h"
#include "geometry/polytope.h"
#include "geometry/polyhedron.h"
#include "geometry/set_interface.h"

#include "python/utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;


template<class R>
void export_polyhedron() 
{
  def("disjoint", (tribool(*)(const Polyhedron<R>&, const Box<R>&))(&disjoint));
  def("disjoint", (tribool(*)(const Box<R>&, const Polyhedron<R>&))(&disjoint));
  def("disjoint", (tribool(*)(const Polyhedron<R>&, const Polytope<R>&))(&disjoint));
  def("disjoint", (tribool(*)(const Polytope<R>&, const Polyhedron<R>&))(&disjoint));
  def("disjoint", (tribool(*)(const Polyhedron<R>&, const Polyhedron<R>&))(&disjoint));
  def("subset", (tribool(*)(const Box<R>&, const Polyhedron<R>&))(&subset));
  def("subset", (tribool(*)(const Polytope<R>&, const Polyhedron<R>&))(&subset));
  def("subset", (tribool(*)(const Polyhedron<R>&, const Polyhedron<R>&))(&subset));
  def("closed_intersection", (Polyhedron<R>(*)(const Polyhedron<R>&, const Polyhedron<R>&))(&closed_intersection));

  class_< Polyhedron<R> >(python_name<R>("Polyhedron").c_str(),init<int>())
    .def(init< Matrix<R>, Vector<R> >())
    .def(init< Polyhedron<R> >())
    .def(init< Box<R> >())
    .def("dimension", &Polyhedron<R>::dimension)
    .def(self_ns::str(self))
  ;
  
}

template<>
void export_polyhedron<Rational>() 
{
  typedef Rational Q;
  
  def("disjoint", (tribool(*)(const Polyhedron<Q>&, const Polyhedron<Q>&))(&disjoint));
  def("disjoint", (tribool(*)(const Polyhedron<Q>&, const Polytope<Q>&))(&disjoint));
  def("disjoint", (tribool(*)(const Polytope<Q>&, const Polyhedron<Q>&))(&disjoint));
  def("subset", (tribool(*)(const Polytope<Q>&, const Polyhedron<Q>&))(&subset));
  def("subset", (tribool(*)(const Polyhedron<Q>&, const Polyhedron<Q>&))(&subset));
  def("closed_intersection", (Polyhedron<Q>(*)(const Polyhedron<Q>&, const Polyhedron<Q>&))(&closed_intersection));

  class_< Polyhedron<Q> >(python_name<Q>("Polyhedron").c_str(),init<int>())
    .def(init< Matrix<Q>, Vector<Q> >())
    .def(init< Polytope<Q> >())
    .def(init< Polyhedron<Q> >())
    .def(init< Box<Q> >())
    .def("dimension", &Polyhedron<Q>::dimension)
    .def(self_ns::str(self))
  ;
  
}

template void export_polyhedron<FloatPy>();
template void export_polyhedron<Rational>();
