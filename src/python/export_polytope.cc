/***************************************************************************
 *            python/export_polytope.cc
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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
#include "python/utilities.h"
#include "python/read_matrix.h"

#include "linear_algebra/matrix.h"

#include "geometry/box.h"
#include "geometry/polytope.h"
#include "geometry/polyhedron.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>  
void
read_polytope(Geometry::Polytope<R>& pltp, const boost::python::object& obj) 
{
  Matrix<R> A;
  read_matrix(A,obj);
  A=transpose(A);
  Matrix<R> G(A.number_of_rows()+1,A.number_of_columns());
  G(slice(0,A.number_of_rows()),slice(0,A.number_of_columns()))=A;
  for(size_type j=0; j!=A.number_of_columns(); ++j) {
    G(A.number_of_rows(),j)=1;
  }
  pltp=Polytope<R>(G);
}

template<class R>
Polytope<R>*
make_polytope(const boost::python::object& obj)
{
  Polytope<R>* pltp=new Polytope<R>;
  read_polytope(*pltp,obj);
  return pltp;
}


template<class R>
void export_polytope() 
{
  def("disjoint", (tribool(*)(const Polytope<R>&, const Polytope<R>&))(&disjoint));
  def("subset", (tribool(*)(const Polytope<R>&, const Polytope<R>&))(&subset));
  def("convex_hull", (Polytope<R>(*)(const Polytope<R>&, const Polytope<R>&))(&convex_hull));

  class_< Polytope<R> >(python_name<R>("Polytope").c_str(),init<int>())
    .def("__init__", make_constructor(&make_polytope<R>) )
    .def(init< PointList<R> >())
    .def(init< Polytope<R> >())
    .def(init< Box<R> >())
    .def("dimension", &Polytope<R>::dimension)
    .def("vertices", &Polytope<R>::vertices)
    .def("bounding_box", &Polytope<R>::bounding_box)
    .def(self_ns::str(self))
  ;

}

template<>
void export_polytope<Rational>() 
{
  typedef Rational Q;

  def("disjoint", (tribool(*)(const Polytope<Q>&, const Polytope<Q>&))(&disjoint));
  def("subset", (tribool(*)(const Polytope<Q>&, const Polytope<Q>&))(&subset));
  def("convex_hull", (Polytope<Q>(*)(const Polytope<Q>&, const Polytope<Q>&))(&convex_hull));

  class_< Polytope<Q> >(python_name<Q>("Polytope").c_str(),init<int>())
    .def("__init__", make_constructor(&make_polytope<Q>) )
    .def(init< Polyhedron<Q> >())
    .def(init< Polytope<Q> >())
    .def(init< Box<Q> >())
    .def("dimension", &Polytope<Q>::dimension)
    .def("vertices", &Polytope<Q>::vertices)
    .def("bounding_box", &Polytope<Q>::bounding_box)
    .def(self_ns::str(self))
  ;

}





template void export_polytope<FloatPy>();
template void export_polytope<Rational>();
