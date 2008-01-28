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


template<class X>
std::string
__str__(const Polytope<X>& p)
{
  std::stringstream ss;
  for(size_type i=0; i!=p.number_of_vertices(); ++i) {
    ss << (i==0?"[":",");
    Point<X> v=p.vertex(i);
    for(size_type j=0; j!=v.dimension(); ++j) {
      ss << (j==0?"(":",") << v[j].get_d();
    }
    ss << ")";
  }
  ss << "]";
  return ss.str();
}

template<class X>
std::string
__repr__(const Polytope<X>& p)
{
  std::stringstream ss;
  ss << p;
  return ss.str();
}


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

  class_< Polytope<R> > polytope_class(python_name<R>("Polytope").c_str(),init<int>());
  polytope_class.def("__init__", make_constructor(&make_polytope<R>) );
  polytope_class.def(init< PointList<R> >());
  polytope_class.def(init< Polytope<R> >());
  polytope_class.def(init< Box<R> >());
  polytope_class.def("dimension", &Polytope<R>::dimension);
  polytope_class.def("number_of_vertices", &Polytope<R>::number_of_vertices);
  polytope_class.def("empty", &Polytope<R>::empty);
  polytope_class.def("vertices", &Polytope<R>::vertices);
  polytope_class.def("bounding_box", &Polytope<R>::bounding_box);
  polytope_class.def(self_ns::str(self));

}

template<>
void export_polytope<Rational>() 
{
  typedef Rational Q;

  def("polytope", (Polytope<Q>(*)(const Polyhedron<Q>&))(&polytope));
  def("disjoint", (tribool(*)(const Polytope<Q>&, const Polytope<Q>&))(&disjoint));
  def("subset", (tribool(*)(const Polytope<Q>&, const Polytope<Q>&))(&subset));
  def("convex_hull", (Polytope<Q>(*)(const Polytope<Q>&, const Polytope<Q>&))(&convex_hull));

  class_< PointList<Q> > point_list_class(python_name<Q>("PointList").c_str(),init<int>());
  point_list_class.def("__getitem__", (Point<Q>(PointList<Q>::*)(size_type)const)&PointList<Q>::operator[]);
  point_list_class.def(self_ns::str(self));

  class_< Polytope<Q> > polytope_class(python_name<Q>("Polytope").c_str(),init<int>());
  polytope_class.def("__init__", make_constructor(&make_polytope<Q>) );
  polytope_class.def(init< Polyhedron<Q> >());
  polytope_class.def(init< Polytope<Q> >());
  polytope_class.def(init< Box<Q> >());
  polytope_class.def("dimension", &Polytope<Q>::dimension);
  polytope_class.def("number_of_vertices", &Polytope<Q>::number_of_vertices);
  polytope_class.def("empty", &Polytope<Q>::empty);
  polytope_class.def("vertices", &Polytope<Q>::vertices);
  polytope_class.def("vertex", &Polytope<Q>::vertex);
  polytope_class.def("bounding_box", &Polytope<Q>::bounding_box);
  polytope_class.def("__str__", &__str__<Q>);
  polytope_class.def("__repr__", &__repr__<Q>);

}






template void export_polytope<FloatPy>();
template void export_polytope<Rational>();
