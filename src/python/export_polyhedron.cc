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
#include "python/read_array.h"
using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class X>
std::string
__str__(const Halfspace<X>& hs)
{
  std::stringstream ss;
  ss << hs;
  return ss.str();
}

template<class X>
void
read_halfspace(Halfspace<X>& hs, const boost::python::object& obj1, const boost::python::object& obj2) 
{
  array< X > a;
  X b;
  read_array(a,obj1);
  read_scalar(b,obj2);
  hs=Halfspace<X>(a.size(),a.begin(),b);
}

template<class X>
Halfspace<X>*
make_halfspace(const boost::python::object& obj1, const boost::python::object& obj2)
{
  Halfspace<X>* hs=new Halfspace<X>;
  read_halfspace(*hs,obj1,obj2);
  return hs;
}



template<class X>
void
read_polyhedron(Polyhedron<X>& plhd, const boost::python::object& obj) 
{
  std::vector< Halfspace<X> > hsv;
  boost::python::list elements=boost::python::extract<boost::python::list>(obj);
  int n=boost::python::len(elements);
  hsv.resize(n);
  for(int i=0; i!=n; ++i) {
    boost::python::extract< Halfspace<X> > hx(elements[i]);
    if(hx.check()) {
      hsv[i]=static_cast< Halfspace<X> >(hx());
    } else {
      boost::python::tuple data=boost::python::extract<boost::python::tuple>(elements[i]);
      ARIADNE_ASSERT(len(data)==2);
      read_halfspace(hsv[i], data[0], data[1]);
    }
  }
  plhd=Polyhedron<X>(hsv);
}

template<class X>
Polyhedron<X>*
make_polyhedron(const boost::python::object& obj)
{
  Polyhedron<X>* plhd=new Polyhedron<X>;
  read_polyhedron(*plhd,obj);
  return plhd;
}


template<class X>
std::string
__str__(const Polyhedron<X>& p)
{
  std::stringstream ss;
  dimension_type d=p.dimension();
  size_type nc=p.number_of_constraints();
  const array<X>& data=p.data();
  for(size_type i=0; i!=nc; ++i) {
    ss << ( i==0 ? "[" : "," );
    for(size_type j=0; j!=d; ++j) {
      ss << ( j==0 ? "(" : ",");
      ss << -data[i*(d+1)+j].get_d(); 
    }
    ss << ";" << data[i*(d+1)+d].get_d() << ")";
  }
  ss << "]";
  return ss.str();
}


template<class X>
std::string
__repr__(const Polyhedron<X>& p)
{
  std::stringstream ss;
  ss << p;
  return ss.str();
}


template<class R>
void export_polyhedron() 
{
  class_< Halfspace<R> > halfspace_class(python_name<R>("Halfspace").c_str(),no_init);
  halfspace_class.def("__init__", make_constructor(&make_halfspace<R>) );
  halfspace_class.def("dimension", &Halfspace<R>::dimension);
  halfspace_class.def("complement", &Halfspace<R>::complement);
  halfspace_class.def("contains", (tribool(Halfspace<R>::*)(const Point<R>&)const) &Halfspace<R>::contains);
  halfspace_class.def(self_ns::str(self));

  def("satisfies", (tribool(*)(const Point<R>&, const Halfspace<R>&)) &satisfies);

  class_< Polyhedron<R> > polyhedron_class(python_name<R>("Polyhedron").c_str(),init<int>());
  polyhedron_class.def("__init__", make_constructor(&make_polyhedron<R>) );
  polyhedron_class.def(init< Matrix<R>, Vector<R> >());
  polyhedron_class.def(init< Halfspace<R> >());
  polyhedron_class.def(init< Polyhedron<R> >());
  polyhedron_class.def(init< Box<R> >());
  polyhedron_class.def("dimension", &Polyhedron<R>::dimension);
  polyhedron_class.def("number_of_constraints", &Polyhedron<R>::number_of_constraints);
  polyhedron_class.def("constraint", &Polyhedron<R>::constraint);
  polyhedron_class.def("empty", &Polyhedron<R>::empty);
  polyhedron_class.def("bounded", &Polyhedron<R>::bounded);
  polyhedron_class.def("contains", (tribool(Polyhedron<R>::*)(const Point<R>&)const) &Polyhedron<R>::contains);
  polyhedron_class.def("__str__", (std::string(*)(const Polyhedron<R>&)) &__str__<R>);
  polyhedron_class.def("__repr__", &__repr__<R>);

  def("polyhedron", (Polyhedron<R>(*)(const Halfspace<R>&))(&polyhedron));
  def("polyhedron", (Polyhedron<R>(*)(const Box<R>&))(&polyhedron));
  def("disjoint", (tribool(*)(const Polyhedron<R>&, const Box<R>&))(&disjoint));
  def("disjoint", (tribool(*)(const Box<R>&, const Polyhedron<R>&))(&disjoint));
  def("disjoint", (tribool(*)(const Polyhedron<R>&, const Polytope<R>&))(&disjoint));
  def("disjoint", (tribool(*)(const Polytope<R>&, const Polyhedron<R>&))(&disjoint));
  def("disjoint", (tribool(*)(const Polyhedron<R>&, const Polyhedron<R>&))(&disjoint));
  def("subset", (tribool(*)(const Box<R>&, const Polyhedron<R>&))(&subset));
  def("subset", (tribool(*)(const Polytope<R>&, const Polyhedron<R>&))(&subset));
  def("subset", (tribool(*)(const Polyhedron<R>&, const Polyhedron<R>&))(&subset));
  def("closed_intersection", (Polyhedron<R>(*)(const Polyhedron<R>&, const Polyhedron<R>&))(&closed_intersection));

}

template<>
void export_polyhedron<Rational>() 
{
  typedef Rational Q;
  
  class_< Halfspace<Q> > halfspace_class(python_name<Q>("Halfspace").c_str(),no_init);
  halfspace_class.def("__init__", make_constructor(&make_halfspace<Q>) );
  halfspace_class.def("dimension", &Halfspace<Q>::dimension);
  halfspace_class.def("complement", &Halfspace<Q>::complement);
  halfspace_class.def("contains", (tribool(Halfspace<Q>::*)(const Point<Q>&)const) &Halfspace<Q>::contains);
  halfspace_class.def(self_ns::str(self));

  def("satisfies", (tribool(*)(const Point<Q>&, const Halfspace<Q>&)) &satisfies);


  class_< Polyhedron<Q> >  polyhedron_class(python_name<Q>("Polyhedron").c_str(),init<int>());
  polyhedron_class.def("__init__", make_constructor(&make_polyhedron<Q>) );
  polyhedron_class.def(init< Matrix<Q>, Vector<Q> >());
  polyhedron_class.def(init< Halfspace<Q> >());
  polyhedron_class.def(init< Polytope<Q> >());
  polyhedron_class.def(init< Polyhedron<Q> >());
  polyhedron_class.def(init< Box<Q> >());
  polyhedron_class.def("dimension", &Polyhedron<Q>::dimension);
  polyhedron_class.def("number_of_constraints", &Polyhedron<Q>::number_of_constraints);
  polyhedron_class.def("constraint", &Polyhedron<Q>::constraint);
  polyhedron_class.def("empty", &Polyhedron<Q>::empty);
  polyhedron_class.def("bounded", &Polyhedron<Q>::bounded);
  polyhedron_class.def("contains", (tribool(Polyhedron<Q>::*)(const Point<Q>&)const) &Polyhedron<Q>::contains);
  polyhedron_class.def("__str__", (std::string(*)(const Polyhedron<Q>&))&__str__<Q>);
  polyhedron_class.def("__repr__", &__repr__<Q>);
  
  def("polyhedron", (Polyhedron<Q>(*)(const Halfspace<Q>&))(&polyhedron));
  def("polyhedron", (Polyhedron<Q>(*)(const Box<Q>&))(&polyhedron));
  def("polyhedron", (Polyhedron<Q>(*)(const Polytope<Q>&))(&polyhedron));
  def("disjoint", (tribool(*)(const Polyhedron<Q>&, const Polyhedron<Q>&))(&disjoint));
  def("disjoint", (tribool(*)(const Polyhedron<Q>&, const Polytope<Q>&))(&disjoint));
  def("disjoint", (tribool(*)(const Polytope<Q>&, const Polyhedron<Q>&))(&disjoint));
  def("subset", (tribool(*)(const Polytope<Q>&, const Polyhedron<Q>&))(&subset));
  def("subset", (tribool(*)(const Polyhedron<Q>&, const Polyhedron<Q>&))(&subset));
  def("intersection", (Polyhedron<Q>(*)(const Polyhedron<Q>&, const Polyhedron<Q>&))(&closed_intersection));

  def("closed_intersection", (Polyhedron<Q>(*)(const Polyhedron<Q>&, const Polyhedron<Q>&))(&closed_intersection));

}

template void export_polyhedron<FloatPy>();
template void export_polyhedron<Rational>();
