/***************************************************************************
 *            python/export_box.cc
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

#include "geometry/box.h"
#include "geometry/list_set.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include "python/utilities.h"
#include "python/read_scalar.h"

using namespace boost::python;


template<class X>
std::string
__str__(const Box<X>& r)
{
  std::stringstream ss;
  ss << r;
  return ss.str();
}

template<class X>
std::string
__repr__(const Box<X>& r)
{
  std::stringstream ss;
  ss << "Box(";
  if(r.empty()) {
    ss << r.dimension();
  } else {
    for(dimension_type i=0; i!=r.dimension(); ++i) {
      ss << (i==0 ? '[' : ',') << '[' << r.lower_bound(i) << ',' << r.upper_bound(i) << ']';
    }
    ss << ']';
  }
  ss << ")";
  return ss.str();
}

template<class R>  
Box<R>* 
make_box(const boost::python::object& obj) 
{
  // See "Extracting C++ objects" in the Boost Python tutorial
  boost::python::list elements=extract<list>(obj);
  int d=boost::python::len(elements);
  Box<R>& r=*new Box<R>(d);
  for(int i=0; i!=d; ++i) {
    extract<list> extract_list(elements[i]);
    if(extract_list.check()) {
      boost::python::list pair=extract_list();
      if(boost::python::len(pair)!=2) {
        throw std::runtime_error("Box must be list of pairs representing intervals");
      }
      R l=read_scalar<R>(pair[0]);
      R u=read_scalar<R>(pair[1]);
      r.set_lower_bound(i,static_cast<R>(l));
      r.set_upper_bound(i,static_cast<R>(u));
    } else {
      extract<std::string> extract_string(elements[i]);
      if(extract_string.check()) {
        Interval<R> interval(extract_string());
        r.set_lower_bound(i,interval.lower());
        r.set_upper_bound(i,interval.upper());
      } else {
        extract< Interval<R> > extract_interval(elements[i]);
        Interval<R> interval=extract_interval();
        r.set_lower_bound(i,interval.lower());
        r.set_upper_bound(i,interval.upper());
      } 
    }
  }
  return &r;
}



template<class R>
void export_box() 
{
  typedef Interval<R> I;
  
  class_< Box<R> >("Box",init<int>())
    .def("__init__", make_constructor(&make_box<R>) )
    .def(init< int >())
    .def(init< Point<R>,Point<R> >())
    .def(init< Box<R> >())
    .def(init< Vector<I> >())
    .def(init< Point<I> >())
    //.def(init<std::string>())
    .def("empty", &Box<R>::empty)
    .def("dimension", &Box<R>::dimension)
    .def("contains", &Box<R>::contains)
    .def("centre", &Box<R>::centre)
    .def("radius", &Box<R>::radius)
    .def("__getitem__", &Box<R>::interval, return_value_policy<copy_const_reference>())
    .def("__setitem__", &Box<R>::set_interval)
    .def("set_lower_bound", &Box<R>::set_lower_bound)
    .def("set_upper_bound", &Box<R>::set_upper_bound)
    .def("bounding_box", &Box<R>::bounding_box)
    .def("lower_corner", &Box<R>::lower_corner)
    .def("upper_corner", &Box<R>::upper_corner)
    .def("lower_bound", (const R&(Box<R>::*)(dimension_type)const)(&Box<R>::lower_bound), return_value_policy<copy_const_reference>())
    .def("upper_bound", (const R&(Box<R>::*)(dimension_type)const)(&Box<R>::upper_bound), return_value_policy<copy_const_reference>())
    .def("neighbourhood", &Box<R>::neighbourhood)
    .def("__str__",&__str__<R>)
    .def("__repr__",&__repr__<R>)
  ;

  def("rectangular_hull", (Box<R>(*)(const Box<R>&, const Box<R>&))(&rectangular_hull));
  def("open_intersection", (Box<R>(*)(const Box<R>&, const Box<R>&))(&open_intersection));
  def("closed_intersection", (Box<R>(*)(const Box<R>&, const Box<R>&))(&closed_intersection));
  def("disjoint", (tribool(*)(const Box<R>&, const Box<R>&))(&disjoint));
  def("subset", (tribool(*)(const Box<R>&, const Box<R>&))(&subset));
  def("equal", (tribool(*)(const Box<R>&, const Box<R>&))(&equal));





}

template void export_box<FloatPy>();
