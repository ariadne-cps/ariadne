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
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "python/python_float.h"

#include "linear_algebra/vector.h"

#include "geometry/rectangle.h"
#include "geometry/list_set.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include "python/python_utilities.h"
#include "python/read_scalar.h"

using namespace boost::python;


template<class X>
std::string
__str__(const Rectangle<X>& r)
{
  std::stringstream ss;
  ss << r;
  return ss.str();
}

template<class X>
std::string
__repr__(const Rectangle<X>& r)
{
  std::stringstream ss;
  ss << "Rectangle(";
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
Rectangle<R>* 
make_rectangle(const boost::python::object& obj) 
{
  // See "Extracting C++ objects" in the Boost Python tutorial
  boost::python::list elements=extract<list>(obj);
  int d=boost::python::len(elements);
  Rectangle<R>& r=*new Rectangle<R>(d);
  for(int i=0; i!=d; ++i) {
    extract<list> extract_list(elements[i]);
    if(extract_list.check()) {
      boost::python::list pair=extract_list();
      if(boost::python::len(pair)!=2) {
        throw std::runtime_error("Rectangle must be list of pairs representing intervals");
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

template<class R> inline
Rectangle<R> over_approximation_of_minkowski_sum(const Rectangle<R>& r1, const Rectangle<R>& r2) {
  return over_approximation(minkowski_sum(r1,r2));
}

template<class R> inline
Rectangle<R> over_approximation_of_minkowski_difference(const Rectangle<R>& r1, const Rectangle<R>& r2) {
  return over_approximation(minkowski_difference(r1,r2));
}


template<class R>
void export_rectangle() 
{
  typedef Interval<R> I;
  
  class_< Rectangle<R> >("Rectangle",init<int>())
    .def("__init__", make_constructor(&make_rectangle<R>) )
    .def(init< int >())
    .def(init< Point<R>,Point<R> >())
    .def(init< Rectangle<R> >())
    .def(init< Vector<I> >())
    .def(init< Point<I> >())
    //.def(init<std::string>())
    .def("empty", &Rectangle<R>::empty)
    .def("dimension", &Rectangle<R>::dimension)
    .def("contains", &Rectangle<R>::contains)
    .def("centre", &Rectangle<R>::centre)
    .def("radius", &Rectangle<R>::radius)
    .def("__getitem__", &Rectangle<R>::interval, return_value_policy<copy_const_reference>())
    .def("__setitem__", &Rectangle<R>::set_interval)
    .def("set_lower_bound", &Rectangle<R>::set_lower_bound)
    .def("set_upper_bound", &Rectangle<R>::set_upper_bound)
    .def("bounding_box", &Rectangle<R>::bounding_box)
    .def("lower_corner", &Rectangle<R>::lower_corner)
    .def("upper_corner", &Rectangle<R>::upper_corner)
    .def("lower_bound", (const R&(Rectangle<R>::*)(dimension_type)const)(&Rectangle<R>::lower_bound), return_value_policy<copy_const_reference>())
    .def("upper_bound", (const R&(Rectangle<R>::*)(dimension_type)const)(&Rectangle<R>::upper_bound), return_value_policy<copy_const_reference>())
    .def("neighbourhood", &Rectangle<R>::neighbourhood)
    .def("__add__", &over_approximation_of_minkowski_sum<R>)
    .def("__sub__", &over_approximation_of_minkowski_difference<R>)
    .def("__str__",&__str__<R>)
    .def("__repr__",&__repr__<R>)
  ;

  def("rectangular_hull", (Rectangle<R>(*)(const Rectangle<R>&, const Rectangle<R>&))(&rectangular_hull));
  def("open_intersection", (Rectangle<R>(*)(const Rectangle<R>&, const Rectangle<R>&))(&open_intersection));
  def("closed_intersection", (Rectangle<R>(*)(const Rectangle<R>&, const Rectangle<R>&))(&closed_intersection));
  def("disjoint", (tribool(*)(const Rectangle<R>&, const Rectangle<R>&))(&disjoint));
  def("subset", (tribool(*)(const Rectangle<R>&, const Rectangle<R>&))(&subset));



  class_< Rectangle<I> >("FuzzyRectangle",init<int>())
    .def(init< Point<I>, Point<I> >())
    .def(init< Rectangle<R> >())
    .def(init< Rectangle<I> >())
    //.def(init<std::string>())
    .def("dimension", &Rectangle<I>::dimension)
    .def("contains", &Rectangle<I>::contains)
    .def("set_lower_bound", &Rectangle<I>::set_lower_bound)
    .def("set_upper_bound", &Rectangle<I>::set_upper_bound)
    .def("bounding_box", &Rectangle<I>::bounding_box)
    .def("lower_corner", &Rectangle<I>::lower_corner)
    .def("upper_corner", &Rectangle<I>::upper_corner)
    .def("lower_bound", &Rectangle<I>::lower_bound, return_value_policy<copy_const_reference>())
    .def("upper_bound", &Rectangle<I>::upper_bound, return_value_policy<copy_const_reference>())
    .def("__str__",&__str__<I>)
    .def("__repr__",&__repr__<I>)
  ;

  def("under_approximation", (Rectangle<R>(*)(const Rectangle<I>&))(&under_approximation));
  def("over_approximation", (Rectangle<R>(*)(const Rectangle<I>&))(&over_approximation));


}

template void export_rectangle<Float>();
