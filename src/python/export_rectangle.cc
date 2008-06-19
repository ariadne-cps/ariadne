/***************************************************************************
 *            python/export_rectangle.cc
 *
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *
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
#include "geometry/rectangle.h"
#include "geometry/list_set.h"

#include "python/utilities.h"
#include "python/read_box.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>


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
  Box<R> bx;
  read_box(bx,obj);
  return new Rectangle<R>(bx);
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




}

template void export_rectangle<FloatPy>();
