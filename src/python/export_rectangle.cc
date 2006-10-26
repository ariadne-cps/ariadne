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

#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"

#include "linear_algebra/vector.h"

#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
#include "python/python_utilities.h"
using namespace boost::python;


template<typename R> inline
Rectangle<R> over_approximation_of_minkowski_sum(const Rectangle<R>& r1, const Rectangle<R>& r2) {
  return over_approximation(minkowski_sum(r1,r2));
}

template<typename R> inline
Rectangle<R> over_approximation_of_minkowski_difference(const Rectangle<R>& r1, const Rectangle<R>& r2) {
  return over_approximation(minkowski_difference(r1,r2));
}


template<typename R>
void export_rectangle() 
{
  typedef Point<R> RPoint;
  typedef Rectangle<R> RRectangle;
  typedef Zonotope<R> RZonotope;
    
  def("rectangular_hull", (RRectangle(*)(const RRectangle&, const RRectangle&))(&rectangular_hull));
  def("open_intersection", (RRectangle(*)(const RRectangle&, const RRectangle&))(&open_intersection));
  def("closed_intersection", (RRectangle(*)(const RRectangle&, const RRectangle&))(&closed_intersection));
  def("disjoint", (tribool(*)(const RRectangle&, const RRectangle&))(&disjoint));
  def("subset", (tribool(*)(const RRectangle&, const RRectangle&))(&subset));

  
  class_<RRectangle>(python_name<R>("Rectangle").c_str(),init<int>())
    .def(init< Point<R>,Point<R> >())
    .def(init<RRectangle>())
    .def(init< Vector<Interval<R> > >())
    .def(init<std::string>())
    .def("empty", &RRectangle::empty)
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
    .def("lower_bound", &Rectangle<R>::lower_bound, return_value_policy<copy_const_reference>())
    .def("upper_bound", &Rectangle<R>::upper_bound, return_value_policy<copy_const_reference>())
    .def("__add__", &over_approximation_of_minkowski_sum<R>)
    .def("__sub__", &over_approximation_of_minkowski_difference<R>)
//    .def("__add__", (RRectangle(*)(const RRectangle&,const RRectangle&))(&minkowski_sum))
//    .def("__add__", (RRectangle(*)(const RRectangle&,const RZonotope&))(&minkowski_sum))
//    .def("__sub__", (RRectangle(*)(const RRectangle&,const RRectangle&))(&minkowski_difference))
//    .def("__sub__", (RRectangle(*)(const RRectangle&,const RZonotope&))(&minkowski_difference))
//    .def("__mul__", (RRectangle(*)(const RRectangle&,const R&))(&scale))
    .def(self_ns::str(self))
  ;
}

template void export_rectangle<MPFloat>();
