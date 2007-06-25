/***************************************************************************
 *            python/export_point.cc
 *
 *  21 October 2005
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

#include "python/python_float.h"

#include "linear_algebra/vector.h"

#include "geometry/point.h"
#include "geometry/point_list.h"
#include "geometry/rectangle.h"

#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

template<class R1, class R2> inline 
void point_setitem(Point<R1>& p, uint i, R2 x) {
  p.at(i)=x;
}

template<class R> inline 
R point_getitem(const Point<R>& p, uint i) {
  return p.at(i);
}

//inline RPoint rpoint_add_rVector(const RPoint& p, const RVector& v) {
//  return p+v;
//}

//inline RVector rpoint_sub_rpoint(const RPoint& p, const RPoint& q) {
//  return p-q;
//}

template<class R> inline 
Rectangle<R> rectangle_expanded(const Point<R>& p, const R& r) {
  return Rectangle<R>(p).expand_by(r);
}

template<class R> inline 
Point<R> point_list_get(const PointList<R>& pl, const size_type& n) {
  return pl[n];
}


template<class R>
void export_point() 
{
  typedef typename Numeric::traits<R>::arithmetic_type A;

  class_< Point<R> >(python_name<R>("Point").c_str(),init<>())
    .def(init< int >())
    .def(init< Point<R> >())
    .def(init<std::string>())
    .def("dimension", &Point<R>::dimension)
    .def("__len__", &Point<R>::dimension)
    .def("__getitem__", &point_getitem<R>)
    .def("__setitem__", &point_setitem<R,R>)
    .def("__setitem__", &point_setitem<R,double>)
    .def("__eq__", &Point<R>::operator==)
    .def("__ne__", &Point<R>::operator!=)
    .def("__add__",(Point<A>(*)(const Point<R>&,const Vector<R>&))&operator+)
    .def("__sub__",(Point<A>(*)(const Point<R>&,const Vector<R>&))&operator-)
    .def("__sub__",(Vector<A>(*)(const Point<R>&,const Point<R>&))&operator-)
    .def(self_ns::str(self))
  ;
}

template<class R>
void export_interval_point() 
{
  typedef Interval<R> I;

  class_< Point< Interval<R> > >(python_name<R>("IntervalPoint").c_str(),init<>())
    .def(init< int >())
    .def(init< Point<R> >())
    .def(init< Point< Interval<R> > >())
    .def(init< Rectangle<R> >())
    .def(init<std::string>())
    .def("dimension", &Point< Interval<R> >::dimension)
    .def("__len__", &Point< Interval<R> >::dimension)
    .def("__getitem__", &point_getitem< Interval<R> >)
    .def("__setitem__", &point_setitem< Interval<R> , Interval<R> >)
    .def("__setitem__", &point_setitem< Interval<R> , double>)
    .def("__eq__", &Point< Interval<R> >::operator==)
    .def("__ne__", &Point< Interval<R> >::operator!=)
    .def("__add__",(Point<I>(*)(const Point<R>&,const Vector<I>&))&operator+)
    .def("__add__",(Point<I>(*)(const Point<I>&,const Vector<R>&))&operator+)
    .def("__add__",(Point<I>(*)(const Point<I>&,const Vector<I>&))&operator+)
    .def("__sub__",(Point<I>(*)(const Point<R>&,const Vector<I>&))&operator-)
    .def("__sub__",(Point<I>(*)(const Point<I>&,const Vector<R>&))&operator-)
    .def("__sub__",(Point<I>(*)(const Point<I>&,const Vector<I>&))&operator-)
    .def("__sub__",(Vector<I>(*)(const Point<I>&,const Point<I>&))&operator-)
    .def(self_ns::str(self))
  ;
}

template<class R>
void export_point_list() 
{
  class_< PointList<R> >("PointList",init<>())
    .def("size", &PointList<R>::size)
    .def("append", &PointList<R>::push_back)
    .def("__getitem__", &point_list_get<R>)
  ;
}




template void export_point<Float>();
template void export_interval_point<Float>();
template void export_point_list<Float>();
