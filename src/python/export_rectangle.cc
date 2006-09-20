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

template<typename R>
void export_rectangle() 
{
  typedef bool (*RectRectBinPred) (const Rectangle<R>&, const Rectangle<R>&);
  typedef Rectangle<R> (*RectRectBinFunc) (const Rectangle<R>&, const Rectangle<R>&);
  typedef Zonotope<R> (*RectZntpBinFunc) (const Rectangle<R>&, const Zonotope<R>&);
  typedef Rectangle<R> (*RectRealBinFunc) (const Rectangle<R>&, const R&);

  typedef bool (Rectangle<R>::*RectPred)(const Rectangle<R> &) const;
  typedef bool (Rectangle<R>::*PointPred) (const Point<R>&) const;

  def("convex_hull", RectRectBinFunc(&rectangular_hull));
  def("regular_intersection", RectRectBinFunc(&regular_intersection));
  def("touching_intersection", RectRectBinFunc(&regular_intersection));
  def("interiors_intersect", RectRectBinPred(&interiors_intersect));
  def("disjoint", RectRectBinPred(&disjoint));
  def("inner_subset", RectRectBinPred(&inner_subset));

  def("intersection", RectRectBinFunc(&intersection));
  def("subset", RectRectBinPred(&subset));
  
  class_< Rectangle<R> >(python_name<R>("Rectangle").c_str(),init<int>())
    .def(init< Point<R>,Point<R> >())
    .def(init< Rectangle<R> >())
    .def(init< Vector<Interval<R> > >())
    .def(init<std::string>())
    .def("empty", & Rectangle<R>::empty)
    .def("empty_interior", &Rectangle<R>::empty_interior)
    .def("dimension", &Rectangle<R>::dimension)
    .def("contains", RectPred(&Rectangle<R>::contains) )
    .def("contains", PointPred(&Rectangle<R>::contains))
    .def("interior_contains", &Rectangle<R>::interior_contains)
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
    .def("__add__", (RectRectBinFunc)(&minkowski_sum))
    .def("__add__", (RectZntpBinFunc)(&minkowski_sum))
    .def("__sub__", (RectRectBinFunc)(&minkowski_difference))
    .def("__sub__", (RectZntpBinFunc)(&minkowski_difference))
    .def("__mul__", (RectRealBinFunc)(&scale))
    .def(self_ns::str(self))
  ;
}

template void export_rectangle<MPFloat>();
