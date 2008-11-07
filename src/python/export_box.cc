/***************************************************************************
 *            python/export_box.cc
 *
 *
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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
#include "geometry/grid_cell.h"
#include "geometry/list_set.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include "python/utilities.h"
#include "python/read_box.h"

using namespace Ariadne;
using namespace Ariadne::Python;

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
  Box<R>* bx=new Box<R>;
  read_box(*bx,obj);
  return bx;
}




template<class R>
void export_box() 
{
  typedef Interval<R> I;
  
  class_< Box<R> > box_class(python_name<R>("Box").c_str(),no_init);
  box_class.def("__init__", make_constructor(&make_box<R>) );
  box_class.def(init< int >());
  box_class.def(init< Point<R> >());
  box_class.def(init< Point<R>,Point<R> >());
  box_class.def(init< Box<R> >());
  box_class.def(init< GridCell<R> >());
  box_class.def(init< Vector<I> >());
  box_class.def(init< Point<I> >());
    //box_class.def(init<std::string>());
  box_class.def("dimension", &Box<R>::dimension);
  box_class.def("contains", (tribool(*)(const Box<R>&,const Point<R>&)) &contains);
  box_class.def("centre", &Box<R>::centre);
  box_class.def("radius", &Box<R>::radius);
  box_class.def("__getitem__", &Box<R>::interval, return_value_policy<copy_const_reference>());
  box_class.def("__setitem__", &Box<R>::set_interval);
  box_class.def("set_lower_bound", &Box<R>::set_lower_bound);
  box_class.def("set_upper_bound", &Box<R>::set_upper_bound);
  box_class.def("lower_corner", &Box<R>::lower_corner);
  box_class.def("upper_corner", &Box<R>::upper_corner);
  box_class.def("lower_bound", (const R&(Box<R>::*)(dimension_type)const)(&Box<R>::lower_bound), return_value_policy<copy_const_reference>());
  box_class.def("upper_bound", (const R&(Box<R>::*)(dimension_type)const)(&Box<R>::upper_bound), return_value_policy<copy_const_reference>());
  box_class.def("neighbourhood", &Box<R>::neighbourhood);
  box_class.def("__str__",&__str__<R>);
  box_class.def("__repr__",&__repr__<R>);

  def("rectangular_hull", (Box<R>(*)(const Box<R>&, const Box<R>&))(&rectangular_hull));
  def("open_intersection", (Box<R>(*)(const Box<R>&, const Box<R>&))(&open_intersection));
  def("closed_intersection", (Box<R>(*)(const Box<R>&, const Box<R>&))(&closed_intersection));
  def("contains", (tribool(*)(const Box<R>&, const Point<R>&))(&contains));
  def("disjoint", (tribool(*)(const Box<R>&, const Box<R>&))(&disjoint));
  def("subset", (tribool(*)(const Box<R>&, const Box<R>&))(&subset));

}

template<>
void export_box<Rational>() 
{
  typedef Rational R;
  typedef Interval<R> I;
  
  class_< Box<R> > box_class(python_name<R>("Box").c_str(),no_init);
  box_class.def("__init__", make_constructor(&make_box<R>) );
  box_class.def(init< int >());
  box_class.def(init< Point<R> >());
  box_class.def(init< Point<R>,Point<R> >());
  box_class.def(init< Box<R> >());
  box_class.def(init< Vector<I> >());
  box_class.def(init< Point<I> >());
    //box_class.def(init<std::string>());
  box_class.def("dimension", &Box<R>::dimension);
  box_class.def("contains", (tribool(*)(const Box<R>&,const Point<R>&)) &contains);
  box_class.def("centre", &Box<R>::centre);
  box_class.def("radius", &Box<R>::radius);
  box_class.def("__getitem__", &Box<R>::interval, return_value_policy<copy_const_reference>());
  box_class.def("__setitem__", &Box<R>::set_interval);
  box_class.def("set_lower_bound", &Box<R>::set_lower_bound);
  box_class.def("set_upper_bound", &Box<R>::set_upper_bound);
  box_class.def("lower_corner", &Box<R>::lower_corner);
  box_class.def("upper_corner", &Box<R>::upper_corner);
  box_class.def("lower_bound", (const R&(Box<R>::*)(dimension_type)const)(&Box<R>::lower_bound), return_value_policy<copy_const_reference>());
  box_class.def("upper_bound", (const R&(Box<R>::*)(dimension_type)const)(&Box<R>::upper_bound), return_value_policy<copy_const_reference>());
  box_class.def("neighbourhood", &Box<R>::neighbourhood);
  box_class.def("__str__",&__str__<R>);
  box_class.def("__repr__",&__repr__<R>);

  def("rectangular_hull", (Box<R>(*)(const Box<R>&, const Box<R>&))(&rectangular_hull));
  def("open_intersection", (Box<R>(*)(const Box<R>&, const Box<R>&))(&open_intersection));
  def("closed_intersection", (Box<R>(*)(const Box<R>&, const Box<R>&))(&closed_intersection));
  def("contains", (tribool(*)(const Box<R>&, const Point<R>&))(&contains));
  def("disjoint", (tribool(*)(const Box<R>&, const Box<R>&))(&disjoint));
  def("subset", (tribool(*)(const Box<R>&, const Box<R>&))(&subset));

}

template void export_box<FloatPy>();
template void export_box<Rational>();
