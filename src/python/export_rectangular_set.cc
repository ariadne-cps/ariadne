/***************************************************************************
 *            python/export_rectangular_set.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

#include "geometry/rectangular_set.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;

template<class R>  
RectangularSet<R>* 
make_rectangular_set(boost::python::object obj) 
{
  // See "Extracting C++ objects" in the Boost Python tutorial
  boost::python::list elements=extract<list>(obj);
  int d=boost::python::len(elements);
  Rectangle<R> r(d);
  for(int i=0; i!=d; ++i) {
    extract< Interval<R> > extract_interval(elements[i]);
    if (extract_interval.check()) {
      Interval<R> interval=extract_interval();
      r.set_lower_bound(i,interval.lower());
      r.set_upper_bound(i,interval.upper());
    } else {
      extract<list> extract_list(elements[i]);
      boost::python::list pair=extract_list();
      if(boost::python::len(pair)!=2) {
        throw std::runtime_error("Box must be list of pairs representing intervals");
      }
      extract<double> l(pair[0]);
      if (l.check()) {
        r.set_lower_bound(i,static_cast<R>(l));
      } else {
        extract<R> l(pair[0]);
        r.set_lower_bound(i,static_cast<R>(l));
      }
      extract<double> u(pair[1]);
      if (u.check()) {
        r.set_upper_bound(i,static_cast<R>(u));
      } else {
        extract<R> u(pair[0]);
        r.set_upper_bound(i,static_cast<R>(u));
      }
    }
  }
  return new RectangularSet<R>(r);
}

template<class R>
void export_rectangular_set() 
{
  typedef Numeric::Interval<R> I;
  
  class_< RectangularSet<R>, bases< SetInterface<R>, Rectangle<R> > >("RectangularSet",no_init)
    .def("__init__", make_constructor(&make_rectangular_set<R>) )
    .def(init< Point<I> >())
    .def(init< Rectangle<R> >())
    .def("rectangle", &RectangularSet<R>::operator const Rectangle<R>&,return_value_policy<copy_const_reference>())
    .def("dimension", &RectangularSet<R>::dimension)
    .def("contains", &RectangularSet<R>::contains)
    .def("superset", &RectangularSet<R>::superset)
    .def("intersects", &RectangularSet<R>::intersects)
    .def("disjoint", &RectangularSet<R>::disjoint)
    .def("subset", &RectangularSet<R>::subset)
    .def("bounding_box", &RectangularSet<R>::bounding_box)
    .def(self_ns::str(self))
  ;
}


template void export_rectangular_set<FloatPy>();
