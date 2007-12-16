/***************************************************************************
 *            python/export_interval_set.cc
 *
 *  Copyright  2007 Pieter Collins
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

#include "numeric/interval.h"
#include "geometry/interval_set.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include "python/utilities.h"
#include "python/read_scalar.h"

using namespace boost::python;


template<class X>
std::string
__str__(const Geometry::IntervalSet<X>& ivl)
{
  std::stringstream ss;
  ss << ivl;
  return ss.str();
}

template<class X>
std::string
__repr__(const Geometry::IntervalSet<X>& ivl)
{
  std::stringstream ss;
  if(empty(ivl)) {
    ss << "IntervalSet()";
  } else {
    ss << "IntervalSet(" << ivl.lower_bound() << "," << ivl.upper_bound() << ")";
  }
  return ss.str();
}



template<class R>
void export_interval_set() 
{
  typedef Numeric::Interval<R> I;

  class_< Geometry::IntervalSet<R> > interval_class("IntervalSet",init<>());
  interval_class.def(init< double >());
  interval_class.def(init< R >());
  interval_class.def(init< double, double >());
  interval_class.def(init< R, R >());
  interval_class.def(init< Geometry::IntervalSet<R> >());
  interval_class.def("lower_bound", &Geometry::IntervalSet<R>::lower_bound, return_value_policy<copy_const_reference>());
  interval_class.def("upper_bound", &Geometry::IntervalSet<R>::upper_bound, return_value_policy<copy_const_reference>());
  interval_class.def("__str__",&__str__<R>);
  interval_class.def("__repr__",&__repr__<R>);
  


  def("centre", (R(*)(const Geometry::IntervalSet<R>&))(&Geometry::centre));
  def("radius", (R(*)(const Geometry::IntervalSet<R>&))(&Geometry::radius));
  def("empty", (tribool(*)(const Geometry::IntervalSet<R>&))(&Geometry::empty));
  def("bounded", (tribool(*)(const Geometry::IntervalSet<R>&))(&Geometry::bounded));
  def("contains", (tribool(*)(const Geometry::IntervalSet<R>&,const R&))(&Geometry::contains));
  def("disjoint", (tribool(*)(const Geometry::IntervalSet<R>&, const Geometry::IntervalSet<R>&))(&Geometry::disjoint));
  def("intersect", (tribool(*)(const Geometry::IntervalSet<R>&, const Geometry::IntervalSet<R>&))(&Geometry::intersect));
  def("subset", (tribool(*)(const Geometry::IntervalSet<R>&, const Geometry::IntervalSet<R>&))(&Geometry::subset));
  def("superset", (tribool(*)(const Geometry::IntervalSet<R>&, const Geometry::IntervalSet<R>&))(&Geometry::superset));
  def("open_intersection", (Geometry::IntervalSet<R>(*)(const Geometry::IntervalSet<R>&, const Geometry::IntervalSet<R>&))(&Geometry::open_intersection));
  def("closed_intersection", (Geometry::IntervalSet<R>(*)(const Geometry::IntervalSet<R>&, const Geometry::IntervalSet<R>&))(&Geometry::closed_intersection));
  def("interval_hull", (Geometry::IntervalSet<R>(*)(const Geometry::IntervalSet<R>&, const Geometry::IntervalSet<R>&))(&Geometry::interval_hull));

 
  class_< Geometry::IntervalSet<I> > fuzzy_interval_class("FuzzyIntervalSet",init<>());
  fuzzy_interval_class.def(init< double >());
  fuzzy_interval_class.def(init< R >());
  fuzzy_interval_class.def(init< double, double >());
  fuzzy_interval_class.def(init< R, R >());
  fuzzy_interval_class.def(init< R, I >());
  fuzzy_interval_class.def(init< I, R >());
  fuzzy_interval_class.def(init< I, I >());
  fuzzy_interval_class.def(init< Geometry::IntervalSet<R> >());
  fuzzy_interval_class.def(init< Geometry::IntervalSet<I> >());
  fuzzy_interval_class.def("lower_bound", &Geometry::IntervalSet<I>::lower_bound, return_value_policy<copy_const_reference>());
  fuzzy_interval_class.def("upper_bound", &Geometry::IntervalSet<I>::upper_bound, return_value_policy<copy_const_reference>());
  fuzzy_interval_class.def("__str__",&__str__<I>);
  fuzzy_interval_class.def("__repr__",&__repr__<I>);

  
}

template void export_interval_set<FloatPy>();
