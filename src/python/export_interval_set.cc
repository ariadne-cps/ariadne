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

#include "numeric/float.h"
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
__str__(const IntervalSet<X>& ivl)
{
  std::stringstream ss;
  ss << ivl;
  return ss.str();
}

template<class X>
std::string
__repr__(const IntervalSet<X>& ivl)
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
  typedef Interval<R> I;

  class_< IntervalSet<R> > interval_class("IntervalSet",init<>());
  interval_class.def(init< double >());
  interval_class.def(init< R >());
  interval_class.def(init< double, double >());
  interval_class.def(init< R, R >());
  interval_class.def(init< IntervalSet<R> >());
  interval_class.def("lower_bound", &IntervalSet<R>::lower_bound, return_value_policy<copy_const_reference>());
  interval_class.def("upper_bound", &IntervalSet<R>::upper_bound, return_value_policy<copy_const_reference>());
  interval_class.def("__str__",&__str__<R>);
  interval_class.def("__repr__",&__repr__<R>);
  


  def("centre", (R(*)(const IntervalSet<R>&))(&centre));
  def("radius", (R(*)(const IntervalSet<R>&))(&radius));
  def("empty", (tribool(*)(const IntervalSet<R>&))(&empty));
  def("bounded", (tribool(*)(const IntervalSet<R>&))(&bounded));
  def("contains", (tribool(*)(const IntervalSet<R>&,const R&))(&contains));
  def("disjoint", (tribool(*)(const IntervalSet<R>&, const IntervalSet<R>&))(&disjoint));
  def("intersect", (tribool(*)(const IntervalSet<R>&, const IntervalSet<R>&))(&intersect));
  def("subset", (tribool(*)(const IntervalSet<R>&, const IntervalSet<R>&))(&subset));
  def("superset", (tribool(*)(const IntervalSet<R>&, const IntervalSet<R>&))(&superset));
  def("open_intersection", (IntervalSet<R>(*)(const IntervalSet<R>&, const IntervalSet<R>&))(&open_intersection));
  def("closed_intersection", (IntervalSet<R>(*)(const IntervalSet<R>&, const IntervalSet<R>&))(&closed_intersection));
  def("interval_hull", (IntervalSet<R>(*)(const IntervalSet<R>&, const IntervalSet<R>&))(&interval_hull));

 
  class_< IntervalSet<I> > fuzzy_interval_class("FuzzyIntervalSet",init<>());
  fuzzy_interval_class.def(init< double >());
  fuzzy_interval_class.def(init< R >());
  fuzzy_interval_class.def(init< double, double >());
  fuzzy_interval_class.def(init< R, R >());
  fuzzy_interval_class.def(init< R, I >());
  fuzzy_interval_class.def(init< I, R >());
  fuzzy_interval_class.def(init< I, I >());
  fuzzy_interval_class.def(init< IntervalSet<R> >());
  fuzzy_interval_class.def(init< IntervalSet<I> >());
  fuzzy_interval_class.def("lower_bound", &IntervalSet<I>::lower_bound, return_value_policy<copy_const_reference>());
  fuzzy_interval_class.def("upper_bound", &IntervalSet<I>::upper_bound, return_value_policy<copy_const_reference>());
  fuzzy_interval_class.def("__str__",&__str__<I>);
  fuzzy_interval_class.def("__repr__",&__repr__<I>);

  
}

template void export_interval_set<FloatPy>();
